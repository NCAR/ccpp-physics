!
!-----------------------------------------------------------------------
!#include "nmm_loop_basemacros.h"
!#include "nmm_loop_macros.h"
!-----------------------------------------------------------------------
!
      MODULE HAFS_radiation
!
      use machine, only : kind_phys
!-----------------------------------------------------------------------
!      USE MODULE_DOMAIN
!      USE MODULE_DM
!      USE MODULE_CONFIGURE
!      USE MODULE_TILES
      USE MODULE_STATE_DESCRIPTION,ONLY : P_QV,P_QC,P_QR,P_QI,P_QS,P_QG,P_QNI,P_QNR
      USE MODULE_MODEL_CONSTANTS
      USE MODULE_RA_GFDLETA,ONLY : CAL_MON_DAY,ZENITH
      USE MODULE_RADIATION_DRIVER
!      USE MODULE_SF_MYJSFC
!      USE MODULE_SURFACE_DRIVER
!      USE MODULE_PBL_DRIVER
!      USE MODULE_GWD
!      USE MODULE_CU_BMJ
!      USE MODULE_CUMULUS_DRIVER
!      USE MODULE_MP_ETANEW
!      USE MODULE_MICROPHYSICS_DRIVER
!      USE MODULE_MICROPHYSICS_ZERO_OUT
!-----------------------------------------------------------------------
!
      public : HAFS_radiation_init, HAFS_radiation_run, HAFS_radiation_finalize

      CONTAINS
     
#if 0
!! \section arg_table_HAFS_radiation_init Argument Table
!! \htmlinclude HAFS_radiation_init.html
!!
#endif
      SUBROUTINE HAFS_radiation_init(ncol, nlev, RTHRATEN,RTHRATENLW,   &
                      RTHRATENSW,CLDFRA,                                &
                      levsiz,XLAT,XLONG,n_ozmixm,                       &
                      alevsiz,no_src_types,                             &
                      cldfra_dp,cldfra_sh,                              & ! Optional for subgrid cloud by ckay
                      cldfra_old,                                       & ! Optional
                      ozmixm,pin,                                       & ! Optional
                      aerodm,pina,                                      & ! Optional
                      ra_lw_physics,ra_sw_physics,                      &
                      rrtmg_lwscheme,rrtmg_swscheme,                    &
                      mpiroot,mpirank,                                  &
                      errmsg, errflg )
!---------------------------------------------------------------------
   USE module_ra_rrtmg_lw  , ONLY : rrtmg_lwinit
   USE module_ra_rrtmg_sw  , ONLY : rrtmg_swinit
   USE module_radiation_gases, ONLY : AEROSOL_IN, oznini
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------
!MZ: change to local
!   LOGICAL,  INTENT(IN)           :: allowed_to_read
    logical, parameter            :: allowed_to_read =.true.

!MZ: change to local
   INTEGER                        :: ids, ide, jds, jde, kds, kde,      &
                                     ims, ime, jms, jme, kms, kme,      &
                                     its, ite, jts, jte, kts, kte


   INTEGER,                   INTENT(IN   )    :: levsiz, n_ozmixm
!   INTEGER,                   INTENT(IN   )    :: paerlev, n_aerosolc
   INTEGER,                   INTENT(IN   )    :: alevsiz, no_src_types
   INTEGER,                   INTENT(IN   )    :: mpirank
   INTEGER,                   INTENT(IN   )    :: mpiroot
   INTEGER,                   INTENT(IN   )    :: ra_lw_physics,        &
                                                  ra_sw_physics,        &
                                                  rrtmg_lwscheme,       &
                                                  rrtmg_swscheme 
   REAL(kind_phys),  DIMENSION(1:ncol) , INTENT(IN) ::  XLAT, XLONG

   REAL(kind_phys),  DIMENSION(1:ncol, levsiz, n_ozmixm ), OPTIONAL,    &
          INTENT(INOUT) ::                                  OZMIXM
   REAL(kind_phys),  DIMENSION(1:ncol, alevsiz, n_ozmixm-1, no_src_types ),
                                                         OPTIONAL,      &
          INTENT(INOUT) ::                                  aerodm


   REAL(kind_phys), DIMENSION(levsiz), OPTIONAL, INTENT(INOUT)  ::  PIN
   REAL(kind_phys), DIMENSION(alevsiz), OPTIONAL, INTENT(INOUT) ::  PINA

   REAL(kind_phys), DIMENSION(1:ncol, 1:nlev) , INTENT(INOUT) ::     &
                                                           RTHRATEN, &
                                                         RTHRATENLW, &
                                                         RTHRATENSW, &
                                                             CLDFRA

   REAL(kind_phys), DIMENSION(1:ncol, 1:nlev) , INTENT(INOUT) ,      & ! for subgrid cloud
                                                         OPTIONAL :: &
                                                          cldfra_dp, &
                                                          cldfra_sh, &
                                                         CLDFRA_OLD

   !REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(INOUT) :: EMISS
   character(len=*),          intent(  out) :: errmsg
   integer,                   intent(  out) :: errflg

   LOGICAL :: acswalloc = .false.
   LOGICAL :: aclwalloc = .false.
   INTEGER :: i, j, k, itf, jtf, ktf


         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (aclwalloc .and. acswalloc) return

         if (mpirank==mpiroot) then
            write(0,*) '--------------------------------------------------------'
            write(0,*) ' ---                   WARNING                       ---'
            write(0,*) ' ---  the CCPP HAFS RRTMG LW/SW scheme is currently  ---'
            write(0,*) ' ---  under development, use at your own risk        ---'
            write(0,*) ' ---                   WARNING                       ---'
            write(0,*) '--------------------------------------------------------'
         end if

         if (ra_lw_physics /= rrtmg_lwscheme .or.                       & 
                             ra_sw_physics /= rrtmg_swscheme ) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of     &
                RADLW/RADSW is different from HAFS RRTMG LWRAD/SWRAD"
            errflg = 1
            return
         end if


         ! Set internal dimensions
         ids = 1
         ims = 1
         its = 1
         ide = ncol
         ime = ncol
         ite = ncol
         jds = 1
         jms = 1
         jts = 1
         jde = 1
         jme = 1
         jte = 1
         kds = 1
         kms = 1
         kts = 1
         kde = Model%levr
         kme = Model%levr
         kte = Model%levr

!---------------------------------------------------------------------

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

!---------------------------------------------------------------------

!-- calculate radiation time step

!    STEPRA = nint(RADT*60./DT)
!    STEPRA = max(STEPRA,1)

!-- initialization


!   IF(start_of_simulation)THEN
     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RTHRATEN(i,k,j)=0.
        RTHRATENLW(i,k,j)=0.
        RTHRATENSW(i,k,j)=0.
        CLDFRA(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

     IF( PRESENT(cldfra_dp) ) THEN
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           cldfra_dp(i,k,j)=0.
           cldfra_sh(i,k,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

     if( present(cldfra_old) ) then
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           cldfra_old(i,k,j) = 0.
        ENDDO
        ENDDO
        ENDDO
     end if
!   ENDIF

!-- ww: attempt to use CAM ozone and some aerosol profiles in all rad schemes
!   note that CAM option will still do the same. 
!   n_ozmixm: no of months; levsiz: = 59, vertical dim
!   Read in CAM ozone data, and interpolate data to model grid
!   Interpolation of aerosols is done on domain 1 only

!#if (EM_CORE==1) 
!   IF ( config_flags%o3input .EQ. 2 .AND. id .EQ. 1 ) THEN
!#else
!   IF ( config_flags%o3input .EQ. 2 ) THEN
!#endif

   IF ( o3input .EQ. 2 ) THEN
      CALL oznini(ozmixm,pin,levsiz,n_ozmixm,XLAT,                &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte)
      if (errflg /= 0) return
   ENDIF

!   IF ( config_flags%aer_opt .EQ. 1 .AND. id .EQ. 1 ) THEN
      IF ( aer_opt .EQ. 1 ) THEN

             CALL aerosol_in(aerodm,pina,alevsiz,n_ozmixm-1,                   &
                     no_src_types,XLAT,XLONG,                           &
                     ids, ide, jds, jde, kds, kde,                      &
                     ims, ime, jms, jme, kms, kme,                      &
                     its, ite, jts, jte, kts, kte)
             if (errflg /= 0) return
      ENDIF


      if (ra_lw_physics == rrtmg_lwscheme) then

             CALL rrtmg_lwinit(                             &
!MZ                           p_top, allowed_to_read ,         &
                           ltp, allowed_to_read,            &
                           ids, ide, jds, jde, kds, kde,    &
                           ims, ime, jms, jme, kms, kme,    &
                           its, ite, jts, jte, kts, kte     )
             if (errflg /= 0) return

            aclwalloc = .true.
       endif

       if (ra_sw_physics == rrtmg_swscheme) then
  
             CALL rrtmg_swinit(                             &
                           allowed_to_read ,                &
                           ids, ide, jds, jde, kds, kde,    &
                           ims, ime, jms, jme, kms, kme,    &
                           its, ite, jts, jte, kts, kte     )
             if (errflg /= 0) return

             acswalloc = .true.

      END SUBROUTINE HAFS_radiation_init

      SUBROUTINE HAFS_radiation_finalize()
      END SUBROUTINE HAFS_radiation_finalize
!

#if 0
!> \section arg_table_HAFS_radiation_run Argument Table
!! \htmlinclude HAFS_radiation_run.html
!!
#endif
      SUBROUTINE HAFS_radiation_run (NTSD,DT,JULDAY,JULYR,XTIME,JULIAN  &
     &                    ,IHRST,NPHS,GLAT,GLON                         &
     &                    ,NRADS,NRADL                                  &
     &                    ,DETA1,DETA2,AETA1,AETA2,ETA1,ETA2,PDTOP,PT   &
     &                    ,PD,RES,PINT,T,Q,MOIST,THS,ALBEDO,EPSR        &
     &                    ,F_ICE,F_RAIN                                 &
     &                    ,GD_CLOUD,GD_CLOUD2                           &
     &                    ,SM,HBM2,CLDFRA,N_MOIST,RESTRT                &
     &                    ,RLWTT,RSWTT,RLWIN,RSWIN,RSWINC,RSWOUT        &
     &                    ,RLWTOA,RSWTOA,CZMEAN                         &
     &                    ,CFRACL,CFRACM,CFRACH,SIGT4                   &
     &                    ,ACFRST,NCFRST,ACFRCV,NCFRCV                  &
     &                    ,CUPPT,VEGFRC,SNOW,HTOP,HBOT                  &
     &                    ,Z,SICE,NUM_AEROSOLC,NUM_OZMIXM,OZMIXM,PIN    &
     &                    ,LEVSIZ,GRID,CONFIG_FLAGS                     &
     &                    ,RTHRATEN                                     &
     &                    ,re_cloud,re_ice,re_snow                      & ! G. Thompson
     &                    ,has_reqc,has_reqi,has_reqs                   & ! G. Thompson
     &                    ,hrswpd,hrlwpd                                & ! J. Henderson AER
     &                    ,SWUPT,SWUPTC,SWDNT,SWDNTC                    &
     &                    ,SWUPB,SWUPBC,SWDNB,SWDNBC                    &
     &                    ,LWUPT,LWUPTC,LWDNT,LWDNTC                    &
     &                    ,LWUPB,LWUPBC,LWDNB,LWDNBC                    &
     &                    ,ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC            &
     &                    ,ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC            &
     &                    ,ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC            &
     &                    ,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC            &
     &                    ,SWVISDIR ,SWVISDIF                           &  !ssib
     &                    ,SWNIRDIR, SWNIRDIF     )!&                            !ssib
!MZ*     &                    ,IDS,IDE,JDS,JDE,KDS,KDE                      &
!MZ*     &                    ,IMS,IME,JMS,JME,KMS,KME                      &
!MZ*     &                    ,ITS,ITE,JTS,JTE,KTS,KTE)
!***  NOTE ***
! RLWIN  - downward longwave at the surface (=TOTLWDN, now a local array)
! RSWIN  - downward shortwave at the surface (=TOTSWDN, now a local array)
! RSWINC - CLEAR-SKY downward shortwave at the surface (=TOTSWDNC, new for AQ)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    RADIATION   RADIATION OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 2002-06-04
!
! ABSTRACT:
!     RADIATION SERVES AS THE INTERFACE BETWEEN THE NCEP NONHYDROSTATIC
!     MESOSCALE MODEL AND THE WRF RADIATION DRIVER.
!
! PROGRAM HISTORY LOG:
!   02-06-04  BLACK      - ORIGINATOR
!   02-09-09  WOLFE      - CONVERTING TO GLOBAL INDEXING
!   04-11-18  BLACK      - THREADED
!   05-12-15  BLACK      - CONVERTED FROM IKJ TO IJK
!
! USAGE: CALL RADIATION FROM SOLVE_NMM
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM
!$$$
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER            :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE 

      INTEGER,INTENT(IN) :: IHRST,JULDAY,JULYR                          &
     &                     ,N_MOIST,NPHS,NRADL,NRADS,NTSD               &
     &                     ,NUM_AEROSOLC,NUM_OZMIXM,LEVSIZ

      INTEGER, INTENT(IN) :: ncol,nlay
     !MZ* REAL, INTENT(IN) :: OZMIXM(ims:ime,LEVSIZ,jms:jme,num_ozmixm), PIN(LEVSIZ)
      REAL(kind_phys), INTENT(IN) :: OZMIXM(1:ncol,LEVSIZ,num_ozmixm), PIN(LEVSIZ)
!
      !MZ* INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCFRCV,NCFRST
      INTEGER,DIMENSION(1:ncol),INTENT(INOUT) :: NCFRCV,NCFRST
!
      REAL(kind_phys),INTENT(IN) :: DT,PDTOP,PT,XTIME,JULIAN
!
!MZ*   REAL,DIMENSION(KMS:KME-1),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2
      REAL(kind_phys),DIMENSION(1:nlay-1),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2
!
!MZ*      REAL,DIMENSION(KMS:KME),INTENT(IN) :: ETA1,ETA2
      REAL(kind_phys),DIMENSION(1:nlay),INTENT(IN) :: ETA1,ETA2
!
      REAL(kind_phys),DIMENSION(1:ncol),INTENT(IN) :: ALBEDO            &
     &                                             ,EPSR,GLAT,GLON      &
     &                                             ,HBM2                &
     &                                             ,PD,RES,SICE,SM      &
     &                                             ,SNOW,THS,VEGFRC
!
      REAL(kind_phys),DIMENSION(1:ncol),INTENT(INOUT) :: CUPPT

!
!MZ*      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: Q,T,Z
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(IN) :: Q,T,Z
!
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(IN) :: F_ICE       &   !<--- Used only with physics (IKJ)
     &                                                     ,F_RAIN
!
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(INOUT) :: RTHRATEN     !<--- Used only with physics (IKJ)
!
!MZ* MOIST need to be seperated
!MZ*      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,N_MOIST)                   &
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay,,N_MOIST)
                                                 ,INTENT(INOUT) :: MOIST
!
      REAL(kind_phys),DIMENSION(1:ncol),INTENT(INOUT) :: ACFRCV,ACFRST  &
     &                                                ,HBOT,HTOP        &
     &                                                ,RLWIN,RLWTOA     &
     &                                                ,RSWIN,RSWOUT     &
     &                                                ,RSWINC,RSWTOA
!
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(INOUT) :: PINT     &
     &                                                        ,RLWTT    &
     &                                                        ,RSWTT
!
      REAL(kind_phys),DIMENSION(1:ncol),INTENT(INOUT) :: CFRACH,CFRACL    &
     &                                                ,CFRACM,CZMEAN    &
     &                                                ,SIGT4

!MZ*
!MZ*      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME ),INTENT(IN) ::            &   !<--- Used only with physics (IKJ)
      REAL(kind_phys),DIMENSION(1:ncol, 1:nlay),INTENT(IN) ::            &   !<--- Used only with physics (IKJ)
     &                              GD_CLOUD,GD_CLOUD2                  
!
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(INOUT) :: CLDFRA

      REAL(kind_phys), DIMENSION(1:ncol), INTENT(INOUT) ::          &
                      ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC,          &
                      ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,          &
                      ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC,          &
                      ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC

! TOA and surface, upward and downward, total and clear fluxes
      REAL(kind_phys), DIMENSION(1:ncol), INTENT(INOUT) ::          &
                        SWUPT,  SWUPTC,  SWDNT,  SWDNTC,          &
                        SWUPB,  SWUPBC,  SWDNB,  SWDNBC,          &
                        LWUPT,  LWUPTC,  LWDNT,  LWDNTC,          &
                        LWUPB,  LWUPBC,  LWDNB,  LWDNBC

      REAL(kind_phys), DIMENSION(1:ncol),                         &
         INTENT(OUT  )  ::                              SWVISDIR, &
                                                        SWVISDIF, &
                                                        SWNIRDIR, &
                                                        SWNIRDIF

!
!..Additions for coupling cloud physics effective radii and radiation.  G. Thompson
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(IN):: re_cloud,   &
                                                            re_ice,     &
                                                            re_snow
      INTEGER, INTENT(INOUT):: has_reqc, has_reqi, has_reqs

!..Output daily longwave and shortwave heating rates: ! J. Henderson AER
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(OUT):: hrlwpd,hrswpd
!
      LOGICAL,INTENT(IN) :: RESTRT
!MZ*
      TYPE(DOMAIN),TARGET :: GRID
!
!MZ* Grid_contrl
      TYPE(GRID_CONFIG_REC_TYPE),INTENT(IN) :: CONFIG_FLAGS
!
!-----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
!-----------------------------------------------------------------------
      INTEGER :: I,IENDX,II,ISTAT,J,JDAY,JMONTH,K,KMNTH,N,NRAD
!
      INTEGER,DIMENSION(3) :: IDAT
      INTEGER,DIMENSION(12) :: MONTH=(/31,28,31,30,31,30,31,31          &
     &                                ,30,31,30,31/)
!
      REAL(kind_phys) :: CAPA,DAYI,DPL,FICE,FRAIN,GMT,HOUR,PLYR,PSFC    &
     &       ,QI,QR,QW,RADT,TIMES,WC,TDUM
!
!MZ*      REAL,DIMENSION(KMS:KME-1) :: QL,TL
      REAL(kind_phys),DIMENSION(1:nlay-1) :: QL,TL
!
      REAL(kind_phys),DIMENSION(1:ncol) :: CUPPTR,CZEN,HBOTR,HTOPR      &
     &                                  ,SWCF, LWCF                     &
     &                                  ,PDSL,REXNSFC,SWNETDN           &
     &                                  ,TOT,TOTLWDN,TOTSWDN,TOTSWDNC   &
     &                                  ,TSFC,XLAND,XLAT,XLON, HT
!
!
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay) :: CLFR,DZ               &   !<--- Used only with physics (IKJ)
     &                                          ,P8W,P_PHY,PI_PHY       &
     &                                          ,RR,T8W                 &
     &                                          ,THRATENLW,THRATENSW    &
     &                                          ,TH_PHY,T_PHY,Z_PHY
!
      REAL(kind_phys),DIMENSION(:,:,:,:),ALLOCATABLE :: MOIST_TRANS
!
      LOGICAL :: WARM_RAIN
      LOGICAL :: IS_CAMMGMP_USED=.FALSE.

      REAL(kind_phys) :: DXKM, DYKM

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!*****
!***** NOTE: THIS IS HARDWIRED FOR CALLS TO LONGWAVE AND SHORTWAVE
!*****       AT EQUAL INTERVALS
!*****
!-----------------------------------------------------------------------
!
      NRAD=NRADS
      RADT=DT*NRADS/60.
!
!-----------------------------------------------------------------------
!
      ALLOCATE(MOIST_TRANS(IMS:IME,KMS:KME,JMS:JME,N_MOIST),STAT=ISTAT)
!
!-----------------------------------------------------------------------
!
      CAPA=R_D/CP

      DXKM=grid%dlmd*0.01745329*6371200. ! numbers from module_initialize_real.F
      DYKM=grid%dy_nmm
!
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=MYJS2,MYJE2
      DO I=MYIS1,MYIE1
!
        PDSL(I,J)=PD(I,J)*RES(I,J)
        P8W(I,KTE+1,J)=PT
        XLAND(I,J)=SM(I,J)+1.
        PSFC=PD(I,J)+PDTOP+PT
        REXNSFC(I,J)=(PSFC*1.E-5)**CAPA
        TSFC(I,J)=THS(I,J)*REXNSFC(I,J)
        T8W(I,KTS,J)=TSFC(I,J)
        P8W(I,KTS,J)=ETA1(KTS)*PDTOP+ETA2(KTS)*PDSL(I,J)+PT
        Z_PHY(I,KTS,J)=Z(I,J,KTS)
        HT(I,J)=Z(I,J,KTS)
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(dpl,i,j,k,plyr,ql,qr,tl)
      DO J=MYJS2,MYJE2
      DO I=MYIS1,MYIE1
        DO K=KTS,KTE
          DPL=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
          QL(K)=MAX(Q(I,J,K),EPSQ)
          PLYR=AETA1(K)*PDTOP+AETA2(K)*PDSL(I,J)+PT
          TL(K)=T(I,J,K)
!
          RR(I,K,J)=PLYR/(R_D*TL(K)*(1.+P608*QL(K)))
          T_PHY(I,K,J)=TL(K)
          TH_PHY(I,K,J)=TL(K)*(1.E5/PLYR)**CAPA
          P8W(I,K+1,J)=ETA1(K+1)*PDTOP+ETA2(K+1)*PDSL(I,J)+PT
          P_PHY(I,K,J)=PLYR
          PI_PHY(I,K,J)=(PLYR*1.E-5)**CAPA
          DZ(I,K,J)=TL(K)*(P608*QL(K)+1.)*R_D                           &
     &                 *(P8W(I,K,J)-P8W(I,K+1,J))                       &
     &                 /(P_PHY(I,K,J)*G)
!
          RTHRATEN(I,K,J)=0.
          THRATENLW(I,K,J)=0.
          THRATENSW(I,K,J)=0.
!         PM2_5_DRY(I,K,J)=0.
!         PM2_5_WATER(I,K,J)=0.

        ENDDO
!
        DO K=KTS+1,KTE
          T8W(I,K,J)=0.5*(TL(K-1)+TL(K))
        ENDDO
!        T8W(I,KTE+1,J)=-1.E20 
! For RRTM 
        T8W(I,KTE+1,J)=T8W(I,KTE,J) + 0.5*(T8W(I,KTE,J)-T8W(I,KTE-1,J))
!
      ENDDO
      ENDDO
!
      GMT=REAL(IHRST)
!
!$omp parallel do                                                       &
!$omp& private(i,j,k)
      DO K=KMS,KME
        DO J=JMS,JME
        DO I=IMS,IME
          CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
      ENDDO
!
!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=JMS,JME
        DO I=IMS,IME
          CFRACH(I,J)=0.
          CFRACL(I,J)=0.
          CFRACM(I,J)=0.
          CZMEAN(I,J)=0.
          SIGT4(I,J)=0.
          TOTSWDN(I,J)=0.   ! TOTAL (clear+cloudy sky) shortwave down at the surface
          TOTSWDNC(I,J)=0.  ! CLEAR SKY shortwave down at the surface
          SWNETDN(I,J)=0.   ! Net (down - up) total (clear+cloudy sky) shortwave at the surface
          TOTLWDN(I,J)=0.   ! Total longwave down at the surface
          CUPPTR(I,J)=CUPPT(I,J)   ! Temporary array set to zero in radiation
          HTOPR(I,J) =0.
          HBOTR(I,J) = REAL(KTE+1)
!          HBOTR(I,J) =0.
          SWCF(I,J) =0.
          LWCF(I,J) =0.
!
!-- NOTE:  HBOTR, HTOPR are passed into radiation and set equal to HBOT, HTOP.  HBOT, HTOP are
!          reset to clear sky values to be used by the ARW.  At the bottom of this subroutine,
!          HBOT, HTOP are re-defined again to values stored in HBOTR, HTOPR.  HBOT, HTOP are
!          reset to clear sky values after the call to radiation and after the top of the hour
!          in subroutine CUCNVC below.
!
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  TRANSPOSE THE MOIST ARRAY (IJK) FOR THE PHYSICS (IKJ).
!***  REMEMBER THAT MOIST AND MOIST_TRANS ARE ONLY USED WITH
!***  THE PHYSICS AND THUS THE P_QV SLOT (=2) IS MIXING RATIO,
!***  NOT SPECIFIC HUMIDITY.
!-----------------------------------------------------------------------
!
      DO N=1,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO K=KMS,KME
        DO J=JMS,JME
        DO I=IMS,IME
          MOIST_TRANS(I,K,J,N)=MOIST(I,J,K,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
!***  CALL THE INNER DRIVER.
!
!-----------------------------------------------------------------------
!
      CALL SET_TILES(GRID,IDS+1,IDE-1,JDS+2,JDE-2,ITS,ITE,JTS,JTE)

      CALL RADIATION_DRIVER(                                            &
     &                  IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=KDS,KDE=KDE &
     &                 ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=KMS,KME=KME &
     &                 ,I_START=GRID%I_START,I_END=GRID%I_END           &
     &                 ,J_START=GRID%J_START,J_END=GRID%J_END           &
     &                 ,KTS=KTS,KTE=KTE,NUM_TILES=GRID%NUM_TILES        &
     &                 ,ITIMESTEP=NTSD,DT=DT                            &
     &                 ,ICLOUD_CU=config_flags%ICLOUD_CU                &
     &                 ,QC_CU=GRID%QC_CU,QI_CU=GRID%QI_CU               &
     &                 ,DX=DXKM,DY=DYKM          &
     &                 ,DXKM=dxkm                                       &
! WRF-Solar variables
                       ,swint_opt=config_flags%swint_opt                &   
     &       ,SWDDIR=grid%swddir,SWDDNI=grid%swddni,SWDDIF=grid%swddif     & ! jararias
     &       ,SWDDIRC=grid%swddirc,SWDDNIC=grid%swddnic                    & ! clearsky sw diags
     &       ,Gx=grid%Gx,Bx=grid%Bx,gg=grid%gg,bb=grid%bb                  & ! for sza-interpolation
     &       ,swdown_ref=grid%swdown_ref,swddir_ref=grid%swddir_ref        & !
     &       ,coszen_ref=grid%coszen_ref                                   & 
     &       ,coszen=grid%coszen                                   & 
     &       ,hrang=grid%hrang                                   & 
     &       ,ht=ht                                   & 
     &       ,aer_type=config_flags%aer_type                                                        &
     &       ,aer_aod550_opt=config_flags%aer_aod550_opt,aer_aod550_val=config_flags%aer_aod550_val &
     &       ,aer_angexp_opt=config_flags%aer_angexp_opt,aer_angexp_val=config_flags%aer_angexp_val &
     &       ,aer_ssa_opt=config_flags%aer_ssa_opt,aer_ssa_val=config_flags%aer_ssa_val             &
     &       ,aer_asy_opt=config_flags%aer_asy_opt,aer_asy_val=config_flags%aer_asy_val             &
     &                 ,RTHRATENLW=THRATENLW,RTHRATENSW=THRATENSW       &
     &                 ,RTHRATEN=RTHRATEN                               &
     &                 ,HRSWPD=HRSWPD, HRLWPD=HRLWPD                    & ! J. Henderson AER
     &                 ,CEN_LAT=grid%cen_lat                            &
     &                 ,GLW=TOTLWDN,GSW=SWNETDN,SWDOWN=TOTSWDN          &
     &                 ,XLAT=grid%HLAT,XLONG=grid%HLON,ALBEDO=ALBEDO,EMISS=EPSR   &
     &                 ,XICE=SICE,XLAND=XLAND,Z=Z,TSK=TSFC              &
     &                 ,N_AEROSOLC=NUM_AEROSOLC,PAERLEV=GRID%PAERLEV    &
     &                 ,ID=grid%id                                      &
     &                 ,CAM_ABS_DIM1=GRID%CAM_ABS_DIM1                  &
     &                 ,CAM_ABS_DIM2=GRID%CAM_ABS_DIM2                  &
     &                 ,CAM_ABS_FREQ_S=GRID%CAM_ABS_FREQ_S              &
     &                 ,ALEVSIZ=grid%alevsiz,no_src_types=grid%no_src_types   &
     &                 ,LEVSIZ=LEVSIZ,N_OZMIXM=NUM_OZMIXM               &
     &                 ,OZMIXM=OZMIXM,PIN=PIN                           &
     &                 ,HTOP=HTOP,HBOT=HBOT,CUPPT=CUPPTR                &
     &                 ,HTOPR=HTOPR,HBOTR=HBOTR                         &
     &                 ,VEGFRA=VEGFRC,SNOW=SNOW                         &
     &                 ,RHO=RR,P8W=P8W,P=P_PHY,PI=PI_PHY                &
     &                 ,DZ8W=DZ,T=T_PHY,T8W=T8W,GMT=GMT                 &
     &                 ,JULDAY=JULDAY,JULYR=JULYR,NPHS=NPHS             &
     &                 ,JULIAN=JULIAN,XTIME=XTIME                       &
     &                 ,YR=JULYR                                        &
     &                 ,LW_PHYSICS=CONFIG_FLAGS%RA_LW_PHYSICS           &
     &                 ,SW_PHYSICS=CONFIG_FLAGS%RA_SW_PHYSICS           &
     &                 ,RADT=RADT,RA_CALL_OFFSET=GRID%RA_CALL_OFFSET    &
     &                 ,STEPRA=NRAD,ICLOUD=config_flags%ICLOUD          &
     &                 ,cldovrlp=config_flags%cldovrlp                  & ! J. Henderson AER: cldovrlp namelist value
     &                 ,WARM_RAIN=WARM_RAIN                             &
     &                 ,SWDOWNC=TOTSWDNC,CLDFRA=CLFR                    &
     &                 ,SWUPT=SWUPT                                     &
     &                 ,SWUPTC=SWUPTC                                   &
     &                 ,SWDNT=SWDNT                                     &
     &                 ,SWDNTC=SWDNTC                                   &
     &                 ,SWUPB=SWUPB                                     &
     &                 ,SWUPBC=SWUPBC                                   &
     &                 ,SWDNB=SWDNB                                     &
     &                 ,SWDNBC=SWDNBC                                   &
     &                 ,LWUPT=LWUPT                                     &
     &                 ,LWUPTC=LWUPTC                                   &
     &                 ,LWDNT=LWDNT                                     &
     &                 ,LWDNTC=LWDNTC                                   &
     &                 ,LWUPB=LWUPB                                     &
     &                 ,LWUPBC=LWUPBC                                   &
     &                 ,LWDNB=LWDNB                                     &
     &                 ,LWDNBC=LWDNBC                                   &
     &                 ,ACSWUPT=ACSWUPT                                 &
     &                 ,ACSWUPTC=ACSWUPTC                               &
     &                 ,ACSWDNT=ACSWDNT                                 &
     &                 ,ACSWDNTC=ACSWDNTC                               &
     &                 ,ACSWUPB=ACSWUPB                                 &
     &                 ,ACSWUPBC=ACSWUPBC                               &
     &                 ,ACSWDNB=ACSWDNB                                 &
     &                 ,ACSWDNBC=ACSWDNBC                               &
     &                 ,ACLWUPT=ACLWUPT                                 &
     &                 ,ACLWUPTC=ACLWUPTC                               &
     &                 ,ACLWDNT=ACLWDNT                                 &
     &                 ,ACLWDNTC=ACLWDNTC                               &
     &                 ,ACLWUPB=ACLWUPB                                 &
     &                 ,ACLWUPBC=ACLWUPBC                               &
     &                 ,ACLWDNB=ACLWDNB                                 &
     &                 ,ACLWDNBC=ACLWDNBC                               &
     &        ,SWVISDIR=swvisdir ,SWVISDIF=swvisdif                     &  !ssib
     &        ,SWNIRDIR=swnirdir ,SWNIRDIF=swnirdif                     &  !ssib
     &                 ,re_cloud=grid%re_cloud                          & ! G. Thompson
     &                 ,re_ice=grid%re_ice                              & ! G. Thompson
     &                 ,re_snow=grid%re_snow                            & ! G. Thompson
     &                 ,has_reqc=has_reqc                               & ! G. Thompson
     &                 ,has_reqi=has_reqi                               & ! G. Thompson
     &                 ,has_reqs=has_reqs                               & ! G. Thompson
     &                 ,RSWTOA=RSWTOA,RLWTOA=RLWTOA                     &
     &                 ,CZMEAN=CZMEAN,CFRACL=CFRACL                     &
     &                 ,CFRACM=CFRACM,CFRACH=CFRACH                     &
     &                 ,ACFRST=ACFRST,NCFRST=NCFRST                     &
     &                 ,ACFRCV=ACFRCV,NCFRCV=NCFRCV                     &
     &                 ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN               &
     &                 ,LWCF=LWCF,SWCF=SWCF                             &
     & ,O3INPUT=config_flags%O3INPUT,AER_OPT=config_flags%AER_OPT,O3RAD=grid%o3rad  &
     &                 ,SF_SURFACE_PHYSICS=CONFIG_FLAGS%SF_SURFACE_PHYSICS &
     &                 ,QV=MOIST_TRANS(IMS,KMS,JMS,P_QV),F_QV=F_QV      &
     &                 ,QC=MOIST_TRANS(IMS,KMS,JMS,P_QC),F_QC=F_QC      &
     &                 ,QR=MOIST_TRANS(IMS,KMS,JMS,P_QR),F_QR=F_QR      &
     &                 ,QI=MOIST_TRANS(IMS,KMS,JMS,P_QI),F_QI=F_QI      &
     &                 ,QS=MOIST_TRANS(IMS,KMS,JMS,P_QS),F_QS=F_QS      &
     &                 ,QG=MOIST_TRANS(IMS,KMS,JMS,P_QG),F_QG=F_QG      &
     &                 ,IS_CAMMGMP_USED=IS_CAMMGMP_USED                 &
     &                 ,EXPLICIT_CONVECTION=config_flags%cu_physics==0  &
     &                 ,AERCU_OPT=config_flags%aercu_opt                &
     &                 ,CU_PHYSICS=config_flags%cu_physics              &
     &                 ,MP_PHYSICS=CONFIG_FLAGS%MP_PHYSICS              &
     &                 ,CALC_CLEAN_ATM_DIAG=CONFIG_FLAGS%CALC_CLEAN_ATM_DIAG)


!
!-----------------------------------------------------------------------
!
!***  UPDATE FLUXES AND TEMPERATURE TENDENCIES.
!
!-----------------------------------------------------------------------
!***  SHORTWAVE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      nrads_block: IF(MOD(NTSD,NRADS)==0)THEN
!-----------------------------------------------------------------------
!
        IF(CONFIG_FLAGS%RA_SW_PHYSICS/=GFDLSWSCHEME)THEN
!
!-----------------------------------------------------------------------
!***  COMPUTE CZMEAN FOR NON-GFDL SHORTWAVE
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,j)
          DO J=MYJS,MYJE
          DO I=MYIS,MYIE
            CZMEAN(I,J)=0.
            TOT(I,J)=0.
          ENDDO
          ENDDO
!
          CALL CAL_MON_DAY(JULDAY,JULYR,JMONTH,JDAY)
          IDAT(1)=JMONTH
          IDAT(2)=JDAY
          IDAT(3)=JULYR
!
          DO II=0,NRADS,NPHS
            TIMES=NTSD*DT+II*DT
            CALL ZENITH(TIMES,DAYI,HOUR,IDAT,IHRST,GLON,GLAT,CZEN       &
     &                 ,MYIS,MYIE,MYJS,MYJE                             &
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                         &
     &                 ,IMS,IME,JMS,JME,KMS,KME                         &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)
!
!$omp parallel do                                                       &
!$omp& private(i,j)
            DO J=MYJS,MYJE
            DO I=MYIS,MYIE
              IF(CZEN(I,J)>0.)THEN
                CZMEAN(I,J)=CZMEAN(I,J)+CZEN(I,J)
                TOT(I,J)=TOT(I,J)+1.
              ENDIF
            ENDDO
            ENDDO
!
          ENDDO
!
!$omp parallel do                                                       &
!$omp& private(i,j)
          DO J=MYJS,MYJE
          DO I=MYIS,MYIE
            IF(TOT(I,J)>0.)CZMEAN(I,J)=CZMEAN(I,J)/TOT(I,J)
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  COMPUTE TOTAL SFC SHORTWAVE DOWN FOR NON-GFDL SCHEMES
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,j)
          DO J=MYJS2,MYJE2
          DO I=MYIS1,MYIE1
!
            IF(HBM2(I,J)>0.5)THEN
              TOTSWDN(I,J)=SWNETDN(I,J)/(1.-ALBEDO(I,J))
!
!--- No value currently available for clear-sky solar fluxes from
!    non GFDL schemes, though it's needed for air quality forecasts.
!    For the time being, set to the total downward solar fluxes.
!
              TOTSWDNC(I,J)=TOTSWDN(I,J)
            ENDIF
!
          ENDDO
          ENDDO
!
        ENDIF   !End non-GFDL block
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,iendx,j)
        DO J=MYJS2,MYJE2
          IENDX=MYIE1
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=MYIS1,IENDX
!
            RSWIN(I,J)=TOTSWDN(I,J)
            RSWINC(I,J)=TOTSWDNC(I,J)
            RSWOUT(I,J)=TOTSWDN(I,J)-SWNETDN(I,J)
!
          ENDDO
        ENDDO
!
!$omp parallel do                                                       &
!$omp& private(i,iendx,j,k)
        DO J=MYJS2,MYJE2
          IENDX=MYIE1
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=MYIS1,IENDX
            DO K=KTS,KTE
              RSWTT(I,J,K)=THRATENSW(I,K,J)*PI_PHY(I,K,J)
            ENDDO
!
          ENDDO
        ENDDO
!
      ENDIF nrads_block
!
!-----------------------------------------------------------------------
!***  LONGWAVE
!-----------------------------------------------------------------------
!
      nradl_block: IF(MOD(NTSD,NRADL)==0)THEN
!
!$omp parallel do                                                       &
!$omp& private(i,iendx,j)
        DO J=MYJS2,MYJE2
          IENDX=MYIE1
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=MYIS1,IENDX
!
            IF(HBM2(I,J)>0.5)THEN
              TDUM=T(I,J,KTS)
              SIGT4(I,J)=STBOLT*TDUM*TDUM*TDUM*TDUM
              RLWIN(I,J)=TOTLWDN(I,J)
            ENDIF
!
          ENDDO
        ENDDO
!
!$omp parallel do                                                       &
!$omp& private(i,iendx,j,k)
        DO J=MYJS2,MYJE2
          IENDX=MYIE1
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
!
          DO K=KTS,KTE
          DO I=MYIS1,IENDX
            IF(HBM2(I,J)>0.5)THEN
                RLWTT(I,J,K)=THRATENLW(I,K,J)*PI_PHY(I,K,J)
            ENDIF
          ENDDO
          ENDDO
!
        ENDDO
!
      ENDIF nradl_block
!
!-----------------------------------------------------------------------
!***  STORE 3D CLOUD FRACTIONS.
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,iendx,j,k)
      DO K=KTS,KTE
        DO J=MYJS2,MYJE2
          IENDX=MYIE1
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=MYIS1,IENDX
            CLDFRA(I,J,K)=CLFR(I,K,J)
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  RESET THE DIAGNOSTIC CONVECTIVE CLOUD TOPS/BOTTOMS AFTER
!***  EACH RADIATION CALL.
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,iendx,j)
!if (config_flags%ra_sw_physics/=hwrfswscheme.and.config_flags%ra_lw_physics/=hwrflwscheme)then
      DO J=MYJS2,MYJE2
        IENDX=MYIE1
        IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
        DO I=MYIS1,IENDX
          HBOT(I,J)=HBOTR(I,J)
          HTOP(I,J)=HTOPR(I,J)
          CUPPT(I,J)=CUPPTR(I,J)
        ENDDO
      ENDDO
!endif
!
!-----------------------------------------------------------------------
!***  ZERO OUT BOUNDARY ROWS.
!-----------------------------------------------------------------------
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        IF(HBM2(I,J)<0.5)THEN
          ACFRST(I,J)=0.
          ACFRCV(I,J)=0.
          CFRACL(I,J)=0.
          CFRACM(I,J)=0.
          CFRACH(I,J)=0.
          RSWTOA(I,J)=0.
          RLWTOA(I,J)=0.
        ENDIF
      ENDDO
      ENDDO
!
!
!-----------------------------------------------------------------------
!***  UPDATE THE PROGNOSTIC MOIST ARRAY.
!-----------------------------------------------------------------------
!
      DO N=2,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO J=JMS,JME
        DO K=KMS,KME
        DO I=IMS,IME
          MOIST(I,J,K,N)=MOIST_TRANS(I,K,J,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO
!
      DEALLOCATE(MOIST_TRANS,STAT=ISTAT)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE HAFS_radiation_run

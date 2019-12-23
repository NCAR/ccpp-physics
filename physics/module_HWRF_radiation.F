!>\file module_HWRF_radiation.F
!! This is the HWRF RRTMG driver wrapper module.

      MODULE HWRF_radiation
!
      use machine, only : kind_phys
!-----------------------------------------------------------------------
      USE MODULE_MODEL_CONSTANTS
      USE module_radiation_astronomy ,ONLY : CAL_MON_DAY,ZENITH
      USE HWRF_rrtmg_driver
!-----------------------------------------------------------------------
!
      public : HWRF_radiation_init, HWRF_radiation_run,                 &
               HWRF_radiation_finalize

      CONTAINS
     
!! \section arg_table_HWRF_radiation_init Argument Table
!! \htmlinclude HWRF_radiation_init.html
!!
      SUBROUTINE HWRF_radiation_init(ncol, nlev, KDT,                   &
                      RTHRATEN,RTHRATENLW,                              &
                      RTHRATENSW,CLDFRA,                                &
                      levsiz,XLAT,XLONG,n_ozmixm,                       &
                      alevsiz,no_src_types,                             &
                      has_reqc,has_reqi,has_reqs,                       &
                      ozmixm,pin,                                       &
                      aerodm,pina,                                      &
                      mpiroot,mpirank,                                  &
                      errmsg, errflg )
!---------------------------------------------------------------------
      USE module_ra_rrtmg_lw  , ONLY : rrtmg_lwinit
      USE module_ra_rrtmg_sw  , ONLY : rrtmg_swinit
      USE module_radiation_gases, ONLY : AEROSOL_IN, oznini
!---------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER,                   INTENT(IN   )    :: levsiz, n_ozmixm
      INTEGER,                   INTENT(IN   )    :: NCOL, NLAY, KDT
      INTEGER,                   INTENT(IN   )    :: alevsiz, no_src_types
      INTEGER,                   INTENT(INOUT)    :: has_reqc, has_reqi,&
                                                     has_reqs
      INTEGER,                   INTENT(IN   )    :: mpirank
      INTEGER,                   INTENT(IN   )    :: mpiroot
      REAL(KIND_PHYS),  DIMENSION(1:ncol) , INTENT(IN) ::  XLAT, XLONG
      REAL(KIND_PHYS),  DIMENSION(1:ncol, levsiz, n_ozmixm ),           &
                                                           OPTIONAL,    &
          INTENT(INOUT) ::                                  OZMIXM
      REAL(KIND_PHYS),  DIMENSION(1:ncol, alevsiz, n_ozmixm-1, no_src_types ),
                                                         OPTIONAL,      &
          INTENT(INOUT) ::                                  aerodm
      REAL(KIND_PHYS), DIMENSION(levsiz), OPTIONAL, INTENT(INOUT)  ::  PIN
      REAL(KIND_PHYS), DIMENSION(alevsiz), OPTIONAL, INTENT(INOUT) ::  PINA
      REAL(KIND_PHYS), DIMENSION(1:ncol, 1:nlev) , INTENT(INOUT) ::     &
                                                              RTHRATEN, &
                                                            RTHRATENLW, &
                                                            RTHRATENSW, &
                                                              CLDFRA
      CHARACTER(LEN=*),          INTENT(  OUT) :: errmsg
      INTEGER,                   INTENT(  OUT) :: errflg

      !local
      LOGICAL, PARAMETER             :: allowed_to_read =.true.
      INTEGER                        :: ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte

      LOGICAL :: acswalloc = .false.
      LOGICAL :: aclwalloc = .false.
      INTEGER :: i, j, k


      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (aclwalloc .and. acswalloc) return

      if (mpirank==mpiroot) then
        write(0,*) '---------------------------------------------------'
        write(0,*) '-                   WARNING                       -'
        write(0,*) '-  the CCPP HWRF RRTMG LW/SW scheme is currently  -'
        write(0,*) '-  under development, use at your own risk        -'
        write(0,*) '-                   WARNING                       -'
        write(0,*) '---------------------------------------------------'
      end if

      !Determin if we will compute and pass radiative effective radii of 
      !cloud water, ice, and snow. Currently not supported in
      !Ferrier-Aligo with RRTMG-LW/SW.
      has_reqc = 0
      has_reqi = 0
      has_reqs = 0


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
      kde = levr
      kme = levr
      kte = levr


!---------------------------------------------------------------------

!-- calculate radiation time step

!      STEPRA = nint(RADT*60./DT)
!      STEPRA = max(STEPRA,1)

!-- initialization
!MZ: rad_reset 
      IF(KDT .LE. 1) THEN
      DO j=jts,jte
      DO k=kts,kte
      DO i=its,ite
        RTHRATEN(i,k,j)=0.
        RTHRATENLW(i,k,j)=0.
        RTHRATENSW(i,k,j)=0.
        CLDFRA(i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO
      ENDIF

!-- use CAM ozone and some aerosol profiles in HWRF rad schemes
!   n_ozmixm: no of months; levsiz: = 59, vertical dim
!   Read in CAM ozone data, and interpolate data to model grid
!   Interpolation of aerosols is done on domain 1 only

      CALL oznini(ozmixm,pin,levsiz,n_ozmixm,XLAT,                      &
                     ids, ide, jds, jde, kds, kde,                      &
                     ims, ime, jms, jme, kms, kme,                      &
                     its, ite, jts, jte, kts, kte,                      &
                     mpirank, mpiroot, errflg, errmsg)
                   
      if (errflg /= 0) return

      CALL aerosol_in(aerodm,pina,alevsiz,n_ozmixm-1,                   &
                     no_src_types,XLAT,XLONG,                           &
                     ids, ide, jds, jde, kds, kde,                      &
                     ims, ime, jms, jme, kms, kme,                      &
                     its, ite, jts, jte, kts, kte)
      if (errflg /= 0) return

      CALL rrtmg_lwinit(                                                &
                     ltp, allowed_to_read,                              &
                     ids, ide, jds, jde, kds, kde,                      &
                     ims, ime, jms, jme, kms, kme,                      &
                     its, ite, jts, jte, kts, kte)
      if (errflg /= 0) return

      aclwalloc = .true.

      CALL rrtmg_swinit(                                                &
                     allowed_to_read ,                                  &
                     ids, ide, jds, jde, kds, kde,                      &
                     ims, ime, jms, jme, kms, kme,                      &
                     its, ite, jts, jte, kts, kte )
      if (errflg /= 0) return

      acswalloc = .true.

      END SUBROUTINE HWRF_radiation_init

      SUBROUTINE HWRF_radiation_finalize()
      END SUBROUTINE HWRF_radiation_finalize
!

!> \section arg_table_HWRF_radiation_run Argument Table
!! \htmlinclude HWRF_radiation_run.html
!!
      SUBROUTINE HWRF_radiation_run (NTSD,DT,JULDAY,JULYR,XTIME,JULIAN  &
     &                    ,IHRST,NPHS,GLAT,GLON                         &
     &                    ,NRADS,NRADL, DX ,p8w,prsl,tsfc,T,Q           &
     &                    ,QV,QC,QI,QR,QS,QG                            & !MOIST: dry mixing ratio
     &                    ,ALBEDO,EPSR                                  &
     &                    ,F_ICE,F_RAIN,SM,CLDFRA                       &
     &                    ,RLWTT,RSWTT,RLWIN,RSWIN,RSWINC,RSWOUT,GSW,   &
     &                    ,RLWTOA,CZMEAN                         &
     &                    ,SNOW,SICE, NUM_OZMIXM,OZMIXM,PIN,LEVSIZ      &
     &                    ,RTHRATEN,RTHRATENLW,RTHRATENSW               &
     &                    ,re_cloud,re_ice,re_snow                      &
     &                    ,has_reqc,has_reqi,has_reqs                   &
     &                    ,swddir, swddni,swddif, swddirc,swddnic       &
     &                    ,gx,bx,gg,bb,swdown_ref, swddir_ref           &
     &                    ,coszen_ref, coszen,hrang,xlat,xlong          &
     &                    ,paerlev, ALEVSIZ,no_src_types,               &   !NML
     &                    ,CAM_ABS_FREQ_S,swint_opt,icloud,cldovrlp     &
     &                    ,RA_CALL_OFFSET, o3input,aer_opt,o3rad        &
     &                    ,sf_surface_physics                           &
     &                    ,hrswpd,hrlwpd                                &
     &                    ,SWUPT,SWUPTC,SWDNT,SWDNTC                    &
     &                    ,SWUPB,SWUPBC,SWDNB,SWDNBC                    &
     &                    ,LWUPT,LWUPTC,LWDNT,LWDNTC                    &
     &                    ,LWUPB,LWUPBC,LWDNB,LWDNBC                    &
     &                    ,ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC            &
     &                    ,ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC            &
     &                    ,ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC            &
     &                    ,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC            &
     &                    ,SWVISDIR ,SWVISDIF                           &
     &                    ,SWNIRDIR, SWNIRDIF                           &
     &                    ,MPIROOT, MPIRANK                             &
     &                    ,ERRMSG, ERRFLG    )

!***********************************************************************
      IMPLICIT NONE

      INTEGER,         INTENT(IN) :: IHRST,JULDAY,JULYR                 &
     &                              ,NPHS,NRADL,NRADS,NTSD              &
     &                              ,NUM_OZMIXM,LEVSIZ

      INTEGER,         INTENT(IN) :: NCOL, NLAY
      REAL(KIND_PHYS), INTENT(IN) :: DX
      REAL(KIND_PHYS), INTENT(IN) :: OZMIXM(1:NCOL,LEVSIZ,NUM_OZMIXM)   &
                                    ,PIN(LEVSIZ)
!MZ: predefine in GFS_typedef
      integer,         intent(in) :: swint_opt
      INTEGER,         INTENT(IN) :: o3input, aer_opt
      INTEGER,         INTENT(IN) :: ICLOUD,ra_call_offset
      INTEGER,         INTENT(IN) :: cldovrlp        
      INTEGER,         INTENT(IN) :: alevsiz, no_src_types
      INTEGER,         INTENT(IN) :: levsiz, n_ozmixm
      INTEGER,         INTENT(IN) :: paerlev, n_aerosolc,
      REAL(KIND_PHYS), INTENT(IN) :: cam_abs_freq_s
      LOGICAL,         INTENT(IN) :: warm_rain
      INTEGER,    OPTIONAL, INTENT(IN   )    ::  sf_surface_physics
!
!      INTEGER,         DIMENSION(1:NCOL),   INTENT(INOUT) :: NCFRCV,    &
!                                                             NCFRST

      REAL(KIND_PHYS), INTENT(IN) :: DT,XTIME,JULIAN

      REAL(KIND_PHYS), DIMENSION(1:ncol),      INTENT(IN)    :: ALBEDO  &
     &                                             ,EPSR,GLAT,GLON      &
!MZ     &                                             ,HBM2                &
     &                                             ,SICE,SM      &
     &                                             ,SNOW
!MZ ,VEGFRC
!MZ      REAL(KIND_PHYS),DIMENSION(1:ncol),       INTENT(INOUT) :: CUPPT
      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlay),INTENT(IN)    :: Q,T,Z,  &
                                                                prsl
      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlay+1),INTENT(IN)  :: p8w 
      REAL(KIND_PHYS),DIMENSION(1:ncol),       INTENT(IN)  :: tsfc

      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlay),INTENT(IN)    :: F_ICE   &
     &                                                         ,F_RAIN
      REAL(KIND_PHYS),DIMENSION(1:NCOL),       INTENT(INOUT) :: SWDDIR  &
                                                               ,SWDDNI  &
                                                               ,SWDDIF  &
                                                               ,SWDDNIC &
                                                               ,SWDDIRC 
      REAL(KIND_PHYS), OPTIONAL, DIMENSION(1:NCOL), INTENT(OUT) :: COSZEN
      REAL(KIND_PHYS), OPTIONAL, DIMENSION(1:NCOL), INTENT(OUT) :: HRANG
      REAL(KIND_PHYS),DIMENSION(1:ncol) , INTENT(IN) ::  XLAT, XLONG
      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlay),INTENT(INOUT) ::RTHRATEN &
                                                            ,RTHRATENSW &
                                                            ,RTHRATENLW
      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlay),INTENT(INOUT) ::    QV   &
                                                                  ,QC   &
                                                                  ,QI   &
                                                                  ,QR   &
                                                                  ,QS   &
                                                                  ,QG
      REAL(kind_phys),DIMENSION(1:ncol),INTENT(INOUT) ::                &!ACFRCV,ACFRST  &
!     &                                                ,HBOT,HTOP        &
     &                                                ,RLWIN,RLWTOA     &
     &                                                ,RSWIN,RSWOUT     &
     &                                                ,GSW              &
     &                                                ,RSWINC
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(INOUT) ::         &
     &                                                        ,RLWTT    &
     &                                                        ,RSWTT
      REAL(kind_phys),DIMENSION(1:ncol),INTENT(INOUT) ::                &
     &                                                        CZMEAN    &
!     &                                                        ,SIGT4
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(INOUT) :: CLDFRA
      REAL(kind_phys), DIMENSION(1:ncol), INTENT(INOUT) ::              &
                            ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC,          &
                            ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,          &
                            ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC,          &
                            ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC

      ! TOA and surface, upward and downward, total and clear fluxes
      REAL(kind_phys), DIMENSION(1:ncol), INTENT(INOUT) ::              &
                              SWUPT,  SWUPTC,  SWDNT,  SWDNTC,          &
                              SWUPB,  SWUPBC,  SWDNB,  SWDNBC,          &
                              LWUPT,  LWUPTC,  LWDNT,  LWDNTC,          &
                              LWUPB,  LWUPBC,  LWDNB,  LWDNBC
      REAL(kind_phys), DIMENSION(1:ncol), INTENT(OUT  )  ::   SWVISDIR, &
                                                              SWVISDIF, &
                                                              SWNIRDIR, &
                                                              SWNIRDIF
      !..Additions for coupling cloud physics effective radii and radiation.  
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(IN):: re_cloud,   &
                                                            re_ice,     &
                                                            re_snow
      INTEGER, INTENT(INOUT):: has_reqc, has_reqi, has_reqs

      !..Output daily longwave and shortwave heating rates: 
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay),INTENT(OUT):: hrlwpd,hrswpd
!mz      LOGICAL,INTENT(IN) :: RESTRT
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER :: I,IENDX,II,ISTAT,J,JDAY,JMONTH,K,KMNTH,N,NRAD,         &
                 NUM_AEROSOLC
      INTEGER,DIMENSION(3) :: IDAT
      INTEGER,DIMENSION(12) :: MONTH=(/31,28,31,30,31,30,31,31          &
     &                                ,30,31,30,31/)
!
      REAL(kind_phys) :: CAPA,DAYI,FICE,FRAIN,GMT,HOUR,PLYR,            &
     &       ,QI,QR,QW,RADT,TIMES,WC,TDUM
!
      REAL(kind_phys),DIMENSION(1:nlay-1) :: QL,TL
!
      REAL(kind_phys),DIMENSION(1:ncol) :: CZEN   !mz,HBOTR,HTOPR      &
     &                                  ,SWCF, LWCF                     &
     &                                   ,REXNSFC,SWNETDN           &
     &                                  ,TOT,TOTLWDN,TOTSWDN,TOTSWDNC   &
     &                                  ,XLAND
      REAL(kind_phys),DIMENSION(1:ncol,1:nlay) :: CLFR,DZ               &
     &                                          ,P_PHY,PI_PHY       &
     &                                          ,RR,T8W                 &
     &                                          ,THRATENLW,THRATENSW    &
     &                                          ,TH_PHY,T_PHY
!

      REAL(kind_phys) :: DXKM, DYKM

      INTEGER            :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE

      integer,                   intent(in)    :: mpirank
      integer,                   intent(in)    :: mpiroot
      character(len=*),          intent(  out) :: errmsg
      integer,                   intent(  out) :: errflg


      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (.not.acswalloc .or. .not. aclwalloc) then
          write(errmsg, fmt='((a))') 'HWRF_radiation_run called before  &
                                      HWRF_radiation_init'
          errflg = 1
          return
      end if

      ! Set internal dimensions
      ids    =    1
      ims    =    1
      its    =    1
      myis   =    1
      myis1  =    1
      ide    =    ncol
      ime    =    ncol
      ite    =    ncol
      myie   =    ncol
      myie1  =    ncol
      jds    =    1
      jms    =    1
      jts    =    1
      myjs   =    1
      myjs2  =    1
      jde    =    1
      jme    =    1
      jte    =    1
      myje   =    1
      myje2  =    1
      kds    =    1
      kms    =    1
      kts    =    1
      kde    =    nlay
      kme    =    nlay
      kte    =    nlay

!----------------------------------------------------------------------
!***  RADIATION
!----------------------------------------------------------------------
!
!***  When allocating CAM radiation 4d arrays (ozmixm, aerosolc), 
!***  the following two scalars are not needed.
!
      NUM_AEROSOLC=1
!
!MZ      IF(grid%ntsd<=0)THEN
!MZ        NTSD_rad=grid%ntsd
!MZ      ELSE
!
!***  Call radiation just BEFORE the top of the hour
!***  so that updated fields are written to history files.
!
      IF (NTSD > 0) then
        NTSD_rad=ntsd+1
      ENDIF

      IF(MOD(NTSD_rad,GRID%NRADS)==0.OR.                               &
     &   MOD(NTSD_rad,GRID%NRADL)==0)THEN



       

     
!
!-----------------------------------------------------------------------
!***** NOTE: THIS IS HARDWIRED FOR CALLS TO LONGWAVE AND SHORTWAVE
!*****       AT EQUAL INTERVALS
!-----------------------------------------------------------------------
!
      NRAD=NRADS
      RADT=DT*NRADS/60.
!
!-----------------------------------------------------------------------
!
      CAPA=R_D/CP
!mz: wrf-nmm
!      DXKM=grid%dlmd*0.01745329*6371200. ! numbers from module_initialize_real.F
!      DYKM=grid%dy_nmm
! In FV3, dx in meter
       DXKM = dx*0.001
       DYKM = dx*0.001
!
!-----------------------------------------------------------------------
!
      DO J=MYJS2,MYJE2
      DO I=MYIS1,MYIE1
!
!mz        PDSL(I,J)=PD(I,J)*RES(I,J)
!mz        P8W(I,KTE+1,J)=PT
        XLAND(I,J)=SM(I,J)+1.
!mz        PSFC=PD(I,J)+PDTOP+PT
!mz        REXNSFC(I,J)=(PSFC*1.E-5)**CAPA
!mz        TSFC(I,J)=THS(I,J)*REXNSFC(I,J)
        T8W(I,KTS,J)=TSFC(I,J)
!MZ:eta1-interface sigma value in pressure domain 
!   eta2-interface sigma value in sigma domain
!  pdtop-mass at i,j in pressure domain, Pa
!  pdsl- sigma-domain pressure at sigma=1
!MZ        P8W(I,KTS,J)=ETA1(KTS)*PDTOP+ETA2(KTS)*PDSL(I,J)+PT
!prsi: air pressure at interface in Pa
!        P8W(I,KTS,J)=prsi(i,kts)
     
!        Z_PHY(I,KTS,J)=Z(I,J,KTS)
!        HT(I,J)=Z(I,J,KTS)
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
      DO J=MYJS2,MYJE2
      DO I=MYIS1,MYIE1
        DO K=KTS,KTE
          !mz DPL=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
          QL(K)=MAX(Q(I,J,K),EPSQ)
          !mz PLYR=AETA1(K)*PDTOP+AETA2(K)*PDSL(I,J)+PT
          PLYR = PRSL(I,J,K)
          TL(K)=T(I,J,K)
!
          RR(I,K,J)=PLYR/(R_D*TL(K)*(1.+P608*QL(K)))
          T_PHY(I,K,J)=TL(K)
          TH_PHY(I,K,J)=TL(K)*(1.E5/PLYR)**CAPA
          !mz P8W(I,K+1,J)=ETA1(K+1)*PDTOP+ETA2(K+1)*PDSL(I,J)+PT
          P_PHY(I,K,J)=PLYR
          PI_PHY(I,K,J)=(PLYR*1.E-5)**CAPA
          DZ(I,K,J)=TL(K)*(P608*QL(K)+1.)*R_D                           &
     &                 *(P8W(I,K,J)-P8W(I,K+1,J))                       &
     &                 /(P_PHY(I,K,J)*G)
!
          RTHRATEN(I,K,J)=0.
          THRATENLW(I,K,J)=0.
          THRATENSW(I,K,J)=0.

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
      DO K=KMS,KME
        DO J=JMS,JME
        DO I=IMS,IME
          CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
      ENDDO
!
      DO J=JMS,JME
        DO I=IMS,IME
!          CFRACH(I,J)=0.
!          CFRACL(I,J)=0.
!          CFRACM(I,J)=0.
          CZMEAN(I,J)=0.
!          SIGT4(I,J)=0.
          TOTSWDN(I,J)=0.   ! TOTAL (clear+cloudy sky) shortwave down at the surface
          TOTSWDNC(I,J)=0.  ! CLEAR SKY shortwave down at the surface
          SWNETDN(I,J)=0.   ! Net (down - up) total (clear+cloudy sky) shortwave at the surface
          TOTLWDN(I,J)=0.   ! Total longwave down at the surface
!mz          CUPPTR(I,J)=CUPPT(I,J)   ! Temporary array set to zero in radiation
!mz          HTOPR(I,J) =0.
!mz          HBOTR(I,J) = REAL(KTE+1)
!          HBOTR(I,J) =0.
          SWCF(I,J) =0.
          LWCF(I,J) =0.
!
!
        ENDDO
      ENDDO
!
!-------------------------------------------------------------------
!
!***  CALL THE INNER DRIVER.
!
!-----------------------------------------------------------------------

      CALL RADIATION_DRIVER(                                            &
     &                  IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=KDS,KDE=KDE &
     &                 ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=KMS,KME=KME &
     &                 ,KTS=KTS,KTE=KTE,                                &
     &                 ,ITIMESTEP=NTSD_rad,DT=DT                        &
     &                 ,DX=DXKM,DY=DYKM                                 &
     &                 ,swint_opt=swint_opt                             &
     &                 ,SWDDIR=swddir,SWDDNI=swddni,SWDDIF=swddif       &
     &                 ,SWDDIRC=swddirc,SWDDNIC=swddnic                 &
     &                 ,Gx=Gx,Bx=Bx,gg=gg,bb=bb                         &
     &                 ,swdown_ref=swdown_ref,swddir_ref=swddir_ref     &
     &                 ,coszen_ref=coszen_ref,coszen=coszen             &
     &                 ,hrang=hrang                                     &
     &                 ,RTHRATENLW=THRATENLW,RTHRATENSW=THRATENSW       &
     &                 ,RTHRATEN=RTHRATEN                               &
     &                 ,HRSWPD=HRSWPD, HRLWPD=HRLWPD                    &
     &                 ,GLW=TOTLWDN,GSW=SWNETDN,SWDOWN=TOTSWDN          &
     &                 ,XLAT=XLAT,XLONG=XLONG,ALBEDO=ALBEDO,EMISS=EPSR  &
     &                 ,XICE=SICE,XLAND=XLAND,Z=Z,TSK=TSFC              &
     &                 ,N_AEROSOLC=NUM_AEROSOLC,PAERLEV=PAERLEV         &
     &                 ,CAM_ABS_FREQ_S=CAM_ABS_FREQ_S                   &
     &                 ,ALEVSIZ=ALEVSIZ,no_src_types=no_src_types       &
     &                 ,LEVSIZ=LEVSIZ,N_OZMIXM=NUM_OZMIXM               &
     &                 ,OZMIXM=OZMIXM,PIN=PIN                           &
     &                 ,SNOW=SNOW                                       &
     &                 ,RHO=RR,P8W=P8W,P=P_PHY,PI=PI_PHY                &
     &                 ,DZ8W=DZ,T=T_PHY,T8W=T8W,GMT=GMT                 &
     &                 ,JULDAY=JULDAY,JULYR=JULYR,NPHS=NPHS             &
     &                 ,JULIAN=JULIAN,XTIME=XTIME                       &
     &                 ,YR=JULYR                                        &
     &                 ,RADT=RADT  ,RA_CALL_OFFSET=RA_CALL_OFFSET       &
     &                 ,STEPRA=NRAD  ,ICLOUD=config_flags%ICLOUD        &
     &                 ,cldovrlp=cldovrlp                               &
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
     &        ,SWVISDIR=swvisdir ,SWVISDIF=swvisdif                     &
     &        ,SWNIRDIR=swnirdir ,SWNIRDIF=swnirdif                     &
     &                 ,re_cloud=grid%re_cloud                          &
     &                 ,re_ice=grid%re_ice                              &
     &                 ,re_snow=grid%re_snow                            &
     &                 ,has_reqc=has_reqc                               &
     &                 ,has_reqi=has_reqi                               &
     &                 ,has_reqs=has_reqs                               &
     &                 ,RLWTOA=RLWTOA                     &
     &                 ,CZMEAN=CZMEAN                                   &
!     &                 ,ACFRST=ACFRST
!     &                 ,NCFRST=NCFRST                     &
!     &                 ,ACFRCV=ACFRCV
!     &                 ,NCFRCV=NCFRCV                     &
     &                 ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN               &
     &                 ,LWCF=LWCF,SWCF=SWCF                             &
     &                 ,O3INPUT=O3INPUT,AER_OPT=AER_OPT,O3RAD=o3rad     &
     &                 ,SF_SURFACE_PHYSICS=SF_SURFACE_PHYSICS           &
     &                 ,QV,F_QV=F_QV      &
     &                 ,QC,F_QC=F_QC      &
     &                 ,QR,F_QR=F_QR      &
     &                 ,QI,F_QI=F_QI      &
     &                 ,QS,F_QS=F_QS      &
     &                 ,QG,F_QG=F_QG      &
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
      nrads_block: IF(MOD(NTSD_rad,NRADS)==0)THEN
!-----------------------------------------------------------------------
!***  COMPUTE CZMEAN FOR NON-GFDL SHORTWAVE
!-----------------------------------------------------------------------
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
            TIMES=NTSD_rad*DT+II*DT
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
!MZ            IF(HBM2(I,J)>0.5)THEN
              TOTSWDN(I,J)=SWNETDN(I,J)/(1.-ALBEDO(I,J))
!
!--- No value currently available for clear-sky solar fluxes from
!    non GFDL schemes, though it's needed for air quality forecasts.
!    For the time being, set to the total downward solar fluxes.
!
              TOTSWDNC(I,J)=TOTSWDN(I,J)
!MZ            ENDIF
!
          ENDDO
          ENDDO
!
!        ENDIF   !End non-GFDL block
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
      nradl_block: IF(MOD(NTSD_rad,NRADL)==0)THEN
!
!$omp parallel do                                                       &
!$omp& private(i,iendx,j)
        DO J=MYJS2,MYJE2
          IENDX=MYIE1
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=MYIS1,IENDX
!
!MZ: hbm2=1. in global FV3
!            IF(HBM2(I,J)>0.5)THEN
              TDUM=T(I,J,KTS)
!              SIGT4(I,J)=STBOLT*TDUM*TDUM*TDUM*TDUM
              RLWIN(I,J)=TOTLWDN(I,J)
!            ENDIF
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
!MZ            IF(HBM2(I,J)>0.5)THEN
                RLWTT(I,J,K)=THRATENLW(I,K,J)*PI_PHY(I,K,J)
!MZ            ENDIF
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
      DO J=MYJS2,MYJE2
        IENDX=MYIE1
        IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
        DO I=MYIS1,IENDX
          HBOT(I,J)=HBOTR(I,J)
!          HTOP(I,J)=HTOPR(I,J)
!          CUPPT(I,J)=CUPPTR(I,J)
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  ZERO OUT BOUNDARY ROWS.
!-----------------------------------------------------------------------
!
!MZ
!      DO J=JTS,JTE
!      DO I=ITS,ITE
!        IF(HBM2(I,J)<0.5)THEN
!          ACFRST(I,J)=0.
!          ACFRCV(I,J)=0.
!!          CFRACL(I,J)=0.
!!          CFRACM(I,J)=0.
!!          CFRACH(I,J)=0.
!          RSWTOA(I,J)=0.
!          RLWTOA(I,J)=0.
!        ENDIF
!      ENDDO
!      ENDDO

!
      DO J=jts,min(jde-1,jte)
      DO I=its,min(ide-1,ite)
          gsw(I,J)=rswin(I,J)-rswout(I,J)
      ENDDO
      ENDDO
     
      ENDIF    !radiation loop

!
!-----------------------------------------------------------------------
!
      END SUBROUTINE HWRF_radiation_run

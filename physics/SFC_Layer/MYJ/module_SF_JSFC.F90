!-----------------------------------------------------------------------
!
      MODULE MODULE_SF_JSFC
!
!-----------------------------------------------------------------------
!
!***  THE J SURFACE SCHEME
!
!-----------------------------------------------------------------------
!
!      USE MODULE_INCLUDE
!
!      USE MODULE_CONSTANTS,ONLY : A2,A3,A4,CP,ELWV                      &
!                                 ,G,P608,PI,PQ0,R_D,R_V,CAPPA
!
!-----------------------------------------------------------------------
!

      USE machine, only: kfpt => kind_phys

      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!       integer,parameter :: isingle=selected_int_kind(r=9)
!       integer,parameter :: idouble=selected_int_kind(r=18)
!       integer,parameter :: single=selected_real_kind(p=6,r=37)
!       integer,parameter :: double=selected_real_kind(p=13,r=200)
!
!       integer,parameter:: &
!       klog=4 &
!      ,kint=isingle &
!      ,kdin=idouble &
!      ,kfpt=single &
!      ,kdbl=double

!       real   (kind=kfpt),parameter :: r4_in=x'ffbfffff'
!       real   (kind=kdbl),parameter :: r8_in=x'fff7ffffffffffff'
!       integer(kind=kint),parameter :: i4_in=-999  ! -huge(1)
!
      ! integer,parameter:: &
      !   klog=4 &                   ! logical variables
      !  ,kint=4 &                   ! integer variables
      !  !,kfpt=4 &                   ! floating point variables
      !  ,kfpt=8 &                   ! floating point variables
      !  ,kdbl=8                     ! double precision
!
      PRIVATE
!
      PUBLIC :: JSFC_INIT,JSFC
!
      INTEGER :: ITRMX=5 ! Iteration count for mixing length computation
!
      REAL(kind=kfpt),PARAMETER :: VKARMAN=0.4

      REAL(kind=kfpt),PARAMETER :: A2=17.2693882,A3=273.15,A4=35.86,CP=1004.6      &
                       ,ELWV=2.501e6,EPSQ2=0.02,G=9.8060226             &
                       ,PQ0=379.90516,R_D=287.04,R_V=461.6              &
                       ,P608=R_V/R_D-1.,CAPPA=R_D/CP                    &
                       ,PI=3.141592653589793

      REAL(kind=kfpt),PARAMETER :: XLV=ELWV
      REAL(kind=kfpt),PARAMETER :: ELOCP=2.72E6/CP
      REAL(kind=kfpt),PARAMETER :: A2S=17.2693882,A3S=273.16,A4S=35.86
      REAL(kind=kfpt),PARAMETER :: GLKBR=10.,GLKBS=30.                             &
                       ,QVISC=2.1E-5,RIC=0.505,SMALL=0.35               &
                       ,SQPR=0.84,SQSC=0.84,SQVISC=258.2                &
                       ,TVISC=2.1E-5                                    &
                       ,USTC=0.7,USTR=0.225,VISC=1.5E-5                 &
                       ,WWST=1.2,ZTFC=1.
      REAL(kind=kfpt),PARAMETER :: SEAFC=0.98,PQ0SEA=PQ0*SEAFC

      REAL(kind=kfpt),PARAMETER :: CZIV=SMALL*GLKBS,GRRS=GLKBR/GLKBS

      REAL(kind=kfpt),PARAMETER :: RTVISC=1./TVISC,RVISC=1./VISC                   &
                       ,ZQRZT=SQSC/SQPR

      REAL(kind=kfpt),PARAMETER :: USTFC=0.018/G                                   &
                       ,FZQ1=RTVISC*QVISC*ZQRZT                         &
                       ,FZQ2=RTVISC*QVISC*ZQRZT                         &
                       ,FZT1=RVISC *TVISC*SQPR                          &
                       ,FZT2=CZIV*GRRS*TVISC*SQPR                       &
                       ,FZU1=CZIV*VISC
      REAL(kind=kfpt),PARAMETER :: WWST2=WWST*WWST                                 &
                       ,RQVISC=1./QVISC

      REAL(kind=kfpt),PARAMETER :: RCAP=1./CAPPA
      REAL(kind=kfpt),PARAMETER :: GOCP02=G/CP*2.,GOCP10=G/CP*10.
      REAL(kind=kfpt),PARAMETER :: EPSU2=1.E-6,EPSUST=1.E-9,EPSZT=1.E-28
      REAL(kind=kfpt),PARAMETER :: CZIL=0.1,EXCML=0.0001,EXCMS=0.0001              &
     &                 ,FH=1.10,TOPOFAC=9.0e-6

      REAL(kind=kfpt),PARAMETER :: ZILFC=-CZIL*VKARMAN*SQVISC
      REAL(kind=kfpt),PARAMETER :: EPSQ=1.e-9
!
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: KZTM=10001,KZTM2=KZTM-2
!
      REAL(kind=kfpt),PRIVATE,SAVE ::      &
         DZETA1,DZETA2,FH01,FH02,ZTMAX1,ZTMAX2,ZTMIN1,ZTMIN2
!
      REAL(kind=kfpt),DIMENSION(KZTM),PRIVATE,SAVE ::    &
                   PSIH1,PSIH2,PSIM1,PSIM2
!
      INTEGER :: IERR
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE JSFC(FLAG_ITER,ITER,ME                                &
     &               ,NTSD,EPSQ2,HT,DZ                                 &
     &               ,PHMID,PHINT,TH,T,Q,QC,U,V,Q2                     &
     &               ,TSK,QSFC,THZ0,QZ0,UZ0,VZ0                        &
     &               ,XLAND                                            &
     &               ,USTAR,Z0,Z0BASE,PBLH,MAVAIL,RMOL                 &
     &               ,AKHS,AKMS,CHKLOWQ,HLFLX,RIB                      &
     &               ,CM,CH,STRESS,FFM,FFH,WIND,FM10,FH2               &
     &               ,A1U,A1T,A1Q                                      &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,LM,errmsg,errflg)
!
!-----------------------------------------------------------------------
!      SUBROUTINE JSFC(NTSD,EPSQ2,HT,DZ                                 &
!     &               ,PHMID,PHINT,TH,T,Q,QC,U,V,Q2                     &
!     &               ,TSK,QSFC,THZ0,QZ0,UZ0,VZ0                        &
!     &               ,XLAND                                            &
!     &               ,VEGFRC,SNOWC                                     & !added 5/17/2013
!     &               ,USTAR,Z0,Z0BASE,PBLH,MAVAIL,RMOL                 &
!     &               ,AKHS,AKMS                                        &
!     &               ,CHS,CHS2,CQS2,HFX,QFX,FLX_LH,FLHC,FLQC           &
!     &               ,QGH,CPM,CT                                       &
!     &               ,U10,V10,T02,TH02,TSHLTR,TH10,Q02,QSHLTR,Q10      &
!     &               ,PSHLTR,RIB                                       & ! Added Bulk Richardson No.
!     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
!     &               ,IMS,IME,JMS,JME,KMS,KME                          &
!     &               ,ITS,ITE,JTS,JTE,KTS,LM)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,LM
!
      INTEGER,INTENT(IN) :: NTSD,ITER,ME
      LOGICAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FLAG_ITER
      real(kind=kfpt),dimension(1:lm),intent(in):: epsq2
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: HT,MAVAIL,TSK      &
     &                                             ,XLAND,Z0BASE
!     &                                             ,VEGFRC,SNOWC
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: DZ,PHMID
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM+1),INTENT(IN) :: PHINT
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: Q,QC,U,V,Q2,T,TH
!
!      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: FLX_LH,HFX,PSHLTR &
!     &                                              ,QFX,Q10,QSHLTR    &
!     &                                              ,TH10,TSHLTR,T02   &
!     &                                              ,U10,V10,TH02,Q02
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME) :: FLX_LH,HFX &
     &                                              ,QFX,Q10,TH10,T02   &
     &                                              ,U10,V10,TH02,Q02
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME) :: PSHLTR,QSHLTR,TSHLTR
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: AKHS,AKMS       &
     &                                                ,PBLH,QSFC,RIB
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: QZ0,RMOL,THZ0   &
     &                                                ,USTAR,UZ0,VZ0   &
     &                                                ,Z0
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: HLFLX,CHKLOWQ
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: CM,CH,STRESS,FFM  &
     &                                              ,FFH,WIND,FM10,FH2 &
     &                                              ,A1U,A1T,A1Q
      character(len=*),intent(out) :: errmsg
      integer,         intent(out) :: errflg
!
!      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: CHS,CHS2,CQS2     &
!     &                                              ,CPM,CT,FLHC,FLQC  &
!     &                                              ,QGH
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME) :: CHS,CHS2,CQS2     &
     &                                              ,FLHC,FLQC
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME) :: QGH,CPM,CT
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: I,J,K,LMH,LPBL
!
      REAL(kind=kfpt) :: A,APESFC,B,BTGX,CWMLOW                                   &
     &       ,DQDT,DTDIF,DTDT,DUDT,DVDT                                &
     &       ,FIS                                                      &
     &       ,P02P,P10P,PLOW,PSFC,PTOP,QLOW,QS02,QS10                  &
     &       ,RAPA,RAPA02,RAPA10,RATIOMX,RDZ,SEAMASK,SM                &
     &       ,T02P,T10P,TEM,TH02P,TH10P,THLOW,THELOW,THM               &
     &       ,TLOW,TZ0,ULOW,VLOW,ZSL
!
      REAL(kind=kfpt),DIMENSION(KTS:LM) :: CWMK,PK,Q2K,QK,THEK,THK,TK,UK,VK
!
      REAL(kind=kfpt),DIMENSION(KTS:LM-1) :: EL,ELM
!
      REAL(kind=kfpt),DIMENSION(KTS:LM+1) :: ZHK
!
      REAL(kind=kfpt),DIMENSION(ITS:ITE,JTS:JTE) :: THSK
!
      REAL(kind=kfpt),DIMENSION(ITS:ITE,JTS:JTE,KTS:LM+1) :: ZINT
!
!----------------------------------------------------------------------
!**********************************************************************
      ! Initialize error-handling
      errflg = 0
      errmsg = ''
!----------------------------------------------------------------------
!
!***  MAKE PREPARATIONS
!
!----------------------------------------------------------------------
      DO J=JTS,JTE
      DO I=ITS,ITE
         IF(FLAG_ITER(I,J))THEN
            DO K=KTS,LM+1
               ZINT(I,J,K)=0.
            ENDDO
         END IF
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        IF(FLAG_ITER(I,J))THEN
        ZINT(I,J,LM+1)=HT(I,J)     ! Z at bottom of lowest sigma layer
        PBLH(I,J)=-1.
!
!!!!!!!!!
!!!!!! UNCOMMENT THESE LINES IF USING ETA COORDINATES
!!!!!!!!!
!!!!!!  ZINT(I,J,LM+1)=1.E-4       ! Z of bottom of lowest eta layer
!!!!!!  ZHK(LM+1)=1.E-4            ! Z of bottom of lowest eta layer
!
        END IF
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
        DO I=ITS,ITE
          IF(FLAG_ITER(I,J))THEN
          DO K=LM,KTS,-1
            ZINT(I,J,K)=ZINT(I,J,K+1)+DZ(I,J,K)
          ENDDO
          END IF
      ENDDO
      ENDDO
!
      IF(NTSD==0.and.iter==1) then
        DO J=JTS,JTE
        DO I=ITS,ITE
          IF(FLAG_ITER(I,J))THEN
          USTAR(I,J)=0.1
          FIS=HT(I,J)*G
          SM=XLAND(I,J)-1.
!!!       Z0(I,J)=SM*Z0SEA+(1.-SM)*(Z0(I,J)*Z0MAX+FIS*FCM+Z0LAND)
!!!       Z0(I,J)=SM*Z0SEA+(1.-SM)*(Z0(I,J)*Z0MAX+FIS*FCM+Z0LAND)
          END IF
        ENDDO
        ENDDO
      ENDIF
!
!!!!  IF(NTSD==1)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          CT(I,J)=0.
        ENDDO
        ENDDO
!!!!  ENDIF
!
!......................................................................
!$omp parallel do &
!$omp private (j,i,lmh,ptop,psfc,seamask,k,thk,tk,ratiomx,qk,pk, &
!$omp          cwmk,thek,q2k,zhk,uk,vk,lpbl,plow,tlow,thlow,thelow,    &
!$omp          qlow,cwmlow,ulow,vlow,zsl,apesfc,tz0,rapa,th02p,th10p,  &
!$omp          rapa02,rapa10,t02p,t10p,p02p,p10p,qs02,qs10)
!......................................................................
!----------------------------------------------------------------------
        setup_integration:  DO J=JTS,JTE
!----------------------------------------------------------------------
!
        DO I=ITS,ITE
          IF(FLAG_ITER(I,J))THEN
!
!***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
!
          LMH=LM
!
          PTOP=PHINT(I,J,1)
          PSFC=PHINT(I,J,LMH+1)
! Define THSK here (for first timestep mostly)
          THSK(I,J)=TSK(I,J)/(PSFC*1.E-5)**CAPPA
!
!***  CONVERT LAND MASK (1 FOR SEA; 0 FOR LAND)
!
          SEAMASK=XLAND(I,J)-1.
!
!***  FILL 1-D VERTICAL ARRAYS
!
          DO K=LM,KTS,-1
            THK(K)=TH(I,J,K)
            TK(K)=T(I,J,K)
            QK(K)=Q(I,J,K)
            PK(K)=PHMID(I,J,K)
            CWMK(K)=QC(I,J,K)
            THEK(K)=(CWMK(K)*(-ELOCP/TK(K))+1.)*THK(K)
            Q2K(K)=Q2(I,J,K)
!
!
!***  COMPUTE THE HEIGHTS OF THE LAYER INTERFACES
!
            ZHK(K)=ZINT(I,J,K)
!
          ENDDO
          ZHK(LM+1)=HT(I,J)          ! Z at bottom of lowest sigma layer
!
          DO K=LM,KTS,-1
            UK(K)=U(I,J,K)
            VK(K)=V(I,J,K)
          ENDDO
!
!***  FIND THE HEIGHT OF THE PBL
!
          LPBL=LMH
          DO K=LMH-1,1,-1
            IF(Q2K(K)<=EPSQ2(K)*FH) THEN
              LPBL=K
              GO TO 110
            ENDIF
          ENDDO
!
          LPBL=1
!
!-----------------------------------------------------------------------
!--------------THE HEIGHT OF THE PBL------------------------------------
!-----------------------------------------------------------------------
!
 110      PBLH(I,J)=ZHK(LPBL)-ZHK(LMH+1)
!
!----------------------------------------------------------------------
          IF(QC(I,J,LM).GT.EPSQ)THEN
             CHKLOWQ(I,J)=0.
          ELSE
             CHKLOWQ(I,J)=1.
          ENDIF
!***
!***  FIND THE SURFACE EXCHANGE COEFFICIENTS
!***
!----------------------------------------------------------------------
          PLOW=PK(LMH)
          TLOW=TK(LMH)
          THLOW=THK(LMH)
          THELOW=THEK(LMH)
          QLOW=QK(LMH)
          CWMLOW=CWMK(LMH)
          ULOW=UK(LMH)
          VLOW=VK(LMH)
          ZSL=(ZHK(LMH)-ZHK(LMH+1))*0.5
!          if(me.eq.0)print*,'ZSL,ZHK(LMH),ZHK(LMH+1,LMH=',ZSL,ZHK(LMH),ZHK(LMH+1),LMH
          APESFC=(PSFC*1.E-5)**CAPPA
       if(NTSD==0) then
          TZ0=TSK(I,J)
       else
          TZ0=THZ0(I,J)*APESFC
       endif
!
         CALL SFCDIF(NTSD,SEAMASK,THSK(I,J),QSFC(I,J),PSFC             &
     &               ,UZ0(I,J),VZ0(I,J),TZ0,THZ0(I,J),QZ0(I,J)         &
     &               ,USTAR(I,J),Z0(I,J),Z0BASE(I,J),CT(I,J),RMOL(I,J) &
     &               ,AKMS(I,J),AKHS(I,J),PBLH(I,J),MAVAIL(I,J)        &
     &               ,CHS(I,J),CHS2(I,J),CQS2(I,J)                     &
     &               ,HFX(I,J),QFX(I,J),FLX_LH(I,J)                    &
     &               ,FLHC(I,J),FLQC(I,J),QGH(I,J),CPM(I,J)            &
     &               ,ULOW,VLOW,TLOW,THLOW,THELOW,QLOW,CWMLOW          &
     &               ,ZSL,PLOW,HLFLX(I,J)                              &
!     &               ,VEGFRC(I,J),SNOWC(I,J)                           &   !added 5/17/2013
     &               ,U10(I,J),V10(I,J),TSHLTR(I,J),TH10(I,J)          &
     &               ,QSHLTR(I,J),Q10(I,J),PSHLTR(I,J)                 &
     &               ,FFM(I,J),FFH(I,J),FM10(I,J),FH2(I,J)             &
     &               ,A1U(I,J),A1T(I,J),A1Q(I,J)                       &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,LM,I,J,ZHK(LMH+1),RIB(I,J)   & ! Added Bulk Richardson No.
     &               ,errmsg, errflg)
!
!***  REMOVE SUPERATURATION AT 2M AND 10M
!
          RAPA=APESFC
          TH02P=TSHLTR(I,J)
          TH10P=TH10(I,J)
          TH02(I,J)=TSHLTR(I,J)
!
          RAPA02=RAPA-GOCP02/TH02P
          RAPA10=RAPA-GOCP10/TH10P
!
          T02P=TH02P*RAPA02
          T10P=TH10P*RAPA10
! 1 may 06 tgs         T02(I,J) = T02P
          T02(I,J) = TH02(I,J)*APESFC
!
          P02P=(RAPA02**RCAP)*1.E5
          P10P=(RAPA10**RCAP)*1.E5
!
          QS02=PQ0/P02P*EXP(A2*(T02P-A3)/(T02P-A4))
          QS10=PQ0/P10P*EXP(A2*(T10P-A3)/(T10P-A4))
!
          IF(QSHLTR(I,J)>QS02)QSHLTR(I,J)=QS02
          IF(Q10   (I,J)>QS10)Q10   (I,J)=QS10
          Q02(I,J)=QSHLTR(I,J)/(1.-QSHLTR(I,J))
!----------------------------------------------------------------------
!          STRESS(I,J)=USTAR(I,J)*USTAR(I,J)
          WIND(I,J)=max(USTAR(I,J)*FFM(I,J)/VKARMAN,1.0)
          CM(I,J)=VKARMAN*VKARMAN/(FFM(I,J)*FFM(I,J))
          CH(I,J)=VKARMAN*VKARMAN/(FFM(I,J)*FFH(I,J))
          TEM=0.00001/DZ(I,J,LM)
          CM(I,J)=max(CM(I,J),tem)
          CH(I,J)=max(CH(I,J),tem)
          STRESS(I,J)=cm(I,J) * wind(I,J) * wind(I,J)
          USTAR(I,J)=sqrt(stress(I,J))
!
          END IF    ! FLAG_ITER
!
        ENDDO
!
!----------------------------------------------------------------------
!
      ENDDO setup_integration
!
!......................................................................
!$omp end parallel do
!......................................................................
!----------------------------------------------------------------------

      END SUBROUTINE JSFC
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
      SUBROUTINE SFCDIF(NTSD,SEAMASK,THS,QS,PSFC                       &
     &                 ,UZ0,VZ0,TZ0,THZ0,QZ0                           &
     &                 ,USTAR,Z0,Z0BASE,CT,RLMO,AKMS,AKHS,PBLH,WETM    &
     &                 ,CHS,CHS2,CQS2,HFX,QFX,FLX_LH,FLHC,FLQC,QGH,CPM &
     &                 ,ULOW,VLOW,TLOW,THLOW,THELOW,QLOW,CWMLOW        &
     &                 ,ZSL,PLOW,HLFLX                                 &
!     &                 ,VEGF,SNOC                                      &  !added 5/17/2013
     &                 ,U10,V10,TH02,TH10,Q02,Q10,PSHLTR               &
     &                 ,FFM,FFH,FM10,FH2,A1U,A1T,A1Q                   &
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                 ,IMS,IME,JMS,JME,KMS,KME                        &
     &                 ,ITS,ITE,JTS,JTE,KTS,LM,I,J,ZSFC,RIB            & ! Added Bulk Richardson No.
     &                 ,errmsg, errflg) 
!     ****************************************************************
!     *                                                              *
!     *                       SURFACE LAYER                          *
!     *                                                              *
!     ****************************************************************
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,LM,i,j
!
      INTEGER,INTENT(IN) :: NTSD
!
      REAL(kind=kfpt),INTENT(IN) :: CWMLOW,PBLH,PLOW,QLOW,PSFC,SEAMASK,ZSFC       &
     &                  ,THELOW,THLOW,THS,TLOW,TZ0,ULOW,VLOW,WETM,ZSL  &
     &                  ,Z0BASE
!                       ,VEGF,SNOC
!
      REAL(kind=kfpt),INTENT(OUT) :: CHS,CHS2,CPM,CQS2,CT,FLHC,FLQC,FLX_LH,HFX    &
     &              ,RIB,PSHLTR,Q02,Q10,QFX,QGH,RLMO,TH02,TH10,U10,V10
      REAL(kind=kfpt),INTENT(OUT) :: FFM,FFH,FM10,FH2,A1U,A1T,A1Q
!
      REAL(kind=kfpt),INTENT(INOUT) :: AKHS,AKMS,QZ0,THZ0,USTAR,UZ0,VZ0,Z0,QS
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: ITR,K
!
      REAL(kind=kfpt) :: A,B,BTGH,BTGX,CXCHL,CXCHS,DTHV,DU2,ELFC,FCT              &
     &       ,HLFLX,HSFLX,HV,PSH02,PSH10,PSHZ,PSHZL,PSM10,PSMZ,PSMZL   &
     &       ,RDZ,RDZT,RLMA,RLMN,RLMP                                  &
     &       ,RLOGT,RLOGU,RWGH,RZ,RZST,RZSU,SIMH,SIMM,TEM,THM          &
     &       ,UMFLX,USTARK,VMFLX,WGHT,WGHTT,WGHTQ,WSTAR2               &
     &       ,X,XLT,XLT4,XLU,XLU4,XT,XT4,XU,XU4,ZETALT,ZETALU          &
     &       ,ZETAT,ZETAU,ZQ,ZSLT,ZSLU,ZT,ZU,TOPOTERM,ZZIL,CZETMAX
!vcw
!
!***  DIAGNOSTICS
!
      REAL(kind=kfpt) :: AKHS02,AKHS10,AKMS02,AKMS10,EKMS10,QSAT10,QSAT2          &
     &       ,RLNT02,RLNT10,RLNU10,SIMH02,SIMH10,SIMM10,T02,T10        &
     &       ,TERM1,RLOW,U10E,V10E,WSTAR,XLT02,XLT024,XLT10            &
     &       ,XLT104,XLU10,XLU104,XU10,XU104,ZT02,ZT10,ZTAT02,ZTAT10   &
     &       ,ZTAU,ZTAU10,ZU10,ZUUZ
!      REAL(kind=kfpt) :: ZILFC1,SNOWZO, Zom_ztmax
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      ! Initialize error-handling
      errflg = 0
      errmsg = ''

      RDZ=1./ZSL
      CXCHL=EXCML*RDZ
      CXCHS=EXCMS*RDZ
!
      BTGX=G/THLOW
      ELFC=VKARMAN*BTGX
!
      IF(PBLH>1000.)THEN
        BTGH=BTGX*PBLH
      ELSE
        BTGH=BTGX*1000.
      ENDIF
!
      WGHT=0.
      WGHTT=0.
      WGHTQ=0.
!----------------------------------------------------------------------
!
!***  SEA POINTS
!
!----------------------------------------------------------------------
!
      IF(SEAMASK>0.5)THEN
!
!----------------------------------------------------------------------
        DO ITR=1,ITRMX
!----------------------------------------------------------------------
          Z0=MAX(USTFC*USTAR*USTAR,1.59E-5)
!
!***  VISCOUS SUBLAYER, JANJIC MWR 1994
!
!----------------------------------------------------------------------
          IF(USTAR<USTC)THEN
!----------------------------------------------------------------------
!
            IF(USTAR<USTR)THEN
!
              IF(NTSD==0) then
!tgs              IF(NTSD==1)THEN
                AKMS=CXCHS
                AKHS=CXCHS
                QS=QLOW
              ENDIF
!
              ZU=FZU1*SQRT(SQRT(Z0*USTAR*RVISC))/USTAR
              WGHT=AKMS*ZU*RVISC
              RWGH=WGHT/(WGHT+1.)
              UZ0=(ULOW*RWGH+UZ0)*0.5
              VZ0=(VLOW*RWGH+VZ0)*0.5
!
              ZT=FZT1*ZU
              ZQ=FZQ1*ZT
              WGHTT=AKHS*ZT*RTVISC
              WGHTQ=AKHS*ZQ*RQVISC
!
              IF(NTSD>0)THEN
                THZ0=((WGHTT*THLOW+THS)/(WGHTT+1.)+THZ0)*0.5
                QZ0=((WGHTQ*QLOW+QS)/(WGHTQ+1.)+QZ0)*0.5
              ELSE
                THZ0=(WGHTT*THLOW+THS)/(WGHTT+1.)
                QZ0=(WGHTQ*QLOW+QS)/(WGHTQ+1.)
              ENDIF
!
            ENDIF
!
            IF(USTAR>=USTR.AND.USTAR<USTC)THEN
              ZU=Z0
              UZ0=0.
              VZ0=0.
!
              ZT=FZT2*SQRT(SQRT(Z0*USTAR*RVISC))/USTAR
              ZQ=FZQ2*ZT
              WGHTT=AKHS*ZT*RTVISC
              WGHTQ=AKHS*ZQ*RQVISC
!
              IF(NTSD>0)THEN
                THZ0=((WGHTT*THLOW+THS)/(WGHTT+1.)+THZ0)*0.5
                QZ0=((WGHTQ*QLOW+QS)/(WGHTQ+1.)+QZ0)*0.5
              ELSE
                THZ0=(WGHTT*THLOW+THS)/(WGHTT+1.)
                QZ0=(WGHTQ*QLOW+QS)/(WGHTQ+1.)
              ENDIF
!
            ENDIF
!----------------------------------------------------------------------
          ELSE
!----------------------------------------------------------------------
            ZU=Z0
            UZ0=0.
            VZ0=0.
!
            ZT=Z0
            THZ0=THS
!
            ZQ=Z0
            QZ0=QS
!----------------------------------------------------------------------
          ENDIF
!----------------------------------------------------------------------
          TEM=(TLOW+TZ0)*0.5
          THM=(THELOW+THZ0)*0.5
!
          A=THM*P608
          B=(ELOCP/TEM-1.-P608)*THM
!
          DTHV=((THELOW-THZ0)*((QLOW+QZ0+CWMLOW)*(0.5*P608)+1.)        &
     &        +(QLOW-QZ0+CWMLOW)*A+CWMLOW*B)
!
          DU2=MAX((ULOW-UZ0)**2+(VLOW-VZ0)**2,EPSU2)
          RIB=BTGX*DTHV*ZSL/DU2
!----------------------------------------------------------------------
!         IF(RIB>=RIC)THEN
!----------------------------------------------------------------------
!           AKMS=MAX( VISC*RDZ,CXCHS)
!           AKHS=MAX(TVISC*RDZ,CXCHS)
!----------------------------------------------------------------------
!         ELSE  !  turbulent branch
!----------------------------------------------------------------------
            ZSLU=ZSL+ZU
            ZSLT=ZSL+ZT
!
            RZSU=ZSLU/ZU
            RZST=ZSLT/ZT
!
            RLOGU=LOG(RZSU)
            RLOGT=LOG(RZST)
!
!----------------------------------------------------------------------
!***  1./MONIN-OBUKHOV LENGTH
!----------------------------------------------------------------------
!
            RLMO=ELFC*AKHS*DTHV/USTAR**3
!
            ZETALU=ZSLU*RLMO
            ZETALT=ZSLT*RLMO
            ZETAU=ZU*RLMO
            ZETAT=ZT*RLMO
!
            ZETALU=MIN(MAX(ZETALU,ZTMIN1),ZTMAX1)
            ZETALT=MIN(MAX(ZETALT,ZTMIN1),ZTMAX1)
            ZETAU=MIN(MAX(ZETAU,ZTMIN1/RZSU),ZTMAX1/RZSU)
            ZETAT=MIN(MAX(ZETAT,ZTMIN1/RZST),ZTMAX1/RZST)
!
!----------------------------------------------------------------------
!***   WATER FUNCTIONS
!----------------------------------------------------------------------
!
            RZ=(ZETAU-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZ=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
            RZ=(ZETALU-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZL=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
            SIMM=PSMZL-PSMZ+RLOGU
!
            RZ=(ZETAT-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZ=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
            RZ=(ZETALT-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZL=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
            SIMH=(PSHZL-PSHZ+RLOGT)*FH01
!----------------------------------------------------------------------
            USTARK=USTAR*VKARMAN
      if(abs(simm)<1.e-10.or.abs(simh)<1.e-10)then
        write(0,*)' simm=',simm,' simh=',simh,' at i=',i,' j=',j
      endif

!            if(abs(SIMM).lt.1.e-5.or.abs(SIMM).gt.1.e5)then
            if(abs(SIMM).lt.1.e-10.or.abs(SIMH).lt.1.e-10)then
              print*,'SIMM,PSMZL,PSMZ,RLOGU=',SIMM,PSMZL,PSMZ,RLOGU
              print*,'SIMH,PSHZL,PSHZ,RLOGT,FH01=',SIMH,PSHZL,PSHZ,RLOGT,FH01
              print*,'USTARK,CXCHS=',USTARK,CXCHS
              print*,'PSIM1(1,2),K=',PSIM1(K+1),PSIM1(K+2),K
              print*,'ZETAU,ZTMIN1,DZETA1=',ZETAU,ZTMIN1,DZETA1
              print*,'PSIH1(1,2),RDZT=',PSIH1(K+1),PSIH1(K+2),RDZT
              print*,'ZSLU,ZSLT,RLMO,ZU,ZT=',ZSLU,ZSLT,RLMO,ZU,ZT
              print*,'A,B,DTHV,DU2,RIB=',A,B,DTHV,DU2,RIB
              errflg = 1
              errmsg = 'ERROR(SFCDIF): '
              return
            end if



            AKMS=MAX(USTARK/SIMM,CXCHS)
            AKHS=MAX(USTARK/SIMH,CXCHS)
!
!----------------------------------------------------------------------
!***  BELJAARS CORRECTION FOR USTAR
!----------------------------------------------------------------------
!
            IF(DTHV<=0.)THEN                                           !zj
              WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)                !zj
            ELSE                                                       !zj
              WSTAR2=0.                                                !zj
            ENDIF                                                      !zj
            USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
!
!----------------------------------------------------------------------
!         ENDIF  !  End of turbulent branch
!----------------------------------------------------------------------
!
        ENDDO  !  End of the iteration loop over sea points
!
!----------------------------------------------------------------------
!
!***  LAND POINTS
!
!----------------------------------------------------------------------
!
      ELSE
!
!----------------------------------------------------------------------
        IF(NTSD==0)THEN
          QS=QLOW
        ENDIF
!
        ZU=Z0
        UZ0=0.
        VZ0=0.
!
        ZT=ZU*ZTFC
        THZ0=THS
!
        ZQ=ZT
        QZ0=QS
!----------------------------------------------------------------------
        TEM=(TLOW+TZ0)*0.5
        THM=(THELOW+THZ0)*0.5
!
        A=THM*P608
        B=(ELOCP/TEM-1.-P608)*THM
!
        DTHV=((THELOW-THZ0)*((QLOW+QZ0+CWMLOW)*(0.5*P608)+1.)          &
     &        +(QLOW-QZ0+CWMLOW)*A+CWMLOW*B)
!
        DU2=MAX((ULOW-UZ0)**2+(VLOW-VZ0)**2,EPSU2)
        RIB=BTGX*DTHV*ZSL/DU2
!----------------------------------------------------------------------
!       IF(RIB>=RIC)THEN
!         AKMS=MAX( VISC*RDZ,CXCHL)
!         AKHS=MAX(TVISC*RDZ,CXCHL)
!----------------------------------------------------------------------
!       ELSE  !  Turbulent branch
!----------------------------------------------------------------------
          ZSLU=ZSL+ZU
!
          RZSU=ZSLU/ZU
!
          RLOGU=LOG(RZSU)

          ZSLT=ZSL+ZU ! u,v and t are at the same level
!----------------------------------------------------------------------
!
!
!mp   Remove Topo modification of ZILFC term
!
!         TOPOTERM=TOPOFAC*ZSFC**2.
!         TOPOTERM=MAX(TOPOTERM,3.0)
!
!vcw
!  RIB modification to ZILFC term
!  7/29/2009 V Wong recommends 5, change pending
!
           CZETMAX = 10.
! stable
          IF(DTHV>0.)THEN
            IF (RIB<RIC) THEN
              ZZIL=ZILFC*(1.0+(RIB/RIC)*(RIB/RIC)*CZETMAX)
            ELSE
              ZZIL=ZILFC*(CZETMAX+1.0)
            ENDIF
! unstable
          ELSE
            ZZIL=ZILFC
          ENDIF

!----------------------------------------------------------------------
!
          land_point_iteration: DO ITR=1,ITRMX
!
!----------------------------------------------------------------------
!***  ZILITINKEVITCH FIX FOR ZT
!----------------------------------------------------------------------
!
            ZT=MAX(EXP(ZZIL*SQRT(USTAR*ZU))*ZU,EPSZT)
! new form  ZT=MAX(EXP(ZZIL*SQRT(USTAR*Z0BASE))*Z0BASE,EPSZT)
            RZST=ZSLT/ZT
            RLOGT=LOG(RZST)
!
!----------------------------------------------------------------------
!***  1./MONIN-OBUKHOV LENGTH-SCALE
!----------------------------------------------------------------------
!
            RLMO=ELFC*AKHS*DTHV/USTAR**3
            ZETALU=ZSLU*RLMO
            ZETALT=ZSLT*RLMO
            ZETAU=ZU*RLMO
            ZETAT=ZT*RLMO
!
            ZETALU=MIN(MAX(ZETALU,ZTMIN2),ZTMAX2)
            ZETALT=MIN(MAX(ZETALT,ZTMIN2),ZTMAX2)
            ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)
            ZETAT=MIN(MAX(ZETAT,ZTMIN2/RZST),ZTMAX2/RZST)
!
!----------------------------------------------------------------------
!***  LAND FUNCTIONS
!----------------------------------------------------------------------
!
            RZ=(ZETAU-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZ=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
            RZ=(ZETALU-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZL=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
            SIMM=PSMZL-PSMZ+RLOGU


!            if(abs(SIMM).lt.1.e-5.or.abs(SIMM).gt.1.e5)then
!              print*,'PSMZL,PSMZ,RLOGU=',PSMZL,PSMZ,RLOGU
!              print*,'RDZT,PSIM2(K+2)=',RDZT,PSIM2(K+2)
!              print*,'ZSL,ZU,Z0,ZSLU,ZSLT,ZT,RLMO=',ZSL,ZU,Z0,ZSLU,ZSLT,ZT,RLMO
!              print*,'RZ,K,ZETAU,ZTMIN2,DZETA2=',RZ,K,ZETAU,ZTMIN2,DZETA2
!              print*,'ELFC,AKHS,DTHV,USTAR=',ELFC,AKHS,DTHV,USTAR
!              stop
!            end if




!
            RZ=(ZETAT-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZ=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
            RZ=(ZETALT-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZL=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
            SIMH=(PSHZL-PSHZ+RLOGT)*FH02
!----------------------------------------------------------------------
            USTARK=USTAR*VKARMAN
            AKMS=MAX(USTARK/SIMM,CXCHL)
            AKHS=MAX(USTARK/SIMH,CXCHL)
!
!----------------------------------------------------------------------
!***  BELJAARS CORRECTION FOR USTAR
!----------------------------------------------------------------------
!
            IF(DTHV<=0.)THEN                                           !zj
              WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)                !zj
            ELSE                                                       !zj
              WSTAR2=0.                                                !zj
            ENDIF                                                      !zj
!
            USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
!
!----------------------------------------------------------------------
          ENDDO land_point_iteration
!----------------------------------------------------------------------
!
!       ENDIF  !  End of turbulant branch over land
!
!----------------------------------------------------------------------
!
      ENDIF  !  End of land/sea branch
!
!----------------------------------------------------------------------
!***  COUNTERGRADIENT FIX
!----------------------------------------------------------------------
!
!     HV=-AKHS*DTHV
!     IF(HV>0.)THEN
!       FCT=-10.*(BTGX)**(-1./3.)
!       CT=FCT*(HV/(PBLH*PBLH))**(2./3.)
!     ELSE
       CT=0.
!     ENDIF
!
!----------------------------------------------------------------------
      A1U=WGHT
      A1T=WGHTT
      A1Q=WGHTQ
      FFM=SIMM*(1.+WGHT)
      FFH=SIMH*(1.+WGHTT)
!      FFQ=SIMH*(1.+WGHTQ)
!----------------------------------------------------------------------
!***  THE FOLLOWING DIAGNOSTIC BLOCK PRODUCES 2-m and 10-m VALUES
!***  FOR TEMPERATURE, MOISTURE, AND WINDS.  IT IS DONE HERE SINCE
!***  THE VARIOUS QUANTITIES NEEDED FOR THE COMPUTATION ARE LOST
!***  UPON EXIT FROM THE ROTUINE.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      WSTAR=SQRT(WSTAR2)/WWST
!
      UMFLX=AKMS*(ULOW -UZ0 )
      VMFLX=AKMS*(VLOW -VZ0 )
      HSFLX=AKHS*(THLOW-THZ0)
      HLFLX=AKHS*(QLOW -QZ0 )
!----------------------------------------------------------------------
!     IF(RIB>=RIC)THEN
!----------------------------------------------------------------------
!       IF(SEAMASK>0.5)THEN
!         AKMS10=MAX( VISC/10.,CXCHS)
!         AKHS02=MAX(TVISC/02.,CXCHS)
!         AKHS10=MAX(TVISC/10.,CXCHS)
!       ELSE
!         AKMS10=MAX( VISC/10.,CXCHL)
!         AKHS02=MAX(TVISC/02.,CXCHL)
!         AKHS10=MAX(TVISC/10.,CXCHL)
!       ENDIF
!----------------------------------------------------------------------
!     ELSE
!----------------------------------------------------------------------
        ZU10=ZU+10.
        ZT02=ZT+02.
        ZT10=ZT+10.
!
        RLNU10=LOG(ZU10/ZU)
        RLNT02=LOG(ZT02/ZT)
        RLNT10=LOG(ZT10/ZT)
!
        ZTAU10=ZU10*RLMO
        ZTAT02=ZT02*RLMO
        ZTAT10=ZT10*RLMO
!
!----------------------------------------------------------------------
!***  SEA
!----------------------------------------------------------------------
!
        IF(SEAMASK>0.5)THEN
!
!----------------------------------------------------------------------
          ZTAU10=MIN(MAX(ZTAU10,ZTMIN1),ZTMAX1)
          ZTAT02=MIN(MAX(ZTAT02,ZTMIN1),ZTMAX1)
          ZTAT10=MIN(MAX(ZTAT10,ZTMIN1),ZTMAX1)
!----------------------------------------------------------------------
          RZ=(ZTAU10-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSM10=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
          SIMM10=PSM10-PSMZ+RLNU10
!
          RZ=(ZTAT02-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH02=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
          SIMH02=(PSH02-PSHZ+RLNT02)*FH01
!
          RZ=(ZTAT10-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH10=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
          SIMH10=(PSH10-PSHZ+RLNT10)*FH01
!
          AKMS10=MAX(USTARK/SIMM10,CXCHS)
          AKHS02=MAX(USTARK/SIMH02,CXCHS)
          AKHS10=MAX(USTARK/SIMH10,CXCHS)
!
!----------------------------------------------------------------------
!***  LAND
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
          ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
          ZTAT02=MIN(MAX(ZTAT02,ZTMIN2),ZTMAX2)
          ZTAT10=MIN(MAX(ZTAT10,ZTMIN2),ZTMAX2)
!----------------------------------------------------------------------
          RZ=(ZTAU10-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
          SIMM10=PSM10-PSMZ+RLNU10
!
          RZ=(ZTAT02-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH02=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
          SIMH02=(PSH02-PSHZ+RLNT02)*FH02
!
          RZ=(ZTAT10-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH10=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
          SIMH10=(PSH10-PSHZ+RLNT10)*FH02
!
          AKMS10=USTARK/SIMM10
          AKHS02=USTARK/SIMH02
          AKHS10=USTARK/SIMH10
!
          IF(AKMS10<=CXCHL) AKMS10=AKMS
          IF(AKHS02<=CXCHL) AKHS02=AKHS
          IF(AKHS10<=CXCHL) AKHS10=AKHS
!
!----------------------------------------------------------------------
        ENDIF
!----------------------------------------------------------------------
!     ENDIF
!----------------------------------------------------------------------
!
      U10 =UMFLX/AKMS10+UZ0
      V10 =VMFLX/AKMS10+VZ0
      TH02=HSFLX/AKHS02+THZ0
      TH10=HSFLX/AKHS10+THZ0
      Q02 =HLFLX/AKHS02+QZ0
      Q10 =HLFLX/AKHS10+QZ0
      TERM1=-0.068283/TLOW
      PSHLTR=PSFC*EXP(TERM1)
!
!----------------------------------------------------------------------
!***  COMPUTE "EQUIVALENT" Z0 TO APPROXIMATE LOCAL SHELTER READINGS.
!----------------------------------------------------------------------
!
      U10E=U10
      V10E=V10
!
      IF(SEAMASK<0.5)THEN

!1st        ZUUZ=MIN(0.5*ZU,0.1)
!1st        ZU=MAX(0.1*ZU,ZUUZ)
!tst        ZUUZ=amin1(ZU*0.50,0.3)
!tst        ZU=amax1(ZU*0.3,ZUUZ)

        ZUUZ=AMIN1(ZU*0.50,0.18)
        ZU=AMAX1(ZU*0.35,ZUUZ)
!
        ZU10=ZU+10.
        RZSU=ZU10/ZU
        RLNU10=LOG(RZSU)

        ZETAU=ZU*RLMO
        ZTAU10=ZU10*RLMO

        ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
        ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)

        RZ=(ZTAU10-ZTMIN2)/DZETA2
        K=INT(RZ)
        RDZT=RZ-REAL(K)
        K=MIN(K,KZTM2)
        K=MAX(K,0)
        PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
        SIMM10=PSM10-PSMZ+RLNU10
        EKMS10=MAX(USTARK/SIMM10,CXCHL)

        U10E=UMFLX/EKMS10+UZ0
        V10E=VMFLX/EKMS10+VZ0

      ENDIF
!
      U10=U10E
      V10=V10E
!
!----------------------------------------------------------------------
!***  SET OTHER WRF DRIVER ARRAYS
!----------------------------------------------------------------------
!
      RLOW=PLOW/(R_D*TLOW)
      CHS=AKHS
      CHS2=AKHS02
      CQS2=AKHS02
      HFX=-RLOW*CP*HSFLX
      QFX=-RLOW*HLFLX*WETM
      FLX_LH=XLV*QFX
      FLHC=RLOW*CP*AKHS
      FLQC=RLOW*AKHS*WETM
!!!   QGH=PQ0/PSHLTR*EXP(A2S*(TSK-A3S)/(TSK-A4S))
      QGH=((1.-SEAMASK)*PQ0+SEAMASK*PQ0SEA)                            &
     &     /PLOW*EXP(A2S*(TLOW-A3S)/(TLOW-A4S))
      CPM=CP*(1.+0.8*QLOW)
!
!***  DO NOT COMPUTE QS OVER LAND POINTS HERE SINCE IT IS
!***  A PROGNOSTIC VARIABLE THERE.  IT IS OKAY TO USE IT
!***  AS A DIAGNOSTIC OVER WATER SINCE IT WILL CAUSE NO
!***  INTERFERENCE BEFORE BEING RECOMPUTED IN MYJPBL.
!
      IF(SEAMASK>0.5)THEN
        QS=QLOW+QFX/(RLOW*AKHS)
        QS=QS/(1.-QS)
      ENDIF
!----------------------------------------------------------------------
      FM10=SIMM10+WGHT*SIMM
      FH2=SIMH02+WGHTT*SIMH
!
      END SUBROUTINE SFCDIF
!
!----------------------------------------------------------------------
      SUBROUTINE JSFC_INIT(USTAR                             &
     &                   ,RESTART                  &
     &                   ,IDS,IDE,JDS,JDE,KDS,KDE                      &
     &                   ,IMS,IME,JMS,JME,KMS,KME                      &
     &                   ,ITS,ITE,JTS,JTE,KTS,LM)
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
      LOGICAL,INTENT(IN) :: RESTART
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,LM
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: USTAR
!
!----------------------------------------------------------------------
!***  LOCAL VARIABLES
!----------------------------------------------------------------------
      INTEGER :: I,J,K,ITF,JTF,KTF
!
      REAL(kind=kfpt) :: X,ZETA1,ZETA2,ZRNG1,ZRNG2
!
      REAL(kind=kfpt) :: PIHF=3.1415926/2.,EPS=1.E-6
!
!----------------------------------------------------------------------
      JTF=MIN0(JTE,JDE-1)
      KTF=MIN0(LM,KDE-1)
      ITF=MIN0(ITE,IDE-1)
!
!
!***  FOR NOW, ASSUME SIGMA MODE FOR LOWEST MODEL LAYER
!
!----------------------------------------------------------------------
      IF(.NOT.RESTART)THEN
        DO J=JTS,JTE
        DO I=ITS,ITF
          USTAR(I,J)=0.1
        ENDDO
        ENDDO
      ENDIF

!----------------------------------------------------------------------
!
!***  COMPUTE SURFACE LAYER INTEGRAL FUNCTIONS
!
!----------------------------------------------------------------------
      FH01=1.
      FH02=1.
!
!      ZTMIN1=-10.0
!      ZTMAX1=2.0
!      ZTMIN2=-10.0
!      ZTMAX2=2.0
!org b
!      ZTMIN1=-5.0
!      ZTMAX1=1.0
!      ZTMIN2=-5.0
!      ZTMAX2=1.0
!org e
      ZTMIN1=-5.0
      ZTMAX1=9.0
      ZTMIN2=-5.0
      ZTMAX2=9.0

      ZRNG1=ZTMAX1-ZTMIN1
      ZRNG2=ZTMAX2-ZTMIN2
!
      DZETA1=ZRNG1/(KZTM-1)
      DZETA2=ZRNG2/(KZTM-1)
!
!----------------------------------------------------------------------
!***  FUNCTION DEFINITION LOOP
!----------------------------------------------------------------------
!
      ZETA1=ZTMIN1
      ZETA2=ZTMIN2
!
      DO K=1,KZTM
!
!----------------------------------------------------------------------
!***  UNSTABLE RANGE
!----------------------------------------------------------------------
!
        IF(ZETA1<0.)THEN
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
          X=SQRT(SQRT(1.-16.*ZETA1))
!
          PSIM1(K)=-2.*LOG((X+1.)/2.)-LOG((X*X+1.)/2.)+2.*ATAN(X)-PIHF
          PSIH1(K)=-2.*LOG((X*X+1.)/2.)
!
!----------------------------------------------------------------------
!***  STABLE RANGE
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
!
!         PSIM1(K)=5.*ZETA1
!         PSIH1(K)=5.*ZETA1
!----------------------------------------------------------------------
!***   HOLTSLAG AND DE BRUIN 1988
!----------------------------------------------------------------------
!
          PSIM1(K)=0.7*ZETA1+0.75*ZETA1*(6.-0.35*ZETA1)*EXP(-0.35*ZETA1)
          PSIH1(K)=0.7*ZETA1+0.75*ZETA1*(6.-0.35*ZETA1)*EXP(-0.35*ZETA1)
!----------------------------------------------------------------------
!
        ENDIF
!
!----------------------------------------------------------------------
!***  UNSTABLE RANGE
!----------------------------------------------------------------------
!
        IF(ZETA2<0.)THEN
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
!
          X=SQRT(SQRT(1.-16.*ZETA2))
!
          PSIM2(K)=-2.*LOG((X+1.)/2.)-LOG((X*X+1.)/2.)+2.*ATAN(X)-PIHF
          PSIH2(K)=-2.*LOG((X*X+1.)/2.)
!----------------------------------------------------------------------
!***  STABLE RANGE
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
!
!         PSIM2(K)=5.*ZETA2
!         PSIH2(K)=5.*ZETA2
!
!----------------------------------------------------------------------
!***  HOLTSLAG AND DE BRUIN 1988
!----------------------------------------------------------------------
!
          PSIM2(K)=0.7*ZETA2+0.75*ZETA2*(6.-0.35*ZETA2)*EXP(-0.35*ZETA2)
          PSIH2(K)=0.7*ZETA2+0.75*ZETA2*(6.-0.35*ZETA2)*EXP(-0.35*ZETA2)
!----------------------------------------------------------------------
!
        ENDIF
!
!----------------------------------------------------------------------
        IF(K==KZTM)THEN
          ZTMAX1=ZETA1
          ZTMAX2=ZETA2
        ENDIF
!
        ZETA1=ZETA1+DZETA1
        ZETA2=ZETA2+DZETA2
!----------------------------------------------------------------------
      ENDDO
!----------------------------------------------------------------------
      ZTMAX1=ZTMAX1-EPS
      ZTMAX2=ZTMAX2-EPS
!----------------------------------------------------------------------
!
      END SUBROUTINE JSFC_INIT
!
!----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_SF_JSFC
!
!-----------------------------------------------------------------------

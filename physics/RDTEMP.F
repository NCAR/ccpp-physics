!
!NCEP_MESO:MODEL_LAYER: PHYSICS
!
!***********************************************************************
      SUBROUTINE RDTEMP(NTSD,DT,JDAT,JULDAY,JULYR,IHRST,GLAT,GLON       &
     &                 ,CZEN,CZMEAN,T,RSWTT,RLWTT                       &
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                         &
     &                 ,IMS,IME,JMS,JME,KMS,KME                         &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    RDTEMP      RADIATIVE TEMPERATURE CHANGE
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 93-12-29
!     
! ABSTRACT:
!     RDTEMP APPLIES THE TEMPERATURE TENDENCIES DUE TO
!     RADIATION AT ALL LAYERS AT EACH ADJUSTMENT TIME STEP
!     
! PROGRAM HISTORY LOG:
!   87-09-??  BLACK      - ORIGINATOR
!   95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!   95-11-20  ABELES     - PARALLEL OPTIMIZATION
!   98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!   02-06-07  BLACK      - WRF CODING STANDARDS
!   02-09-09  WOLFE      - CONVERTING TO GLOBAL INDEXING
!   05-12-19  BLACK      - CONVERTED FROM IKJ TO IJK
!     
! USAGE: CALL RDTEMP FROM SUBROUTINE SOLVE_RUNSTREAM
!  
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!$$$  
!-----------------------------------------------------------------------
!mz      USE MODULE_MPP
!mz      USE MODULE_RA_GFDLETA,ONLY : CAL_MON_DAY,ZENITH
      USE module_radiation_astronomy ,ONLY : ZENITH
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: IHRST,JULDAY,JULYR,NTSD,JDAT(1:8)
!
      REAL,INTENT(IN) :: DT
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CZMEAN,GLAT,GLON    
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: RLWTT       &
     &                                                     ,RSWTT
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: T
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(in) :: CZEN
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,JDAY,JMONTH,K
      !MZ: Height boundary mask
      REAL,DIMENSION(IMS:IME,JMS:JME) :: HBM2
!
      INTEGER,DIMENSION(3) :: IDAT
!
      REAL :: DAYI,HOUR,TIMES,TTNDKL
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: CZEN2,XLAT2,XLON2
!
      REAL,DIMENSION(ITS:ITE,JTS:JTE) :: FACTR
!
      REAL :: DEGRAD=3.1415926/180.
      real :: xlat1,xlon1

      INTEGER :: MYIS, MYIE, MYJS, MYJE
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      MYIS=MAX(IDS,ITS)
      MYIE=MIN(IDE,ITE)
      MYJS=MAX(JDS,JTS)
      MYJE=MIN(JDE,JTE)
!-----------------------------------------------------------------------
!
!***  GET CURRENT VALUE OF COS(ZENITH ANGLE)
!
!      TIMES=NTSD*DT
!MZ      TIMES=XTIME*60.
!
!      DO J=MYJS,MYJE
!      DO I=MYIS,MYIE
!        XLAT2(I,J)=GLAT(I,J)
!        XLON2(I,J)=GLON(I,J)
!!!!!!!!!!!!Remove the following lines after bit-correct answers
!!!!!!!!!!!!are established with the control
!       xlat1=glat(i,j)/degrad
!       xlat2(i,j)=xlat1*degrad
!       xlon1=glon(i,j)/degrad
!       xlon2(i,j)=xlon1*degrad
!!!!!!!!!!!!
!!!!!!!!!!!!
!      ENDDO
!      ENDDO
!
!mz WRF
!      CALL CAL_MON_DAY(JULDAY,JULYR,JMONTH,JDAY)

!      IDAT(1)=JMONTH
!      IDAT(2)=JDAY
!      IDAT(3)=JULYR

!mz* using existing interstitials
!      IDAT(1)=jdat(2)
!      IDAT(2)=jdat(3)
!      IDAT(3)=jdat(1)
!
!mz   
      write(0,*)'rdtemp:max/min(czen) =', maxval(czen),minval(czen)
      write(0,*)'rdtemp:max/min(czmean) =',maxval(czmean),minval(czmean)
      write(0,*)'rdtemp:ntsd,dt,times,ihrst =',                  &
!     &       ntsd,dt, times,ihrst
!      CALL ZENITH(TIMES,DAYI,HOUR,IDAT,IHRST,XLON2,XLAT2,CZEN2          &
!     &           ,MYIS,MYIE,MYJS,MYJE                                   &
!     &           ,IDS,IDE,JDS,JDE,KDS,KDE                               &
!     &           ,IMS,IME,JMS,JME,KMS,KME                               &
!     &           ,ITS,ITE,JTS,JTE,KTS,KTE)
!mz
!      write(0,*)'rdtemp: max/min(czen2) =', maxval(czen2),minval(czen2)
!
      DO J=MYJS,MYJE
      DO I=MYIS,MYIE
!        CZEN(I,J)=CZEN2(I,J)
        IF(CZMEAN(I,J)>0.)THEN 
          FACTR(I,J)=CZEN(I,J)/CZMEAN(I,J)
        ELSE
          FACTR(I,J)=0.
        ENDIF
      ENDDO
      ENDDO
!
!mz*
      write(0,*)'RDTEMP: max/min(factr) = ',maxval(factr),minval(factr) 
      write(0,*)'RDTEMP: max/min(rswtt) = ',maxval(rswtt),minval(rswtt) 
      write(0,*)'RDTEMP: max/min(rlwtt) = ',maxval(rlwtt),minval(rlwtt)


      DO K=KTS,KTE
        DO J=MYJS,MYJE
        DO I=MYIS,MYIE
!MZ* see module_initialize_real.F
          HBM2(I,J) =1.
          TTNDKL=RSWTT(I,J,K)*FACTR(I,J)+RLWTT(I,J,K)
          T(I,J,K)=T(I,J,K)+TTNDKL*DT*HBM2(I,J)
        ENDDO
        ENDDO
      ENDDO
      write(0,*)'RDTEMP: max/min(T) = ', maxval(T), minval(T) 

!-----------------------------------------------------------------------
      END SUBROUTINE RDTEMP
!-----------------------------------------------------------------------

!>\file HWRF_radiation_post.F90
!! This file applies temperature tendency due to HWRF radiation.
      module HWRF_radiation_post
      contains

      subroutine HWRF_radiation_post_init ()
      end subroutine HWRF_radiation_post_init

      subroutine HWRF_radiation_post_finalize ()
      end subroutine HWRF_radiation_post_finalize


!> \section arg_table_HWRF_radiation_post_run Argument Table
!! \htmlinclude HWRF_radiation_post_run.html
!!
      subroutine HWRF_radiation_post_run (ncol, nlay,JDAT,              &
                            ntsd, dt,  julday, julyr,  ihrst,           &
                            glat, glon, czmean, rswtt, rlwtt,           &
                            t, czen, errmsg, errflg)
      
      USE MACHINE , only : kind_phys

      IMPLICIT NONE

      !-- interface variables
      INTEGER,               INTENT(IN) :: NCOL, NLAY, IHRST,JULDAY,    &
     &                                     JULYR,NTSD,JDAT(1:8)
      !MZ* dt-time step for physics?
      REAL(KIND_PHYS),       INTENT(IN) :: DT
      REAL(KIND_PHYS),DIMENSION(1:NCOL),INTENT(IN) :: CZMEAN,GLAT,GLON  

      REAL(KIND_PHYS),DIMENSION(1:NCOL,1:NLAY),INTENT(IN) :: RLWTT      &
     &                                                     ,RSWTT

      REAL(KIND_PHYS),DIMENSION(1:NCOL,1:NLAY),INTENT(INOUT) :: T

      REAL(KIND_PHYS),DIMENSION(1:NCOL),INTENT(in) :: CZEN

      !-- local variables
      INTEGER            :: I,K
      INTEGER            :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE

      REAL(KIND_PHYS)        :: XTIME
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0


      ! Set internal dimension
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
         kde = nlay
         kme = nlay
         kte = nlay

!----------------------------------------------------------------------
!***  APPLY TEMPERATURE TENDENCY DUE TO RADIATION
!----------------------------------------------------------------------
!
!mz* 
      xtime = ntsd*dt/60.

      CALL RDTEMP(ntsd,DT,JDAT,JULDAY,JULYR                             &
     &           ,XTIME,IHRST,glat,glon                                 &
     &           ,czen,czmean,t                                         &
     &           ,rswtt,rlwtt                                           &
     &           ,IDS,IDE,JDS,JDE,KDS,KDE                               &
     &           ,IMS,IME,JMS,JME,KMS,KME                               &
     &           ,ITS,ITE,JTS,JTE,KTS,KTE)
!

      end subroutine HWRF_radiation_post_run

      end module HWRF_radiation_post


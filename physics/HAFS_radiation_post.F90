!>\file HAFS_radiation_post.F90
!! This file applies temperature tendency due to HAFS radiation.
      module HAFS_radiation_post

      implicit none
     
      private
      
      public :: HAFS_radiation_post_init, HAFS_radiation_post_run,      &
                HAFS_radiation_post_finalize


      contains

      subroutine HAFS_radiation_post_init ()
      end subroutine HAFS_radiation_pre_init

      subroutine HAFS_radiation_post_finalize ()
      end subroutine HAFS_radiation_post_finalize


!> \section arg_table_HAFS_radiation_post_run Argument Table
!! \htmlinclude HAFS_radiation_post_run.html
!!
      subroutine HAFS_radiation_post_run (ncol, nlay, rswin, rswout,    &
                            ntsd, dt,  julday, julyr, xtime, ihrst,     &
                            glat, glon, czmean, rswtt, rlwtt, hbm2,     &
                            t, gsw, czen, errmsg, errflg)
      
      USE MACHINE , only : kind_phys

      IMPLICIT NONE

      !-- interface variables
      INTEGER,               INTENT(IN) :: IHRST,JULDAY,JULYR,NTSD
      REAL(KIND_PHYS),       INTENT(IN) :: DT,XTIME
      REAL(KIND_PHYS),DIMENSION(1:NCOL),INTENT(IN) :: CZMEAN,GLAT,GLON  &
     &                                             ,HBM2

      REAL(KIND_PHYS),DIMENSION(1:NCOL,1:NLAY),INTENT(IN) :: RLWTT       &
     &                                                     ,RSWTT

      REAL(KIND_PHYS),DIMENSION(1:NCOL,1:NLAY),INTENT(INOUT) :: T

      REAL(KIND_PHYS),DIMENSION(1:NCOL),INTENT(OUT) :: CZEN

      !-- local variables
      INTEGER            :: I,K
      INTEGER,           :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE

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

 

        DO J=jts,jte
        DO I=its,ite
          gsw(I,J)=rswin(I,J)-rswout(I,J)
        ENDDO
        ENDDO


!----------------------------------------------------------------------
!***  APPLY TEMPERATURE TENDENCY DUE TO RADIATION
!----------------------------------------------------------------------
!

      CALL RDTEMP(ntsd,DT,JULDAY,JULYR                                  &
     &           ,XTIME,IHRST,glat,glon                                 &
     &           ,czen,czmean,t                                         &
     &           ,rswtt,rlwtt,hbm2                                      &
     &           ,IDS,IDF,JDS,JDF,KDS,KDE                               &
     &           ,IMS,IME,JMS,JME,KMS,KME                               &
     &           ,ITS,ITE,JTS,JTE,KTS,KTE)
!

      end subroutine HAFS_radiation_post_run
      end module HAFS_radiation_post


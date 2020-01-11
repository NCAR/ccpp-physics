!>\file HWRF_radiation_pre.F90
!! This file contains CCPP-compliant ETAMP_TO_MOIST() in HWRF.  
!!
!> This module contains the CCPP-compliant ETAMP_TO_MOIST in HWRF .
     module HWRF_radiation_pre

     contains

     subroutine HWRF_radiation_pre_init ()
     end subroutine HWRF_radiation_pre_init

     subroutine HWRF_radiation_pre_finalize ()
     end subroutine HWRF_radiation_pre_finalize

!> \defgroup hwrf_etamp_to_moist HWRF Update Moist Module
!! This subroutine is to update water array with CWM, F_RAIN, and F_ICE. 
!! \section arg_table_HWRF_radiation_pre_run Argument Table
!! \htmlinclude HWRF_radiation_pre_run.html
!!
     subroutine HWRF_radiation_pre_run (CWM,QV, F_ICE,F_RAIN            &
                                        ,qc,qi,qr                       &
                                        ,qv_r, qc_r,qr_r,qi_r           & 
                                        ,qs_r, qg_r                     &      
                                        ,LM,IME,errmsg,errflg )

       USE MACHINE , only : kind_phys
       IMPLICIT NONE
     
!----------------------
!-- Argument Variables
!----------------------
!
      INTEGER,INTENT(IN) :: LM,IME
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(IN) ::  CWM,QV, &
                                                                QC,QI,  &
                                                                QR
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(IN) ::          &
                                                               F_ICE,   &
                                                               F_RAIN 

      !dry mixing ratio used in HWRF RRTMG/TURBL/FER
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(OUT) :: qv_r,   &
                                                                qc_r,   &
                                                                qi_r,   &
                                                                qr_r,   &
                                                                qs_r,   &
                                                                qg_r

!
!--------------------
!--  Local Variables
!--------------------
!
      INTEGER :: I,K
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

            if(size(f_ice,1)*size(f_ice,2) <= 1) then
              return
            endif

            DO K=1,LM
              DO I=1,IME
                qv_r (i,k)=qv(i,k)/(1.-qv(i,k))
                qi_r(I,K) =0.
                qr_r(I,K) =0.
                qc_r(I,K) =0.
                IF(F_ICE(I,K)>=1.) THEN
                   qi_r(I,K) = CWM(i,k)
                ELSEIF(F_ICE(I,K)<=0.) THEN
                   qc_r(I,K) = CWM(I,K)
                ELSE
                   qi_r(I,K) = F_ICE(I,K)*CWM(I,K)
                   qc_r(I,K) = CWM(I,K)-qi_r(I,K)
                ENDIF

                IF(qc_r(I,K)>0. .AND. F_RAIN(I,K)>0.) THEN
                   IF(F_RAIN(I,K)>=1.)THEN
                      qr_r(I,K)=qc_r(I,K)
                      qc_r(I,K)=0.
                   ELSE
                      qr_r(I,K)=F_RAIN(I,K)*qc_r(I,K)
                      qc_r(I,K)=qc_r(I,K)-qr_r(I,K)
                   ENDIF
                ENDIF
                 
                !no snow and graupel mixing ratio in HWRF physics
                qg_r (I,K) = 0.  
                qs_r (I,K) = 0.


              ENDDO
            ENDDO

      end subroutine HWRF_radiation_pre_run

      end module HWRF_radiation_pre

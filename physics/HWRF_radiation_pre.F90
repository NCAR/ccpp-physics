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
      subroutine HWRF_radiation_pre_run (qv,qc,qi,qr,qs,qg              &
                                        ,qv_r, qc_r,qr_r,qi_r           &
                                        ,qs_r, qg_r                     &
                                        ,imp_physics                    &
                                        ,imp_physics_fer_hires          &
                                        ,KDT,LM,IME,mpirank, mpiroot    &
                                        ,errmsg,errflg )

       USE MACHINE , only : kind_phys
       IMPLICIT NONE
     
!----------------------
!-- Argument Variables
!----------------------
!
      INTEGER,INTENT(IN) :: LM,IME
      INTEGER,INTENT(IN) :: KDT, imp_physics, imp_physics_fer_hires
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(IN) ::  QV,     &
                                                                QC,QI,  &
                                                                QR,QS,  &
                                                                QG

      !dry mixing ratio used in HWRF RRTMG/TURBL/FER
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(OUT) :: qv_r,   &
                                                                qc_r,   &
                                                                qi_r,   &
                                                                qr_r,   &
                                                                qs_r,   &
                                                                qg_r
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      integer,                   intent(in)    :: mpirank
      integer,                   intent(in)    :: mpiroot


!
!--------------------
!--  Local Variables
!--------------------
!
      INTEGER :: I,K
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0


      if (mpirank == mpiroot )  then
         write(0,*)'------- kdt = ', kdt
         write(0,*)'max/min(qv) = ',maxval(qv),minval(qv)
         write(0,*)'max/min(qc) = ',maxval(qc),minval(qc)
         write(0,*)'max/min(qi) = ',maxval(qi),minval(qi)
         write(0,*)'max/min(qr) = ',maxval(qr),minval(qr)
      endif

      DO K=1,LM
        DO I=1,IME
           qv_r(i,k) = qv(i,k)/(1.0_kind_phys-qv(i,k))
           qc_r(i,k) = qc(i,k)/(1.0_kind_phys-qv(i,k))
           qi_r(i,k) = qi(i,k)/(1.0_kind_phys-qv(i,k))
           qr_r(i,k) = qr(i,k)/(1.0_kind_phys-qv(i,k))
           if (imp_physics == imp_physics_fer_hires) then
              qg_r (I,K) = 0.  
              qs_r (I,K) = 0.
           else
              qg_r(i,k) = qg(i,k)/(1.0_kind_phys-qv(i,k))
              qs_r(i,k) = qs(i,k)/(1.0_kind_phys-qv(i,k))
           endif
        ENDDO
      ENDDO

      if (mpirank == mpiroot )  then
         write(0,*)'------------'
         write(0,*)'max/min(qv_r) = ',maxval(qv_r),minval(qv_r)
         write(0,*)'max/min(qc_r) = ',maxval(qc_r),minval(qc_r)
         write(0,*)'max/min(qi_r) = ',maxval(qi_r),minval(qi_r)
         write(0,*)'max/min(qr_r) = ',maxval(qr_r),minval(qr_r)
         write(0,*)'max/min(qg_r) = ',maxval(qg_r),minval(qg_r) 
         write(0,*)'max/min(qs_r) = ',maxval(qs_r),minval(qs_r)
      endif


      end subroutine HWRF_radiation_pre_run

      end module HWRF_radiation_pre

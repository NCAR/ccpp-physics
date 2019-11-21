!>\file HAFS_update_moist.F90
!! This file contains CCPP-compliant UPDATE_MOIST() in HWRF.  
!!
!! In HWRF , this subroutine should be called before:
!! - radiation()
!! - CUCNVC()
!! - TURBLE()

!> This module contains the CCPP-compliant UPDATE_MOIST for Ferrier-Aligo MP.
     module HAFS_update_moist

     implicit none

     private

     public :: HAFS_update_moist_init, HAFS_update_moist_run,           &
               HAFS_update_moist_finalize

     contains

     subroutine HAFS_update_moist_init ()
     end subroutine HAFS_update_moist_init

     subroutine HAFS_update_moist_finalize ()
     end subroutine HAFS_update_moist_finalize

!> \defgroup hafs_update_moist HAFS Update Moist Module
!! This subroutine is to update water array with CWM, F_RAIN, and F_ICE 
!! and convert moist mixing ratio to dry mixing ratio
!! \section arg_table_HAFS_update_moist_run Argument Table
!! \htmlinclude HAFS_update_moist_run.html
!!
     subroutine HAFS_update_moist_run (CWM,QV, F_ICE,F_RAIN             &
                                        ,qc,qi,qr                       &
                                        ,imp_physics                    &
                                        ,imp_physics_fer_hires          &
                                        ,qv_r, qc_r,qr_r,qi_r           & !-output: dry mixing ratioes for HWRF physics
                                        ,qs_r, qg_r                     &      
                                        ,spec_adv                       &
                                        ,LM,IME,errmsg,errflg )

       USE MACHINE , only : kind_phys
       IMPLICIT NONE
     
!----------------------
!-- Argument Variables
!----------------------
!
       INTEGER,INTENT(IN) :: LM,IME
!

      LOGICAL,INTENT(IN) :: SPEC_ADV
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(IN) ::  CWM,QV, &
                                                                QC,QI,  &
                                                                QR
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(IN) ::          &
                                                               F_ICE,   &
                                                               F_RAIN 
      integer,                        intent(in)    :: imp_physics
      integer,                        intent(in)    :: imp_physics_fer_hires

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

!-- Update WATER arrays when advecting only total condensate (spec_adv=F)
!-- and F_*   or at the initial time step
      if (imp_physics == imp_physics_fer_hires) then
        if (spec_adv) then
            qv_r (i,k)=qv(i,k)/(1.-qv(i,k))
            qc_r (i,k)=qc(i,k)/(1.-qv(i,k))
            qi_r (i,k)=qi(i,k)/(1.-qv(i,k))
            qr_r (i,k)=qr(i,k)/(1.-qv(i,k))
            qs_r (i,k)= 0.
            qg_r (i,k)= 0.

        else   ! .not.spec_adv
            DO K=1,LM
              DO I=1,IME
                !!MZ: HWRF::UPDATE_MOIST() solution
                !calculate dry mixing ratio of all q
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
                qs_r (I,K) = 0.
                qg_r (I,K) = 0.
              ENDDO
            ENDDO
        end if
      else    !.not. fer_hires
        write(errmsg,'(*(a))') "Logic error: HAFS_update_moist only works for HWRF physics"
        errflg = 1
        return
      end if

      end subroutine HAFS_update_moist_run

  end module HAFS_update_moist

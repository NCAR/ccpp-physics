!>\file mp_FA_update_moist.F90
!! This file contains CCPP-compliant UPDATE_MOIST() in HWRF.  
!!
!! In HWRF (.not.spec_adv), this subroutine should be called before:
!! - radiation()
!! - CUCNVC()
!! - TURBLE()

!> This module contains the CCPP-compliant UPDATE_MOIST for Ferrier-Aligo MP.
     module mp_FA_update_moist

     implicit none

     private

     public :: mp_FA_update_moist_init, mp_FA_update_moist_run, &
               mp_FA_update_moist_finalize

     contains

     subroutine mp_FA_update_moist_init ()
     end subroutine mp_FA_update_moist_init

     subroutine mp_FA_update_moist_finalize ()
     end subroutine mp_FA_update_moist_finalize

!> \defgroup hafs_fa_update HAFS Ferrier-Aligo MP Scheme Update Water Module
!! \ingroup hafs_famp
!! This subroutine is to update water array with CWM, F_RAIN, and F_ICE 
!! for Ferrier-Aligo MP scheme. 
!! \section arg_table_mp_FA_update_moist_run Argument Table
!! | local_name     | standard_name                                         | long_name                                                                                  | units   | rank | type      | kind      | intent | optional |
!! |----------------|-------------------------------------------------------|--------------------------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | cwm            | total_cloud_condensate_mixing_ratio_updated_by_physics| total cloud condensate mixing ratio (except water vapor) updated by physics                | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | f_ice          | fraction_of_ice_water_cloud                           | fraction of ice water cloud                                                                | frac    |    2 | real      | kind_phys | in     | F        |
!! | f_rain         | fraction_of_rain_water_cloud                          | fraction of rain water cloud                                                               | frac    |    2 | real      | kind_phys | in     | F        |
!! | qc             | cloud_condensed_water_mixing_ratio_updated_by_physics | moist (dry+vapor, no condensates) mixing ratio of cloud condensed water updated by physics | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | qr             | rain_water_mixing_ratio_updated_by_physics            | moist (dry+vapor, no condensates) mixing ratio of rain water updated by physics            | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | qi             | ice_water_mixing_ratio_updated_by_physics             | moist (dry+vapor, no condensates) mixing ratio of ice water updated by physics             | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | qs             | snow_water_mixing_ratio_updated_by_physics            | moist (dry+vapor, no condensates) mixing ratio of snow water updated by physics            | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | lm             | vertical_dimension                                    | number of vertical levels                                                                  | count   |    0 | integer   |           | in     | F        |
!! | ime            | horizontal_dimension                                  | horizontal dimension                                                                       | count   |    0 | integer   |           | in     | F        |
!! | errmsg         | ccpp_error_message                                    | error message for error handling in CCPP                                                   | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                       | error flag for error handling in CCPP                                                      | flag    |    0 | integer   |           | out    | F        |
!!
     subroutine mp_FA_update_moist_run (CWM,F_ICE,F_RAIN,               &
                                        ,QC,QR,QI,QS,                   &
                                        ,LM,IME,errmsg,errflg )

       USE MACHINE , only : kind_phys
       IMPLICIT NONE
     
!----------------------
!-- Argument Variables
!----------------------
!
       INTEGER,INTENT(IN) :: LM,IME
!
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(IN) :: CWM,     &
                                                               F_ICE,   &
                                                               F_RAIN 
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(OUT) :: QC,QR,  &
                                                           ,QI,QS    ! qs is total ice in F-A
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

!MZ
        IF (SPEC_ADV) then
          DO K=1,LM
            DO I=1,IME
             !MZ: need to check
             CWM(I,K)=QC(I,K)+QR(I,K)+QI(I,K)+QS(I,K)
            ENDDO
          ENDDO
        ELSE

!-- Update WATER arrays when advecting only total condensate (spec_adv=F)
!-- and F_*   or at the initial time step
            DO K=1,LM
              DO I=1,IME
                !!MZ: HWRF::UPDATE_MOIST() solution
                !mixing ratio of qv
                !MOIST (i,k)=q(i,k)/(1.-q(i,k))
                QI(I,K) =0.
                QR(I,K) =0.
                QC(I,K) =0.
                IF(F_ICE(I,K)>=1.) THEN
                   QI(I,K) = CWM(i,k)
                ELSEIF(F_ICE(I,K)<=0.) THEN
                   QC(I,K) = CWM(I,K)
                ELSE
                   QI(I,K) = F_ICE(I,K)*CWM(I,K)
                   QC(I,K) = CWM(I,K)-QI(I,K)
                ENDIF

                IF(QC(I,K)>0. .AND. F_RAIN(I,K)>0.) THEN
                   IF(F_RAIN(I,K)>=1.)THEN
                      QR(I,K)=QC(I,K)
                      QC(I,K)=0.
                   ELSE
                      QR(I,K)=F_RAIN(I,K)*QC(I,K)
                      QC(I,K)=QC(I,K)-QR(I,K)
                   ENDIF
                ENDIF
                qi(i,k) = 0.
                qs(i,k) = qi(i,k)    !MZ: qs contains total ice in HWRF
              ENDDO
            ENDDO
        ENDIF
      end subroutine mp_FA_update_moist_run

  end module mp_FA_update_moist

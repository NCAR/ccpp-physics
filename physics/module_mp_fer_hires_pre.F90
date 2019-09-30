!>\file module_mp_fer_hires_pre.F90
!! This file contains a wrap around UPDATE_WATER in NAM/module_MICROPHYSICS.F90
!! which should be called the first in the physics suite using Ferrier-Aligo MP
!! scheme.

!> This module contains the CCPP-compliant UPDATE_WATER for Ferrier-Aligo MP.
     module mp_fer_hires_pre

     implicit none

     private

     public :: mp_fer_hires_pre_init, mp_fer_hires_pre_run, &
               mp_fer_hires_pre_finalize

     contains

     subroutine mp_fer_hires_pre_init ()
     end subroutine mp_fer_hires_pre_init

     subroutine mp_fer_hires_pre_finalize ()
     end subroutine mp_fer_hires_pre_finalize

!> \defgroup hafs_fa_update HAFS Ferrier-Aligo MP Scheme Update Water Module
!! \ingroup hafs_famp
!! This subroutine is to update water array with CWM, F_RAIN, and F_ICE 
!! for Ferrier-Aligo MP scheme. 
!! \section arg_table_mp_fer_hires_pre_run Argument Table
!! | local_name     | standard_name                                         | long_name                                                                                  | units   | rank | type      | kind      | intent | optional |
!! |----------------|-------------------------------------------------------|--------------------------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | cwm            | total_cloud_condensate_mixing_ratio_updated_by_physics| total cloud condensate mixing ratio (except water vapor) updated by physics                | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | f_ice          | fraction_of_ice_water_cloud                           | fraction of ice water cloud                                                                | frac    |    2 | real      | kind_phys | inout  | F        |
!! | f_rain         | fraction_of_rain_water_cloud                          | fraction of rain water cloud                                                               | frac    |    2 | real      | kind_phys | inout  | F        |
!! | f_rimef        | rime_factor                                           | rime factor                                                                                | frac    |    2 | real      | kind_phys | inout  | F        |
!! | epsq           | minimum_value_of_specific_humidity                    | floor value for specific humidity                                                          | kg kg-1 |    0 | real      | kind_phys | in     | F        |
!! | t              | air_temperature                                       | model layer mean temperature                                                               | K       |    2 | real      | kind_phys | inout  | F        |
!! | qc             | cloud_condensed_water_mixing_ratio_updated_by_physics | moist (dry+vapor, no condensates) mixing ratio of cloud condensed water updated by physics | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | qr             | rain_water_mixing_ratio_updated_by_physics            | moist (dry+vapor, no condensates) mixing ratio of rain water updated by physics            | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | qi             | ice_water_mixing_ratio_updated_by_physics             | moist (dry+vapor, no condensates) mixing ratio of ice water updated by physics             | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | qg             | mass_weighted_rime_factor_updated_by_physics          | mass weighted rime factor updaed by physics                                                | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | spec_adv       | flag_for_individual_cloud_species_advected            | flag for individual cloud species advected                                                 | flag    |    0 | logical   |           | in     | F        |
!! | kdt            | index_of_time_step                                    | current forecast interation                                                                | index   |    0 | integer   |           | in     | F        |
!! | lm             | vertical_dimension                                    | number of vertical levels                                                                  | count   |    0 | integer   |           | in     | F        |
!! | ime            | horizontal_dimension                                  | horizontal dimension                                                                       | count   |    0 | integer   |           | in     | F        |
!! | errmsg         | ccpp_error_message                                    | error message for error handling in CCPP                                                   | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                       | error flag for error handling in CCPP                                                      | flag    |    0 | integer   |           | out    | F        |
!!
     subroutine mp_fer_hires_pre_run (CWM,F_ICE,F_RAIN,F_RIMEF          &
                             ,EPSQ                                      &
                             ,T,QC,QR,QI,QG                             &
                             ,SPEC_ADV,kdt                              &
                             ,LM,IME,errmsg,errflg )

       USE MACHINE , only : kind_phys
       IMPLICIT NONE
     
       REAL(KIND=KIND_PHYS) :: EPSQ

!----------------------
!-- Argument Variables
!----------------------
!
       INTEGER,INTENT(IN) :: KDT,LM,IME
!
      LOGICAL,INTENT(IN) :: SPEC_ADV
!
      REAL(kind=kind_phys),DIMENSION(1:IME,1:LM),INTENT(INOUT) :: CWM   &
                                                           ,F_ICE       &
                                                           ,F_RAIN      &
                                                           ,F_RIMEF     &
                                                           ,T,QC,QR     &
                                                           ,QI,QG        ! QS
!
!--------------------
!--  Local Variables
!--------------------
!
      INTEGER :: I,J,K
      REAL(kind=kind_phys) :: FRACTION, LIQW, OLDCWM
      LOGICAL :: CLD_INIT
      LOGICAL :: deep_ice
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!
!zm      IF(NTIMESTEP<=1)THEN
      IF(kdt<=2)THEN
        CLD_INIT=.TRUE.
      ELSE
        CLD_INIT=.FALSE.
      ENDIF


      IF (.NOT.SPEC_ADV .OR. CLD_INIT) THEN
!-- Update WATER arrays when advecting only total condensate (spec_adv=F)
!-- and F_*   or at the initial time step
            DO K=1,LM
              DO I=1,IME
                !!CWM(I,K)=QC(I,K)+QR(I,K)+QI(I,K)
                IF (CWM(I,K)>EPSQ) THEN
                  LIQW=(1.-F_ice(I,K))*CWM(I,K)
                  QC(I,K)=(1.-F_rain(I,K))*LIQW
                  QR(I,K)=F_rain(I,K)*LIQW
                  QI(I,K)=F_ice(I,K)*CWM(I,K)
                  QG(I,K)=F_RIMEF(I,K)*QI(I,K)
                ELSE
                  QC(I,K)=0.
                  QR(I,K)=0.
                  QI(I,K)=0.
                  QG(I,K)=0.
                ENDIF
              ENDDO
            ENDDO

      ELSE 
!-- Update CWM,F_ICE,F_RAIN arrays from separate species advection (spec_adv=T)
            DO K=1,LM
              DO I=1,IME
                CWM(I,K)=QC(I,K)+QR(I,K)+QI(I,K)
                IF (QI(I,K)>EPSQ) THEN
                  !F_ICE(I,K)=QI(I,K)/CWM(I,K)
                  F_ICE(I,K) = MAX(0., MIN(1., QI(I,K)/CWM(I,K)))
                  !MZ: f_rimef= qrimef/qi
                  F_RIMEF(I,K) = QG(I,K)/QI(I,K)
                ELSE
                  F_ICE(I,K)=0.0
                  F_RIMEF(I,K) = 1.0
                ENDIF
                IF (QR(I,K)>EPSQ) THEN
                  F_RAIN(I,K)=QR(I,K)/(QC(I,K)+QR(I,K))
                ELSE
                  F_RAIN(I,K)=0.
                ENDIF
              ENDDO
            ENDDO
      ENDIF 

      end subroutine mp_fer_hires_pre_run

  end module mp_fer_hires_pre

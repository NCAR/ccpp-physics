!>\file module_soil_pre.F90
!! This file contains subroutines that initialize RUC LSM levels, soil
!! temperature/moisture.
module module_soil_pre

!tgs Initialize RUC LSM levels, soil temp/moisture

      implicit none

      private

      public init_soil_depth_3, init_soil_3_real

contains

!>\ingroup lsm_ruc_group
!> This subroutine defines level depth in soil and thickness of soil 
!! layers RUC LSM.
!!
!! In RUC LSM ZS-soil levels, and DZS-soil layer thicknesses, not used
!! ZS is specified in the namelist: num_soil_levels =6 or 9.
!! other options with number of levels are possible, but 
!! WRF users should change consistently the namelist entry with the 
!! ZS array in this subroutine.
   SUBROUTINE init_soil_depth_3 ( zs , dzs , num_soil_levels )

      INTEGER, INTENT(IN) :: num_soil_levels

      REAL, DIMENSION(1:num_soil_levels), INTENT(OUT)  ::  zs, dzs
      REAL, DIMENSION(1:num_soil_levels)  :: zs2

      INTEGER                   ::      l

      CHARACTER (LEN=132) :: message

! in RUC LSM ZS - soil levels, and DZS - soil layer thicknesses, not used
! ZS is specified in the namelist: num_soil_levels = 6 or 9.
! Other options with number of levels are possible, but
! WRF users should change consistently the namelist entry with the
!    ZS array in this subroutine.

     IF ( num_soil_levels .EQ. 6) THEN
      zs  = (/ 0.00 , 0.05 , 0.20 , 0.40 , 1.60 , 3.00 /)
     ELSEIF ( num_soil_levels .EQ. 9) THEN
      zs  = (/ 0.00 , 0.01 , 0.04 , 0.10 , 0.30, 0.60, 1.00 , 1.60, 3.00 /)
      !zs  = (/ 0.00 , 0.05 , 0.20 , 0.40 , 0.60, 1.00, 1.60 , 2.20, 3.00 /)
     ENDIF

      zs2(1) = 0.
      zs2(2) = (zs(2) + zs(1))*0.5
      dzs(1) = zs2(2) - zs2(1)
     do l = 2, num_soil_levels - 1
      zs2(l) = (zs(l+1) + zs(l)) * 0.5
      dzs(l) = zs2(l) - zs2(l-1)
     enddo
      zs2(num_soil_levels) = zs(num_soil_levels)
      dzs(num_soil_levels) = zs2(num_soil_levels) - zs2(num_soil_levels-1)

      IF ( num_soil_levels .EQ. 4 .OR. num_soil_levels .EQ. 5 ) THEN
         WRITE(message,FMT= '(A)')'Usually, the RUC LSM uses 6, 9 or more levels.  Change this in the namelist.'
!         CALL wrf_error_fatal ( message )
      END IF

   END SUBROUTINE init_soil_depth_3

!>\ingroup lsm_ruc_group
!! This subroutine initializes soil moisture and temperature at RUC vertical levels
!! from the Noah layers. RUC has 3 levels in the top Noah layer, therefore, initialization of 
!! soil moisture at these top levels is questionable.
   SUBROUTINE init_soil_3_real ( tsk , tmn , smois , tslb , &
                                 st_input , sm_input , landmask, sst, &
                                 zs , dzs , &
                                 st_levels_input , sm_levels_input , &
                                 num_soil_layers , num_st_levels_input , num_sm_levels_input ,  &
                                 num_st_levels_alloc , num_sm_levels_alloc , &
                                 flag_sst , flag_soil_layers , flag_soil_levels , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )

      INTEGER , INTENT(IN) :: num_soil_layers , &
                              num_st_levels_input , num_sm_levels_input , &
                              num_st_levels_alloc , num_sm_levels_alloc , &
                              ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte

      INTEGER , INTENT(IN) :: flag_sst, flag_soil_layers, flag_soil_levels

      INTEGER , DIMENSION(1:num_st_levels_input) , INTENT(INOUT) :: st_levels_input
      INTEGER , DIMENSION(1:num_sm_levels_input) , INTENT(INOUT) :: sm_levels_input

      REAL , DIMENSION(ims:ime,1:num_st_levels_alloc,jms:jme) , INTENT(INOUT) :: st_input
      REAL , DIMENSION(ims:ime,1:num_sm_levels_alloc,jms:jme) , INTENT(INOUT) :: sm_input
      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(IN) :: landmask , sst

      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(IN) :: tmn
      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(INOUT) :: tsk
      REAL , DIMENSION(num_soil_layers) :: zs , dzs

      REAL , DIMENSION(ims:ime,num_soil_layers,jms:jme) , INTENT(OUT) :: tslb , smois

      REAL , ALLOCATABLE , DIMENSION(:) :: zhave

      logical :: debug_print = .false.
      INTEGER :: i , j , l , lout , lin , lwant , lhave, k
      REAL :: temp

      !  Allocate the soil layer array used for interpolating.      

      IF ( ( num_st_levels_input .LE. 0 ) .OR. &
           ( num_sm_levels_input .LE. 0 ) ) THEN
         write (0, FMT='(A)')&
'No input soil level data (either temperature or moisture, or both are missing).  Required for RUC LSM.'
      ELSE
         IF ( flag_soil_levels == 1 ) THEN
           if (debug_print) write(0, FMT='(A)') ' Assume RUC LSM input'
           ALLOCATE ( zhave( MAX(num_st_levels_input,num_sm_levels_input)  ) )
         ELSE
           if (debug_print) write(0, FMT='(A)') ' Assume non-RUC LSM input'
           ALLOCATE ( zhave( MAX(num_st_levels_input,num_soil_layers)  ) )
         END IF
      END IF
      
      !  Sort the levels for temperature.

      outert : DO lout = 1 , num_st_levels_input-1
         innert : DO lin = lout+1 , num_st_levels_input
            IF ( st_levels_input(lout) .GT. st_levels_input(lin) ) THEN
               temp = st_levels_input(lout)
               st_levels_input(lout) = st_levels_input(lin)
               st_levels_input(lin) = NINT(temp)
               DO j = jts , jte
                  DO i = its ,ite
                     temp = st_input(i,lout,j)
                     st_input(i,lout,j) = st_input(i,lin,j)
                     st_input(i,lin,j) = temp
                  END DO
               END DO
            END IF
         END DO innert
      END DO outert

      IF ( flag_soil_layers == 1 ) THEN
      DO j = jts , jte
         DO i = its , ite
            st_input(i,1,j) = tsk(i,j)
            st_input(i,num_st_levels_input+2,j) = tmn(i,j)
         END DO
      END DO
      END IF

      !  Sort the levels for moisture.      

      outerm: DO lout = 1 , num_sm_levels_input-1
         innerm : DO lin = lout+1 , num_sm_levels_input
            IF ( sm_levels_input(lout) .GT. sm_levels_input(lin) ) THEN
               temp = sm_levels_input(lout)
               sm_levels_input(lout) = sm_levels_input(lin)
               sm_levels_input(lin) = NINT(temp)
               DO j = jts ,jte
                  DO i = its , ite
                     temp = sm_input(i,lout,j)
                     sm_input(i,lout,j) = sm_input(i,lin,j)
                     sm_input(i,lin,j) = temp
                  END DO
               END DO
            END IF
         END DO innerm
      END DO outerm

      IF ( flag_soil_layers == 1 ) THEN
      DO j = jts , jte
         DO i = its , ite
            sm_input(i,1,j) = (sm_input(i,2,j)-sm_input(i,3,j))/   &
                              (st_levels_input(2)-st_levels_input(1))*st_levels_input(1)+  &
                              sm_input(i,2,j)

            sm_input(i,num_sm_levels_input+2,j) = sm_input(i,num_sm_levels_input+1,j)
         END DO
      END DO
      END IF

      !  Here are the levels that we have from the input for temperature.      

      IF ( flag_soil_levels == 1 ) THEN
         DO l = 1 , num_st_levels_input
            zhave(l) = st_levels_input(l) / 100.
         END DO

      
      !  Interpolate between the layers we have (zhave) and those that we want
      !  (zs).

      z_wantt : DO lwant = 1 , num_soil_layers
         z_havet : DO lhave = 1 , num_st_levels_input -1
            IF ( ( zs(lwant) .GE. zhave(lhave  ) ) .AND. &
                 ( zs(lwant) .LE. zhave(lhave+1) ) ) THEN
               DO j = jts , jte
                  DO i = its , ite
                     tslb(i,lwant,j)= ( st_input(i,lhave,j ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                        st_input(i,lhave+1,j) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                ( zhave(lhave+1) - zhave(lhave) )
                  END DO
               END DO
               EXIT z_havet
            END IF
         END DO z_havet
      END DO z_wantt

      ELSE

         zhave(1) = 0.
         DO l = 1 , num_st_levels_input
            zhave(l+1) = st_levels_input(l) / 100.
         END DO
         zhave(num_st_levels_input+2) = 300. / 100.

      !  Interpolate between the layers we have (zhave) and those that we want
      !  (zs).      

      z_wantt_2 : DO lwant = 1 , num_soil_layers
         z_havet_2 : DO lhave = 1 , num_st_levels_input +2
            IF ( ( zs(lwant) .GE. zhave(lhave  ) ) .AND. &
                 ( zs(lwant) .LE. zhave(lhave+1) ) ) THEN
               DO j = jts , jte
                  DO i = its , ite
                     tslb(i,lwant,j)= ( st_input(i,lhave,j ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                        st_input(i,lhave+1,j) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                ( zhave(lhave+1) - zhave(lhave) )
                  END DO
               END DO
               EXIT z_havet_2
            END IF
         END DO z_havet_2
      END DO z_wantt_2

      END IF

      !  Here are the levels that we have from the input for moisture.      

      IF ( flag_soil_levels .EQ. 1 ) THEN
         DO l = 1 , num_sm_levels_input
            zhave(l) = sm_levels_input(l) / 100.
         END DO

           !  Interpolate between the layers we have (zhave) and those that we
           !  want (zs).      

      z_wantm : DO lwant = 1 , num_soil_layers
         z_havem : DO lhave = 1 , num_sm_levels_input -1
            IF ( ( zs(lwant) .GE. zhave(lhave  ) ) .AND. &
                 ( zs(lwant) .LE. zhave(lhave+1) ) ) THEN
               DO j = jts , jte
                  DO i = its , ite
                     smois(i,lwant,j)= ( sm_input(i,lhave,j ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                         sm_input(i,lhave+1,j) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                 ( zhave(lhave+1) - zhave(lhave) )
                  END DO
               END DO
               EXIT z_havem
            END IF
         END DO z_havem
      END DO z_wantm

      ELSE

         zhave(1) = 0.
         DO l = 1 , num_sm_levels_input
            zhave(l+1) = sm_levels_input(l) / 100.
         END DO
         zhave(num_sm_levels_input+2) = 300. / 100.

      z_wantm_2 : DO lwant = 1 , num_soil_layers
         z_havem_2 : DO lhave = 1 , num_sm_levels_input +2
            IF ( ( zs(lwant) .GE. zhave(lhave  ) ) .AND. &
                 ( zs(lwant) .LE. zhave(lhave+1) ) ) THEN
               DO j = jts , jte
                  DO i = its , ite
                     smois(i,lwant,j)= ( sm_input(i,lhave,j ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                         sm_input(i,lhave+1,j) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                 ( zhave(lhave+1) - zhave(lhave) )
                  END DO
               END DO
               EXIT z_havem_2
            END IF
         END DO z_havem_2
      END DO z_wantm_2

      END IF

      

      IF ( flag_sst .EQ. 1 ) THEN
         DO j = jts , jte
            DO i = its , ite
               IF ( landmask(i,j) .LT. 0.5 ) THEN
                  DO l = 1 , num_soil_layers
                     tslb(i,l,j) = sst(i,j)
                     tsk(i,j)    = sst(i,j)
                     smois(i,l,j)= 1.0
                  END DO
               END IF
            END DO
         END DO
      ELSE
         DO j = jts , jte
            DO i = its , ite
               IF ( landmask(i,j) .LT. 0.5 ) THEN
                  DO l = 1 , num_soil_layers
                     tslb(i,l,j)= tsk(i,j)
                     smois(i,l,j)= 1.0
                  END DO
               END IF
            END DO
         END DO
      END IF

      DEALLOCATE (zhave)

   END SUBROUTINE init_soil_3_real

end module module_soil_pre

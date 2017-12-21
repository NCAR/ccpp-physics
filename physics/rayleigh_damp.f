      module rayleigh_damp
      contains

!! \section arg_table_rayleigh_damp_init Argument Table
!!
      subroutine rayleigh_damp_init ()
      end subroutine rayleigh_damp_init


!! \section arg_table_rayleigh_damp_run Argument Table
!! | local var name | longname                                             | description                                          | units      | rank | type    | kind      | intent | optional |
!! |----------------|------------------------------------------------------|------------------------------------------------------|------------|------|---------|-----------|--------|----------|
!! | lsidea         | flag_idealized_physics                               | flag for idealized physics                           | flag       | 0    | logical | default   | in     | F        |
!! | im             | horizontal_loop_extent                               | horizontal loop extent                               | index      | 0    | integer | default   | in     | F        |
!! | ix             | horizontal_dimension                                 | horizontal dimension                                 | index      | 0    | integer | default   | in     | F        |
!! | iy             | horizontal_loop_extent                               | horizontal dimension                                 | index      | 0    | integer | default   | in     | F        |
!! | km             | vertical_dimension                                   | number of vertical layers                            | index      | 0    | integer | default   | in     | F        |
!! | A              | tendency_of_y_wind_due_to_model_physics              | meridional wind tendency due to model physics        | m s-2      | 2    | real    | kind_phys | inout  | F        |
!! | B              | tendency_of_x_wind_due_to_model_physics              | zonal wind tendency due to model physics             | m s-2      | 2    | real    | kind_phys | inout  | F        |
!! | C              | tendency_of_air_temperature_due_to_model_physics     | air temperature tendency due to model physics        | K s-1      | 2    | real    | kind_phys | inout  | F        |
!! | u1             | x_wind                                               | zonal wind                                           | m s-1      | 2    | real    | kind_phys | in     | F        |
!! | v1             | y_wind                                               | meridional wind                                      | m s-1      | 2    | real    | kind_phys | in     | F        |
!! | dt             | time_step_for_physics                                | physics time step                                    | s          | 0    | real    | kind_phys | in     | F        |
!! | cp             | specific_heat_of_dry_air_at_constant_pressure        | specific heat of dry air at constant pressure        | J kg-1 K-1 | 0    | real    | kind_phys | in     | F        |
!! | levr           | number_of_vertical_layers_for_radiation_calculations | number of vertical layers for radiation calculations | index      | 0    | integer | default   | in     | F        |
!! | pgr            | surface_air_pressure                                 | surface pressure                                     | Pa         | 1    | real    | kind_phys | in     | F        |
!! | prsl           | air_pressure                                         | mid-layer pressure                                   | Pa         | 2    | real    | kind_phys | in     | F        |
!! | prslrd0        | pressure_cutoff_for_rayleigh_damping                 | pressure level above which to apply Rayleigh damping | Pa         | 0    | real    | kind_phys | in     | F        |
!! | ral_ts         | time_scale_for_rayleigh_damping                      | time scale for Rayleigh damping                      | d          | 0    | real    | kind_phys | in     | F        |
!!
      subroutine rayleigh_damp_run (
     &           lsidea,IM,IX,IY,KM,A,B,C,U1,V1,DT,CP,
     &           LEVR,pgr,PRSL,PRSLRD0,ral_ts)
!
!   ********************************************************************
! ----->  I M P L E M E N T A T I O N    V E R S I O N   <----------
!
!          --- rayleigh friction with total energy conservation ---
!              ----------------     -----------------------
!
!------ friction coefficient is based on deldif ----
!----------------------------------------------------------------------C
!    USE
!        ROUTINE IS CALLED FROM GBPHYS  (AFTER CALL TO GWDPS)
!
!    PURPOSE
!        USING THE GWD PARAMETERIZATIONS OF PS-GLAS AND PH-
!        GFDL TECHNIQUE.  THE TIME TENDENCIES OF U V ARE 
!        ALTERED TO INCLUDE/MIMIC THE EFFECT OF NON-STATIONARY 
!        GRAVITY WAVE DRAG FROM CONVECTION, FRONTGENOSIS,
!        WIND SHEAR ETC.  LOSS OF KINETIC ENERGY FORM GWD DRAG
!        IS CONVERTED INTO INTERNAL ENERGY.   
!
!  INPUT
!        A(IY,KM)  NON-LIN TENDENCY FOR V WIND COMPONENT
!        B(IY,KM)  NON-LIN TENDENCY FOR U WIND COMPONENT
!        C(IY,KM)  NON-LIN TENDENCY FOR TEMPERATURE
!        U1(IX,KM) ZONAL WIND M/SEC  AT T0-DT
!        V1(IX,KM) MERIDIONAL WIND M/SEC AT T0-DT
!        T1(IX,KM) TEMPERATURE DEG K AT T0-DT
!
!        DT  TIME STEP    SECS
!        pgr(im)          surface pressure (Pa)
!        prsl(IX,KM)      PRESSURE AT MIDDLE OF LAYER (Pa)
!        prslrd0          pressure level above which to apply Rayleigh damping
!        ral_ts           timescale in days for Rayleigh damping
!
!  OUTPUT
!        A, B, C AS AUGMENTED BY TENDENCY DUE TO RAYLEIGH FRICTION
!   ********************************************************************
      USE MACHINE , ONLY : kind_phys
      implicit none
!
      logical,intent(in)                 :: lsidea
      integer,intent(in)                 :: im, ix, iy, km,levr
      real(kind=kind_phys),intent(in)    :: DT, CP, PRSLRD0, ral_ts
      real(kind=kind_phys),intent(in)    :: pgr(im), PRSL(IX,KM)
      real(kind=kind_phys),intent(in)    :: U1(IX,KM), V1(IX,KM)
      real(kind=kind_phys),intent(inout) :: A(IY,KM), B(IY,KM), C(IY,KM)

!--- local variables
      real(kind=kind_phys), parameter :: cons1=1.0, cons2=2.0, half=0.5
      real(kind=kind_phys) DTAUX, DTAUY, wrk1, rtrd1, rfactrd, wrk2
     &,                    ENG0, ENG1, tem1, tem2, dti, hfbcpdt, rtrd
      real(kind=kind_phys) tx1(im)
      integer              i, k
!
      if (lsidea .or. ral_ts <= 0.0 .or. prslrd0 == 0.0) return
!
      RTRD1 = 1.0/(ral_ts*86400) ! RECIPROCAL OF TIME SCALE PER SCALE HEIGHT
                                 ! ABOVE BEGINNING SIGMA LEVEL FOR RAYLEIGH DAMPING
      dti = cons1 / dt
      hfbcpdt = half / (cp*dt)
!
      DO K=1,km
        IF(PRSL(1,K) < PRSLRD0) THEN    ! applied only on constant pressure surfaces
          wrk1 = LOG(PRSLRD0/PRSL(1,K))
          if (k > levr) then
            RTRD = RTRD1 * wrk1 * wrk1
          else
            RTRD = RTRD1 * wrk1
          endif
        ELSE
          RTRD = 0
        ENDIF
        DO I = 1,IM
          RFACTRD = CONS1 / (CONS1+DT*RTRD) - cons1
          DTAUX   = U1(I,k) * RFACTRD
          DTAUY   = V1(I,k) * RFACTRD
          ENG0    = U1(I,K)*U1(I,K) + V1(I,K)*V1(I,K)
          tem1    = U1(I,K) + DTAUX
          tem2    = V1(I,K) + DTAUY
          ENG1    = tem1*tem1 + tem2*tem2
          A(I,K)  = A(I,K) + DTAUY * dti
          B(I,K)  = B(I,K) + DTAUX * dti
          C(I,K)  = C(I,K) + max((ENG0-ENG1),0.0) * hfbcpdt
        ENDDO
      ENDDO


      RETURN
      end subroutine rayleigh_damp_run


!! \section arg_table_rayleigh_damp_finalize Argument Table
!!
      subroutine rayleigh_damp_finalize ()
      end subroutine rayleigh_damp_finalize


      end module rayleigh_damp

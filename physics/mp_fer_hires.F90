!>\file  mp_fer_hires.F90
!! This file contains 

!
module mp_fer_hires
     
      use machine, only : kind_phys

      use module_mp_fer_hires, only : ferrier_init_hr, FER_HIRES

      implicit none

      public :: mp_fer_hires_init, mp_fer_hires_run, mp_fer_hires_finalize
    
      private

      logical :: is_initialized = .False.

   contains

!> \section arg_table_mp_fer_hires_init Argument Table
!! | local_name            | standard_name                             | long_name                                                 | units    | rank |  type                 |   kind    | intent | optional |
!! |-----------------------|-------------------------------------------|-----------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model                 | GFS_control_type_instance                 | Fortran DDT containing FV3-GFS model control parameters   | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | imp_physics           | flag_for_microphysics_scheme              | choice of microphysics scheme                             | flag     |    0 | integer               |           | in     | F        |
!! | imp_physics_fer_hires | flag_for_fer_hires_microphysics_scheme    | choice of Ferrier-Aligo microphysics scheme               | flag     |    0 | integer               |           | in     | F        |
!! | mpicomm               | mpi_comm                                  | MPI communicator                                          | index    |    0 | integer               |           | in     | F        |
!! | mpirank               | mpi_rank                                  | current MPI-rank                                          | index    |    0 | integer               |           | in     | F        |
!! | mpiroot               | mpi_root                                  | master MPI-rank                                           | index    |    0 | integer               |           | in     | F        |
!! | threads               | omp_threads                               | number of OpenMP threads available to scheme              | count    |    0 | integer               |           | in     | F        | 
!! | errmsg                | ccpp_error_message                        | error message for error handling in CCPP                  | none     |    0 | character             | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                           | error flag for error handling in CCPP                     | flag     |    0 | integer               |           | out    | F        |
!!
     subroutine mp_fer_hires_init(Model, imp_physics,                   &
                                  imp_physics_fer_hires,                &
                                  mpicomm, mpirank,mpiroot,             &
                                  threads, errmsg, errflg)

      USE machine,             ONLY : kind_phys
      USE MODULE_MP_FER_HIRES, ONLY : FERRIER_INIT_HR
      USE GFS_typedefs,        ONLY : GFS_control_type
      implicit none

      type(GFS_control_type),         intent(in)    :: Model
      integer,                        intent(in)    :: imp_physics
      integer,                        intent(in)    :: imp_physics_fer_hires
      integer,                        intent(in)    :: mpicomm
      integer,                        intent(in)    :: mpirank
      integer,                        intent(in)    :: mpiroot
      integer,                        intent(in)    :: threads
      character(len=*),               intent(out)   :: errmsg
      integer,                        intent(out)   :: errflg

      ! Local variables
      real(kind=kind_phys)                          :: DT_MICRO

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0
     
      if (is_initialized) return


       ! MZ* temporary 
       if (mpirank==mpiroot) then
         write(0,*) ' ---------------------------------------------------------------------------------------------------------------------'
         write(0,*) ' --- WARNING --- the CCPP Ferrier-Aligo MP scheme is currently under development, use at your own risk --- WARNING ---'
         write(0,*) ' ---------------------------------------------------------------------------------------------------------------------'
       end if
       ! MZ* temporary

       if (imp_physics /= imp_physics_fer_hires ) then
          write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Ferrier-Aligo MP"
          errflg = 1
          return
       end if
          
       !MZ: FERRIER_INIT_HR in NAM_typedefs.F90
       !Model%F_QC=.TRUE.
       !Model%F_QR=.TRUE.
       !Model%F_QS=.TRUE.
       !Model%F_QG=.TRUE.
       !MZ NPRECIP- number of dynamics timesteps between calls to 
       !convection and microphysics
       ! DT_MICRO=Model%NPRECIP*Model%dtp
        DT_MICRO=Model%dtp

       CALL FERRIER_INIT_HR(DT_MICRO,mpicomm,mpirank,threads)

       if (errflg /= 0 ) return
       
       is_initialized = .true.

 
     end subroutine mp_fer_hires_init


!> This is the CCPP-compliant FER_HIRES driver module.
!> \section arg_table_mp_fer_hires_run Argument Table
!! | local_name  | standard_name                                                             | long_name                                                                                          | units      | rank | type     |   kind    | intent | optional |
!! |-------------|---------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|------------|------|----------|-----------|--------|----------|
!! |  ncol       | horizontal_loop_extent                                                    | horizontal loop extent                                                                             | count      |   0  | integer  |           | in     | F        |
!! |  nlev       | vertical_dimension                                                        | number of vertical levels                                                                          | count      |   0  | integer  |           | in     | F        |
!! |  dt         | time_step_for_physics                                                     | physics timestep                                                                                   | s          |   0  | real     | kind_phys | in     | F        |
!! |  slmsk      | sea_land_ice_mask_real                                                    | landmask: sea/land/ice=0/1/2                                                                       | flag       |   1  | real     | kind_phys | in     | F        |
!! |  prsi       | air_pressure_at_interface                                                 | air pressure at model layer interfaces                                                             | Pa         |   2  | real     | kind_phys | in     | F        |
!! |  p_phy      | air_pressure                                                              | mean layer pressure                                                                                | Pa         |   2  | real     | kind_phys | in     | F        |
!! |  t          | air_temperature_updated_by_physics                                        | temperature updated by physics                                                                     | K          |   2  | real     | kind_phys | inout  | F        |
!! |  q          | water_vapor_specific_humidity_updated_by_physics                          | water vapor specific humidity updated by physics                                                   | kg kg-1    |   2  | real     | kind_phys | inout  | F        |
!! |  cwm        | total_cloud_condensate_mixing_ratio_updated_by_physics                    | total cloud condensate mixing ratio (except water vapor) updated by physics                        | kg kg-1    |   2  | real     | kind_phys | inout  | F        |
!! |  train      | accumulated_tendency_of_air_temperature_due_to_FA_scheme                  | accumulated tendency of air temperature due to FA MP scheme                                        | K          |   2  | real     | kind_phys | inout  | F        |
!! |  sr         | ratio_of_snowfall_to_rainfall                                             | snow ratio: ratio of snow to total precipitation  (explicit only)                                  | frac       |   1  | real     | kind_phys | out    | F        |
!! |  f_ice      | mass_fraction_of_ice_water_cloud                                          | mass fraction of ice water cloud                                                                   | frac       |   2  | real     | kind_phys | inout  | F        |
!! |  f_rain     | mass_fraction_of_rain_water_cloud                                         | mass fraction of rain water cloud                                                                  | frac       |   2  | real     | kind_phys | inout  | F        |
!! |  f_rimef    | mass_fraction_of_rime_factor                                              | mass fraction of rime factor                                                                       | frac       |   2  | real     | kind_phys | inout  | F        |
!! |  qc         | cloud_condensed_water_mixing_ratio_updated_by_physics                     | moist (dry+vapor, no condensates) mixing ratio of cloud condensed water updated by physics         | kg kg-1    |   2  | real     | kind_phys | inout  | F        |
!! |  qr         | rain_water_mixing_ratio_updated_by_physics                                | moist (dry+vapor, no condensates) mixing ratio of rain water updated by physics                    | kg kg-1    |   2  | real     | kind_phys | inout  | F        |
!! |  qi         | ice_water_mixing_ratio_updated_by_physics                                 | moist (dry+vapor, no condensates) mixing ratio of ice water updated by physics                     | kg kg-1    |   2  | real     | kind_phys | inout  | F        |
!! |  qs         | snow_water_mixing_ratio_updated_by_physics                                | moist (dry+vapor, no condensates) mixing ratio of snow water updated by physics                    | kg kg-1    |   2  | real     | kind_phys | inout  | F        |
!! |  qg         | graupel_mixing_ratio_updated_by_physics                                   | moist (dry+vapor, no condensates) mixing ratio of graupel updated by physics                       | kg kg-1    |   2  | real     | kind_phys | inout  | F        |
!! |  prec       | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep    | total precipitation amount in each time step                                                       | m          |   1  | real     | kind_phys | inout  | F        |
!! |  acprec     | accumulated_lwe_thickness_of_precipitation_amount                         | accumulated total precipitation                                                                    | m          |   1  | real     | kind_phys | inout  | F        |
!! |  threads    | omp_threads                                                               | number of OpenMP threads available to scheme                                                       | count      |   0  | integer  |           | in     | F        | 
!! |  refl_10cm  | radar_reflectivity_10cm                                                   | instantaneous refl_10cm                                                                            | dBZ        |   2  | real     | kind_phys | inout  | F        |
!! |  rhgrd      | fa_threshold_relative_humidity_for_onset_of_condensation                  | relative humidity threshold parameter for condensation for FA scheme                               | none       |   0  | real     | kind_phys | in     | F        |
!! |  EPSQ       | minimum_value_of_specific_humidity                                        | floor value for specific humidity                                                                  | kg kg-1    |   0  | real     | kind_phys | in     | F        |    
!! |  R_D        | gas_constant_dry_air                                                      | ideal gas constant for dry air                                                                     | J kg-1 K-1 |   0  | real     | kind_phys | in     | F        |
!! |  P608       | ratio_of_vapor_to_dry_air_gas_constants_minus_one                         | (rv/rd) - 1 (rv = ideal gas constant for water vapor)                                              | none       |   0  | real     | kind_phys | in     | F        | 
!! |  CP         | specific_heat_of_dry_air_at_constant_pressure                             | specific heat of dry air at constant pressure                                                      | J kg-1 K-1 |   0  | real     | kind_phys | in     | F        |
!! |  G          | gravitational_acceleration                                                | gravitational acceleration                                                                         | m s-2      |   0  | real     | kind_phys | in     | F        | 
!! |  errmsg     | ccpp_error_message                                                        | error message for error handling in CCPP                                                           | none       |   0  | character| len=*     | out    | F        |
!! |  errflg     | ccpp_error_flag                                                           | error flag for error handling in CCPP                                                              | flag       |   0  | integer  |           | out    | F        |
!!
       SUBROUTINE mp_fer_hires_run(NCOL, NLEV, DT                       &
                         ,SLMSK                                         &
                         ,PRSI,P_PHY                                    &
                         ,T,Q,CWM                                       &
                         ,TRAIN,SR                                      &
                         ,F_ICE,F_RAIN,F_RIMEF                          &
                         ,QC,QR,QI,QS,QG                                &
                         ,PREC,ACPREC                                   &
                         ,threads                                       &
                         ,refl_10cm                                     &
                         ,RHGRD                                         &
                         ,EPSQ,R_D,P608,CP,G                            &
                         ,errmsg,errflg)

!-----------------------------------------------------------------------
      USE MACHINE,    ONLY: kind_phys
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: D_SS=1
!
!------------------------
!***  Argument Variables
!------------------------

      integer,           intent(in   ) :: ncol
      integer,           intent(in   ) :: nlev
      real(kind_phys),   intent(in   ) :: dt
      integer,           intent(in   ) :: threads
      real(kind_phys),   intent(in   ) :: slmsk(1:ncol)
      real(kind_phys),   intent(in   ) :: prsi(1:ncol,1:nlev+1)
      real(kind_phys),   intent(in   ) :: p_phy(1:ncol,1:nlev)
      real(kind_phys),   intent(in   ) :: epsq,r_d,p608,cp,g
      real(kind_phys),   intent(inout) :: t(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: q(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: cwm(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: train(1:ncol,1:nlev)
      real(kind_phys),   intent(out  ) :: sr(1:ncol)
      real(kind_phys),   intent(inout) :: f_ice(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: f_rain(1:ncol,1:nlev) 
      real(kind_phys),   intent(inout) :: f_rimef(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qc(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qr(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qi(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qs(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qg(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: prec(1:ncol)
      real(kind_phys),   intent(inout) :: acprec(1:ncol)
      real(kind_phys),   intent(inout) :: refl_10cm(1:ncol,1:nlev)
      real(kind_phys),   intent(in   ) :: rhgrd
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg
!
!---------------------
!***  Local Variables
!---------------------
!
      integer            :: I,J,K,N
      integer            :: lowlyr(1:ncol)
      real(kind_phys)    :: mprates(1:ncol,1:nlev,d_ss)
      real(kind_phys)    :: sm(1:ncol), xland(1:ncol)
      real(kind_phys)    :: DTPHS,PCPCOL,RDTPHS,TNEW  
      real(kind_phys)    :: ql(1:ncol),tl(1:ncol)
      real(kind_phys)    :: rainnc(1:ncol),rainncv(1:ncol)
      real(kind_phys)    :: snownc(1:ncol),snowncv(1:ncol)
      real(kind_phys)    :: graupelncv(1:ncol)
      real(kind_phys)    :: dz(1:ncol,1:nlev)
      real(kind_phys)    :: pi_phy(1:ncol,1:nlev)
      real(kind_phys)    :: rr(1:ncol,1:nlev)
      real(kind_phys)    :: th_phy(1:ncol,1:nlev)
      real(kind_phys)    :: R_G, CAPPA

! Dimension
      integer            :: ims, ime, jms, jme, lm

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      R_G=1./G
      CAPPA=R_D/CP
      
      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (.not. is_initialized) then
         write(errmsg, fmt='((a))') 'mp_fer_hires_run called before mp_fer_hires_init'
         errflg = 1
         return
      end if 

   
!ZM      NTSD=ITIMESTEP
!ZM presume nphs=1     DTPHS=NPHS*DT
      DTPHS=DT
      RDTPHS=1./DTPHS
!ZM      AVRAIN=AVRAIN+1.

! Set internal dimensions
      ims = 1
      ime = ncol
      jms = 1
      jme = 1
      lm  = nlev


!ZM: module_SOLVER_GRID_COMP.F90 
      DO  I = IMS, IME
          !Sfcprop%sm(i)=1.; if(Sfcprop%slmsk(i) > 0.5 ) Sfcprop%sm(i)=0.
          sm(i) = 1.; if(slmsk(i) > 0.5) sm(i)=0.
      ENDDO

!-----------------------------------------------------------------------
!***  NOTE:  THE NMMB HAS IJK STORAGE WITH LAYER 1 AT THE TOP.
!***         THE WRF PHYSICS DRIVERS HAVE IKJ STORAGE WITH LAYER 1
!***         AT THE BOTTOM.
!-----------------------------------------------------------------------
!
!.......................................................................
!$OMP PARALLEL DO SCHEDULE(dynamic) num_threads(threads) &
!$omp& private(i,k,ql,tl)
!.......................................................................
      DO I=IMS,IME
!
        LOWLYR(I)=1
        XLAND(I)=SM(I)+1.
!
!-----------------------------------------------------------------------
!***   FILL RAINNC WITH ZERO (NORMALLY CONTAINS THE NONCONVECTIVE
!***                          ACCUMULATED RAIN BUT NOT YET USED BY NMM)
!***   COULD BE OBTAINED FROM ACPREC AND CUPREC (ACPREC-CUPREC)
!-----------------------------------------------------------------------
!..The NC variables were designed to hold simulation total accumulations
!.. whereas the NCV variables hold timestep only values, so change below
!.. to zero out only the timestep amount preparing to go into each
!.. micro routine while allowing NC vars to accumulate continually.
!.. But, the fact is, the total accum variables are local, never saved
!.. nor written so they go nowhere at the moment.
!
        RAINNC (I)=0. ! NOT YET USED BY NMM
        RAINNCv(I)=0.
        SNOWNCv(I)=0.
        graupelncv(i) = 0.0
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
        DO K=LM,1,-1   ! We are moving down from the top in the flipped arrays
!
          TL(K)=T(I,K)
          QL(K)=AMAX1(Q(I,K),EPSQ)
!
          RR(I,K)=P_PHY(I,K)/(R_D*TL(K)*(P608*QL(K)+1.))
          PI_PHY(I,K)=(P_PHY(I,K)*1.E-5)**CAPPA
          TH_PHY(I,K)=TL(K)/PI_PHY(I,K)
          DZ(I,K)=(PRSI(I,K+1)-PRSI(I,K))*R_G/RR(I,K)
!
        ENDDO    !- DO K=LM,1,-1
!
      ENDDO    !- DO I=IMS,IME
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!***  CALL MICROPHYSICS
!
!---------------------------------------------------------------------
!*** Update the rime factor array after 3d advection
!---------------------------------------------------------------------
              DO K=1,LM
              DO I=IMS,IME
                IF (QG(I,K)>EPSQ .AND. QS(I,K)>EPSQ) THEN
                  F_RIMEF(I,K)=MIN(50.,MAX(1.,QG(I,K)/QS(I,K)))
                ELSE
                  F_RIMEF(I,K)=1.
                ENDIF
              ENDDO
              ENDDO
!---------------------------------------------------------------------

            CALL FER_HIRES(                                             &
                   DT=dtphs,RHgrd=RHGRD                                 &
                  ,DZ8W=dz,RHO_PHY=rr,P_PHY=p_phy,PI_PHY=pi_phy         &
                  ,TH_PHY=th_phy                                        &
                  ,Q=Q,QC=QC,QS=QS,QR=QR,QT=cwm                         &
                  ,LOWLYR=LOWLYR,SR=SR                                  &
                  ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN                    &
                  ,F_RIMEF_PHY=F_RIMEF                                  &
                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
                  ,threads=threads                                      &
                  ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,LM=LM                &
                  ,D_SS=d_ss,MPRATES=mprates                            &
                  ,refl_10cm=refl_10cm)

!---------------------------------------------------------------------
!*** Calculate graupel from snow array and rime factor
!---------------------------------------------------------------------
              DO K=1,LM
              DO I=IMS,IME
                QG(I,K)=QS(I,K)*F_RIMEF(I,K)
              ENDDO
              ENDDO
!

!.......................................................................
!$OMP PARALLEL DO SCHEDULE(dynamic) num_threads(threads) &
!$omp& private(i,k,TNEW)
!.......................................................................
      DO K=1,LM
        DO I=IMS,IME
!
!-----------------------------------------------------------------------
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, CLOUD WATER, AND HEATING.
!-----------------------------------------------------------------------
!
          TNEW=TH_PHY(I,K)*PI_PHY(I,K)
          TRAIN(I,K)=TRAIN(I,K)+(TNEW-T(I,K))*RDTPHS
          T(I,K)=TNEW
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  UPDATE PRECIPITATION
!-----------------------------------------------------------------------
!
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j,pcpcol)
      DO I=IMS,IME
        PCPCOL=RAINNCV(I)*1.E-3
        PREC(I)=PREC(I)+PCPCOL
        ACPREC(I)=ACPREC(I)+PCPCOL
!
! NOTE: RAINNC IS ACCUMULATED INSIDE MICROPHYSICS BUT NMM ZEROES IT OUT ABOVE
!       SINCE IT IS ONLY A LOCAL ARRAY FOR NOW
!
      ENDDO
!
!-----------------------------------------------------------------------
!
       end subroutine mp_fer_hires_run 


!> \section arg_table_mp_fer_hires_finalize Argument Table
!!
       subroutine mp_fer_hires_finalize ()
       end subroutine mp_fer_hires_finalize

end module mp_fer_hires
     

                  

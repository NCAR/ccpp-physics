!>\file  mp_fer_hires.F90
!! This file contains 

!
module mp_fer_hires
     
      use machine, only : kind_phys, kind_evod

      use module_mp_fer_hires, only : fer_hires_init, FER_HIRES, FER_HIRES_ADVECT

      implicit none

      public :: mp_fer_hires_init, mp_fer_hires_run, mp_fer_hires_finalize
    
      private

      logical :: is_initialized = .False.

   contains

!> \section arg_table_mp_fer_hires_init Argument Table
!! | local_name            | standard_name                             | long_name                                                 | units    | rank |  type                 |   kind    | intent | optional |
!! |-----------------------|-------------------------------------------|-----------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model                 | GFS_control_type_instance                 | Fortran DDT containing FV3-GFS model control parameters   | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | mpirank               | mpi_rank                                  | current MPI-rank                                          | index    |    0 | integer               |           | in     | F        |
!! | mpiroot               | mpi_root                                  | master MPI-rank                                           | index    |    0 | integer               |           | in     | F        |
!! | imp_physics           | flag_for_microphysics_scheme              | choice of microphysics scheme                             | flag     |    0 | integer               |           | in     | F        |
!! | imp_physics_fer_hires | flag_for_fer_hires_microphysics_scheme    | choice of Ferrier-Aligo microphysics scheme               | flag     |    0 | integer               |           | in     | F        |
!! | errmsg                | ccpp_error_message                        | error message for error handling in CCPP                  | none     |    0 | character             | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                           | error flag for error handling in CCPP                     | flag     |    0 | integer               |           | out    | F        |
!!
     subroutine mp_fer_hires_init(Model, mpirank, mpiroot, imp_physics, &
                                  imp_physics_fer_hires, errmsg, errflg)

     ! USE machine,             ONLY : kind_evod
     ! USE MODULE_MP_FER_HIRES, ONLY : FERRIER_INIT_HR
      USE GFS_typedefs,        ONLY : GFS_control_type
      implicit none

      type(GFS_control_type),         intent(in)    :: Model
      integer,                        intent(in)    :: mpirank
      integer,                        intent(in)    :: mpiroot
      integer,                        intent(in)    :: imp_physics
      integer,                        intent(in)    :: imp_physics_fer_hires
      character(len=*),               intent(out)   :: errmsg
      integer,                        intent(out)   :: errflg

      ! Local variables
      real(kind=kind_evod)                          :: DT_MICRO

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

       CALL FERRIER_INIT_HR(DT_MICRO)

       if (errflg /= 0 ) return
       
       is_initialized = .true

 
     end subroutine mp_fer_hires_init


!> This is the CCPP-compliant FER_HIRES driver module.
!> \section arg_table_mp_fer_hires_run Argument Table
!! | local_name  | standard_name                                                             | long_name                                                                                          | units    | rank | type     |   kind    | intent | optional |
!! |-------------|---------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|----------|------|----------|-----------|--------|----------|
!! |  ncol       | horizontal_loop_extent                                                    | horizontal loop extent                                                                             | count    |   0  | integer  |           | in     | F        |   
!! |  nlev       | vertical_dimension                                                        | number of vertical levels                                                                          | count    |   0  | integer  |           | in     | F        |
!! |  dt         | time_step_for_physics                                                     | physics timestep                                                                                   | s        |   0  | real     | kind_phys | in     | F        |
!! |  slmsk      | sea_land_ice_mask_real                                                    | landmask: sea/land/ice=0/1/2                                                                       | flag     |   1  | real     | kind_phys | in     | F        |
!! |  prsi       | air_pressure_at_interface                                                 | air pressure at model layer interfaces                                                             | Pa       |   2  | real     | kind_phys | in     | F        |
!! |  p_phy      | air_pressure                                                              | mean layer pressure                                                                                | Pa       |   2  | real     | kind_phys | in     | F        |
!! |  t          | air_temperature_updated_by_physics                                        | temperature updated by physics                                                                     | K        |   2  | real     | kind_phys | inout  | F        |
!! |  q          | water_vapor_specific_humidity_updated_by_physics                          | water vapor specific humidity updated by physics                                                   | kg kg-1  |   2  | real     | kind_phys | inout  | F        |
!> |  cwm        | total_cloud_condensate_mixing_ratio                                       | total cloud condensate mixing ratio (except water vapor) in NAM                                    | kg kg-1  |   2  | real     | kind_phys | inout  | F        |
!> |  train      | accumulated_tendency_of_air_temperature_due_to_FA_scheme                  |
!> |  sr         | fraction_of_surface_precipitation_associated_with_snow 
!> |  f_ice      | mass_fraction_of_ice_water_cloud                                          | mass fraction of ice water cloud for Ferrier-Aligo MP scheme                                       | frac     |   2  | real     | kind_phys | inout  | F        |
!> |  f_rain     | mass_fraction_of_rain_water_cloud                                         | mass fraction of rain water cloud for Ferrier-Aligo MP scheme                                      | frac     |   2  | real     | kind_phys | inout  | F        |
!> |  f_rimef    | FA_scheme_rime_factor                         
!! |  qc         | cloud_condensed_water_mixing_ratio_updated_by_physics                     | moist (dry+vapor, no condensates) mixing ratio of cloud condensed water updated by physics         | kg kg-1  |   2  | real     | kind_phys | inout  | F        |
!! |  qr         | rain_water_mixing_ratio_updated_by_physics                                | moist (dry+vapor, no condensates) mixing ratio of rain water updated by physics                    | kg kg-1  |   2  | real     | kind_phys | inout  | F        | 
!! |  qi         | ice_water_mixing_ratio_updated_by_physics                                 | moist (dry+vapor, no condensates) mixing ratio of ice water updated by physics                     | kg kg-1  |   2  | real     | kind_phys | inout  | F        | 
!! |  qs         | snow_water_mixing_ratio_updated_by_physics                                | moist (dry+vapor, no condensates) mixing ratio of snow water updated by physics                    | kg kg-1  |   2  | real     | kind_phys | inout  | F        |
!! |  qg         | graupel_mixing_ratio_updated_by_physics                                   | moist (dry+vapor, no condensates) mixing ratio of graupel updated by physics                       | kg kg-1  |   2  | real     | kind_phys | inout  | F        |
!! |  prec       | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep    | total precipitation amount in each time step                                                       | m        |   1  | real     | kind_phys | inout  | F        |
!! |  acprec     | accumulated_lwe_thickness_of_precipitation_amount                         | accumulated total precipitation                                                                    | m        |   1  | real     | kind_phys | inout  | F        |
!! |  refl_10cm  | radar_reflectivity_10cm                                                   | instantaneous refl_10cm                                                                            | dBZ      |   2  | real     | kind_phys | inout  | F        |
!> |  rhgrd      | fa_threshold_relative_humidity_for_onset_of_condensation                  | relative humidity threshold parameter for condensation for FA scheme                               | none     |   0  | real     | kind_phys | in     | F        |           
!! |  errmsg     | ccpp_error_message                                                        | error message for error handling in CCPP                                                           | none     |   0  | character| len=*     | out    | F        |
!! |  errflg     | ccpp_error_flag                                                           | error flag for error handling in CCPP                                                              | flag     |   0  | integer  |           | out    | F        |
!!
       SUBROUTINE mp_fer_hires_run(NCOL, NLEV, DT,                      &
                         ,SLMSK                                         &
                         ,PRSI,P_PHY                                    &
                         ,T,Q,CWM                                       &
                         ,TRAIN,SR                                      &
                         ,F_ICE,F_RAIN,F_RIMEF                          &
                         ,QC,QR,QI,QS,QG                                &  ,NI,NR  
!ZM                         ,F_QC,F_QR,F_QI,F_QS,F_QG,F_NI,F_NR         &
!                         ,has_reqc, has_reqi, has_reqs                 &
                         ,PREC,ACPREC                                   &
                         ,refl_10cm                                     &
                         ,RHGRD                                         &
                         ,errmsg,errflg)
!                         ,IMS,IME,LM,errmsg,errflg)

!          CALL GSMDRIVE(Model%dtp,Model%NPRECIP               &
!                       ,Sfcprop%sm(:),Sfcprop%oro(:)                    &
!                       ,Statein%prsi(:,:),Statein%prsl(:,:)             &
!                       ,Stateout%gt0(:,:)                               &
!                       ,Stateout%gq0(:,:,1),Stateout%gq0(:,:,2)         & !rv CW vs. QC
!                       ,Diag%TRAIN(:,:),Diag%sr(:)                      &
!                       ,Statein%f_ice(:,:),Statein%f_rain(:,:)          &
!                       ,Statein%f_rimef(:,:)                            &
!                       ,Stateout%gq0(:,:,Model%ntcw)                    &
!                       ,Stateout%gq0(:,:,Model%ntrw)                    &
!                       ,Stateout%gq0(:,:,Model%ntiw)                    &
!                       ,Stateout%gq0(:,:,Model%ntsw)                    &
!                       ,Stateout%gq0(:,:,Model%ntgl)                    &
!                       ,Stateout%gq0(:,:,Model%ntinc)                   &
!                       ,Stateout%gq0(:,:,Model%ntrnc)                   &
!                       ,Model%F_QC,Model%F_QR,Model%F_QI,Model%F_QS     &
!                       ,Model%F_QG,Model%F_NI,Model%F_NR                &
!                       ,Model%has_reqc,Model%has_reqi,Model%has_reqs    &
!                       ,Sfcprop%tprcp(:),Diag%totprcp(:)                &
!                       ,Diag%AVRAIN,Statein%refl_10cm(:,:)              &
!                       ,Statein%re_cloud(:,:)                           &
!                       ,Statein%re_ice(:,:),Statein%re_snow(:,:)        &
!                       ,Model%MICROPHYSICS,Model%RHGRD,Diag%TP1(:,:)    &
!                       ,Diag%QP1(:,:),Diag%PSP1(:)                      &
!                       ,ims,ime,LM)


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
      integer,           intent(in   ) :: dt
      real(kind_phys),   intent(in   ) :: slmsk(1:ncol)
      real(kind_phys),   intent(in   ) :: prsi(1:ncol,1:nlev+1)
      real(kind_phys),   intent(in   ) :: p_phy(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: t(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: q(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: cwm(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: train(1:ncol,1:nlev)
      real(kind_phys),   intent(out)   :: sr(1:ncol)
      real(kind_phys),   intent(inout) :: f_ice(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: f_rain(1:ncol,1:nlev) 
      real(kind_phys),   intent(inout) :: f_rimef(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qc(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qr(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qi(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qs(1:ncol,1:nlev)
      real(kind_phys),   intent(inout) :: qg(1:ncol,1:nlev)
!      real(kind_phys),   intent(inout) :: ni(1:ncol,1:nlev)
!      real(kind_phys),   intent(inout) :: nr(1:ncol,1:nlev)
!      logical,           intent(in   ) :: f_qc
!      logical,           intent(in   ) :: f_qr
!      logical,           intent(in   ) :: f_qi
!      logical,           intent(in   ) :: f_qs
!      logical,           intent(in   ) :: f_qg
!      logical,           intent(in   ) :: f_ni
!      logical,           intent(in   ) :: f_nr
!      integer,           intent(in   ) :: has_reqc, has_reqi,has_reqs
      real(kind_phys),   intent(inout) :: prec(1:ncol)
      real(kind_phys),   intent(inout) :: acprec(1:ncol)
      real(kind_phys),   intent(inout) :: refl_10cm(1:ncol,1:nlev)
      real(kind_phys),   intent(in   ) :: rhgrd
!      integer,           intent(in   ) :: IMS,IME,LM      
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
      real(kind_phys)    :: DTPHS,PCPCOL,RDTPHS,TNEW   !qw
      real(kind_phys)    :: ql(1:ncol),tl(1:ncol)
      real(kind_phys)    :: rainnc(1:ncol),rainncv(1:ncol)
      real(kind_phys)    :: snownc(1:ncol),snowncv(1:ncol)
      real(kind_phys)    :: graupelncv(1:ncol)
      real(kind_phys)    :: dz(1:ncol,1:nlev)
      real(kind_phys)    :: pi_phy(1:ncol,1:nlev)
      real(kind_phys)    :: rr(1:ncol,1:nlev)
      real(kind_phys)    :: th_phy(1:ncol,1:nlev)

! Dimension
      integer            :: ims, ime, jms, jme, lm

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      
      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (.not. is_initialized) then
         write(errmsg, fmt='((a))') 'mp_fer_hires_run called before mp_fer_hires_init'
         errflag = 1
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
!$omp parallel do                                                       &
!$omp& private(j,i,k,ql,tl)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
!
        LOWLYR(I,J)=1
        XLAND(I,J)=SM(I,J)+1.
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
        RAINNC (I,J)=0. ! NOT YET USED BY NMM
        RAINNCv(I,J)=0.
        SNOWNCv(I,J)=0.
        graupelncv(i,j) = 0.0
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
        DO K=LM,1,-1   ! We are moving down from the top in the flipped arrays
!
          TL(K)=T(I,J,K)
          QL(K)=AMAX1(Q(I,J,K),EPSQ)
!
          RR(I,J,K)=P_PHY(I,K)/(R_D*TL(K)*(P608*QL(K)+1.))
          PI_PHY(I,J,K)=(P_PHY(I,K)*1.E-5)**CAPPA
          TH_PHY(I,J,K)=TL(K)/PI_PHY(I,J,K)
          DZ(I,J,K)=(PRSI(I,K+1)-PRSI(I,K))*R_G/RR(I,J,K)
!
        ENDDO    !- DO K=LM,1,-1
!
      ENDDO    !- DO I=IMS,IME
      ENDDO    !- DO J=JMS,JME
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
              DO J=JMS,JME
              DO I=IMS,IME
                IF (QG(I,J,K)>EPSQ .AND. QS(I,J,K)>EPSQ) THEN
                  F_RIMEF(I,J,K)=MIN(50.,MAX(1.,QG(I,J,K)/QS(I,J,K)))
                ELSE
                  F_RIMEF(I,J,K)=1.
                ENDIF
              ENDDO
              ENDDO
              ENDDO
!---------------------------------------------------------------------

            CALL FER_HIRES(                                             &
                   DT=dtphs,RHgrd=RHGRD                  &
                  ,DZ8W=dz,RHO_PHY=rr,P_PHY=p_phy,PI_PHY=pi_phy         &
                  ,TH_PHY=th_phy                                        &
                  ,Q=Q,QC=QC,QS=QS,QR=QR,QT=cwm                         &
                  ,LOWLYR=LOWLYR,SR=SR                                  &
                  ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN                    &
                  ,F_RIMEF_PHY=F_RIMEF                                  &
                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
                  ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,LM=LM                &
                  ,D_SS=d_ss,MPRATES=mprates                            &
                  ,refl_10cm=refl_10cm)

!---------------------------------------------------------------------
!*** Calculate graupel from snow array and rime factor
!---------------------------------------------------------------------
              DO K=1,LM
              DO J=JMS,JME
              DO I=IMS,IME
                QG(I,J,K)=QS(I,J,K)*F_RIMEF(I,J,K)
              ENDDO
              ENDDO
              ENDDO
!

!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k,TNEW)
!.......................................................................
      DO K=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
!
!-----------------------------------------------------------------------
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, CLOUD WATER, AND HEATING.
!-----------------------------------------------------------------------
!
          TNEW=TH_PHY(I,J,K)*PI_PHY(I,J,K)
          TRAIN(I,J,K)=TRAIN(I,J,K)+(TNEW-T(I,J,K))*RDTPHS
          T(I,J,K)=TNEW
        ENDDO
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
      DO J=JMS,JME
      DO I=IMS,IME
        PCPCOL=RAINNCV(I,J)*1.E-3
        PREC(I,J)=PREC(I,J)+PCPCOL
        ACPREC(I,J)=ACPREC(I,J)+PCPCOL
!
! NOTE: RAINNC IS ACCUMULATED INSIDE MICROPHYSICS BUT NMM ZEROES IT OUT ABOVE
!       SINCE IT IS ONLY A LOCAL ARRAY FOR NOW
!
      ENDDO
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
     

                  

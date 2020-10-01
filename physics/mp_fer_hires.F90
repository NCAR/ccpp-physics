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
      
      ! * T_ICE - temperature (C) threshold at which all remaining liquid water
      !           is glaciated to ice
      ! * T_ICE_init - maximum temperature (C) at which ice nucleation occurs
      REAL, PUBLIC, PARAMETER :: T_ICE=-40.,                             &
                                 T0C=273.15,                             &
                                 T_ICEK=T0C+T_ICE

   contains

!> This subroutine initialize constants & lookup tables for Ferrier-Aligao MP
!! scheme.
!> \section arg_table_mp_fer_hires_init Argument Table
!! \htmlinclude mp_fer_hires_init.html
!!
     subroutine mp_fer_hires_init(ncol, nlev, dtp, imp_physics,         &
                                  imp_physics_fer_hires,                &
                                  restart,                              &
                                  f_ice,f_rain,f_rimef,                 &
                                  mpicomm, mpirank,mpiroot,             &
                                  threads, errmsg, errflg)

      USE machine,             ONLY : kind_phys
      USE MODULE_MP_FER_HIRES, ONLY : FERRIER_INIT_HR
      implicit none

      integer,                        intent(in)    :: ncol
      integer,                        intent(in)    :: nlev
      real(kind_phys),                intent(in)    :: dtp
      integer,                        intent(in)    :: imp_physics
      integer,                        intent(in)    :: imp_physics_fer_hires
      integer,                        intent(in)    :: mpicomm
      integer,                        intent(in)    :: mpirank
      integer,                        intent(in)    :: mpiroot
      integer,                        intent(in)    :: threads
      logical,                        intent(in)    :: restart
      character(len=*),               intent(out)   :: errmsg
      integer,                        intent(out)   :: errflg
      real(kind_phys),                intent(out), optional  :: f_ice(1:ncol,1:nlev)
      real(kind_phys),                intent(out), optional  :: f_rain(1:ncol,1:nlev)
      real(kind_phys),                intent(out), optional  :: f_rimef(1:ncol,1:nlev)


      ! Local variables
      integer                                       :: ims, ime, lm,i,k
      !real(kind=kind_phys)                          :: DT_MICRO

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0
     
      if (is_initialized) return

      ! Set internal dimensions
      ims = 1
      ime = ncol
      lm  = nlev

       ! MZ* temporary 
       if (mpirank==mpiroot) then
         write(0,*) ' -----------------------------------------------'
         write(0,*) ' ---            !!! WARNING !!!              ---'
         write(0,*) ' --- the CCPP Ferrier-Aligo MP scheme is     ---'
         write(0,*) ' --- currently under development, use at     ---'
         write(0,*) ' --- your own risk .                         ---'
         write(0,*) ' -----------------------------------------------'
       end if
       ! MZ* temporary

       if (imp_physics /= imp_physics_fer_hires ) then
          write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Ferrier-Aligo MP"
          errflg = 1
          return
       end if
          
     !MZ: fer_hires_init() in HWRF
      IF(.NOT.RESTART .AND.  present(F_ICE)) THEN  !HWRF
        write(errmsg,'(*(a))') " WARNING: F_ICE,F_RAIN AND F_RIMEF IS REINITIALIZED "
        DO K = 1,lm
        DO I= ims,ime
          F_ICE(i,k)=0.
          F_RAIN(i,k)=0.
          F_RIMEF(i,k)=1.
        ENDDO
        ENDDO
      ENDIF
      !MZ: fer_hires_init() in HWRF

       CALL FERRIER_INIT_HR(dtp,mpicomm,mpirank,mpiroot,threads,errmsg,errflg)

       if (mpirank==mpiroot) write (0,*)'F-A: FERRIER_INIT_HR finished ...'
       if (errflg /= 0 ) return
       
       is_initialized = .true.

 
     end subroutine mp_fer_hires_init

!>\defgroup hafs_famp HAFS Ferrier-Aligo Cloud Microphysics Scheme
!> This is the CCPP-compliant FER_HIRES driver module.
!> \section arg_table_mp_fer_hires_run Argument Table
!! \htmlinclude mp_fer_hires_run.html
!!
       SUBROUTINE mp_fer_hires_run(NCOL, NLEV, DT ,SPEC_ADV             &
                         ,SLMSK                                         &
                         ,PRSI,P_PHY                                    &
                         ,T,Q,CWM                                       &
                         ,TRAIN,SR                                      &
                         ,F_ICE,F_RAIN,F_RIMEF                          &
                         ,QC,QR,QI,QG                                   & ! wet mixing ratio
                         !,qc_m,qi_m,qr_m                               & 
                         ,PREC                                          &!,ACPREC  -MZ:not used 
                         ,mpirank, mpiroot, threads                     &
                         ,refl_10cm                                     &
                         ,RHGRD,dx                                      &
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
      logical,           intent(in   ) :: spec_adv
      integer,           intent(in   ) :: mpirank
      integer,           intent(in   ) :: mpiroot
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
      real(kind_phys),   intent(inout) :: qg(1:ncol,1:nlev) ! QRIMEF

      real(kind_phys),   intent(inout) :: prec(1:ncol)
!      real(kind_phys)                  :: acprec(1:ncol)   !MZ: change to local
      real(kind_phys),   intent(inout) :: refl_10cm(1:ncol,1:nlev)
      real(kind_phys),   intent(in   ) :: rhgrd
      real(kind_phys),   intent(in   ) :: dx(1:ncol)
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg
!
!---------------------
!***  Local Variables
!---------------------
!
      integer            :: I,J,K,N
      integer            :: lowlyr(1:ncol)
      integer            :: dx1
      !real(kind_phys)    :: mprates(1:ncol,1:nlev,d_ss)
      real(kind_phys)    :: DTPHS,PCPCOL,RDTPHS,TNEW  
      real(kind_phys)    :: ql(1:nlev),tl(1:nlev)
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

! Use the dx of the 1st i point to set an integer value of dx to be used for
! determining where RHgrd should be set to 0.98 in the coarse domain when running HAFS.
      DX1=NINT(DX(1))

!-----------------------------------------------------------------------
!***  NOTE:  THE NMMB HAS IJK STORAGE WITH LAYER 1 AT THE TOP.
!***         THE WRF PHYSICS DRIVERS HAVE IKJ STORAGE WITH LAYER 1
!***         AT THE BOTTOM.
!-----------------------------------------------------------------------
!.......................................................................
      DO I=IMS,IME
!
        LOWLYR(I)=1
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
!          TL(K)=T(I,K)
!          QL(K)=AMAX1(Q(I,K),EPSQ)
!
          RR(I,K)=P_PHY(I,K)/(R_D*T(I,K)*(P608*AMAX1(Q(I,K),EPSQ)+1.))
          PI_PHY(I,K)=(P_PHY(I,K)*1.E-5)**CAPPA
          TH_PHY(I,K)=T(I,K)/PI_PHY(I,K)
          DZ(I,K)=(PRSI(I,K)-PRSI(I,K+1))*R_G/RR(I,K)

!
!***  CALL MICROPHYSICS

!MZ* in HWRF
!-- 6/11/2010: Update cwm, F_ice, F_rain and F_rimef arrays
         cwm(I,K)=QC(I,K)+QR(I,K)+QI(I,K)
         IF (QI(I,K) <= EPSQ) THEN
            F_ICE(I,K)=0.
            F_RIMEF(I,K)=1.
            IF (T(I,K) < T_ICEK) F_ICE(I,K)=1.
         ELSE
            F_ICE(I,K)=MAX( 0., MIN(1., QI(I,K)/cwm(I,K) ) )
            F_RIMEF(I,K)=QG(I,K)/QI(I,K)
         ENDIF
         IF (QR(I,K) <= EPSQ) THEN
            F_RAIN(I,K)=0.
         ELSE
            F_RAIN(I,K)=QR(I,K)/(QR(I,K)+QC(I,K))
         ENDIF

        end do
      enddo

!---------------------------------------------------------------------
!*** Update the rime factor array after 3d advection
!---------------------------------------------------------------------
!MZ* in namphysics
!              DO K=1,LM
!              DO I=IMS,IME
!                IF (QG(I,K)>EPSQ .AND. QI(I,K)>EPSQ) THEN
!                  F_RIMEF(I,K)=MIN(50.,MAX(1.,QG(I,K)/QI(I,K)))
!                ELSE
!                  F_RIMEF(I,K)=1.
!                ENDIF
!              ENDDO
!              ENDDO


!---------------------------------------------------------------------
        
            CALL FER_HIRES(                                             &
                   DT=dtphs,RHgrd=RHGRD                                 &
                  ,DZ8W=dz,RHO_PHY=rr,P_PHY=p_phy,PI_PHY=pi_phy         &
                  ,TH_PHY=th_phy,T_PHY=t                                &
                  ,Q=Q,QT=cwm                                           &
                  ,LOWLYR=LOWLYR,SR=SR                                  &
                  ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN                    &
                  ,F_RIMEF_PHY=F_RIMEF                                  &
                  ,QC=QC,QR=QR,QS=QI                                    &
                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
                  ,threads=threads                                      &
                  ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,LM=LM                &
                  ,D_SS=d_ss                                            &
                  ,refl_10cm=refl_10cm,DX1=DX1)


!.......................................................................

!MZ*
!Aligo Oct-23-2019 
! - Convert dry qc,qr,qi back to wet mixing ratio
!    DO K = 1, LM
!     DO I= IMS, IME
!       qc_m(i,k) = qc(i,k)/(1.0_kind_phys+q(i,k))
!       qi_m(i,k) = qi(i,k)/(1.0_kind_phys+q(i,k))
!       qr_m(i,k) = qr(i,k)/(1.0_kind_phys+q(i,k))
!     ENDDO
!    ENDDO
    


!-----------------------------------------------------------
      DO K=1,LM
        DO I=IMS,IME

!---------------------------------------------------------------------
!*** Calculate graupel from total ice array and rime factor
!---------------------------------------------------------------------

!MZ            
            IF (SPEC_ADV) then
                QG(I,K)=QI(I,K)*F_RIMEF(I,K)
            ENDIF

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

!
!-----------------------------------------------------------------------
!***  UPDATE PRECIPITATION
!-----------------------------------------------------------------------
!
      DO I=IMS,IME
        PCPCOL=RAINNCV(I)*1.E-3        !MZ:unit:m
        PREC(I)=PREC(I)+PCPCOL
!MZ        ACPREC(I)=ACPREC(I)+PCPCOL     !MZ: not used 
!
! NOTE: RAINNC IS ACCUMULATED INSIDE MICROPHYSICS BUT NMM ZEROES IT OUT ABOVE
!       SINCE IT IS ONLY A LOCAL ARRAY FOR NOW
!
      ENDDO
!-----------------------------------------------------------------------
!
       end subroutine mp_fer_hires_run 


!> \section arg_table_mp_fer_hires_finalize Argument Table
!!
       subroutine mp_fer_hires_finalize ()
       end subroutine mp_fer_hires_finalize

end module mp_fer_hires

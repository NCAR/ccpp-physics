!>\file  mp_fer_hires.F90
!! This file contains the Ferrier-Aligo microphysics scheme driver.

!
module mp_fer_hires

      use machine, only : kind_phys

      use module_mp_fer_hires, only : ferrier_init_hr, FER_HIRES,       &
                                      ferhires_finalize

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
                                  imp_physics_fer_hires, restart,       &
                                  mpicomm, mpirank,mpiroot,             &
                                  threads, errmsg, errflg)

      USE mpi_f08
      USE machine,             ONLY : kind_phys
      USE MODULE_MP_FER_HIRES, ONLY : FERRIER_INIT_HR
      implicit none

      integer,                        intent(in)    :: ncol
      integer,                        intent(in)    :: nlev
      real(kind_phys),                intent(in)    :: dtp
      integer,                        intent(in)    :: imp_physics
      integer,                        intent(in)    :: imp_physics_fer_hires
      type(MPI_Comm),                 intent(in)    :: mpicomm
      integer,                        intent(in)    :: mpirank
      integer,                        intent(in)    :: mpiroot
      integer,                        intent(in)    :: threads
      logical,                        intent(in)    :: restart
      character(len=*),               intent(out)   :: errmsg
      integer,                        intent(out)   :: errflg

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

        if (mpirank==mpiroot) write (0,*) 'F-A: calling FERRIER_INIT_HR ...'
       CALL FERRIER_INIT_HR(dtp,mpicomm,mpirank,mpiroot,threads,errmsg,errflg)

       if (mpirank==mpiroot) write (0,*)'F-A: FERRIER_INIT_HR finished ...'
       if (errflg /= 0 ) return

       is_initialized = .true.

     end subroutine mp_fer_hires_init

!>\defgroup hafs_famp HWRF Ferrier-Aligo Microphysics Scheme
!> This is the CCPP-compliant FER_HIRES driver module.
!> \section arg_table_mp_fer_hires_run Argument Table
!! \htmlinclude mp_fer_hires_run.html
!!
       SUBROUTINE mp_fer_hires_run(NCOL, NLEV, DT ,SPEC_ADV             &
                         ,SLMSK                                         &
                         ,PRSI,P_PHY                                    &
                         ,T,Q                                           &
                         ,TRAIN,SR                                      &
                         ,QC,QR,QI,QG                                   &
                         ,PREC                                          &
                         ,mpirank, mpiroot, threads                     &
                         ,refl_10cm                                     &
                         ,RHGRD,dx                                      &
                         ,EPSQ,R_D,P608,CP,G                            &
                         ,ten_t,ten_qv,ten_ql,ten_qr,ten_qi,ten_qg      &
                         ,ten_q,errmsg,errflg)

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
      real(kind_phys),   intent(in   ) :: slmsk(:)
      real(kind_phys),   intent(in   ) :: prsi(:,:)
      real(kind_phys),   intent(in   ) :: p_phy(:,:)
      real(kind_phys),   intent(in   ) :: epsq,r_d,p608,cp,g
      real(kind_phys),   intent(in   ) :: t(:,:)
      real(kind_phys),   intent(in   ) :: q(:,:)
      real(kind_phys),   intent(inout), optional :: train(:,:)
      real(kind_phys),   intent(out  ) :: sr(:)
      real(kind_phys),   intent(in   ) :: qc(:,:)
      real(kind_phys),   intent(in   ) :: qr(:,:)
      real(kind_phys),   intent(in   ) :: qi(:,:)
      real(kind_phys),   intent(in   ) :: qg(:,:) ! QRIMEF

      real(kind_phys),   intent(inout) :: prec(:)
      real(kind_phys),   intent(inout) :: refl_10cm(:,:)
      real(kind_phys),   intent(in   ) :: rhgrd
      real(kind_phys),   intent(in   ) :: dx(:)
      real(kind_phys),   intent(  out) :: ten_t(:,:)
      real(kind_phys),   intent(  out) :: ten_qv(:,:)
      real(kind_phys),   intent(  out) :: ten_ql(:,:)
      real(kind_phys),   intent(  out) :: ten_qr(:,:)
      real(kind_phys),   intent(  out) :: ten_qi(:,:)
      real(kind_phys),   intent(  out) :: ten_qg(:,:)
      real(kind_phys),   intent(  out) :: ten_q(:,:,:)
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
      real(kind_phys)    :: PCPCOL
      real(kind_phys)    :: ql(1:nlev),tl(1:nlev)
      real(kind_phys)    :: rainnc(1:ncol),rainncv(1:ncol)
      real(kind_phys)    :: snownc(1:ncol),snowncv(1:ncol)
      real(kind_phys)    :: graupelncv(1:ncol)
      real(kind_phys)    :: train_phy(1:ncol,1:nlev)
      real(kind_phys)    :: f_ice(1:ncol,1:nlev)
      real(kind_phys)    :: f_rain(1:ncol,1:nlev)
      real(kind_phys)    :: f_rimef(1:ncol,1:nlev)
      real(kind_phys)    :: cwm(1:ncol,1:nlev)
      real(kind_phys)    :: new_t(1:ncol,1:nlev)
      real(kind_phys)    :: new_qv(1:ncol,1:nlev)
      real(kind_phys)    :: new_ql(1:ncol,1:nlev)
      real(kind_phys)    :: new_qr(1:ncol,1:nlev)
      real(kind_phys)    :: new_qi(1:ncol,1:nlev)
      real(kind_phys)    :: new_qg(1:ncol,1:nlev)

! Dimension
      integer            :: ims, ime, lm

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      new_t = t
      new_qv = q
      new_ql = qc
      new_qr = qr
      new_qi = qi
      new_qg = qg
      
      ten_t = 0.0
      ten_q = 0.0 !set tendency of entire tracer array to zero to make sure that those tracers not affected by this scheme do not change when tendencies are applied
      ten_qv = 0.0
      ten_ql = 0.0
      ten_qr = 0.0
      ten_qi = 0.0
      ten_qg = 0.0
      
      ! Check initialization state
      if (.not. is_initialized) then
         write(errmsg, fmt='((a))') 'mp_fer_hires_run called before mp_fer_hires_init'
         errflg = 1
         return
      end if

! Set internal dimensions
      ims = 1
      ime = ncol
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
        DO K=LM,1,-1   !mz* We are moving down from the top in the flipped arrays

!***  CALL MICROPHYSICS

!MZ* in HWRF
!-- 6/11/2010: Update cwm, F_ice, F_rain and F_rimef arrays
         cwm(I,K)=new_ql(I,K)+new_qr(I,K)+new_qi(I,K)
         IF (new_qi(I,K) <= EPSQ) THEN
            F_ICE(I,K)=0.
            F_RIMEF(I,K)=1.
            IF (new_t(I,K) < T_ICEK) F_ICE(I,K)=1.
         ELSE
            F_ICE(I,K)=MAX( 0., MIN(1., new_qi(I,K)/cwm(I,K) ) )
            F_RIMEF(I,K)=new_QG(I,K)!/QI(I,K)
         ENDIF
         IF (new_qr(I,K) <= EPSQ) THEN
            F_RAIN(I,K)=0.
         ELSE
            F_RAIN(I,K)=new_qr(I,K)/(new_qr(I,K)+new_ql(I,K))
         ENDIF

        ENDDO

      ENDDO

!---------------------------------------------------------------------
!aligo
       DO K = 1, LM
       DO I= IMS, IME
         cwm(i,k) = cwm(i,k)/(1.0_kind_phys-new_qv(i,k))
         new_qr(i,k) = new_qr(i,k)/(1.0_kind_phys-new_qv(i,k))
         new_qi(i,k) = new_qi(i,k)/(1.0_kind_phys-new_qv(i,k))
         new_ql(i,k) = new_ql(i,k)/(1.0_kind_phys-new_qv(i,k))
       ENDDO
       ENDDO
!aligo
!---------------------------------------------------------------------

            CALL FER_HIRES(                                             &
                   DT=DT,RHgrd=RHGRD                                    &
                  ,PRSI=prsi,P_PHY=p_phy,T_PHY=new_t                    &
                  ,Q=new_qv,QT=cwm                                      &
                  ,LOWLYR=LOWLYR,SR=SR,TRAIN_PHY=train_phy              &
                  ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN                    &
                  ,F_RIMEF_PHY=F_RIMEF                                  &
                  ,QC=new_ql,QR=new_qr,QS=new_qi                        &
                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
                  ,threads=threads                                      &
                  ,IMS=IMS,IME=IME,LM=LM                                &
                  ,D_SS=d_ss                                            &
                  ,refl_10cm=refl_10cm,DX1=DX1)


!.......................................................................

!MZ*
!Aligo Oct-23-2019
! - Convert dry qc,qr,qi back to wet mixing ratio
    DO K = 1, LM
     DO I= IMS, IME
        new_ql(i,k) = new_ql(i,k)/(1.0_kind_phys+new_qv(i,k))
        new_qi(i,k) = new_qi(i,k)/(1.0_kind_phys+new_qv(i,k))
        new_qr(i,k) = new_qr(i,k)/(1.0_kind_phys+new_qv(i,k))
     ENDDO
    ENDDO

!-----------------------------------------------------------
      DO K=1,LM
        DO I=IMS,IME

!---------------------------------------------------------------------
!*** Calculate graupel from total ice array and rime factor
!---------------------------------------------------------------------

!MZ
            IF (SPEC_ADV) then
                new_QG(I,K)=new_qi(I,K)*F_RIMEF(I,K)
                ten_qg(i,k) = (new_qg(i,k) - qg(i,k))/dt
            ENDIF

!
!-----------------------------------------------------------------------
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, CLOUD WATER, AND HEATING.
!-----------------------------------------------------------------------
!
          ten_t(i,k) = (new_t(i,k) - t(i,k))/dt
          ten_qv(i,k) = (new_qv(i,k) - q(i,k))/dt
          ten_ql(i,k) = (new_ql(i,k) - qc(i,k))/dt
          ten_qr(i,k) = (new_qr(i,k) - qr(i,k))/dt
          ten_qi(i,k) = (new_qi(i,k) - qi(i,k))/dt
          TRAIN(I,K)=TRAIN(I,K)+TRAIN_PHY(I,K)
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
!! \htmlinclude mp_fer_hires_finalize.html
!!
       subroutine mp_fer_hires_finalize (errmsg,errflg)
         implicit none

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (.not.is_initialized) return

         call ferhires_finalize()

         is_initialized = .false.
       end subroutine mp_fer_hires_finalize

end module mp_fer_hires

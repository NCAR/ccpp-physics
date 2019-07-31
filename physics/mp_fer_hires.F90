!>\file  mp_fer_hires.F90
!! This file contains 

!
module mp_fer_hires
     
      use machine, only : kind_phys

      use module_mp_fer_hires, only : fer_hires_init, FER_HIRES, FER_HIRES_ADVECT

      implicit none

      public :: mp_fer_hires_init, mp_fer_hires_run, mp_fer_hires_finalize
    
      private

      logical :: is_initialized = .False.

   contains

!>\section arg_table_mp_fer_hires_init Argument Table
!!
     subroutine mp_fer_hires_init(ncol, nlev, F_ICE_PHY, F_RAIN_PHY,     &)
 
       implicit none
       ! Interface variables
       integer,       intent(in)       :: ncol
       integer,       intent(in)       :: nlev

       logical,       intent(in)       :: allowed_to_read
       integer,       intent(in)       ::

       ! Local variables: dimensions used in fer_hires_init
       integer                         :: ids, ide, jds, jde, kds, kde, &
                                          ims, ime, jms, jme, kms, kme, &
                                          its, ite, jts, jte, kts, kte


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
          
       ! Set internal dimensions
       
       ids = 1
       ims = 1
       its = 1
       ide = ncol
       ime = ncol
       ite = ncol
       jds = 1
       jms = 1
       jts = 1
       jde = 1
       jme = 1
       jte = 1
       kds = 1
       kms = 1
       kts = 1
       kde = nlev
       kme = nlev
       kte = nlev

!MZ: FERRIER_INIT_HR in NAM_typedefs.F90
       Model%F_QC=.TRUE.
       Model%F_QR=.TRUE.
       Model%F_QS=.TRUE.
       Model%F_QG=.TRUE.
       !MZ NPRECIP- number of dynamics timesteps between calls to 
       !convection and microphysics
       DT_MICRO=Model%NPRECIP*Model%dtp

       CALL FERRIER_INIT_HR(DT_MICRO)

       

       if (errflg /= 0 ) return
       
       is_initialized = .true

 
     end subroutine mp_fer_hires_init


!> This is the CCPP-compliant FER_HIRES driver module.
      subroutine mp_fer_hires_run ()
   
      implicit none
      
      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (.not. is_initialized) then
         write(errmsg, fmt='((a))') 'mp_fer_hires_run' called before mp_fer_hires_init'
         errflag = 1
         return
      end if 

!HWRF v4.0 version
!           CALL FER_HIRES(                                              &
!                   ITIMESTEP=itimestep,DT=dt, GID=id        &
!                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
!                  ,DZ8W=dz8w,RHO_PHY=rho,P_PHY=p                        &
!                  ,PI_PHY=pi_phy,TH_PHY=th                              &
!                  ,QV=qv_curr                                           &
!                  ,QT=qt_curr                                           &
!                  ,LOWLYR=LOWLYR,SR=SR                                  &
!                  ,F_ICE_PHY=F_ICE_PHY,F_RAIN_PHY=F_RAIN_PHY            &
!                  ,F_RIMEF_PHY=F_RIMEF_PHY                              &
!                  ,QC=qc_curr,QR=Qr_curr,QI=Qi_curr                     &
!                  ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde    &
!                  ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme    &
!                  ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte    &
!                  ,errmsg=errmsg, errflg=errflg )

!-------------------------------------------------------------------
!*** NAMPHYSICS below
!-------------------------------------------------------------------
   
        NTSD=ITIMESTEP
        DTPHS=NPHS*DT
!       RDTPHS=1./DTPHS
!       AVRAIN=AVRAIN+1.
       

!zm:......................................................................
!$omp parallel do                                                       &
!$omp& private(j,i,k,ql,tl)
!.......................................................................

        DO J=JMS,JME
        DO I=IMS,IME
!
           LOWLYR(I,J)=1
!          XLAND(I,J)=SM(I,J)+1.   


!------------------------------------------------------------------
!*** FILL RAINNC WITH ZERO (NORMALLY CONTAINS THE NONCONVECTIVE
!***                        ACCUMULATED RAIN BUT NOT YET USED BY NMM)
!*** COULD BE OBTAINED FROM ACPREC AND CUPREC (ACPREC-CUPREC)
!------------------------------------------------------------------
!.. The NC variables were designed to hold simulation total accumulations
!.. whereas the NCV variables hold timestep only values, so change below
!.. to zero out only the timestep amount preparing to go into each 
!.. micro routine while allowing NC vars to accumulate continually. 
!.. But, the fact is, the total accum variables are local, never saved 
!.. nor written so they go nowhere at the moment.
!
           RAINNC (I,J)=0. ! NOT YET USED BY NMM
           RAINNCv(I,J)=0.
           SNOWNCv(I,J)=0.
           graupelncv(i,j) =0.0

!-------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-------------------------------------------------------------------
!
           DO K=LM,1,-1           ! We are moving down from the top in the flipped arrays
!        
!zm             TL(K)=T(I,J,K)
!             QL(K)=AMAX1(Q(I,J,K),EPSQ)
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

!zm: OMP.......................................................................
!$omp end parallel do
!.......................................................................

             


!-------------------------------------------------------------------
!*** NOTE: THE NMMB HAS IJK STORAGE WITH LAYER 1 AT THE TOP.
!*** THE WRF PHYSICS DRIVER HAVE IKJ STORAGE WITH LAYER 1
!*** AT THE BOTTOM.
!-------------------------------------------------------------------
!
! module_MICROPHYSICS.F90 in namephysics
! Update the rime factor array after 3d advection
        Do K =1, LM
        DO J=JMS,JME
        DO I=IMS,IME
          IF (QG(I,J,K)>EPSQ  .AND. QS(I,J,K)>EPSQ) THEN
             F_RIMEF(I,J,K)=MIN(50., MAX(1.,QG(I,J,K)/QS(I,J,K)))
          ELSE 
             F_RIMEF(I,J,K)=1.
          ENDIF 
        ENDDO
        ENDDO
        ENDDO

!---------------------------------------------------------------------
! *** CALL MICROPHYSICS

            CALL FER_HIRES(                                             &
                   ITIMESTEP=ntsd,DT=dtphs,RHgrd=RHGRD                  &
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
!---------------------------------------------------------------------

!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k,TNEW,MP_TTEN)
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
!         IF (USE_RADAR) THEN
!           MP_TTEN=(TNEW-T(I,J,K))*RDTPHS
!           IF(DFI_TTEN(I,J,K)>MP_TTEN.AND.DFI_TTEN(I,J,K)<0.01        &
!                                     .AND.MP_TTEN<0.0018)THEN
!             MP_TTEN=DFI_TTEN(I,J,K)
!           END IF
!           T(I,J,K)=T(I,J,K)+MP_TTEN/RDTPHS
!         ELSE
            T(I,J,K)=TNEW
!         ENDIF
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




       end subroutine mp_fer_hires_run 


!> \section arg_table_mp_fer_hires_finalize Argument Table
!!
       subroutine mp_fer_hires_finalize ()
       end subroutine mp_fer_hires_finalize

end module mp_fer_hires
     

                  

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

! FERRIER_INIT_HR in NAM_typedefs.F90

       
       ! Call init to initialize constans & lookup tables for microphysics
       CALL fer_hires_init (MPDT,DT,DX,DY,LOWLYR,restart,         &
                            allowed_to_read,                        &
                            ids, ide, jds, jde, kds, kde,           &
                            ims, ime, jms, jme, kms, kme,           &
                            its, ite, jts, jte, kts, kte,           &
                            F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY)

       
       is_initialized = .true

 
     end subroutine mp_fer_hires_init



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

      
!      if (is_fer_mp_hires_advect) then
!           CALL FER_HIRES_ADVECT(                                       &
!                   ITIMESTEP=itimestep,DT=dt,DX=dx,DY=dy, GID=id        &
!                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
!                  ,DZ8W=dz8w,RHO_PHY=rho,P_PHY=p                        &
!                  ,PI_PHY=pi_phy,TH_PHY=th                              &
!                  ,QV=qv_curr                                           &
!                  ,LOWLYR=LOWLYR,SR=SR                                  &
!                  ,QC=qc_curr,QR=Qr_curr,QI=Qi_curr,QRIMEF=qrimef_curr  &
!                  ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde    &
!                  ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme    &
!                  ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte    &
!                  ,errmsg=errmsg, errflg=errflg )
!      else 

!-------------------------------------------------------------------
!NAM: 'fer_hires': module_MICROPHYSICS.F90
!*** Update the rime factor array after 3d advection
!-------------------------------------------------------------------
!           DO K=1, LM
!           DO J=JMS,JME
!           DO I=IMS,IME
!             IF  (QG(I,J,K)>EPSQ .AND. QS(I,J,K)>EPSQ) THEN
!                F_RIMEF(I,J,K)=MIN(50.,MAX(1.,QG(I,J,K)/QS(I,J,K)))
!             ELSE
!                F_RIMEF(I,J,K)=1.
!             ENDIF
!           ENDDO
!           ENDDO
!           ENDDO
          

!HWRF
           CALL FER_HIRES(                                              &
!MZ                   ITIMESTEP=itimestep,DT=dt,DX=dx,DY=dy, GID=id        &
                   ITIMESTEP=itimestep,DT=dt, GID=id        &
                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
                  ,DZ8W=dz8w,RHO_PHY=rho,P_PHY=p                        &
                  ,PI_PHY=pi_phy,TH_PHY=th                              &
                  ,QV=qv_curr                                           &
                  ,QT=qt_curr                                           &
                  ,LOWLYR=LOWLYR,SR=SR                                  &
                  ,F_ICE_PHY=F_ICE_PHY,F_RAIN_PHY=F_RAIN_PHY            &
                  ,F_RIMEF_PHY=F_RIMEF_PHY                              &
                  ,QC=qc_curr,QR=Qr_curr,QI=Qi_curr                     &
                  ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde    &
                  ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme    &
                  ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte    &
                  ,errmsg=errmsg, errflg=errflg )
         endif

       end subroutine mp_fer_hires_run 


!> \section arg_table_mp_fer_hires_finalize Argument Table
!!
       subroutine mp_fer_hires_finalize ()
       end subroutine mp_fer_hires_finalize

end module mp_fer_hires
     

                  

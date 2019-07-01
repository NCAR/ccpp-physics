!see mp_thompson.F90

      subroutine mp_fer_hires_run ()




             if (is_fer_mp_hires_advect) then
                  CALL FER_HIRES_ADVECT(                                &
                         ITIMESTEP=itimestep,DT=dt,DX=dx,DY=dy, GID=id  &
                  ,RAINNC=rainnc,RAINNCV=rainncv                             &
                  ,DZ8W=dz8w,RHO_PHY=rho,P_PHY=p,PI_PHY=pi_phy,TH_PHY=th     &
                  ,QV=qv_curr                                        &
                  ,LOWLYR=LOWLYR,SR=SR                               &
                  ,QC=qc_curr,QR=Qr_curr,QI=Qi_curr,QRIMEF=qrimef_curr   &
                  ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
                  ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
                  ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
                                                                     )




             else 
                CALL FER_HIRES(                                      &
                   ITIMESTEP=itimestep,DT=dt,DX=dx,DY=dy, GID=id &
                  ,RAINNC=rainnc,RAINNCV=rainncv                     &
                  ,DZ8W=dz8w,RHO_PHY=rho,P_PHY=p,PI_PHY=pi_phy,TH_PHY=th &
                  ,QV=qv_curr                                        &
                  ,QT=qt_curr                                        &
                  ,LOWLYR=LOWLYR,SR=SR                               &
                  ,F_ICE_PHY=F_ICE_PHY,F_RAIN_PHY=F_RAIN_PHY         &
                  ,F_RIMEF_PHY=F_RIMEF_PHY                           &
                  ,QC=qc_curr,QR=Qr_curr,QI=Qi_curr                  &
                  ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
                  ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
                  ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
                                                                     )
              endif

                  

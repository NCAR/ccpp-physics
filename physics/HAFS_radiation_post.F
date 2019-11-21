      btimx=now_time()
!
!     Pass in XTIME (elapsed time from start of parent) to compute
!     the time passed into the zenith angle code consistently between
!     RDTEMP and RADIATION.

      call check_grid(grid,config_flags,'before RDTEMP', &
            ids,ide,jds,jde,kds,kde, &
            ims,ime,jms,jme,kms,kme, &
            its,ite,jts,jte,kts,kte)
      CALL RDTEMP(grid%ntsd,GRID%DT,GRID%JULDAY,GRID%JULYR              &
     &           ,GRID%XTIME,IHRST,grid%glat,grid%glon                  &
     &           ,grid%czen,grid%czmean,grid%t                          &
     &           ,grid%rswtt,grid%rlwtt,grid%hbm2                       &
     &           ,IDS,IDF,JDS,JDF,KDS,KDE                               &
     &           ,IMS,IME,JMS,JME,KMS,KME                               &
     &           ,ITS,ITE,JTS,JTE,KTS,KTE)
!
      rdtemp_tim=rdtemp_tim+now_time()-btimx
!
!



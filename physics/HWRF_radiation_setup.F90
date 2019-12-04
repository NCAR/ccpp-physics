!>\file 
!! This file include HWRF  RRTMG variable setup in module_initialize_real
!albedo
!sm

     ! julday,julyr in interger
     ! main/real_nmm
      current_date = start_date

      !  After current_date has been set, fill in the julgmt stuff.

      CALL geth_julgmt ( config_flags%julyr , config_flags%julday , &
                                              config_flags%gmt )
      !xtime
      
      ! ihrst: solve_nmm
      ihrst = gmt
    
      !phy_init
      !- calculate sm - sea mask
       DO j = jts, MIN(jte,jde-1)
         DO i = its, MIN(ite,ide-1)
           if (grid%landmask(I,J) .gt. 0.5) grid%sm(I,J)=0.
           if (grid%landmask(I,J) .le. 0.5) grid%sm(I,J)=1.
           IF(grid%sm(I,J).GT.0.9) THEN
             IF (grid%xice(I,J) .gt. 0) then
                grid%si(I,J)=1.0
             ENDIF


         ENDDO
       ENDDO

     !- calculate sice - sea ice mask

     
       




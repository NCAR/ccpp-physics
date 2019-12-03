
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
     
       




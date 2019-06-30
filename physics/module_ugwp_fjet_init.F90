!=========================================================================
!
!    ugwp_fjet_init
!    
!=========================================================================  
     
   module ugwp_fjet_init 
      implicit none
      real    ::  eff_fj                      ! scale factors for conv GWs
      integer ::  nwfj                       ! number of waves      
      integer ::  nazfj                      ! number of azimuths
      integer ::  nstfj                      ! flag for stochastic choice of launch level above Conv-cloud  
!        
      real, parameter    ::  fjet_trig=0.       ! if ( abs(frgf) > fjet_trig ) launch GW-packet
      
      
      real, parameter    :: cmin =  2.5
      real, parameter    :: cmax = 67.5  
      real               :: dc
      real, allocatable  :: ch_fjet(:) , spf_fjet(:) 
      real, allocatable  :: xaz_fjet(:), yaz_fjet(:)
     contains    
     subroutine init_fjet_gws(nwaves, nazdir, nstoch, effac, &
                              lonr, kxw )
     use ugwp_common,  only : pi2, arad
     use ugwp_conv_init, only: init_nazdir

     implicit none
     
      integer :: nwaves, nazdir, nstoch
      integer :: lonr 
      real    :: kxw,  effac , chk
      
      integer :: k
              
      nwfj  =  nwaves
      nazfj =  nazdir
      nstfj =  nstoch      
      eff_fj  =  effac
      
       if (.not. allocated(ch_fjet))  allocate (ch_fjet(nwaves))
       if (.not. allocated(spf_fjet)) allocate (spf_fjet(nwaves))
       if (.not. allocated(xaz_fjet)) allocate (xaz_fjet(nazdir))
       if (.not. allocated(yaz_fjet)) allocate (yaz_fjet(nazdir)) 
           
      dc = (cmax-cmin)/float(nwaves-1)
      do k = 1,nwaves
         
         chk        = cmin + (k-1)*dc
         ch_fjet(k) =  chk
         spf_fjet(k) =  1.0	 
      enddo	
      call init_nazdir(nazdir,  xaz_fjet,  yaz_fjet)
           
     end subroutine init_fjet_gws

    end module ugwp_fjet_init      
!

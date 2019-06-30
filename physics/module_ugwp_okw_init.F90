      module ugwp_okw_init 
!=========================================================================    
      implicit none
     
      real    ::  eff_okw                     ! scale factors for conv GWs
      integer ::  nwokw                       ! number of waves      
      integer ::  nazokw                      ! number of azimuths
      integer ::  nstokw                      ! flag for stochastic choice of launch level above Conv-cloud  
!        
      real, parameter    ::  okw_trig=0.       ! if ( abs(okwp) > okw_trig ) launch GW-packet
      
      
      real, parameter    :: cmin =  2.5
      real, parameter    :: cmax = 67.5  
      real               :: dc
      real, allocatable  :: ch_okwp(:),   spf_okwp(:) 
      real, allocatable  :: xaz_okwp(:),   yaz_okwp(:)

     contains
!              
      subroutine init_okw_gws(nwaves, nazdir, nstoch, effac, &
                              lonr, kxw )
			      
     use ugwp_common,  only : pi2, arad
     use ugwp_conv_init, only: init_nazdir
 
      implicit none
     
      integer :: nwaves, nazdir, nstoch
      integer :: lonr 
      real    :: kxw,  effac , chk
      
      integer :: k
              
      nwokw  =  nwaves
      nazokw =  nazdir
      nstokw =  nstoch      
      eff_okw  =  effac
      
       if (.not. allocated(ch_okwp))    allocate (ch_okwp(nwaves))
       if (.not. allocated(spf_okwp))   allocate (spf_okwp(nwaves))
       if (.not. allocated(xaz_okwp)) allocate (xaz_okwp(nazdir))
       if (.not. allocated(yaz_okwp)) allocate (yaz_okwp(nazdir))      
      dc = (cmax-cmin)/float(nwaves-1)
      do k = 1,nwaves
         chk =  cmin + (k-1)*dc
         ch_okwp(k) = chk
	 spf_okwp(k) = 1. 
      enddo
      
      call init_nazdir(nazdir,  xaz_okwp,  yaz_okwp)           
         
      end subroutine init_okw_gws
          
      end module ugwp_okw_init     
          
!=============================== end of GW  sources    

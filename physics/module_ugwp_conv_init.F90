! ========================================================================= 
!
!    ugwp_conv_init
!    
!=========================================================================     
    module ugwp_conv_init      

     implicit none
      real    ::  eff_con                   ! scale factors for conv GWs
      integer ::  nwcon                     ! number of waves      
      integer ::  nazcon                    ! number of azimuths
      integer ::  nstcon                    ! flag for stochastic choice of launch level above Conv-cloud    
      real    ::  con_dlength 
      real    ::  con_cldf
      
      real, parameter    :: cmin =  5  !2.5
      real, parameter    :: cmax = 95. !82.5
      real, parameter    :: cmid = 22.5  
      real, parameter    :: cwid = cmid  
      real, parameter    :: bns = 2.e-2, bns2 = bns*bns, bns4=bns2*bns2
      real, parameter    :: mstar = 6.28e-3/2.  !  2km        
      real               :: dc
      
      real, allocatable  :: ch_conv(:), spf_conv(:)    
      real, allocatable  :: xaz_conv(:),   yaz_conv(:)                         
     contains 
!     
     subroutine init_conv_gws(nwaves, nazdir, nstoch, effac, &
                              lonr, kxw, cgwf )
     use ugwp_common,  only : pi2, arad     
     implicit none
     
      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real    :: cgwf(2)    
      real    :: kxw,  effac 
      real    :: work1 = 0.5
      real    :: chk, tn4, snorm      
      integer :: k   
         
      nwcon  =  nwaves
      nazcon =  nazdir
      nstcon =  nstoch      
      eff_con  =  effac
      
      con_dlength = pi2*arad/float(lonr)
      con_cldf   = cgwf(1)    * work1 + cgwf(2) *(1.-work1)     
!
! allocate & define spectra in "selected direction": "dc" "ch(nwaves)"
!
       if (.not. allocated(ch_conv))  allocate (ch_conv(nwaves))
       if (.not. allocated(spf_conv)) allocate (spf_conv(nwaves)) 
       if (.not. allocated(xaz_conv)) allocate (xaz_conv(nazdir))
       if (.not. allocated(yaz_conv)) allocate (yaz_conv(nazdir)) 
          
      dc = (cmax-cmin)/float(nwaves-1)
!
! we may use different spectral "shapes"
! for example FVS-93 "Desabeius"
!  E(s=1, t=3,m, w, k) ~ m^s/(m*^4 + m^4) ~ m^-3 saturated tail
!  
      do k = 1,nwaves
         chk         = cmin + (k-1)*dc
	 tn4 = (mstar*chk)**4
         ch_conv(k)  =  chk
         spf_conv(k) =  bns4*chk/(bns4+tn4)	 
      enddo
         
	 snorm = sum(spf_conv)
	 spf_conv = spf_conv/snorm*1.5
	 
      call init_nazdir(nazdir,  xaz_conv,  yaz_conv)
     end subroutine init_conv_gws  

      subroutine init_nazdir(naz,  xaz,  yaz)
      use ugwp_common , only : pi2
      implicit none
      integer :: naz
      real, dimension(naz) :: xaz,  yaz
      integer :: idir
      real    :: phic, drad
      drad  = pi2/float(naz)
      if (naz.ne.4) then
        do idir =1, naz
         Phic = drad*(float(idir)-1.0)
         xaz(idir) = cos(Phic)
         yaz(idir) = sin(Phic)
        enddo
      else
!                         if (naz.eq.4) then
          xaz(1) = 1.0     !E
          yaz(1) = 0.0
          xaz(2) = 0.0
          yaz(2) = 1.0     !N
          xaz(3) =-1.0     !W
          yaz(3) = 0.0
          xaz(4) = 0.0
          yaz(4) =-1.0     !S
      endif
      end  subroutine init_nazdir     
        
    end module ugwp_conv_init 
!=========================================================================

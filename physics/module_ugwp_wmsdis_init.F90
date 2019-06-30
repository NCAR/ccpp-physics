  module ugwp_wmsdis_init 
 
    implicit none
    
      real,      parameter  :: maxdudt = 250.e-5	    
            
      real,      parameter  :: hpscale= 7000., rhp2   =  0.5/hpscale     
      real,      parameter  :: omega2 = 2.*6.28/86400
      real,      parameter  :: gptwo=2.0       
        
      real ,    parameter   :: dked_min =0.01
      real ,    parameter   :: gssec = (6.28/30.)**2        ! max-value for bn2
      real ,    parameter   :: bv2min = (6.28/60./120.)**2  ! min-value for bn2  7.6(-7)  2 hrs
      real, parameter       :: minvel = 0.5     
      
! 
! make parameter list that will be passed to SOLVER
!
!      
        
      real, parameter       :: v_kxw  = 6.28e-3/200.
      real, parameter       :: v_kxw2 = v_kxw*v_kxw
      real, parameter       :: tamp_mpa = 30.e-3
      real, parameter       :: zfluxglob= 3.75e-3
               
      real ,     parameter  :: nslope=1                ! the GW sprctral slope at small-m     
      integer, parameter    :: klaunch=55       ! 32 - ~ 1km ;55 - 5.5 km ; 52 4.7km ; 60-7km index for selecting launch level 
      integer, parameter    :: ilaunch=klaunch                  
      
      integer  , parameter  :: iazidim=4       ! number of azimuths
      integer  , parameter  :: incdim=25       !number of discrete cx - spectral elements in launch spectrum
      real ,     parameter  :: ucrit2  =0.5   
      
      real ,     parameter  :: zcimin = ucrit2  
      real ,     parameter  :: zcimax = 125.0
      real ,     parameter  :: zgam =   0.25              
      real ,     parameter  :: zms_l  = 2000.0        
        
      real                  :: gw_eff
      
!===========================================================================
      integer  :: nwav, nazd, nst
      real     :: eff  
      
      real                :: zaz_fct , zms 
      real, allocatable   :: zci(:),  zci4(:), zci3(:),zci2(:), zdci(:) 
      real, allocatable   :: zcosang(:), zsinang(:)      
      contains
!============================================================================
     subroutine initsolv_wmsdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw)
     
!        call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
!         knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw)
!
     use ugwp_common, only :   pi, pi2  
     implicit none
!     
!input -control for solvers:
!      nwaves, nazdir, nstoch, effac, do_physb, kxw
!     
!
     integer  :: me, master
     integer  :: nwaves, nazdir
     integer  :: nstoch
     real     :: effac  
     logical  :: do_physb
     real     :: kxw
!     
!locals     
!
      integer   :: inc, jk, jl, iazi  
!        


      real :: zang, zang1, znorm
      real :: zx1, zx2, ztx, zdx, zxran, zxmin, zxmax, zx
      real :: zpexp     
       
     if( nwaves == 0) then
!
!     redefine from the deafault
!       
       nwav = incdim
       nazd = iazidim
       nst  = 0
       eff  = 1.0
       gw_eff = eff
     else
!
! from input.nml
!     
       nwav = nwaves
       nazd = nazdir
       nst  = nstoch
       gw_eff  = effac  
     endif 
     
      allocate ( zci(nwav),  zci4(nwav), zci3(nwav),zci2(nwav), zdci(nwav)  )   
      allocate ( zcosang(nazd), zsinang(nazd) ) 
            
 
      
      if (me == master) then
         print *, 'ugwp_v0: init_gw_wmsdis_control '         
         print *, 'ugwp_v0: WMSDIS launch layer ',      klaunch 
         print *, 'ugwp_v0: WMSDID tot_mflux in mpa',   tamp_mpa*1000.
       endif     
          
        
        zpexp=gptwo/2.0                    ! gptwo=2 , zpexp = 1.
                  	 	   
!
!      set up azimuth directions and some trig factors
!
!
           zang=pi2/float(nazd)
      
! get normalization factor to ensure that the same amount of momentum
! flux is directed (n,s,e,w) no mater how many azimuths are selected.
!
        znorm=0.0 
        do iazi=1, nazd
         zang1= (iazi-1)*zang
         zcosang(iazi)=cos(zang1)
         zsinang(iazi)=sin(zang1)
         znorm=znorm+abs(zcosang(iazi))
        enddo
        zaz_fct= 1.0	
        zaz_fct= 2./znorm            ! correction factot for azimuthal sums
	 
!       define coordinate transform for "Ch"   ....x = 1/c stretching transform
!       -----------------------------------------------     
! note that this is expresed in terms of the intrinsic phase speed
! at launch ci=c-u_o so that the transformation is identical
! see eq. 28-30 of scinocca 2003.   x = 1/c stretching transform
!
          zxmax=1.0 /zcimin
          zxmin=1.0 /zcimax
          zxran=zxmax-zxmin
          zdx=zxran/real(nwav-1)                            ! dkz  
!
          zx1=zxran/(exp(zxran/zgam)-1.0 )                    !zgam =1./4.
          zx2=zxmin-zx1
  
! 
! computations for zci =1/zx
!                                   if(lgacalc) zgam=(zxmax-zxmin)/log(zxmax/zxmin)
!                                   zx1=zxran/(exp(zxran/zgam)-1.0_jprb)
!                                   zx2=zxmin-zx1
        zms=2.*pi/zms_l    
        do inc=1, nwav
          ztx=real(inc-1)*zdx+zxmin
          zx = zx1*exp((ztx-zxmin)/zgam)+zx2                          !eq. 29 of scinocca 2003
          zci(inc)=1.0 /zx                                            !eq. 28 of scinocca 2003
        zdci(inc)=zci(inc)**2*(zx1/zgam)*exp((ztx-zxmin)/zgam)*zdx    !eq. 30 of scinocca 2003
          zci4(inc)=(zms*zci(inc))**4	
	  zci2(inc)=(zms*zci(inc))**2
	  zci3(inc)=(zms*zci(inc))**3 
        enddo
!	
!
!  all done and print-out 
!	
!	
	
         if (me == master) then  
         print *  
         print *,  'ugwp_v0: zcimin=' , zcimin
         print *,  'ugwp_v0: zcimax=' , zcimax
         print *,  'ugwp_v0: cd_crit=',  zgam                        ! m/s precision for crit-level 
         print *,  'ugwp_v0: launch-level',  ilaunch
         print *, ' ugwp_v0 zms_l=', zms_l
         print *, ' ugwp_vgw  nslope=', nslope

         print *
         endif      
 
        
     end subroutine initsolv_wmsdis
!
! make a list of  all-initilized parameters needed for "gw_solver_wmsdis"
!     
     
  end module ugwp_wmsdis_init  
!=========================================================================   

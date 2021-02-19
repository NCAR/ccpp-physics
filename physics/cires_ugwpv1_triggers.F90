module cires_ugwpv1_triggers

        use machine,            only: kind_phys

contains

!
!
!      
!>\ingroup cires_ugwp_run
!> @{
!!
!!
      subroutine slat_geos5_tamp_v1(im, tau_amp, xlatdeg, tau_gw)
!=================
! V1: GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch) =rho*<u'w'>
!=================
      implicit none
      integer :: im     
      real(kind=kind_phys)    :: tau_amp, xlatdeg(im), tau_gw(im)
      real(kind=kind_phys)    :: latdeg, flat_gw, tem
      integer :: i
      
!
! if-lat
!
      do i=1, im
        latdeg = abs(xlatdeg(i))    
        if (latdeg < 15.3) then
          tem = (latdeg-3.0) / 8.0
          flat_gw = 0.75 * exp(-tem * tem)
          if (flat_gw < 1.2 .and. latdeg <= 3.0) flat_gw = 0.75
        elseif (latdeg <  31.0 .and. latdeg >=  15.3) then
           flat_gw =  0.10
        elseif (latdeg <  60.0 .and. latdeg >=  31.0) then
          tem = (latdeg-60.0) / 23.0
          flat_gw =  0.50 * exp(- tem * tem)
        elseif (latdeg >=  60.0) then
          tem = (latdeg-60.0) / 70.0
          flat_gw =  0.50 * exp(- tem * tem)
        endif
        tau_gw(i) = tau_amp*flat_gw 
      enddo
!      
      end subroutine slat_geos5_tamp_v1
!      
      subroutine slat_geos5_2020(im, tau_amp, xlatdeg, tau_gw)
!=================================================================
! modified for FV3GFS-127L/C96 QBO-experiments
! GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch)
!================================================================
      implicit none
      integer :: im     
      real(kind=kind_phys)    :: tau_amp, xlatdeg(im), tau_gw(im)
      real(kind=kind_phys)    :: latdeg, flat_gw, tem
      real(kind=kind_phys), parameter    :: fampqbo = 1.25    ! 1.5     
      real(kind=kind_phys), parameter    :: famp60S = 1.0     ! 1.5
      real(kind=kind_phys), parameter    :: famp60N = 1.0     ! 1.0
      real(kind=kind_phys), parameter    :: famp30  = 0.25    ! 0.4
        
      real(kind=kind_phys), parameter    :: swid15  = 12.5  
      real(kind=kind_phys), parameter    :: swid60S  = 30.0   ! 40
      real(kind=kind_phys), parameter    :: swid60N  = 25.0   ! 30                 
      integer :: i
!    
!
!
      do i=1, im
      
        latdeg = abs(xlatdeg(i))    
        if (latdeg < 15.3) then
          tem = (latdeg-3.0) / swid15
          flat_gw = fampqbo * exp(-tem * tem)
          if (latdeg <= 3.0) flat_gw = fampqbo
        elseif (latdeg <  31.0 .and. latdeg >=  15.3) then
           flat_gw =  famp30
        elseif (latdeg <  60.0 .and. latdeg >=  31.0) then
          tem = (latdeg-60.0) / 23.0
          flat_gw =   famp60N* exp(- tem * tem)
        elseif (latdeg >=  60.0) then
          tem = (latdeg-60.0) /swid60N
          flat_gw =  famp60N * exp(- tem * tem)
        endif

          if (xlatdeg(i) <= -31.0) then   
!                 
            if (latdeg <  60.0 .and. latdeg >=  31.0) then
              tem = (latdeg-60.0) / 23.0
              flat_gw =  famp60S * exp(- tem * tem)           
            endif           
            if (latdeg >=  60.0) then
               tem = (latdeg-60.0) /swid60S
               flat_gw =  famp60S * exp(- tem * tem)
            endif
            
          endif
        tau_gw(i) = tau_amp*flat_gw 
      enddo
!      
      end subroutine slat_geos5_2020   
      
           
      subroutine slat_geos5(im, xlatdeg, tau_gw)
      
!=================
!
! WAM: GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch) =rho*<u'w'>
! 
!=================
      implicit none
      integer :: im     
      real(kind=kind_phys)  :: xlatdeg(im)          
      real(kind=kind_phys)  :: tau_gw(im)
      real(kind=kind_phys)  :: latdeg
      real(kind=kind_phys), parameter  :: tau_amp = 3.5e-3    ! 3.5 mPa 
      real(kind=kind_phys)             :: trop_gw,  flat_gw 
      integer :: i
!
! if-lat
!
      trop_gw = 0.75
      do i=1, im
      latdeg = xlatdeg(i)    
      if (-15.3 < latdeg .and. latdeg < 15.3) then
          flat_gw = trop_gw*exp(-( (abs(latdeg)-3.)/8.0)**2)
          if (flat_gw < 1.2 .and. abs(latdeg) <= 3.) flat_gw = trop_gw
      else if (latdeg > -31. .and. latdeg <= -15.3) then
          flat_gw =  0.10
      else if (latdeg <  31. .and. latdeg >=  15.3) then
           flat_gw =  0.10
      else if (latdeg > -60. .and. latdeg <= -31.) then
        flat_gw =  0.50*exp(-((abs(latdeg)-60.)/23.)**2)
      else if (latdeg <  60. .and. latdeg >=  31.) then
        flat_gw =  0.50*exp(-((abs(latdeg)-60.)/23.)**2)
      else if (latdeg <= -60.) then
         flat_gw =  0.50*exp(-((abs(latdeg)-60.)/70.)**2)
      else if (latdeg >=  60.) then
         flat_gw =  0.50*exp(-((abs(latdeg)-60.)/70.)**2)
      end if
      tau_gw(i) = tau_amp*flat_gw 
      enddo
!      
      end subroutine slat_geos5      
      
!===============================================
!
!   Spontaneous GW triggers by dynamical inbalances (OKW, fronts/jets, and convection)
!       not activated due to "limited" set of GFS-physics 
!       statein-type ( needs horizontal gradients of winds and temperature, humodity)
!
!===============================================       
      subroutine get_spectra_tau_convgw &
      (nw, im, levs, dcheat,  scheat, precip, icld, xlatd, sinlat, coslat,taub, klev, if_src, nf_src)
!
! temporarily can put GEOS-5/MERRA-2 GW-lat dependent function
!      
      integer                   :: nw, im, levs 
      integer,dimension(im,3)   :: icld   
      real(kind=kind_phys), dimension(im, levs) :: dcheat,  scheat   
      real(kind=kind_phys), dimension(im)       :: precip, xlatd, sinlat, coslat      
      real(kind=kind_phys), dimension(im)       :: taub
      integer, dimension(im)    :: klev, if_src
      integer ::  nf_src
!
! locals
      real(kind=kind_phys), parameter           ::  precip_max = 100.   ! mm/day
      real(kind=kind_phys), parameter           ::  tau_amp    = 3.5e-3  ! 3.5 mPa   
       
      integer  :: i, k, klow, ktop, kmid
      real(kind=kind_phys)     :: dtot, dmax, daver
!      
      nf_src = 0
      if_src(1:im) = 0 
      taub(1:im)   = 0.0
      do i=1, im
        klow = icld(i,1)
        ktop = icld(i,2)
        kmid=  icld(i,3)
        if (klow == -99 .and.  ktop == -99) then 
          cycle
        else 
          klev(i) = ktop
          k = klow
          klev(i) = k
          dmax = abs(dcheat(i,k) + scheat(i,k))
          do k=klow+1, ktop
            dtot =abs(dcheat(i,k) + scheat(i,k))
            if ( dtot >  dmax) then
              klev(i) = k
              dmax =  dtot
            endif  
          enddo
!	  
! klev  as  max( dcheat(i,k) + scheat)
!	vertical width of conv-heating
!  
! counts/triiger=1   & taub(i) 
!
          nf_src = nf_src +1
          if_src(i) = 1          
          taub(i) = tau_amp* precip(i)/precip_max*coslat(i)
        endif 

      enddo
!	
!  
!      
      call Slat_geos5(im, xlatd, taub)
      nf_src =im
      do i=1, im
        if_src(i) = 1
        klev(i)   = 127-45
      enddo
    
! with info on precip/clouds/dc_heat create Bulk 
!        taub(im), klev(im)      
!   
!      print *, '   get_spectra_tau_convgw '  
      end subroutine get_spectra_tau_convgw
!
      subroutine get_spectra_tau_nstgw(nw, im, levs, trig_fgf, xlatd, sinlat, coslat, taub, klev, if_src, nf_src)
      integer                      :: nw, im, levs 
      real(kind=kind_phys), dimension(im, levs)    :: trig_fgf
!      real(kind=kind_phys), dimension(im, levs+1)  :: pint      
      real(kind=kind_phys), dimension(im)          :: xlatd, sinlat, coslat        
      real(kind=kind_phys), dimension(im)          :: taub
      integer, dimension(im)          :: klev, if_src
      integer ::  nf_src
! locals
      real(kind=kind_phys), parameter           ::  tlim_fgf = 100.     ! trig_fgf > tlim_fgf, launch waves should scale-dependent
      real(kind=kind_phys), parameter           ::  tau_amp   = 3.5e-3  ! 3.5 mPa
      real(kind=kind_phys), parameter           ::  pmax = 750.e2,  pmin = 100.e2  
      integer, parameter        ::  klow =127-92,   ktop=127-45
      integer, parameter        ::  kwidth = ktop-klow+1   
      integer  :: i, k,  kex
      real(kind=kind_phys)     :: dtot, dmax, daver      
      real(kind=kind_phys)     :: fnorm, tau_min     
      nf_src = 0      
      if_src(1:im) = 0 
      taub(1:im)   = 0.0
      fnorm = 1.0 / float(kwidth)
      tau_min = tau_amp*fnorm
      do i=1, im
!
! only trop-c fjets so find max(trig_fgf) => klev
!                      use abs-values to scale tau_amp
!
  
        k = klow
        klev(i) = k
        dmax = abs(trig_fgf(i,k))
        kex  = 0
        if (dmax >= tlim_fgf) kex = kex+1 
        do k=klow+1, ktop
          dtot = abs(trig_fgf(i,k))
          if (dtot >= tlim_fgf) kex = kex+1 
          if ( dtot >  dmax) then
            klev(i) = k
            dmax    =  dtot
          endif  
        enddo

        if (dmax .ge. tlim_fgf) then
          nf_src = nf_src +1
          if_src(i) = 1          
          taub(i) = tau_min*float(kex)          !* precip(i)/precip_max*coslat(i)
        endif 

      enddo
!     
!      print *, '   get_spectra_tau_nstgw '
      call Slat_geos5(im, xlatd, taub)
      nf_src =im
      do i=1, im
        if_src(i) = 1
        klev(i) = 127-45                  ! FV3-127L
      enddo
!                 
      end subroutine get_spectra_tau_nstgw   
! 
      subroutine get_spectra_tau_okw(nw, im, levs,  trig_okw, xlatd, sinlat, coslat, taub, klev, if_src, nf_src)
      integer                      :: nw, im, levs 
      real(kind=kind_phys), dimension(im, levs)    :: trig_okw
!      real(kind=kind_phys), dimension(im, levs+1)  :: pint    
      real(kind=kind_phys), dimension(im)          :: xlatd, sinlat, coslat     
      real(kind=kind_phys), dimension(im)          :: taub
      integer, dimension(im)       :: klev, if_src
      integer ::  nf_src
! locals
      real(kind=kind_phys), parameter           ::  tlim_okw = 100.                  ! trig_fgf > tlim_fgf, launch GWs should scale-dependent
      real(kind=kind_phys), parameter           ::  tau_amp   = 35.e-3               ! 35 mPa  
      real(kind=kind_phys), parameter           ::  pmax = 750.e2,  pmin = 100.e2  
      integer, parameter        ::  klow =127-92,   ktop=127-45      ! for FV3-127L
      integer, parameter        ::  kwidth = ktop-klow+1   
      integer  :: i, k,  kex
      real(kind=kind_phys)     :: dtot, dmax, daver      
      real(kind=kind_phys)     :: fnorm, tau_min       
      
      nf_src = 0  
      if_src(1:im) = 0    
      taub(1:im) = 0.0    
       fnorm = 1./float(kwidth)
      tau_min = tau_amp*fnorm 
      print *, '   get_spectra_tau_okwgw '   
      do i=1, im
        k = klow
        klev(i) = k	     
        dmax = abs(trig_okw(i,k))
        kex = 0
        if (dmax >= tlim_okw) kex = kex+1 
        do k=klow+1, ktop
          dtot = abs(trig_okw(i,k))
          if (dtot >= tlim_fgf ) kex = kex+1 
          if ( dtot >  dmax) then
            klev(i) = k
            dmax =  dtot
          endif  
        enddo
!	  
        if (dmax >= tlim_okw) then
          nf_src = nf_src + 1
          if_src(i) = 1          
          taub(i) = tau_min*float(kex)          !* precip(i)/precip_max*coslat(i)
        endif 

      enddo 
      print *, '   get_spectra_tau_okwgw '           
      end subroutine get_spectra_tau_okw   

end module cires_ugwpv1_triggers

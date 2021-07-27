!
      subroutine slat_geos5_tamp_v0(im, tau_amp, xlatdeg, tau_gw)
!=================
! GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch) =rho*<u'w'>
!=================
      implicit none
      integer :: im     
      real    :: tau_amp, xlatdeg(im), tau_gw(im)
      real    :: latdeg, flat_gw, tem
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
      end subroutine slat_geos5_tamp_v0
      
      subroutine slat_geos5_v0(im, xlatdeg, tau_gw)
!=================
! GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch) =rho*<u'w'>
!=================
      implicit none
      integer :: im     
      real  :: xlatdeg(im)          
      real  :: tau_gw(im)
      real  :: latdeg
      real, parameter  :: tau_amp = 100.e-3
      real             :: trop_gw,  flat_gw 
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
      end subroutine slat_geos5_v0
!      
      subroutine init_nazdir_v0(naz,  xaz,  yaz)
      use ugwp_common_v0 , only : pi2
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
      end  subroutine init_nazdir_v0

        subroutine ugwp_triggers
        implicit none
        write(6,*) ' physics-based triggers for UGWP ' 
        end subroutine ugwp_triggers
!	
     SUBROUTINE  subs_diag_geo(nx, ny,  lat, lon, rlat, rlon, dy, dx,  &
                               cosv, rlatc, brcos, brcos2, dlam1, dlam2, dlat, divJp, divJm)
      use ugwp_common , only : deg_to_rad
      
      implicit none
       integer :: nx, ny
       real    :: lon(nx), lat(ny)
       real    :: rlon(nx), rlat(ny) , cosv(ny), tanlat(ny)  
       real    :: rlatc(ny-1), brcos(ny), brcos2(ny)              
       real    :: earth_r, ra1, ra2, dx, dy, dlat
       real    :: dlam1(ny), dlam2(ny),  divJp(ny), divJm(ny)
       integer :: j
!
!    specify common constants and
!    geometric factors to compute deriv-es etc ...
!    coriolis coslat tan etc...
!
      earth_r = 6370.e3
      ra1     = 1.0 / earth_r
      ra2     = ra1*ra1
!
      rlat   = lat*deg_to_rad
      rlon   = lon*deg_to_rad
      tanlat = atan(rlat)
      cosv   =  cos(rlat)     
      dy     = rlat(2)-rlat(1)
      dx     = rlon(2)-rlon(1)
!
      do j=1, ny-1
        rlatc(j) = 0.5 * (rlat(j)+rlat(j+1))
      enddo
!
      do j=2, ny-1
        brcos(j) = 1.0 / cos(rlat(j))*ra1
      enddo       
       
      brcos(1)  = brcos(2)
      brcos(ny) = brcos(ny-1)
      brcos2    = brcos*brcos
!
      dlam1 = brcos  / (dx+dx)
      dlam2 = brcos2 / (dx*dx)
  
      dlat = ra1 / (dy+dy)
  
      divJp = dlat*cosv
      divJM = dlat*cosv
!
      do  j=2, ny-1 
        divJp(j) = dlat*cosv(j+1)/cosv(j) 
        divJM(j) = dlat*cosv(j-1)/cosv(j) 
      enddo
      divJp(1)  = divjp(2)   !*divjp(1)/divjp(2)
      divJp(ny) = divjp(1)
      divJM(1)  = divjM(2)   !*divjM(1)/divjM(2)
      divJM(ny) = divjM(1)
!
      return
      end SUBROUTINE  subs_diag_geo
!      
      subroutine get_xy_pt(V, Vx, Vy, nx, ny, dlam1, dlat)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! compute for each Vert-column: grad(V)
! periodic in X and central diff ...
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none
      integer :: nx, ny
      real    :: V(nx, ny), dlam1(ny), dlat
      real    :: Vx(nx, ny), Vy(nx, ny)      
      integer :: i, j
      do i=2, nx-1 
        Vx(i,:) = dlam1(:)*(V(i+1,:)-V(i-1,:))
      enddo
      Vx(1,:)  =  dlam1(:)*(V(2,:)-V(nx,:))
      Vx(nx,:) =  dlam1(:)*(V(1,:)-V(nx-1,:))

      do j=2, ny-1
        Vy(:,j) = dlat*(V(:,j+1)-V(:, j-1))
      enddo 
      Vy(:, 1) = dlat*2.*(V(:,2)-V(:,1))
      Vy(:,ny) = dlat*2.*(V(:,ny)-V(:,ny-1))

    end subroutine get_xy_pt
    
    subroutine get_xyd_wind( V, Vx, Vy, Vyd, nx, ny, dlam1, dlat, divJp, divJm)
!    
! compute for each Vert-column: grad(V)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      implicit none
      integer :: nx, ny
      real    :: V(nx, ny), dlam1(ny), dlat
      real    :: Divjp(ny), Divjm(ny)
      real    :: Vx(nx, ny), Vy(nx, ny), Vyd(nx, ny)      
      integer :: i, j
      do i=2, nx-1 
        Vx(i,:) = dlam1(:)*(V(i+1,:)-V(i-1,:))
      enddo
      Vx(1,:)  =  dlam1(:)*(V(2,:)-V(nx,:))
      Vx(nx,:) =  dlam1(:)*(V(1,:)-V(nx-1,:))

      do j=2, ny-1
        Vy(:,j) = dlat*(V(:,j+1)-V(:, j-1))
      enddo 
       Vy(:, 1) = dlat*2.*(V(:,2)-V(:,1))
       Vy(:,ny) = dlat*2.*(V(:,ny)-V(:,ny-1))
!~~~~~~~~~~~~~~~~~~~~
! 1/cos*d(vcos)/dy
!~~~~~~~~~~~~~~~~~~~~
      do j=2, ny-1      
        Vyd(:,j) = divJP(j)*V(:,j+1)-V(:, j-1)*divJM(j)
      enddo
      Vyd(:, 1) =  Vyd(:,2) 
      Vyd(:,ny) = Vyd(:,ny-1) 

      end  subroutine get_xyd_wind

      subroutine trig3d_fjets( nx, ny, nz, U, V, T, Q, P3D, PS, delp, delz, lon, lat, pmid, trig3d_fgf)
      implicit none
      integer ::  nx, ny, nz
      real    ::  lon(nx), lat(ny)
!      
      real, dimension(nz)          :: pmid
      real, dimension(nx, ny, nz)  :: U, V, T, Q, delp, delz, p3d
      real, dimension(nx, ny    )  :: PS                
      real, dimension(nx, ny, nz)  :: trig3d_fgf
!
! locals
!      
      real, dimension(nx, ny)  ::   ux, uy, uyd, vy, vx, vyd, ptx, pty     
      integer :: k, i, j
      
      real, parameter :: cappa=2./7., pref=1.e5
      real, dimension(nx, ny)  :: pt, w1, w2      
      
      real    :: rlon(nx), rlat(ny) , cosv(ny), tanlat(ny)  
      real    :: rlatc(ny-1), brcos(ny), brcos2(ny)              

      real    ::  dx, dy, dlat
       real    :: dlam1(ny), dlam2(ny),  divJp(ny), divJm(ny)  
       
            
      call subs_diag_geo(nx, ny,  lat, lon, rlat, rlon, dy, dx,  &
           cosv, rlatc, brcos, brcos2, dlam1, dlam2, dlat, divJp, divJm)

      do k=1, nz 
        w1(:,:)  = P3d(:,:,k)
        w2(:,:)  = T(:,:,k)

        pt = w2*(pref/w1)**cappa
        call get_xy_pt(Pt, ptx, pty, nx, ny, dlam1, dlat)
        w1(:,:) = V(:,:, K)
        call get_xyd_wind( w1, Vx, Vy, Vyd, nx, ny, dlam1, dlat, divJp, divJm)
        w1(:,:) = U(:,:, K)
        call get_xyd_wind( w1, Ux, Uy, Uyd, nx, ny, dlam1, dlat, divJp, divJm)

        trig3d_fgf(:,:,k) =   -ptx*ptx*ux - pty*pty*vy -(vx+uyd)*ptx*pty

      enddo
      end  subroutine trig3d_fjets
      
      subroutine trig3d_okubo( nx, ny, nz, U, V, T, Q, P3d, PS, delp, delz, lon, lat, pmid, trig3d_okw)
      implicit none
      integer ::  nx, ny, nz
      real    ::  lon(nx), lat(ny)
!      
      real, dimension(nz)          :: pmid
      real, dimension(nx, ny, nz)  :: U, V, T, Q, delp, delz, p3d
      real, dimension(nx, ny    )  :: PS                
      real, dimension(nx, ny, nz)  :: trig3d_okw
!
! locals
!      
      real, dimension(nx, ny)  ::   ux, uy, uyd, vy, vx, vyd, ptx, pty     
      integer :: k, i, j
      
      real, parameter :: cappa=2./7., pref=1.e5
      real, dimension(nx, ny)  :: pt, w1, w2, d1     
      
      real    :: rlon(nx), rlat(ny) , cosv(ny), tanlat(ny)  
      real    :: rlatc(ny-1), brcos(ny), brcos2(ny)              

      real    ::  dx, dy, dlat
      real    :: dlam1(ny), dlam2(ny),  divJp(ny), divJm(ny) 
             
      call subs_diag_geo(nx, ny,  lat, lon, rlat, rlon, dy, dx,  &
           cosv, rlatc, brcos, brcos2, dlam1, dlam2, dlat, divJp, divJm)

      do k=1, nz 
        w1(:,:)  = P3d(:,:,k)
        w2(:,:)  = T(:,:,k)

        pt = w2*(pref/w1)**cappa
        call get_xy_pt(Pt, ptx, pty, nx, ny, dlam1, dlat)
        w1(:,:) = V(:,:, K)
        call get_xyd_wind( w1, Vx, Vy, Vyd, nx, ny, dlam1, dlat, divJp, divJm)
        w1(:,:) = U(:,:, K)
        call get_xyd_wind( w1, Ux, Uy, Uyd, nx, ny, dlam1, dlat, divJp, divJm)

        trig3d_okw(:,:,k) =   -ptx*ptx*ux - pty*pty*vy -(vx+uyd)*ptx*pty
        w1 = (Ux -Vy)*(Ux-Vy) + (Vx +Uy)*(Vx+Uy)      ! S2
        W2 =  (Vx - Uyd)*(Vx - Uyd)
        D1 =  Ux + Vyd
     trig3d_okw(:,:,k) = W1 -W2
!    trig3d_okw(:, :, k)  =S2 -W2
!    trig3d_okw(:, :, k)  =D1*D1 + 4*(Vx*Uyd -Ux*Vyd)                         !  ocean
!    trig3d_okw(:, :, k) = trig3d_okw(:,:,k) + D1*D1 + 2.*D1*sqrt(abs(W1-W2)) !  S2 =W1Ted-luk	
      enddo   
      end  subroutine trig3d_okubo     
!      
      subroutine trig3d_dconv(nx, ny, nz, U, V, T, Q, P3d, PS, delp, delz, lon, lat, pmid, trig3d_conv, &
                              dcheat3d, precip2d, cld_klevs2d, scheat3d)

      implicit none
      integer ::  nx, ny, nz
      real    ::  lon(nx), lat(ny)
!      
      real, dimension(nz)          :: pmid
      real, dimension(nx, ny, nz)  :: U, V, T, Q, delp, delz, p3d
      real, dimension(nx, ny    )  :: PS                
      real, dimension(nx, ny, nz)  :: trig3d_conv
      
      real, dimension(nx, ny, nz)  :: dcheat3d, scheat3d                       
      real, dimension(nx, ny    )  :: precip2d
      integer,dimension(nx, ny, 3 ):: cld_klevs2d
      integer :: k
      end subroutine trig3d_dconv

      subroutine cires_3d_triggers( nx, ny, nz, lon, lat, pmid,          &
        U, V, W, T, Q, delp, delz, p3d, PS, HS, Hyam, Hybm, Hyai, Hybi,  &
        trig3d_okw, trig3d_fgf, trig3d_conv,                             &
        dcheat3d, precip2d, cld_klevs2d, scheat3d) 

      implicit none
      integer ::  nx, ny, nz
      real    ::  lon(nx), lat(ny)
!
! reversed ???  Hyai, Hybi , pmid
!      
      real, dimension(nz+2)        :: Hyai, Hybi
      real, dimension(nz+1)        :: Hyam, Hybm
!      
      real, dimension(nz)          :: pmid
      real, dimension(nx, ny, nz)  :: U, V, W, T, Q, delp, delz, p3d
      real, dimension(nx, ny    )  :: PS, HS                 
      real, dimension(nx, ny, nz)  :: trig3d_okw, trig3d_fgf, trig3d_conv
      real, dimension(nx, ny, nz)  :: dcheat3d, scheat3d                       
      real, dimension(nx, ny    )  :: precip2d
      integer,dimension(nx, ny, 3 ):: cld_klevs2d
      real    :: dzkm, zkm
      integer :: k
!==================================================================================      
! fgf and OW-triggers
! read PRECIP + SH/DC conv heating  + cloud-top-bot-middle from "separate" file !!!
!
!===================================================================================  

      call trig3d_fjets( nx, ny, nz, U, V, T, Q, P3D, PS, delp, delz, lon, lat, pmid, trig3d_fgf) 
      call trig3d_okubo( nx, ny, nz, U, V, T, Q, P3D, PS, delp, delz, lon, lat, pmid, trig3d_okw)
      call trig3d_dconv(nx, ny, nz, U, V, T, Q,   P3D, PS, delp, delz, lon, lat, pmid, trig3d_conv, &
         dcheat3d, precip2d, cld_klevs2d, scheat3d)  
!=====================================================================================================
! output of triggers: trig3d_fgf, trig3d_okw, trig3d_conv, cheat3d, precip2d, cld_klevs2d, scheat3d
!
! Bulk momentum flux=/ 0 and levels for launches
!
!=====================================================================================================	                
 111  format(i6, 4(3x, F8.3),  '  trigger-grid ')     
 
      do k=1, nz-1
       zkm  = -7.*alog(pmid(k)*1.e-3)
       dzkm = zkm +7.*alog(pmid(k+1)*1.e-3)
       write(6,111) k, hybi(k), pmid(k), zkm, dzkm !' triggers '
      enddo
      
      end subroutine cires_3d_triggers
!==================================================================================      
!                              tot-flux  launch  0 or 1 # of Launches 
! specify time-dep bulk sources:   taub, klev,  if_src, nf_src
!
!==================================================================================
      subroutine get_spectra_tau_convgw &
      (nw, im, levs, dcheat,  scheat, precip, icld, xlatd, sinlat, coslat,taub, klev, if_src, nf_src)
!
! temporarily can put GEOS-5/MERRA-2 GW-lat dependent function
!      
      integer                   :: nw, im, levs 
      integer,dimension(im,3)   :: icld   
      real, dimension(im, levs) :: dcheat,  scheat   
      real, dimension(im)       :: precip, xlatd, sinlat, coslat      
      real, dimension(im)       :: taub
      integer, dimension(im)    :: klev, if_src
      integer ::  nf_src
!
! locals
      real, parameter           ::  precip_max = 100.   ! mm/day
      real, parameter           ::  tau_amp    = 35.e-3  ! 35 mPa   
       
      integer  :: i, k, klow, ktop, kmid
      real     :: dtot, dmax, daver
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
!  100 mb launch  and MERRA-2 slat-forcing
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
      real, dimension(im, levs)    :: trig_fgf
!      real, dimension(im, levs+1)  :: pint      
      real, dimension(im)          :: xlatd, sinlat, coslat        
      real, dimension(im)          :: taub
      integer, dimension(im)          :: klev, if_src
      integer ::  nf_src
! locals
      real, parameter           ::  tlim_fgf = 100.     ! trig_fgf > tlim_fgf, launch waves should scale-dependent
      real, parameter           ::  tau_amp   = 35.e-3  ! 35 mPa
      real, parameter           ::  pmax = 750.e2,  pmin = 100.e2  
      integer, parameter        ::  klow =127-92,   ktop=127-45
      integer, parameter        ::  kwidth = ktop-klow+1   
      integer  :: i, k,  kex
      real     :: dtot, dmax, daver      
      real     :: fnorm, tau_min     
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
        klev(i) = 127-45     
      enddo
!                 
      end subroutine get_spectra_tau_nstgw   
! 
      subroutine get_spectra_tau_okw(nw, im, levs,  trig_okw, xlatd, sinlat, coslat, taub, klev, if_src, nf_src)
      integer                      :: nw, im, levs 
      real, dimension(im, levs)    :: trig_okw
!      real, dimension(im, levs+1)  :: pint    
      real, dimension(im)          :: xlatd, sinlat, coslat     
      real, dimension(im)          :: taub
      integer, dimension(im)       :: klev, if_src
      integer ::  nf_src
! locals
      real, parameter           ::  tlim_okw = 100.     ! trig_fgf > tlim_fgf, launch waves should scale-dependent
      real, parameter           ::  tau_amp   = 35.e-3  ! 35 mPa  
      real, parameter           ::  pmax = 750.e2,  pmin = 100.e2  
      integer, parameter        ::  klow =127-92,   ktop=127-45
      integer, parameter        ::  kwidth = ktop-klow+1   
      integer  :: i, k,  kex
      real     :: dtot, dmax, daver      
      real     :: fnorm, tau_min       
      
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
!
!
!      
!>\ingroup cires_ugwp_run
!> @{
!!
!!
      subroutine slat_geos5_tamp(im, tau_amp, xlatdeg, tau_gw)
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
      end subroutine slat_geos5_tamp
      
      subroutine slat_geos5(im, xlatdeg, tau_gw)
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
      end subroutine slat_geos5
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

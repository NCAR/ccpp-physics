  subroutine ugwp_wmsdis_naz(levs, ksrc, nw, naz, kxw, taub_lat, ch, xaz, yaz,    &
              fcor, c2f2, dp, zmid, zint, pmid, pint, rho, ui, vi, ti,            &
              kvg, ktg, krad, kion, bn2i, bvi, rhoi, ax, ay, eps, ked, tau1)
!
!
!   use para_taub,  only         :  tau_ex
   use    ugwp_common, only      : rcpd, grav, rgrav     
   implicit none 
!   
   integer                       ::  levs
   integer                       ::  nw, naz       ! # - waves for each azimuth (naz)
   integer                       ::  ksrc          ! source level
   real                          ::  kxw           ! horizontal wn
   real                          ::  taub_lat      ! lat-dep tau_bulk N/m2      
!
   real, dimension(nw)           ::   ch, dch, taub_spect
   real, dimension(naz)          ::   xaz, yaz
   real, dimension(levs+1)       ::   ui, vi, ti, bn2i,  bvi, rhoi, zint, pint
   real, dimension(levs  )       ::   dp, rho, pmid, zmid
   real                          ::   fcor, c2f2
   real, dimension(levs+1)       ::   kvg, ktg,  kion, krad, kmol 
   
! output/locals        
   real, dimension(levs  )       ::  ax, ay, eps
   real, dimension(levs+1)       ::  ked , tau1
   real, dimension(levs+1  )     ::  uaz
   
   real, dimension(levs,  naz  ) ::  epsd
   real, dimension(levs+1, naz ) ::  atau, kedd
   
   real, dimension(levs+1      ) ::  taux, tauy, bnrho  
   real, dimension(levs  )       ::  dzirho , dzpi 

!  
   integer     :: iaz,  k , inc
   real, parameter               ::   gcstar=1.0
   integer , parameter           ::   nslope=1
   real                          ::   spnorm        ! source level normalization factor for the Broad Spectra
   real                          ::   bnrhos        ! sum(taub_spect*dc) = spnorm   taub_sect_norm  = taub_spect/spnorm
!     
   atau=0.0 ;  epsd=0.0 ; kedd=0.0  
   bnrhos = bvi(ksrc)/rhoi(ksrc)
   do k=1,levs
     dzpi(k)   = zint(k+1)-zint(k)   
     dzirho(k) = 1.0 / (rho(k)*dzpi(k))                 !   grav/abs(dp(k)) still hydrostatic "ugwp"
     bnrho(k)  = (rhoi(k)/bvi(k)) !*bnrhos * gcstar     !   gcstar=1.0 and   bnrho(k=ksrc) =1. 
   enddo
   k = levs+1
   bnrho(k)  = (rhoi(k)/bvi(k))*bnrhos    
!
! re-define    ch,  dch, taub_spect,  this portion can be moved to "ugwp_init"
!              
!  
!   
   call FVS93_ugwps(nw,  ch,  dch, taub_spect, spnorm,  nslope, bn2i(ksrc), bvi(ksrc), bnrho(ksrc)) 
    
    
!   print *, ' after FVS93_ugwp ', nw,  maxval(ch), minval(ch)  
!
!    do normaalization for the spectral element of the saturated flux
!
   bnrho = bnrho *spnorm

!  print *	
!  do inc=1, nw	
!    write(6,221) inc, ch(INC),taub_lat*taub_spect(inc),  spnorm, dch(inc)
!221 FORMAT( i6, 2x, F8.2, 3(2x, E10.3))
!  enddo
!  pause

   loop_iaz:  do iaz =1, naz
       
                do k=1,levs+1
                  uaz(k) =ui(k)*xaz(iaz) +vi(k)*yaz(iaz)  
                enddo
!
!
! multi-wave broad spectrum of FVS-93 with ~scheme of WMS-IFS 2010
!
!       print *, ' iaz before ugwp_wmsdis_az1 ', iaz
!
       
       call  ugwp_wmsdis_az1(levs, ksrc, nw, kxw, ch, dch, taub_spect, taub_lat,                &
          spnorm, fcor, c2f2, zmid, zint, rho, uaz, ti, bn2i, bvi, rhoi, bnrho, dzirho, dzpi,   &
          kvg, ktg, krad, kion, kmol, epsd(:, iaz), kedd(:,iaz), atau(:, iaz) )

!       print *, ' iaz after ugwp_wmsdis_az1 ', iaz        

!                               
   enddo loop_iaz                    ! azimuth of gw propagation directions
!
!     sum over azimuth and project atau(z, iza) =>(taux and tauy)
!         for scalars              for "wave-drag vector"
! 
      eps =0. ; ked =0.
      do k=ksrc, levs 
         eps(k) = sum(epsd(k,:))*rcpd  
      enddo  
 
      do k=ksrc, levs+1
        taux(k) = sum( atau(k,:)*xaz(:))
        tauy(k) = sum( atau(k,:)*yaz(:))
        ked(k)  = sum( kedd(k,:))
      enddo   
!
      tau1(ksrc:levs) = taux(ksrc:levs)         
      tau1(1:ksrc-1) = tau1(ksrc)

! end  solver: gw_azimuth_solver_ls81    
!      sign ax in rho*du/dt  = -d(rho*tau)/dz
!                             [(k) - (k+1)]
!      du/dt = ax = -1/rho*d( tau) /dz
!
      ax =0. ; ay = 0.
  
      do k=ksrc, levs
        ax(k) = dzirho(k)*(taux(k)-taux(k+1))
        ay(k) = dzirho(k)*(tauy(k)-tauy(k+1))
      enddo  
      call ugwp_limit_1d(ax, ay, eps, ked, levs) 

      return 
      end subroutine ugwp_wmsdis_naz


! =======================================================================
   subroutine ugwp_wmsdis_az1(levs, ksrc, nw, kxw, ch, dch, taub_sp, tau_bulk,  &
          spnorm, fcor, c2f2, zm, zi, rho, um, tm, bn2, bn, rhoi,  bnrho,       &
          dzirho, dzpi, kvg, ktg, krad, kion, kmol, eps, ked, tau ) 
!	  
!      use para_taub,  only        :  tau_ex, xlatdeg !for exchange src-tau
!
      use cires_ugwp_module, only : f_coriol, f_nonhyd, f_kds, linsat                 
      use cires_ugwp_module, only : ipr_ktgw, ipr_spgw, ipr_turb, ipr_mol
      use cires_ugwp_module, only : rhp4, rhp2, rhp1, khp, cd_ulim 
! =======================================================================
      integer  ::  levs, ksrc, nw
      real     ::  fcor, c2f2, kxw
! 
      real, dimension(nw)         :: taub_sp, ch, dch
      real                        :: tau_bulk, spnorm   
      real, dimension(levs)       :: zm, rho, dzirho, dzpi
      real, dimension(levs+1)     :: zi, um, tm, bn2, bn, rhoi, bnrho    
      real, dimension(levs+1)     :: kvg, ktg, krad, kion, kmol
      real, dimension(levs+1)     :: ked, tau
      real, dimension(levs  )     :: eps
!
!locals
      integer                     ::  k, inc
      real, dimension(levs+1)     ::  umi
      real                        ::  zcin, zci_min, ztmp, zcinc
      real                        ::  zcimin=0.5                ! crit-level precision, 0.5 and start of Ch_MIN
      real,  parameter            ::  Keff = 0.2

      real, dimension(nw)         ::  zflux                                       !
      real, dimension(nw)         ::  wzact, zacc                                 ! =1  ..crit level change it

      real, dimension(levs)       ::  zcrt                                        !
      real, dimension(nw, levs)   ::  zflux_z, zact   
      
      real                        ::  zdelp, kxw2 
      real                        ::  vu_eff, vu_lin, v_kzw, v_cdp, v_wdp, v_kzi
      real                        ::  dfsat,  fdis, fsat, fmode, expdis
      real                        ::  vc_zflx_mode, vm_zflx_mode
      real                        ::  tau_g5     
! =======================================================================
!eps, ked, tau 

           eps (:) =0;  ked = 0.0 ;
           kxw2    = kxw*kxw
!                 
           zcrt(1:levs)    = 0.0
           umi(1:levs+1) = um
!           umi(1:levs+1) = um(1:levs+1) -um(ksrc) 

           zci_min = zcimin
	   
!          CALL slat_geos5(1, xlatdeg(1), tau_g5) 	   
!	   tau_bulk = tau_g5                       !tau_bulk*0.75               !3.75e-2   
!	     
           zflux(:)  = taub_sp(:)*tau_bulk         ! includes tau_bulk(x,y) and spectral normalization
	   
           zflux_z(1:nw,ksrc)=zflux(:)                   
           
	   tau(1:levs+1) = tau_bulk            ! constant flux for all layers k<ksrc & sum(taub_sp*dc)
	   
!           print *, ksrc, tau_bulk*1.e3, nw, maxval(ch), minval(ch), sum(taub_sp*dch), ' 1.00= NORM-N bf z-loop in wmsdis '
!	   pause
!	 
           wzact(1:nw)      = 1.0
	   zact(1:nw,:)     = 1.0
        do k=ksrc+1, levs 
	   umi(k) = um(k) -um(ksrc) 
!                   
! first check for critical levels
!
	    
            zci_min=  max(zcimin,   abs(umi(k)))
	    zci_min = zci_min*sign(1.0, umi(k))              
! ----------------------------------------------          
! set zact to zero if critical level encountered
! ----------------------------------------------
           do inc=1, nw	       
	       if (wzact(inc) == 0.) cycle
               zcin=ch(inc)                                     !zci -transformed	       
               ztmp= 0.5 + sign(0.5,zcin-zci_min)
               zacc(inc)  = wzact(inc)-ztmp
               wzact(inc) = ztmp                                 !  o -critical layer absorption
	       if ((zcin-zci_min) .le. 0.0)  zact(inc, k)=0      ! flag for critical layer
           enddo
	   
               tau(k) = sum(zact(:,k)*zflux(:)*dch(:))          ! critical "band"-zact filtering

!
! integrate to get critical-level contribution to mom deposition on this level, i.e. zacc=1
!
! --------------------------------------------------------
! get weighted average of phase speed in layer/azimuth
!  diagnostics to help  "abrupt"; do not apply "zcrt" in V1
! ----------------------------------------------------------
!            if(tau(k)>0.0 ) then
!                 ztmp = sum( ch(:)*zacc(:)*zflux(:)*dch(:) )
!                 zcrt(k)=ztmp/tau(k)
!             else
!                  zcrt( k )=zcrt(k-1)
!            endif
! ---------------------------------------------------------
! do saturation (eq. (26) and (27) of scinocca 2003)
!  +  add molecular/eddy dissipation od gw-spectra vay-2015
!      for each mode & direction
!      x by exp(-mi*zdelp)  x introduce ....... mi(nw)
!
! mode-loop +  add molecular/eddy dissipation od gw-spectra vay-2015
!
            do inc=1,nw
	       if (zact(inc,k) == 0.0) then
		zflux(inc)      =  0.0
		zflux_z(inc,k)  = zflux(inc) 
            else            
	    vu_eff  =  kvg(k)                  ! + ktg (k)        !* ipr_ktgw                
            vu_lin  =  kion(k)                 ! + krad(k)        !* ipr_ktgw   	    
            vu_eff = 2.e-5*exp(zi(k)/7000.)+.01
                zcin= ch(inc)
               
!=======================================================================
! saturated limit    wfit = kzw*kzw*kt; wfdt = wfit/(kxw*cx)*betat
! & dissipative      kzi = 2.*kzw*(wfdm+wfdt)*dzpi(k)
!           define   kxw = 
!=======================================================================  
               v_cdp =  zcin-umi(k) 
               v_wdp = kxw*v_cdp
               if (v_wdp.gt.0) then 
                  v_kzw =  bn(k)/v_cdp                        !can be non-hydrostatic 
                  v_kzi  = abs(( v_kzw*v_kzw*vu_eff + vu_lin) /v_wdp*v_kzw)
                  expdis = exp(-2.*v_kzi*dzpi(k) )
               else
                  v_kzi  = 0.
                  expdis = 1.0  
               endif
               fmode =  zflux(inc) 
               fdis  =  fmode*expdis                             ! only dissipation/crit_lev degrades it
!------------------------
! includes rho/bn /(rhos/bns) *spnorm
!------------------------
               fsat  = bnrho(k)* v_cdp*v_cdp /zcin               !  expression for saturated flux <u'w'> 
                                                                 !  zfluxs=gcstar*zfct(  k)*(zcin-zui(  k ))**2/zcin
! flux_tot - sat.flux
!
                   dfsat= fdis-fsat
                   if( dfsat > 0.0 ) then
! put sat-n limit
                      zflux(inc)    = fsat
                   else
! assign dis-ve flux
                      zflux(inc)    =fdis
                   endif  
                     zflux_z(inc,k)=zflux(inc)
                  
                    if (zflux_z(inc,k) > zflux_z(inc,k-1) ) zflux_z(inc,k) = zflux_z(inc,k-1)
                 
	       endif
		   
             enddo                 
!
! integrate over spectral modes  zpu(y, z, azimuth)    zact(  inc, )*zflux(  inc, )*[d("zcinc")]
!
                    tau(k)  = sum( zflux_z(:,k)*dch(:))
!------------------------------------------------------------------------------
! define expressions for eps-heat + Ked,  needs more work for the broad spectra
!                                         formulation especially for Ked
!              after defining Ked .....GW-eddy cooling needs to be added
!                                      for now "only" heating here
!==============================================================================
                    eps(k) =0.     
               do inc=1, nw
	       if (zact(inc,k) == 0.0) cycle                                           ! dc-integration  + dtau/dz
                    vc_zflx_mode =  zflux(inc)
                 
                    zdelp= abs(ch(inc)-umi(k)) * dch(inc) /dzpi(k)
                    vm_zflx_mode=zflux_z(inc,k-1)
                    eps(k) =eps( k ) + (vm_zflx_mode-vc_zflx_mode)*zdelp     !  heating >0                        
                               
       
                enddo                               !inc=1, nw
                ked(k) = Keff*eps(k)/bn2(k)
!
! --------------
!
           enddo                              ! end k do-loop vertical loop        do k=ksrc+1, levs 

!top lid	   
	   k =levs+1
	   ked(k) = ked(k-1)
!	   eps(k) = eps(k-1)
	   tau(k) =tau(k-1)*0.933

! from surface to ksrc-1   
!	   tau(1:ksrc) = tau(ksrc)
	   ked(1:ksrc)  = 0.
	   eps(	1:ksrc) = 0.    	    
	   
!      
!  output:    eps, ked, tau   for given azimuth
!    
      end subroutine ugwp_wmsdis_az1
!
!
      subroutine FVS93_ugwps(nw,  ch,  dch, taub_sp, spnorm, nslope, bn2, bn, bnrhos)
      implicit none
      integer :: nw, nslope
      real    :: bn2, bn, bnrhos
!!      real    :: taub_lat                                                ! bulk - lat-dep momentum flux
      real, dimension (nw) :: ch,  dch, taub_sp
! locals
      integer :: i, inc
      real, parameter :: zcimin = 0.5, zcimax = 95.0, zgam =1./4.
      real, parameter :: zms    = 6.28e-3/2.                             ! mstar  Lz ~ 2km
      real            :: zxran, zxmax, zxmin, zx1, zx2, zdx, ztx, rch
      real            :: bn3, bn4, zcin, tn4, tn3, tn2, cstar
      real            :: spnorm                                          ! needs to be passed for saturation flux norm-n
      real            :: tau_bulk  
!--------------------------------------------------------------------
!
!     transforms ch -uniform =>  1/ch and back to non-uniform ch, dch
!
!-------------------------------------------------------------------     
! note that this is expresed in terms of the intrinsic ch or vertical wn=N/cd
! at launch cd=ch-um(ksrc), the transformation is identical for all
!    levels, azimuths and horizontal pixels
! see eq. 28-30 of scinocca 2003.   x = 1/c stretching transform
!
          zxmax=1.0 /zcimin
          zxmin=1.0 /zcimax
          zxran=zxmax-zxmin
          zdx=zxran/float(nw-1)     ! d_kz or d_mi
! 
!
          zx1=zxran/(exp(zxran/zgam)-1.0 )           !zgam =1./4.
          zx2=zxmin-zx1
!
!                            add idl computations for zci =1/zx
!                            x = 1/c stretching transform to look at final ch(i), dch(i)
! 

        do i=1, nw
          ztx=float(i-1)*zdx+zxmin
          rch=zx1*exp((ztx-zxmin)/zgam)+zx2                                 !eq. 29 of scinocca 2003
          ch(i)=1.0 /rch                                                    !eq. 28 of scinocca 2003
          dch(i)=ch(i)*ch(i)*(zx1/zgam)*exp((ztx-zxmin)/zgam)*zdx           !eq. 30 of scinocca 2003
        enddo   
!
! nslope-dependent flux taub_spect(nw)  momentum flux spectral density
!                  need to check math....expressions
! eq. (25) of scinocca 2003 with u-uo=0 it is identical to all azimuths
!
!
            cstar=bn/zms
            bn4=bn2*bn2  ! four times
            bn3=bn2*bn 
   if(nslope==1) then
! s=1 case
      do inc=1, nw 
           zcin=ch(inc)
           tn4=(zms*zcin)**4
           taub_sp(inc) =bnrhos  * zcin*bn4/(bn4+tn4)
      enddo
!
   elseif(nslope==2) then
! s=2 case
      do inc=1, nw 
            zcin=ch(inc)
            tn4=(zms*zcin)**4           
        taub_sp(inc)= bnrhos*zcin*bn4/(bn4+tn4*zcin/cstar)
       enddo
!
    elseif(nslope==-1) then
! s=-1 case
       do inc=1, nw 
         zcin=ch(inc)
         tn2=(zms*zcin)**2
        taub_sp(inc)=bnrhos*zcin*bn2/(bn2+tn2)
       enddo
! s=0 case
    elseif(nslope==0) then

       do inc=1, nw 
           zcin=ch(inc)
           tn3=(zms*zcin)**3
        taub_sp(inc)=bnrhos*zcin*bn3/(bn3+tn3)
       enddo
     endif                                     ! for n-slopes
!=============================================
!       normalize launch momentum flux
!    ------------------------------------
! (rho x f^h = rho_o x f_p^total)   integrate (zflux x dx)

        tau_bulk= sum(taub_sp(:)*dch(:))
        spnorm= 1./tau_bulk
       
         do inc=1, nw 
          taub_sp(inc)=spnorm*taub_sp(inc)
         enddo

      end  subroutine FVS93_ugwps


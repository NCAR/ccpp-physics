module cires_ugwp_solv2_v1_mod


contains


!---------------------------------------------------
!  Broad spectrum FVS-1993, mkz^nSlope with nSlope = 0, 1,2 
!  dissipative solver with NonHyd/ROT-effects
!  reflected GWs treated as waves with "negligible" flux,
!  they are out of given column
!---------------------------------------------------
      subroutine cires_ugwp_solv2_v1(im, levs, dtp    ,          &
                 tm , um, vm, qm, prsl, prsi, zmet,  zmeti,      &
                 prslk, xlatd, sinlat, coslat,                   &
                 grav, cpd, rd, rv, omega, pi, fv,               &
                 pdudt, pdvdt, pdtdt, dked, tauabs, wrms, trms,  &
                 tau_ngw, mpi_id, master, kdt)
!
!--------------------------------------------------------------------------------
!      nov 2015 alternative gw-solver for nggps-wam
!      nov 2017 nh/rotational gw-modes for nh-fv3gfs
!      oct 2019 adding empirical satellite-based
!          source function and *F90 CIRES-style of the code
! --------------------------------------------------------------------------------
!  
 
      use machine,          only : kind_phys  

      use cires_ugwp_module_v1,only :  krad, kvg, kion, ktg 
      
      use cires_ugwp_module_v1,only :  knob_ugwp_doheat, knob_ugwp_dokdis, idebug_gwrms
      
      use ugwp_common_v1 ,     only : dw2min, velmin, hpscale, rhp, rh4
!
      use ugwp_wmsdis_init_v1, only : v_kxw,  rv_kxw,   v_kxw2, tamp_mpa, tau_min, ucrit, &    
                                   maxdudt, gw_eff,  dked_min,  dked_max, maxdtdt,     &
                                   nslope,  ilaunch, zms,                              &
                                   zci,     zdci,    zci4, zci3, zci2,                 &
                                   zaz_fct, zcosang, zsinang,  nwav,    nazd,          &
                                   zcimin, zcimax, rimin, sc2, sc2u, ric
! 
      implicit none
!23456 

      integer, intent(in)  :: levs              ! vertical level
      integer, intent(in)  :: im              ! horiz tiles

      real ,intent(in)   :: dtp               ! model time step       
      real ,intent(in)   :: vm(im,levs)       ! meridional wind 
      real ,intent(in)   :: um(im,levs)       ! zonal wind
      real ,intent(in)   :: qm(im,levs)       ! spec. humidity
      real ,intent(in)   :: tm(im,levs)       ! kinetic temperature 

      real ,intent(in)   :: prsl(im,levs)     ! mid-layer pressure
      real ,intent(in)   :: prslk(im,levs)    ! mid-layer exner function     
      real ,intent(in)   :: zmet(im,levs)     ! meters now !!!!!       phil =philg/grav
      real ,intent(in)   :: prsi(im,levs+1)   !  interface pressure
      real ,intent(in)   :: zmeti(im,levs+1)  !  interface geopi/meters      
      real ,intent(in)   :: xlatd(im)         ! lat was in radians, now with xlat_d in degrees
      real ,intent(in)   :: sinlat(im)
      real ,intent(in)   :: coslat(im)
      real ,intent(in)   :: tau_ngw(im)

      integer, intent(in):: mpi_id, master, kdt

      real ,intent(in)   :: grav, cpd, rd, rv, omega, pi, fv
!  
!
! out-gw effects
!
      real ,intent(out) :: pdudt(im,levs)     ! zonal momentum tendency
      real ,intent(out) :: pdvdt(im,levs)     ! meridional momentum tendency
      real ,intent(out) :: pdtdt(im,levs)     ! gw-heating (u*ax+v*ay)/cp
      real ,intent(out) :: dked(im,levs)      ! gw-eddy diffusion
!
! GW diagnostics => next move it to "module_gw_diag"
!      
      real ,intent(out) :: tauabs(im,levs)    !
      real ,intent(out) :: wrms(im,levs)      ! 
      real ,intent(out) :: trms(im,levs)      !   
                                          	    
      real              :: zwrms(nwav,nazd), wrk1(levs), wrk2(levs)
      real              :: atrms(nazd, levs),awrms(nazd, levs), akzw(nwav,nazd, levs+1)                   
!
! local ===========================================================================================  
      real              :: taux(levs+1)         ! EW component of vertical momentum flux (pa)
      real              :: tauy(levs+1)         ! NS component of vertical momentum flux (pa)
      real              :: fpu(nazd, levs+1)    ! az-momentum flux
      real              :: ui(nazd, levs+1)     ! azimuthal wind
           
      real              :: fden_bn(levs+1)         ! density/brent
      real              :: flux_z(nwav,levs+1)   
      real              :: flux(nwav, nazd)                   	
!
! ===============================================================================================
!                ilaunch:levs ....... MOORTHI's improvements
!      all computations of GW-effects include interface layers from ilaunch+1 to levs +1
!      at k=levs+1, extrapolation of MF-state has been made, "ideally" all spectral modes should
!      be absorbed; 2-options for this "ideal" requirement
!      a) properly truncate GW-spectra ; b) dissipate all GW-energy in the top layers ( GW-sponge)
!=====================================================================================================    
!   	 
       real  :: bn(levs+1)                  ! interface BV-frequency 
       real  :: bn2(levs+1)                 ! interface BV*BV-frequency 
       real  :: rhoint(levs+1)               ! interface density 
       real  :: uint(levs+1)                 ! interface zonal wind
       real  :: vint(levs+1)                 ! meridional wind
       
       real  :: irhodz_mid(levs), dzdt(levs+1), bnk(levs+1), rhobnk(levs+1)       
      
       real  :: v_zmet(levs+1)
       real  :: vueff(levs+1)     
       real  :: dfdz_v(nazd, levs)           ! axj = -df*rho/dz       directional momentum deposition
       

       real  :: suprf(levs+1)                      ! RF-super linear dissipation  
 
       real, dimension(levs)   ::  atm , aum, avm, aqm, aprsl, azmet
       real, dimension(levs+1) ::  aprsi, azmeti 

       real  :: wrk3(levs)    
       real, dimension(levs) :: uold, vold, told, unew, vnew, tnew
       real, dimension(levs) :: dktur, rho, rhomid, adif, cdif 

       real  :: rdci(nwav),  rci(nwav)
       real  :: wave_act(nwav, nazd)           ! active waves at given vert-level
       real  :: ul(nazd)                       ! velocity in azimuthal direction at launch level
       real  :: bvi, bvi2, bvi3, bvi4, rcms    ! BV at launch level	 
       real  :: c2f2, cf1   
           

       real  ::  flux_norm                      ! norm-factor
       real  ::  taub_src, rho_src
!
! scalars
!
       real  :: zthm, dtau, cgz, ucrit_maxdc     
       real  :: vm_zflx_mode, vc_zflx_mode
       real  :: kzw2, kzw3, kdsat, cdf2, cdf1, wdop2,v_cdp2
       real  :: ucrit_max
       real  :: pwrms, ptrms
       real  :: zu, zcin, zcin2, zcin3, zcin4, zcinc
       real  :: zatmp, fluxs, zdep,  ze1, ze2
!
       real  :: rcpdl, grav2cpd, rcpd, rcpd2, pi2, rad_to_deg
       real  :: deg_to_rad, rdi, gor, grcp, gocp, bnv2min, bnv2max, gr2
       real  :: grav2, rgrav, rgrav2, mkzmin, mkz2min
!  
       real  :: zdelp, zdelm, taud_min
       real  :: tvc,  tvm, ptc, ptm
       real  :: umfp, umfm, umfc, ucrit3   
       real  :: fmode, expdis, fdis
       real  :: v_kzi, v_kzw, v_cdp, v_wdp, tx1, fcorsat, dzcrit
       real  :: v_wdi, v_wdpc
       real  :: ugw, vgw, ek1, ek2, rdtp, rdtp2 
      
       integer :: j, jj, k, kk, inc, jk, jkp, jl, iaz
       integer ::  ksrc, km2, km1, kp1, ktop      
!
! Kturb-part
!  
           
      real     :: uz, vz, shr2 , ritur, ktur      
            
      real     :: kamp, zmetk, zgrow
      real     :: stab, stab_dt, dtstab
      integer  :: nstab, ist, anstab(levs)   
      real     :: w1, w2, w3, dtdif
 
      real     ::  dzmetm, dzmetp, dzmetf, bdif, kturp   
      real     ::  bnrh_src   
!--------------------------------------------------------------------------
!

        if (mpi_id == master .and. kdt < 2) then
              print *, im, levs, dtp, kdt,    ' vay-solv2-v1'
              print *,  minval(tm), maxval(tm), ' min-max-tm '
	      print *,  minval(vm), maxval(vm), ' min-max-vm '
	      print *,  minval(um), maxval(um), ' min-max-um '
	      print *,  minval(qm), maxval(qm), ' min-max-qm '	
	      print *,  minval(prsl), maxval(prsl), ' min-max-Pmid '  
	      print *,  minval(prsi), maxval(prsi), ' min-max-Pint ' 
	      print *,  minval(zmet), maxval(zmet), ' min-max-Zmid '
	      print *,  minval(zmeti), maxval(zmeti), ' min-max-Zint ' 	  
	      print *,  minval(prslk), maxval(prslk), ' min-max-Exner ' 
	      print *,  minval(tau_ngw), maxval(tau_ngw), ' min-max-taungw '
	      print *, 	 tau_min,  ' tau_min ',   tamp_mpa, ' tamp_mpa '    	              	      	      
! 
	endif
        	
        if (idebug_gwrms == 1) then
	  tauabs=0.0; wrms =0.0 ; trms =0.0 
	endif
	
      
       grav2 = grav + grav
       rgrav = 1.0/grav
       rgrav2 = rgrav*rgrav
       rdi = 1.0/rd
       gor = grav/rd
       gr2 = grav*gor
       rcpd = 1.0/cpd
       rcpd2 = 0.5/cpd
       rcpdl = cpd*rgrav             ! 1/[g/cp]  == cp/g
       pi2 = 2.0*pi
       grcp = grav*rcpd
       gocp = grcp
       grav2cpd = grav*grcp          !  g*(g/cp)= g^2/cp
       rad_to_deg=180.0/pi
       deg_to_rad=pi/180.0
       bnv2min = (pi2/1800.)*(pi2/1800.)
       bnv2max = (pi2/30.)*(pi2/30.) 
       mkzmin = pi2/80.0e3
       mkz2min = mkzmin*mkzmin
 
       rci(:) = 1./zci(:)
       rdci(:) = 1./zdci(:)
          
       rdtp = 1./dtp
       rdtp2 = 0.5*rdtp 
!
! launch level control ksrc > 2
!       
   
       ksrc= max(ilaunch, 3)
       km2 = ksrc - 2       
       km1 = ksrc - 1
       kp1 = ksrc + 1
       ktop= levs+1

       do k=1,levs
	   suprf(k) = kion(k)    ! approximate 1-st order damping with Fast super-RF of FV3
            pdvdt(:,k) = 0.0
            pdudt(:,k) = 0.0
            pdtdt(:,k) = 0.0
            dked(: ,k) = 0.0
        enddo
	
!-----------------------------------------------------------	
! column-based j=1,im pjysics with 1D-arrays
!-----------------------------------------------------------
       DO j=1, im   
           
         jl =j
	 tx1           = 2*omega * sinlat(j) *rv_kxw        
	 cf1 = abs(tx1)
         c2f2      = tx1 * tx1
	 ucrit_max = max(ucrit, cf1)
	 ucrit3 = ucrit_max*ucrit_max*ucrit_max
!
! ngw-fluxes at all gridpoints (with tau_min at least)
!	 	 	 
         taub_src = max(tau_ngw(jl), tau_min)  
	 aum(km2:levs) = um(jl,km2:levs)
	 avm(km2:levs) = vm(jl,km2:levs)	 	 
	 atm(km2:levs) = tm(jl,km2:levs)	
	 aqm(km2:levs) = qm(jl,km2:levs)
	 aprsl(km2:levs) = prsl(jl,km2:levs)
	 azmet(km2:levs) =	 zmet(jl,km2:levs)
	 aprsi(km2:levs+1) =	 prsi(jl,km2:levs+1)
	 azmeti(km2:levs+1) =	 zmeti(jl,km2:levs+1)
	 
	 rho_src = aprsl(ksrc)*rdi/atm(ksrc)	 
	 	 		 
!       ---------------------------------------------  
!   interface mean flow parameters launch -> levs+1
!       ---------------------------------------------
       do jk= km1,levs
           tvc = atm(jk)   * (1. +fv*aqm(jk))
           tvm = atm(jk-1) * (1. +fv*aqm(jk-1))
	   ptc =  tvc/ prslk(jl, jk)
	   ptm =  tvm/prslk(jl,jk-1)	   
!
           zthm          = 2.0 / (tvc+tvm)
!	   
           uint(jk)   = 0.5 *(aum(jk-1)+aum(jk))
           vint(jk)   = 0.5 *(avm(jk-1)+avm(jk))
           rhomid(jk) = aprsl(jk)*rdi/atm(jk)
           rhoint(jk) = aprsi(jk)*rdi*zthm                  !  rho = p/(RTv) 
           zdelp          = azmeti(jk+1)-azmeti(jk)         !  >0 ...... dz-meters
           zdelm          = 1./(azmet(jk)-azmet(jk-1))      ! 1/dz  ...... 1/meters	  
	   dzdt(jk) = dtp/zdelp 
!	   
!          bvf2 = grav2 * zdelm * (ptc-ptm)/ (ptc + ptm) ! N2=[g/PT]*(dPT/dz)
!	   
           bn2(jk)    = grav2cpd*zthm  * (1.0+rcpdl*(tvc-tvm)*zdelm)
           bn2(jk)    = max(min(bn2(jk), bnv2max), bnv2min)
           bn(jk)     = sqrt(bn2(jk))   
	   bnk(jk) = bn(jk)*v_kxw
	   rhobnk(jk)=rhoint(jk)/bnk(jk)*v_kxw 
	   wrk3(jk)= 1./zdelp/rhomid(jk)                           ! 1/rho_mid(k)/[Z_int(k+1)-Z_int(k)]
           irhodz_mid(jk) = rdtp*zdelp*rhomid(jk)/rho_src 

           v_zmet(jk)  = 2.*zdelp                                  ! 2*kzi*[Z_int(k+1)-Z_int(k)]
! 	          
!
! diagnostics -Kzz above PBL
!
           uz = aum(jk) - aum(jk-1)
           vz = avm(jk) - avm(jk-1)
           shr2 = (max(uz*uz+vz*vz, dw2min)) * zdelm *zdelm 
                    
           zmetk  =  azmet(jk)* rh4                     ! mid-layer height k_int => k_int+1
           zgrow = exp(zmetk)
           ritur = bn2(jk)/shr2         
           kamp = sqrt(shr2)*sc2 *zgrow
           w1 = 1./(1. + 5*ritur)
           ktur= min(max(kamp * w1 * w1, dked_min)+kvg(k), dked_max) 
	   zmetk =  azmet(jk)* rhp       
           vueff(jk)  = ktur*0. + 2.e-5*exp( zmetk)
          enddo

        if (idebug_gwrms == 1) then          
         do jk= km1,levs
	   wrk1(jk) = rv_kxw/rhoint(jk) 
	   wrk2(jk)= rgrav2*zthm*zthm*bn2(jk)	! dimension [K*K]*(c2/m2)	
          enddo
        endif     

!
! extrapolating values for ktop = levs+1 (lev-interface for prsi(levs+1) =/= 0)
!     
         jk = levs
	 
         suprf(ktop) = kion(jk)
        
           rhoint(ktop) = aprsi(ktop)*rdi/atm(jk)
	   
           uint(ktop)  = aum(jk)
           vint(ktop)  = avm(jk)
           
           v_zmet(ktop) =  v_zmet(jk)
           vueff(ktop) = vueff(jk)
           bn2(ktop)   = bn2(jk) 
	    bn(ktop)   = bn(jk)
	    bnk(ktop) = bn(ktop)*v_kxw
	    
           rhobnk(ktop) = rhoint(ktop)/bnk(ktop)*v_kxw	 
	            
           bvi = bn(ksrc);  bvi2 = bvi * bvi;  
	   bvi3 = bvi2*bvi; bvi4 = bvi2 * bvi2; rcms = zms/bvi
	   bnrh_src =  bvi/rhoint(ksrc)
!	
!        define intrinsic velocity (relative to ilaunch) u(z)-u(zo), and coefficinets
!       ------------------------------------------------------------------------------------------        
        do iaz=1, nazd
           ul(iaz) = zcosang(iaz) *uint(ksrc) + zsinang(iaz) *vint(ksrc)
        enddo
!
         do jk=ksrc, ktop   	 
           do iaz=1, nazd
               zu = zcosang(iaz)*uint(jk) + zsinang(iaz)*vint(jk)
               ui(iaz, jk) =  zu   !- ul(iaz)*0.
             enddo
          enddo
!      ----------------------------------------- 
!       set launch momentum flux spectral density
!       ----------------------------------------- 

         fpu(1, ksrc) =0.
        do inc=1,nwav
           zcin  = zci(inc)
           zcin4 = zci4(inc)/bvi4
!           
          if(nslope == 0) then
	    zcin3 = zci3(inc)/bvi3
	    flux(inc,1) = zcin/(1.+zcin3)
          endif
	  
          if(nslope == 1) flux(inc,1) = zcin/(1.+zcin4)       
          if(nslope == 2) flux(inc,1)=  zcin/(1.+zcin4*zcin*rcms)  
        
! integrate (flux x dx)
            fpu(1,ksrc) = fpu(1,ksrc) + flux(inc,1)*zdci(inc)
	    
           do iaz=1,nazd   	    
             akzw(inc, iaz,ksrc:ktop) = bvi*rci(inc)
	   enddo
	    
         enddo
!                
         flux_norm  = taub_src / fpu(1, ksrc)   
!
      do iaz=1,nazd        
         fpu(iaz, ksrc) = taub_src
      enddo

! adjust rho/bn vertical factors for saturated fluxes (E(m) ~m^-3)
      bnrh_src=bnrh_src*flux_norm
      do jk=ksrc, ktop
          fden_bn(jk) = bnrh_src*rhoint(jk) / bn(jk)               !*bvi/rhoint(ksrc) 
      enddo
         
!                           
      do inc=1, nwav
          flux(inc,1) = flux_norm*flux(inc,1)
      enddo

       if (idebug_gwrms == 1) then      
	  pwrms =0.
	  ptrms =0.	  	  
	  tx1 = real(nazd)/rhoint(ksrc)*rv_kxw
	  ze2 = wrk2(ksrc)    !  (bvi*atm(ksrc)*rgrav)**2  
	 do inc=1, nwav
	  v_kzw = bvi*rci(inc)
	  ze1 = flux(inc,1)*zdci(inc)*tx1*v_kzw
	  pwrms = pwrms + ze1
	  ptrms = ptrms + ze1 * ze2	  
	 enddo
	   wrms(jl, ksrc) = pwrms 
	   trms(jl, ksrc) = ptrms 	 
        endif
  
!     copy flux-1 into other azimuths
!     --------------------------------
      do iaz=2, nazd
        do inc=1,nwav
            flux(inc,iaz) = flux(inc,1)
          enddo
      enddo	
      	
! constant flux below  ilaunch
        do jk=km1, ksrc
         do inc=1, nwav	
	     flux_z(inc,jk)=flux(inc,1)
	  enddo   
	enddo
	
	wave_act(:,:) = 1.0	
!                                        vertical do-loop
        do jk=ksrc, levs
	   jkp = jk+1
!                                        azimuth do-loop
         do iaz=1, nazd	

 	     umfp = ui(iaz, jkp)
             umfm = ui(iaz, jk)
             umfc = .5*(umfm + umfp)
!                                        wave-cin loop
          do inc=1, nwav
          
            zcin  = zci(inc)          ! zcin =/0  by definition
            zcinc = rci(inc) 
            
             if(wave_act(inc,iaz) == 1.0) then 
!=======================================================================
! discrete mode
! saturated limit    wfit = kzw*kzw*kt; wfdt = wfit/(kxw*cx)*betat
! & dissipative      kzi = 2.*kzw*(wfdm+wfdt)*dzpi(k)     
!=======================================================================  
        
               v_cdp =  zcin - umfp
	       	                    
	     if (v_cdp .le. ucrit_max) then 
!
! between layer [k-1,k or jk-jkp]  (Chi - Uk) -> ucrit_max ; wave's absorption
!
	        wave_act(inc,iaz) =0. 
                akzw(inc, iaz, jkp) = pi/v_zmet(jk)   ! pi2/dzmet
		fluxs = 0.0                           !max(0., rhobnk(jkp)*ucrit3)*rdci(inc)                            
		flux(inc,iaz)   = fluxs
                flux_z(inc,jkp) = fluxs
!		ucrit_maxdc =0.
	     else 
	           
              v_wdp = v_kxw*v_cdp
              wdop2 = v_wdp* v_wdp
              v_cdp2=v_cdp*v_cdp
!	      
! rotational cut-off
!	      
              cdf2  = v_cdp2 - c2f2
	      
              if (cdf2 > 0.0) then
                kzw2 = (bn2(jkp)-wdop2)/Cdf2
              else
                kzw2 = mkz2min
              endif
	      
              if ( kzw2 > mkz2min ) then
                v_kzw = sqrt(kzw2)
                akzw(inc, iaz, jkp) = v_kzw
!	       
!linsatdis:  kzw2, kzw3, kdsat, c2f2,  cdf2, cdf1
!	             
!kzw2 = (bn2(k)-wdop2)/Cdf2  - rhp4 - v_kx2w  ! full lin DS-NGW (N2-wd2)*k2=(m2+k2+[1/2H]^2)*(wd2-f2) 
!              Kds = kxw*Cdf1*rhp2/kzw3
!
                v_cdp  = sqrt( cdf2 )
                v_wdp  = v_kxw *  v_cdp
		v_wdi = kzw2*vueff(jk) + supRF(jk)      ! supRF - diss due to FRF-FV3dycore for "all" vars
		v_wdpc = sqrt(v_wdp*v_wdp +v_wdi*v_wdi)
		v_kzi  = v_kzw*v_wdi/v_wdpc
!       
                ze1 = v_kzi*v_zmet(jk)  
		 
                if (ze1 .ge. 1.e-2) then 
		    expdis = max(exp(-ze1), 0.01)
		else
		    expdis = 1./(1.+ ze1) 
		endif     
				
!
		wave_act(inc,iaz) = 1.0
                fmode =  flux(inc,iaz)
						
              else                                    ! kzw2 <= mkz2min large "Lz"-reflection
                
                expdis = 1.0
                v_kzw  = mkzmin
		
                v_cdp  = 0.                          ! no effects of reflected waves
                wave_act(inc,iaz) = 0.0
                akzw(inc, iaz, jkp) = v_kzw
		fmode = 0.
              endif
       
              fdis  =  fmode*expdis
!
! saturated flux + wave dissipation - Keddy_gwsat in UGWP-V1
!  linsatdis = 1.0 , here:   u'^2 ~ linsatdis* [v_cdp*v_cdp]
!
!                     fluxs= fden_bn(jkp)*cdf2*zcinc
                    fluxs= fden_bn(jkp)*sqrt(cdf2)
		   
!                                     
! S2003             fluxs= fden_bn(jk)*(zcin-ui(jk,iaz))**2/zcin
! WM2001            fluxs= fden_bn(jk)*(zcin-ui(jk,iaz))
!
              zdep = wave_act(inc,iaz)* (fdis-fluxs)
              if(zdep > 0.0 ) then
! subs on sat-limit
                 flux(inc,iaz)   = fluxs
                 flux_z(inc,jkp)  = fluxs
              else
! assign dis-ve flux
                flux(inc,iaz) =   fdis
                flux_z(inc,jkp) = fdis
              endif
	      
!             cgz = bnk(jk)/max(mkz2min, kzw2)	
     
	     dtau = flux_z(inc,jk)-flux_z(inc,jkp)
	     if (dtau .lt. 0) flux_z(inc,jkp) = flux_z(inc,jk)
	     
!	     if (dtau .ge. ucrit_maxdc) then	     
!	     flux_z(inc,jkp) = max(flux_z(inc,jk)-ucrit_maxdc, 0.)
!	       ze1 = zci(inc)-umfc-ucrit_maxdc
!	     write(6,287) dzdt(jk)/cgz, dtau/ucrit_maxdc, flux_z(inc,jkp)*1.e3, fluxs*1.e3, jk, zci(inc), ze1
!	    	     	           
!	     endif
! 287         format(' dtau >ucrit_max', 4(2x, F12.7), I4, 2x, 2(2x,F8.3))		      
!

	      endif  !  coriolis or CL condition-checkif  => (v_cdp .le. ucrit_max) then 
             endif   !  only for waves w/o CL-absorption   wave_act=1
	     
		      
!	    
          enddo      ! wave-inc-loop 
!
! integrate over spectral modes  fpu(y, z, azimuth)    wave_act(jl,inc,iaz)*flux(jl,inc,iaz)*[d("zcinc")]
!
        if (idebug_gwrms == 1) then            
	   pwrms =0.
	   ptrms =0.	  
! new arrays
	  	  	  
	 do inc=1, nwav	
	  if (wave_act(inc,iaz) > 0.) then 
           v_kzw =akzw(inc, iaz, jk) 
           ze1 = flux(inc,iaz)*v_kzw*zdci(inc)*wrk1(jk)
	   pwrms = pwrms + ze1
	   ptrms = ptrms + ze1*wrk2(jk)	  
	  endif 	  
	 enddo	  
           Awrms(iaz, jk) = pwrms 
	   Atrms(iaz, jk) = ptrms 	  
	endif


          dfdz_v(iaz, jk) = 0.0
          fpu(iaz, jkp)    = 0.0
	
          do inc=1, nwav
	    if (wave_act(inc,iaz) > 0.) then 
             
	      zcinc =zdci(inc)            
              vc_zflx_mode      = flux(inc,iaz)
              fpu(iaz, jkp) = fpu(iaz,jkp) + vc_zflx_mode*zcinc
              
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! (heat deposition integration over spectral mode for each azimuth
!      later sum over selected azimuths as "non-negative" scalars)
!      cdf1 = sqrt( (zci(inc)-umfc)**2-c2f2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             zdelp = wrk3(jk)*cdf1 *zcinc        
             zdelp = wrk3(jk)*abs(zci(inc)-umfc) *zcinc
             vm_zflx_mode = flux_z(inc,jk)	     
             dfdz_v(iaz, jk) = dfdz_v(iaz,jk) +(vm_zflx_mode-vc_zflx_mode)*zdelp    ! heating >0
            endif
          enddo                            !waves inc=1,nwav

            ze1 =fpu(iaz, jk)
            if (fpu(iaz, jkp) > ze1 ) fpu(iaz, jkp) = ze1
! --------------
        enddo                              ! end  Azimuth do-loop
	
! 
! extra- eddy wave dissipation  to  limit GW-rms        
!	    tx1 = sum(abs(dfdz_v(jk,1:nazd)))/bn2(jk)
!	    ze1=max(dked_min, tx1)
!	    ze2=min(dked_max, ze1)
!	    vueff(jkp) = ze2 + vueff(jkp)	       
!

	 
       enddo                               ! end Vertical do-loop
!
! top-layers constant interface-fluxes and zero-heat
!           
            fpu(1:nazd,ktop) = fpu(1:nazd, levs) 
            dfdz_v(1:nazd, levs) = 0.0
             
! ---------------------------------------------------------------------    
!       sum contribution for total zonal and meridional fluxes +
!           energy dissipation
!       ---------------------------------------------------
!       
!========================================================================
! at the source level and below taux = 0 (taux_E=-taux_W by assumption)
!========================================================================
   
    
       
         do jk=ksrc,  levs
            taux(jk) = 0.0 
            tauy(jk) = 0.0     
           do iaz=1,nazd
             taux(jk)  = taux(jk)  + fpu(iaz,jk)*zcosang(iaz)     
             tauy(jk)  = tauy(jk)  + fpu(iaz,jk)*zsinang(iaz)  
             pdtdt(jl,jk) = pdtdt(jl,jk)+ dfdz_v(iaz,jk)
          enddo	  
         enddo
             jk = ktop; taux(jk)=0.; tauy(jk)=0.
          do iaz=1,nazd           
             taux(jk)  = taux(jk)  + fpu(iaz,jk)*zcosang(iaz)     
             tauy(jk)  = tauy(jk)  + fpu(iaz,jk)*zsinang(iaz)           
          enddo
       
       if (idebug_gwrms == 1) then 

         do jk=kp1,  levs  
            do iaz=1,nazd           
            wrms(jl,jk)  =wrms(jl,jk)   + Awrms(iaz,jk)
            trms(jl,jk)  =trms(jl,jk)   + Atrms(iaz,jk)	
	    tauabs(jl,jk)=tauabs(jl,jk) + fpu(iaz,jk)
          enddo
	 enddo     
 
       endif	           
!

       do jk=ksrc,levs
          jkp = jk + 1
           zdelp = wrk3(jk)*gw_eff           
	   ze1   = (taux(jkp)-taux(jk))* zdelp
           ze2   = (tauy(jkp)-tauy(jk))* zdelp
	   
           if (abs(ze1) >= maxdudt ) then
             ze1 = sign(maxdudt, ze1)
           endif 
           if (abs(ze2) >= maxdudt ) then
             ze2 = sign(maxdudt, ze2)
           endif
           pdudt(jl,jk) = -ze1
           pdvdt(jl,jk) = -ze2
!
! Cx =0 based Cx=/= 0. above
!
! 
        if (knob_ugwp_doheat == 1) then
!	
!	  ek1 =aum(jk)*aum(jk) +avm(jk)*avm(jk)
!	  ugw = aum(jk)- ze1*dtp; vgw = avm(jk)- ze2*dtp
!	  ek2 = ugw*ugw +vgw*vgw
!         pdtdt(jl,jk) = rdtp2*max(ek1-ek2, 0.0)                !=ze1*um + 0.5*ze1^2*dtp
!         pdtdt(jl,jk) = max(ze1*aum(jk) + ze2*avm(jk), 0.)     ! gw_eff => in "ze1 and ze2"
             pdtdt(jl,jk) = max(pdtdt(jl,jk) , 0.)*gw_eff
        endif

           if (abs(pdtdt(jl,jk)) >= maxdtdt ) pdtdt(jl,jk) = maxdtdt
           ze1  = max(dked_min, pdtdt(jl,jk)/bn2(jk))
           dked(jl,jk)  = min(dked_max, ze1)	   
     
       enddo
!	
! add limiters/efficiency for "unbalanced ics" if it is needed
!     
       do jk=ksrc,levs         
           pdtdt(jl,jk) = pdtdt(jl,jk)*rcpd
         enddo          
! 	    	 
	  dktur(1:levs) = dked(jl,1:levs)
!
       do ist= 1, 3	  
         do jk=ksrc,levs-1
           adif(jk)  =  .25*(dktur(jk-1)+ dktur(jk+1)) + .5*dktur(jk)     
          enddo
          dktur(ksrc:levs-1) = adif(ksrc:levs-1)
       enddo
         
!	 dked(jl, ksrc:levs-1) = dktur(ksrc:levs-1)
!	 dked(jl, levs) =dked(jl, levs-1)
	 	        
!
! perform "diffusive" 3-point smoothing of "u-v-t"
! from the surface to the "top"
!      
      if (knob_ugwp_dokdis == 2) then
             
	  uold(1:levs) = aum(1:levs)+pdudt(jl,1:levs)*dtp
          vold(1:levs) = avm(1:levs)+pdvdt(jl,1:levs)*dtp
          told(1:levs) = atm(1:levs)+pdtdt(jl,1:levs)*dtp
	  
	  do jk=1,levs 
	  zmetk= azmet(jk)*rhp
	  ktur  = kvg(k) + 2.e-5*exp( zmetk)
	  dktur(jk) = dked(jl,jk) + ktur	  	  	  
	  enddo
	  
	  dzmetm= azmet(ksrc)-  azmet(ksrc-1)
	  	       
          do jk=2,levs-1  
	    dzmetf = (azmeti(jk+1)-  azmeti(jk))*rhomid(jk)	  
	    ktur = .5*(dktur(jk-1)+dktur(jk)) *rhoint(jk)/dzmetf 
	    kturp = .5*(dktur(jk+1)+dktur(jk))*rhoint(jk+1)/dzmetf 
	    	    
	    dzmetp =  azmet(jk+1)-azmet(jk)
	    Adif(jk) = ktur/dzmetm
	    Cdif(jk) = kturp/dzmetp
	    bdif = adif(jk)+cdif(jk)
	    if (rdtp  <  bdif ) then
	      Anstab(jk) = nint( bdif/rdtp + 1) 
	    else
	     Anstab(jk) = 1
	    endif
	    dzmetm = dzmetp
	   enddo  
	     
	    nstab = maxval( Anstab(ksrc:levs-1))
	    if (nstab .ge. 2) print *, 'nstab ', nstab
	    dtdif = dtp/real(nstab)
	    do ist= 1, nstab
             do k=ksrc,levs-1  	     
	      Bdif = nstab*rdtp-Adif(k)-Cdif(k)
	      unew(k) = uold(k)*Bdif+ uold(k-1)*Adif(k) + uold(k)*Cdif(k)
	      vnew(k) = vold(k)*Bdif+ vold(k-1)*Adif(k) + vold(k)*Cdif(k)
	      tnew(k) = told(k)*Bdif+ told(k-1)*Adif(k) + told(k)*Cdif(k)  
	     enddo
	      uold = unew*dtdif
	      vold = vnew*dtdif
	      told = tnew*dtdif
            enddo
!
! create "smoothed" tendencies by molecular + GW-eddy diffusion
!	    
             do k=ksrc,levs-1 	    
	      pdtdt(jl,jk)= rdtp*(told(k) - tm(jl,k))
	      ze2 = rdtp*(uold(k) - aum(k))
	      ze1 = rdtp*(vold(k) - avm(k))
	      if (abs(pdtdt(jl,jk)) >= maxdtdt ) pdtdt(jl,jk) = maxdtdt
           if (abs(ze1) >= maxdudt ) then
             ze1 = sign(maxdudt, ze1)
           endif 
           if (abs(ze2) >= maxdudt ) then
             ze2 = sign(maxdudt, ze2)
           endif
	       pdudt(jl, k) = ze2
	       pdvdt(jl, k) = ze1
!
! add eddy viscosity heating
!            pdtdt(jl,jk) = pdtdt(jl,jk) - max(ze1*aum(jk) + ze2*avm(jk), 0.) *rcpd 
!	       	      
	     enddo
	     
	     
	  ENDIF    !  dissipative IF-loop for "abrupt" tendencies
	     
       enddo        ! J-loop
!    


       RETURN    
       
!  
! Print/Debugging  ----------------------------------------------------------------------- 
!
 239   continue
       if (kdt ==1 .and. mpi_id == master) then
!      
          print *, 'ugwp-vay: nazd-nw-ilaunch=', nazd, nwav,ilaunch, maxval(kvg), ' kvg '
          print *,  'ugwp-vay: zdci(inc)=' ,  maxval(zdci), minval(zdci)
          print *,  'ugwp-vay: zcimax=' , maxval(zci) ,' zcimin=' , minval(zci)  
!         print *,  'ugwp-vay: tau_ngw=' , maxval(taub_src)*1.e3,  minval(taub_src)*1.e3, tau_min

          print *
	
        endif
 
       if (kdt == 1 .and. mpi_id == master) then
         print *, 'vgw done nstab ', nstab
!
         print *, maxval(pdudt)*86400.,  minval(pdudt)*86400, 'vgw ax  ugwp'
         print *, maxval(pdvdt)*86400.,  minval(pdvdt)*86400, 'vgw ay  ugwp'
         print *, maxval(dked)*1.,  minval(dked)*1,  'vgw keddy m2/sec ugwp'
         print *, maxval(pdtdt)*86400.,  minval(pdtdt)*86400,'vgw eps  ugwp'
!
!        print *, ' ugwp -heating rates '
       endif



       return
       end subroutine cires_ugwp_solv2_v1


end module cires_ugwp_solv2_v1_mod

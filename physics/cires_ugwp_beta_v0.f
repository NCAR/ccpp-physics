!===============================================================================
!    use fv3gfs-v0
!    first beta version of ugwp for fv3gfs-128l
!    cires/swpc - jan 2018
!    non-tested wam ugwp-solvers in fv3gfs: "lsatdis", "dspdis", "ado99dis"
!    they reqiure extra-work to put them in with intializtion and namelists
!          next will be lsatdis for both fv3wam & fv3gfs-128l implementations
!    with (a) stochastic-deterministic propagation solvers for wave packet/spectra
!         (b) gw-sources: oro/convection/dyn-instability (fronts/jets/pv-anomalies)
!         (c) guidance from high-res runs for gw sources and res-aware tune-ups  
!23456 
!================================2017:    new wam scheme ======================
!    scinocca-solver adapted for "surface" -> top like in wam
!=========================================================================  
!
!      call gwdrag_wam(1,  im,   ix,  levs,   ksrc, dtp,
!     & xlat, precip, gw_dudt, gw_dvdt,  taux,  tauy)
!        call fv3_ugwp_wms17(kid1, im,  ix,  levs,  ksrc_ifs, dtp,
!     &  adt,adu,adv,prsl,prsi,phil,xlat,precip,gw_dudt, gw_dvdt, gw_dtdt, gw_ked, taux,tauy,
!     &  grav, amol_i, me, lstep_first )
!
!
!23456==============================================================================
      subroutine fv3_ugwp_betav0(klon, klev, dtime,
     & tm1 , um1, vm1,  prsl,  prsi, phil981, xlatd, sinlat, coslat,
     & ptenu, ptenv, ptent, dked, mpi_id, master, kdt)
!
!=======================================================
!
!      nov 2015 alternative gw-solver for nggps-wam
!      nov 2017 nh/rotational gw-modes for nh-fv3gfs
! ---------------------------------------------------------------------------------
!  
!       use physcons,      only :  rg1  => con_rog,  rd=>con_rd
!       use physcons,      only :  grav => con_g, cpd => con_cp
!       use physcons,      only :  rcpd1 => con_rocp, pi=>con_pi   
      use ugwp_common ,      only : rgrav, grav, cpd, rd, rv, rcpd2      
      use ugwp_common ,      only : pi, rad_to_deg, deg_to_rad, pi2       
      use ugwp_common ,      only : rdi, gor,  grcp,  gocp,  fv,  gr2
      use ugwp_common ,      only : bnv2min, dw2min, velmin                
!
! 
       implicit none
!23456 
       
        integer  ,intent(in)  :: klev              ! vertical level
        integer  ,intent(in)  :: klon              ! horiz tiles
	
        real ,intent(in)   :: dtime              !model time step       
        real ,intent(in)   :: vm1(klon,klev)      !full model level meridional velocity 
        real ,intent(in)   :: um1(klon,klev)      !full model level meridional velocity

        real ,intent(in)   :: tm1(klon,klev)      !full model level pot-temperature 

        real ,intent(in) :: prsl(klon,klev)       ! prsl full model level pressure
        real ,intent(in) :: phil981(klon,klev)    ! phil full model level geopotential dt) => meters !!!!!       phil =phil981/g981
        real ,intent(in) :: prsi(klon,klev+1)     !   prsi interface pressure
        real ,intent(in) :: xlatd(klon)           ! lat was in radians, now with xlat_d in degrees
        real ,intent(in) :: sinlat(klon)
        real ,intent(in) :: coslat(klon)          		
	 	
        integer, intent(in)   :: mpi_id, master, kdt
!       logical, intent(in)   :: lstep_first   
!
! out-gw effects
!
        real ,intent(out):: ptenu(klon,klev)     !full-model level zonal momentum tendency
        real ,intent(out):: ptenv(klon,klev)     !full-model level meridional momentum tendency
        real ,intent(out):: ptent(klon,klev)     ! gw-heating (u*ax+v*ay)/cp
        real ,intent(out):: dked(klon,klev)      ! gw-eddy diffusion
	
	
!         
      real, parameter      :: maxdudt = 250.e-5, gw_eff = 1.0	    
      integer, parameter   :: nslope = 1
      integer, parameter   :: klaunch=55       ! 32 - ~ 1km ;55 - 5.5 km ; 52 4.7km ; 60-7km index for selecting launch level     
      
        real,      parameter  :: omega2 = 2.*6.28/86400
        real,      parameter  :: gptwo=2.0       
        integer  , parameter  :: iazidim=4       ! number of azimuths
        integer  , parameter  :: incdim=25       !number of discrete cx - spectral elements in launch spectrum
        
        real ,    parameter   :: dked_min =0.01
        real ,    parameter   :: gssec = (6.28/30.)**2        ! max-value for bn2
        real ,    parameter   :: bv2min = (6.28/60./120.)**2  ! min-value for bn2  7.6(-7)  2 hrs
        real, parameter       :: minvel = 0.5       
        real, parameter       :: v_kxw = 6.28e-3/200.
        real, parameter       :: v_kxw2 = v_kxw*v_kxw
        real, parameter       :: tamp_mpa = 30.e-3
        real, parameter       :: zfluxglob= 3.75e-3
        real, parameter       :: zms_l    = 2000.0                     !m* (given as length, in m)       
              
!vay-2018
        real, dimension(klon) :: tau_merra  
	
        real             :: pprecip(klon)        !total surface precipitation
        real             :: agrav(klon, klev)    ! wam gravity(z,x) m/s2 < 9.81
        real             :: kmol(klon,klev+1)    ! molec-dissipation (visc) m2/s
        real             :: pfluxv(klon,klev+1)  != meridional component of vertical momentum flux (pa)
        real             :: pfluxu(klon,klev+1)  != meridional component of vertical momentum flux (pa)
        real             :: phil(klon,klev)      !gphil/9.81	
!
!work==========================================================================================================
!
	 
         real  :: zbvfhm1(klon,klev), zbn2(klon,klev)	 
	 real  :: kzw2, kzw3, kdsat, cdf2, cdf1, wdop2
	 real  :: c2f2(klon)	 
	 
         real  :: zuhm1(klon,klev)             !half-model level zonal velocity
         real  :: zvhm1(klon,klev)             !half-model level meridional velocity
         real  :: zrhohm1(klon,klev)           !half-model level densityaisalla frequency
         real  :: zx(incdim)                   !coordinate transformation
         real  :: zci(incdim)                  !phase speed element
         real  :: zdci(incdim)
!23456
       real  :: zul(klon,iazidim)            !velocity in azimuthal direction at launch level
       real  :: zbvfl(klon)                  !buoyancy at launch level
       real  :: zcosang(iazidim)             !cos of azimuth angle
       real  :: zsinang(iazidim)             !sin of azimuth angle
       real  :: zfct(klon,klev)
       real  :: zfnorm(klon)                 !normalisation factor (a)
       real  :: zci_min(klon,iazidim)
       real  :: zthm1(klon,klev)             !temperature on half-model levels
!

       real  :: zdfdz_v(klon,klev,iazidim)   ! axj = -df*rho/dz       directional momentum depositiom
       real  :: zflux(klon,incdim,iazidim)   ! momentum flux at each level   stored as ( ix, mode, iazdim)

       real  :: zflux_z (klon,incdim,klev)    !momentum flux at each azimuth stored as ( ix, mode, klev)
!
       real  :: wave_flux_z(klon,incdim,klev,iazidim)    !momentum flux at each azimuth stored as ( ix, mode, klev, iazi)
!
       real  :: vm_zflx_mode, vc_zflx_mode
       real  :: zact(klon,incdim,iazidim)    !if =1 then critical level encountered => c-u
       real  :: zacc(klon,incdim,iazidim)
!
       real  :: zpu(klon,klev,iazidim)       !momentum flux
       real  :: zdfl(klon,klev,iazidim)

       real  :: zcrt(klon,klev,iazidim)
       integer   :: ilaunch                   !model level from which gw spectrum is launched

       
       real  :: zcimin, zcimax
       real  :: zgam, zpexp, zxmax, zxmin, zxran, zdx, zx1, zx2
       real  :: zang, zaz_fct, znorm, zang1, ztx
       real  :: zu, zcin, zcpeak, zcin4, zbvfl4
       real  :: zcin2, zbvfl2, zcin3, zbvfl3, zcinc
       real  :: zatmp, zfluxs, zdep, zfluxsq, zulm, zdft, ze1, ze2

!       real  :: zgauss(klon), zfluxlaun(klon), zcngl(klon)fluxlaun(klon), zcngl(klon)

       real  :: zcons1,zcons2,zdelp,zrgpts
       real  :: zthstd,zrhostd,zbvfstd
       
       real  :: zhook_handle
       real  ::  zms, zfluxlaun(klon)
       real  ::  zui(klon, klev,iazidim), zcngl(klon)

       real  ::  rcpd, grav2cpd
       
       real :: v_zmet(klon, klev), fmode, expdis, fdis
       real :: vueff(klon, klev)
       real :: v_kzi, v_kzw, v_cdp, v_wdp, sc
       

           
       integer :: j, k, inc, jk, jl, iazi
              
        
        zms=2.*pi/zms_l    
        zpexp=gptwo/2.0                    ! gptwo=2 , zpexp = 1.   
!       
!--------------------------------------------------------------------------
!
      
        phil =phil981/grav
       
                       	
         rcpd = 1./(grav/cpd)        ! 1/[g/cp]
         grav2cpd=grav*grav/cpd      ! g*(g/cp)= g^2/cp
	 
         zcimin= 0.50               !     zcimax=100.0
         zgam=0.25                  !     m/s precision for crit-level  
         zcimax=  125.0             !75-95 -based
         ilaunch= klaunch
!   
!
! vay-2017 lat-dep forcing !!!!
!                                                       
!         pprecip(1:klon) = 0.0
!         agrav(:,:)     = grav
!         kmol(:, :)     = 1.e-5
!
         if (kdt ==1 .and. mpi_id == master) then

         print *,  maxval(tm1), minval(tm1), 'vgw: temp-res '
         print *,  'ugwp-v0: zcimin=' , zcimin
         print *,  'ugwp-v0: zcimax=' , zcimax
         print *,  'ugwp-v0: cd_crit=',  zgam           ! m/s precision for crit-level 
         print *,  'ugwp-v0: klaunch/ilaunch',  klaunch, ilaunch
         print *, ' ugwp-v0: vgw zms_l=', zms_l
         print *, ' ugwp-v0: vgw zpexp=', zpexp
         print *, ' ugwp-v0: zfluxglob=' , zfluxglob
         print *, ' ugwp-v0: nslope=', nslope
         print *, 'ugwp-v0 int-ers', klon,  klev 
         print *
	 
         endif	 
!
!=================================================
       do iazi=1,iazidim
          do jk=1,klev
          do jl=1,klon                       
         zpu(jl,jk,iazi) =0.0 
         zcrt(jl,jk,iazi)=0.0 
         zdfl(jl,jk,iazi)=0.0 
         enddo
         enddo      
       enddo  

!vay-ok   
! set initial min ci in each column and azimuth (used for critical levels)
       do iazi=1,iazidim
        do jl=1,klon
            zci_min(jl,iazi)=zcimin
        enddo
       enddo
!       define half model level winds and temperature
!       -----------------------------------
       do jk=2,klev
        do jl=1,klon                                        
        zthm1(jl,jk) =0.5 *(tm1(jl,jk-1)+tm1(jl,jk)) 
        zuhm1(jl,jk) =0.5 *(um1(jl,jk-1)+um1(jl,jk))
        zvhm1(jl,jk) =0.5 *(vm1(jl,jk-1)+vm1(jl,jk))
        enddo
       enddo

         jk=1
        do jl=1,klon                                        
        zthm1(jl,jk)=tm1(jl,jk) 
        zuhm1(jl,jk)=um1(jl,jk)
        zvhm1(jl,jk)=vm1(jl,jk)
        enddo
!ok
!
       do jk=2, klev
         do jl=1,klon
       zdelp=phil(jl,jk)-phil(jl,jk-1)      !>0 ...... dz-meters
!
       v_zmet(jl,jk) = 2.*zdelp
       vueff(jl,jk) = 
     &  2.e-5*exp(.5*(phil(jl,jk)+phil(jl,jk-1))/7000.)+dked_min


       zrhohm1(jl,jk)=prsi(jl,jk)*rdi/zthm1(jl,jk)                !  rho = p/(RT) 
       zbn2(jl,jk)= grav2cpd/zthm1(jl,jk)*
     &    (1.0 + rcpd*(tm1(jl,jk)-tm1(jl,jk-1))/zdelp)
        
!
!    (r*r*g*g)((1/cp)/t = [g/t]*(g/cp)*(r)*[1.+dt/z/(g/cp)]
!                         (g/t)*[ g/cp + dt/dz]   
!
!       bn2 = g/t[ g/cp - g/cp*cp*dt/dp] = g/t*(ga +dt/dz=> -dt/dp*rho*g)
! make more sense .......
! 
         zbn2(jl,jk)=min(zbn2(jl,jk),gssec)
         zbn2(jl,jk)=max(zbn2(jl,jk),bv2min)
!
!  ifs-error       zbvfhm1(jl,jk)= max (zbvfhm1(jl,jk),gssec)
!
         zbvfhm1(jl,jk)=sqrt(zbn2(jl,jk))       ! bn = sqrt(bn2)
!
       enddo
       enddo
    

       DO JL=1,klon
         ZBVFHM1(JL,1)  = ZBVFHM1(JL,2)
          V_ZMET(JL,1)  = V_ZMET(JL,2)
          VUEFF(JL,1)   = KMOL(JL,1) + DKED_MIN
	   ZBN2(JL,1)     = ZBN2(JL,2)
	 C2F2(JL) = (OMEGA2*SINLAT(JL))**2/V_KXW2
       ENDDO
!
!      set up azimuth directions and some trig factors
!
!
           zang=2.*pi/iazidim
           zaz_fct=1.0 
! get normalization factor to ensure that the same amount of momentum
! flux is directed (n,s,e,w) no mater how many azimuths are selected.
! note, however, the code below assumes a symmetric distribution of
! of azimuthal directions (ie 4,8,16,32,...)
          znorm=0.0 
        do iazi=1,iazidim
         zang1=(iazi-1)*zang
         zcosang(iazi)=cos(zang1)
         zsinang(iazi)=sin(zang1)
         znorm=znorm+abs(zcosang(iazi))
        enddo
         zaz_fct=2. *zaz_fct/znorm
!       define coordinate transform
!       -----------------------------------------------     
! note that this is expresed in terms of the intrinsic phase speed
! at launch ci=c-u_o so that the transformation is identical at every
! launch site.
! see eq. 28-30 of scinocca 2003.   x = 1/c stretching transform
!
          zxmax=1.0 /zcimin
          zxmin=1.0 /zcimax
          zxran=zxmax-zxmin
          zdx=zxran/real(incdim-1)     ! dkz
! vay-error
!
          zx1=zxran/(exp(zxran/zgam)-1.0 )           !zgam =1./4.
          zx2=zxmin-zx1
!
!                                   add idl computations for zci =1/zx
!                                                   x = 1/c stretching transform
!                                
!                                   zx1=zxran/(exp(zxran/zgam)-1.0_jprb)
!                                   zx2=zxmin-zx1

        do inc=1,incdim
          ztx=real(inc-1)*zdx+zxmin
          zx(inc)=zx1*exp((ztx-zxmin)/zgam)+zx2                       !eq. 29 of scinocca 2003
          zci(inc)=1.0 /zx(inc)                                       !eq. 28 of scinocca 2003
          zdci(inc)=zci(inc)**2*(zx1/zgam)*exp((ztx-zxmin)/zgam)*zdx  !eq. 30 of scinocca 2003
        enddo
	
!        define intrinsic velocity (relative to launch level velocity) u(z)-u(zo), and coefficinets
!       ------------------------------------------------------------------------------------------        
        do iazi=1,iazidim
         do jl=1,klon
         zul(jl,iazi)=zcosang(iazi)*zuhm1(jl,ilaunch)+
     &                zsinang(iazi)*zvhm1(jl,ilaunch)
         enddo
        enddo
	
	
        do jl=1,klon
         zbvfl(jl)=zbvfhm1(jl,ilaunch)
        enddo
!
!======================================= ok 
!
         do jk=ilaunch, klev-1     ! from z-launch up   model level from which gw spectrum is launched
          do iazi=1,iazidim
           do jl=1,klon
          zu=zcosang(iazi)*zuhm1(jl,jk)+zsinang(iazi)*zvhm1(jl,jk)
          zui(jl,jk,iazi)=zu-zul(jl,iazi)
          enddo
        enddo

        enddo
!                                         define rho(zo)/n(zo)
!       ------------------- 
      do jk=ilaunch, klev-1     ! from z-launch up     2,ilaunch
        do jl=1,klon
        zfct(jl,jk)=zrhohm1(jl,jk)/zbvfhm1(jl,jk)
        enddo
      enddo

!      ----------------------------------------- 
!       set launch momentum flux spectral density
!       ----------------------------------------- 
! do this for only one azimuth and copy for all other azimuths
!
!
      if(nslope==1) then
! s=1 case
       do inc=1,incdim
             zcin=zci(inc)
             zcin4=(zms*zcin)**4
       do jl=1,klon
!n4
         zbvfl4=zbvfl(jl)**4
         zflux(jl,inc,1)=zfct(jl,ilaunch)*zbvfl4*zcin/(zbvfl4+zcin4)
         zact(jl,inc,1)=1.0 
      enddo
      enddo
           elseif(nslope==2) then
! s=2 case
      do inc=1,incdim
      zcin=zci(inc)
      zcin4=(zms*zcin)**4
        do jl=1,klon
         zbvfl4=zbvfl(jl)**4
         zcpeak=zbvfl(jl)/zms
        zflux(jl,inc,1)=zfct(jl,ilaunch)*
     &     zbvfl4*zcin*zcpeak/(zbvfl4*zcpeak+zcin4*zcin)

         zact(jl,inc,1)=1.0 
         enddo
       enddo
          elseif(nslope==-1) then
! s=-1 case
       do inc=1,incdim
         zcin=zci(inc)
         zcin2=(zms*zcin)**2
       do jl=1,klon
        zbvfl2=zbvfl(jl)**2
        zflux(jl,inc,1)=zfct(jl,ilaunch)*zbvfl2*zcin/(zbvfl2+zcin2)
        zact(jl,inc,1)=1.0 
       enddo
       enddo
!s=0 case
           elseif(nslope==0) then

       do inc=1,incdim
           zcin=zci(inc)
           zcin3=(zms*zcin)**3
       do jl=1,klon
        zbvfl3=zbvfl(jl)**3
        zflux(jl,inc,1)=zfct(jl,ilaunch)*zbvfl3*zcin/(zbvfl3+zcin3)
        zact(jl,inc,1)=1.0 
        zacc(jl,inc,1)=1.0 
       enddo
       enddo

       endif  ! for slopes
!
! normalize momentum flux at the src-level
!       ------------------------------
! (rho x f^h = rho_o x f_p^total)
! integrate (zflux x dx)
      do inc=1,incdim
         zcinc=zdci(inc)
      do jl=1,klon
      zpu(jl,ilaunch,1)=zpu(jl,ilaunch,1)+zflux(jl,inc,1)*zcinc
      enddo    
      enddo
!       normalize and include lat-dep  (precip or merra-2)
!       -----------------------------------------------------------
! also other options to alter tropical values
!
! vay-2017   zfluxglob=> lat-dep here from geos-5/merra-2  tamp = 100.e-3*1.e3 = 100 mpa
!-----------------------------------------------------------
       call slat_geos5_tamp(klon, tamp_mpa, xlatd, tau_merra)
       do jl=1,klon       
       zfluxlaun(jl)=tau_merra(jl)     !*(.5+.75*coslat(JL))      !zfluxglob/2  oon poles
        zfnorm(jl)= zfluxlaun(jl)/zpu(jl,ilaunch,1)
       enddo
!
      do iazi=1,iazidim
      do jl=1,klon
      zpu(jl,ilaunch,iazi)=zfluxlaun(jl)
      enddo
      enddo

!       adjust constant zfct

         do jk=ilaunch, klev-1
        do jl=1,klon
            zfct(jl,jk)=zfnorm(jl)*zfct(jl,jk)
        enddo
        enddo
!               renormalize each spectral element in first azimuth

        do inc=1,incdim
         do jl=1,klon
          zflux(jl,inc,1)=zfnorm(jl)*zflux(jl,inc,1)
         enddo
         enddo

!       copy zflux into all other azimuths
!       --------------------------------
          do iazi=2,iazidim
            do inc=1,incdim
            do jl=1,klon
          zflux(jl,inc,iazi)=zflux(jl,inc,1)
          zact(jl,inc,iazi)=1.0 
          zacc(jl,inc,iazi)=1.0 
           enddo
          enddo
          enddo

! -------------------------------------------------------------
!       
!                                       iazidim do-loop
! --------------------
        do iazi=1,iazidim
!                                       vertical do-loop
! ----------------
          do jk=ilaunch, klev-1                     
! first check for critical levels
! ------------------------
           do jl=1,klon
            zci_min(jl,iazi)=max(zci_min(jl,iazi),zui(jl,jk,iazi))               
           enddo
! set zact to zero if critical level encountered
! ----------------------------------------------
           do inc=1,incdim
            zcin=zci(inc)
            do jl=1,klon
               zatmp=minvel+sign(minvel,zcin-zci_min(jl,iazi))
               zacc(jl,inc,iazi)=zact(jl,inc,iazi)-zatmp
               zact(jl,inc,iazi)=zatmp
            enddo
           enddo
!
! integrate to get critical-level contribution to mom deposition
! ---------------------------------------------------------------
           do inc=1,incdim
            zcinc=zdci(inc)
            do jl=1,klon
               zdfl(jl,jk,iazi)=zdfl(jl,jk,iazi)+
     &                zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zcinc
            enddo
           enddo
! --------------------------------------------
! get weighted average of phase speed in layer
! --------------------------------------------
          do jl=1,klon
            if(zdfl(jl,jk,iazi)>0.0 ) then
               zatmp=zcrt(jl,jk,iazi)
             do inc=1,incdim
                  zatmp=zatmp+zci(inc)*
     &                   zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zdci(inc)
             enddo
!
               zcrt(jl,jk,iazi)=zatmp/zdfl(jl,jk,iazi)
            else
               zcrt(jl,jk,iazi)=zcrt(jl,jk-1,iazi)
            endif
         enddo

!======================================================= default case
!        if(gptwo==2.0 ) then
!
! mode-loop +  add molecular/eddy dissipation of gw-spectra vay-2015
!
             do inc=1,incdim
                zcin=zci(inc)
                zcinc=1.0 /zcin
                do jl=1,klon
!=======================================================================
! saturated limit    wfit = kzw*kzw*kt; wfdt = wfit/(kxw*cx)*betat
! & dissipative      kzi = 2.*kzw*(wfdm+wfdt)*dzpi(k)
!           define   kxw = 
!=======================================================================  
               v_cdp =  abs(zcin-zui(jL,jk,iazi))
               v_wdp = v_kxw*v_cdp
	       wdop2 = v_wdp* v_wdp
	       cdf2 = v_cdp*v_cdp - c2f2(jL) 
	       if (cdf2 .gt. 0) then
	           kzw2 = (zBn2(jL,jk)-wdop2)/Cdf2  - v_kxw2  
		else
		 kzw2 = 0.0
	       endif	        	       
               if ( kzw2 .gt. 0 ) then 
	       v_kzw = sqrt(kzw2)
!	       
!linsatdis:  kzw2, kzw3, kdsat, c2f2,  cdf2, cdf1
!	             
!	       kzw2 = (zBn2(k)-wdop2)/Cdf2  - rhp4 - v_kx2w  ! full lin DS-NiGW (N2-wd2)*k2=(m2+k2+[1/2H]^2)*(wd2-f2) 
!              Kds = kxw*Cdf1*rhp2/kzw3
!
                v_cdp = sqrt( cdf2 ) 
	        v_wdp = v_kxw *  v_cdp    
                v_kzi  = abs(v_kzw*v_kzw*vueff(jl,jk)/v_wdp*v_kzw)
                  expdis = exp(-v_kzi*v_zmet(jl,jk))
               else
                  v_kzi = 0.
                  expdis = 1.0  
		  v_kzw = 0.
		  v_cdp = 0.   ! no effects of reflected waves
               endif
	       
               fmode =  zflux(jl,inc,iazi) 
               fdis  =  fmode*expdis
!
! saturated flux + wave dissipation - Keddy_gwsat in UGWP-V1
!  linsatdis = 1.0 , here:   u'^2 ~ linsatdis* [v_cdp*v_cdp]
!
               zfluxs= zfct(jl,jk)*v_cdp*v_cdp*zcinc
!                                     
!               zfluxs= zfct(jl,jk)*(zcin-zui(jl,jk,iazi))**2/zcin
! flux_tot - sat.flux
! 

               zdep=zact(jl,inc,iazi)* (fdis-zfluxs)		   
                  if(zdep>0.0 ) then
! subs on sat-limit
                      zflux(jl,inc,iazi)=zfluxs
                      zflux_z(jl,inc,jk)=zfluxs
                   else
! assign dis-ve flux
                      zflux(jl,inc,iazi)=fdis 
                      zflux_z(jl,inc,jk)=fdis 
                   endif  
               enddo                 
             enddo
!
! integrate over spectral modes  zpu(y, z, azimuth)    zact(jl,inc,iazi)*zflux(jl,inc,iazi)*[d("zcinc")]
!
           zdfdz_v(:,jk, iazi) =0.
        
          do inc=1,incdim
                 zcinc=zdci(inc)                    ! dc-integration
            do jl=1,klon
!orig       zpu(jl,jk,iazi)=zpu(jl,jk,iazi)+zact(jl,inc,iazi)*zflux(jl,inc,iazi)*zcinc

               vc_zflx_mode = zact(jl,inc,iazi)*zflux(jl,inc,iazi)
               zpu(jl,jk,iazi)=zpu(jl,jk,iazi) + vc_zflx_mode*zcinc
               wave_flux_z(jl,inc, jk,iazi) = vc_zflx_mode*zcinc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! check monotonic decrease
!     (heat deposition integration over spectral mode for each azimuth
!      later sum over selected azimuths as "non-negative" scalars)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (jk.gt.ilaunch)then
       zdelp= grav/(prsi(jl,jk-1)-prsi(jl,jk))*
     &        abs(zcin-zui(jl,jk,iazi)) *zcinc
       vm_zflx_mode=zact(jl,inc,iazi)* zflux_z(jl,inc,jk-1)

      if (vc_zflx_mode.gt.vm_zflx_mode) vc_zflx_mode =vm_zflx_mode ! no-flux increase
      zdfdz_v( jl,jk,iazi) =zdfdz_v( jl,jk,iazi) +
     &           (vm_zflx_mode-vc_zflx_mode)*zdelp                 ! heating >0            
!
!              
       endif              
            enddo                          !jl=1,klon             
         enddo                             !waves inc=1,incdim

! --------------
        enddo                              ! end jk do-loop vertical loop
! ---------------
       enddo                               ! end iazidim do-loop
! ----------------------------------------------------------------------------    
!       sum contribution for total zonal and meridional flux +
!           energy dissipation
!       ---------------------------------------------------
! !vay-ok       
       do jk=1,klev+1
        do jl=1,klon
         pfluxu(jl,jk)=0.0 
         pfluxv(jl,jk)=0.0 
       enddo
       enddo     
        ptenu(:,:)  =0.0
        ptenv(:,:)  =0.0
	ptent(:,:)  =0.0
        dked(:,:)   =0.0      
      do iazi=1,iazidim
        do jk=ilaunch,  klev-1      
        do jl=1,klon
      pfluxu(jl,jk)=pfluxu(jl,jk)+zpu(jl,jk,iazi)*zaz_fct*zcosang(iazi)       ! zaz_fct - "azimuth"-norm-n factor
      pfluxv(jl,jk)=pfluxv(jl,jk)+zpu(jl,jk,iazi)*zaz_fct*zsinang(iazi)
      ptent(jl,jk) = ptent(jl,jk)+ zdfdz_v(jl,jk,iazi)*zaz_fct/cpd         ! eps_dis =sum( +d(flux_e)/dz) > 0.
        enddo
        enddo

      enddo
!
!    update du/dt and dv/dt tendencies   ..... no contribution to heating => keddy/tracer-mom-heat
!    ----------------------------   
! 
 

        do jk=ilaunch,klev
        do jl=1, klon
        zdelp= grav/(prsi(jl,jk-1)-prsi(jl,jk))            
        ze1=(pfluxu(jl,jk)-pfluxu(jl,jk-1))*zdelp
        ze2=(pfluxv(jl,jk)-pfluxv(jl,jk-1))*zdelp  
	if (abs(ze1) .ge. maxdudt ) then
	    ze1 = sign(maxdudt, ze1)
	endif 
	if (abs(ze2) .ge. maxdudt ) then
	    ze2 = sign(maxdudt, ze2)
	endif 	
        ptenu(jl,jk)=-ze1
        ptenv(jl,jk)=-ze2
!	
! Cx =0 based Cx=/= 0. above
!
!        ptent(jl,jk) = (ze1*um1(jl,jk) + ze2*vm1(jl,jk))/cp1003
!
        dked(jl,jk)=ptent(jl,jk)/zbn2(jl,jk)
        ptent(jl,jk) =  0.0
        if (dked(jl,jk).lt.0)  dked(jl,jk) = dked_min
        enddo
        enddo
!	
! add limiters/efficiency for "unbalanced ics"
!       
	ptenu = gw_eff*ptenu
	ptenv = gw_eff*ptenv
	ptent = gw_eff*ptent		
	dked =  gw_eff* dked	        	
!  
!--------------------------------------------------------------------------- 
!
       if (kdt == 1 .and. mpi_id == master) then
             print *, 'vgw done gwdrag_wms in ifs '
!
        print *, maxval(ptenu)*86400.,  minval(ptenu)*86400, 'vgw ax'
        print *, maxval(ptenv)*86400.,  minval(ptenv)*86400, 'vgw ay'
        print *, maxval(dked)*1.,  minval(dked)*1,  'vgw keddy m2/sec'
        print *, maxval(ptent)*86400.,  minval(ptent)*86400,'vgw eps'
!
        print *, ' ugwp -heating rates '
        endif

        return
        end subroutine fv3_ugwp_betav0
	
	
!=================================================================
!      call gwdps_2017(im, ix, im, levs, model%lonr, me, 
!        kdt, model%nmtvr, kpbl, dtp,    dvdt, dudt, dtdt,   &    
!                 statein%ugrs, statein%vgrs, statein%tgrs,  &
!                 statein%qgrs, statein%prsi, del,           &
!                 statein%prsl, statein%prslk, statein%phii, &
!                 statein%phil,                              &
!                 sfcprop%hprime(1,1), oc, oa4, clx, theta,  &
!                 sigma, gamma, elvmax, dusfcg, dvsfcg,      &
!                 con_g, con_cp, con_rd, con_rv,             &
!                 model%cdmbgwd, lprnt,ipr)
!==========================================================================================================
      subroutine gwdps_2017(im,ix,iy,km, imx, me, 
     &               kdt, nmtvr,kpbl,deltim, 
     &               a,b,c,            u1,v1,t1,q1,              
     &               prsi,del,prsl,prslk,phii, phil,         
     &               hprime,oc,oa4,clx4,theta,sigma,gamma,elvmax,       
     &               dusfc, dvsfc,   g, cp, rd, rv,                    
     &               cdmbgwd, lprnt, ipr)
!
      use machine , only : kind_phys
      implicit none
!
      logical lprnt
      integer ::  im, ix, iy, km
      integer ::  nmtvr
      integer            ::  imx, kdt, ipr, me
      integer            ::  kpbl(im)                 ! index for the pbl top layer!
!
      real(kind=kind_phys) deltim, g, cp, rd, rv,  cdmbgwd(2)
      real(kind=kind_phys) a(iy,km),    b(iy,km),      c(iy,km)
      real(kind=kind_phys)
     &                     u1(ix,km),   v1(ix,km),     t1(ix,km),       
     &                     q1(ix,km),   prsi(ix,km+1), del(ix,km),      
     &                     prsl(ix,km), prslk(ix,km),  phil(ix,km),     
     &                     phii(ix,km+1)
      real(kind=kind_phys) oc(im),     oa4(iy,4), clx4(iy,4), hprime(im)
! 
      real(kind=kind_phys) elvmax(im),theta(im),sigma(im),gamma(im)
      real(kind=kind_phys) wk(im)
      real(kind=kind_phys) bnv2lm(im,km),pe(im),ek(im),zbk(im),up(im)
      real(kind=kind_phys) db(im,km),ang(im,km),uds(im,km)
      real(kind=kind_phys) zlen, dbtmp, r, phiang, cdmb, dbim
      real(kind=kind_phys) eng0, eng1
!
!     some constants
!
      real(kind=kind_phys) pi, dw2min, rimin, ric, bnv2min, efmin
     &,                 efmax,hpmax,hpmin, rad_to_deg, deg_to_rad

      real(kind=kind_phys) frc,    ce,  ceofrc, frmax, cg, gmax
     &,                    veleps, factop, rlolev, rdi
      real(kind=kind_phys) dpmin,hminmt,hncrit,minwnd,sigfac
!
      parameter (pi=3.1415926535897931)
      parameter (rad_to_deg=180.0/pi, deg_to_rad=pi/180.0)
      parameter (dw2min=1., rimin=-100., ric=0.25, bnv2min=1.0e-5) 
      parameter (frc=0.7, ce=0.8, ceofrc=ce/frc, frmax=10., cg=0.5)
      parameter (gmax=1.0, veleps=1.0, factop=0.5)
      parameter (rlolev=50000.0)
!vay
      parameter (efmin=0.0, efmax=5.0)
!vay   
      real, parameter  :: xlinvmax = 6.28e-3/42.                    ! ~4*13 km max-kx  
       
      parameter (hpmax=2400.0, hpmin=75.0)   ! sgo-hprime range
      parameter (hncrit=8000.)   ! max value in meters for elvmax (*j*)
      parameter (sigfac=4.0)     ! mb3a expt test for elvmax factor (*j*)  elvmax+sigfac*sgo_hprime
      parameter (hminmt=50.)     ! min mtn height (*j*)
      parameter (minwnd=0.1)     ! min wind component (*j*)
      parameter (dpmin=5000.0)   ! minimum thickness of the reference layer in pa
!

      integer, parameter :: mdir=8
      integer nwdir(mdir)
      data nwdir/6,7,5,8,2,3,1,4/
      save nwdir
      real(kind=kind_phys) fdir, r_cpdeltim   
!
      logical icrilv(im)
!
!----   mountain induced gravity wave drag
!
      real(kind=kind_phys) taub(im),  xn(im),     yn(im),    ubar(im)     &
     &,                    vbar(im),  ulow(im),   oa(im),    clx(im)      &
     &,                    roll(im),  uloi(im),   dusfc(im), dvsfc(im)    &
     &,                    dtfac(im), xlinv(im),  delks(im), delks1(im)
!
      real(kind=kind_phys) bnv2(im,km),  taup(im,km+1), ri_n(im,km)       &
     &,                    taud(im,km),  ro(im,km),     vtk(im,km)        &
     &,                    vtj(im,km),   scor(im),      velco(im,km-1)    &
     &,                    bnv2bar(im)
!
!
      integer   kref(im), kint(im), iwk(im), ipt(im)
!
      integer   kreflm(im), iwklm(im)
      
      integer   idxzb(im), ktrial, klevm1
!
      real(kind=kind_phys) gor,    gocp,  fv,    gr2,  bnv,  fr           
     &,                    brvf,   cleff, tem,   tem1,  tem2, temc, temv  
     &,                    wdir,   ti,    rdz,   dw2,   shr2, bvf2        
     &,                    rdelks, efact, coefm, gfobnv                   
     &,                    scork,  rscor, hd,    fro,   rim,  sira        
     &,                    dtaux,  dtauy, pkp1log, pklog
      integer kmm1, kmm2, lcap, lcapp1, kbps, kbpsp1,kbpsm1               
     &, kmps, idir, nwd, i, j, k, klcap, kp1, kmpbl, npt            
     &, kmll                                

      real(kind=kind_phys) rgg, mbv1, mbu1, mbu2, mbv2  
     
      real :: hpeff, hdsat
!========================================= start
      rgg= 1.0/g/g
      rdi  = 1.0 / rd
      gor  = g/rd
      gr2  = g*gor
      gocp = g/cp
      fv   = rv/rd - 1
      r_cpdeltim =1./cp/deltim
!
      kmm1   = km - 1
      kmm2   = km - 2
      lcap   = km                ! all range of levels from (1, top_mt) => (km, top lid)
      lcapp1 = lcap + 1
!    
      fdir=mdir/(pi+pi)
!
! res-awareness of gfs         =========================  cdmb = 192.0/float(imx)
!
      cdmb  = 4.0 * 192.0/float(imx)
      cleff = 0.5e-5 / sqrt(float(imx)/192.0) !  this is inverse of cleff!
      if (cdmbgwd(1) >= 0.0) cdmb = cdmb   * cdmbgwd(1)
      if (cdmbgwd(2) >= 0.0) cleff = cleff * cdmbgwd(2)
!=========================                                ====== res-awareness of gfs
     
      do i = 1, im
         dusfc(i) = 0.
         dvsfc(i) = 0.
      enddo
!
      do k = 1, km
        do i = 1, im
          db(i,k)  = 0.
          ang(i,k) = 0.
          uds(i,k) = 0.
!	  
! zero tendencym w/o accumulation
!	  
	  a(i,k) = 0.
	  b(i,k) = 0.
	  c(i,k) = 0.
        enddo
      enddo
!
!
      if ( nmtvr .eq. 14) then 
! ----  for lm and gwd calculation points
        ipt = 0
        npt = 0
        do i = 1,im
          if ( (elvmax(i) .gt. hminmt) 
     &       .and. (hprime(i) .gt. hpmin) )  then
             npt      = npt + 1
             ipt(npt) = i
          endif
        enddo

        if (npt .eq. 0) return     ! no gwd/mb calculation done!
!
! --- iwklm is the level above the height of the of the mountain.
! --- idxzb is the level of the dividing streamline.
! initialize dividing streamline (ds) control vector
!==========================================================================
! part -i    ---  oro mtb 2017   with  options to change  zblk-definition
!==========================================================================
        do i=1,npt
          iwklm(i)  = 2
          idxzb(i)  = 0 
          kreflm(i) = 0
        enddo
        kmll = kmm1
! --- subgrid mountain blocking section
!  (*j*)  11/03:  test upper limit on kmll=km - 1
!      then do not need hncrit -- test with large hncrit first.
!       kmll  = km / 2 ! maximum mtnlm height : # of vertical levels / 2
! --- no mtn should be as high as kmll (so we do not have to start at 
! --- the top of the model but could do calc for all levels).
!
          do i = 1, npt
            j = ipt(i)
            elvmax(j) = min (elvmax(j) + sigfac * hprime(j), hncrit)
          enddo
!
        do k = 1,kmll
          do i = 1, npt
            j = ipt(i)
! --- interpolate to max mtn height for index, iwklm(i) wk[gz]
! --- elvmax is limited to hncrit because to hi res topo30 orog.
            pkp1log =  phil(j,k+1) / g
            pklog =  phil(j,k)   / g
!!!-------     elvmax(j) = min (elvmax(j) + sigfac * hprime(j), hncrit)
            if ( ( elvmax(j) .le.  pkp1log ) .and. 
     &           ( elvmax(j) .ge.   pklog  ) ) then
!  
               wk(i)  = g * elvmax(j) / ( phil(j,k+1) - phil(j,k) )
               iwklm(i)  =  max(iwklm(i), k+1 )    
            endif
! ---        find at prsl levels large scale environment variables   
            vtj(i,k)  = t1(j,k)  * (1.+fv*q1(j,k))
            vtk(i,k)  = vtj(i,k) / prslk(j,k)
            ro(i,k)   = rdi * prsl(j,k) / vtj(i,k)                 ! density kg/m**3
          enddo
        enddo
!
         
        klevm1 = kmll - 1
        do k = 1, klevm1  
          do i = 1, npt
           j   = ipt(i)
            rdz  = g   / ( phil(j,k+1) - phil(j,k) )
            bnv2lm(i,k) = (g+g) * rdz * ( vtk(i,k+1)-vtk(i,k) )
     &                     / ( vtk(i,k+1)+vtk(i,k) )
            bnv2lm(i,k) = max( bnv2lm(i,k), bnv2min )
          enddo
        enddo
!   
!
        do i = 1, npt
          j   = ipt(i)
          delks(i)  = 1.0 / (prsi(j,1) - prsi(j,iwklm(i)))
          delks1(i) = 1.0 / (prsl(j,1) - prsl(j,iwklm(i)))
          ubar (i)  = 0.0
          vbar (i)  = 0.0
          roll (i)  = 0.0
          pe   (i)  = 0.0
          ek   (i)  = 0.0
          bnv2bar(i) = (prsl(j,1)-prsl(j,2)) * delks1(i) * bnv2lm(i,1)
        enddo

! --- find the dividing stream line height 
! --- starting from the level above the max mtn downward
! --- iwklm(i) is the k-index of mtn elvmax elevation
!
!! the maximum mountain height and processing downward.
        do ktrial = kmll, 1, -1
          do i = 1, npt
             if ( ktrial .lt. iwklm(i) .and. kreflm(i) .eq. 0 ) then
                kreflm(i) = ktrial
             endif
          enddo
        enddo
!     print *,' in gwdps_lm.f 4 npt=',npt,kreflm(npt),me
!
! --- in the layer kreflm(i) to 1 find pe (which needs n, elvmax)
! ---  make averages, guess dividing stream (ds) line layer.
! ---  this is not used in the first cut except for testing and
! --- is the vert ave of quantities from the surface to mtn top.
!   
        do i = 1, npt
          do k = 1, kreflm(i)
            j        = ipt(i)
            rdelks     = del(j,k) * delks(i)
            ubar(i)    = ubar(i)  + rdelks * u1(j,k) ! trial mean u below 
            vbar(i)    = vbar(i)  + rdelks * v1(j,k) ! trial mean v below 
            roll(i)    = roll(i)  + rdelks * ro(i,k) ! trial mean ro below 
            rdelks     = (prsl(j,k)-prsl(j,k+1)) * delks1(i)
            bnv2bar(i) = bnv2bar(i) + bnv2lm(i,k) * rdelks
          enddo
        enddo
!     
! --- need the first layer where pe>ek - as soon as 
! --- idxzb is not 0 we have a hit and zb is found.
!
        do i = 1, npt
          j = ipt(i)
          do k = iwklm(i), 1, -1
            phiang   =  atan2(v1(j,k),u1(j,k))*rad_to_deg
            ang(i,k) = theta(j) - phiang 
            if ( ang(i,k) .gt.  90. ) ang(i,k) = ang(i,k) - 180.
            if ( ang(i,k) .lt. -90. ) ang(i,k) = ang(i,k) + 180.
            ang(i,k) = ang(i,k) * deg_to_rad
!
            uds(i,k) = 
     &          max(sqrt(u1(j,k)*u1(j,k) + v1(j,k)*v1(j,k)), minwnd)
! --- test to see if we found zb previously
            if (idxzb(i) .eq. 0 ) then
              pe(i) = pe(i) + bnv2lm(i,k) * 
     &           ( g * elvmax(j) - phil(j,k) ) * 
     &           ( phii(j,k+1) - phii(j,k) ) * rgg
! 
! --- kinetic energy is at the layer zb
! --- theta ranges from -+90deg |_ to the mtn "largest topo variations"
              up(i)  =  uds(i,k) * cos(ang(i,k))
              ek(i)  = 0.5 *  up(i) * up(i) 

! --- dividing stream lime  is found when pe =exceeds ek.
              if ( pe(i) .ge.  ek(i) ) idxzb(i) = k

! --- then mtn blocked flow is between zb=k(idxzb(i)) and surface
            endif
          enddo
        enddo
!
        do i = 1, npt
          j    = ipt(i)
! --- calc if n constant in layers (zb guess) - a diagnostic only.
          zbk(i) = elvmax(j)
     &           - sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i))/bnv2bar(i)
        enddo
!
! --- the drag for mtn blocked flow
! 
        do i = 1, npt
          j = ipt(i)
          if ( idxzb(i) .gt. 0 ) then 
            do k = idxzb(i), 1, -1
!
              if ( phil(j,idxzb(i)) .gt.  phil(j,k) ) then
!
              zlen = sqrt( ( phil(j,idxzb(i)) - phil(j,k) ) / 
     &                       ( phil(j,k ) + g * hprime(j) ) )
! --- lm eq 14:
                r = (cos(ang(i,k))**2 + gamma(j) * sin(ang(i,k))**2) / 
     &              (gamma(j) * cos(ang(i,k))**2 + sin(ang(i,k))**2)
!
!! the drag per unit area and per unit height is written (eq.15 in 
!! lott and miller (1997)
!
                dbtmp = 0.25 *  cdmb *
     &                  max( 2. - 1. / r, 0. ) * sigma(j) *
     &                  max(cos(ang(i,k)), gamma(j)*sin(ang(i,k))) *
     &                  zlen / hprime(j) 
                db(i,k) =  dbtmp * uds(i,k)    
!
              endif
            enddo
          endif
        enddo
! 
!.............................
!.............................
! end  mtn blocking section
!
      elseif ( nmtvr .ne. 14) then 
! ----  for mb not present and  gwd (nmtvr .ne .14) 
        ipt     = 0
        npt     = 0
        do i = 1,im
          if ( hprime(i) .gt. hpmin )  then
             npt      = npt + 1
             ipt(npt) = i
          endif
        enddo
	
        if (npt .eq. 0) return     ! no gwd/mb calculation done!
        do i=1,npt
          idxzb(i) = 0
        enddo
      endif

!=============================================================================
! part-ii    --- orographic gravity wave drag section with 
!
!            traditional "linear saturation" w/o rim and palmer's 1987 solver
!
!=============================================================================

      kmpbl  = km / 2 ! maximum pbl height : # of vertical levels / 2
!
      do k = 1,km
        do i =1,npt
          j         = ipt(i)
          vtj(i,k)  = t1(j,k)  * (1.+fv*q1(j,k))
          vtk(i,k)  = vtj(i,k) / prslk(j,k)
          ro(i,k)   = rdi * prsl(j,k) / vtj(i,k) ! density tons/m**3
          taup(i,k) = 0.0
        enddo
      enddo

      do k = 1,kmm1
        do i =1,npt
          j         = ipt(i)
          ti        = 2.0 / (t1(j,k)+t1(j,k+1))
          tem       = ti  / (prsl(j,k)-prsl(j,k+1))
          rdz       = g   / (phil(j,k+1) - phil(j,k))
          tem1      = u1(j,k) - u1(j,k+1)
          tem2      = v1(j,k) - v1(j,k+1)
          dw2       = tem1*tem1 + tem2*tem2
          shr2      = max(dw2,dw2min) * rdz * rdz
          bvf2      = g*(gocp+rdz*(vtj(i,k+1)-vtj(i,k))) * ti
          ri_n(i,k) = max(bvf2/shr2,rimin)   ! richardson number
!                                              brunt-vaisala frequency
!
          bnv2(i,k) = (g+g) * rdz * (vtk(i,k+1)-vtk(i,k))
     &                            / (vtk(i,k+1)+vtk(i,k))
          bnv2(i,k) = max( bnv2(i,k), bnv2min )
        enddo
      enddo
! 
!
      do i=1,npt
        iwk(i) = 2
      enddo
      do k=3,kmpbl
        do i=1,npt
          j   = ipt(i)
          tem = (prsi(j,1) - prsi(j,k))
          if (tem .lt. dpmin) iwk(i) = k
        enddo
      enddo
!
!! kpbl is the index for the pbl top layer.
      kbps = 1
      kmps = km
      do i=1,npt
        j         = ipt(i)
        kref(i)   = max(iwk(i), kpbl(j)+1 ) ! reference level 
        delks(i)  = 1.0 / (prsi(j,1) - prsi(j,kref(i)))
        delks1(i) = 1.0 / (prsl(j,1) - prsl(j,kref(i)))
        ubar (i)  = 0.0
        vbar (i)  = 0.0
        roll (i)  = 0.0
        kbps      = max(kbps,  kref(i))
        kmps      = min(kmps,  kref(i))
!
        bnv2bar(i) = (prsl(j,1)-prsl(j,2)) * delks1(i) * bnv2(i,1)
      enddo
!      print *,' in gwdps_lm.f gwd:15  =',kbps,kmps
      kbpsp1 = kbps + 1
      kbpsm1 = kbps - 1
!
      do k = 1,kbps
        do i = 1,npt
          if (k .lt. kref(i)) then
            j          = ipt(i)
            rdelks     = del(j,k) * delks(i)
            ubar(i)    = ubar(i)  + rdelks * u1(j,k)   ! mean u below kref
            vbar(i)    = vbar(i)  + rdelks * v1(j,k)   ! mean v below kref
!
            roll(i)    = roll(i)  + rdelks * ro(i,k)   ! mean ro below kref
            rdelks     = (prsl(j,k)-prsl(j,k+1)) * delks1(i)
            bnv2bar(i) = bnv2bar(i) + bnv2(i,k) * rdelks
          endif
        enddo
      enddo
!      print *,' in gwdps_lm.f gwd:15b =',bnv2bar(npt)
!
!     figure out low-level horizontal wind direction and find 'oa'
!
!             nwd  1   2   3   4   5   6   7   8
!              wd  w   s  sw  nw   e   n  ne  se
!
!> - calculate low-level horizontal wind direction, the derived 
!! orographic asymmetry parameter (oa), and the derived lx (clx). 
      do i = 1,npt
        j      = ipt(i)
        wdir   = atan2(ubar(i),vbar(i)) + pi
        idir   = mod(nint(fdir*wdir),mdir) + 1
        nwd    = nwdir(idir)
        oa(i)  = (1-2*int( (nwd-1)/4 )) * oa4(j,mod(nwd-1,4)+1)
        clx(i) = clx4(j,mod(nwd-1,4)+1)
      enddo
!
!-----xn,yn            "low-level" wind projections in zonal
!                                    & meridional directions
!-----ulow             "low-level" wind magnitude -        (= u)
!-----bnv2             bnv2 = n**2
!-----taub             base momentum flux
!-----= -(ro * u**3/(n*xl)*gf(fr) for n**2 > 0
!-----= 0.                        for n**2 < 0
!-----fr               froude    =   n*hprime / u
!-----g_ph_factor      gmax*fr**2/(fr**2+cg/oc)
!
!-----initialize some arrays
!
      do i = 1,npt
        xn(i)     = 0.0
        yn(i)     = 0.0
        taub (i)  = 0.0
        ulow (i)  = 0.0
        dtfac(i)  = 1.0
        icrilv(i) = .false. ! initialize critical level control vector
        
!
!----compute the "low level" wind magnitude (m/s)
!
        ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0)
        uloi(i) = 1.0 / ulow(i)
      enddo
!
      do  k = 1,kmm1
        do  i = 1,npt
          j            = ipt(i)
          velco(i,k)   = 0.5 * ((u1(j,k)+u1(j,k+1))*ubar(i)
     &                       +  (v1(j,k)+v1(j,k+1))*vbar(i))
          velco(i,k)   = velco(i,k) * uloi(i)
!     
        enddo
      enddo
!      
!
!
!  warning  kint = kref !!!!!!!!!
      do i=1,npt
        kint(i) = kref(i)
      enddo
!
      do i = 1,npt
        j      = ipt(i)
        bnv    = sqrt( bnv2bar(i) )
	hpeff = min(hprime(j),hpmax)
        fr     = bnv     * uloi(i) * hpeff
        fr     = min(fr, frmax)
        xn(i)  = ubar(i) * uloi(i)
        yn(i)  = vbar(i) * uloi(i)
!
!     compute the base level stress and store it in taub
!     calculate enhancement factor, number of mountains & aspect
!     ratio const. use simplified relationship between standard
!     deviation & critical hgt
!
! - calculate enhancement factor (e),number of mountans (m') and 
!! aspect ratio constant.
! as in eq.(4.9),(4.10),(4.11) in kim and arakawa (1995) 
!
        efact    = (oa(i) + 2.) ** (ceofrc*fr)
        efact    = min( max(efact,efmin), efmax )
        coefm    = (1. + clx(i)) ** (oa(i)+1.)
        xlinv(i) = min (xlinvmax, coefm * cleff)     ! does not exceed  42km ~4*dx
        tem      = fr    * fr * oc(j)
        gfobnv   = gmax  * tem / ((tem + cg)*bnv)  ! g/n0
!	
! source flux tau0
        taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)
     &           * ulow(i)  * gfobnv  * efact 
        if(fr.le.frc)
     &	 taub(i)= xlinv(i)*roll(i)*bnv*ulow(i)*hpeff*hpeff        
!
        k        = max(1, kref(i)-1)
        tem      = max(velco(i,k)*velco(i,k), 0.1)
        scor(i)  = bnv2(i,k) / tem                           ! scorer parameter below ref level
      enddo
!                                   
!----set up bottom values of stress
!
      do k = 1, kbps
        do i = 1,npt
          if (k .le. kref(i)) taup(i,k) = taub(i)
        enddo
      enddo
!======================    from here compute taup(z) ====================
!
!    along main-ridge
!
!=======================================================================
!   now compute vertical structure of the stress.
!
      do k = kmps, kmm1                   ! vertical level k loop!

        kp1 = k + 1
        do i = 1, npt
!
!-----unstable layer if ri < ric
!---- at (u-c)=0. crit layer exists and bit vector should be set (.le.)
!
          if (k .ge. kref(i)) then
            icrilv(i) = icrilv(i) .or. ( ri_n(i,k) .lt. ric)
     &                            .or. (velco(i,k) .le. 0.0)
           endif
        enddo
!
        do i = 1,npt
          if (k .ge. kref(i))   then
            if (.not.icrilv(i) .and. taup(i,k) .gt. 0.0 ) then
              temv = 1.0 / max(velco(i,k), 0.01)
!        
              if (oa(i).gt.0. .and. kp1 .lt. kint(i)) then
                scork   = bnv2(i,k) * temv * temv
                rscor   = min(1.0, scork / scor(i))
                scor(i) = scork
              else 
                rscor   = 1.
              endif
!

              brvf = sqrt(bnv2(i,k))                                 ! brunt-vaisala frequency
              tem1 = xlinv(i)*(ro(i,kp1)+ro(i,k))*brvf*0.5           ! linear flux-mult w/o hp*hp
     &                       * max(velco(i,k),0.01)
!
! ??? hp(z)
!     
              hd   = sqrt(taup(i,k) / tem1)
              fro  = brvf * hd * temv                        !reverse froude
	      if (fro .le. frc) then
	        hdsat = frc/brvf/temv
	        taup(i,kp1) =	tem1 *hdsat *hdsat           ! effective dissip-n kdis can be added later 
	      else
	        taup(i,kp1) =taup(i,k)
	      endif
!	      
	      goto 888                                       ! skip palmer/shutts ri-solvers
!
!    rim is the  minimum-richardson number by shutts (1985), palmer et al.(1986)
! see eq.(4.6) in kim and arakawa (1995) \cite kim_and_arakawa_1995.
!
              tem2   = sqrt(ri_n(i,k))
              tem    = 1. + tem2 * fro
              rim    = ri_n(i,k) * (1.-fro) / (tem * tem)
!
!    check stability to employ non-standard 'saturation hypothesis'
!    -  this is a  nwp approach of palmer et al. (1986) 
!
      if (rim .le. ric .and. (oa(i) .le. 0. .or.  kp1.ge.kint(i))) then
                 temc = 2.0 + 1.0 / tem2
                 hd   = velco(i,k) * (2.*sqrt(temc)-temc) / brvf
                 taup(i,kp1) = tem1 * hd * hd
      else 
                 taup(i,kp1) = taup(i,k) * rscor
      endif
      
888   continue       !    skip palmer/shutts ri-solvers    for ogw-propagation
      
              taup(i,kp1) = min(taup(i,kp1), taup(i,k))
            endif
          endif
        enddo
      enddo                   ! vertical level k loop!
!
!     do i=1,im
!       taup(i,km+1) = taup(i,km)
!     enddo
!
      if(lcap .le. km) then
         do klcap = lcapp1, km+1
            do i = 1,npt
              sira          = prsi(ipt(i),klcap) / prsi(ipt(i),lcap)
              taup(i,klcap) = sira * taup(i,lcap)   ! reduction by p(k)/p(lcap) < 1
            enddo
         enddo
      endif

!======================================================================
!     calculate - (g/p*)*d(tau)/d(sigma) and decel terms dtaux, dtauy
!      indep-nt on solver type 
!======================================================================
      do k = 1,km
        do i = 1,npt
          taud(i,k) = g * (taup(i,k+1) - taup(i,k)) / del(ipt(i),k)
        enddo
      enddo
!
!------limit de-acceleration (momentum deposition ) at top to 1/2 value
!
      do klcap = lcap, km
         do i = 1,npt
            taud(i,klcap) = taud(i,klcap) * factop    ! factop=0.5
         enddo
      enddo
!
!------if the gravity wave drag would force a critical line in the
!------layers below sigma=rlolev during the next deltim timestep,
!------then only apply drag until that critical line is reached.
!
      do k = 1,kmm1
        do i = 1,npt
           if (k .gt. kref(i) .and. prsi(ipt(i),k) .ge. rlolev) then
             if(taud(i,k).ne.0.) then
               tem = deltim * taud(i,k)
               dtfac(i) = min(dtfac(i),abs(velco(i,k)/tem))
             endif
           endif
        enddo
      enddo
!
! -orodrag effects: a, b, c, +  dusfc, dvsfc (see parameter description).
! 
      do k = 1,km
        do i = 1,npt
          j          = ipt(i)
          taud(i,k)  = taud(i,k) * dtfac(i)
          dtaux      = taud(i,k) * xn(i)
          dtauy      = taud(i,k) * yn(i)
          eng0       = 0.5*(u1(j,k)*u1(j,k)+v1(j,k)*v1(j,k))
!
! ---  lm mb (*j*)  changes overwrite gwd
!
          if ( k .lt. idxzb(i) .and. idxzb(i) .ne. 0 ) then
!
! here dbim*u = krayleigh*u => compute updated drag along "main"-axis
!      and project it back on 2-directions (x-y)
!
            dbim = db(i,k) / (1.+db(i,k)*deltim)
            mbv1 = - dbim * v1(j,k)
            mbu1 = - dbim * u1(j,k)
            a(j,k)  = mbv1 + a(j,k)
            b(j,k)  = mbu1 + b(j,k)
!
!explicit bettet to use implicit eng0/(1.0+dbim*deltim)/(1.0+dbim*deltim)
!
            eng1    = eng0*(1.0-dbim*deltim)*(1.0-dbim*deltim)
! 
            dusfc(j)   = dusfc(j) + mbu1 * del(j,k)
            dvsfc(j)   = dvsfc(j) + mbv1 * del(j,k)
          else
!
            a(j,k)     = dtauy     + a(j,k)
            b(j,k)     = dtaux     + b(j,k)
            eng1       = 0.5*(
     &                   (u1(j,k)+dtaux*deltim)*(u1(j,k)+dtaux*deltim)
     &                 + (v1(j,k)+dtauy*deltim)*(v1(j,k)+dtauy*deltim))
            dusfc(j)   = dusfc(j)  + dtaux * del(j,k)
            dvsfc(j)   = dvsfc(j)  + dtauy * del(j,k)
          endif
!
! compute "heating"-only depositions
!
          c(j,k) = c(j,k) + max(eng0-eng1,0.)*r_cpdeltim
        enddo
      enddo
! ========================> -1/gravity factor    tauxgw tauygw n/m/m
      tem    = -1.0/g
      do i = 1,npt
        j          = ipt(i)
        dusfc(j) = tem * dusfc(j)
        dvsfc(j) = tem * dvsfc(j)
      enddo
!                                                                       
!                   line-763
!
      return
      end

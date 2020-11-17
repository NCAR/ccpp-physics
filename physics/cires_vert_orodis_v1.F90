module cires_vert_orodis_v1


contains


! subroutine ugwp_drag_mtb
! subroutine ugwp_taub_oro  
! subroutine ugwp_oro_lsatdis  
!
     subroutine ugwp_drag_mtb( iemax,  nz,                              &
        elvpd, elvp, hprime , sigma, theta, oc,  oa4, clx4, gam, zpbl,  &
        up, vp, tp, qp, dp, zpm, zpi, pmid, pint, idxzb, drmtb,taumtb)

      use ugwp_common_v1,  only  : bnv2min, grav, grcp, fv, rad_to_deg, dw2min, velmin, rdi
      use ugwp_oro_init_v1, only : nridge, cdmb, fcrit_mtb, frmax, frmin, strver
      
      implicit none   
!========================
! several versions for drmtb => high froude mountain blocking
! version 1 => vay_2018 ; 
! version 2 => kdn_2005 ; Kim & Doyle in NRL-2005
! version 3 => ncep/gfs-2017 -gfs_2017 with lm1997
!======================== 
!     real, parameter   ::  Fcrit_mtb  = 0.7   
      
      integer, intent(in)              ::  nz     
      integer, intent(in)              ::  iemax    ! standard ktop  z=elvpd  + 4 * hprime
      real  , intent(out)              ::  taumtb    
          
      integer , intent(out)            ::  idxzb           
      real, dimension(nz), intent(out) ::  drmtb

      real, intent(in)  ::  elvp, elvpd              !elvp  = min (elvpd  + sigfac * hprime(j), hncrit=10000meters)
      real, intent(in)  ::  hprime , sigma, theta, oc,  oa4(4), clx4(4), gam
      real, intent(in)  ::  zpbl     

      real, dimension(nz),    intent(in) :: up, vp, tp, qp, dp, zpm, pmid
      real, dimension(nz+1),  intent(in) :: zpi, pint

      ! character(len=*), intent(out) :: errmsg
      ! integer,          intent(out) :: errflg
!
      real, dimension(nz+1)              :: zpi_zero
      real, dimension(nz)                :: zpm_zero
      real ::  vtj, rhok, bnv2, rdz, vtkp, vtk, dzp
      
      real, dimension(nz) ::  bn2, uds, umf, cosang, sinang
      
      integer :: k,  klow, ktop, kpbl
      real    :: uhm, vhm, bn2hm, rhohm,                  &
                 mtb_fix, umag, bnmag, frd_src,           &
                 zblk, who_iz_normal, rlm97,              &
                 phiang, ang, pe, ek,                     &
                 cang, sang, ss2, cs2, zlen,  dbtmp,      &
                 hamp, bgamm, cgamm      
 

    ! Initialize CCPP error handling variables
    ! errmsg = ''
    ! errflg = 0
 
!==================================================
!
!     elvp + hprime <=>elvp + nridge*hprime, ns =2
!                 ns = sigfac
!     tau_parel & tau_normal along major "axes" 
!
!  options to block the "flow", choices for [klow, ktop]
!
!  1-directional (normal)  & 2-directional "blocking"
!
!==================================================
! no - blocking: drmtb(1:nz) = 0.0
!=================
      idxzb = -1
      drmtb(1:nz) = 0.0
      taumtb      = 0.0
      klow = 2

      ktop = iemax
      hamp = nridge*hprime
      
! reminder:        cdmb = 4.0 * 192.0/float(imx)*cdmbgwd(1)   Lellipse= a/2=sigma/hprime

      mtb_fix = cdmb*sigma/hamp        !hamp  ~ 2*hprime and 1/sigfac = 0.25 is inside 1/hamp

      ! if (mtb_fix == 0.) then
      !    write(errmsg,'(*(a))') cdmb, sigma, hamp, ' MTB == 0'
      !    errflg = 1
      !    return
      ! endif

      if (strver == 'vay_2018')  then    

        zpm_zero = zpm - zpi(1)
        zpi_zero = zpi - zpi(1)

         do k=1, nz-1
           if (hamp .le. zpi_zero(k+1) .and. (hamp .gt. zpi_zero(k)  ) ) then
               ktop  =  k+1      !......simply k+1 next interface level
               exit
            endif
          enddo
!        print *,  klow, ktop,   ' klow-ktop '
        call um_flow(nz, klow, ktop, up, vp, tp, qp, dp, zpm, zpi, pmid, pint, &
             bn2, uhm, vhm, bn2hm, rhohm)   
     
        umag = max(sqrt(uhm*uhm + vhm*vhm), velmin)       !velmin=dw2min =1.0 m/s 
        ! if (bn2hm .le. 0.0) then
        !    write(errmsg,'(*(a))') 'unstable MF for MTB  - RETURN '
        !    errflg = 1
        !    return    ! unstable PBL
        ! end if

        bnmag =sqrt(bn2hm)

        frd_src  =  min(hamp*bnmag/umag, frmax)        ! frmax =10.

!       print *, frd_src,  Fcrit_mtb/frd_src,  ' no-Blocking >  1 '        	          
!      
        if ( frd_src .le. Fcrit_mtb)  RETURN   ! no-blocking, although on small ridges  with weak winds can be blocking
!
!  zblk > 0
!  Fcrit_mtb > Fcrit_ogw    h_clip = Fr_mtb*U/N ! h_hill  minus h_clip = zblk
!
        zblk = hamp*(1. - Fcrit_mtb/frd_src)   
        idxzb =1
        do k = 2, ktop

          if ( zblk < zpm_zero(k) .and. zblk >= zpm_zero(k-1)) then 
             idxzb = k 
             exit
          endif
        enddo  
!   
        if (idxzb == 1)  RETURN  ! first surface level block is not "important"

        if (idxzb > 1) then      ! let start with idxzb = 2....and up with LM1997
!
! several options to compute MTB-drag: a) IFS_1997 ; b) WRF_KD05 ; c) SJM_2000
!   
          bgamm = 1.0 - 0.18*gam -0.04*gam*gam
          cgamm = 0.48*gam +0.3*gam*gam   
        
          do k = 1, idxzb-1  
            zlen = sqrt( (zblk - zpm_zero(k) ) /  ( zpm_zero(k) +hprime ))

            umag =  max(sqrt(up(k)*up(k) + vp(k)*vp(k)), velmin)

            phiang   =  atan(vp(k)/umag)
! theta -90/90
            ang  = theta - phiang
            cang = cos(ang)  ; sang = sin(ang)  
              
            who_iz_normal = max(cang, gam*sang )         !gfs-2018

            cs2 = cang* cang ; ss2 = 1.-cs2

            rlm97 =(gam * cs2 + ss2)/  (cs2 + gam * ss2) ! ... (cs2 + gam * ss2) / (gam * cs2 + ss2)   ! check it
!
            if (rlm97 > 2.0 ) rlm97 = 2.0                ! zero mtb-friction at this level
!

            who_iz_normal = bgamm*cs2 + cgamm*ss2     ! LM1997/IFS

            dbtmp = mtb_fix* max(0.,  2.- rlm97)*zlen*who_iz_normal 
            if (dbtmp < 0) dbtmp  = 0.0
!
! several approximation can be made to implement MTB-drag 
!   as a "nonlinear level dependent"-drag or "constant"-drag
!   uds(k) == umag = const between the 1-layer and idxzb
!
               
            drmtb(k) =  dbtmp * abs(umag)  ! full mtb-drag =   -drmtb(k) * uds = -kr*u 
            taumtb = taumtb - drmtb(k)*umag *rdi * pmid(k)/tp(k)*(zpi(k+1)-zpi(k))
! 
! 2-wave appr for anisotropic drmtb_Bellipse(k) and drmtb_Aell(k) can be used 
!   with Umag-projections on A & B ellipse axes
!   mtb_fix =0.25*cdmb*sigma/hprime,  
!   in SM-2000 mtb_fix~ 1/8*[cdmb_A, cdmb_B]*sigma/hprimesum ( A+B)  = 1/4.
!          
!333       format(i4, 7(2x, F10.3))         
!                write(6,333) , k, zpm_zero(k),  zblk, hamp*Fcrit_mtb/frd_src, taumtb*1.e3,  drmtb(k) ,  -drmtb(k)*up(k)*1.e5
          enddo
!   
        endif     
      endif          ! strver=='vay_2018'
!
!
!
      if (strver == 'kdn_2005' .or. strver == 'wrf_2018' ) then
      
        print *, ' kdn_2005  with # of hills '
!
! compute flow-blocking stress based on WRF  'gwdo2d'
!    
      endif
!
!
      if (strver == 'gfs_2018') then

        ktop = iemax;        klow = 2 

        call um_flow(nz, klow, ktop, up, vp, tp, qp, dp, zpm, zpi, pmid, pint,  &
                     bn2, uhm, vhm, bn2hm, rhohm)
        if (bn2hm <= 0.0) RETURN                        ! unstable PBL   
!--------------------------------------------- 	   
!
!'gfs_2018'  .... does not rely on Fr_crit
!            and Fr-regimes 
!----gfs17 for mtn ignores "averaging of the flow"
!   for MTB-part  it is only works with "angles"
!   no projections on [uhm, vhm] -direction
!   kpbl can be used for getting high values of iemax-hill
!----------------------------------------------------------- 
        zpm_zero = zpm - zpi(1)
        zpi_zero = zpi - zpi(1)
        do k=1, nz-1
          if (zpbl .le. zpm_zero(k+1) .and. (zpbl .ge. zpm_zero(k)  ) ) then
            kpbl  =  k+1   
            exit
          endif
        enddo

        do k = iemax, 1, -1
            
          uds(k)    =  max(sqrt(up(k)*up(k) + vp(k)*vp(k)), velmin)
          phiang    =  atan(vp(k)/uds(k))
          ang       = theta - phiang
          cosang(k) = cos(ang)
          sinang(k) = sin(ang)

          if (idxzb == 0) then
            pe     = pe + bn2(k) * (elvp - zpm(k)) *(zpi(k+1) - zpi(k)) 
            umf(k) = uds(k) * cosang(k)                  ! normal to main axis
            ek     = 0.5 * umf(k) * umf(k) 
!
! --- dividing stream lime  is found when pe =>exceeds ek first from the "top"
!
            if (pe >= ek) idxzb = k
            exit
          endif
        enddo

!       idxzb = min(kpbl, idxzb)
!
!  
!
! last:   mtb-drag
!
        if (idxzb > 1) then
          zblk = zpm(idxzb)
          print *, zpm(idxzb)*1.e-3, ' mtb-gfs18 block-lev km ', idxzb, iemax, int(elvp)
          do k = idxzb-1, 1, -1      
!
            zlen = sqrt( (zblk - zpm_zero(k) ) /  ( zpm_zero(k) +hprime ))
            cs2 = cosang(k)* cosang(k)
            ss2 = 1.-cs2
            rlm97 =(gam * cs2 + ss2)/  (cs2 + gam * ss2) ! (cs2 + gam * ss2) / (gam * cs2 + ss2)   ! check it

            who_iz_normal = max(cosang(k), gam*sinang(k))
!
! high res-n higher mtb  0.125 => 3.5 ; (negative of db -- see sign at tendency)
!
            dbtmp = mtb_fix* max(0.,  2.- rlm97)*zlen*who_iz_normal 

            drmtb(k) =  dbtmp * abs(uds(k))  ! full mtb-drag =   -drmtb(k) * uds = -kr*u 
!
            taumtb = taumtb - drmtb(k) * uds(k) *rdi * pmid(k)/tp(k)*(zpi(k+1)-zpi(k))
!            
          enddo
        endif   
      endif                          !  strver=='gfs17'
!                        
!	     
      end subroutine ugwp_drag_mtb
! 
!
! ugwp_taub_oro  - Computes  [taulin,  taufrb,  drlee(levs) ]
!     
! 
      subroutine ugwp_taub_oro(levs, izb, kxw, tau_izb,  fcor,                 &
          hprime , sigma, theta, oc,  oa4, clx4, gamm,                         &
          elvp, up, vp, tp, qp, dp, zpm, zpi, pmid, pint, xn, yn, umag,        &
          tautot, tauogw,  taulee, drlee, tau_src, kxridge, kdswj, krefj, kotr)
!	  
      use ugwp_common_v1,       only : bnv2min, grav, pi, pi2, dw2min, velmin
      use ugwp_common_v1,       only : mkz2min, mkzmin     
      use cires_ugwp_module_v1, only : frcrit, ricrit, linsat
      use ugwp_oro_init_v1,     only :  hpmax, cleff, frmax
      use ugwp_oro_init_v1,     only :  nwdir, mdir, fdir    
      use ugwp_oro_init_v1,     only :  efmin, efmax , gmax, cg, ceofrc  
      use ugwp_oro_init_v1,     only :  fcrit_sm, fcrit_gfs, frmin, frmax
      use ugwp_oro_init_v1,     only :  coro, nridge, odmin, odmax
      use ugwp_oro_init_v1,     only :  strver
! 
      use ugwp_oro_init_v1,     only : zbr_pi
! ---
! 
! define oro-GW fluxes: taulin,  taufrb amd if kdswj > 0 (LWB-lee wave breaking)
!                                approximate for drlee-momentum tendency
! --- 
      implicit none
!      
      integer, intent(in)  ::       levs,  izb
      real   , intent(in)  ::       tau_izb       ! integrated (1:izb) drag -Kr_mtb*U, or Zero  
      integer, intent(out) ::       kdswj, krefj, kotr
      integer              ::       klwb
      real, intent(in)     ::       kxw, fcor
      real, intent(in)     ::       hprime, sigma, theta, oc, gamm, elvp

!
      real, intent(in)     ::   oa4(4), clx4(4)  

      real, dimension(levs),   intent(in)     ::  up, vp, tp, qp, dp
      real, dimension(levs+1), intent(in)     ::  zpi,  pint
      real, dimension(levs  ), intent(in)     ::  zpm,  pmid  
!         
      real,dimension(levs),   intent(out)     ::  drlee
      real,dimension(levs+1), intent(out)     ::  tau_src            
!      
      real, intent(out)                       ::  tauogw, tautot, taulee
      real                                    ::  taulin, tauhcr, taumtb
      real, intent(out)                       ::  xn, yn, umag, kxridge 
! 
!  
! locals
! four possible versions to compute "taubase as  a function of Fr-number"
!   character :: strver='smc_2000'  ! 'kd_2005', 'gfs_2017', 'vay_2018'
!
      real, dimension(levs+1) :: zpi_zero

      real    ::   oa,    clx,  odir, cl4p(4),   clxp 
   
      real    ::   uhm, vhm, bn2hm, rhohm, bnv  

      real    ::   elvpMTB, wdir
      real    ::   tem, efact, coefm, kxlinv, gfobnv
      
      real    ::   fr,  frlin,  frlin2, frlin3, frlocal, dfr
      real    ::   betamax, betaf, frlwb, frmtb
      integer ::   klow, ktop, kph
 
      integer ::   i, j, k, nwd, ind4, idir
      
      real                 ::   sg_ridge, kx2, umd2
      real                 ::   mkz, mkz2, zbr_mkz, mkzi

      real                 ::   hamp      ! clipped hprime*elvmax/elv_clip > hprime
      real                 ::   hogw      ! hprime or hamp for free-prop OGWs z > z(krefj)
      real                 ::   hdsw      ! empirical like DNS amplitudes for Lee-dsw trapped waves
      real                 ::   hcrit
      real                 ::   hblk      ! blocking div-stream height
 
      real                 ::   coef_h2, frnorm
     
 
      real, dimension(levs)   ::  bn2
      real                    ::  rho(levs)
      real, dimension(levs+1) ::  ui, vi, ti,  bn2i, bvi, rhoi
      real, dimension(levs+1) ::  umd, phmkz
      real                    ::  c2f2, umag2, dzwidth, udir
      real                    ::  hogwi, hdswi, hogwz, hdswz         ! height*height wave-amp  
      real                    ::  uogwi, udswi, uogwz, udswz         ! wind2   wave-rms
      real, dimension(levs+1) ::  dtrans, deff
      real                    ::  pdtrans
      logical                 ::  do_klwb_phase = .false.   ! phase-crireria for LLWB of SM00
      logical                 ::  do_dtrans     = .true.   ! dissipative saturation to deposit momentum
                                                           ! between  ZMTB => ZHILL
!-----------------------------------------------------------------------------
!
!  downslope/lee/GW wave regimes kdswj: between ZMTB and ZOGW(krefj)
!                        ZMTB <  ZOGW = ns*HPRIME  < ELVP
!  define krefj as a level for OGWs above ZMTB and "2-3-4*hprime" + ZMTB
!  we rely on the concept of the "CLIPPED-SG" mountain above ZMTB & new 
!  inverse Froude number for the "mean flow" averaged from ZMTB to ZOGW
!  here we can use "elvp" as only for hprime adjustment ...elvp/elvp_MTB
!
!"empirical" specification of tauwave = taulee+tauogw in [ZMTB : ns*HPRIME]
! can be based on numerical runs like WRF-model 
!  for  Frc < Fr< [Frc : 2.5-3 Frc]
!  see suggestions proposed in SM-2000 and Eckermann et al. (2010) 
!-----------------------------------------------------------------------------
      tautot = 0.  ; taulin = 0.  ; taulee = 0.  ;  drlee(1:levs) = 0. ; tau_src = 0.0  
      krefj  = 1   ; kotr = levs+1; kdswj = 1
      xn = 1.0     ; yn = 0.      ; umag = velmin;  kxridge = kxw 

      dtrans  = 0.  ; deff =0.
      klow    = 2 
      elvpMTB = elvp
!
! clipped mountain  H-zmtb for estimating wave-regimes new Fr and MF  above ZMTB 
!
      if (izb > 0 ) then 
         klow = izb
         elvpMTB = max(elvp - zpi(izb), 0.0) 
      endif
      if (elvpMTB <=0 ) print *, ' blocked flow '
      if (elvpMTB <=0 ) return  ! "blocked flow" from the surface to elvMAX	 
     
      zpi_zero(:) = zpi(:) - zpi(1) 
      hblk =  zpi_zero(klow)

      sg_ridge = max( nridge*hprime * (elvp/elvpMTB), hblk+hprime*0.333)

!
! enhance sg_ridge by  elvp/elvpMTB >1 and H_clip = H-hiilnew - zblk  later  for hamp 
!      
      sg_ridge = min(sg_ridge, hpmax)

!       print *, 'sg_ridge ', sg_ridge
     
       do k=1, levs
         if (sg_ridge .gt. zpi_zero(k) .and. ( sg_ridge .le. zpi_zero(k+1)  ) ) then
           ktop  = k+1 
           exit
          endif
       enddo
        
        krefj = ktop              ! the mountain top index for sg_ridge = ns*hprime 

!        if ( izb > 0 .and. krefj .le. izb) then 
!          print *, izb, krefj,  sg_ridge,  zpi_zero(izb), ' izb >ktop '   
!        endif
 
!
! here ktop displays sg_ridge-position not elvP !!!! klow =2  to avoid for 127-126L
! instability due to extreme "thin" layer...128L-model needs cruder vertical resolution
!	     
      call um_flow(levs, klow, ktop, up, vp, tp, qp, dp, zpm, zpi, pmid, pint, &
          bn2, uhm, vhm, bn2hm, rhohm)

      call get_unit_vector(uhm, vhm, xn, yn, umag)
      
      if (bn2hm <= 0.0) RETURN        !  "unstable/neutral" hill need different treatment
      bnv  = sqrt(bn2hm)
      hamp = sg_ridge-zpi_zero(klow)     ! hamp >= nridge*hprime due higher SG-elevations - zblk or first layer
      hogw = hamp
      hdsw = hamp
      

      fr  = bnv * hamp /umag
      fr  = min(fr, frmax)     
      kxridge = max(sigma/hamp, kxw)   ! to get rid from "SSO-errors" kxw-provides max-value for kx
      kx2     = kxridge*kxridge
      umag    = max( umag,  velmin)
      c2f2    = fcor*fcor/kx2
      umag2   = umag*umag - c2f2

      if (umag2 <= 0.0)  RETURN      ! Coriolis cut-off at high-mid latitudes for low kx

      mkz2 = bn2hm/umag2 - kx2         ! we add Coriolis corrections for crude model resolutions "low-kx"
                                       ! and non-stationary waves coro, fcor for small umag
                                       ! bn2hm/[(coro-umag)^2 -fc2/kx2] - kx2, cf = fc/kx => 2 m/s to 11 m/s for 60deg 
      IF (mkz2 <  mkz2min .or. krefj <= 2 )  THEN
!
! case then no effects of wave-orography
!
        krefj  = 1  ; kdswj = 1; kotr = levs ; klwb = 1
        tautot = 0. 
        tauogw = 0. 
        taulee = 0.        
        drlee  = 0. ; tau_src(1:levs+1) =  0.  
        return
      ENDIF
!=========================================================================	  
!  find orographic asymmetry and convexity :'oa/clx' for clipped SG-hill
!             nwd  1   2   3   4   5   6   7   8
!              wd  w   s  sw  nw   e   n  ne  se
! make sure that SM_00 and KD_05 oro-characteristics can match each other
! OD-KDO5    = Gamma=a/b [0:2]   ; hsg = 2.*hprime
! OC-KD05      mount sharpness sigma^4   "height to half-width"[0:1]
! alph-SM00    fraction of h2d contributed to hprime           [0:1]
!
! OA-KDO5      OA > dwstream  OA=0 sym   OA < 0 upstram [-1. 0. 1]
! delt-SM00    dw/up asymmetry -1 < delta < 1
! Gamma-LM97   anisotropy of the orography   g2 =(dh/dx)^2/(dh/dy)^2
!..
!A parametrization of low-level wave breaking which includes a dependence on
!the degree of 2-dimensionality of SG; it is active over a finite range of Fr
!=========================================================================
      wdir   = atan2(uhm,vhm) + pi
      idir   = mod( int(fdir*wdir),mdir) + 1

      nwd    = nwdir(idir)
      ind4 =   mod(nwd-1,4) + 1
      if (ind4 < 1 ) ind4 = 1
      if (ind4 > 4 ) ind4 = 4

      oa      = ( 1-2*int( (nwd-1)/4 )) * oa4(ind4)
      clx     = clx4(ind4)
      cl4p(1) = clx4(2)
      cl4p(2) = clx4(1)
      cl4p(3) = clx4(4)
      cl4p(4) = clx4(3)
      clxp   = cl4p(ind4) 

      odir = clxp/max(clx, 1.e-5)    ! WRF-based definition for "odir" 

      odir = min(odmax, odir)  
      odir = max(odmin, odir)  


      if  (strver == 'smc_2000' .or. strver == 'vay_2018') then
!=========================================================================
!
!  thrree-piece def-n for    tautot(Fr): 0-Fr_lin - Fr_lee -Fr_mtb
!                                     taulin/tauogw taulee    taumtb
!  here    tau_src(levs+1): approximate wave flux from surface to LLWB
!  Following attempts of Scinocca +McFarlane, 2000 & Eckermann etal.(2010)
!=========================================================================  
! 
!  if (mkz2 < 0)... mkzi = sqrt(-mkz2) trapped wave regime don't a case in UGWP-V1
!  wave flux ~ rho_src*kx_src/mkz_src*wind_rms
!              bn2, uhm, vhm, bn2hm, rhohm
!
!      IF (mkz2.ge. mkz2min .and. krefj > 2 ) THEN
!
! wave regimes
! 
        mkz  = sqrt(mkz2)
        frlwb   = fcrit_sm          ! should  be higher than LOGW to get  zblk < zlwb
        frlin   = fcrit_sm
        frlin2  = 1.5*fcrit_sm
        frlin3  = 3.0*fcrit_sm

        hcrit   = fcrit_sm*umag/bnv      
        hogw    = min(hamp, hcrit)
        hdsw    = min(hamp, frlwb*umag/bnv)      ! no trapped-wave solution

        coef_h2 = kxridge * rhohm * bnv * umag 

        taulin  = coef_h2 * hamp*hamp 
        tauhcr  = coef_h2 * hcrit*hcrit 

        IF (fr < frlin ) then 
           tauogw =  taulin  
           taulee =  0.0
           taumtb =  0.0
        else if (fr .ge. frlin ) then 
            tauogw = tauhcr 
            taulin = coef_h2 * hamp*hamp 
            taumtb = tau_izb                 ! integrated form MTB 
!
! SM-2000 approach for taulee, shall we put limits on BetaMax_max  ~ 20 or  Betaf   ??
! 
          frnorm   = fr/fcrit_sm      ! frnorm below [1.0 to 3.0]
          BetaMax  = 1.0 + 2.0*OC     ! alpha of SM00 or OC-mountain sharphess KD05 OC=[10, 0]

          if ( fr <= frlin2 ) then 
            Betaf= 2.*BetaMax*(frNorm-1.0)
            taulee = (1. + Betaf )*taulin - tauhcr 
          else if ( (fr > frlin2).and.(fr <= frlin3))then
            Betaf=-1.+ 1./frnorm/frnorm + & 
             (BetaMax + 0.555556)*(2.0 - 0.666*frnorm)* (2.0 - 0.666*frnorm)         
            taulee = (1. + Betaf )*taulin - tauhcr    
!==============
! Eck-2010 WRF-alternatve through Dp_surf = P'*grad(h(x,y))
! 1 < Fr < 2.5  tauwave = taulee+tauogw = tau_dp*(fr)**(-0.9)
!    Fr  > 2.5  tauwave = tau_dp*(2.5)**(-0.9)
! to apply it need tabulated Dp(fr, Dlin) Dp=function(Dlin, U, N, h)
!
!============== 
          else
            taulee = 0.0
            hdsw   = 0.0
          endif
        ENDIF

        tautot = tauogw + taulee + taumtb*0.

        IF  (taulee > 0.0 )  THEN

          hdsw = sqrt(tautot/coef_h2)   ! averaged value for hdsw - mixture of lee+ogw with mkz/kxridge
!
! compute  vertical profile "drlee" with the low-level wave breaking & "locally" trapped waves
!  make "empirical" height above elvp that may represent DSW-wave breaking & trapping
!  here  we will assign tau_sso(z) profile between:   zblk(zsurf) - zlwb - ztop_sso = ns*sridge
!
          call mflow_tauz(levs, up, vp, tp, qp, dp,   zpm, zpi, &
                 pmid, pint, rho, ui, vi, ti,  bn2i, bvi,  rhoi)

          kph  = max(izb, 2)    !  kph marks  the low-level of wave solutions
          klwb = kph            !  klwb above blocking marks  wave-breaking
          kotr = levs+1         !  kotr marks mkz2(z) <= 0., reflection level

          if (do_dtrans) pdtrans = log(tautot/tauogw)/(zpi(krefj) - zpi(kph))

          udir =  max(ui(krefj)*xn +vi(krefj)*yn, velmin)
          hogwi = hogw*hogw* rhohm/rhoi(krefj) * umag/udir * bnv/bvi(krefj)
          umd(krefj) = udir

          udir =  max(ui(kph)*xn +vi(kph)*yn, velmin)
          hdswi = hdsw*hdsw* rhohm/rhoi(kph) * umag/udir * bnv/bvi(kph)
          umd(kph) = udir
                                 ! what we can put between k =[kph:krefj]
          phmkz(:) = 0.0                                !               
          phmkz(kph-1) = fr     ! initial Phase of the low-level wave
!
!  now transfer tau_layer => tau_level assuming   tau_layer = tau_level
!   kx*rho_layer*bn_layer*u_layer* HL*HL = kx*rho_top*bn_top*u_top * HT*HT        
!   apply it for both hdsw  & hogw  with linear saturation-solver for Cx =0
!
   loop_lwb_otr:  do k=kph+1, krefj                   ! levs

            umd(k) = max(ui(k)*xn +vi(k)*yn, velmin)
            umd2   =(coro- umd(k))*(coro- umd(k))
            umd2   = max(umd2, dw2min) -c2f2 
 
           
            if (umd2 <= 0.0) then 
!
! critical layer
!
              klwb = k
              kotr = k
              exit   loop_lwb_otr
            endif
           
            mkz2 = bn2i(k)/umd2 - kx2
          
            if ( mkz2 >= mkz2min ) then
!
! find klwb having some "kinematic" phase "break-down" crireria SM00 or LM97
! at finest  vertical resolution  we can meet "abrupt" mkz
!  mkzmax = 6.28/(2*dz), mkzmin = 6.28/ztrop=18km
! to regularize SG-solution  mkz = max(mkzmax, min(mkz,in, mkz)) 
!
              mkz  = sqrt(mkz2)
              hdswz = hdswi* rhoi(k-1)/rhoi(k) * umd(k-1)/umd(k) * bvi(k-1)/bvi(k)
              udswz = hdswz *bn2i(k)
!===========================================================================================
!linsat wave ampl.:     mkz*sqrt(hdswz) <= 1.0 or     udswz <= linsat2*umd2
!
!                tautot = tausat = rhoi(k) *udswz_sat * kxridge/mkz
! by k = krefj   tautot = tauogw(krefj)
!===========================================================================================
              if  (do_klwb_phase)  then                    
                 phmkz(k) = phmkz(k-1) + mkz*(zpm(k)-zpm(k-1))
                 if( ( phmkz(k) .ge. zbr_pi).and.(klwb == kph)) then
                   klwb = min(k, krefj)     
                   exit loop_lwb_otr
                 endif 
              endif 
            else                !  mkz2 < mkz2min
               kotr = k     ! trapped/reflected  waves /
               exit   loop_lwb_otr
            endif
          enddo     loop_lwb_otr
!
! define tau_src(1:zblk:klwb) = sum(tau_oro+tau_dsw+tau_ogw) and define drlee
!        tau_trapped ???
! 
          if (do_klwb_phase) then
            do k=kph, kotr-1

             if (klwb > kph .and. k < klwb)  then
               drlee(k)   =  (tautot -tauogw)/(zpi(kph) - zpi(klwb))   !  negative Ax*rho
               tau_src(k) = tautot + (zpi(k) - zpi(klwb))*drlee(k)
               drlee(k)   = drlee(k)/rho(k)         
             else if ( k >=  klwb .and. k < kotr) then
               tau_src(k) = tauogw
               drlee(k)   = 0.0
             endif
           enddo  
           kdswj = klwb   ! assign to the "low-level" wave breaking
         endif
!
! simplest exponential transmittance d(tau)/dz = - pdtrans *tau(z)
! more  complicated is dissipative saturation pdtrans =/= constant
!
         if (do_dtrans) then
           do k=kph, krefj
             tau_src(k)= tautot*exp(-pdtrans*(zpi(k)-zpi(kph)))
             drlee(k)  = -tau_src(k)/rho(k) * pdtrans 
           enddo
         endif


       ENDIF      !taulee > 0.0 


     endif       !strver 
!

!=========================================================================       
     if  (strver == 'gfs_2018' .or. strver == 'kd_2005') then
!=========================================================================
!
! orowaves: OGW+DSW/Lee
!
       efact    = (oa + 2.0) ** (ceofrc*fr)
       efact    = min( max(efact,efmin), efmax )
       coefm    = (1. + clx) ** (oa+1.)

       kxlinv    = min (kxw,  coefm * cleff)                ! does not exceed  42km ~4*dx
       kxlinv    = coefm * cleff        
       tem       = fr  * fr * oc
       gfobnv   = gmax * tem / ((tem + cg)*bnv)             ! g/n0
!=========================================================================	
! source fluxes:  taulin,  taufrb
!=========================================================================
       tautot  = kxlinv * rhohm * umag * umag *umag* gfobnv  * efact

       coef_h2 = kxlinv *rhohm * bnv*umag
       taulin  = coef_h2 *hamp*hamp
       hcrit   = fcrit_gfs*umag/bnv
       tauhcr  = coef_h2 *hcrit*hcrit

       IF (fr <= fcrit_gfs) then
         tauogw  = taulin   
         tautot  = taulin   
         taulee  =  0.
         drlee(:) = 0.
       ELSE                      !fr > fcrit_gfs
         tauogw  = tauhcr
         taulee  = max(tautot - tauogw, 0.0)
         if (taulee > 0.0 ) hdsw = sqrt(taulee/coef_h2) 
!  approximate drlee(k) between  [izb, klwb]
! find klwb and decrease taulee(izb) => taulee(klwb) = 0.
! above izb  tau
         if (mkz2 > mkz2min.and. krefj > 2 .and. taulee > 0.0) then

           mkz = sqrt(mkz2)
           call mflow_tauz(levs, up, vp, tp, qp, dp,   zpm, zpi, &
                pmid, pint, rho, ui, vi, ti,  bn2i, bvi,  rhoi)
          
           kph  = max(izb, 2)
           phmkz(:) = 0.0
           klwb = max(izb, 1)
           kotr = levs+1
           phmkz(kph-1) = fr  ! initial Phase of the Lee-OGW

           loop_lwb_gfs18:  do k=kph, levs

             umd(k) = max(ui(k)*xn +vi(k)*yn, velmin)
             umd2   =(coro- umd(k))*(coro- umd(k))
             umd2   = max(umd2, velmin*velmin)
             mkz2   = bn2i(k)/umd2 - kx2
             if ( mkz2 > mkz2min ) then
               mkz      = sqrt(mkz2)
               frlocal  = max(hdsw*bvi(k)/umd(k), frlwb) 
               phmkz(k) = phmkz(k-1) + mkz*(zpm(k)-zpm(k-1))
               if( ( phmkz(k) >= zbr_pi ) .and. (frlocal > frlin)) klwb = k
             else 
               kotr = k
               exit loop_lwb_gfs18
             endif
           enddo   loop_lwb_gfs18
!
!
           do k=kph, kotr-1

             if (klwb > kph .and. k < klwb) then
               drlee(k)   = -(tautot -tauogw)/(zpi(kph) - zpi(klwb))
               tau_src(k) = tautot + (zpi(k) - zpi(klwb))*drlee(k)
               drlee(k)   = drlee(k)/rho(k)
             else if ( k >= klwb .and. k < kotr) then
               tau_src(k) = tauogw
               drlee(k)   = 0.0
             endif
           enddo  
           kdswj = klwb   ! assign to the "low-level" wave breaking
         endif              !  mkz2 > mkz2min.and. krefj > 2 .and. taulee > 0.0
       ENDIF               !fr > fcrit_gfs


     ENDIF     !strbase='gfs2017' .or. strbase='kd_2005'		 
 

! output : taulin,  taufrb, taulee, xn, yn, umag, kxw/kxridge	
!       print *,  krefj, levs,  tauogw,  tautot , ' ugwp_taub_oro ' 
!
      end subroutine ugwp_taub_oro
!
!--------------------------------------
! 
!     call ugwp_oro_lsatdis( krefj, levs,  tauogw(j),  tautot(j), tau_src,      &
!           con_pi, con_g, kxw, fcor(j), c2f2(j), up, vp, tp, qp, dp, zpm, zpi, &
!           pmid1, pint1, xn, yn, umag, drtau, kdis_oro)
 
      subroutine ugwp_oro_lsatdis( krefj, levs,  tauogw,  tautot,  tau_src,     &
           pi, grav, kxw, fcor,  kxridge, up, vp, tp, qp, dp, zpm, zpi,        &
           pmid, pint, xn, yn, umag, drtau, kdis)

      use ugwp_common_v1,       only : dw2min, velmin
      use cires_ugwp_module_v1, only : frcrit, ricrit, linsat, hps, rhp1, rhp2
      use cires_ugwp_module_v1, only : kvg, ktg, krad, kion
      use ugwp_oro_init_v1,     only : coro , fcrit_sm , fcrit_sm2
      implicit none 
!      	   
      integer, intent(in)        ::   krefj,   levs
      real , intent(in)          ::   tauogw,  tautot, kxw
      real , intent(in)          ::   fcor

      real , dimension(levs+1)   ::  tau_src

      real, intent(in)           ::  pi, grav

      real, dimension(levs) , intent(in)     ::  up, vp, tp, qp, dp, zpm
      real, dimension(levs+1), intent(in)    ::  zpi, pmid, pint
      real , intent(in)                      ::  xn, yn, umag
      real , intent(in)                      ::  kxridge


      real, dimension(levs), intent(out)      ::  drtau, kdis
!
! locals
!
      real  :: bnv2min, pi2, rgrav
      real                       ::  uref, udir, uf2, ufd, uf2p
      real, dimension(levs+1)    ::  tauz
      real, dimension(levs)      ::  rho
      real, dimension(levs+1)    ::  ui, vi, ti, bn2i, bvi, rhoi
         
      integer                    :: i, j, k, kcrit, kref
      real                       :: kx2,  kx2w, kxs
      real                       :: mkzm, mkz, dkz, mkz2, ch, kzw3
      real                       :: wfdM, wfdT, wfiM, wfiT
      real                       :: fdis, mkzi, keff_m, keff_t
      real                       :: betadis, betam, betat, cdfm, cdft
      real                       :: fsat, hsat, hsat2, kds , c2f2

      pi2 = 2.0*pi
      bnv2min = (pi2/1800.)*(pi2/1800.)
      rgrav = 1.0/grav

      drtau(1:levs) =  0.0
      kdis (1:levs) =  0.0

         ch = coro

         kx2w = kxw*kxw
         kx2 = kxridge*kxridge
         if( kx2 < kx2w ) kx2 = kx2w
         kxs = sqrt(kx2)
         c2f2 = fcor*fcor/kx2
!
! non-hydrostatic LinSatDis for Ch = 0 (with set of horizontal wavenumber kxw)
!  
!      print *,  krefj, levs,  tauogw,  tautot , ' orolsatdis '
      call mflow_tauz(levs, up, vp, tp, qp, dp, zpm, zpi, &
                 pmid, pint, rho, ui, vi, ti, bn2i, bvi, rhoi) 
!===============================================================================
! for stationary  oro-GWs only "single"-azimuth cd = 0 -(-Udir) = Udir > 0
!  rotational/non-hyrostatic effects  are important only for high-res runs
!  Udir = 0, Udir < 0  are not  
!   future"revisions"  shear effects for d mkz /dt  = -kxw*dU/dz
!                      horizontal wavelength spectra mkz2 = l2 -kxw(n)*kxw(n) 
!                      stochastic "tauogw'-setup+ sigma_tau ; 
!                      3D-wave effects 1+ (k/l)^2      and NS vs EW orowaves
!                      target is to get "multiple"-saturation levels for OGWs
!===============================================================================
         tauz(1:krefj)  =  tauogw          ! constant flux for OGW-packet or single mode
                                           ! sign of tauz > 0...and its attenuate with Z
         k = krefj
         uref = ui(k)*xn +vi(k)*yn - ch    ! stationary waves
         uf2  =  uref*uref - c2f2
         if (uf2 > 0) then
           mkz2 = bn2i(k)/uf2 -kx2
           if (mkz2.gt.0) then
             mkzm = sqrt(mkz2)
           else
             return                        ! wave reflection mkz2 <=0.
           endif
         else
           return                          ! wave absorption uf2 <= 0.
         endif
!
! upward solver for single "mode" with tauz(levs+1) =0. at the top 
! 
         kds   =   0.1*  kvg(krefj)   ! eddy wave diffusion from the previous layer
         kcrit =   levs
         do k= krefj+1, levs  
!
! 2D-wave propagation along reference-wind direction
!           udir = 0 critical wind for coro =0
!           cdop = -uref .... upwind waves travel against MF
!
           udir =  ui(k)*xn +vi(k)*yn
           uf2  =  udir*udir - c2f2
        
    
           if (uf2 < dw2min .or. udir <= 0.0) then
             kcrit =K
             tauz(kcrit:levs) = 0. 
             exit                    ! vert-level loop
           endif
!
! wave-based solution
!
           mkz2 = bn2i(k)/uf2 -kx2
           if (mkz2 > 0) then 
             mkzm = sqrt(mkz2) 
!
! do dissipative flux vs saturation: kvg, ktg, krad, kion 
!
             kzw3 = mkzm*mkz2
!
             keff_m =  kvg(k)*mkz2 + kion(k)
!            keff_t = kturb(k)*iPr_turb + kmol(k)*iPr_mol
             keff_t = ktg(k)*mkz2  + krad(k)
!       
!
             uf2p    = uf2 + 2.0*c2f2
             betadis = uf2/uf2p
             betaM   = 1.0 / (1.0+betadis)       ! if c2f2 = 0. betaM = betaT  =0.5 ekw = epw
             betaT   = 1.0- BetaM

!
!imaginary frequencies of momentum and heat with "kds at (k-1) level"
!
             wfiM    = kds*mkz2  + keff_m
             wfiT    = kds*mkz2  + keff_t
!
             cdfm    = sqrt(uf2)*kxs
             cdft    = abs(udir)*kxs
             wfdM    = wfiM/cdfm  *BetaM
             wfdT    = wfiT/Cdft  *BetaT
             mkzi    = 2.0*mkzm*(wfdM+wfdT)

             fdis    = tauz(k-1)*exp(-mkzi*(zpi(k)-zpi(k-1)) )
             tauz(k) = fdis
             hsat2   = fcrit_sm2 * uf2 *bn2i(k)
             fsat    = rhoi(k)* hsat2 * sqrt(uf2) * bvi(k)
             if (fdis > fsat) then 
               tauz(k)  = min(fsat, tauz(k-1))
!=================================================================
! two definitions for eddy mixing of MF:
!  a) wave damping-Lindzen :  Ked ~ kx/(2H)*(u-c)^4/N^3
!  b) heat-based turbulence: 4/3 Richardson Ked ~eps^1/3 *Lt^4/3 
!=================================================================
               kds     = rhp2*kxs*uf2*uf2/bn2i(k)/bvi(k)
               kdis(k) = kds 
             endif         
           else
             tauz(k:levs) = 0.     ! wave is reflected above
             kds  = 0.
           endif  
         enddo

         do k=krefj+1, kcrit
           drtau(k) = rgrav*(tauz(k+1)-tauz(k))/dp(k)
         enddo         
!
!
      end  subroutine ugwp_oro_lsatdis
!
!
     subroutine ugwp_tofd(im, levs, con_cp, sigflt, elvmax, zpbl, u, v, zmid, &
                          utofd, vtofd, epstofd, krf_tofd)
       use machine ,      only : kind_phys 
       use ugwp_oro_init_v1, only : n_tofd, const_tofd, ze_tofd, a12_tofd, ztop_tofd
!       
     implicit none
!     
     integer ::  im, levs
     real(kind_phys) :: con_cp
     real(kind_phys), dimension(im, levs)  ::  u, v, zmid
     real(kind_phys), dimension(im)        ::  sigflt, elvmax, zpbl
     real(kind_phys), dimension(im, levs)  ::  utofd, vtofd, epstofd, krf_tofd
!
! locals
!
     integer :: i, k
     real    :: rcpd2
     real    :: sgh = 30.
     real    :: sgh2, ekin, zdec, rzdec, umag, zmet, zarg, zexp, krf
!     
     utofd =0.0 ; vtofd = 0.0 ;  epstofd =0.0 ; krf_tofd =0.0    
     rcpd2 = 0.5/con_cp
!      
     
     do i=1, im 
       
       zdec = max(n_tofd*sigflt(i), zpbl(i))
       zdec = min(ze_tofd, zdec)
       rzdec = 1.0/zdec
       sgh2  = max(sigflt(i)*sigflt(i), sgh*sgh)
       
      do k=1, levs
         zmet = zmid(i,k)   
      if (zmet > ztop_tofd) cycle
       ekin = u(i,k)*u(i,k) + v(i,k)*v(i,k)
       umag = sqrt(ekin)
       zarg = zmet*rzdec
       zexp = exp(-zarg*sqrt(zarg))
       krf         = const_tofd* a12_tofd *sgh2* zmet ** (-1.2) *zexp
       utofd(i,k)  = -krf*u(i,k)
       vtofd(i,k)  = -krf*v(i,k)
       epstofd(i,k)= rcpd2*krf*ekin  ! more accurate heat/mom form using "implicit tend-solver"
                                     ! to update momentum and temp-re
       krf_tofd(i,k) = krf
      enddo
     enddo
!                
     end subroutine ugwp_tofd  
!
!     
     subroutine ugwp_tofd1d(levs, con_cp, sigflt, elvmax, zsurf, zpbl,  u, v, &
                            zmid, utofd, vtofd, epstofd, krf_tofd)
       use machine ,      only : kind_phys 
       use ugwp_oro_init_v1, only : n_tofd, const_tofd, ze_tofd, a12_tofd, ztop_tofd
!       
     implicit none
      integer ::  levs
      real(kind_phys) :: con_cp
      real(kind_phys), dimension(levs)  ::   u, v, zmid
      real(kind_phys)                   ::  sigflt, elvmax, zpbl, zsurf
      real(kind_phys), dimension(levs)  ::  utofd, vtofd, epstofd, krf_tofd
!
! locals
!
     integer :: i, k
     real    :: rcpd2
     real    :: sghmax = 5.
     real    :: sgh2, ekin, zdec, rzdec, umag, zmet, zarg, ztexp, krf
!     
     utofd =0.0 ; vtofd = 0.0 ;  epstofd =0.0 ; krf_tofd =0.0    
     rcpd2 = 0.5/con_cp
!         
       zdec = max(n_tofd*sigflt, zpbl)          ! ntimes*sgh_turb or Zpbl
       zdec = min(ze_tofd, zdec)                ! cannot exceed 18 km
       rzdec = 1.0/zdec
       sgh2 = max(sigflt*sigflt, sghmax*sghmax) ! 25 meters dz-of the first layer
       
      do k=1, levs  
         zmet = zmid(k)-zsurf   
      if (zmet > ztop_tofd) cycle
         ekin = u(k)*u(k) + v(k)*v(k)
         umag = sqrt(ekin)
         zarg = zmet*rzdec
         ztexp = exp(-zarg*sqrt(zarg))
         krf   = const_tofd* a12_tofd *sgh2* zmet ** (-1.2) *ztexp

         utofd(k)    = -krf*u(k)
         vtofd(k)    = -krf*v(k)
         epstofd(k)  = rcpd2*krf*ekin  ! more accurate heat/mom form using "implicit tend-solver"
                                       ! to update momentum and temp-re; epstofd(k) can be skipped
         krf_tofd(k) = krf
     enddo
!                
     end subroutine ugwp_tofd1d    


end module cires_vert_orodis_v1

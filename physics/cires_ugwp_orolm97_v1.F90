module cires_ugwp_orolm97_v1


contains



     subroutine gwdps_oro_v1(im,  km,    imx, do_tofd,           &
         pdvdt, pdudt, pdtdt, pkdis, u1,v1,t1,q1,kpbl,           &  
         prsi,del,prsl,prslk, zmeti, zmet, dtp, kdt, hprime,     &
         oc, oa4, clx4, theta, sigmad, gammad, elvmaxd,          &
         grav, con_omega, rd, cpd, rv, pi, arad, fv, sgh30,      &
         dusfc, dvsfc,  xlatd, sinlat, coslat, sparea,           &
         cdmbgwd, me, master, rdxzb,                             &
         zmtb, zogw, tau_mtb, tau_ogw, tau_tofd,                 &
         dudt_mtb, dudt_ogw, dudt_tms)
!----------------------------------------
! ugwp_v1: gwdps_oro_v1 following recent updates of Lott & Miller 1997
!          eventually will be replaced with more "advanced"LLWB
!          and multi-wave solver that produce competitive FV3GFS-skills  
! 
!   computation of kref for ogw + coorde diagnostics
!   all constants/parameters inside cires_ugwp_initialize.f90 
!----------------------------------------    

      use machine ,      only : kind_phys
      use ugwp_common_v1,  only : dw2min, velmin
 
      use ugwp_oro_init_v1, only : rimin,  ric,     efmin,     efmax   ,  &
                                hpmax,  hpmin,   sigfaci => sigfac  ,  &
                                dpmin,  minwnd,  hminmt,    hncrit  ,  &
                                rlolev, gmax,    veleps,    factop  ,  &
                                frc,    ce,      ceofrc,    frmax, cg, &
                                fdir,   mdir,    nwdir,                &
                                cdmb,   cleff,   fcrit_gfs, fcrit_mtb, &
                                n_tofd, ze_tofd, ztop_tofd

      use cires_ugwp_module_v1, only : kxw,  max_kdis, max_axyz

      use cires_orowam2017, only : oro_wam_2017

      use cires_vert_orodis_v1, only : ugwp_tofd1d
      
      
!      use sso_coorde,        only : pgwd, pgwd4
!----------------------------------------
      implicit none
      real(kind=kind_phys), parameter :: pgwd=1, pgwd4= pgwd    
      real(kind=kind_phys), parameter :: sigfac = 3, sigfacs = 0.5      
      character(len=8)                :: strsolver='pss-1986'  ! current operational solver or  'wam-2017'
      real(kind=kind_phys)            :: gammin = 0.00999999
      real(kind=kind_phys), parameter :: nhilmax = 25.
      real(kind=kind_phys), parameter :: sso_min = 3000.
      logical, parameter              :: do_adjoro = .false.
!----------------------------------------      
      
      integer, intent(in) :: im, km, imx, kdt
      integer, intent(in) :: me, master
      logical, intent(in) :: do_tofd
      
      
 
      integer, intent(in)              :: kpbl(im)    ! index for the pbl top layer!
      real(kind=kind_phys), intent(in) :: dtp         !  time step
      real(kind=kind_phys), intent(in) :: cdmbgwd(2)
      
      real(kind=kind_phys), intent(in) :: hprime(im), oc(im), oa4(im,4),   &
                                  clx4(im,4), theta(im), sigmad(im),        &
                                  gammad(im), elvmaxd(im)

      real(kind=kind_phys), intent(in) :: grav, con_omega, rd, cpd, rv,     &
                                          pi, arad, fv
      real(kind=kind_phys), intent(in) :: sgh30(im)       
      real(kind=kind_phys), intent(in), dimension(im,km) ::   &
                                  u1,  v1,   t1, q1,del, prsl, prslk, zmet
     
      real(kind=kind_phys), intent(in),dimension(im,km+1):: prsi, zmeti
      real(kind=kind_phys), intent(in) :: xlatd(im),sinlat(im), coslat(im)
      real(kind=kind_phys), intent(in) :: sparea(im)

!    
!output -phys-tend
      real(kind=kind_phys),dimension(im,km),intent(out) ::   &
                           pdvdt,    pdudt,    pkdis, pdtdt
! output - diag-coorde
      real(kind=kind_phys),dimension(im,km),intent(out) ::   &
                          dudt_mtb, dudt_ogw, dudt_tms
!                     
      real(kind=kind_phys),dimension(im) :: rdxzb,  zmtb,  zogw ,    &
                    tau_ogw, tau_mtb, tau_tofd, dusfc,   dvsfc

!
!---------------------------------------------------------------------
! # of permissible sub-grid orography hills for "any" resolution  < 25
!    correction for "elliptical" hills based on shilmin-area =sgrid/25 
!     4.*gamma*b_ell*b_ell  >=  shilmin
!     give us limits on [b_ell & gamma *b_ell] > 5 km =sso_min
!     gamma_min = 1/4*shilmin/sso_min/sso_min
!23.01.2019:  cdmb = 4.*192/768_c192=1 x 0.5
!     192: cdmbgwd        = 0.5, 2.5
!     cleff = 2.5*0.5e-5 * sqrt(192./768.) => lh_eff = 1004. km
!      6*dx = 240 km 8*dx = 320. ~ 3-5 more effective OGW-lin
!---------------------------------------------------------------------
! 
! locals SSO
!
      real(kind=kind_phys) :: vsigma(im),  vgamma(im)

      real(kind=kind_phys)            :: ztoph,zlowh,ph_blk, dz_blk
      real(kind=kind_phys)            :: shilmin, sgrmax, sgrmin
      real(kind=kind_phys)            :: belpmin, dsmin,  dsmax
!     real(kind=kind_phys)            :: arhills(im)              ! not used why do we need?
      real(kind=kind_phys)            :: xlingfs
         
! 
! locals        mean flow  ...etc
!
      real(kind=kind_phys), dimension(im,km) :: ri_n, bnv2, ro
      real(kind=kind_phys), dimension(im,km) :: vtk, vtj, velco
!mtb     
      real(kind=kind_phys), dimension(im)    :: oa,  clx , sigma, gamma,  &
                                                elvmax, wk
      real(kind=kind_phys), dimension(im)    :: pe, ek, up
      
      real(kind=kind_phys), dimension(im,km) :: db, ang, uds

      real(kind=kind_phys) :: zlen, dbtmp, r, phiang, dbim, zr
      real(kind=kind_phys) :: eng0, eng1, cosang2, sinang2
      real(kind=kind_phys) :: bgam, cgam, gam2, rnom, rdem       
!
! tofd
!     some constants now in "use ugwp_oro_init" +   "use ugwp_common"
!
!==================
      real(kind=kind_phys)   :: unew, vnew,  zpbl,  sigflt, zsurf
      real(kind=kind_phys), dimension(km)    :: utofd1, vtofd1
      real(kind=kind_phys), dimension(km)    :: epstofd1, krf_tofd1
      real(kind=kind_phys), dimension(km)    :: up1, vp1, zpm
      
      real(kind=kind_phys),dimension(im, km) :: axtms, aytms
! 
! ogw
!
      logical icrilv(im)
!
      real(kind=kind_phys), dimension(im) :: xn, yn, ubar, vbar, ulow, &
                     roll,  bnv2bar, scor, dtfac, xlinv, delks, delks1
!
      real(kind=kind_phys) :: taup(im,km+1), taud(im,km)
      real(kind=kind_phys) :: taub(im), taulin(im), heff, hsat, hdis

      integer, dimension(im) :: kref, idxzb, ipt, kreflm, iwklm, iwk, izlow
   
!      
!check what we need
!
      real(kind=kind_phys) ::   bnv,  fr, ri_gw, brvf
      real(kind=kind_phys) ::   tem,   tem1,  tem2, temc, temv
      real(kind=kind_phys) ::   ti,     rdz,   dw2,   shr2, bvf2
      real(kind=kind_phys) ::   rdelks, efact, coefm, gfobnv
      real(kind=kind_phys) ::   scork,  rscor, hd,    fro,  sira
      real(kind=kind_phys) ::   dtaux,  dtauy, zmetp, zmetk
      
      real(kind=kind_phys) ::   grav2, rcpdt, windik, wdir
      real(kind=kind_phys) ::   sigmin, dxres,sigres,hdxres, cdmb4, mtbridge

      real(kind=kind_phys) ::   kxridge, inv_b2eff, zw1, zw2
      real(kind=kind_phys) ::   belps, aelps, nhills, selps

      real(kind=kind_phys) ::   rgrav, rcpd, rcpd2, rad_to_deg, deg_to_rad
      real(kind=kind_phys) ::   pi2, rdi, gor, grcp, gocp, gr2, bnv2min

!      
! various integers
!     
      integer ::   kmm1, kmm2, lcap, lcapp1
      integer ::   npt,   kbps, kbpsp1,kbpsm1
      integer ::   kmps,  idir, nwd,  klcap, kp1, kmpbl, kmll
      integer ::   k_mtb, k_zlow, ktrial, klevm1
      integer ::   i, j, k
! 
! initialize gamma and sigma
      gamma(:) = gammad(:)
      sigma(:) = sigmad(:)
!
      rcpdt = 1.0 / (cpd*dtp)
      grav2 = grav + grav
!
      rgrav = 1.0/grav
      rcpd = 1.0/cpd
      rcpd2 = 0.5/cpd
      rad_to_deg=180.0/pi
      deg_to_rad=pi/180.0
      pi2 = 2.*pi
      rdi = 1.0/rd
      gor = grav/rd
      grcp = grav*rcpd
      gocp = grcp
      gr2  = grav*gor
      bnv2min = (pi2/1800.)*(pi2/1800.)
!       
! mtb-blocking  sigma_min and dxres => cires_initialize
!  
      sgrmax = maxval(sparea) ; sgrmin = minval(sparea)
      dsmax  = sqrt(sgrmax)   ; dsmin  = sqrt(sgrmin)

      dxres   = pi2*arad/float(imx)
      hdxres  = 0.5*dxres
!     shilmin = sgrmin/nhilmax            ! not used - moorthi

!     gammin = min(sso_min/dsmax, 1.)     ! moorthi - with this results are not reproducible
      gammin = min(sso_min/dxres, 1.)     ! moorthi

!     sigmin = 2.*hpmin/dsmax      !dxres ! moorthi - this will not reproduce
      sigmin = 2.*hpmin/dxres      !dxres

!     if (kdt == 1) then
!       print *, sgrmax, sgrmin , ' min-max sparea '
!       print *, 'sigmin-hpmin-dsmax', sigmin, hpmin, dsmax
!       print *, 'dxres/dsmax ', dxres, dsmax
!       print *, ' shilmin gammin ', shilmin, gammin
!     endif

      kxridge = float(imx)/arad * cdmbgwd(2)

      if (me == master .and. kdt == 1) then
        print *, ' gwdps_v0 kxridge ', kxridge
        print *, ' gwdps_v0 scale2 ', cdmbgwd(2)
        print *, ' gwdps_v0 imx ', imx
        print *, ' gwdps_v0 gam_min ', gammin
        print *, ' gwdps_v0 sso_min ', sso_min
      endif

      do i=1,im
        idxzb(i)    = 0
        zmtb(i)     = 0.0
        zogw(i)     = 0.0
        rdxzb(i)    = 0.0      
        tau_ogw(i)  = 0.0
        tau_mtb(i)  = 0.0 
        dusfc(i)    = 0.0
        dvsfc(i)    = 0.0
        tau_tofd(i) = 0.0
!
        ipt(i) = 0
! 
      enddo

      do k=1,km
        do i=1,im
          pdvdt(i,k)    = 0.0
          pdudt(i,k)    = 0.0
          pdtdt(i,k)    = 0.0
          pkdis(i,k)    = 0.0
          dudt_mtb(i,k) = 0.0
          dudt_ogw(i,k) = 0.0
          dudt_tms(i,k) = 0.0
        enddo
      enddo
 
! ----  for lm and gwd calculation points
!cires_ugwp_initialize.F90:      real, parameter :: hpmax=2400.0, hpmin=25.0  
!cires_ugwp_initialize.F90:       real,      parameter :: hminmt=50.     ! min mtn height (*j*) 
!----  for lm and gwd calculation points  


      npt = 0
      
      do i = 1,im
        if ( elvmaxd(i) >= hminmt .and. hprime(i)  >= hpmin ) then          
          npt      = npt + 1
          ipt(npt) = i
	endif  
      enddo
      
      if (npt == 0) then
!         print *,  'oro-npt = 0 elvmax ', maxval(elvmaxd), hminmt
!         print *,  'oro-npt = 0 hprime ', maxval(hprime), hpmin	     	    
        return      ! no gwd/mb calculation done
      endif
!========================================================

!      
      if (do_adjoro ) then 
      	  
        do i = 1,im	  
!         arhills(i) = 1.0
!
          sigres = max(sigmin, sigma(i))
!         if (sigma(i) < sigmin) sigma(i)=  sigmin
          dxres = sqrt(sparea(i))
          if (2.*hprime(i)/sigres > dxres) sigres=2.*hprime(i)/dxres
          aelps = min(2.*hprime(i)/sigres, 0.5*dxres)
          if (gamma(i) > 0.0 ) belps = min(aelps/gamma(i),.5*dxres)
!
! small-scale "turbulent" oro-scales < sso_min
!
          if( aelps < sso_min ) then
 
! a, b > sso_min upscale ellipse  a/b > 0.1 a>sso_min & h/b=>new_sigm
!
            aelps = sso_min 
            if (belps < sso_min ) then
              gamma(i) = 1.0
              belps = aelps*gamma(i)
            else
              gamma(i) = min(aelps/belps, 1.0)
            endif
            
	    sigma(i) = 2.*hprime(i)/aelps
            gamma(i) = min(aelps/belps, 1.0)
	    
          endif
 
          selps      = belps*belps*gamma(i)*4.    ! ellipse area of the el-c hill 
          nhills     = min(nhilmax, sparea(i)/selps)
!         arhills(i) = max(nhills, 1.0)

!333   format( ' nhil: ', i6, 4(2x, f9.3), 2(2x, e9.3))	    
!      if (kdt==1 )
!     & write(6,333) nint(nhills)+1,xlatd(i), hprime(i),aelps*1.e-3,
!     &   belps*1.e-3, sigma(i),gamma(i)

        
      enddo
      endif        !(do_adjoro )



      do i=1,npt
        iwklm(i)  = 2
        idxzb(i)  = 0 
        kreflm(i) = 0
      enddo
 
      do k=1,km
        do i=1,im
          db(i,k)  = 0.0
          ang(i,k) = 0.0
          uds(i,k) = 0.0 
        enddo
      enddo

      kmm1 = km - 1 ;  kmm2   = km - 2 ; kmll   = kmm1
      lcap = km     ;  lcapp1 = lcap + 1 
      
      cdmb4 = 0.25*cdmb 
       
      do i = 1, npt
        j = ipt(i)
        elvmax(j) = min (elvmaxd(j)*0. + sigfac * hprime(j), hncrit)
        izlow(i)  = 1          ! surface-level
      enddo
!
      do k = 1, kmm1
        do i = 1, npt
          j = ipt(i)
          ztoph   = sigfac * hprime(j)
          zlowh   = sigfacs* hprime(j) 
          zmetp   =  zmet(j,k+1) 
          zmetk   =  zmet(j,k)  
!         if (( elvmax(j) <= zmetp) .and. (elvmax(j).ge.zmetk) )
!     &      iwklm(i)  =  max(iwklm(i), k+1 ) 
          if (( ztoph <= zmetp) .and. (ztoph >= zmetk) ) iwklm(i)  =  max(iwklm(i), k+1 )
          if (zlowh <= zmetp .and. zlowh >= zmetk)       izlow(i)  =  max(izlow(i),k)
    
        enddo
      enddo
!
      do k = 1,km
        do i =1,npt
          j         = ipt(i)
          vtj(i,k)  = t1(j,k)  * (1.+fv*q1(j,k))
          vtk(i,k)  = vtj(i,k) / prslk(j,k)
          ro(i,k)   = rdi * prsl(j,k) / vtj(i,k)       ! density mid-levels
          taup(i,k) = 0.0
        enddo
      enddo
!
! check ri_n or ri_mf computation
!
      do k = 1,kmm1
        do i =1,npt
          j         = ipt(i)
          rdz       = 1.   / (zmet(j,k+1) - zmet(j,k))
          tem1      = u1(j,k) - u1(j,k+1)
          tem2      = v1(j,k) - v1(j,k+1)
          dw2       = tem1*tem1 + tem2*tem2
          shr2      = max(dw2,dw2min) * rdz * rdz
!          ti        = 2.0 / (t1(j,k)+t1(j,k+1))
!          bvf2      = grav*(gocp+rdz*(vtj(i,k+1)-vtj(i,k)))* ti
!          ri_n(i,k) = max(bvf2/shr2,rimin)   ! richardson number
!
          bvf2 = grav2 * rdz * (vtk(i,k+1)-vtk(i,k))/ (vtk(i,k+1)+vtk(i,k))
     
          bnv2(i,k+1) = max( bvf2, bnv2min )
          ri_n(i,k+1) = bnv2(i,k)/shr2        ! richardson number consistent with bnv2	
!
! add here computation for ktur and ogw-dissipation fro ve-gfs
!	    
        enddo
      enddo
      k = 1
      do i = 1, npt
        bnv2(i,k) = bnv2(i,k+1)
      enddo
!		
! level iwklm => zmet(j,k) < sigfac * hprime(j) < zmet(j,k+1) 
!
      do i = 1, npt
        j   = ipt(i)
        k_zlow = izlow(i)
        if (k_zlow == iwklm(i)) k_zlow = 1
        delks(i)   = 1.0 / (prsi(j,k_zlow) - prsi(j,iwklm(i)))
!       delks1(i)  = 1.0 /(prsl(j,k_zlow) - prsl(j,iwklm(i)))
        ubar (i)   = 0.0
        vbar (i)   = 0.0
        roll (i)   = 0.0
        pe   (i)   = 0.0
        ek   (i)   = 0.0
        bnv2bar(i) = 0.0   
      enddo
!
      do i = 1, npt
        k_zlow = izlow(i)
        if (k_zlow == iwklm(i)) k_zlow = 1
        do k = k_zlow, iwklm(i)-1                        ! kreflm(i)= iwklm(i)-1 
          j       = ipt(i)                               ! laye-aver rho, u, v
          rdelks  = del(j,k) * delks(i)
          ubar(i) = ubar(i)  + rdelks * u1(j,k)          ! trial mean u below 
          vbar(i) = vbar(i)  + rdelks * v1(j,k)          ! trial mean v below 
          roll(i) = roll(i)  + rdelks * ro(i,k)          ! trial mean ro below 
!   
          bnv2bar(i) = bnv2bar(i) + .5*(bnv2(i,k)+bnv2(i,k+1))* rdelks
        enddo
      enddo
!
      do i = 1, npt
        j = ipt(i)
!
! integrate from ztoph = sigfac*hprime  down to zblk if exists
! find ph_blk, dz_blk like in LM-97 and ifs
!	
        ph_blk =0.  
        do k = iwklm(i), 1, -1
          phiang   =  atan2(v1(j,k),u1(j,k))*rad_to_deg
          ang(i,k) = ( theta(j) - phiang )
          if ( ang(i,k) >  90. ) ang(i,k) = ang(i,k) - 180.
          if ( ang(i,k) < -90. ) ang(i,k) = ang(i,k) + 180.
          ang(i,k) = ang(i,k) * deg_to_rad
          uds(i,k) = max(sqrt(u1(j,k)*u1(j,k) + v1(j,k)*v1(j,k)), velmin)
!
          if (idxzb(i) == 0 ) then
            dz_blk = zmeti(j,k+1) - zmeti(j,k)
            pe(i)  =  pe(i) + bnv2(i,k) *( elvmax(j) - zmet(j,k) ) * dz_blk

            up(i)  =  max(uds(i,k) * cos(ang(i,k)), velmin)  
            ek(i)  = 0.5 *  up(i) * up(i) 

            ph_blk = ph_blk + dz_blk*sqrt(bnv2(i,k))/up(i)

! --- dividing stream lime  is found when pe =exceeds ek. oper-l gfs
!           if ( pe(i) >=  ek(i) ) then
            if ( ph_blk >=  fcrit_gfs ) then
               idxzb(i) = k
               zmtb (j) = zmet(j, k)
               rdxzb(j) = real(k, kind=kind_phys)
            endif

          endif
        enddo
!
! alternative expression: zmtb = max(heff*(1. -fcrit_gfs/fr), 0)
! fcrit_gfs/fr	 
!
        goto 788

        bnv     = sqrt( bnv2bar(i) )
        heff    = 2.*min(hprime(j),hpmax)
        zw2     = ubar(i)*ubar(i)+vbar(i)*vbar(i)
        ulow(i) = sqrt(max(zw2,dw2min))
        fr      = heff*bnv/ulow(i)
        zw1     = max(heff*(1. -fcrit_gfs/fr), 0.0)
        zw2     = zmet(j,2)
        if (fr > fcrit_gfs .and. zw1 > zw2 ) then 
          do k=2, kmm1
            zmetp =  zmet(j,k+1) 
            zmetk   =  zmet(j,k)   
            if (zw1 <= zmetp .and. zw1 >= zmetk)  exit
          enddo
            idxzb(i) = k
            zmtb (j) = zmet(j, k)
        else
           zmtb (j) = 0.
           idxzb(i) = 0
        endif
	
788     continue
!
! --- the drag for mtn blocked flow
!
        if ( idxzb(i) > 0 ) then
	
! (4.16)-ifs	  
          gam2 = gamma(j)*gamma(j)
          bgam = 1.0 - 0.18*gamma(j) - 0.04*gam2
          cgam =       0.48*gamma(j) + 0.30*gam2
	  
          do k = idxzb(i)-1, 1, -1
            zlen = sqrt( (zmtb(j)-zmet(j,k) )/(zmet(j,k ) + hprime(j)) )
            tem     = cos(ang(i,k))
            cosang2 = tem * tem
            sinang2 = 1.0 - cosang2 
!	      
!  cos =1 sin =0 =>   1/r= gam     zr = 2.-gam 
!  cos =0 sin =1 =>   1/r= 1/gam   zr = 2.- 1/gam
!
            rdem = cosang2      +  gam2 * sinang2
            rnom = cosang2*gam2 +         sinang2
!	       
! metoffice dec 2010
! correction of h. wells & a. zadra for the
! aspect ratio  of the hill seen by mean flow
! (1/r , r-inverse below: 2-r)

            rdem = max(rdem, 1.e-6)       
            r    = sqrt(rnom/rdem)
            zr   =  max( 2. - r, 0. )

            sigres = max(sigmin, sigma(j))
            if (hprime(j)/sigres > dxres) sigres = hprime(j)/dxres
            mtbridge = zr * sigres*zlen / hprime(j)
! (4.15)-ifs 	   
!           dbtmp = cdmb4 * mtbridge *                               &
!     &           max(cos(ang(i,k)), gamma(j)*sin(ang(i,k)))
! (4.16)-ifs
            dbtmp  = cdmb4*mtbridge*(bgam* cosang2 +cgam* sinang2)
            db(i,k)= dbtmp * uds(i,k)
          enddo
!                  
        endif
      enddo
!.............................
!.............................
! end  mtn blocking section
!.............................
!.............................
!
!--- orographic gravity wave drag section
!     
!  scale cleff between im=384*2 and 192*2 for t126/t170 and t62
!  inside "cires_ugwp_initialize.f90" now
!
      kmpbl  = km / 2 
      iwk(1:npt) = 2
!
! meto/UK-scheme: 
! k_mtb = max(k_zmtb, k_n*hprime/2] to reduce diurnal variations taub_ogw 
!     
      do k=3,kmpbl
        do i=1,npt
          j   = ipt(i)
          tem = (prsi(j,1) - prsi(j,k))
          if (tem < dpmin) iwk(i) = k           ! dpmin=50 mb
  
!===============================================================	  
! lev=111      t=311.749     hkm=0.430522     ps-p(iwk)=52.8958 
!           below "hprime" - source of ogws  and below zblk !!!
!           27           2  kpbl ~ 1-2 km   < hprime
!===============================================================	  
        enddo
      enddo
!
! iwk - adhoc gfs-parameter to select ogw-launch level between
!      level ~0.4-0.5 km from surface or/and  pbl-top
! in ugwp-v1: options to modify as  htop ~ (2-3)*hprime > zmtb
! in ugwp-v0 we ensured that : zogw > zmtb
!

      kbps  = 1
      kmps  = km
      k_mtb = 1
      do i=1,npt
        j         = ipt(i)
        k_mtb     = max(1, idxzb(i))

        kref(i)   = max(iwk(i),  kpbl(j)+1 )            ! reference level pbl or smt-else ????
        kref(i)   = max(kref(i), iwklm(i) )             ! iwklm => sigfac*hprime

        if (kref(i) <= k_mtb)  kref(i) = k_mtb + 1      ! layer above zmtb
        kbps      = max(kbps,  kref(i))
        kmps      = min(kmps,  kref(i))
!
        delks(i)  = 1.0 / (prsi(j,k_mtb) - prsi(j,kref(i)))
        ubar (i)  = 0.0
        vbar (i)  = 0.0
        roll (i)  = 0.0
        bnv2bar(i)= 0.0
      enddo
!
      kbpsp1 = kbps + 1
      kbpsm1 = kbps - 1
      k_mtb  = 1
!
      do i = 1,npt
        k_mtb = max(1, idxzb(i))
        do k = k_mtb,kbps            !kbps = max(kref) ;kmps= min(kref)
          if (k < kref(i)) then
            j          = ipt(i)
            rdelks     = del(j,k) * delks(i)
            ubar(i)    = ubar(i)  + rdelks * u1(j,k)   ! mean u below kref
            vbar(i)    = vbar(i)  + rdelks * v1(j,k)   ! mean v below kref
            roll(i)    = roll(i)  + rdelks * ro(i,k)   ! mean ro below kref
            bnv2bar(i) = bnv2bar(i) + .5*(bnv2(i,k)+bnv2(i,k+1))* rdelks
          endif
        enddo
      enddo
!
! orographic asymmetry parameter (oa), and (clx)
      do i = 1,npt
        j      = ipt(i)
        wdir   = atan2(ubar(i),vbar(i)) + pi
        idir   = mod(nint(fdir*wdir),mdir) + 1
        nwd    = nwdir(idir)
        oa(i)  = (1-2*int( (nwd-1)/4 )) * oa4(j,mod(nwd-1,4)+1)
        clx(i) = clx4(j,mod(nwd-1,4)+1)
      enddo
!
      do i = 1,npt
       dtfac(i)  = 1.0
       icrilv(i) = .false.                      ! initialize critical level control vector
       ulow(i) = max(sqrt(ubar(i)*ubar(i)+vbar(i)*vbar(i)),velmin)
       xn(i)  = ubar(i) / ulow(i)
       yn(i)  = vbar(i) / ulow(i)
      enddo
!
      do  k = 1, kmm1
        do  i = 1,npt
          j            = ipt(i)
          velco(i,k)   = 0.5 * ((u1(j,k)+u1(j,k+1))*xn(i)+  (v1(j,k)+v1(j,k+1))*yn(i))
          
        enddo
      enddo
!
!------------------
! v0: incorporates latest modifications for kxridge and heff/hsat
!             and taulin for fr <=fcrit_gfs 
!             and concept of "clipped" hill if zmtb > 0. to make
! the integrated "tau_sso = tau_ogw +tau_mtb" close to reanalysis data
!      it is still used the "single-orowave"-approach along ulow-upwind
!
!      in contrast to the 2-orthogonal wave (2otw) schemes of ifs/meto/e-canada  
! 2otw scheme requires "aver angle" and wind projections on 2 axes of ellipse a-b
!     with 2-stresses:  taub_a & taub_b as of Phillips  (1984)
!------------------
      taub(:)  = 0. ; taulin(:)= 0.
      do i = 1,npt
        j    = ipt(i)
        bnv  = sqrt( bnv2bar(i) )
        heff = min(hprime(j),hpmax)

        if( zmtb(j) > 0.) heff = max(sigfac*heff-zmtb(j), 0.)/sigfac  
        if (heff <= 0) cycle

        hsat = fcrit_gfs*ulow(i)/bnv
        heff = min(heff, hsat)

        fr   = min(bnv  * heff /ulow(i), frmax)
!
        efact    = (oa(i) + 2.) ** (ceofrc*fr)
        efact    = min( max(efact,efmin), efmax )
!
        coefm    = (1. + clx(i)) ** (oa(i)+1.)
!
        xlinv(i) = coefm * cleff           ! effective kxw for lin-wave
        xlingfs  = coefm * cleff 
!
        tem      = fr    * fr * oc(j)
        gfobnv   = gmax  * tem / ((tem + cg)*bnv) 
!
!new specification of xlinv(i) & taulin(i)

        sigres = max(sigmin, sigma(j))
        if (heff/sigres > hdxres) sigres = heff/hdxres
        inv_b2eff =  0.5*sigres/heff
        kxridge   =  1.0 / sqrt(sparea(j))
        xlinv(i)  = xlingfs    !or max(kxridge, inv_b2eff)  ! 6.28/lx ..0.5*sigma(j)/heff = 1./lridge  
        taulin(i) = 0.5*roll(i)*xlinv(i)*bnv*ulow(i)*heff*heff*pgwd4

        if ( fr > fcrit_gfs ) then
          taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)  &
                     * ulow(i)  * gfobnv  * efact          ! nonlinear flux tau0...xlinv(i)
!
        else
!         taub(i)  = taulin(i)                             !  linear flux for fr <= fcrit_gfs	 
          taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)   &
                     * ulow(i)  * gfobnv  * efact                    
!   
        endif
!		 
!
        k       = max(1, kref(i)-1)
        tem     = max(velco(i,k)*velco(i,k), dw2min)
        scor(i) = bnv2(i,k) / tem           ! scorer parameter below kref level
!
! diagnostics for zogw > zmtb
!
        zogw(j) = zmeti(j, kref(i) )
      enddo
!                                                                       
!----set up bottom values of stress
!
      do k = 1, kbps
        do i = 1,npt
          if (k <= kref(i)) taup(i,k) = taub(i)
        enddo
      enddo
      
      if (strsolver == 'pss-1986') then     
       
!======================================================
!   v0-gfs orogw-solver of palmer et al 1986 -"pss-1986"
!   in v1-orogw  linsatdis of                 "wam-2017"
!     with llwb-mechanism for
!     rotational/non-hydrostat ogws important for 
!     highres-fv3gfs with dx < 10 km
!======================================================
      
        do k = kmps, kmm1                   ! vertical level loop from min(kref)
          kp1 = k + 1
          do i = 1, npt
!
            if (k >= kref(i)) then
              icrilv(i) = icrilv(i) .or. ( ri_n(i,k) < ric).or. (velco(i,k) <= 0. )
            endif
          enddo
!
          do i = 1,npt
            if (k >= kref(i))   then
              if (.not.icrilv(i) .and. taup(i,k) > 0.0 ) then
                temv = 1.0 / max(velco(i,k), velmin)
!
                if (oa(i) > 0. .and. kp1 < kref(i)) then
!				
                  scork   = bnv2(i,k) * temv * temv
                  rscor   = min(1.0, scork / scor(i))
                  scor(i) = scork
                else 
                  rscor   = 1.
                endif
!
                brvf = sqrt(bnv2(i,k))        ! brent-vaisala frequency interface
!               tem1 = xlinv(i)*(ro(i,kp1)+ro(i,k))*brvf*velco(i,k)*0.5

                tem1 = xlinv(i)*(ro(i,kp1)+ro(i,k))*brvf*0.5        &
                              * max(velco(i,k), velmin)
                hd   = sqrt(taup(i,k) / tem1)
                fro  = brvf * hd * temv
!
!    rim is the  "wave"-richardson number by palmer,shutts, swinbank 1986
!

                tem2   = sqrt(ri_n(i,k))
                tem    = 1. + tem2 * fro
                ri_gw  = ri_n(i,k) * (1.0-fro) / (tem * tem)
!
!    check stability to employ the 'dynamical saturation hypothesis'
!    of palmer,shutts, swinbank 1986
!                                       
             if (ri_gw <= ric .and.(oa(i) <= 0. .or.  kp1 >= kref(i) )) then
                   temc = 2.0 + 1.0 / tem2
                   hd   = velco(i,k) * (2.*sqrt(temc)-temc) / brvf
                   taup(i,kp1) = tem1 * hd * hd
              else 
                   taup(i,kp1) = taup(i,k) * rscor
              endif
!	      
                taup(i,kp1) = min(taup(i,kp1), taup(i,k))
              endif
            endif
          enddo
        enddo
!     
!  zero momentum deposition at the top model layer
!      
        taup(1:npt,km+1) = taup(1:npt,km)      
!
!     calculate wave acc-n: - (grav)*d(tau)/d(p) = taud
!
        do k = 1,km
          do i = 1,npt
            taud(i,k) = grav*(taup(i,k+1) - taup(i,k))/del(ipt(i),k)
          enddo
        enddo
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!------if the gravity wave drag would force a critical line in the
!------layers below sigma=rlolev during the next deltim timestep,
!------then only apply drag until that critical line is reached.
! empirical implementation of the llwb-mechanism: lower level wave breaking
! by limiting "ax = dtfac*ax" due to possible llwb around kref and 500 mb
! critical line [v - ax*dtp = 0.] is smt like "llwb" for stationary ogws
!2019:  this option limits sensitivity of taux/tauy to  increase/decrease of taub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        do k = 1,kmm1
          do i = 1,npt
           if (k  >= kref(i) .and. prsi(ipt(i),k) >= rlolev) then

               if(taud(i,k) /= 0.) then
                 tem = dtp * taud(i,k)                          ! tem = du/dt-oro*dt => U/dU vs 1
                 dtfac(i) = min(dtfac(i),abs(velco(i,k)/tem))   ! reduce Ax= Ax*(1, or U/dU <=1)
!	         dtfac(i) = 1.0
               endif
            endif
          enddo
        enddo
!
!--------------------------- orogw-solver of gfs pss-1986
!  
      else 
!
!-----------Unified orogw-solver of wam2017
!
!       sigres = max(sigmin, sigma(j))
!	if (heff/sigres.gt.dxres) sigres=heff/dxres
!         inv_b2eff =  0.5*sigres/heff 	
!       xlinv(i) = max(kxridge, inv_b2eff)           ! 0.5*sigma(j)/heff = 1./lridge      

        dtfac(:) =  1.0  
	     
        call oro_wam_2017(im, km, npt, ipt, kref, kdt, me, master,             &
           dtp, dxres, taub, u1, v1, t1, xn, yn, bnv2, ro, prsi,prsl,          &
           grav, con_omega, rd,                                                &
           del, sigma, hprime, gamma, theta, sinlat, xlatd, taup, taud, pkdis)
     
      endif            !  oro_wam_2017 - linsatdis-solver of wam-2017
!      
!---- above orogw-solver of wam2017
! 
! tofd as in beljaars-2004
!
! ---------------------------    
      if( do_tofd ) then
        axtms(:,:) = 0.0 ; aytms(:,:) = 0.0 
        if ( kdt == 1 .and. me == 0) then
          print *, 'vay do_tofd  from surface to ', ztop_tofd 
        endif
        do i = 1,npt          
          j = ipt(i)
          zpbl  = zmet( j, kpbl(j) )
 
          sigflt = min(sgh30(j), 0.3*hprime(j)) ! cannot exceed 30% of ls-sso
 
          zsurf = zmeti(j,1)
          do k=1,km
            zpm(k) = zmet(j,k)
            up1(k) = u1(j,k)
            vp1(k) = v1(j,k)
          enddo
 
          call ugwp_tofd1d(km, cpd, sigflt, elvmaxd(j), zsurf, zpbl,  &
               up1, vp1, zpm,  utofd1, vtofd1, epstofd1, krf_tofd1)
     
          do k=1,km
            axtms(j,k) = utofd1(k)
            aytms(j,k) = vtofd1(k)
!	 
! add tofd to gw-tendencies
!	 
            pdvdt(j,k)  = pdvdt(j,k) + aytms(j,k)
            pdudt(j,k)  = pdudt(j,k) + axtms(j,k)
          enddo
!2018-diag
          tau_tofd(j) = sum( utofd1(1:km)* del(j,1:km))
        enddo
      endif             ! do_tofd 

!--------------------------------------------
! combine oro-drag effects MB +TOFD + OGWs
!--------------------------------------------  
! +  diag-3d

      dudt_tms = axtms 
      tau_ogw  = 0.
      tau_mtb  = 0.
 
      do k = 1,km
        do i = 1,npt
          j    = ipt(i)
!
          eng0 = 0.5*(u1(j,k)*u1(j,k)+v1(j,k)*v1(j,k))
!
          if ( k < idxzb(i) .and. idxzb(i) /= 0 ) then
!
! if blocking layers -- no ogws
!	  
            dbim       = db(i,k) / (1.+db(i,k)*dtp)
            pdvdt(j,k) = - dbim * v1(j,k) +pdvdt(j,k)
            pdudt(j,k) = - dbim * u1(j,k) +pdudt(j,k)
            eng1       = eng0*(1.0-dbim*dtp)*(1.-dbim*dtp)
            
            dusfc(j)   = dusfc(j) - dbim * u1(j,k) * del(j,k)
            dvsfc(j)   = dvsfc(j) - dbim * v1(j,k) * del(j,k)
!2018-diag 
            dudt_mtb(j,k) = -dbim * u1(j,k)
            tau_mtb(j)    = tau_mtb(j) + dudt_mtb(j,k)* del(j,k)

          else
!
! ogw-s above blocking height
!	  
            taud(i,k)  = taud(i,k) * dtfac(i)
            dtaux      = taud(i,k) * xn(i) * pgwd
            dtauy      = taud(i,k) * yn(i) * pgwd
    
            pdvdt(j,k)   = dtauy  +pdvdt(j,k)
            pdudt(j,k)   = dtaux  +pdudt(j,k)
    
            unew   = u1(j,k) + dtaux*dtp       !   pdudt(j,k)*dtp
            vnew   = v1(j,k) + dtauy*dtp       !   pdvdt(j,k)*dtp
            eng1   = 0.5*(unew*unew + vnew*vnew)
!
            dusfc(j) = dusfc(j)  + dtaux * del(j,k)
            dvsfc(j) = dvsfc(j)  + dtauy * del(j,k)
!2018-diag
            dudt_ogw(j,k) = dtaux
            tau_ogw(j)    = tau_ogw(j) +dtaux*del(j,k)
          endif
!	  
! local energy deposition sso-heat
!	
          pdtdt(j,k) = max(eng0-eng1,0.)*rcpdt  
        enddo
      enddo
! dusfc w/o tofd  sign as in the era-i, merra  and cfsr
      do i = 1,npt
         j           = ipt(i)
         dusfc(j)    = -rgrav * dusfc(j)
         dvsfc(j)    = -rgrav * dvsfc(j)
         tau_mtb(j)  = -rgrav * tau_mtb(j)
         tau_ogw(j)  = -rgrav * tau_ogw(j)
         tau_tofd(j) = -rgrav * tau_tofd(j)
       enddo

       return


!============ debug ------------------------------------------------
       if (kdt <= 2 .and. me == 0) then
        print *, 'vgw-oro done gwdps_v0 in ugwp-v0 step-proc ', kdt, me
!
        print *, maxval(pdudt)*86400.,  minval(pdudt)*86400, 'vgw_axoro'
        print *, maxval(pdvdt)*86400.,  minval(pdvdt)*86400, 'vgw_ayoro'
!       print *, maxval(kdis),  minval(kdis),  'vgw_kdispro m2/sec'
        print *, maxval(pdtdt)*86400.,  minval(pdtdt)*86400,'vgw_epsoro'
        print *, maxval(zmtb), ' z_mtb ',  maxval(tau_mtb), ' tau_mtb '
        print *, maxval(zogw), ' z_ogw ',  maxval(tau_ogw), ' tau_ogw '
!       print *, maxval(tau_tofd),  ' tau_tofd '
!       print *, maxval(axtms)*86400.,  minval(axtms)*86400, 'vgw_axtms'
!       print *,maxval(dudt_mtb)*86400.,minval(dudt_mtb)*86400,'vgw_axmtb'
        if (maxval(abs(pdudt))*86400. > 100.) then

          print *, maxval(u1),  minval(u1),  ' u1 gwdps-v0 '
          print *, maxval(v1),  minval(v1),  ' v1 gwdps-v0 '
          print *, maxval(t1),  minval(t1),  ' t1 gwdps-v0 '
          print *, maxval(q1),  minval(q1),  ' q1 gwdps-v0 '
          print *, maxval(del), minval(del), ' del gwdps-v0 '
          print *, maxval(zmet),minval(zmet), 'zmet'
          print *, maxval(zmeti),minval(zmeti), 'zmeti'
          print *, maxval(prsi), minval(prsi), ' prsi '
          print *, maxval(prsl), minval(prsl), ' prsl '
          print *, maxval(ro), minval(ro), ' ro-dens '
          print *, maxval(bnv2(1:npt,:)), minval(bnv2(1:npt,:)),' bnv2 '
          print *, maxval(kpbl), minval(kpbl), ' kpbl '
          print *, maxval(sgh30), maxval(hprime), maxval(elvmax),'oro-d'
          print *
          do i =1, npt
            j= ipt(i)
            print *,zogw(j)/hprime(j), zmtb(j)/hprime(j), &
                   zmet(j,1)*1.e-3, nint(hprime(j)/sigma(j))
!
!....................................................................
!
!   zogw/hp=5.9   zblk/hp=10.7    zm=11.1m   ridge/2=2,489m/9,000m
!   from 5 to 20 km , we need to count for "ridges" > dx/4 ~ 15 km
!   we must exclude blocking by small ridges
!   vay-kref < iblk           zogw-lev 15        block-level:  39
!
! velmin => 1.0, 0.01, 0.1 etc.....unification of wind limiters
! max(sqrt(u1(j,k)*u1(j,k) + v1(j,k)*v1(j,k)), minwnd)
! max(dw2,dw2min) * rdz * rdz
! ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0)
! tem      = max(velco(i,k)*velco(i,k), 0.1)
!              temv = 1.0 / max(velco(i,k), 0.01)
!     &                   * max(velco(i,k),0.01) 
!....................................................................     
          enddo
          print *
          stop
        endif
       endif
      
!cires_ugwp_solv2_v1.f90 
      return
      end subroutine gwdps_oro_v1 


end module cires_ugwp_orolm97_v1

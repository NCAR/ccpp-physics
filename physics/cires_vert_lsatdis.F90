   subroutine ugwp_lsatdis_naz(levs, ksrc, nw, naz, kxw, taub_spect, ch, xaz, yaz,  &                                 
              fcor, c2f2, dp, zmid, zint, pmid, pint, rho, ui, vi, ti,              &
              kvg, ktg, krad, kion, bn2i, bvi, rhoi, ax, ay, eps, ked, tau1)
!
!         call ugwp_lsatdis_naz(levs, ksrc, nw, naz, kxw, taub_spect, ch, xaz, yaz, &                                
!              fcor(j), c2f2(j), dp, zmid, zint, pmid, pint, rho, ui, vi, ti,       &
!	      kvg, ktg, krad, kion, bn2i, bvi, rhoi, ax1, ay1, eps1, ked1)
   use    ugwp_common, only  : rcpd, grav, rgrav    
   implicit none 
!   
   integer ::  levs, nw, naz, ksrc
   real    ::  kxw
   real, dimension(nw)     ::   taub_spect, ch
   real, dimension(naz)    ::   xaz, yaz
   real, dimension(levs+1) ::  ui, vi, ti, bn2i,  bvi, rhoi, zint, pint
   real, dimension(levs  ) ::  dp, rho, pmid, zmid
   real                    ::  fcor, c2f2
   real, dimension(levs+1) ::  kvg, ktg,  kion, krad, kmol 
   
! output/locals        
   real, dimension(levs  )       ::  ax, ay, eps
   real, dimension(levs+1)       ::  ked , tau1

   real, dimension(levs+1  )     ::  uaz
   real, dimension(levs,  naz  ) ::  epsd
   real, dimension(levs+1, naz ) ::  atau, kedd
   real, dimension(levs+1      ) ::  taux, tauy  
   real, dimension(levs  )       ::  dzirho , dzpi 
   real                          ::  usrc
!  
   integer     :: iaz,  k 
!     
   atau=0.0 ;  epsd=0.0 ; kedd=0.0  
   
   do k=1,levs
    dzpi(k)   = -(pint(k+1)-pint(k))/rho(k)*rgrav 
    dzirho(k) = 1./rho(k)/dzpi(k)          !   grav/abs(dp(k)) still hydrostatic "UGWP"
   enddo
   
   LOOP_IAZ:  do iaz =1, naz
                usrc = ui(ksrc)*xaz(iaz) +vi(ksrc)*yaz(iaz)  
                do k=1,levs+1
                  uaz(k) =ui(k)*xaz(iaz) +vi(k)*yaz(iaz) -usrc  
                enddo
!
!      if (nw .le. 4) call stochastic ..ugwp_lsatdis_az1  only 4-waves ch_ngw1, fuw_ngw1, eff_ngw1=1
!
! multi-wave scheme
!
      if (nw .gt. 4) then 
       call  ugwp_lsatdis_az1(levs, ksrc, nw, kxw, ch, taub_spect,               &
          fcor, c2f2, zmid, zint, rho, uaz, ti, bn2i, bvi, rhoi, dzirho, dzpi,   &
          kvg, ktg, krad, kion, kmol, epsd(:, iaz), kedd(:,iaz), atau(:, iaz) )
       
      endif
!                               
      ENDDO LOOP_IAZ                    ! Azimuth of GW propagation directions
!
!     sum over azimuth and project aTau(z, iza) =>(taux and tauy)
!         for scalars              for "wave-drag vector"
!
         eps =0. ; ked =0.
          do k=ksrc, levs
            eps(k) = sum(epsd(k,:))*rcpd
         enddo

         do k=ksrc, levs+1
            taux(k) = sum( atau(k,:)*xaz(:))
            tauy(k) = sum( atau(k,:)*yaz(:))
            ked(k)  = sum(kedd(k,:))
         enddo  

         tau1(ksrc:levs) = taux(ksrc:levs) 
         tau1(1:ksrc-1) = tau1(ksrc)
!
! end  solver: gw_azimuth_solver_LS81    
!      sign Ax in rho*dU/dt  = -d(rho*tau)/dz
!                             [(k) - (k+1)]
         ax =0. ; ay = 0.
         do k=ksrc, levs
           ax(k) = dzirho(k)*(taux(k)-taux(k+1))
           ay(k) = dzirho(k)*(tauy(k)-tauy(k+1))
         enddo  
       call ugwp_limit_1d(ax, ay, eps, ked, levs) 
      return 

!
      print *
      print *, ' Ax: ', maxval(Ax(ksrc:levs))*86400.,  minval(Ax(ksrc:levs))*86400.
      print *, ' Ay: ', maxval(Ay(ksrc:levs))*86400.,  minval(Ay(ksrc:levs))*86400.
      print *, 'Eps: ', maxval(Eps(ksrc:levs))*86400.,  minval(Eps(ksrc:levs))*86400.
      print *, 'Ked: ', maxval(Ked(ksrc:levs))*1.,  minval(Ked(ksrc:levs))*1.
!      print *, 'Atau  ', maxval(atau(ksrc:levs, 1:Naz))*1.e3, minval(atau(ksrc:levs, 1:Naz))*1.e3
!      print *, 'taux_gw: ', maxval(taux( ksrc:levs))*1.e3,  minval(taux( ksrc:levs))*1.e3
      print *     
!-----------------------------------------------------------------------
! Here we can apply "ad-hoc" or/and "stability-based" limiters on
! (axy_gw, ked_gw and eps_gw) and check vert-inegrated conservation laws:
! energy and momentum and after that => final update gw-phys tendencies
!-----------------------------------------------------------------------
         
    end subroutine ugwp_lsatdis_naz
!
    subroutine ugwp_lsatdis_az1(levs, ksrc, nw, kxw, ch, taub_sp,               &
          fcor, c2f2, zm, zi, rho, um, tm, bn2, bn, rhoi,                       &
          dzirho, dzpi, kvg, ktg, krad, kion, kmol, eps, ked, tau )
  
!       call  ugwp_lsatdis_az1(levs, ksrc, nw, kxw, ch, taub_spect,               &
!          fcor, c2f2, zmid, zint, rho, uaz, ti, bn2i, bvi, rhoi, dzirho, dzpi,   &
!          kvg, ktg, krad, kion, kmol, epsd(:, iaz), kedd(:,iaz), atau(:, iaz) )

      use cires_ugwp_module, only : F_coriol, F_nonhyd, F_kds, linsat, linsat2
      use cires_ugwp_module, only : iPr_ktgw, iPr_spgw, iPr_turb, iPr_mol
      use cires_ugwp_module, only : rhp4, rhp2, rhp1, khp, cd_ulim
!                                                              
      implicit NONE
!
      integer, intent(in) :: nw                         ! number of GW modes in given direction
      integer, intent(in) :: levs                       ! vertical layers
      integer, intent(in) :: ksrc                       ! level of GW-launch layer
      
      real , intent(in)  :: kxw                         ! horizontal wavelength
      real , intent(in)  :: ch(nw)                      ! horizontal phase velocities
      real , intent(in)  :: taub_sp(nw)                 ! spectral distribution of the mom-flux
!
      real,  intent(in)   :: fcor, c2f2                 ! Corilois factors

      real , intent(in)   :: um(levs+1)
      real , intent(in)   :: tm(levs+1)
!in
      real, intent(in), dimension(levs)     :: rho,  zm
      real, intent(in), dimension(levs+1)   :: rhoi, zi
      real, intent(in), dimension(levs+1)   :: bn2, bn
      real, intent(in), dimension(levs)     :: dzpi, dzirho
      real, intent(in), dimension(levs+1)   :: kvg, ktg, krad, kion, kmol
!========================================================================
!out
      real, dimension(levs+1) :: tau, ked
      real, dimension(levs)   :: eps 

!=========================================================================
!local
      real                       :: Fd1, Fd2
      real, dimension(levs)      :: a_mkz
      real, dimension(levs+1,nw) :: sp_tau, sp_ked, sp_kth
      real, dimension(levs,nw)   :: sp_eps

      real, dimension(levs,nw)   :: sp_mkz, sp_etot
      real, dimension(levs,nw)   :: sp_ek, sp_ep


      real, dimension(levs)      :: swg_ep, swg_ek, swg_et, swg_kz

      real, dimension(nw)       :: rtaus         ! spectral distribution at ksrc
      real                      :: sum_rtaus     ! total flux in iaz-azimuth
      real :: Chnorm, Cx, Cs, Cxs, Cx2sat
      real :: Fdis, Fdisat
      real :: Cdf2, Cdf1                         ! (Cd*cd-f*f) and sqrt   
!
! two-level => upward integration for wave-filtering (dissip + breaking)
!
      real :: taus, tauk, tau_lin
      real :: etws, etwk, etw_lin
      real :: epss, epsk
      real :: kds, kdk
      real :: kzw, kzw2, kzw3, kzi, kzs
      real :: wfd, wfi                           !
!
! for GW dissipation on the rotational sphere
!
      real :: Betadis                            ! Ep/Ek ratio
      real :: BetaM, BetaT                       ! 0.5 or 1./1+b and 1-1/(1+b)
      real :: wfdM, wfdT, wfiM, wfiT, wdop2

      real :: dzi, keff, keff_m, keff_t, keffs

      real ::  sf2k2, cf2
      real ::  Lzkm, Lzsat                 
      
      integer :: i, k, igw
      integer :: ksat1, ksat2
      
      real    :: zsat1, zsat2
      real    :: kx2_nh
      
      real :: rab1, rab2, rab3, rab4, cd_ulim2
      
      integer :: Ind_out(nw,    levs+1)
      
!
      logical, parameter ::  dbg_print = .false.    
!
!===================================================================             
!     Nullify arrays
!     tau, eps, ked
!====================================================================
     
      tau = 0.0
      eps = 0.0
      ked = 0.0
      Ind_out(1:nw,:) = 0
!
! GW-spectral arrays ..... sp_etot ....sp_tau 
!
      sp_tau  = 0.
      sp_eps  = 0.
      sp_ked  = 0.
      sp_mkz  = -99.
      sp_etot = 0.
      sp_ek   = 0.
      sp_ep   = 0.       
      sp_kth  = 0.
!
      swg_et = 0.
      swg_ep = 0.
      swg_ek = 0.
      swg_kz = 0.
      cd_ulim2 = cd_ulim*cd_ulim
      cf2      = F_coriol*c2f2
      kx2_nh   = F_nonhyd*kxw*kxw

      if (dbg_print) then
        write(6,*) linsat , ' eff-linsat & kx ', kxw
        write(6,*) maxval(ch), minval(ch), ' ch '
        write(6,*)
        write(6,*) maxval(rhoi), minval(rhoi), 'rhoi '
        write(6,*) zi(ksrc) , ' zi(ksrc) '
        write(6,*) cd_ulim,     ' crit-level  cd_ulim '
        write(6,*) F_coriol,    ' F_coriol'
        write(6,*) F_nonhyd,    ' F_nonhyd '
        write(6,*) maxval(Bn), minval(BN), ' BN-BV '
        write(6,*) Um(ksrc),    ' Um-ksrc ', cd_ulim2 ,  'cd_ulim2 ', c2f2, ' c2f2 '
        !pause
      endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Loop_GW:  over GW-spectra 
!               of individual non-interactive modes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    
  Loop_GW:  do i=1, nw
!
            Kds = 0.0
!
! src-level
!
         Cx       = ch(i) - Um(ksrc)
         Cdf2     = Cx*Cx - cf2
         taus     = taub_sp(i)                 ! momentum flux for i-mode w/o  rhoi(ksrc)
         kzw      = Bn(ksrc) / Ch(i)           ! ch(i) > 0.  Cx(i) < 0. critica
         etws     = taus*kzw / kxw
         rtaus(i) = taus*rhoi(ksrc)
!         
        IF( Cx <= cd_ulim .or. Cdf2 <= cd_ulim2) THEN
          Ind_out(i, ksrc) =-1                ! -1 - diagnostic index for critical levels
          cycle Loop_GW                       !  got to the next mode of GW-spectra
        ELSE
!
          kzw2 = Bn2(ksrc)/Cdf2 - rhp4 - kx2_nh
!
          if (kzw2 <= 0.) then
            Ind_out(i, ksrc) =-2               ! -2 - diagnostic index for reflected waves
            cycle  Loop_GW                     ! no wave reflection in GW-LSD scheme
          endif
     
          kzw  = sqrt(kzw2)
          kzw3 = kzw2*kzw
          etws = taus*kzw/kxw
!
! Here Linsat == Fr_critical
!
          Cx2sat = Linsat2*Cdf2
          if (etws >= cx2sat) then
               Kds  = kxw*Cx*rhp2/kzw3
               etws = cx2sat
               taus = etws*kxw/kzw
               Ind_out(i, ksrc) =-3           ! -3 - dignostic index for saturated waves
          endif
!        
          betadis = cdf2/(Cx*Cx+cf2)
          betaM   = 1.0 /(1.0+betadis)
          betaT   = 1.0 - BetaM
!
          Cxs = Cx
          kzs = kzw
!         keffs =  (kvg(ksrc)+kds)*iPr_turb*.5*khp  	                
!         sp_kth(ksrc, i)  = rhoi(ksrc)*keffs*(Tm(ksrc)+Tm(ksrc-1))
          rtaus(i) = taus*rhoi(ksrc)
          sp_tau(ksrc, i)  = rtaus(i)
          sp_etot(ksrc, i) = etws
          sp_mkz(ksrc, i)  = kzw
          sp_ek(ksrc, i)   = etws*betam
          sp_ep(ksrc, i)   = etws*betaT         ! can be transferred to (T'**2) T-rms         
             
!
        ENDIF   ! vertical propagation of i-mode to the next upper layer = (ksrc+1)
!
!  Loop_Zint  ..................................   VERTICAL "INTERFACE" LOOP from ksrc => ktop_GW
!               
  Loop_Zi: do k=ksrc+1, levs
!
           Cx   = ch(i)-Um(k)                      ! Um(k) is defined at the interface pressure levels
           Cdf2 = Cx*Cx -cf2
           if( Cx <= cd_ulim .or. Cdf2 <= 0.) then
             Ind_out(i, k) =-1                     ! 1 - diagnostic index for critical levels       
                                                   ! print*,'crit level C-U ',int(Cx),int(sqrt(cf2)),' Um ',Um(k)
             cycle Loop_GW
           endif

           cdf1 =sqrt(Cdf2)
           wdop2 = (kxw*Cx)* (kxw*Cx)
           kzw2 = (Bn2(k)-wdop2)/Cdf2  - rhp4 - kx2_nh    ! full lin DS-NIGW (N2-wd2)*k2=(m2+k2+[1/2H]^2)*(wd2-f2) 
 
              if (kzw2 < 0.) then  
               Ind_out(i, k) =-2                          ! 2 - diagnostic index for reflected waves      
               cycle Loop_GW 
              endif
            kzw = sqrt(kzw2)
            kzw3 =kzw2*kzw
!
            keff_m =  kvg(k)*kzw2 + kion(k)
!           keff_t = kturb(k)*iPr_turb + kmol(k)*iPr_mol
            keff_t = ktg(k)*kzw2  + krad(k)
!       
!
           betadis = cdf2 / (Cx*Cx+cf2)
           betaM   = 1.0 / (1.0+betadis)
           betaT   = 1.0 - BetaM

!
!imaginary frequencies of momentum and heat with "kds at (k-1) level"
!
         wfiM = kds*kzw2*F_kds             + keff_m
         wfiT = kds*iPr_ktgw*F_kds * kzw2  + keff_t
!
            wfdM = wfiM/(kxw*Cdf1)*BetaM
            wfdT = wfiT/(kxw*Cx)*BetaT
! exp-l: "kzi*dz"
            kzi = 2.*kzw*(wfdM+wfdT)*dzpi(k)             ! 2-factor energy-momentum (U')^2
!-------------------------------------------------------
! dissipative factor: Fdis
! we can replace WKB-solver by Numerical integration of
! tau_gw == etot_gw/kzw*kxw
! d(rho*tau_gw) = -kdis*rho*tau_gw
!     |tau_gw| <= |tau_gwsat|  
!                  linear limit for single mode
!     generalization for the "broad" spectra
!     or treating single mode breaking
!        over finite "vertical"-depth with "efficiency"
! Now: time-step + hor-l scale
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Fdis  = exp(-kzi) 
!          
!
! dissipative "wave rms" by WKB 
!
            etwk = etws*rhoi(k-1)/rhoi(k)*Fdis*kzw/kzs
! 
            Cx2sat = Linsat2*Cdf2
!
! Linear saturation
!
            if (etwk.ge.cx2sat) then 
                 
                Ind_out(i, k) =-3                    ! 3 - dignostic index for saturated waves
!                                                    !    saturate energy and "trigger" keddy  <kds>                
                etw_lin = etwk     
                etwk = cx2sat  
                Kds  = kxw*Cdf1*rhp2/kzw3 
                tauk = etwk*kxw/kzw
               
!===================================================================================     
!  WAM/case with high Kds   tau_lin = (etw_lin-etwk)*kxw/kzw !tau_loss by sat theory        
!      Lzsat = 6,28/kzw       Zsat1 = Zi(k)-.5*Lzsat
!                             Zsat2 = Zi(k)+.5*Lzsat   
! in WAM triggering from "kds = 0 m2/s" => "200 m2/s" for Lzw ~ 10 km
!
!             call sat_domain(zi, Zsat1, Zsat2, pver, ksat1, ksat2)
!
! to avoid it do the new diss-n factor with eddy "kds"  added to the
! background keff_m and keff_t
!
! can be taken out for the strato-mesosphere in GFS
!             wfiM = kds*kzw2 + keff_m
!             wfiT = kds*iPr_ktgw * kzw2  +keff_t                  
!             wfdM = wfiM/(kxw*Cdf1)*BetaM
!             wfdT = wfiT/(kxw*Cx)*BetaT
!             kzi = 2.*kzw*(wfdM+wfdT)*dzpi(k)
!             Fdisat  = exp(-kzi)    
!            etwk = etws*rhoi(k-1)/rhoi(k)*Fdis*(kzw/Kzs)
!  updated breaking in the Lzsat-domain:  zsat1 < zi < zsat2
! =================================================================================       
            else
              kds = 0.0  
              tauk = etwk*kxw/kzw               ! <u'w'> = Ekin*kx/kz
            ENDIF
!--------------------------------------
!
! Fill in spectral arrays(levs, nw)
!
!--------------------------------------
            sp_ked(k,i)  = kds                            ! defined at interfaces
            sp_tau(k, i) = tauk*rhoi(k)                   ! defined at interfaces   
       
!            keff = (kds + kvg(k))*iPr_turb*0.5*KHP      
!            sp_kth(k, i)  = rhoi(k)*keff*(Tm(k)+Tm(k-1))  ! defined at mid-layers  

            sp_etot(k, i) = etwk                          ! defined at interfaces
            sp_mkz(k, i)  = kzw                           ! defined at interfaces
            sp_ek(k, i)   = etwk*betam                    ! defined at interfaces
            sp_ep(k, i)   = etwk*betaT                    ! can be transferred to (T'**2)
!
!
            if (sp_tau(k,i) > sp_tau(k-1,i)) then
                sp_tau(k,i) = sp_tau(k-1,i)               ! prevent "possible" numerical "noise"
            endif
!
! updates for "eps and keff" from 
!
            rab1 =.5*(cx+cxs)*dzirho(k)                    
! heating
! due to wave dissipation
!
            sp_eps(k,i) = rab1*(sp_tau(k-1,i)- sp_tau(k,i)) ! defined at mid-layers 
!                                              
! cooling term due to eddy heat conduction =0 if Keff_cond =>0, 
!  usually updated by 1D-heat implict tridiagonal solver
!  explicit local solver ---->sp_kth(k,i) = Kt*(dT/dz+ R/Cp*T/Hp~>g/cp)
!
!        sp_eps(k,i)=sp_eps(k,i)+dzirho(k)*(sp_kth(k,i)- sp_kth(k-1,i)) 
!
            kzs  = kzw
            cxs  = cX
            taus = tauk
            etws = etwk             
!           keffs = keff             

        enddo Loop_Zi              ! ++++++++++++++ vertical layer
!
!  ................................!    stop ' in solver single-mode'
!       
      enddo Loop_GW               ! i-mode of GW-spectra
!
       sum_rtaus =sum(rtaus)       ! total momentum flux at k=ksrc

!       print *, sum_rtaus, ' tau-src ', nint(zi(ksrc)*1.e-3)
!       print *, maxval(ch), minval(ch), ' Ch ', ngwv, ' N-modes '
!
!==============================================================================
! Perform spectral integartion (sum) & apply "efficiency/inremittency" factors
!
!       eff_factor:  ~ 1./[number of modes in 1-direction of model columns]
!       
!==============================================================================
          do k=ksrc, levs

           ked(k) =0.
           Eps(k) = 0.
           Tau(k) = 0.
           swg_et(k) =0.
           swg_ep(k) =0.
           swg_ek(k) =0.

           do i=1,nw
           Ked(k) = Ked(k)+sp_ked(k,i) 
           Eps(k) = Eps(k)+sp_eps(k,i)  
           Tau(k) = Tau(k)+sp_tau(k,i)  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! GW-energy + GW-en flux  ~ Cgz*E, diagnostics-only
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           swg_et(k) = swg_et(k)+sp_etot(k,i)           !*eff_fact 
           swg_ep(k) = swg_ep(k)+sp_ep(k,i)             !*eff_fact 
           swg_ek(k) = swg_ek(k)+sp_ek(k,i)             !*eff_fact          
           enddo

          enddo  
! fill in below the "source" level ..... [1:ksrc-1]
!
         do k=1, ksrc-1
! no loss of the total momentum flux
           ked(k)  =0.
           eps(k) = 0.
           tau(k) = tau(ksrc)
! lin-theory   diagnostics-only
           swg_et(k) =swg_et(ksrc)*rhoi(ksrc)/rhoi(k) 
           swg_ep(k) =swg_ep(ksrc)*rhoi(ksrc)/rhoi(k) 
           swg_ek(k) =swg_ek(ksrc)*rhoi(ksrc)/rhoi(k) 
         enddo
!
         RETURN
!
! diagnostics below
!
345      FORMAT(2x, F8.2, 4(2x, F10.3), 2x, F8.2)
        if (dbg_print) then           
            print *
            print *,  ' Zkm        EK m2/s2      Ked m2/s    Eps m2/s3      tau-Mpa '
            do k=ksrc, levs
!                Fd1 = maxval(Fdis_modes(1:nw,k))
!                Fd2 = minval(Fdis_modes(1:nw,k))
            write(6, 345) Zi(k)*1.e-3, sqrt(swg_ek(k)), Ked(k), Eps(k), Tau(k)*1.e3, Um(k)  !, Fd1, Fd2
            enddo
           print *
           write(6,*) nw , ' nwaves-linsat '
           write(6,*) maxval(sp_ked), minval(sp_ked), 'ked '
           write(6,*) maxval(sp_tau), minval(sp_tau), 'sp_tau '
          !pause
        endif
          
!
         end subroutine ugwp_lsatdis_az1
!	 
         subroutine ugwp_limit_1d(ax, ay,eps, ked, levs) 
           use cires_ugwp_module, only : max_kdis, max_eps, max_axyz
           implicit none
           integer  :: levs
	   real, dimension(levs)   :: ax, ay,eps
	   real, dimension(levs+1) :: ked
	   real, parameter         :: xtiny = 1.e-30
           where (abs(ax)  > max_axyz ) ax  = ax/abs(ax+xtiny)*max_axyz
           where (abs(ay)  > max_axyz ) ay  = ay/abs(ay+xtiny)*max_axyz
           where (abs(eps) > max_eps  ) eps = eps/abs(eps+xtiny)*max_eps
           where (ked > max_kdis  )     ked =  max_kdis
         end subroutine ugwp_limit_1d

!>  \file sfc_diff.f
!! this file contains the surface roughness length formulation based on
!! the surface sublayer scheme in
!! zeng and dickinson (1998) \cite zeng_and_dickinson_1998.

!> this module contains the ccpp-compliant gfs surface layer scheme for
!!hafs.
      module module_sfc_diff

      use machine , only : kind_phys
      
      use physcons, grav => con_g
      real (kind=kind_phys), parameter :: ca=.4  ! ca - von karman constant

      contains

      subroutine sfc_diff_init
      end subroutine sfc_diff_init

      subroutine sfc_diff_finalize
      end subroutine sfc_diff_finalize

      subroutine sfc_diff_hafs(im,ps,u1,v1,t1,q1,z1,                    &!intent(in)
     &                    prsl1,prslki,ddvel,                           &!intent(in)
     &                    sigmaf,vegtype,shdmax,ivegsrc,                &!intent(in)
     &                    z0pert,ztpert,                                &! mg, sfc-perts !intent(in)
     &                    flag_iter,redrag,                             &!intent(in)
     &                    u10m,v10m,sfc_z0_type,                        &!wang,z0 type !intent(in)
     &                    wet,dry,icy,                                  &!intent(in)
     &                    tskin, tsurf, snwdph, z0rl, ustar,
     &                    cm, ch, rb, stress, fm, fh, fm10, fh2,
     &                    wind)                                         &!intent(out)
!
      use funcphys, only : fpvs
      use physcons, rvrdm1 => con_fvirt
     &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
! 1 - land, 2 - ice, 3 - water
! --------  -------- ---------
      integer, intent(in) :: im, ivegsrc
      integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean

      integer, dimension(im), intent(in) :: vegtype

      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
      logical, dimension(im), intent(in) :: flag_iter, dry, wet, icy

      real(kind=kind_phys), dimension(im), intent(in)    :: u10m,v10m
      real(kind=kind_phys), dimension(im), intent(in)    ::
     &                    ps,u1,v1,t1,q1,z1,prsl1,prslki,ddvel,         &
     &                    sigmaf,shdmax,                                &
     &                    z0pert,ztpert ! mg, sfc-perts       
      real(kind=kind_phys), dimension(im,3), intent(in)    ::           &
     &                    tskin, tsurf, snwdph

      real(kind=kind_phys), dimension(im,3), intent(inout) ::           &
     &                       z0rl, ustar

! 1 - land, 2 - ice, 3 - water
! --------  -------- ---------
      real(kind=kind_phys), dimension(im,3), intent(out)   ::           &
     &                       cm, ch, rb, stress, fm, fh, fm10, fh2
      real(kind=kind_phys), dimension(im),   intent(out)   :: wind
!
!     locals
!
      real(kind=kind_phys), dimension(im)                :: wind10m

      integer   i
!
      real(kind=kind_phys) :: qs1,  rat, thv1, restar,                  &
     &                      czilc, tem1, tem2

      real(kind=kind_phys) :: tvs_ocn,  tvs_lnd,  tvs_ice,              &
     &                         z0_ocn,   z0_lnd,   z0_ice,              &
     &                      z0max_ocn,z0max_lnd,z0max_ice,              &
     &                      ztmax_ocn,ztmax_lnd,ztmax_ice
!
      real(kind=kind_phys), parameter ::
     &              charnock=.014, z0s_max=.317e-2,                     & ! a limiting value at high winds over sea
     &              vis=1.4e-5, rnu=1.51e-5, visi=1.0/vis,              &
     &              log01=log(0.01), log05=log(0.05), log07=log(0.07)

!     parameter (charnock=.014,ca=.4)!c ca is the von karman constant
!     parameter (alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041)
!     parameter (a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5)

!     real(kind=kind_phys) aa1,bb1,bb2,cc,cc1,cc2,arnu
!     parameter (aa1=-1.076,bb1=.7045,cc1=-.05808)
!     parameter (bb2=-.1954,cc2=.009999)
!     parameter (arnu=.135*rnu)
!
!    z0s_max=.196e-2 for u10_crit=25 m/s
!    z0s_max=.317e-2 for u10_crit=30 m/s
!    z0s_max=.479e-2 for u10_crit=35 m/s
!
! mbek -- toga-coare flux algorithm
!     parameter (rnu=1.51e-5,arnu=0.11*rnu)
!
!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals, wind is wind speed, 
!  surface roughness length is converted to m from cm
!

!       write(0,*)'in sfc_diff, sfc_z0_type=',sfc_z0_type

      do i=1,im
        ztmax_ocn = 0.0 ; ztmax_lnd = 0.0 ; ztmax_ice = 0.0
        tvs_lnd   = 0.0 ; tvs_ice   = 0.0 ; tvs_ocn   = 0.0

        wind10m(i) = max(sqrt( u10m(i)*u10m(i) + v10m(i)*v10m(i)),      &
     &                 1.0)


        if(flag_iter(i)) then 
          wind(i) = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))                 &
     &                + max(0.0, min(ddvel(i), 30.0)), 1.0)
          tem1    = 1.0 + rvrdm1 * max(q1(i),1.e-8)
          thv1    = t1(i) * prslki(i) * tem1
          if (dry(i)) tvs_lnd = 0.5 * (tsurf(i,1)+tskin(i,1)) * tem1
          if (icy(i)) tvs_ice = 0.5 * (tsurf(i,2)+tskin(i,2)) * tem1
          if (wet(i)) tvs_ocn = 0.5 * (tsurf(i,3)+tskin(i,3)) * tem1

          qs1     = fpvs(t1(i))
          qs1     = max(1.0e-8, eps * qs1 / (prsl1(i) + epsm1 * qs1))

          z0_lnd    = 0.01 * z0rl(i,1)
          z0max_lnd = max(1.0e-6, min(z0_lnd,z1(i)))
          z0_ice    = 0.01 * z0rl(i,2)
          z0max_ice = max(1.0e-6, min(z0_ice,z1(i)))
          z0_ocn    = 0.01 * z0rl(i,3)
          z0max_ocn = max(1.0e-6, min(z0_ocn,z1(i)))

!  compute stability dependent exchange coefficients
!  this portion of the code is presently suppressed
!

          if (wet(i)) then ! some open ocean
            ustar(i,3) = sqrt(grav * z0_ocn / charnock)

!**  test xubin's new z0

!           ztmax  = z0max

            restar = max(ustar(i,3)*z0max_ocn*visi, 0.000001)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat    = min(7.0, 2.67 * sqrt(sqrt(restar)) - 2.57)
            ztmax_ocn  = z0max_ocn * exp(-rat)

           if (sfc_z0_type == 6) then
             call znot_t_v6(wind10m(i),ztmax_ocn)   ! 10-m wind,m/s, ztmax(m)
           else if (sfc_z0_type == 7) then
             call znot_t_v7(wind10m(i),ztmax_ocn)   ! 10-m wind,m/s, ztmax(m)
           else if (sfc_z0_type .ne. 0) then
             write(0,*)'no option for sfc_z0_type=',sfc_z0_type
             stop
           endif

          endif ! open ocean

          if (dry(i) .or. icy(i)) then ! over land or sea ice
!** xubin's new z0  over land and sea ice
            tem1 = 1.0 - shdmax(i)
            tem2 = tem1 * tem1
            tem1 = 1.0  - tem2
          
            if( ivegsrc == 1 ) then

              if (vegtype(i) == 10) then
                z0max_lnd = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype(i) == 6) then
                z0max_lnd = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype(i) == 7) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max_lnd = 0.01
              elseif (vegtype(i) == 16) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max_lnd = 0.01
              else
                z0max_lnd = exp( tem2*log01 + tem1*log(z0max_lnd) )
              endif

            elseif (ivegsrc == 2 ) then

                if (vegtype(i) == 7) then
                  z0max_lnd = exp( tem2*log01 + tem1*log07 )
                elseif (vegtype(i) == 8) then
                  z0max_lnd = exp( tem2*log01 + tem1*log05 )
                elseif (vegtype(i) == 9) then
!                 z0max = exp( tem2*log01 + tem1*log01 )
                  z0max_lnd = 0.01
                elseif (vegtype(i) == 11) then
!                 z0max = exp( tem2*log01 + tem1*log01 )
                  z0max_lnd = 0.01
                else
                  z0max_lnd = exp( tem2*log01 + tem1*log(z0max_lnd) )
                endif

            endif ! over land or sea ice

            z0max_ice = z0max_lnd

! mg, sfc-perts: add surface perturbations to z0max over land
            if (dry(i) .and. z0pert(i) /= 0.0 ) then
              z0max_lnd = z0max_lnd * (10.**z0pert(i))
            endif
 
            z0max_lnd = max(z0max_lnd,1.0e-6)
            z0max_ice = max(z0max_ice,1.0e-6)

!           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
            czilc = 0.8

            tem1 = 1.0 - sigmaf(i)
            ztmax_lnd = z0max_lnd*exp( - tem1*tem1
     &                     * czilc*ca*sqrt(ustar(i,1)*(0.01/1.5e-05)))
            ztmax_ice = z0max_ice*exp( - tem1*tem1
     &                     * czilc*ca*sqrt(ustar(i,2)*(0.01/1.5e-05)))


! mg, sfc-perts: add surface perturbations to ztmax/z0max ratio over land
            if (dry(i) .and. ztpert(i) /= 0.0) then
              ztmax_lnd = ztmax_lnd * (10.**ztpert(i))
            endif


          endif       ! end of if(sfctype flags) then

          ztmax_ocn  = max(ztmax_ocn,1.0e-6)
          ztmax_lnd  = max(ztmax_lnd,1.0e-6)
          ztmax_ice  = max(ztmax_ice,1.0e-6)

! bwg begin "stability" block, 2019-03-23
          if (wet(i)) then ! some open ocean
            call stability                                              &
!  ---  inputs:                                                  
     &       (z1(i),snwdph(i,3),thv1,wind(i),                           &
     &        z0max_ocn,ztmax_ocn,tvs_ocn,                              &
!  ---  outputs:
     &        rb(i,3), fm(i,3), fh(i,3), fm10(i,3), fh2(i,3),           &
     &        cm(i,3), ch(i,3), stress(i,3), ustar(i,3))               
          endif ! open ocean points

          if (dry(i)) then ! some land
            call stability                                              &
!  ---  inputs:                                                  
     &       (z1(i),snwdph(i,1),thv1,wind(i),                           &
     &        z0max_lnd,ztmax_lnd,tvs_lnd,                              &
!  ---  outputs:
     &        rb(i,1), fm(i,1), fh(i,1), fm10(i,1), fh2(i,1),           &
     &        cm(i,1), ch(i,1), stress(i,1), ustar(i,1))
          endif ! dry points

          if (icy(i)) then ! some ice
            call stability                                              &
!  ---  inputs:                                                  
     &       (z1(i),snwdph(i,2),thv1,wind(i),                           &
     &        z0max_ice,ztmax_ice,tvs_ice,
!  ---  outputs:
     &        rb(i,2), fm(i,2), fh(i,2), fm10(i,2), fh2(i,2),           &
     &        cm(i,2), ch(i,2), stress(i,2), ustar(i,2))                
          endif ! icy points

! bwg: everything from here to end of subroutine was after
!      the stuff now put into "stability"

!
!  update z0 over ocean
!
          if (wet(i)) then
            z0_ocn = (charnock / grav) * ustar(i,3) * ustar(i,3)

! mbek -- toga-coare flux algorithm
!           z0 = (charnock / grav) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!           cc = ustar(i) * z0 / rnu
!           pp = cc / (1. + cc)
!           ff = grav * arnu / (charnock * ustar(i) ** 3)
!           z0 = arnu / (ustar(i) * ff ** pp)

            if (redrag) then
              z0rl(i,3) = 100.0 * max(min(z0_ocn, z0s_max), 1.e-7)
            else
              z0rl(i,3) = 100.0 * max(min(z0_ocn,.1), 1.e-7)
            endif

             if (sfc_z0_type == 6) then   ! wang
               call znot_m_v6(wind10m(i),z0_ocn)   ! wind, m/s, z0, m
               z0rl(i,3) = 100.0 * z0_ocn          ! cm
             endif              !wang
             if (sfc_z0_type == 7) then   ! wang
               call znot_m_v7(wind10m(i),z0_ocn)   ! wind, m/s, z0, m
               z0rl(i,3) = 100.0 * z0_ocn          ! cm
             endif              !wang


          endif              ! end of if(open ocean)
        endif                ! end of if(flagiter) loop
      enddo

      return
      end subroutine sfc_diff


!----------------------------------------
      subroutine stability
!........................................
!  ---  inputs:
     &     ( z1, snwdph, thv1, wind, z0max, ztmax, tvs,                 &
!  ---  outputs:
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar)
!-----

!  ---  inputs:
      real(kind=kind_phys), intent(in) ::                               &
     &       z1, snwdph, thv1, wind, z0max, ztmax, tvs

!  ---  outputs:
      real(kind=kind_phys), intent(out) ::                              &
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar

!  ---  locals:
      real(kind=kind_phys), parameter :: alpha=5., a0=-3.975,           &
     &              a1=12.32, alpha4=4.0*alpha,                         &
     &              b1=-7.755,  b2=6.041,  alpha2=alpha+alpha, beta=1.0,&
     &              a0p=-7.941, a1p=24.75, b1p=-8.705, b2p=7.899,       &
     &              ztmin1=-999.0

      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,    &
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2,     &
     &                     z1i,                                         &
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,       &
     &                     hl110,  hlt,    hltinf, olinf,               &
     &                     tem1,   tem2, ztmax1

          z1i = 1.0 / z1

          tem1   = z0max/z1
          if (abs(1.0-tem1) > 1.0e-6) then
            ztmax1 = - beta*log(tem1)/(alpha2*(1.-tem1))
          else
            ztmax1 = 99.0
          endif
          if( z0max < 0.05 .and. snwdph < 10.0 ) ztmax1 = 99.0

!  compute stability indices (rb and hlinf)

          dtv     = thv1 - tvs
          adtv    = max(abs(dtv),0.001)
          dtv     = sign(1.,dtv) * adtv
          rb      = max(-5000.0, (grav+grav) * dtv * z1
     &            / ((thv1 + tvs) * wind * wind))
          tem1    = 1.0 / z0max
          tem2    = 1.0 / ztmax
          fm      = log((z0max+z1) * tem1)
          fh      = log((ztmax+z1) * tem2)
          fm10    = log((z0max+10.)* tem1)
          fh2     = log((ztmax+2.) * tem2)
          hlinf   = rb * fm * fm / fh
          hlinf   = min(max(hlinf,ztmin1),ztmax1)
!
!  stable case
!
          if (dtv >= 0.0) then
            hl1 = hlinf
            if(hlinf > .25) then
              tem1   = hlinf * z1i
              hl0inf = z0max * tem1
              hltinf = ztmax * tem1
              aa     = sqrt(1. + alpha4 * hlinf)
              aa0    = sqrt(1. + alpha4 * hl0inf)
              bb     = aa
              bb0    = sqrt(1. + alpha4 * hltinf)
              pm     = aa0 - aa + log( (aa + 1.)/(aa0 + 1.) )
              ph     = bb0 - bb + log( (bb + 1.)/(bb0 + 1.) )
              fms    = fm - pm
              fhs    = fh - ph
              hl1    = fms * fms * rb / fhs
              hl1    = min(max(hl1, ztmin1), ztmax1)
            endif
!
!  second iteration
!
            tem1  = hl1 * z1i
            hl0   = z0max * tem1
            hlt   = ztmax * tem1
            aa    = sqrt(1. + alpha4 * hl1)
            aa0   = sqrt(1. + alpha4 * hl0)
            bb    = aa
            bb0   = sqrt(1. + alpha4 * hlt)
            pm    = aa0 - aa + log( (1.0+aa)/(1.0+aa0) )
            ph    = bb0 - bb + log( (1.0+bb)/(1.0+bb0) )
            hl110 = hl1 * 10. * z1i
            hl110 = min(max(hl110, ztmin1), ztmax1)
            aa    = sqrt(1. + alpha4 * hl110)
            pm10  = aa0 - aa + log( (1.0+aa)/(1.0+aa0) )
            hl12  = (hl1+hl1) * z1i
            hl12  = min(max(hl12,ztmin1),ztmax1)
!           aa    = sqrt(1. + alpha4 * hl12)
            bb    = sqrt(1. + alpha4 * hl12)
            ph2   = bb0 - bb + log( (1.0+bb)/(1.0+bb0) )
!
!  unstable case - check for unphysical obukhov length
!
          else                          ! dtv < 0 case
            olinf = z1 / hlinf
            tem1  = 50.0 * z0max
            if(abs(olinf) <= tem1) then
              hlinf = -z1 / tem1
              hlinf = min(max(hlinf,ztmin1),ztmax1)
            endif
!
!  get pm and ph
!
            if (hlinf >= -0.5) then
              hl1   = hlinf
              pm    = (a0  + a1*hl1)  * hl1   / (1.+ (b1+b2*hl1)  *hl1)
              ph    = (a0p + a1p*hl1) * hl1   / (1.+ (b1p+b2p*hl1)*hl1)
              hl110 = hl1 * 10. * z1i
              hl110 = min(max(hl110, ztmin1), ztmax1)
              pm10  = (a0 + a1*hl110) * hl110 / (1.+(b1+b2*hl110)*hl110)
              hl12  = (hl1+hl1) * z1i
              hl12  = min(max(hl12, ztmin1), ztmax1)
              ph2   = (a0p + a1p*hl12) * hl12 / (1.+(b1p+b2p*hl12)*hl12)
            else                       ! hlinf < 0.05
              hl1   = -hlinf
              tem1  = 1.0 / sqrt(hl1)
              pm    = log(hl1) + 2. * sqrt(tem1) - .8776
              ph    = log(hl1) + .5 * tem1 + 1.386
!             pm    = log(hl1) + 2.0 * hl1 ** (-.25) - .8776
!             ph    = log(hl1) + 0.5 * hl1 ** (-.5) + 1.386
              hl110 = hl1 * 10. * z1i
              hl110 = min(max(hl110, ztmin1), ztmax1)
              pm10  = log(hl110) + 2.0 / sqrt(sqrt(hl110)) - .8776
!             pm10  = log(hl110) + 2. * hl110 ** (-.25) - .8776
              hl12  = (hl1+hl1) * z1i
              hl12  = min(max(hl12, ztmin1), ztmax1)
              ph2   = log(hl12) + 0.5 / sqrt(hl12) + 1.386
!             ph2   = log(hl12) + .5 * hl12 ** (-.5) + 1.386
            endif

          endif          ! end of if (dtv >= 0 ) then loop
!
!  finish the exchange coefficient computation to provide fm and fh
!
          fm        = fm - pm
          fh        = fh - ph
          fm10      = fm10 - pm10
          fh2       = fh2 - ph2
          cm        = ca * ca / (fm * fm)
          ch        = ca * ca / (fm * fh)
          tem1      = 0.00001/z1
          cm        = max(cm, tem1)
          ch        = max(ch, tem1)
          stress    = cm * wind * wind
          ustar     = sqrt(stress)

      return
!.................................
      end subroutine stability
!---------------------------------


!! add fitted z0,zt curves for hurricane application (used in hwrf/hmon)
!! weiguo wang, 2019-0425

       subroutine znot_m_v6(uref,znotm)
       implicit none
! calculate areodynamical roughness over water with input 10-m wind
! for low-to-moderate winds, try to match the cd-u10 relationship from coare v3.5 (edson et al. 2013)
! for high winds, try to fit available observational data
!
! bin liu, noaa/ncep/emc 2017
! 
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!

       real (kind=kind_phys), intent(in) :: uref
       real (kind=kind_phys), intent(out):: znotm
       real             :: p13, p12, p11, p10
       real             :: p25, p24, p23, p22, p21, p20
       real             :: p35, p34, p33, p32, p31, p30
       real             :: p40

       real(kind=kind_phys), parameter :: p13 = -1.296521881682694e-02, &
     & p12 =  2.855780863283819e-01 p11 = -1.597898515251717e+00,       &
     & p10 = -8.396975715683501e+00,                                    & 

     & p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,      &
     & p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,      &
     & p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,      &

     & p35 =  1.840430200185075e-07, p34 = -2.793849676757154e-05,      &
     & p33 =  1.735308193700643e-03, p32 = -6.139315534216305e-02,      &
     & p31 =  1.255457892775006e+00, p30 = -1.663993561652530e+01,      &

     & p40 =  4.579369142033410e-04

       if (uref >= 0.0 .and.  uref <= 6.5 ) then
        znotm = exp( p10 + p11*uref + p12*uref**2 +                     &
     &          p13*uref**3)
       elseif (uref > 6.5 .and. uref <= 15.7) then
        znotm = p25*uref**5 + p24*uref**4 + p23*uref**3 +               & 
     &          p22*uref**2 + p21*uref + p20
       elseif (uref > 15.7 .and. uref <= 53.0) then
        znotm = exp( p35*uref**5 + p34*uref**4 +                        &
     &         p33*uref**3 + p32*uref**2 + p31*uref + p30 )
       elseif ( uref > 53.0) then
         znotm = p40
       else
        print*, 'wrong input uref value:',uref
       endif

       end subroutine znot_m_v6

       subroutine znot_t_v6(uref,znott)
       implicit none
! calculate scalar roughness over water with input 10-m wind
! for low-to-moderate winds, try to match the ck-u10 relationship from coare algorithm
! for high winds, try to retain the ck-u10 relationship of fy2015 hwrf
!
! bin liu, noaa/ncep/emc 2017
!
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
!

       real (kind=kind_phys), intent(in) :: uref
       real (kind=kind_phys), intent(out):: znott

       real             :: p00
       real             :: p15, p14, p13, p12, p11, p10
       real             :: p25, p24, p23, p22, p21, p20
       real             :: p35, p34, p33, p32, p31, p30
       real             :: p45, p44, p43, p42, p41, p40
       real             :: p56, p55, p54, p53, p52, p51, p50
       real             :: p60

       real(kind=kind_phys), parameter :: p00 =  1.100000000000000e-04, &

     & p15 = -9.144581627678278e-10, p14 =  7.020346616456421e-08,      &
     & p13 = -2.155602086883837e-06, p12 =  3.333848806567684e-05,      &
     & p11 = -2.628501274963990e-04, p10 =  8.634221567969181e-04,      &

     & p25 = -8.654513012535990e-12, p24 =  1.232380050058077e-09,      &
     & p23 = -6.837922749505057e-08, p22 =  1.871407733439947e-06,      &
     & p21 = -2.552246987137160e-05, p20 =  1.428968311457630e-04,      &

     & p35 =  3.207515102100162e-12, p34 = -2.945761895342535e-10,      &
     & p33 =  8.788972147364181e-09, p32 = -3.814457439412957e-08,      &
     & p31 = -2.448983648874671e-06, p30 =  3.436721779020359e-05,      &

     & p45 = -3.530687797132211e-11, p44 =  3.939867958963747e-09,      &
     & p43 = -1.227668406985956e-08, p42 = -1.367469811838390e-05,      &
     & p41 =  5.988240863928883e-04, p40 = -7.746288511324971e-03,      &

     & p56 = -1.187982453329086e-13, p55 =  4.801984186231693e-11,      &
     & p54 = -8.049200462388188e-09, p53 =  7.169872601310186e-07,      &
     & p52 = -3.581694433758150e-05, p51 =  9.503919224192534e-04,      &
     & p50 = -1.036679430885215e-02,

     & p60 =  4.751256171799112e-05

       if (uref >= 0.0 .and. uref < 5.9 ) then
         znott = p00
        elseif (uref >= 5.9 .and. uref <= 15.4) then
         znott = p15*uref**5 + p14*uref**4 + p13*uref**3                &
     &          + p12*uref**2 + p11*uref + p10
        elseif (uref > 15.4 .and. uref <= 21.6) then
         znott = p25*uref**5 + p24*uref**4 + p23*uref**3                &
     &          + p22*uref**2 + p21*uref + p20
        elseif (uref > 21.6 .and. uref <= 42.2) then
         znott = p35*uref**5 + p34*uref**4 + p33*uref**3                &
     &          + p32*uref**2 + p31*uref + p30
        elseif ( uref > 42.2 .and. uref <= 53.3) then
         znott = p45*uref**5 + p44*uref**4 + p43*uref**3                &
     &          + p42*uref**2 + p41*uref + p40
        elseif ( uref > 53.3 .and. uref <= 80.0) then
         znott = p56*uref**6 + p55*uref**5 + p54*uref**4                &
     &      + p53*uref**3 + p52*uref**2 + p51*uref + p50
        elseif ( uref > 80.0) then
         znott = p60
       else
         print*, 'wrong input uref value:',uref
       endif

       end subroutine znot_t_v6


       subroutine znot_m_v7(uref,znotm)
        implicit none
! calculate areodynamical roughness over water with input 10-m wind
! for low-to-moderate winds, try to match the cd-u10 relationship from coare v3.5 (edson et al. 2013)
! for high winds, try to fit available observational data
! comparing to znot_t_v6, slightly decrease cd for higher wind speed
!
! bin liu, noaa/ncep/emc 2018
!
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!

      real (kind=kind_phys), intent(in) :: uref
      real (kind=kind_phys), intent(out):: znotm
      real             :: p13, p12, p11, p10
      real             :: p25, p24, p23, p22, p21, p20
      real             :: p35, p34, p33, p32, p31, p30
      real             :: p40

       real(kind=kind_phys), parameter :: p13 = -1.296521881682694e-02, &
     & p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,      &
     & p10 = -8.396975715683501e+00,                                    &

     & p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,      &
     & p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,      &
     & p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,      &


     & p35 =  1.897534489606422e-07, p34 = -3.019495980684978e-05,      &
     & p33 =  1.931392924987349e-03, p32 = -6.797293095862357e-02,      &
     & p31 =  1.346757797103756e+00, p30 = -1.707846930193362e+01,      &

     & p40 =  3.371427455376717e-04

       if (uref >= 0.0 .and.  uref <= 6.5 ) then
         znotm = exp( p10 + p11*uref + p12*uref**2 + p13*uref**3)
         elseif (uref > 6.5 .and. uref <= 15.7) then
         znotm = p25*uref**5 + p24*uref**4 + p23*uref**3 +              &
     &          p22*uref**2 + p21*uref + p20
         elseif (uref > 15.7 .and. uref <= 53.0) then
         znotm = exp( p35*uref**5 + p34*uref**4 + p33*uref**3 +         & 
     &          p32*uref**2 + p31*uref + p30 )
         elseif ( uref > 53.0) then
         znotm = p40
       else
        print*, 'wrong input uref value:',uref
       endif

      end subroutine znot_m_v7
      subroutine znot_t_v7(uref,znott)
       implicit none
! calculate scalar roughness over water with input 10-m wind
! for low-to-moderate winds, try to match the ck-u10 relationship from coare algorithm
! for high winds, try to retain the ck-u10 relationship of fy2015 hwrf
! to be compatible with the slightly decreased cd for higher wind speed
!
! bin liu, noaa/ncep/emc 2018
!
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
!

      real (kind=kind_phys), intent(in) :: uref
      real (kind=kind_phys), intent(out):: znott

      real             :: p00
      real             :: p15, p14, p13, p12, p11, p10
      real             :: p25, p24, p23, p22, p21, p20
      real             :: p35, p34, p33, p32, p31, p30
      real             :: p45, p44, p43, p42, p41, p40
      real             :: p56, p55, p54, p53, p52, p51, p50
      real             :: p60

      real(kind=kind_phys), parameter ::p00 =  1.100000000000000e-04,   &

     &  p15 = -9.193764479895316e-10, p14 =  7.052217518653943e-08,     &
     &  p13 = -2.163419217747114e-06, p12 =  3.342963077911962e-05,     &
     &  p11 = -2.633566691328004e-04, p10 =  8.644979973037803e-04,     &

     &  p25 = -9.402722450219142e-12, p24 =  1.325396583616614e-09,     &
     &  p23 = -7.299148051141852e-08, p22 =  1.982901461144764e-06,     &
     &  p21 = -2.680293455916390e-05, p20 =  1.484341646128200e-04,     &

     &  p35 =  7.921446674311864e-12, p34 = -1.019028029546602e-09,     &
     &  p33 =  5.251986927351103e-08, p32 = -1.337841892062716e-06,     &
     &  p31 =  1.659454106237737e-05, p30 = -7.558911792344770e-05,     &

     &  p45 = -2.694370426850801e-10, p44 =  5.817362913967911e-08,     &
     &  p43 = -5.000813324746342e-06, p42 =  2.143803523428029e-04,     &
     &  p41 = -4.588070983722060e-03, p40 =  3.924356617245624e-02,     &
      
     &  p56 = -1.663918773476178e-13, p55 =  6.724854483077447e-11,     &
     &  p54 = -1.127030176632823e-08, p53 =  1.003683177025925e-06,     &
     &  p52 = -5.012618091180904e-05, p51 =  1.329762020689302e-03,     &
     &  p50 = -1.450062148367566e-02,

     &  p60 =  6.840803042788488e-05

        if (uref >= 0.0 .and. uref < 5.9 ) then
            znott = p00
         elseif (uref >= 5.9 .and. uref <= 15.4) then
           znott = p15*uref**5 + p14*uref**4 + p13*uref**3 + 
     &             p12*uref**2 + p11*uref + p10
         elseif (uref > 15.4 .and. uref <= 21.6) then
           znott = p25*uref**5 + p24*uref**4 + p23*uref**3 + 
     &             p22*uref**2 + p21*uref + p20
         elseif (uref > 21.6 .and. uref <= 42.6) then
           znott = p35*uref**5 + p34*uref**4 + p33*uref**3 + 
     &             p32*uref**2 + p31*uref + p30
         elseif ( uref > 42.6 .and. uref <= 53.0) then
           znott = p45*uref**5 + p44*uref**4 + p43*uref**3 + 
     &             p42*uref**2 + p41*uref + p40
         elseif ( uref > 53.0 .and. uref <= 80.0) then
           znott = p56*uref**6 + p55*uref**5 + p54*uref**4 + 
     &             p53*uref**3 + p52*uref**2 + p51*uref + p50
         elseif ( uref > 80.0) then
           znott = p60
        else
           print*, 'wrong input uref value:',uref
         endif

        end subroutine znot_t_v7


!---------------------------------
      end module module_sfc_diff

!> \File drag_suite.F90
!! This file is the  parameterization of orographic gravity wave
!! drag, mountain blocking, and form drag.

!> This module contains the CCPP-compliant orographic gravity wave dray scheme.
      module drag_suite

      contains

      subroutine drag_suite_init()
      end subroutine drag_suite_init

! \defgroup GFS_ogwd GFS Orographic Gravity Wave Drag
!> \defgroup gfs_drag_suite GFS drag_suite Main
!! \brief This subroutine includes orographic gravity wave drag,  mountain
!! blocking, and form drag.
!!
!> The time tendencies of zonal and meridional wind are altered to
!! include the effect of mountain induced gravity wave drag from
!! subgrid scale orography including convective breaking, shear
!! breaking and the presence of critical levels.
!!
!> \section arg_table_drag_suite_run Argument Table
!! \htmlinclude drag_suite_run.html
!!
!> \section gen_drag_suite GFS Orographic GWD Scheme General Algorithm
!! -# Calculate subgrid mountain blocking
!! -# Calculate orographic wave drag
!!
!! The NWP model gravity wave drag (GWD) scheme in the GFS has two
!! main components: how the surface stress is computed, and then how
!! that stress is distributed over a vertical column where it may
!! interact with the models momentum.  Each of these depends on the
!! large scale environmental atmospheric state and assumptions about
!! the sub-grid scale processes.  In Alpert GWD (1987) based on linear,
!! two-dimensional non-rotating, stably stratified flow over a mountain ridge,
!! sub-grid scale gravity wave motions are assumed which propagate away
!! from the mountain.  Described in Alpert (1987), the flux measured over
!! a "low level" vertically averaged layer, in the atmosphere defines a base
!! level flux.  "Low level" was taken to be the first 1/3 of the troposphere
!! in the 1987 implementation.  This choice was meant to encompass a thick
!! low layer for vertical averages of the environmental (large scale) flow
!! quantities.  The vertical momentum flux or gravity wave stress in a
!! grid box due to a single mountain is given as in Pierrehumbert, (1987) (PH):
!!
!! \f$ \tau =  \frac {\rho \: U^{3}\: G(F_{r})} {\Delta X \; N } \f$
!!
!! emetic \f$ \Delta X \f$ is a grid increment, N is the Brunt Viasala frequency
!!
!!
!! \f$ N(\sigma) = \frac{-g \: \sigma \:
!!  \frac{\partial\Theta}{\partial\sigma}}{\Theta \:R \:T} \f$
!!
!! The environmental variables are calculated from a mass weighted vertical
!! average over a base layer.  G(Fr) is a monotonically increasing
!! function of Froude number,
!!
!! \f$ F_{r} = \frac{N h^{'}}{U} \f$
!!
!! where U is the wind speed calculated as a mass weighted vertical average in
!! the base layer, and  h', is the vertical displacement caused by the orography
!! variance.  An effective mountain length for the gravity wave processes,
!!
!! \f$ l^{*} =  \frac{\Delta X}{m} \f$
!!
!! where m is the number of mountains in a grid box, can then
!! be defined to obtain the form of the base level stress
!!
!!
!! \f$ \tau =  \frac {\rho \: U^{3} \: G(F_{r})} {N \;l^{*}} \f$
!!
!! giving the stress induced from the surface in a model grid box.
!!   PH gives the form for the function G(Fr) as
!!
!!
!! \f$ G(F_{r}) = \bar{G}\frac{F^{2}_{r}}{F^{2}_{r}\: + \:a^{2}} \f$
!!
!! Where \f$ \bar{G}  \f$  is an order unity non-dimensional saturation
!! flux set to 1  and 'a' is a function of the mountain aspect ratio also
!!set to 1 in the 1987 implementation of the GFS GWD.  Typical values of
!! U=10m/s, N=0.01 1/s, l*=100km, and a=1, gives a flux of 1 Pascal and
!! if this flux is made to go to zero linearly with height then the
!! decelerations would be about 10/m/s/day which is consistent with
!! observations in PH.
!!
!!
!! In Kim, Moorthi, Alpert's (1998, 2001) GWD currently in GFS operations,
!! the GWD scheme has the same physical basis as in Alpert (1987) with the addition
!! of enhancement factors for the amplitude, G, and mountain shape details
!! in G(Fr) to account for effects from the mountain blocking.  A factor,
!! E m', is an enhancement factor on the stress in the Alpert '87 scheme.
!!  The E ranges from no enhancement to an upper limit of 3, E=E(OA)[1-3],
!!  and is a function of OA, the Orographic Asymmetry defined in KA (1995) as
!!
!! Orographic Asymmetry (OA) = \f$  \frac{ \bar{x} \; - \;
!!  \sum\limits_{j=1}^{N_{b}} x_{j} \; n_{j} }{\sigma_{x}} \f$
!!
!! where Nb is the total number of bottom blocks in the mountain barrier,
!! \f$ \sigma_{x} \f$ is the standard deviation of the horizontal distance defined by
!!
!! \f$ \sigma_{x} = \sqrt{ \frac{\sum\limits_{j=1}^{N_{b}}
!! \; (x_{j} \; - \; \bar{x} )^2}{N_{x}} } \f$
!!
!!
!! where Nx is the number of grid intervals for the large scale domain being
!! considered. So the term, E(OA)m'/  \f$ \Delta X \f$ in Kim's scheme represents
!! a multiplier on G shown in Alpert's eq (1), where m' is the number of mountains
!! in a sub-grid scale box. Kim increased the complexity of m' making it a
!! function of the fractional area of the sub-grid mountain and the asymmetry
!! and convexity statistics which are found from running a gravity wave
!!  model for a large number of cases:
!!
!! \f$ m^{'} = C_{m} \Delta X \left[  \frac{1 \; + \;
!!  \sum\limits_{x} L_{h} }{\Delta X}  \right]^{OA+1}   \f$
!!
!! Where, according to Kim,  \f$ \sum \frac{L_{h}}{\Delta X} \f$  is
!! the fractional area covered by the subgrid-scale orography higher than
!! a critical height  \f$ h_{c} = Fr_{c} U_{0}/N_{0} \f$ , over the
!! "low level" vertically averaged layer, for a grid box with the interval
!! \f$ \Delta X \f$.  Each \f$ L_{n}\f$  is the width of a segment of
!! orography intersection at the critical height:
!!
!! \f$  Fr_{0} = \frac{N_{0} \; h^{'}}{U_{0}}  \f$
!!
!! \f$ G^{'}(OC,Fr_{0}) = \frac{Fr_{0}^{2}}{Fr_{0}^{2} \; + \; a^{2}}  \f$
!!
!! \f$  a^{2} = \frac{C_{G}}{OC}  \f$
!!
!! \f$  E(OA, Fr_{0}) = (OA \; + \; 2)^{\delta} \f$ and \f$  \delta
!! \; = \; \frac{C_{E} \; Fr_{0}}{Fr_{c}}  \f$  where \f$ Fr_{c} \f$
!! is as in Alpert.
!!
!!
!! This represents a closed scheme, somewhat empirical adjustments
!! to the original scheme to calculate the surface stress.
!!
!! Momentum is deposited by the sub-grid scale gravity waves break due
!! to the presence of convective mixing assumed to occur when the
!! minimum Richardson number:
!!
!! Orographic Convexity (OC) = \f$  \frac{ \sum\limits_{j=1}^{N_{x}}
!!  \; (h_{j} \; - \; \bar{h})^4 }{N_{x} \;\sigma_{h}^4} \f$  ,
!!   and where  \f$ \sigma_{h} = \sqrt{ \frac{\sum\limits_{j=1}^{N_{x}}
!!  \; (h_{j} \; - \; \bar{h} )^2}{N_{x}} } \f$
!!
!! This represents a closed scheme, somewhat empirical adjustments
!! to the original scheme to calculate the surface stress.
!!
!! Momentum is deposited by the sub-grid scale gravity waves break due
!!  to the presence of convective mixing assumed to occur when
!!  the minimum Richardson number:
!!
!! \f$ Ri_{m} = \frac{Ri(1 \; - \; Fr)}{(1 \; + \; \sqrt{Ri}Fr)^2} \f$
!!
!! Is less than 1/4  Or if critical layers are encountered in a layer
!! the the momentum flux will vanish.  The critical layer is defined
!! when the base layer wind becomes perpendicular to the environmental
!!  wind.  Otherwise, wave breaking occurs at a level where the amplification
!!  of the wave causes the local Froude number or similarly a truncated
!!  (first term of the) Scorer parameter, to be reduced below a critical
!!  value by the saturation hypothesis (Lindzen,).  This is done through
!!  eq 1 which can be written as
!!
!! \f$ \tau = \rho U N k h^{'2} \f$
!!
!! For small Froude number this is discretized in the vertical so at each
!!  level the stress is reduced by ratio of the Froude or truncated Scorer
!!  parameter, \f$ \frac{U^{2}}{N^{2}} = \frac{N \tau_{l-1}}{\rho U^{3} k} \f$ ,
!!  where the stress is from the layer below beginning with that found near
!!  the surface.   The respective change in momentum is applied in
!! that layer building up from below.
!!
!! An amplitude factor is part of the calibration of this scheme which is
!! a function of the model resolution and the vertical diffusion.  This
!! is because the vertical diffusion and the GWD account encompass
!! similar physical processes. Thus, one needs to run the model over
!! and over for various amplitude factors for GWD and vertical diffusion.
!!
!! In addition, there is also mountain blocking from lift and frictional
!!  forces.  Improved integration between how the GWD is calculated and
!! the mountain blocking of wind flow around sub-grid scale orography
!! is underway at NCEP.  The GFS already has convectively forced GWD
!!  an independent process.  The next step is to test
!!
!> \section det_drag_suite GFS Orographic GWD Scheme Detailed Algorithm
!> @{
   subroutine drag_suite_run(                                           &
     &           IM,KM,dvdt,dudt,dtdt,U1,V1,T1,Q1,KPBL,                 &
     &           PRSI,DEL,PRSL,PRSLK,PHII,PHIL,DELTIM,KDT,              &
     &           var,oc1,oa4,ol4,                                       &
     &           varss,oc1ss,oa4ss,ol4ss,                               &
     &           THETA,SIGMA,GAMMA,ELVMAX,                              &
     &           dtaux2d_ls,dtauy2d_ls,dtaux2d_bl,dtauy2d_bl,           &
     &           dtaux2d_ss,dtauy2d_ss,dtaux2d_fd,dtauy2d_fd,           &
     &           dusfc,dvsfc,                                           &
     &           dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,                   &
     &           dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd,                   &
     &           slmsk,br1,hpbl,                                        &
     &           g, cp, rd, rv, fv, pi, imx, cdmbgwd, me, master,       &
     &           lprnt, ipr, rdxzb, dx, gwd_opt,                        &
     &           do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd,   &
     &           errmsg, errflg     )

!   ********************************************************************
! ----->  I M P L E M E N T A T I O N    V E R S I O N   <----------
!
!                ----- This code -----
!begin WRF code

!  this code handles the time tendencies of u v due to the effect of  mountain
!  induced gravity wave drag from sub-grid scale orography. this routine
!  not only treats the traditional upper-level wave breaking due to mountain
!  variance (alpert 1988), but also the enhanced lower-tropospheric wave
!  breaking due to mountain convexity and asymmetry (kim and arakawa 1995).
!  thus, in addition to the terrain height data in a model grid box,
!  additional 10-2d topographic statistics files are needed, including
!  orographic standard  deviation (var), convexity (oc1), asymmetry (oa4)
!  and ol (ol4). these data sets are prepared based on the 30 sec usgs orography
!  hong (1999). the current scheme was implmented as in hong et al.(2008)
!
!  Originally coded by song-you hong and young-joon kim and implemented by song-you hong
!
!  program history log:
!    2014-10-01  Hyun-Joo Choi (from KIAPS)  flow-blocking drag of kim and doyle
!                              with blocked height by dividing streamline theory
!    2017-04-06  Joseph Olson (from Gert-Jan Steeneveld) added small-scale
!                    orographic grabity wave drag:
!    2017-09-15  Joseph Olson, with some bug fixes from Michael Toy: added the
!                    topographic form drag of Beljaars et al. (2004, QJRMS)
!           Activation of each component is done by specifying the integer-parameters
!           (defined below) to 0: inactive or 1: active
!                    gwd_opt_ls = 0 or 1: large-scale  
!                    gwd_opt_bl = 0 or 1: blocking drag
!                    gwd_opt_ss = 0 or 1: small-scale gravity wave drag
!                    gwd_opt_fd = 0 or 1: topographic form drag
!    2017-09-25  Michael Toy (from NCEP GFS model) added dissipation heating
!                    gsd_diss_ht_opt = 0: dissipation heating off
!                    gsd_diss_ht_opt = 1: dissipation heating on
!    2020-08-25  Michael Toy changed logic control for drag component selection
!                    for CCPP.
!                    Namelist options:
!                    do_gsl_drag_ls_bl - logical flag for large-scale GWD + blocking
!                    do_gsl_drag_ss - logical flag for small-scale GWD
!                    do_gsl_drag_tofd - logical flag for turbulent form drag
!                    Compile-time options (same as before):
!                    gwd_opt_ls = 0 or 1: large-scale GWD
!                    gwd_opt_bl = 0 or 1: blocking drag
!
!  References:
!        Hong et al. (2008), wea. and forecasting
!        Kim and Doyle (2005), Q. J. R. Meteor. Soc.
!        Kim and Arakawa (1995), j. atmos. sci.
!        Alpert et al. (1988), NWP conference.
!        Hong (1999), NCEP office note 424.
!        Steeneveld et al (2008), JAMC
!        Tsiringakis et al. (2017), Q. J. R. Meteor. Soc.
!        Beljaars et al. (2004), Q. J. R. Meteor. Soc.
!
!  notice : comparible or lower resolution orography files than model resolution
!           are desirable in preprocess (wps) to prevent weakening of the drag
!-------------------------------------------------------------------------------
!
!  input
!        dudt (im,km)  non-lin tendency for u wind component
!        dvdt (im,km)  non-lin tendency for v wind component
!        u1(im,km) zonal wind / sqrt(rcl)  m/sec  at t0-dt
!        v1(im,km) meridional wind / sqrt(rcl) m/sec at t0-dt
!        t1(im,km) temperature deg k at t0-dt
!        q1(im,km) specific humidity at t0-dt
!        deltim  time step    secs
!        del(km)  positive increment of pressure across layer (pa)
!        KPBL(IM) is the index of the top layer of the PBL
!        ipr & lprnt for diagnostics 
!
!  output
!        dudt, dvdt    wind tendency due to gwdo
!        dTdt
!
!-------------------------------------------------------------------------------

!end wrf code
!----------------------------------------------------------------------C
!    USE
!        ROUTINE IS CALLED FROM CCPP (AFTER CALLING PBL SCHEMES)
!
!    PURPOSE
!        USING THE GWD PARAMETERIZATIONS OF PS-GLAS AND PH-
!        GFDL TECHNIQUE.  THE TIME TENDENCIES OF U V
!        ARE ALTERED TO INCLUDE THE EFFECT OF MOUNTAIN INDUCED
!        GRAVITY WAVE DRAG FROM SUB-GRID SCALE OROGRAPHY INCLUDING
!        CONVECTIVE BREAKING, SHEAR BREAKING AND THE PRESENCE OF
!        CRITICAL LEVELS
!
!
!   ********************************************************************
   USE MACHINE , ONLY : kind_phys
   implicit none

   ! Interface variables
   integer, intent(in) :: im, km, imx, kdt, ipr, me, master
   integer, intent(in) :: gwd_opt
   logical, intent(in) :: lprnt
   integer, intent(in) :: KPBL(im)
   real(kind=kind_phys), intent(in) :: deltim, G, CP, RD, RV, cdmbgwd(2)

   integer              ::  kpblmax
   integer, parameter   ::  ims=1, kms=1, its=1, kts=1
   real(kind=kind_phys), intent(in) ::  fv, pi
   real(kind=kind_phys) ::  rcl, cdmb
   real(kind=kind_phys) ::  g_inv

   real(kind=kind_phys), intent(inout) ::                        &
     &                   dudt(im,km),dvdt(im,km),                &
     &                   dtdt(im,km)
   real(kind=kind_phys), intent(out) :: rdxzb(im)
   real(kind=kind_phys), intent(in) ::                           &
     &                            u1(im,km),v1(im,km),           &
     &                            t1(im,km),q1(im,km),           &
     &                            PHII(im,km+1),prsl(im,km),     &
     &                            prslk(im,km),PHIL(im,km)
   real(kind=kind_phys), intent(in) ::  prsi(im,km+1),           &
     &                                  del(im,km)
   real(kind=kind_phys), intent(in) ::   var(im),oc1(im),        &
     &                                   oa4(im,4),ol4(im,4),    &
     &                                   dx(im)
   real(kind=kind_phys), intent(in) ::   varss(im),oc1ss(im),    &
     &                              oa4ss(im,4),ol4ss(im,4)
   real(kind=kind_phys), intent(in) :: THETA(im),SIGMA(im),      &
     &                                 GAMMA(im),ELVMAX(im)

! added for small-scale orographic wave drag
   real(kind=kind_phys), dimension(im,km) :: utendwave,vtendwave,thx,thvx
   real(kind=kind_phys), intent(in) ::     br1(im),              &
     &                                     hpbl(im),             &
     &                                     slmsk(im)
   real(kind=kind_phys), dimension(im)    :: govrth,xland
   !real(kind=kind_phys), dimension(im,km) :: dz2
   real(kind=kind_phys)                   :: tauwavex0,tauwavey0,  &
     &                                     XNBV,density,tvcon,hpbl2
   integer                          ::     kpbl2,kvar
   !real(kind=kind_phys), dimension(im,km+1)         ::     zq      ! = PHII/g
   real(kind=kind_phys), dimension(im,km)           ::     zl      ! = PHIL/g

!SPP
   real(kind=kind_phys), dimension(im)              :: rstoch

!Output:
   real(kind=kind_phys), intent(out) ::                          &
     &                      dusfc(im),   dvsfc(im)
!Output (optional):
   real(kind=kind_phys), intent(out) ::                          &
     &                      dusfc_ls(:),dvsfc_ls(:),             &
     &                      dusfc_bl(:),dvsfc_bl(:),             &
     &                      dusfc_ss(:),dvsfc_ss(:),             &
     &                      dusfc_fd(:),dvsfc_fd(:)
   real(kind=kind_phys), intent(out) ::                          &
     &         dtaux2d_ls(:,:),dtauy2d_ls(:,:),                  &
     &         dtaux2d_bl(:,:),dtauy2d_bl(:,:),                  &
     &         dtaux2d_ss(:,:),dtauy2d_ss(:,:),                  &
     &         dtaux2d_fd(:,:),dtauy2d_fd(:,:)

!Misc arrays
   real(kind=kind_phys), dimension(im,km)     :: dtaux2d, dtauy2d

!-------------------------------------------------------------------------
! Flags to regulate the activation of specific components of drag suite:
! Each component is tapered off automatically as a function of dx, so best to
! keep them activated (.true.).
      logical, intent(in) ::   &
      do_gsl_drag_ls_bl,       & ! large-scale gravity wave drag and blocking
      do_gsl_drag_ss,          & ! small-scale gravity wave drag (Steeneveld et al. 2008)
      do_gsl_drag_tofd           ! form drag (Beljaars et al. 2004, QJRMS)

! Additional flags
      integer, parameter ::    &
      gwd_opt_ls      = 1,     & ! large-scale gravity wave drag
      gwd_opt_bl      = 1,     & ! blocking drag
      gsd_diss_ht_opt = 0

! Parameters for bounding the scale-adaptive variability:
! Small-scale GWD + turbulent form drag
   real(kind=kind_phys), parameter      :: dxmin_ss = 1000.,                     &
     &                     dxmax_ss = 12000.  ! min,max range of tapering (m)
! Large-scale GWD + blocking
   real(kind=kind_phys), parameter      :: dxmin_ls = 3000.,                     &
     &                     dxmax_ls = 13000.  ! min,max range of tapering (m)
   real(kind=kind_phys)                 :: ss_taper, ls_taper ! small- and large-scale tapering factors (-)
!
! Variables for limiting topographic standard deviation (var)
   real(kind=kind_phys), parameter      :: varmax_ss = 50.,                      &
                           varmax_fd = 150.,                     &
                           beta_ss = 0.1,                        &
                           beta_fd = 0.2
   real(kind=kind_phys)                 :: var_temp, var_temp2

! added Beljaars orographic form drag
   real(kind=kind_phys), dimension(im,km) :: utendform,vtendform
   real(kind=kind_phys)                 :: a1,a2,wsp
   real(kind=kind_phys)                 :: H_efold

! critical richardson number for wave breaking : ! larger drag with larger value
   real(kind=kind_phys), parameter       ::  ric     = 0.25
   real(kind=kind_phys), parameter       ::  dw2min  = 1.
   real(kind=kind_phys), parameter       ::  rimin   = -100.
   real(kind=kind_phys), parameter       ::  bnv2min = 1.0e-5
   real(kind=kind_phys), parameter       ::  efmin   = 0.0
   real(kind=kind_phys), parameter       ::  efmax   = 10.0
   real(kind=kind_phys), parameter       ::  xl      = 4.0e4
   real(kind=kind_phys), parameter       ::  critac  = 1.0e-5
   real(kind=kind_phys), parameter       ::  gmax    = 1.
   real(kind=kind_phys), parameter       ::  veleps  = 1.0
   real(kind=kind_phys), parameter       ::  factop  = 0.5
   real(kind=kind_phys), parameter       ::  frc     = 1.0
   real(kind=kind_phys), parameter       ::  ce      = 0.8
   real(kind=kind_phys), parameter       ::  cg      = 0.5
   integer,parameter    ::  kpblmin = 2

!
!  local variables
!
   integer              ::  i,j,k,lcap,lcapp1,nwd,idir,           &
                            klcap,kp1
!
   real(kind=kind_phys) ::  rcs,csg,fdir,cleff,cleff_ss,cs,       &
                            rcsks,wdir,ti,rdz,tem2,dw2,shr2,      &
                            bvf2,rdelks,wtkbj,tem,gfobnv,hd,fro,  &
                            rim,temc,tem1,efact,temv,dtaux,dtauy, &
                            dtauxb,dtauyb,eng0,eng1
!
   logical              ::  ldrag(im),icrilv(im),                 &
                            flag(im),kloop1(im)
!
   real(kind=kind_phys) ::  taub(im),taup(im,km+1),               &
                            xn(im),yn(im),                        &
                            ubar(im),vbar(im),                    &
                            fr(im),ulow(im),                      &
                            rulow(im),bnv(im),                    &
                            oa(im),ol(im),                        &
                            oass(im),olss(im),                    &
                            roll(im),dtfac(im),                   &
                            brvf(im),xlinv(im),                   &
                            delks(im),delks1(im),                 &
                            bnv2(im,km),usqj(im,km),              &
                            taud_ls(im,km),taud_bl(im,km),        &
                            ro(im,km),                            &
                            vtk(im,km),vtj(im,km),                &
                            zlowtop(im),velco(im,km-1),           &
                            coefm(im),coefm_ss(im)
!
   integer              ::  kbl(im),klowtop(im)
   integer,parameter    ::  mdir=8
   !integer              ::  nwdir(mdir)
   !data nwdir/6,7,5,8,2,3,1,4/
   integer, parameter :: nwdir(8) = (/6,7,5,8,2,3,1,4/)
!
!  variables for flow-blocking drag
!
   real(kind=kind_phys),parameter       :: frmax  = 10.
   real(kind=kind_phys),parameter       :: olmin  = 1.0e-5
   real(kind=kind_phys),parameter       :: odmin  = 0.1
   real(kind=kind_phys),parameter       :: odmax  = 10.
   real(kind=kind_phys),parameter       :: erad   = 6371.315e+3
   integer              :: komax(im)
   integer              :: kblk
   real(kind=kind_phys)                 :: cd
   real(kind=kind_phys)                 :: zblk,tautem
   real(kind=kind_phys)                 :: pe,ke
   real(kind=kind_phys)                 :: delx,dely,dxy4(4),dxy4p(4)
   real(kind=kind_phys)                 :: dxy(im),dxyp(im)
   real(kind=kind_phys)                 :: ol4p(4),olp(im),od(im)
   real(kind=kind_phys)                 :: taufb(im,km+1)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! Calculate inverse of gravitational acceleration
   g_inv = 1./G

   ! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0
   var_temp2 = 0.


!--------------------------------------------------------------------
! SCALE-ADPTIVE PARAMETER FROM GFS GWD SCHEME
!--------------------------------------------------------------------
!     parameter (cdmb = 1.0)     ! non-dim sub grid mtn drag Amp (*j*)
! non-dim sub grid mtn drag Amp (*j*)
!     cdmb = 1.0/float(IMX/192)
!     cdmb = 192.0/float(IMX)
      cdmb = 4.0 * 192.0/float(IMX)
      if (cdmbgwd(1) >= 0.0) cdmb = cdmb * cdmbgwd(1)

!>-# Orographic Gravity Wave Drag Section
      kpblmax = km / 2 ! maximum pbl height : # of vertical levels / 2
!
!  Scale cleff between IM=384*2 and 192*2 for T126/T170 and T62
!
      if (imx > 0) then
!       cleff = 1.0E-5 * SQRT(FLOAT(IMX)/384.0)
!       cleff = 1.0E-5 * SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
!       cleff = 0.5E-5 * SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
!       cleff = 1.0E-5 * SQRT(FLOAT(IMX)/192)/float(IMX/192)
!       cleff = 1.0E-5 / SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
        cleff = 0.5E-5 / SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
! hmhj for ndsl
! jw        cleff = 0.1E-5 / SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
!       cleff = 2.0E-5 * SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
!       cleff = 2.5E-5 * SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
      endif
      if (cdmbgwd(2) >= 0.0) cleff = cleff * cdmbgwd(2)
!--------------------------------------------------------------------
! END SCALE-ADPTIVE PARAMETER SECTION
!--------------------------------------------------------------------
!
!---- constants
!
   rcl    = 1.
   rcs    = sqrt(rcl)
   cs     = 1. / sqrt(rcl)
   csg    = cs * g
   lcap   = km
   lcapp1 = lcap + 1
   fdir   = mdir / (2.0*pi)

   do i=1,im
      if (slmsk(i)==1. .or. slmsk(i)==2.) then !sea/land/ice mask (=0/1/2) in FV3
         xland(i)=1.0                          !but land/water = (1/2) in this module
      else
         xland(i)=2.0
      endif
      RDXZB(i)  = 0.0
   enddo

!temporary use of large-scale data:
!   do i=1,im
!      varss(i)=var(i)
!      oc1ss(i)=oc1(i)
!      do j=1,4
!         oa4ss(i,j)=oa4(i,j)
!         ol4ss(i,j)=ol4(i,j)
!      enddo
!   enddo
!
!--- calculate scale-aware tapering factors
!NOTE: if dx(1) is not representative of most/all dx, this needs to change...
if ( dx(1) .ge. dxmax_ls ) then
   ls_taper = 1.
else
   if ( dx(1) .le. dxmin_ls) then
      ls_taper = 0.
   else
      ls_taper = 0.5 * ( SIN(pi*(dx(1)-0.5*(dxmax_ls+dxmin_ls))/    &
                                (dxmax_ls-dxmin_ls)) + 1. )
   end if
end if
! if (me==master) print *,"in Drag Suite, dx(1:2):",dx(1),dx(2)
if ( dx(1) .ge. dxmax_ss ) then
   ss_taper = 1.
else
   if ( dx(1) .le. dxmin_ss) then
      ss_taper = 0.
   else
      ss_taper = dxmax_ss * (1. - dxmin_ss/dx(1))/(dxmax_ss-dxmin_ss)
   end if
end if
! if (me==master) print *,"in Drag Suite, ss_taper:",ss_taper

!--- calculate length of grid for flow-blocking drag
!
   delx   = dx(1)
   dely   = dx(1)
   dxy4(1)  = delx
   dxy4(2)  = dely
   dxy4(3)  = sqrt(delx*delx + dely*dely)
   dxy4(4)  = dxy4(3)
   dxy4p(1) = dxy4(2)
   dxy4p(2) = dxy4(1)
   dxy4p(3) = dxy4(4)
   dxy4p(4) = dxy4(3)
!
!-----initialize arrays
!
   dtaux = 0.0
   dtauy = 0.0
   do i = its,im
     klowtop(i)    = 0
     kbl(i)        = 0
   enddo
!
   do i = its,im
     xn(i)         = 0.0
     yn(i)         = 0.0
     ubar (i)      = 0.0
     vbar (i)      = 0.0
     roll (i)      = 0.0
     taub (i)      = 0.0
     oa(i)         = 0.0
     ol(i)         = 0.0
     oass(i)       = 0.0
     olss(i)       = 0.0
     ulow (i)      = 0.0
     dtfac(i)      = 1.0
     rstoch(i)     = 0.0
     ldrag(i)      = .false.
     icrilv(i)     = .false.
     flag(i)       = .true.
   enddo

   do k = kts,km
     do i = its,im
       usqj(i,k) = 0.0
       bnv2(i,k) = 0.0
       vtj(i,k)  = 0.0
       vtk(i,k)  = 0.0
       taup(i,k) = 0.0
       taud_ls(i,k) = 0.0
       taud_bl(i,k) = 0.0
       dtaux2d(i,k) = 0.0
       dtauy2d(i,k) = 0.0
     enddo
   enddo
!
   if ( (gwd_opt == 33).or.(gwd_opt == 22) ) then
     do i = its,im
       dusfc_ls(i) = 0.0
       dvsfc_ls(i) = 0.0
       dusfc_bl(i) = 0.0
       dvsfc_bl(i) = 0.0
       dusfc_ss(i) = 0.0
       dvsfc_ss(i) = 0.0
       dusfc_fd(i) = 0.0
       dvsfc_fd(i) = 0.0
     enddo
     do k = kts,km
       do i = its,im
         dtaux2d_ls(i,k)= 0.0
         dtauy2d_ls(i,k)= 0.0
         dtaux2d_bl(i,k)= 0.0
         dtauy2d_bl(i,k)= 0.0
         dtaux2d_ss(i,k)= 0.0
         dtauy2d_ss(i,k)= 0.0
         dtaux2d_fd(i,k)= 0.0
         dtauy2d_fd(i,k)= 0.0
       enddo
     enddo
   endif

   do i = its,im
     taup(i,km+1) = 0.0
     xlinv(i)     = 1.0/xl
     dusfc(i) = 0.0
     dvsfc(i) = 0.0
   enddo
!
!  initialize array for flow-blocking drag
!
   taufb(1:im,1:km+1) = 0.0
   komax(1:im) = 0
!
   do k = kts,km
     do i = its,im
       vtj(i,k)  = t1(i,k)  * (1.+fv*q1(i,k))
       vtk(i,k)  = vtj(i,k) / prslk(i,k)
       ro(i,k)   = 1./rd * prsl(i,k) / vtj(i,k) ! density kg/m**3
     enddo
   enddo
!
!  calculate mid-layer height (zl), interface height (zq), and layer depth (dz2).
!
   !zq=0.
   do k = kts,km
     do i = its,im
       !zq(i,k+1) = PHII(i,k+1)*g_inv
       !dz2(i,k)  = (PHII(i,k+1)-PHII(i,k))*g_inv
       zl(i,k)   = PHIL(i,k)*g_inv
     enddo
   enddo
!
!  determine reference level: maximum of 2*var and pbl heights
!
   do i = its,im
     zlowtop(i) = 2. * var(i)
   enddo
!
   do i = its,im
     kloop1(i) = .true.
   enddo
!
   do k = kts+1,km
     do i = its,im
       if(kloop1(i).and.zl(i,k)-zl(i,1).ge.zlowtop(i)) then
         klowtop(i) = k+1
         kloop1(i)  = .false.
       endif
     enddo
   enddo
!
   do i = its,im
     kbl(i)   = max(kpbl(i), klowtop(i))
     kbl(i)   = max(min(kbl(i),kpblmax),kpblmin)
   enddo
!
!  determine the level of maximum orographic height
!
   ! komax(:) = kbl(:)
   komax(:) = klowtop(:) - 1    ! modification by NOAA/GSD March 2018
!
   do i = its,im
     delks(i)  = 1.0 / (prsi(i,1) - prsi(i,kbl(i)))
     delks1(i) = 1.0 / (prsl(i,1) - prsl(i,kbl(i)))
   enddo
!
!  compute low level averages within pbl
!
   do k = kts,kpblmax
     do i = its,im
       if (k.lt.kbl(i)) then
         rcsks   = rcs     * del(i,k) * delks(i)
         rdelks  = del(i,k)  * delks(i)
         ubar(i) = ubar(i) + rcsks  * u1(i,k)      ! pbl u  mean
         vbar(i) = vbar(i) + rcsks  * v1(i,k)      ! pbl v  mean
         roll(i) = roll(i) + rdelks * ro(i,k)      ! ro mean
       endif
     enddo
   enddo
!
!     figure out low-level horizontal wind direction
!
!             nwd  1   2   3   4   5   6   7   8
!              wd  w   s  sw  nw   e   n  ne  se
!
   do i = its,im
     wdir   = atan2(ubar(i),vbar(i)) + pi
     idir   = mod(nint(fdir*wdir),mdir) + 1
     nwd    = nwdir(idir)
     oa(i)  = (1-2*int( (nwd-1)/4 )) * oa4(i,mod(nwd-1,4)+1)
     ol(i) = ol4(i,mod(nwd-1,4)+1)
     ! Repeat for small-scale gwd
     oass(i)  = (1-2*int( (nwd-1)/4 )) * oa4ss(i,mod(nwd-1,4)+1)
     olss(i) = ol4ss(i,mod(nwd-1,4)+1)

!
!----- compute orographic width along (ol) and perpendicular (olp)
!----- the direction of wind
!
     ol4p(1) = ol4(i,2)
     ol4p(2) = ol4(i,1)
     ol4p(3) = ol4(i,4)
     ol4p(4) = ol4(i,3)
     olp(i)  = ol4p(mod(nwd-1,4)+1)
!
!----- compute orographic direction (horizontal orographic aspect ratio)
!
     od(i) = olp(i)/max(ol(i),olmin)
     od(i) = min(od(i),odmax)
     od(i) = max(od(i),odmin)
!
!----- compute length of grid in the along(dxy) and cross(dxyp) wind directions
!
     dxy(i)  = dxy4(MOD(nwd-1,4)+1)
     dxyp(i) = dxy4p(MOD(nwd-1,4)+1)
   enddo
!
! END INITIALIZATION; BEGIN GWD CALCULATIONS:
!
IF ( (do_gsl_drag_ls_bl).and.                            &
     ((gwd_opt_ls .EQ. 1).or.(gwd_opt_bl .EQ. 1)).and.   &
               (ls_taper .GT. 1.E-02) ) THEN   !====
!
!---  saving richardson number in usqj for migwdi
!
   do k = kts,km-1
     do i = its,im
       ti        = 2.0 / (t1(i,k)+t1(i,k+1))
       rdz       = 1./(zl(i,k+1) - zl(i,k))
       tem1      = u1(i,k) - u1(i,k+1)
       tem2      = v1(i,k) - v1(i,k+1)
       dw2       = rcl*(tem1*tem1 + tem2*tem2)
       shr2      = max(dw2,dw2min) * rdz * rdz
       bvf2      = g*(g/cp+rdz*(vtj(i,k+1)-vtj(i,k))) * ti
       usqj(i,k) = max(bvf2/shr2,rimin)
       bnv2(i,k) = 2.0*g*rdz*(vtk(i,k+1)-vtk(i,k))/(vtk(i,k+1)+vtk(i,k))
       bnv2(i,k) = max( bnv2(i,k), bnv2min )
     enddo
   enddo
!
!----compute the "low level" or 1/3 wind magnitude (m/s)
!
   do i = its,im
     ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0)
     rulow(i) = 1./ulow(i)
   enddo
!
   do k = kts,km-1
     do i = its,im
       velco(i,k)  = (0.5*rcs) * ((u1(i,k)+u1(i,k+1)) * ubar(i)                &
                                + (v1(i,k)+v1(i,k+1)) * vbar(i))
       velco(i,k)  = velco(i,k) * rulow(i)
       if ((velco(i,k).lt.veleps) .and. (velco(i,k).gt.0.)) then
         velco(i,k) = veleps
       endif
     enddo
   enddo
!
!  no drag when critical level in the base layer
!
   do i = its,im
     ldrag(i) = velco(i,1).le.0.
   enddo
!
!  no drag when velco.lt.0
!
   do k = kpblmin,kpblmax
     do i = its,im
       if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. velco(i,k).le.0.
     enddo
   enddo
!
!  no drag when bnv2.lt.0
!
   do k = kts,kpblmax
     do i = its,im
       if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. bnv2(i,k).lt.0.
     enddo
   enddo
!
!-----the low level weighted average ri is stored in usqj(1,1; im)
!-----the low level weighted average n**2 is stored in bnv2(1,1; im)
!---- this is called bnvl2 in phys_gwd_alpert_sub not bnv2
!---- rdelks (del(k)/delks) vert ave factor so we can * instead of /
!
   do i = its,im
     wtkbj     = (prsl(i,1)-prsl(i,2)) * delks1(i)
     bnv2(i,1) = wtkbj * bnv2(i,1)
     usqj(i,1) = wtkbj * usqj(i,1)
   enddo
!
   do k = kpblmin,kpblmax
     do i = its,im
       if (k .lt. kbl(i)) then
         rdelks    = (prsl(i,k)-prsl(i,k+1)) * delks1(i)
         bnv2(i,1) = bnv2(i,1) + bnv2(i,k) * rdelks
         usqj(i,1) = usqj(i,1) + usqj(i,k) * rdelks
       endif
     enddo
   enddo
!
   do i = its,im
     ldrag(i) = ldrag(i) .or. bnv2(i,1).le.0.0
     ldrag(i) = ldrag(i) .or. ulow(i).eq.1.0
     ldrag(i) = ldrag(i) .or. var(i) .le. 0.0
   enddo
!
!  set all ri low level values to the low level value
!
   do k = kpblmin,kpblmax
     do i = its,im
       if (k .lt. kbl(i)) usqj(i,k) = usqj(i,1)
     enddo
   enddo
!
   do i = its,im
     if (.not.ldrag(i))   then
       bnv(i) = sqrt( bnv2(i,1) )
       fr(i) = bnv(i)  * rulow(i) * 2. * var(i) * od(i)
       fr(i) = min(fr(i),frmax)
       xn(i)  = ubar(i) * rulow(i)
       yn(i)  = vbar(i) * rulow(i)
     endif
   enddo
!
!  compute the base level stress and store it in taub
!  calculate enhancement factor, number of mountains & aspect
!  ratio const. use simplified relationship between standard
!  deviation & critical hgt

   do i = its,im
     if (.not. ldrag(i))   then
       efact    = (oa(i) + 2.) ** (ce*fr(i)/frc)
       efact    = min( max(efact,efmin), efmax )
!!!!!!! cleff (effective grid length) is highly tunable parameter
!!!!!!! the bigger (smaller) value produce weaker (stronger) wave drag
!WRF       cleff    = sqrt(dxy(i)**2. + dxyp(i)**2.)
!WRF       cleff    = 3. * max(dx(i),cleff)
       coefm(i) = (1. + ol(i)) ** (oa(i)+1.)
!WRF       xlinv(i) = coefm(i) / cleff
       xlinv(i) = coefm(i) * cleff
       tem      = fr(i) * fr(i) * oc1(i)
       gfobnv   = gmax * tem / ((tem + cg)*bnv(i))
       if ( gwd_opt_ls .NE. 0 ) then
          taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)                       &
                   * ulow(i) * gfobnv * efact
       else     ! We've gotten what we need for the blocking scheme
          taub(i) = 0.0
       end if
     else
       taub(i) = 0.0
       xn(i)   = 0.0
       yn(i)   = 0.0
     endif
   enddo

ENDIF   ! (do_gsl_drag_ls_bl).and.((gwd_opt_ls .EQ. 1).or.(gwd_opt_bl .EQ. 1))

!=========================================================
! add small-scale wavedrag for stable boundary layer
!=========================================================
  XNBV=0.
  tauwavex0=0.
  tauwavey0=0.
  density=1.2
  utendwave=0.
  vtendwave=0.
!
  IF ( (do_gsl_drag_ss).and.(ss_taper.GT.1.E-02) ) THEN
    ! if (me==master) print *,"in Drag Suite: Running small-scale gravity wave drag"
!
! declaring potential temperature
!
    do k = kts,km
      do i = its,im
        thx(i,k) = t1(i,k)/prslk(i,k)
      enddo
    enddo
!
    do k = kts,km
      do i = its,im
        tvcon = (1.+fv*q1(i,k))
        thvx(i,k) = thx(i,k)*tvcon
      enddo
    enddo

    do i=its,im
       hpbl2 = hpbl(i)+10.
       kpbl2 = kpbl(i)
       !kvar = MIN(kpbl, k-level of var)
       kvar = 1
       do k=kts+1,MAX(kpbl(i),kts+1)
!          IF (zl(i,k)>2.*var(i) .or. zl(i,k)>2*varmax) then
          IF (zl(i,k)>300.) then
             kpbl2 = k
             IF (k == kpbl(i)) then
                hpbl2 = hpbl(i)+10.
             ELSE
                hpbl2 = zl(i,k)+10.
             ENDIF
             exit
          ENDIF
       enddo
       if((xland(i)-1.5).le.0. .and. 2.*varss(i).le.hpbl(i))then
          if(br1(i).gt.0. .and. thvx(i,kpbl2)-thvx(i,kts) > 0.)then
            cleff_ss    = sqrt(dxy(i)**2 + dxyp(i)**2)   ! WRF
!            cleff_ss    = 3. * max(dx(i),cleff_ss)
!            cleff_ss    = 10. * max(dxmax_ss,cleff_ss)
            cleff_ss    = 0.1 * max(dxmax_ss,cleff_ss)  ! WRF
!             cleff_ss    = 0.1 * 12000.
            coefm_ss(i) = (1. + olss(i)) ** (oass(i)+1.)
            xlinv(i) = coefm_ss(i) / cleff_ss
            !govrth(i)=g/(0.5*(thvx(i,kpbl(i))+thvx(i,kts)))
            govrth(i)=g/(0.5*(thvx(i,kpbl2)+thvx(i,kts)))
            !XNBV=sqrt(govrth(i)*(thvx(i,kpbl(i))-thvx(i,kts))/hpbl(i))
            XNBV=sqrt(govrth(i)*(thvx(i,kpbl2)-thvx(i,kts))/hpbl2)
!
            !if(abs(XNBV/u1(i,kpbl(i))).gt.xlinv(i))then
            if(abs(XNBV/u1(i,kpbl2)).gt.xlinv(i))then
              !tauwavex0=0.5*XNBV*xlinv(i)*(2*MIN(varss(i),75.))**2*ro(i,kts)*u1(i,kpbl(i))
              !tauwavex0=0.5*XNBV*xlinv(i)*(2.*MIN(varss(i),40.))**2*ro(i,kts)*u1(i,kpbl2)
              !tauwavex0=0.5*XNBV*xlinv(i)*(2.*MIN(varss(i),40.))**2*ro(i,kts)*u1(i,3)
              var_temp = MIN(varss(i),varmax_ss) +                                   &
                            MAX(0.,beta_ss*(varss(i)-varmax_ss))
              ! Note:  This is a semi-implicit treatment of the time differencing
              var_temp2 = 0.5*XNBV*xlinv(i)*(2.*var_temp)**2*ro(i,kvar)  ! this is greater than zero
              tauwavex0=-var_temp2*u1(i,kvar)/(1.+var_temp2*deltim)
              tauwavex0=tauwavex0*ss_taper
            else
              tauwavex0=0.
            endif
!
            !if(abs(XNBV/v1(i,kpbl(i))).gt.xlinv(i))then
            if(abs(XNBV/v1(i,kpbl2)).gt.xlinv(i))then
              !tauwavey0=0.5*XNBV*xlinv(i)*(2*MIN(varss(i),75.))**2*ro(i,kts)*v1(i,kpbl(i))
              !tauwavey0=0.5*XNBV*xlinv(i)*(2.*MIN(varss(i),40.))**2*ro(i,kts)*v1(i,kpbl2)
              !tauwavey0=0.5*XNBV*xlinv(i)*(2.*MIN(varss(i),40.))**2*ro(i,kts)*v1(i,3)
              var_temp = MIN(varss(i),varmax_ss) +                                   &
                            MAX(0.,beta_ss*(varss(i)-varmax_ss))
              ! Note:  This is a semi-implicit treatment of the time differencing
              tauwavey0=-var_temp2*v1(i,kvar)/(1.+var_temp2*deltim)
              tauwavey0=tauwavey0*ss_taper
            else
              tauwavey0=0.
            endif

            do k=kts,kpbl(i) !MIN(kpbl2+1,km-1)
!original
              !utendwave(i,k)=-1.*tauwavex0*2.*max((1.-zl(i,k)/hpbl(i)),0.)/hpbl(i)
              !vtendwave(i,k)=-1.*tauwavey0*2.*max((1.-zl(i,k)/hpbl(i)),0.)/hpbl(i)
!new
              utendwave(i,k)=-1.*tauwavex0*2.*max((1.-zl(i,k)/hpbl2),0.)/hpbl2
              vtendwave(i,k)=-1.*tauwavey0*2.*max((1.-zl(i,k)/hpbl2),0.)/hpbl2
!mod-to be used in HRRRv3/RAPv4
              !utendwave(i,k)=-1.*tauwavex0 * max((1.-zl(i,k)/hpbl2),0.)**2
              !vtendwave(i,k)=-1.*tauwavey0 * max((1.-zl(i,k)/hpbl2),0.)**2
            enddo
          endif
       endif
    enddo ! end i loop

    do k = kts,km
       do i = its,im
         dudt(i,k)  = dudt(i,k) + utendwave(i,k)
         dvdt(i,k)  = dvdt(i,k) + vtendwave(i,k)
         dusfc(i)   = dusfc(i) + utendwave(i,k) * del(i,k)
         dvsfc(i)   = dvsfc(i) + vtendwave(i,k) * del(i,k)
       enddo
    enddo
    if ( (gwd_opt == 33).or.(gwd_opt == 22) ) then
      do k = kts,km
        do i = its,im
          dusfc_ss(i) = dusfc_ss(i) + utendwave(i,k) * del(i,k)
          dvsfc_ss(i) = dvsfc_ss(i) + vtendwave(i,k) * del(i,k)
          dtaux2d_ss(i,k) = utendwave(i,k)
          dtauy2d_ss(i,k) = vtendwave(i,k)
        enddo
      enddo
    endif

ENDIF  ! if (do_gsl_drag_ss)

!================================================================
! Topographic Form Drag from Beljaars et al. (2004, QJRMS, equ. 16):
!================================================================
IF ( (do_gsl_drag_tofd).and.(ss_taper.GT.1.E-02) ) THEN
    ! if (me==master) print *,"in Drag Suite: Running form drag"

   utendform=0.
   vtendform=0.

   DO i=its,im
      IF ((xland(i)-1.5) .le. 0.) then
         !(IH*kflt**n1)**-1 = (0.00102*0.00035**-1.9)**-1 = 0.00026615161
          var_temp = MIN(varss(i),varmax_fd) +                            &
                     MAX(0.,beta_fd*(varss(i)-varmax_fd))
          var_temp = MIN(var_temp, 250.)
          a1=0.00026615161*var_temp**2
!         a1=0.00026615161*MIN(varss(i),varmax)**2
!         a1=0.00026615161*(0.5*varss(i))**2
         ! k1**(n1-n2) = 0.003**(-1.9 - -2.8) = 0.003**0.9 = 0.005363
         a2=a1*0.005363
         ! Revise e-folding height based on PBL height and topographic std. dev. -- M. Toy 3/12/2018
         H_efold = max(2*varss(i),hpbl(i))
         H_efold = min(H_efold,1500.)
         DO k=kts,km
            wsp=SQRT(u1(i,k)**2 + v1(i,k)**2)
            ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759
            var_temp = 0.0759*EXP(-(zl(i,k)/H_efold)**1.5)*a2*       &
                              zl(i,k)**(-1.2)*ss_taper   ! this is greater than zero
            !  Note:  This is a semi-implicit treatment of the time differencing
            !  per Beljaars et al. (2004, QJRMS)
            utendform(i,k) = - var_temp*wsp*u1(i,k)/(1. + var_temp*deltim*wsp)
            vtendform(i,k) = - var_temp*wsp*v1(i,k)/(1. + var_temp*deltim*wsp)
            !IF(zl(i,k) > 4000.) exit
         ENDDO
      ENDIF
   ENDDO

   do k = kts,km
      do i = its,im
         dudt(i,k)  = dudt(i,k) + utendform(i,k)
         dvdt(i,k)  = dvdt(i,k) + vtendform(i,k)
         dusfc(i)   = dusfc(i) + utendform(i,k) * del(i,k)
         dvsfc(i)   = dvsfc(i) + vtendform(i,k) * del(i,k)
      enddo
   enddo
   if ( (gwd_opt == 33).or.(gwd_opt == 22) ) then
     do k = kts,km
       do i = its,im
         dtaux2d_fd(i,k) = utendform(i,k)
         dtauy2d_fd(i,k) = vtendform(i,k)
         dusfc_fd(i) = dusfc_fd(i) + utendform(i,k) * del(i,k)
         dvsfc_fd(i) = dvsfc_fd(i) + vtendform(i,k) * del(i,k)
       enddo
     enddo
   endif

ENDIF  ! if (do_gsl_drag_tofd)
!=======================================================
! More for the large-scale gwd component
IF ( (do_gsl_drag_ls_bl).and.                                        &
     (gwd_opt_ls .EQ. 1).and.(ls_taper.GT.1.E-02) ) THEN
    ! if (me==master) print *,"in Drag Suite: Running large-scale gravity wave drag"
!
!   now compute vertical structure of the stress.
   do k = kts,kpblmax
      do i = its,im
         if (k .le. kbl(i)) taup(i,k) = taub(i)
      enddo
   enddo
!
   do k = kpblmin, km-1                   ! vertical level k loop!
      kp1 = k + 1
      do i = its,im
!
!   unstablelayer if ri < ric
!   unstable layer if upper air vel comp along surf vel <=0 (crit lay)
!   at (u-c)=0. crit layer exists and bit vector should be set (.le.)
!
         if (k .ge. kbl(i)) then
           icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric)                  &
                                 .or. (velco(i,k) .le. 0.0)
           brvf(i)  = max(bnv2(i,k),bnv2min) ! brunt-vaisala frequency squared
           brvf(i)  = sqrt(brvf(i))          ! brunt-vaisala frequency
         endif
      enddo
!
      do i = its,im
        if (k .ge. kbl(i) .and. (.not. ldrag(i)))   then
          if (.not.icrilv(i) .and. taup(i,k) .gt. 0.0 ) then
            temv = 1.0 / velco(i,k)
            tem1 = coefm(i)/dxy(i)*(ro(i,kp1)+ro(i,k))*brvf(i)*velco(i,k)*0.5
            hd   = sqrt(taup(i,k) / tem1)
            fro  = brvf(i) * hd * temv
!
!  rim is the minimum-richardson number by shutts (1985)
            tem2   = sqrt(usqj(i,k))
            tem    = 1. + tem2 * fro
            rim    = usqj(i,k) * (1.-fro) / (tem * tem)
!
!  check stability to employ the 'saturation hypothesis'
!  of lindzen (1981) except at tropospheric downstream regions
!
            if (rim .le. ric) then  ! saturation hypothesis!
              if ((oa(i) .le. 0.).or.(kp1 .ge. kpblmin )) then
                temc = 2.0 + 1.0 / tem2
                hd   = velco(i,k) * (2.*sqrt(temc)-temc) / brvf(i)
                taup(i,kp1) = tem1 * hd * hd
              endif
            else                    ! no wavebreaking!
              taup(i,kp1) = taup(i,k)
            endif
          endif
        endif
      enddo
   enddo
!
   if(lcap.lt.km) then
      do klcap = lcapp1,km
         do i = its,im
           taup(i,klcap) = prsi(i,klcap) / prsi(i,lcap) * taup(i,lcap)
         enddo
      enddo
   endif

ENDIF !END LARGE-SCALE TAU CALCULATION
!===============================================================
!COMPUTE BLOCKING COMPONENT                                     
!===============================================================
IF ( (do_gsl_drag_ls_bl) .and.                                       &
     (gwd_opt_bl .EQ. 1) .and. (ls_taper .GT. 1.E-02) ) THEN
   ! if (me==master) print *,"in Drag Suite: Running blocking drag"

   do i = its,im
      if(.not.ldrag(i)) then
!
!------- determine the height of flow-blocking layer
!
        kblk = 0
        pe = 0.0
        do k = km, kpblmin, -1
          if(kblk.eq.0 .and. k.le.komax(i)) then
            pe = pe + bnv2(i,k)*(zl(i,komax(i))-zl(i,k))*del(i,k)/g/ro(i,k)
            ke = 0.5*((rcs*u1(i,k))**2.+(rcs*v1(i,k))**2.)
!
!---------- apply flow-blocking drag when pe >= ke
!
            if(pe.ge.ke) then
              kblk = k
              kblk = min(kblk,kbl(i))
              zblk = zl(i,kblk)-zl(i,kts)
              RDXZB(i) = real(k,kind=kind_phys)
            endif
          endif
        enddo
        if(kblk.ne.0) then
!
!--------- compute flow-blocking stress
!
          cd = max(2.0-1.0/od(i),0.0)
          taufb(i,kts) = 0.5 * roll(i) * coefm(i) / max(dxmax_ls,dxy(i))**2 * cd * dxyp(i)   &
                         * olp(i) * zblk * ulow(i)**2
          tautem = taufb(i,kts)/float(kblk-kts)
          do k = kts+1, kblk
            taufb(i,k) = taufb(i,k-1) - tautem
          enddo
!
!----------sum orographic GW stress and flow-blocking stress
!
          ! taup(i,:) = taup(i,:) + taufb(i,:)   ! Keep taup and taufb separate for now
        endif
      endif
   enddo

ENDIF   ! end blocking drag
!===========================================================
IF ( (do_gsl_drag_ls_bl) .and.                                       &
     (gwd_opt_ls .EQ. 1 .OR. gwd_opt_bl .EQ. 1) .and. (ls_taper .GT. 1.E-02) ) THEN
!
!  calculate - (g)*d(tau)/d(pressure) and deceleration terms dtaux, dtauy
!
   do k = kts,km
     do i = its,im
       taud_ls(i,k) = 1. * (taup(i,k+1) - taup(i,k)) * csg / del(i,k)
       taud_bl(i,k) = 1. * (taufb(i,k+1) - taufb(i,k)) * csg / del(i,k)
     enddo
   enddo
!
!  limit de-acceleration (momentum deposition ) at top to 1/2 value
!  the idea is some stuff must go out the 'top'
   do klcap = lcap,km
     do i = its,im
       taud_ls(i,klcap) = taud_ls(i,klcap) * factop
       taud_bl(i,klcap) = taud_bl(i,klcap) * factop
     enddo
   enddo
!
!  if the gravity wave drag would force a critical line
!  in the lower ksmm1 layers during the next deltim timestep,
!  then only apply drag until that critical line is reached.
!
   do k = kts,kpblmax-1
      do i = its,im
         if (k .le. kbl(i)) then
           if((taud_ls(i,k)+taud_bl(i,k)).ne.0.)                         &
              dtfac(i) = min(dtfac(i),abs(velco(i,k)                     &
                   /(deltim*rcs*(taud_ls(i,k)+taud_bl(i,k)))))
         endif
      enddo
   enddo
!
   do k = kts,km
      do i = its,im
         taud_ls(i,k)  = taud_ls(i,k) * dtfac(i) * ls_taper *(1.-rstoch(i))
         taud_bl(i,k)  = taud_bl(i,k) * dtfac(i) * ls_taper *(1.-rstoch(i))

         dtaux  = taud_ls(i,k) * xn(i)
         dtauy  = taud_ls(i,k) * yn(i)
         dtauxb = taud_bl(i,k) * xn(i)
         dtauyb = taud_bl(i,k) * yn(i)

         !add blocking and large-scale contributions to tendencies 
         dudt(i,k)  = dtaux + dtauxb + dudt(i,k)
         dvdt(i,k)  = dtauy + dtauyb + dvdt(i,k)

         if ( gsd_diss_ht_opt .EQ. 1 ) then
            ! Calculate dissipation heating
            ! Initial kinetic energy (at t0-dt)
            eng0 = 0.5*( (rcs*u1(i,k))**2. + (rcs*v1(i,k))**2. )
            ! Kinetic energy after wave-breaking/flow-blocking
            eng1 = 0.5*( (rcs*(u1(i,k)+(dtaux+dtauxb)*deltim))**2 + &
                         (rcs*(v1(i,k)+(dtauy+dtauyb)*deltim))**2 )
            ! Modify theta tendency
            dtdt(i,k) = dtdt(i,k) + max((eng0-eng1),0.0)/cp/deltim
         end if

         dusfc(i)   = dusfc(i) + taud_ls(i,k)*xn(i)*del(i,k) + taud_bl(i,k)*xn(i)*del(i,k)
         dvsfc(i)   = dvsfc(i) + taud_ls(i,k)*yn(i)*del(i,k) + taud_bl(i,k)*yn(i)*del(i,k)
      enddo
   enddo

   !  Finalize dusfc and dvsfc diagnostics
   do i = its,im
    dusfc(i) = (-1./g*rcs) * dusfc(i)
    dvsfc(i) = (-1./g*rcs) * dvsfc(i)
   enddo

   if ( (gwd_opt == 33).or.(gwd_opt == 22) ) then
     do k = kts,km
       do i = its,im
         dtaux2d_ls(i,k) = taud_ls(i,k) * xn(i)
         dtauy2d_ls(i,k) = taud_ls(i,k) * yn(i)
         dtaux2d_bl(i,k) = taud_bl(i,k) * xn(i)
         dtauy2d_bl(i,k) = taud_bl(i,k) * yn(i)
         dusfc_ls(i)  = dusfc_ls(i) + dtaux2d_ls(i,k) * del(i,k)
         dvsfc_ls(i)  = dvsfc_ls(i) + dtauy2d_ls(i,k) * del(i,k)
         dusfc_bl(i)  = dusfc_bl(i) + dtaux2d_bl(i,k) * del(i,k)
         dvsfc_bl(i)  = dvsfc_bl(i) + dtauy2d_bl(i,k) * del(i,k)
       enddo
     enddo
   endif

ENDIF  ! (do_gsl_drag_ls_bl).and.(gwd_opt_ls.EQ.1 .OR. gwd_opt_bl.EQ.1)

if ( (gwd_opt == 33).or.(gwd_opt == 22) ) then
  !  Finalize dusfc and dvsfc diagnostics
  do i = its,im
    dusfc_ls(i) = (-1./g*rcs) * dusfc_ls(i)
    dvsfc_ls(i) = (-1./g*rcs) * dvsfc_ls(i)
    dusfc_bl(i) = (-1./g*rcs) * dusfc_bl(i)
    dvsfc_bl(i) = (-1./g*rcs) * dvsfc_bl(i)
    dusfc_ss(i) = (-1./g*rcs) * dusfc_ss(i)
    dvsfc_ss(i) = (-1./g*rcs) * dvsfc_ss(i)
    dusfc_fd(i) = (-1./g*rcs) * dusfc_fd(i)
    dvsfc_fd(i) = (-1./g*rcs) * dvsfc_fd(i)
  enddo
endif
!
   return
   end subroutine drag_suite_run
!-------------------------------------------------------------------
!

      subroutine drag_suite_finalize()
      end subroutine drag_suite_finalize

      end module drag_suite

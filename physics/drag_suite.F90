!> \file drag_suite.F90
!! This file is the  parameterization of orographic gravity wave
!! drag, mountain blocking, and form drag.

      module drag_suite

      contains

!> \defgroup gfs_drag_suite_mod GSL drag_suite Module
!> This module contains the CCPP-compliant GSL orographic gravity wave drag scheme.
!> @{
!!
!> \brief This subroutine initializes the orographic gravity wave drag scheme.
!!
!> \section arg_table_drag_suite_init Argument Table
!! \htmlinclude drag_suite_init.html
!!
      subroutine drag_suite_init(gwd_opt, errmsg, errflg)

      integer,          intent(in)  :: gwd_opt
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg


      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Consistency checks
      if (gwd_opt/=3 .and. gwd_opt/=33) then
        write(errmsg,'(*(a))') "Logic error: namelist choice of gravity wave &
          & drag is different from drag_suite scheme"
        errflg = 1
        return
      end if        
      end subroutine drag_suite_init

!> \brief This subroutine includes orographic gravity wave drag, mountain
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
     &           dtaux2d_ms,dtauy2d_ms,dtaux2d_bl,dtauy2d_bl,           &
     &           dtaux2d_ss,dtauy2d_ss,dtaux2d_fd,dtauy2d_fd,           &
     &           dusfc,dvsfc,                                           &
     &           dusfc_ms,dvsfc_ms,dusfc_bl,dvsfc_bl,                   &
     &           dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd,                   &
     &           slmsk,br1,hpbl,                                        &
     &           g, cp, rd, rv, fv, pi, imx, cdmbgwd, me, master,       &
     &           lprnt, ipr, rdxzb, dx, gwd_opt,                        &
     &           do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd,   &
     &           dtend, dtidx, index_of_process_orographic_gwd,         &
     &           index_of_temperature, index_of_x_wind,                 &
     &           index_of_y_wind, ldiag3d, ldiag_ugwp, ugwp_seq_update, & 
     &           spp_wts_gwd, spp_gwd, errmsg, errflg)

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
!                    gwd_opt_ms = 0 or 1: mesoscale    (changed to logical flag)
!                    gwd_opt_bl = 0 or 1: blocking drag  (changed to logical flag)
!                    gwd_opt_ss = 0 or 1: small-scale gravity wave drag   (removed)
!                    gwd_opt_fd = 0 or 1: topographic form drag      (removed)
!    2017-09-25  Michael Toy (from NCEP GFS model) added dissipation heating (logical flags)
!                    gsd_diss_ht_opt = .false. : dissipation heating off
!                    gsd_diss_ht_opt = .true.  : dissipation heating on
!    2020-08-25  Michael Toy changed logic control for drag component selection
!                    for CCPP.
!                    Namelist options:
!                    do_gsl_drag_ls_bl - logical flag for mesoscale GWD + blocking
!                    do_gsl_drag_ss - logical flag for small-scale GWD
!                    do_gsl_drag_tofd - logical flag for turbulent form drag
!                    Compile-time options (changed from integer switches to logical flags):
!                    gwd_opt_ms = : mesoscale GWD (active if == .true.)
!                    gwd_opt_bl = : blocking drag (active if == .true.)
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
   integer, intent(in) :: KPBL(:)
   real(kind=kind_phys), intent(in) :: deltim, G, CP, RD, RV, cdmbgwd(:)
   real(kind=kind_phys), intent(inout) :: dtend(:,:,:)
   logical, intent(in) :: ldiag3d
   integer, intent(in) :: dtidx(:,:), index_of_temperature,      &
     &  index_of_process_orographic_gwd, index_of_x_wind, index_of_y_wind

   integer              ::  kpblmax
   integer, parameter   ::  ims=1, kms=1, its=1, kts=1
   real(kind=kind_phys), intent(in) ::  fv, pi
   real(kind=kind_phys) ::  rcl, cdmb
   real(kind=kind_phys) ::  g_inv

   real(kind=kind_phys), intent(inout) ::                        &
     &                   dudt(:,:),dvdt(:,:),                &
     &                   dtdt(:,:)
   real(kind=kind_phys), intent(out) :: rdxzb(:)
   real(kind=kind_phys), intent(in) ::                           &
     &                            u1(:,:),v1(:,:),           &
     &                            t1(:,:),q1(:,:),           &
     &                            PHII(:,:),prsl(:,:),     &
     &                            prslk(:,:),PHIL(:,:)
   real(kind=kind_phys), intent(in) ::  prsi(:,:),           &
     &                                  del(:,:)
   real(kind=kind_phys), intent(in) ::   var(:),oc1(:),        &
     &                                   oa4(:,:),ol4(:,:),    &
     &                                   dx(:)
   real(kind=kind_phys), intent(in) ::   varss(:),oc1ss(:),    &
     &                              oa4ss(:,:),ol4ss(:,:)
   real(kind=kind_phys), intent(in) :: THETA(:),SIGMA(:),      &
     &                                 GAMMA(:),ELVMAX(:)

! added for small-scale orographic wave drag
   real(kind=kind_phys), dimension(im,km) :: utendwave,vtendwave,thx,thvx
   real(kind=kind_phys), intent(in) ::     br1(:),              &
     &                                     hpbl(:),             &
     &                                     slmsk(:)
   real(kind=kind_phys), dimension(im)    :: govrth,xland
   !real(kind=kind_phys), dimension(im,km) :: dz2
   real(kind=kind_phys)                   :: tauwavex0,tauwavey0,  &
     &                                     XNBV,density,tvcon,hpbl2
   integer                          ::     kpbl2,kvar
   !real(kind=kind_phys), dimension(im,km+1)         ::     zq      ! = PHII/g
   real(kind=kind_phys), dimension(im,km)           ::     zl      ! = PHIL/g

!SPP
   real(kind=kind_phys), dimension(im) :: var_stoch, varss_stoch, &
                                       varmax_fd_stoch
   real(kind=kind_phys), intent(in) :: spp_wts_gwd(:,:)
   integer, intent(in) :: spp_gwd

   real(kind=kind_phys), dimension(im)              :: rstoch

!Output:
   real(kind=kind_phys), intent(inout) ::                        &
     &                      dusfc(:),   dvsfc(:)
!Output (optional):
   real(kind=kind_phys), intent(inout) ::                        &
     &                      dusfc_ms(:),dvsfc_ms(:),             &
     &                      dusfc_bl(:),dvsfc_bl(:),             &
     &                      dusfc_ss(:),dvsfc_ss(:),             &
     &                      dusfc_fd(:),dvsfc_fd(:)
   real(kind=kind_phys), intent(inout) ::                        &
     &         dtaux2d_ms(:,:),dtauy2d_ms(:,:),                  &
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
      do_gsl_drag_ls_bl,       & ! mesoscale gravity wave drag and blocking
      do_gsl_drag_ss,          & ! small-scale gravity wave drag (Steeneveld et al. 2008)
      do_gsl_drag_tofd           ! form drag (Beljaars et al. 2004, QJRMS)

! Flag for diagnostic outputs
      logical, intent(in) :: ldiag_ugwp

! Flag for sequential update of u and v between
! LSGWD + BLOCKING and SSGWD + TOFD calculations
      logical, intent(in) :: ugwp_seq_update

! More variables for sequential updating of winds
      ! Updated winds
      real(kind=kind_phys), dimension(im,km) :: uwnd1, vwnd1
      real(kind=kind_phys) :: tmp1, tmp2   ! temporary variables

! Additional flags
      logical, parameter ::    &
      gwd_opt_ms      = .true.,     & ! mesoscale gravity wave drag
      gwd_opt_bl      = .true.,     & ! blocking drag
      gsd_diss_ht_opt = .false.       ! dissipation heating

! Parameters for bounding the scale-adaptive variability:
! Small-scale GWD + turbulent form drag
   real(kind=kind_phys), parameter      :: dxmin_ss = 1000.,                     &
     &                     dxmax_ss = 12000.  ! min,max range of tapering (m)
! Mesoscale GWD + blocking
   real(kind=kind_phys), parameter      :: dxmin_ms = 3000.,                     &
     &                     dxmax_ms = 13000.  ! min,max range of tapering (m)
   real(kind=kind_phys), dimension(im)  :: ss_taper, ls_taper ! small- and meso-scale tapering factors (-)
!
! Variables for limiting topographic standard deviation (var)
   real(kind=kind_phys), parameter      :: varmax_ss = 50.,      &  ! varmax_ss not used
                           varmax_fd = 500.,                     &
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
   real(kind=kind_phys), parameter       ::  sgmalolev  = 0.5  ! max sigma lvl for dtfac
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
   logical              ::  prop_test
!
   real(kind=kind_phys) ::  onebgrcs
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
                            taud_ms(im,km),taud_bl(im,km),        &
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
   real(kind=kind_phys)                 :: delx,dely
   real(kind=kind_phys)                 :: dxy4(im,4),dxy4p(im,4)
   real(kind=kind_phys)                 :: dxy(im),dxyp(im)
   real(kind=kind_phys)                 :: ol4p(4),olp(im),od(im)
   real(kind=kind_phys)                 :: taufb(im,km+1)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   integer :: udtend, vdtend, Tdtend

   ! Calculate inverse of gravitational acceleration
   g_inv = 1./G

   ! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   ! Initialize local variables
   var_temp2 = 0.
   udtend = -1
   vdtend = -1
   Tdtend = -1

   if(ldiag3d) then
      udtend = dtidx(index_of_x_wind,index_of_process_orographic_gwd)
      vdtend = dtidx(index_of_y_wind,index_of_process_orographic_gwd)
      Tdtend = dtidx(index_of_temperature,index_of_process_orographic_gwd)
   endif


   ! Initialize winds for sequential updating.
   ! These are for optional sequential updating of the wind between the
   ! LSGWD+BLOCKING and SSGWD+TOFD steps.  They are placeholders
   ! for the u1,v1 winds that are updated within the scheme if
   ! ugwp_seq_update == T, otherwise, they retain the values
   ! passed in to the subroutine.
   uwnd1(:,:) = u1(:,:)
   vwnd1(:,:) = v1(:,:)

!--------------------------------------------------------------------
! SCALE-ADPTIVE PARAMETER FROM GFS GWD SCHEME
!--------------------------------------------------------------------
!     parameter (cdmb = 1.0)     ! non-dim sub grid mtn drag Amp (*j*)
! non-dim sub grid mtn drag Amp (*j*)
!     cdmb = 1.0/float(IMX/192)
!     cdmb = 192.0/float(IMX)
      ! New cdmbgwd addition for GSL blocking drag
      cdmb = 1.0
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
   onebgrcs = 1._kind_phys/g*rcs

   do i=1,im
      if (slmsk(i)==1. .or. slmsk(i)==2.) then !sea/land/ice mask (=0/1/2) in FV3
         xland(i)=1.0                          !but land/water = (1/2) in this module
      else
         xland(i)=2.0
      endif
      RDXZB(i)  = 0.0
   enddo

!--- calculate scale-aware tapering factors
do i=1,im
   if ( dx(i) .ge. dxmax_ms ) then
      ls_taper(i) = 1.
   else
      if ( dx(i) .le. dxmin_ms) then
         ls_taper(i) = 0.
      else
         ls_taper(i) = 0.5 * ( SIN(pi*(dx(i)-0.5*(dxmax_ms+dxmin_ms))/  &
                                  (dxmax_ms-dxmin_ms)) + 1. )
      endif
   endif
enddo

! Remove ss_tapering
ss_taper(:) = 1.

! SPP, if spp_gwd is 0, no perturbations are applied.
if ( spp_gwd==1 ) then
  do i = its,im
    var_stoch(i)   = var(i)   + var(i)*0.75*spp_wts_gwd(i,1)
    varss_stoch(i) = varss(i) + varss(i)*0.75*spp_wts_gwd(i,1)
    varmax_fd_stoch(i) = varmax_fd + varmax_fd*0.75*spp_wts_gwd(i,1)
  enddo
else
  do i = its,im
    var_stoch(i)   = var(i)
    varss_stoch(i) = varss(i)
    varmax_fd_stoch(i) = varmax_fd
  enddo
endif

!--- calculate length of grid for flow-blocking drag
!
do i=1,im
   delx   = dx(i)
   dely   = dx(i)
   dxy4(i,1)  = delx
   dxy4(i,2)  = dely
   dxy4(i,3)  = sqrt(delx*delx + dely*dely)
   dxy4(i,4)  = dxy4(i,3)
   dxy4p(i,1) = dxy4(i,2)
   dxy4p(i,2) = dxy4(i,1)
   dxy4p(i,3) = dxy4(i,4)
   dxy4p(i,4) = dxy4(i,3)
enddo
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
       taud_ms(i,k) = 0.0
       taud_bl(i,k) = 0.0
       dtaux2d(i,k) = 0.0
       dtauy2d(i,k) = 0.0
     enddo
   enddo
!
   do i = its,im
     xlinv(i)     = 1.0/xl
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
     zlowtop(i) = 2. * var_stoch(i)
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
     dxy(i)  = dxy4(i,MOD(nwd-1,4)+1)
     dxyp(i) = dxy4p(i,MOD(nwd-1,4)+1)
   enddo
!
! END INITIALIZATION; BEGIN GWD CALCULATIONS:
!

IF ( (do_gsl_drag_ls_bl).and.                            &
     (gwd_opt_ms.or.gwd_opt_bl) ) then

   do i=its,im

      if ( ls_taper(i).GT.1.E-02 ) then

!
!---  saving richardson number in usqj for migwdi
!
         do k = kts,km-1
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
!
!----compute the "low level" or 1/3 wind magnitude (m/s)
!
         ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0)
         rulow(i) = 1./ulow(i)
!
         do k = kts,km-1
            velco(i,k)  = (0.5*rcs) * ((u1(i,k)+u1(i,k+1)) * ubar(i)       &
                                     + (v1(i,k)+v1(i,k+1)) * vbar(i))
            velco(i,k)  = velco(i,k) * rulow(i)
            if ((velco(i,k).lt.veleps) .and. (velco(i,k).gt.0.)) then
               velco(i,k) = veleps
            endif
         enddo
!
!  no drag when critical level in the base layer
!
         ldrag(i) = velco(i,1).le.0.
!
!  no drag when velco.lt.0
!
         do k = kpblmin,kpblmax
            if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. velco(i,k).le.0.
         enddo
!
!  no drag when bnv2.lt.0
!
         do k = kts,kpblmax
            if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. bnv2(i,k).lt.0.
         enddo
!
!-----the low level weighted average ri is stored in usqj(1,1; im)
!-----the low level weighted average n**2 is stored in bnv2(1,1; im)
!---- this is called bnvl2 in phys_gwd_alpert_sub not bnv2
!---- rdelks (del(k)/delks) vert ave factor so we can * instead of /
!
         wtkbj     = (prsl(i,1)-prsl(i,2)) * delks1(i)
         bnv2(i,1) = wtkbj * bnv2(i,1)
         usqj(i,1) = wtkbj * usqj(i,1)
!
         do k = kpblmin,kpblmax
            if (k .lt. kbl(i)) then
               rdelks    = (prsl(i,k)-prsl(i,k+1)) * delks1(i)
               bnv2(i,1) = bnv2(i,1) + bnv2(i,k) * rdelks
               usqj(i,1) = usqj(i,1) + usqj(i,k) * rdelks
            endif
         enddo
!
         ldrag(i) = ldrag(i) .or. bnv2(i,1).le.0.0
         ldrag(i) = ldrag(i) .or. ulow(i).eq.1.0
         ldrag(i) = ldrag(i) .or. var_stoch(i) .le. 0.0
!
!  set all ri low level values to the low level value
!
         do k = kpblmin,kpblmax
            if (k .lt. kbl(i)) usqj(i,k) = usqj(i,1)
         enddo
!
         if (.not.ldrag(i))   then
            bnv(i) = sqrt( bnv2(i,1) )
            fr(i) = bnv(i)  * rulow(i) * 2. * var_stoch(i) * od(i)
            fr(i) = min(fr(i),frmax)
            xn(i)  = ubar(i) * rulow(i)
            yn(i)  = vbar(i) * rulow(i)
         endif
!
!  compute the base level stress and store it in taub
!  calculate enhancement factor, number of mountains & aspect
!  ratio const. use simplified relationship between standard
!  deviation & critical hgt

         if (.not. ldrag(i))   then
            efact    = (oa(i) + 2.) ** (ce*fr(i)/frc)
            efact    = min( max(efact,efmin), efmax )
!!!!!!! cleff (effective grid length) is highly tunable parameter
!!!!!!! the bigger (smaller) value produce weaker (stronger) wave drag
!WRF         cleff    = sqrt(dxy(i)**2. + dxyp(i)**2.)
!WRF         cleff    = 3. * max(dx(i),cleff)
            coefm(i) = (1. + ol(i)) ** (oa(i)+1.)
!WRF         xlinv(i) = coefm(i) / cleff
            xlinv(i) = coefm(i) * cleff
            tem      = fr(i) * fr(i) * oc1(i)
            gfobnv   = gmax * tem / ((tem + cg)*bnv(i))
            if ( gwd_opt_ms ) then
               taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)           &
                           * ulow(i) * gfobnv * efact
            else     ! We've gotten what we need for the blocking scheme
               taub(i) = 0.0
            end if
         else
            taub(i) = 0.0
            xn(i)   = 0.0
            yn(i)   = 0.0
         endif

      endif   ! (ls_taper(i).GT.1.E-02)

   enddo  ! do i=its,im

ENDIF   ! (do_gsl_drag_ls_bl).and.(gwd_opt_ms.or.gwd_opt_bl)




!=======================================================
! Mesoscale GWD + blocking
!=======================================================
IF ( (do_gsl_drag_ls_bl).and.(gwd_opt_ms) ) THEN

   do i=its,im

      if ( ls_taper(i).GT.1.E-02 ) then

!
!   now compute vertical structure of the stress.
         do k = kts,kpblmax
            if (k .le. kbl(i)) taup(i,k) = taub(i)
         enddo
!
         do k = kpblmin, km-1                   ! vertical level k loop!
            kp1 = k + 1
!
!   unstablelayer if ri < ric
!   unstable layer if upper air vel comp along surf vel <=0 (crit lay)
!   at (u-c)=0. crit layer exists and bit vector should be set (.le.)
!
            if (k .ge. kbl(i)) then
               icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric)            &
                                     .or. (velco(i,k) .le. 0.0)
               brvf(i)  = max(bnv2(i,k),bnv2min) ! brunt-vaisala frequency squared
               brvf(i)  = sqrt(brvf(i))          ! brunt-vaisala frequency
            endif
!
            if (k .ge. kbl(i) .and. (.not. ldrag(i)))   then
               if (.not.icrilv(i) .and. taup(i,k) .gt. 0.0 ) then
                  temv = 1.0 / velco(i,k)
                  tem1 = coefm(i)/dxy(i)*(ro(i,kp1)+ro(i,k))*brvf(i)*   &
                                                      velco(i,k)*0.5
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
!
         if(lcap.lt.km) then
            do klcap = lcapp1,km
               taup(i,klcap) = prsi(i,klcap) / prsi(i,lcap) * taup(i,lcap)
            enddo
         endif

      endif  ! if ( ls_taper(i).GT.1.E-02 )

   enddo  ! do i=its,im

ENDIF  ! (do_gsl_drag_ls_bl).and.(gwd_opt_ms)
!===============================================================
!COMPUTE BLOCKING COMPONENT                                     
!===============================================================
IF ( do_gsl_drag_ls_bl .and. gwd_opt_bl ) THEN

   do i=its,im

      if ( ls_taper(i).GT.1.E-02 ) then

         if (.not.ldrag(i)) then
!
!------- determine the height of flow-blocking layer
!
            kblk = 0
            pe = 0.0
            do k = km, kpblmin, -1
               if(kblk.eq.0 .and. k.le.komax(i)) then
                  pe = pe + bnv2(i,k)*(zl(i,komax(i))-zl(i,k))*      &
                            del(i,k)/g/ro(i,k)
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
               ! New cdmbgwd addition for GSL blocking drag
               taufb(i,kts) = cdmb * 0.5 * roll(i) * coefm(i) /            &
                                 max(dxmax_ms,dxy(i))**2 * cd * dxyp(i) *  &
                                 olp(i) * zblk * ulow(i)**2
               tautem = taufb(i,kts)/float(kblk-kts)
               do k = kts+1, kblk
                  taufb(i,k) = taufb(i,k-1) - tautem
               enddo
!
!----------sum orographic GW stress and flow-blocking stress
!
              ! taup(i,:) = taup(i,:) + taufb(i,:)   ! Keep taup and taufb separate for now
            endif

         endif   ! if (.not.ldrag(i))

      endif   ! if ( ls_taper(i).GT.1.E-02 )

   enddo   ! do i=its,im

ENDIF   ! IF ( do_gsl_drag_ls_bl .and. gwd_opt_bl )
!===========================================================
IF ( (do_gsl_drag_ls_bl) .and.                                       &
     (gwd_opt_ms .OR. gwd_opt_bl) ) THEN 

   do i=its,im

      if ( ls_taper(i) .GT. 1.E-02 ) then

!
!  calculate - (g)*d(tau)/d(pressure) and deceleration terms dtaux, dtauy
!
!  First, set taup (momentum flux) at model top equal to that of the layer
!  interface just below the top, i.e., taup(km)
!  The idea is to allow the momentum flux to go out the 'top'.  This
!  ensures there is no GWD force at the top layer.
!
         taup(i,km+1) = taup(i,km)
         do k = kts,km
            taud_ms(i,k) = (taup(i,k+1) - taup(i,k)) * csg / del(i,k)
            taud_bl(i,k) = (taufb(i,k+1) - taufb(i,k)) * csg / del(i,k)
         enddo
!
!
!  if the gravity wave drag + blocking would force a critical line
!  in the layers below pressure-based 'sigma' level = sgmalolev during the next deltim
!  timestep, then only apply drag until that critical line is reached, i.e.,
!  reduce drag to limit resulting wind components to zero
!  Note: 'sigma' = prsi(k)/prsi(k=1), where prsi(k=1) is the surface pressure
!
         do k = kts,kpblmax-1
            if (prsi(i,k).ge.sgmalolev*prsi(i,1)) then
               if ((taud_ms(i,k)+taud_bl(i,k)).ne.0.)                   &
                  dtfac(i) = min(dtfac(i),abs(velco(i,k)                &
                       /(deltim*rcs*(taud_ms(i,k)+taud_bl(i,k)))))
            else
               exit
            endif
         enddo
!
         do k = kts,km
            taud_ms(i,k)  = taud_ms(i,k)*dtfac(i)* ls_taper(i) *(1.-rstoch(i))
            taud_bl(i,k)  = taud_bl(i,k)*dtfac(i)* ls_taper(i) *(1.-rstoch(i))

            dtaux  = taud_ms(i,k) * xn(i)
            dtauy  = taud_ms(i,k) * yn(i)
            dtauxb = taud_bl(i,k) * xn(i)
            dtauyb = taud_bl(i,k) * yn(i)

            !add blocking and mesoscale contributions to tendencies 
            tmp1 = dtaux + dtauxb
            tmp2 = dtauy + dtauyb
            dudt(i,k)  = tmp1 + dudt(i,k)
            dvdt(i,k)  = tmp2 + dvdt(i,k)

            ! Update winds if sequential updating is selected
            ! and SSGWD and TOFD will be calculated
            ! Note:  uwnd1 and vwnd1 replace u1 and u2,respectively,
            ! for the SSGWD and TOFD calculations
            if ( ugwp_seq_update .and. (do_gsl_drag_ss.or.do_gsl_drag_tofd) ) then
               uwnd1(i,k) = uwnd1(i,k) + tmp1*deltim
               vwnd1(i,k) = vwnd1(i,k) + tmp2*deltim
            endif

            if ( gsd_diss_ht_opt ) then
               ! Calculate dissipation heating
               ! Initial kinetic energy (at t0-dt)
               eng0 = 0.5*( (rcs*u1(i,k))**2. + (rcs*v1(i,k))**2. )
               ! Kinetic energy after wave-breaking/flow-blocking
               eng1 = 0.5*( (rcs*(u1(i,k)+(dtaux+dtauxb)*deltim))**2 +  &
                            (rcs*(v1(i,k)+(dtauy+dtauyb)*deltim))**2 )
               ! Modify theta tendency
               dtdt(i,k) = dtdt(i,k) + max((eng0-eng1),0.0)/cp/deltim
               if ( Tdtend>0 ) then
                  dtend(i,k,Tdtend) = dtend(i,k,Tdtend) + max((eng0-eng1),0.0)/cp
               endif
            endif

            dusfc(i) = dusfc(i) - onebgrcs * ( taud_ms(i,k)*xn(i)*del(i,k) + &
                                    taud_bl(i,k)*xn(i)*del(i,k) )
            dvsfc(i) = dvsfc(i) - onebgrcs * ( taud_ms(i,k)*yn(i)*del(i,k) + &
                                    taud_bl(i,k)*yn(i)*del(i,k) )
            if(udtend>0) then
               dtend(i,k,udtend) = dtend(i,k,udtend) + (taud_ms(i,k) *  &
                    xn(i) + taud_bl(i,k) * xn(i)) * deltim
            endif
            if(vdtend>0) then
               dtend(i,k,vdtend) = dtend(i,k,vdtend) + (taud_ms(i,k) *  &
                    yn(i) + taud_bl(i,k) * yn(i)) * deltim
            endif

         enddo

         if ( ldiag_ugwp ) then
            do k = kts,km
               dtaux2d_ms(i,k) = taud_ms(i,k) * xn(i)
               dtauy2d_ms(i,k) = taud_ms(i,k) * yn(i)
               dtaux2d_bl(i,k) = taud_bl(i,k) * xn(i)
               dtauy2d_bl(i,k) = taud_bl(i,k) * yn(i)
               dusfc_ms(i)  = dusfc_ms(i) - onebgrcs * dtaux2d_ms(i,k) * del(i,k)
               dvsfc_ms(i)  = dvsfc_ms(i) - onebgrcs * dtauy2d_ms(i,k) * del(i,k)
               dusfc_bl(i)  = dusfc_bl(i) - onebgrcs * dtaux2d_bl(i,k) * del(i,k)
               dvsfc_bl(i)  = dvsfc_bl(i) - onebgrcs * dtauy2d_bl(i,k) * del(i,k)
            enddo
         endif

      endif    ! if ( ls_taper(i) .GT. 1.E-02 )

   enddo   ! do i=its,im

ENDIF  ! (do_gsl_drag_ls_bl).and.(gwd_opt_ms .OR. gwd_opt_bl)


!====================================================================
! Calculate small-scale gravity wave drag for stable boundary layer
!====================================================================
  XNBV=0.
  tauwavex0=0.
  tauwavey0=0.
  density=1.2
  utendwave=0.
  vtendwave=0.
!
IF ( do_gsl_drag_ss ) THEN

   do i=its,im

      if ( ss_taper(i).GT.1.E-02 ) then
   !
   ! calculating potential temperature
   !
         do k = kts,km
            thx(i,k) = t1(i,k)/prslk(i,k)
         enddo
   !
         do k = kts,km
            tvcon = (1.+fv*q1(i,k))
            thvx(i,k) = thx(i,k)*tvcon
         enddo

         hpbl2 = hpbl(i)+10.
         kpbl2 = kpbl(i)
         !kvar = MIN(kpbl, k-level of var)
         kvar = 1
         do k=kts+1,MAX(kpbl(i),kts+1)
!            IF (zl(i,k)>2.*var(i) .or. zl(i,k)>2*varmax) then
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
         if((xland(i)-1.5).le.0. .and. 2.*varss_stoch(i).le.hpbl(i))then
            if(br1(i).gt.0. .and. thvx(i,kpbl2)-thvx(i,kts) > 0.)then
              ! Modify xlinv to represent wave number of "typical" small-scale topography 
!              cleff_ss    = 3. * max(dx(i),cleff_ss)
!              cleff_ss    = 10. * max(dxmax_ss,cleff_ss)
!               cleff_ss    = 0.1 * 12000.
              xlinv(i) = 0.001*pi   ! 2km horizontal wavelength
              !govrth(i)=g/(0.5*(thvx(i,kpbl(i))+thvx(i,kts)))
              govrth(i)=g/(0.5*(thvx(i,kpbl2)+thvx(i,kts)))
              !XNBV=sqrt(govrth(i)*(thvx(i,kpbl(i))-thvx(i,kts))/hpbl(i))
              XNBV=sqrt(govrth(i)*(thvx(i,kpbl2)-thvx(i,kts))/hpbl2)
!
              ! check for possibility of vertical wave propagation
              ! (avoids division by zero if uwnd1(i,kpbl2).eq.0.)
              if (uwnd1(i,kpbl2).eq.0.) then
                 prop_test = .true.
              elseif (abs(XNBV/uwnd1(i,kpbl2)).gt.xlinv(i)) then
                 prop_test = .true.
              else
                 prop_test = .false.
              endif
              if (prop_test) then
                ! Remove limit on varss_stoch
                var_temp = varss_stoch(i)
                ! Note:  This is a semi-implicit treatment of the time differencing
                var_temp2 = 0.5*XNBV*xlinv(i)*(2.*var_temp)**2*ro(i,kvar)  ! this is greater than zero
                tauwavex0=var_temp2*uwnd1(i,kvar)/(1.+var_temp2*deltim)
                tauwavex0=tauwavex0*ss_taper(i)
              else
                tauwavex0=0.
              endif
!
              ! check for possibility of vertical wave propagation
              ! (avoids division by zero if vwnd1(i,kpbl2).eq.0.)
              if (vwnd1(i,kpbl2).eq.0.) then
                 prop_test = .true.
              elseif (abs(XNBV/vwnd1(i,kpbl2)).gt.xlinv(i)) then
                 prop_test = .true.
              else
                 prop_test = .false.
              endif
              if (prop_test) then
                ! Remove limit on varss_stoch
                var_temp = varss_stoch(i)
                ! Note:  This is a semi-implicit treatment of the time differencing
                var_temp2 = 0.5*XNBV*xlinv(i)*(2.*var_temp)**2*ro(i,kvar)  ! this is greater than zero
                tauwavey0=var_temp2*vwnd1(i,kvar)/(1.+var_temp2*deltim)
                tauwavey0=tauwavey0*ss_taper(i)
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

         do k = kts,km
            dudt(i,k)  = dudt(i,k) + utendwave(i,k)
            dvdt(i,k)  = dvdt(i,k) + vtendwave(i,k)
            dusfc(i)   = dusfc(i) - onebgrcs * utendwave(i,k) * del(i,k)
            dvsfc(i)   = dvsfc(i) - onebgrcs * vtendwave(i,k) * del(i,k)
         enddo
         if(udtend>0) then
            dtend(i,kts:km,udtend) = dtend(i,kts:km,udtend) + utendwave(i,kts:km)*deltim
         endif
         if(vdtend>0) then
            dtend(i,kts:km,vdtend) = dtend(i,kts:km,vdtend) + vtendwave(i,kts:km)*deltim
         endif
         if ( ldiag_ugwp ) then
            do k = kts,km
               dusfc_ss(i) = dusfc_ss(i) + utendwave(i,k) * del(i,k)
               dvsfc_ss(i) = dvsfc_ss(i) + vtendwave(i,k) * del(i,k)
               dtaux2d_ss(i,k) = utendwave(i,k)
               dtauy2d_ss(i,k) = vtendwave(i,k)
            enddo
         endif

      endif  ! if (ss_taper(i).GT.1.E-02)

   enddo  ! i=its,im

ENDIF  ! if (do_gsl_drag_ss)


!===================================================================
! Topographic Form Drag from Beljaars et al. (2004, QJRMS, equ. 16):
!===================================================================
IF ( do_gsl_drag_tofd ) THEN

   do i=its,im

      if ( ss_taper(i).GT.1.E-02 ) then

         utendform=0.
         vtendform=0.

         IF ((xland(i)-1.5) .le. 0.) then
            !(IH*kflt**n1)**-1 = (0.00102*0.00035**-1.9)**-1 = 0.00026615161
            var_temp = MIN(varss_stoch(i),varmax_fd_stoch(i)) +                &
                       MAX(0.,beta_fd*(varss_stoch(i)-varmax_fd_stoch(i)))
            a1=0.00026615161*var_temp**2
!           a1=0.00026615161*MIN(varss(i),varmax)**2
!           a1=0.00026615161*(0.5*varss(i))**2
           ! k1**(n1-n2) = 0.003**(-1.9 - -2.8) = 0.003**0.9 = 0.005363
            a2=a1*0.005363
            ! Beljaars H_efold
            H_efold = 1500.
            DO k=kts,km
               wsp=SQRT(uwnd1(i,k)**2 + vwnd1(i,k)**2)
               ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759
               var_temp = 0.0759*EXP(-(zl(i,k)/H_efold)**1.5)*a2*          &
                                 zl(i,k)**(-1.2)*ss_taper(i) ! this is greater than zero
               !  Note:  This is a semi-implicit treatment of the time differencing
               !  per Beljaars et al. (2004, QJRMS)
               utendform(i,k) = - var_temp*wsp*uwnd1(i,k)/(1. + var_temp*deltim*wsp)
               vtendform(i,k) = - var_temp*wsp*vwnd1(i,k)/(1. + var_temp*deltim*wsp)
               !IF(zl(i,k) > 4000.) exit
            ENDDO
         ENDIF

         do k = kts,km
            dudt(i,k)  = dudt(i,k) + utendform(i,k)
            dvdt(i,k)  = dvdt(i,k) + vtendform(i,k)
            dusfc(i)   = dusfc(i) - onebgrcs * utendform(i,k) * del(i,k)
            dvsfc(i)   = dvsfc(i) - onebgrcs * vtendform(i,k) * del(i,k)
         enddo
         if(udtend>0) then
            dtend(i,kts:km,udtend) = dtend(i,kts:km,udtend) + utendform(i,kts:km)*deltim
         endif
         if(vdtend>0) then
            dtend(i,kts:km,vdtend) = dtend(i,kts:km,vdtend) + vtendform(i,kts:km)*deltim
         endif
         if ( ldiag_ugwp ) then
            do k = kts,km
               dtaux2d_fd(i,k) = utendform(i,k)
               dtauy2d_fd(i,k) = vtendform(i,k)
               dusfc_fd(i) = dusfc_fd(i) + utendform(i,k) * del(i,k)
               dvsfc_fd(i) = dvsfc_fd(i) + vtendform(i,k) * del(i,k)
            enddo
         endif

      endif   ! if (ss_taper(i).GT.1.E-02)

   enddo  ! i=its,im

ENDIF  ! if (do_gsl_drag_tofd)



if ( ldiag_ugwp ) then
   !  Finalize dusfc and dvsfc diagnostics for gsl small-scale drag components
   dusfc_ss(:) = -onebgrcs * dusfc_ss(:)
   dvsfc_ss(:) = -onebgrcs * dvsfc_ss(:)
   dusfc_fd(:) = -onebgrcs * dusfc_fd(:)
   dvsfc_fd(:) = -onebgrcs * dvsfc_fd(:)
endif
!
   return
   end subroutine drag_suite_run
!-------------------------------------------------------------------
!
!> @}

      end module drag_suite

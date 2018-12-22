!> \file module_mynnsfc_wrapper.F90
!!  Contains all of the code related to running the MYNN surface layer scheme 

      MODULE mynnsfc_wrapper

      contains

      subroutine mynnsfc_wrapper_init ()
      end subroutine mynnsfc_wrapper_init

      subroutine mynnsfc_wrapper_finalize ()
      end subroutine mynnsfc_wrapper_finalize

!!
!> \brief This scheme (1) performs pre-mynnsfc work, (20 runs the mynn sfc layer scheme, and (3) performs post-mynnsfc work
#if 0
!! \section arg_table_mynnsfc_wrapper_run Argument Table
!! | local_name          | standard_name                                                               | long_name                                             | units         | rank | type      |    kind   | intent | optional |
!! |---------------------|-----------------------------------------------------------------------------|-------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | ix                  | horizontal_dimension                                                        | horizontal dimension                                  | count         |    0 | integer   |           | in     | F        |
!! | im                  | horizontal_loop_extent                                                      | horizontal loop extent                                | count         |    0 | integer   |           | in     | F        |
!! | levs                | vertical_dimension                                                          | vertical layer dimension                              | count         |    0 | integer   |           | in     | F        |
!! | iter                | ccpp_loop_counter                                                           | loop counter for subcycling loops in CCPP             | index         |    0 | integer   |           | in     | F        |
!! | flag_init           | flag_for_first_time_step                                                 | flag signaling first time step for time integration loop | flag          |    0 | logical   |           | in     | F        |
!! | flag_restart        | flag_for_restart                                                            | flag for restart (warmstart) or coldstart             | flag          |    0 | logical   |           | in     | F        |
!! | delt                | time_step_for_physics                                                       | time step for physics                                 | s             |    0 | real      | kind_phys | in     | F        |
!! | dx                  | cell_size                                                                   | size of the grid cell                                 | m             |    1 | real      | kind_phys | in     | F        |
!! | u                   | x_wind                                                                      | x component of layer wind                             | m s-1         |    2 | real      | kind_phys | in     | F        |
!! | v                   | y_wind                                                                      | y component of layer wind                             | m s-1         |    2 | real      | kind_phys | in     | F        |
!! | t3d                 | air_temperature                                                             | layer mean air temperature                            | K             |    2 | real      | kind_phys | in     | F        |
!! | qvsh                | water_vapor_specific_humidity                                               | water vapor specific humidity                         | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | qc                  | cloud_condensed_water_mixing_ratio                                          | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate)   | kg kg-1   |    2 | real   | kind_phys | in     | F        |
!! | prsl                | air_pressure                                                                | mean layer pressure                                   | Pa            |    2 | real      | kind_phys | in     | F        |
!! | phii                | geopotential_at_interface                                                   | geopotential at model layer interfaces                | m2 s-2        |    2 | real      | kind_phys | in     | F        |
!! | exner               | dimensionless_exner_function_at_model_layers                                | Exner function at layers                              | none          |    2 | real      | kind_phys | in     | F        |
!! | tsq                 | t_prime_squared                                                             | temperature fluctuation squared                       | K2            |    2 | real      | kind_phys | in     | F        |
!! | qsq                 | q_prime_squared                                                             | water vapor fluctuation squared                       | kg2 kg-2      |    2 | real      | kind_phys | in     | F        |
!! | cov                 | t_prime_q_prime                                                             | covariance of temperature and moisture                | K kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | el_pbl              | mixing_length                                                               | mixing length in meters                               | m             |    2 | real      | kind_phys | in     | F        |
!! | Sh3D                | stability_function_for_heat                                                 | stability function for heat                           | none          |    2 | real      | kind_phys | in     | F        |
!! | QC_BL               | subgrid_cloud_mixing_ratio_pbl                                              | subgrid cloud cloud mixing ratio from PBL scheme      | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | CLDFRA_BL           | subgrid_cloud_fraction_pbl                                                  | subgrid cloud fraction from PBL scheme                | frac          |    2 | real      | kind_phys | in     | F        |
!! | ps                  | surface_air_pressure                                                        | surface pressure                                      | Pa            |    1 | real      | kind_phys | in     | F        |
!! | PBLH                | atmosphere_boundary_layer_thickness                                         | PBL thickness                                         | m             |    1 | real      | kind_phys | in     | F        |
!! | slmsk               | sea_land_ice_mask_real                                                      | landmask: sea/land/ice=0/1/2                          | flag          |    1 | real      | kind_phys | in     | F        |
!! | tsk                 | surface_skin_temperature                                                    | surface temperature                                   | K             |    1 | real      | kind_phys | in     | F        |
!! | qsfc                | surface_specific_humidity                                                   | surface air saturation specific humidity              | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | snowd               | surface_snow_thickness_water_equivalent                                     | water equivalent snow depth over land                 | mm            |    1 | real      | kind_phys | in     | F        |
!! | zorl                | surface_roughness_length                                                    | surface roughness length in cm                        | cm            |    1 | real      | kind_phys | inout  | F        |
!! | ust                 | surface_friction_velocity                                                   | boundary layer parameter                              | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | ustm                | surface_friction_velocity_drag                                              | friction velocity isolated for momentum only          | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | zol                 | surface_stability_parameter                                                 | monin obukhov surface stability parameter             | none          |    1 | real      | kind_phys | inout  | F        |
!! | mol                 | theta_star                                                                  | temperature flux divided by ustar (temperature scale) | K             |    1 | real      | kind_phys | inout  | F        |
!! | rmol                | reciprocal_of_obukhov_length                                                | one over obukhov length                               | m-1           |    1 | real      | kind_phys | inout  | F        |
!! | fm                  | Monin-Obukhov_similarity_function_for_momentum                              | Monin-Obukhov similarity parameter for momentum       | none          |    1 | real      | kind_phys | inout  | F        |
!! | fh                  | Monin-Obukhov_similarity_function_for_heat                                  | Monin-Obukhov similarity parameter for heat           | none          |    1 | real      | kind_phys | inout  | F        |
!! | fm10                | Monin-Obukhov_similarity_function_for_momentum_at_10m                       | Monin-Obukhov similarity parameter for momentum       | none          |    1 | real      | kind_phys | inout  | F        |
!! | fh2                 | Monin-Obukhov_similarity_function_for_heat_at_2m                            | Monin-Obukhov similarity parameter for heat           | none          |    1 | real      | kind_phys | inout  | F        |
!! | wspd                | wind_speed_at_lowest_model_layer                                            | wind speed at lowest model level                      | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | br                  | bulk_richardson_number_at_lowest_model_level                                | bulk Richardson number at the surface                 | none          |    1 | real      | kind_phys | inout  | F        |
!! | ch                  | surface_drag_wind_speed_for_momentum_in_air                                 | momentum exchange coefficient                         | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | hflx                | kinematic_surface_upward_sensible_heat_flux                                 | kinematic surface upward sensible heat flux           | K m s-1       |    1 | real      | kind_phys | inout  | F        |
!! | QFX                 | kinematic_surface_upward_latent_heat_flux                                   | kinematic surface upward latent heat flux             | kg kg-1 m s-1 |    1 | real      | kind_phys | inout  | F        |
!! | lh                  | surface_latent_heat                                                         | latent heating at the surface (pos = up)              | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | flhc                | surface_exchange_coefficient_for_heat                                       | surface exchange coefficient for heat                 | W m-2 K-1     |    1 | real      | kind_phys | inout  | F        |
!! | flqc                | surface_exchange_coefficient_for_moisture                                   | surface exchange coefficient for moisture             | kg m-2 s-1    |    1 | real      | kind_phys | inout  | F        |
!! | u10                 | x_wind_at_10m                                                               | 10 meter u wind speed                                 | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | v10                 | y_wind_at_10m                                                               | 10 meter v wind speed                                 | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | th2                 | potential_temperature_at_2m                                                 | 2 meter potential temperature                         | K             |    1 | real      | kind_phys | inout  | F        |
!! | t2                  | temperature_at_2m                                                           | 2 meter temperature                                   | K             |    1 | real      | kind_phys | inout  | F        |
!! | q2                  | specific_humidity_at_2m                                                     | 2 meter specific humidity                             | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | wstar               | surface_wind_enhancement_due_to_convection                                  | surface wind enhancement due to convection            | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | chs2                | surface_exchange_coefficient_for_heat_at_2m                                 | exchange coefficient for heat at 2 meters             | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | cqs2                | surface_exchange_coefficient_for_moisture_at_2m                             | exchange coefficient for moisture at 2 meters         | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | cda                 | surface_drag_coefficient_for_momentum_in_air                                | surface exchange coeff for momentum                   | none          |    1 | real      | kind_phys | inout  | F        |
!! | cka                 | surface_drag_coefficient_for_heat_and_moisture_in_air                       | surface exchange coeff heat & moisture                | none          |    1 | real      | kind_phys | inout  | F        |
!! | stress              | surface_wind_stress                                                         | surface wind stress                                   | m2 s-2        |    1 | real      | kind_phys | inout  | F        |
!! | bl_mynn_cloudpdf    | cloudpdf                                                                    | flag to determine which cloud PDF to use              | flag          |    0 | integer   |           | in     | F        |
!! | icloud_bl           | couple_sgs_clouds_to_radiation_flag                                         | flag for coupling sgs clouds to radiation             | flag          |    0 | integer   |           | in     | F        |
!! | lprnt               | flag_print                                                                  | control flag for diagnostic print out                 | flag          |    0 | logical   |           | none   | F        |
!! | errmsg              | ccpp_error_message                                                          | error message for error handling in CCPP              | none          |    0 | character | len=*     | out    | F        |
!! | errflg              | ccpp_error_flag                                                             | error flag for error handling in CCPP                 | flag          |    0 | integer   |           | out    | F        |
!!
#endif
!###===================================================================
SUBROUTINE mynnsfc_wrapper_run(         &
     &  ix,im,levs,                     &
     &  iter,flag_init,flag_restart,    &
     &  delt,dx,                        &
     &  u, v, t3d, qvsh, qc, prsl, phii,&
     &  exner, tsq, qsq, cov, sh3d,     &
     &  el_pbl, qc_bl, cldfra_bl,       &
     &  ps, PBLH, slmsk, TSK,           &
     &  QSFC, snowd,                    &
     &  zorl,UST,USTM, ZOL,MOL,RMOL,    &
     &  fm, fh, fm10, fh2, WSPD, br, ch,&
     &  HFLX, QFX, LH, FLHC, FLQC,      &
     &  U10, V10, TH2, T2, Q2,          &
     &  wstar, CHS2, CQS2,              &
     &  cda, cka, stress,               &
!     &  CP, G, ROVCP, R, XLV,           &
!     &  SVP1, SVP2, SVP3, SVPT0,        &
!     &  EP1,EP2,KARMAN,                 &
     &  icloud_bl, bl_mynn_cloudpdf,    &
     &  lprnt, errmsg, errflg           )


! should be moved to inside the mynn:
      use machine , only : kind_phys
!      use funcphys, only : fpvs

      use physcons, only : cp     => con_cp,              &
     &                     g      => con_g,               &
     &                     r_d    => con_rd,              &
     &                     r_v    => con_rv,              &
     &                     cpv    => con_cvap,            &
     &                     cliq   => con_cliq,            &
     &                     Cice   => con_csol,            &
     &                     rcp    => con_rocp,            &
     &                     XLV    => con_hvap,            &
     &                     XLF    => con_hfus,            &
     &                     EP_1   => con_fvirt,           &
     &                     EP_2   => con_eps

      USE module_sf_mynn, only : SFCLAY_mynn 

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 
!  ---  constant parameters:
!      real(kind=kind_phys), parameter :: rvovrd  = r_v/r_d
      real(kind=kind_phys), parameter :: karman  = 0.4
!      real(kind=kind_phys), parameter :: XLS     = 2.85E6
!      real(kind=kind_phys), parameter :: p1000mb = 100000.
      real(kind=kind_phys), parameter :: SVP1    = 0.6112
      real(kind=kind_phys), parameter :: SVP2    = 17.67
      real(kind=kind_phys), parameter :: SVP3    = 29.65
      real(kind=kind_phys), parameter :: SVPT0   = 273.15

!-------------------------------------------------------------------
!For WRF:
!-------------------------------------------------------------------
!  USE module_model_constants, only: &
!       &karman, g, p1000mb, &
!       &cp, r_d, r_v, rcp, xlv, xlf, xls, &
!      &svp1, svp2, svp3, svpt0, ep_1, ep_2, rvovrd, &
!       &cpv, cliq, cice

!-------------------------------------------------------------------
!For reference
!   REAL    , PARAMETER :: karman       = 0.4
!   REAL    , PARAMETER :: g            = 9.81
!   REAL    , PARAMETER :: r_d          = 287.
!   REAL    , PARAMETER :: cp           = 7.*r_d/2.
!   REAL    , PARAMETER :: r_v          = 461.6
!   REAL    , PARAMETER :: cpv          = 4.*r_v
!   REAL    , PARAMETER :: cliq         = 4190.
!   REAL    , PARAMETER :: Cice         = 2106.
!   REAL    , PARAMETER :: rcp          = r_d/cp
!   REAL    , PARAMETER :: XLS          = 2.85E6
!   REAL    , PARAMETER :: XLV          = 2.5E6
!   REAL    , PARAMETER :: XLF          = 3.50E5
!   REAL    , PARAMETER :: p1000mb      = 100000.
!   REAL    , PARAMETER :: rvovrd       = r_v/r_d
!   REAL    , PARAMETER :: SVP1         = 0.6112
!   REAL    , PARAMETER :: SVP2         = 17.67
!   REAL    , PARAMETER :: SVP3         = 29.65
!   REAL    , PARAMETER :: SVPT0        = 273.15
!   REAL    , PARAMETER :: EP_1         = R_v/R_d-1.
!   REAL    , PARAMETER :: EP_2         = R_d/R_v

  REAL, PARAMETER :: xlvcp=xlv/cp, xlscp=(xlv+xlf)/cp, ev=xlv, rd=r_d, &
       &rk=cp/rd, svp11=svp1*1.e3, p608=ep_1, ep_3=1.-ep_2, g_inv=1/g


  character(len=*), intent(out) :: errmsg
  integer, intent(out) :: errflg

! NAMELIST OPTIONS (INPUT):
      INTEGER, INTENT(IN) ::                                &
     &       bl_mynn_cloudpdf,                              &
     &       icloud_bl

!MISC CONFIGURATION OPTIONS
      INTEGER, PARAMETER ::                                 &
     &       spp_pbl  = 0,                                  &
     &       isftcflx = 0,                                  &
     &       iz0tlnd  = 0,                                  &
     &       isfflx   = 1

!MYNN-1D
      REAL    :: delt
      INTEGER :: im, ix, levs
      INTEGER :: iter, k, i, itimestep
      LOGICAL :: flag_init,flag_restart,lprnt
      INTEGER :: IDS,IDE,JDS,JDE,KDS,KDE,                   &
     &            IMS,IME,JMS,JME,KMS,KME,                  &
     &            ITS,ITE,JTS,JTE,KTS,KTE

!MYNN-3D
      real(kind=kind_phys), dimension(im,levs+1) :: phii
      real(kind=kind_phys), dimension(im,levs) ::           &
     &        exner, PRSL,                                  &
     &        u, v, t3d, qvsh, qc,                          &
     &        Sh3D, EL_PBL, EXCH_H,                         &
     &        qc_bl, cldfra_bl,                             &
     &        Tsq, Qsq, Cov
     !LOCAL
      real(kind=kind_phys), dimension(im,levs) ::           &
     &        dz, rho, th, qv,                              &
     &        pattern_spp_pbl

!MYNN-2D
      real(kind=kind_phys), dimension(im) ::                &
     &        dx, pblh, slmsk, tsk, qsfc, ps,               &
     &        zorl, ust, ustm, hflx, qfx, br, wspd, snowd,  &
     &        FLHC, FLQC, U10, V10, TH2, T2, Q2,            &
     &        CHS2, CQS2, rmol, zol, mol, ch,               &
     &        fm, fh, fm10, fh2,                            &
     &        lh, cda, cka, stress, wstar
     !LOCAL
      real, dimension(im) ::                                &
     &        qcg, hfx, znt, ts, snowh, psim, psih,         &
     &        chs, ck, cd, mavail, regime, xland, GZ1OZ0       

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lprnt) then
         write(0,*)"=============================================="
         write(0,*)"in mynn surface layer wrapper..."
         write(0,*)"flag_init=",flag_init
         write(0,*)"flag_restart=",flag_restart
         write(0,*)"iter=",iter
      endif

      ! If initialization is needed and mynnsfc_wrapper is called
      ! in a subcycling loop, then test for (flag_init==.T. .and. iter==1);
      ! initialization in sfclay_mynn is triggered by itimestep == 1
      ! DH* TODO: Use flag_restart to distinguish which fields need
      ! to be initialized and which are read from restart files
      if (flag_init.and.iter==1) then
          itimestep = 1
      else
          itimestep = 2
      endif

      !prep MYNN-only variables
            do k=1,levs
              do i=1,im
                 dz(i,k)=(phii(i,k+1) - phii(i,k))*g_inv
                 th(i,k)=t3d(i,k)/exner(i,k)
                 !qc(i,k)=MAX(qgrs(i,k,ntcw),0.0)
                 qv(i,k)=qvsh(i,k)/(1.0 - qvsh(i,k))
                 rho(i,k)=prsl(i,k)/(r_d*t3d(i,k)) !gt0(i,k))
                 pattern_spp_pbl(i,k)=0.0
              enddo
            enddo
            do i=1,im
                if (slmsk(i)==1. .or. slmsk(i)==2.)then !sea/land/ice mask (=0/1/2) in FV3
                  xland(i)=1.0                          !but land/water = (1/2) in SFCLAY_mynn
                else
                  xland(i)=2.0
                endif
!                ust(i) = sqrt(stress(i))
                !ch(i)=0.0
                HFX(i)=hflx(i)*rho(i,1)*cp
                !QFX(i)=evap(i)
                !wstar(i)=0.0
                qcg(i)=0.0
                snowh(i)=snowd(i)*800. !mm -> m
                znt(i)=zorl(i)*0.01    !cm -> m?
                ts(i)=tsk(i)/exner(i,1)  !theta
!                qsfc(i)=qss(i)
!                ps(i)=pgr(i)
!                wspd(i)=wind(i)
                mavail(i)=1.0  !????
            enddo

      if (lprnt) then
          write(0,*)"CALLING SFCLAY_mynn; input:"
          print*,"T:",t3d(1,1),t3d(1,2),t3d(1,3)
          print*,"TH:",th(1,1),th(1,2),th(1,3)
          print*,"rho:",rho(1,1),rho(1,2),rho(1,3)
          print*,"u:",u(1,1:3)
          !print*,"qv:",qv(1,1:3,1)
          print*,"p:",prsl(1,1)," snowh=",snowh(1)
          print*,"dz:",dz(1,1)," qsfc=",qsfc(1)
          print*,"rmol:",rmol(1)," ust:",ust(1)
          print*,"Tsk:",tsk(1)," Thetasurf:",ts(1)
          print*,"HFX:",hfx(1)," qfx",qfx(1)
          print*,"qsfc:",qsfc(1)," ps:",ps(1)
          print*,"wspd:",wspd(1),"br=",br(1)
          print*,"znt:",znt(1)," delt=",delt
          print*,"im=",im," levs=",levs
          print*,"flag_init=",flag_init !," ntcw=",ntcw!," ntk=",ntk
          print*,"flag_restart=",flag_restart !," ntcw=",ntcw!," ntk=",ntk
          print*,"iter=",iter
          !print*,"ncld=",ncld," ntrac(gq0)=",ntrac
          print*,"zlvl(1)=",dz(1,1)*0.5
          print*,"PBLH=",pblh(1)," xland=",xland(1)
       endif


        CALL SFCLAY_mynn(                                                 &
                     u3d=u,v3d=v,t3d=t3d,qv3d=qv,p3d=prsl,dz8w=dz,        &
                     CP=cp,G=g,ROVCP=rcp,R=r_d,XLV=xlv,                   &
                     PSFCPA=ps,CHS=chs,CHS2=chs2,CQS2=cqs2,               &
                     ZNT=znt,UST=ust,PBLH=pblh,MAVAIL=mavail,             &
                     ZOL=zol,MOL=mol,REGIME=regime,psim=psim,psih=psih,   &
                     psix=fm,psit=fh,psix10=fm10,psit2=fh2,               &
!                     fm=psix,fh=psit,fm10=psix10,fh2=psit2,               &
                     XLAND=xland,HFX=hfx,QFX=qfx,LH=lh,TSK=tsk,           &
                     FLHC=flhc,FLQC=flqc,QSFC=qsfc,RMOL=rmol,             &
                     U10=u10,V10=v10,TH2=th2,T2=t2,Q2=q2,SNOWH=snowh,     &
                     GZ1OZ0=GZ1OZ0,WSPD=wspd,BR=br,ISFFLX=isfflx,DX=dx,   &
                     SVP1=svp1,SVP2=svp2,SVP3=svp3,SVPT0=svpt0,           &
                     EP1=ep_1,EP2=ep_2,KARMAN=karman,                     &
                     itimestep=itimestep,ch=ch,                           &
                     th3d=th,pi3d=exner,qc3d=qc,rho3d=rho,                &
                     tsq=tsq,qsq=qsq,cov=cov,sh3d=sh3d,el_pbl=el_pbl,     &
                     qcg=qcg,wstar=wstar,                                 &
                     icloud_bl=icloud_bl,qc_bl=qc_bl,cldfra_bl=cldfra_bl, &
                     spp_pbl=spp_pbl,pattern_spp_pbl=pattern_spp_pbl,     &
                     ids=1,ide=im, jds=1,jde=1, kds=1,kde=levs,           &
                     ims=1,ime=im, jms=1,jme=1, kms=1,kme=levs,           &
                     its=1,ite=im, jts=1,jte=1, kts=1,kte=levs,           &
                     ustm=ustm, ck=ck, cka=cka, cd=cd, cda=cda,           &
                     isftcflx=isftcflx, iz0tlnd=iz0tlnd,                  &
                     bl_mynn_cloudpdf=bl_mynn_cloudpdf      )


     ! POST MYNN SURFACE LAYER (INTERSTITIAL) WORK:
        do i = 1, im
           hflx(i)=hfx(i)/(rho(i,1)*cp)
           !QFX(i)=evap(i)                                                                                                                                                                                                      
           zorl(i)=znt(i)*100.             !m -> cm
           stress(i) = ust(i)**2
        enddo


      if (lprnt) then
         print*
         print*,"finished with mynn_surface layer; output:"
         print*,"xland=",xland(1)," cda=",cda(1)
         print*,"rmol:",rmol(1)," ust:",ust(1)
         print*,"Tsk:",tsk(1)," Thetasurf:",ts(1)
         print*,"HFX:",hfx(1)," qfx",qfx(1)
         print*,"qsfc:",qsfc(1)," ps:",ps(1)
         print*,"wspd:",wspd(1)," br=",br(1)
         print*,"znt:",znt(1),"pblh:",pblh(1)
         print*,"FLHC=",FLHC(1)," CHS=",CHS(1)
         print*
      endif


  END SUBROUTINE mynnsfc_wrapper_run

!###=================================================================

END MODULE mynnsfc_wrapper

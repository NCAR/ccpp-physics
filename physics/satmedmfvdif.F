!> \file satmedmfvdif.F
!! This file  contains the CCPP-compliant SATMEDMF scheme which computes
!! subgrid vertical turbulence mixing using scale-aware TKE-based moist
!! eddy-diffusion mass-flux (EDMF) parameterization (by Jongil Han).
!!
!! Subroutine 'satmedmfvdif_run' computes subgrid vertical turbulence mixing
!! using scale-aware TKE-based moist eddy-diffusion mass-flux (EDMF) parameterization
!! (by Jongil Han).
!!
!! - For the convective boundary layer, the scheme adopts
!! EDMF parameterization (Siebesma et al., 2007) to take
!! into account nonlocal transport by large eddies (mfpblt.f).
!!
!! - A new mass-flux parameterization for stratocumulus-top-induced turbulence
!! mixing has been introduced (previously, it was eddy diffusion form)
!! [mfscu.f].
!!
!! - For local turbulence mixing, a TKE closure model is used.
!
!----------------------------------------------------------------------

      module satmedmfvdif
      contains

      subroutine satmedmfvdif_init ()
      end subroutine satmedmfvdif_init

      subroutine satmedmfvdif_finalize ()
      end subroutine satmedmfvdif_finalize

!> \defgroup SATMEDMF FV3GFS satmedmfvdif_run Main
!! \brief This subroutine contains all of the logic for the 
!! scale-aware TKE-based moist eddy-diffusion mass-flux (EDMF) scheme.
!!
!> \section arg_table_satmedmfvdif_run Argument Table
!! | local_name     | standard_name                                                               | long_name                                             | units         | rank | type      |    kind   | intent | optional |
!! |----------------|-----------------------------------------------------------------------------|-------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | ix             | horizontal_dimension                                                        | horizontal dimension                                  | count         |    0 | integer   |           | in     | F        |
!! | im             | horizontal_loop_extent                                                      | horizontal loop extent                                | count         |    0 | integer   |           | in     | F        |
!! | km             | vertical_dimension                                                          | vertical layer dimension                              | count         |    0 | integer   |           | in     | F        |
!! | ntrac          | number_of_vertical_diffusion_tracers                                        | number of tracers to diffuse vertically               | count         |    0 | integer   |           | in     | F        |
!! | ntcw           | index_for_liquid_cloud_condensate                                           | tracer index for cloud condensate (or liquid water)   | index         |    0 | integer   |           | in     | F        |
!! | ntiw           | index_for_ice_cloud_condensate_vertical_diffusion_tracer             | tracer index for ice water in the vertically diffused tracer array | index         |    0 | integer   |           | in     | F        |
!! | ntke           | index_for_turbulent_kinetic_energy_vertical_diffusion_tracer | index for turbulent kinetic energy in the vertically diffused tracer array | index   |    0 | integer   |           | in     | F        |
!! | grav           | gravitational_acceleration                                                  | gravitational acceleration                            | m s-2         |    0 | real      | kind_phys | in     | F        |
!! | rd             | gas_constant_dry_air                                                        | ideal gas constant for dry air                        | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | cp             | specific_heat_of_dry_air_at_constant_pressure                               | specific heat of dry air at constant pressure         | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | rv             | gas_constant_water_vapor                                                    | ideal gas constant for water vapor                    | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | hvap           | latent_heat_of_vaporization_of_water_at_0C                                  | latent heat of evaporation/sublimation                | J kg-1        |    0 | real      | kind_phys | in     | F        |
!! | hfus           | latent_heat_of_fusion_of_water_at_0C                                        | latent heat of fusion                                 | J kg-1        |    0 | real      | kind_phys | in     | F        |
!! | fv             | ratio_of_vapor_to_dry_air_gas_constants_minus_one                           | (rv/rd) - 1 (rv = ideal gas constant for water vapor) | none          |    0 | real      | kind_phys | in     | F        |
!! | eps            | ratio_of_dry_air_to_water_vapor_gas_constants                               | rd/rv                                                 | none          |    0 | real      | kind_phys | in     | F        |
!! | epsm1          | ratio_of_dry_air_to_water_vapor_gas_constants_minus_one                     | (rd/rv) - 1                                           | none          |    0 | real      | kind_phys | in     | F        |
!! | dv             | tendency_of_y_wind_due_to_model_physics                                     | updated tendency of the y wind                        | m s-2         |    2 | real      | kind_phys | inout  | F        |
!! | du             | tendency_of_x_wind_due_to_model_physics                                     | updated tendency of the x wind                        | m s-2         |    2 | real      | kind_phys | inout  | F        |
!! | tdt            | tendency_of_air_temperature_due_to_model_physics                            | updated tendency of the temperature                   | K s-1         |    2 | real      | kind_phys | inout  | F        |
!! | rtg            | tendency_of_vertically_diffused_tracer_concentration                        | updated tendency of the tracers due to vertical diffusion in PBL scheme | kg kg-1 s-1 | 3 | real | kind_phys | inout | F |
!! | u1             | x_wind                                                                      | x component of layer wind                             | m s-1         |    2 | real      | kind_phys | in     | F        |
!! | v1             | y_wind                                                                      | y component of layer wind                             | m s-1         |    2 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature                                                             | layer mean air temperature                            | K             |    2 | real      | kind_phys | in     | F        |
!! | q1             | vertically_diffused_tracer_concentration                                    | tracer concentration diffused by PBL scheme           | kg kg-1       |    3 | real      | kind_phys | in     | F        |
!! | swh            | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step | total sky shortwave heating rate                      | K s-1         |    2 | real      | kind_phys | in     | F        |
!! | hlw            | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step  | total sky longwave heating rate                       | K s-1         |    2 | real      | kind_phys | in     | F        |
!! | xmu            | zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes                | zenith angle temporal adjustment factor for shortwave | none          |    1 | real      | kind_phys | in     | F        |
!! | garea          | cell_area                                                                   | area of the grid cell                                 | m2            |    1 | real      | kind_phys | in     | F        |
!! | psk            | dimensionless_exner_function_at_lowest_model_interface                      | dimensionless Exner function at the surface interface | none          |    1 | real      | kind_phys | in     | F        |
!! | rbsoil         | bulk_richardson_number_at_lowest_model_level                                | bulk Richardson number at the surface                 | none          |    1 | real      | kind_phys | in     | F        |
!! | zorl           | surface_roughness_length                                                    | surface roughness length in cm                        | cm            |    1 | real      | kind_phys | in     | F        |
!! | u10m           | x_wind_at_10m                                                               | x component of wind at 10 m                           | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v10m           | y_wind_at_10m                                                               | y component of wind at 10 m                           | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | fm             | Monin-Obukhov_similarity_function_for_momentum                              | Monin-Obukhov similarity function for momentum        | none          |    1 | real      | kind_phys | in     | F        |
!! | fh             | Monin-Obukhov_similarity_function_for_heat                                  | Monin-Obukhov similarity function for heat            | none          |    1 | real      | kind_phys | in     | F        |
!! | tsea           | surface_skin_temperature                                                    | surface skin temperature                              | K             |    1 | real      | kind_phys | in     | F        |
!! | heat           | kinematic_surface_upward_sensible_heat_flux                                 | kinematic surface upward sensible heat flux           | K m s-1       |    1 | real      | kind_phys | in     | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux                                   | kinematic surface upward latent heat flux             | kg kg-1 m s-1 |    1 | real      | kind_phys | in     | F        |
!! | stress         | surface_wind_stress                                                         | surface wind stress                                   | m2 s-2        |    1 | real      | kind_phys | in     | F        |
!! | spd1           | wind_speed_at_lowest_model_layer                                            | wind speed at lowest model level                      | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | kpbl           | vertical_index_at_top_of_atmosphere_boundary_layer                          | PBL top model level index                             | index         |    1 | integer   |           | out    | F        |
!! | prsi           | air_pressure_at_interface                                                   | air pressure at model layer interfaces                | Pa            |    2 | real      | kind_phys | in     | F        |
!! | del            | air_pressure_difference_between_midlayers                                   | pres(k) - pres(k+1)                                   | Pa            |    2 | real      | kind_phys | in     | F        |
!! | prsl           | air_pressure                                                                | mean layer pressure                                   | Pa            |    2 | real      | kind_phys | in     | F        |
!! | prslk          | dimensionless_exner_function_at_model_layers                                | Exner function at layers                              | none          |    2 | real      | kind_phys | in     | F        |
!! | phii           | geopotential_at_interface                                                   | geopotential at model layer interfaces                | m2 s-2        |    2 | real      | kind_phys | in     | F        |
!! | phil           | geopotential                                                                | geopotential at model layer centers                   | m2 s-2        |    2 | real      | kind_phys | in     | F        |
!! | delt           | time_step_for_physics                                                       | time step for physics                                 | s             |    0 | real      | kind_phys | in     | F        |
!! | dspheat        | flag_TKE_dissipation_heating                                                | flag for using TKE dissipation heating                | flag          |    0 | logical   |           | in     | F        |
!! | dusfc          | instantaneous_surface_x_momentum_flux                                       | x momentum flux                                       | Pa            |    1 | real      | kind_phys | out    | F        |
!! | dvsfc          | instantaneous_surface_y_momentum_flux                                       | y momentum flux                                       | Pa            |    1 | real      | kind_phys | out    | F        |
!! | dtsfc          | instantaneous_surface_upward_sensible_heat_flux                             | surface upward sensible heat flux                     | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | dqsfc          | instantaneous_surface_upward_latent_heat_flux                               | surface upward latent heat flux                       | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | hpbl           | atmosphere_boundary_layer_thickness                                         | PBL thickness                                         | m             |    1 | real      | kind_phys | out    | F        |
!! | kinver         | index_of_highest_temperature_inversion                                      | index of highest temperature inversion                | index         |    1 | integer   |           | in     | F        |
!! | xkzm_m         | atmosphere_momentum_diffusivity_background                                  | background value of momentum diffusivity              | m2 s-1        |    0 | real      | kind_phys | in     | F        |
!! | xkzm_h         | atmosphere_heat_diffusivity_background                                      | background value of heat diffusivity                  | m2 s-1        |    0 | real      | kind_phys | in     | F        |
!! | xkzm_s         | diffusivity_background_sigma_level                                          | sigma level threshold for background diffusivity      | none          |    0 | real      | kind_phys | in     | F        |
!! | errmsg         | ccpp_error_message                                                          | error message for error handling in CCPP              | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                             | error flag for error handling in CCPP                 | flag          |    0 | integer   |           | out    | F        |
!!
      subroutine satmedmfvdif_run(ix,im,km,ntrac,ntcw,ntiw,ntke,
     &     grav,rd,cp,rv,hvap,hfus,fv,eps,epsm1,
     &     dv,du,tdt,rtg,u1,v1,t1,q1,swh,hlw,xmu,garea,
     &     psk,rbsoil,zorl,u10m,v10m,fm,fh,
     &     tsea,heat,evap,stress,spd1,kpbl,
     &     prsi,del,prsl,prslk,phii,phil,delt,
     &     dspheat,dusfc,dvsfc,dtsfc,dqsfc,hpbl,
     &     kinver,xkzm_m,xkzm_h,xkzm_s,errmsg,errflg)
!
      use machine  , only : kind_phys
      use funcphys , only : fpvs
!
      implicit none
!
!----------------------------------------------------------------------
      integer, intent(in)  :: ix, im, km, ntrac, ntcw, ntiw, ntke
      integer, intent(in)  :: kinver(im)
      integer, intent(out) :: kpbl(im)
!
      real(kind=kind_phys), intent(in) :: grav,rd,cp,rv,hvap,hfus,fv,
     &                                    eps,epsm1
      real(kind=kind_phys), intent(in) :: delt, xkzm_m, xkzm_h, xkzm_s
      real(kind=kind_phys), intent(inout) :: dv(im,km),     du(im,km),
     &                     tdt(im,km),    rtg(im,km,ntrac)
      real(kind=kind_phys), intent(in) ::
     &                     u1(ix,km),     v1(ix,km),
     &                     t1(ix,km),     q1(ix,km,ntrac),
     &                     swh(ix,km),    hlw(ix,km),
     &                     xmu(im),       garea(im),
     &                     psk(ix),       rbsoil(im),
     &                     zorl(im),      tsea(im),
     &                     u10m(im),      v10m(im),
     &                     fm(im),        fh(im),
     &                     evap(im),      heat(im),
     &                     stress(im),    spd1(im),
     &                     prsi(ix,km+1), del(ix,km),
     &                     prsl(ix,km),   prslk(ix,km),
     &                     phii(ix,km+1), phil(ix,km)
      real(kind=kind_phys), intent(out) ::
     &                     dusfc(im),     dvsfc(im),
     &                     dtsfc(im),     dqsfc(im),
     &                     hpbl(im) 
!
      logical, intent(out)  :: dspheat
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!          flag for tke dissipative heating
!
!----------------------------------------------------------------------
!***
!***  local variables
!***
      integer i,is,k,kk,n,km1,kmpbl,kmscu,ntrac1
      integer lcld(im),kcld(im),krad(im),mrad(im)
      integer kx1(im), kpblx(im)
!
      real(kind=kind_phys) tke(im,km),  tkeh(im,km-1)
!
      real(kind=kind_phys) theta(im,km),thvx(im,km),  thlvx(im,km),
     &                     qlx(im,km),  thetae(im,km),thlx(im,km),
     &                     slx(im,km),  svx(im,km),   qtx(im,km),
     &                     tvx(im,km),  pix(im,km),   radx(im,km-1),
     &                     dku(im,km-1),dkt(im,km-1), dkq(im,km-1),
     &                     cku(im,km-1),ckt(im,km-1)
!
      real(kind=kind_phys) plyr(im,km), rhly(im,km),  cfly(im,km),
     &                     qstl(im,km)
!
      real(kind=kind_phys) dtdz1(im), gdx(im),
     &                     phih(im),  phim(im),    prn(im,km-1),
     &                     rbdn(im),  rbup(im),    thermal(im),
     &                     ustar(im), wstar(im),   hpblx(im),
     &                     ust3(im),  wst3(im),
     &                     z0(im),    crb(im),
     &                     hgamt(im), hgamq(im),
     &                     wscale(im),vpert(im),
     &                     zol(im),   sflux(im),   radj(im),
     &                     tx1(im),   tx2(im)
!
      real(kind=kind_phys) radmin(im)
!
      real(kind=kind_phys) zi(im,km+1),  zl(im,km),   zm(im,km),
     &                     xkzo(im,km-1),xkzmo(im,km-1),
     &                     xkzm_hx(im),  xkzm_mx(im),
     &                     rdzt(im,km-1),
     &                     al(im,km-1),  ad(im,km),   au(im,km-1),
     &                     f1(im,km),    f2(im,km*(ntrac-1))
!
      real(kind=kind_phys) elm(im,km),   ele(im,km),  rle(im,km-1),
     &                     ckz(im,km),   chz(im,km),
     &                     diss(im,km-1),prod(im,km-1),
     &                     bf(im,km-1),  shr2(im,km-1),
     &                     xlamue(im,km-1), xlamde(im,km-1),
     &                     gotvx(im,km), rlam(im,km-1)
!
!   variables for updrafts (thermals)
!
      real(kind=kind_phys) tcko(im,km),  qcko(im,km,ntrac),
     &                     ucko(im,km),  vcko(im,km),
     &                     buou(im,km),  xmf(im,km)
!
!   variables for stratocumulus-top induced downdrafts
!
      real(kind=kind_phys) tcdo(im,km),  qcdo(im,km,ntrac),
     &                     ucdo(im,km),  vcdo(im,km),
     &                     buod(im,km),  xmfd(im,km)
!
      logical  pblflg(im), sfcflg(im), flg(im)
      logical  scuflg(im), pcnvflg(im)
      logical  mlenflg
!
!  pcnvflg: true for unstable pbl
!
      real(kind=kind_phys) aphi16,  aphi5,
     &                     wfac,    cfac,
     &                     gamcrt,  gamcrq, sfcfrac,
     &                     conq,    cont,   conw,
     &                     dsdz2,   dsdzt,  dkmax,
     &                     dsig,    dt2,    dtodsd,
     &                     dtodsu,  g,      factor, dz,
     &                     gocp,    gravi,  zol1,   zolcru,
     &                     buop,    shrp,   dtn,    cdtn,
     &                     prnum,   prmax,  prmin,  prtke,
     &                     prscu,   dw2,    dw2min, zk,
     &                     elmfac,  elefac, dspmax,
     &                     alp,     clwt,   cql,
     &                     f0,      robn,   crbmin, crbmax,
     &                     es,      qs,     value,  onemrh,
     &                     cfh,     gamma,  elocp,  el2orc,
     &                     epsi,    beta,   chx,    cqx,
     &                     rdt,     rdz,    qmin,   qlmin,
     &                     ri,      rimin,
     &                     rbcr,    rbint,  tdzmin,
     &                     rlmn,    rlmx,   elmx,
     &                     ttend,   utend,  vtend,  qtend,
     &                     zfac,    zfmin,  vk,     spdk2,
     &                     tkmin,   xkzinv, dspfac, xkgdx,
     &                     zlup,    zldn,   bsum,
     &                     tem,     tem1,   tem2,
     &                     ptem,    ptem0,  ptem1,  ptem2
!
      real(kind=kind_phys) ck0, ck1, ch0, ch1, ce0, rchck
!
      real(kind=kind_phys) qlcr, zstblmax
!
      real(kind=kind_phys) h1 
!!
      parameter(wfac=7.0,cfac=4.5)
      parameter(gamcrt=3.,gamcrq=0.,sfcfrac=0.1)
      parameter(vk=0.4,rimin=-100.)
      parameter(rbcr=0.25,zolcru=-0.02,tdzmin=1.e-3)
      parameter(rlmn=30.,rlmx=500.,elmx=500.)
      parameter(prmin=0.25,prmax=4.0,prtke=1.0,prscu=0.67)
      parameter(f0=1.e-4,crbmin=0.15,crbmax=0.35)
      parameter(tkmin=1.e-9,dspfac=0.5,dspmax=10.0)
      parameter(qmin=1.e-8,qlmin=1.e-12,zfmin=1.e-8)
      parameter(aphi5=5.,aphi16=16.)
      parameter(elmfac=1.0,elefac=1.0,cql=100.)
      parameter(dw2min=1.e-4,dkmax=1000.,xkgdx=25000.)
      parameter(qlcr=3.5e-5,zstblmax=2500.,xkzinv=0.15)
      parameter(h1=0.33333333)
      parameter(ck0=0.4,ck1=0.15,ch0=0.4,ch1=0.15,ce0=0.4)
      parameter(rchck=1.5,cdtn=25.)

      gravi=1.0/grav
      g=grav
      gocp=g/cp
      cont=cp/g
      conq=hvap/g
      conw=1.0/g  ! for del in pa
!     parameter(cont=1000.*cp/g,conq=1000.*hvap/g,conw=1000./g) !kpa
      elocp=hvap/cp
      el2orc=hvap*hvap/(rv*cp)

!
!************************************************************************
! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      dt2  = delt
      rdt = 1. / dt2
!
! the code is written assuming ntke=ntrac
! if ntrac > ntke, the code needs to be modified
!
      ntrac1 = ntrac - 1
      km1 = km - 1
      kmpbl = km / 2
      kmscu = km / 2
!
      do k=1,km
        do i=1,im
          zi(i,k) = phii(i,k) * gravi
          zl(i,k) = phil(i,k) * gravi
          xmf(i,k) = 0.
          xmfd(i,k) = 0.
          buou(i,k) = 0.
          buod(i,k) = 0.
          ckz(i,k) = ck1
          chz(i,k) = ch1
        enddo
      enddo
      do i=1,im
        zi(i,km+1) = phii(i,km+1) * gravi
      enddo
      do k=1,km
        do i=1,im
          zm(i,k) = zi(i,k+1)
        enddo
      enddo
! horizontal grid size
      do i=1,im
        gdx(i) = sqrt(garea(i))
      enddo
!
      do k=1,km
        do i=1,im
          tke(i,k) = max(q1(i,k,ntke), tkmin)
        enddo
      enddo
      do k=1,km1
        do i=1,im
          tkeh(i,k) = 0.5 * (tke(i,k) + tke(i,k+1))
        enddo
      enddo
!
      do k = 1,km1
        do i=1,im
          rdzt(i,k) = 1.0 / (zl(i,k+1) - zl(i,k))
          prn(i,k)  = 1.0
        enddo
      enddo
!
! set background diffusivities as a function of
!  horizontal grid size with xkzm_h & xkzm_m for gdx >= 25km
!  and 0.01 for gdx=5m, i.e.,
!  xkzm_hx = 0.01 + (xkzm_h - 0.01)/(xkgdx-5.) * (gdx-5.)
!  xkzm_mx = 0.01 + (xkzm_h - 0.01)/(xkgdx-5.) * (gdx-5.)
!
      do i=1,im
        kx1(i) = 1
        tx1(i) = 1.0 / prsi(i,1)
        tx2(i) = tx1(i)
        if(gdx(i) >= xkgdx) then
          xkzm_hx(i) = xkzm_h
          xkzm_mx(i) = xkzm_m
        else
          tem  = 1. / (xkgdx - 5.)
          tem1 = (xkzm_h - 0.01) * tem
          tem2 = (xkzm_m - 0.01) * tem
          ptem = gdx(i) - 5.
          xkzm_hx(i) = 0.01 + tem1 * ptem
          xkzm_mx(i) = 0.01 + tem2 * ptem
        endif
      enddo
      do k = 1,km1
        do i=1,im
          xkzo(i,k)  = 0.0
          xkzmo(i,k) = 0.0
          if (k < kinver(i)) then
!                                  vertical background diffusivity
            ptem      = prsi(i,k+1) * tx1(i)
            tem1      = 1.0 - ptem
            tem1      = tem1 * tem1 * 10.0
            xkzo(i,k) = xkzm_hx(i) * min(1.0, exp(-tem1))
!                                 vertical background diffusivity for momentum
            if (ptem >= xkzm_s) then
              xkzmo(i,k) = xkzm_mx(i)
              kx1(i)     = k + 1
            else
              if (k == kx1(i) .and. k > 1) tx2(i) = 1.0 / prsi(i,k)
              tem1 = 1.0 - prsi(i,k+1) * tx2(i)
              tem1 = tem1 * tem1 * 5.0
              xkzmo(i,k) = xkzm_mx(i) * min(1.0, exp(-tem1))
            endif
          endif
        enddo
      enddo
!
      do i = 1,im
         z0(i)    = 0.01 * zorl(i)
         dusfc(i) = 0.
         dvsfc(i) = 0.
         dtsfc(i) = 0.
         dqsfc(i) = 0.
         kpbl(i) = 1
         hpbl(i) = 0.
         kpblx(i) = 1
         hpblx(i) = 0.
         pblflg(i)= .true.
         sfcflg(i)= .true.
         if(rbsoil(i) > 0.) sfcflg(i) = .false.
         pcnvflg(i)= .false.
         scuflg(i)= .true.
         if(scuflg(i)) then
           radmin(i)= 0.
           mrad(i)  = km1
           krad(i)  = 1
           lcld(i)  = km1
           kcld(i)  = km1
         endif
      enddo
!
      do k=1,km
        do i=1,im
          pix(i,k)   = psk(i) / prslk(i,k)
          theta(i,k) = t1(i,k) * pix(i,k)
          if(ntiw > 0) then
            tem = max(q1(i,k,ntcw),qlmin)
            tem1 = max(q1(i,k,ntiw),qlmin)
            qlx(i,k) = tem + tem1
            ptem = hvap*tem + (hvap+hfus)*tem1
            slx(i,k)   = cp * t1(i,k) + phil(i,k) - ptem
          else
            qlx(i,k) = max(q1(i,k,ntcw),qlmin)
            slx(i,k)   = cp * t1(i,k) + phil(i,k) - hvap*qlx(i,k)
          endif
          tem2       = 1.+fv*max(q1(i,k,1),qmin)-qlx(i,k)
          thvx(i,k)  = theta(i,k) * tem2
          tvx(i,k)   = t1(i,k) * tem2
          qtx(i,k) = max(q1(i,k,1),qmin)+qlx(i,k)
          thlx(i,k)  = theta(i,k) - pix(i,k)*elocp*qlx(i,k)
          thlvx(i,k) = thlx(i,k) * (1. + fv * qtx(i,k))
          svx(i,k)   = cp * tvx(i,k)
          ptem1      = elocp * pix(i,k) * max(q1(i,k,1),qmin)
          thetae(i,k)= theta(i,k) +  ptem1
          gotvx(i,k) = g / tvx(i,k)
        enddo
      enddo
!
! The background vertical diffusivities in the inversion layers are limited 
!    to be less than or equal to xkzminv
!
      do k = 1,km1
        do i=1,im
          tem1 = (tvx(i,k+1)-tvx(i,k)) * rdzt(i,k)
          if(tem1 > 1.e-5) then
             xkzo(i,k)  = min(xkzo(i,k),xkzinv)
             xkzmo(i,k) = min(xkzmo(i,k),xkzinv)
          endif
        enddo
      enddo
!
! compute an empirical cloud fraction based on 
!    Xu & Randall's (1996,JAS) study
!
      do k = 1, km
        do i = 1, im
          plyr(i,k)   = 0.01 * prsl(i,k)   ! pa to mb (hpa)
!  --- ...  compute relative humidity
          es  = 0.01 * fpvs(t1(i,k))       ! fpvs in pa
          qs  = max(qmin, eps * es / (plyr(i,k) + epsm1*es))
          rhly(i,k) = max(0.0, min(1.0, max(qmin, q1(i,k,1))/qs))
          qstl(i,k) = qs
        enddo
      enddo
!
      do k = 1, km
        do i = 1, im
          cfly(i,k) = 0.
          clwt = 1.0e-6 * (plyr(i,k)*0.001)
          if (qlx(i,k) > clwt) then
            onemrh= max(1.e-10, 1.0-rhly(i,k))
            tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)
            tem1  = cql / tem1
            value = max(min( tem1*qlx(i,k), 50.0), 0.0)
            tem2  = sqrt(sqrt(rhly(i,k)))
            cfly(i,k) = min(max(tem2*(1.0-exp(-value)), 0.0), 1.0)
          endif
        enddo
      enddo
!
!  compute buoyancy modified by clouds
!
      do k = 1, km1
        do i = 1, im
          tem  = 0.5 * (svx(i,k) + svx(i,k+1))
          tem1 = 0.5 * (t1(i,k) + t1(i,k+1))
          tem2 = 0.5 * (qstl(i,k) + qstl(i,k+1))
          cfh  = min(cfly(i,k+1),0.5*(cfly(i,k)+cfly(i,k+1)))
          alp  = g / tem
          gamma = el2orc * tem2 / (tem1**2)
          epsi  = tem1 / elocp
          beta  = (1. + gamma*epsi*(1.+fv)) / (1. + gamma)
          chx   = cfh * alp * beta + (1. - cfh) * alp
          cqx   = cfh * alp * hvap * (beta - epsi)
          cqx   = cqx + (1. - cfh) * fv * g
          ptem1 = (slx(i,k+1)-slx(i,k))*rdzt(i,k)
          ptem2 = (qtx(i,k+1)-qtx(i,k))*rdzt(i,k)
          bf(i,k) = chx * ptem1 + cqx * ptem2
        enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k=1,km1
        do i=1,im
          dku(i,k)  = 0.
          dkt(i,k)  = 0.
          dkq(i,k)  = 0.
          cku(i,k)  = 0.
          ckt(i,k)  = 0.
          tem       = zi(i,k+1)-zi(i,k)
          radx(i,k) = tem*(swh(i,k)*xmu(i)+hlw(i,k))
        enddo
      enddo
!
      do i = 1,im
         sflux(i)  = heat(i) + evap(i)*fv*theta(i,1)
         if(.not.sfcflg(i) .or. sflux(i) <= 0.) pblflg(i)=.false.
      enddo
!
!  compute critical bulk richardson number 
!
      do i = 1,im
        if(pblflg(i)) then
!         thermal(i) = thvx(i,1)
          thermal(i) = thlvx(i,1)
          crb(i) = rbcr
        else
          thermal(i) = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
          tem = sqrt(u10m(i)**2+v10m(i)**2)
          tem = max(tem, 1.)
          robn = tem / (f0 * z0(i))
          tem1 = 1.e-7 * robn
          crb(i) = 0.16 * (tem1 ** (-0.18))
          crb(i) = max(min(crb(i), crbmax), crbmin)
        endif
      enddo
!
      do i=1,im
         dtdz1(i)  = dt2 / (zi(i,2)-zi(i,1))
      enddo
!
      do i=1,im
         ustar(i) = sqrt(stress(i))
      enddo
!
!  compute buoyancy (bf) and winshear square
!
      do k = 1, km1
      do i = 1, im
         rdz  = rdzt(i,k)
!        bf(i,k) = gotvx(i,k)*(thvx(i,k+1)-thvx(i,k))*rdz
         dw2  = (u1(i,k)-u1(i,k+1))**2
     &        + (v1(i,k)-v1(i,k+1))**2
         shr2(i,k) = max(dw2,dw2min)*rdz*rdz
      enddo
      enddo
!
! find pbl height based on bulk richardson number (mrf pbl scheme)
!   and also for diagnostic purpose
!
      do i=1,im
         flg(i) = .false.
         rbup(i) = rbsoil(i)
      enddo
!
      do k = 1, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          spdk2   = max((u1(i,k)**2+v1(i,k)**2),1.)
!         rbup(i) = (thvx(i,k)-thermal(i))*
!    &              (g*zl(i,k)/thvx(i,1))/spdk2
          rbup(i) = (thlvx(i,k)-thermal(i))*
     &              (g*zl(i,k)/thlvx(i,1))/spdk2
          kpblx(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
      do i = 1,im
        if(kpblx(i) > 1) then
          k = kpblx(i)
          if(rbdn(i) >= crb(i)) then
            rbint = 0.
          elseif(rbup(i) <= crb(i)) then
            rbint = 1.
          else
            rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
          endif
          hpblx(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
          if(hpblx(i) < zi(i,kpblx(i))) kpblx(i)=kpblx(i)-1
        else
          hpblx(i) = zl(i,1)
          kpblx(i) = 1
        endif
        hpbl(i) = hpblx(i)
        kpbl(i) = kpblx(i)
        if(kpbl(i) <= 1) pblflg(i)=.false.
      enddo
!
!     compute similarity parameters
!
      do i=1,im
         zol(i) = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
         if(sfcflg(i)) then
           zol(i) = min(zol(i),-zfmin)
         else
           zol(i) = max(zol(i),zfmin)
         endif
!
         zol1 = zol(i)*sfcfrac*hpbl(i)/zl(i,1)
         if(sfcflg(i)) then
           tem     = 1.0 / (1. - aphi16*zol1)
           phih(i) = sqrt(tem)
           phim(i) = sqrt(phih(i))
         else
           phim(i) = 1. + aphi5*zol1
           phih(i) = phim(i)
         endif
      enddo
!
      do i=1,im
        if(pblflg(i)) then
          if(zol(i) < zolcru) then
            pcnvflg(i) = .true.
          endif
          wst3(i) = gotvx(i,1)*sflux(i)*hpbl(i)
          wstar(i)= wst3(i)**h1
          ust3(i) = ustar(i)**3.
          wscale(i)=(ust3(i)+wfac*vk*wst3(i)*sfcfrac)**h1
          ptem = ustar(i)/aphi5
          wscale(i) = max(wscale(i),ptem)
        endif
      enddo
!
! compute a thermal excess
!
      do i = 1,im
         if(pcnvflg(i)) then
           hgamt(i) = heat(i)/wscale(i)
           hgamq(i) = evap(i)/wscale(i)
           vpert(i) = hgamt(i) + hgamq(i)*fv*theta(i,1)
           vpert(i) = max(vpert(i),0.)
           tem = min(cfac*vpert(i),gamcrt)
           thermal(i)= thermal(i) + tem
         endif
      enddo
!
!  enhance the pbl height by considering the thermal excess
!     (overshoot pbl top)
!
      do i=1,im
         flg(i)  = .true.
         if(pcnvflg(i)) then
           flg(i)  = .false.
           rbup(i) = rbsoil(i)
         endif
      enddo
      do k = 2, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          spdk2   = max((u1(i,k)**2+v1(i,k)**2),1.)
          rbup(i) = (thlvx(i,k)-thermal(i))*
     &              (g*zl(i,k)/thlvx(i,1))/spdk2
          kpbl(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
      do i = 1,im
        if(pcnvflg(i)) then
           k = kpbl(i)
           if(rbdn(i) >= crb(i)) then
             rbint = 0.
           elseif(rbup(i) <= crb(i)) then
             rbint = 1.
           else
             rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
           endif
           hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
           if(hpbl(i) < zi(i,kpbl(i))) then 
             kpbl(i) = kpbl(i) - 1
           endif
           if(kpbl(i) <= 1) then
              pcnvflg(i) = .false.
              pblflg(i) = .false.
           endif
        endif
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  look for stratocumulus
!
      do i=1,im
         flg(i)  = scuflg(i)
      enddo
      do k = 1, km1
        do i=1,im
          if(flg(i).and.zl(i,k) >= zstblmax) then
             lcld(i)=k
             flg(i)=.false.
          endif
      enddo
      enddo
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmscu,1,-1
      do i = 1, im
        if(flg(i) .and. k <= lcld(i)) then
          if(qlx(i,k) >= qlcr) then
             kcld(i)=k
             flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. kcld(i)==km1) scuflg(i)=.false.
      enddo
!
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmscu,1,-1
      do i = 1, im
        if(flg(i) .and. k <= kcld(i)) then
          if(qlx(i,k) >= qlcr) then
            if(radx(i,k) < radmin(i)) then
              radmin(i)=radx(i,k)
              krad(i)=k
            endif
          else
            flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. krad(i) <= 1) scuflg(i)=.false.
        if(scuflg(i) .and. radmin(i)>=0.) scuflg(i)=.false.
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute components for mass flux mixing by large thermals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
            tcko(i,k) = t1(i,k)
            ucko(i,k) = u1(i,k)
            vcko(i,k) = v1(i,k)
          endif
          if(scuflg(i)) then
            tcdo(i,k) = t1(i,k)
            ucdo(i,k) = u1(i,k)
            vcdo(i,k) = v1(i,k)
          endif
        enddo
      enddo
      do kk = 1, ntrac1
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
            qcko(i,k,kk) = q1(i,k,kk)
          endif
          if(scuflg(i)) then
            qcdo(i,k,kk) = q1(i,k,kk)
          endif
        enddo
      enddo
      enddo
!
      call mfpblt(im,ix,km,kmpbl,ntcw,ntrac1,dt2,
     &    pcnvflg,zl,zm,q1,t1,u1,v1,plyr,pix,thlx,thvx,
     &    gdx,hpbl,kpbl,vpert,buou,xmf,
     &    tcko,qcko,ucko,vcko,xlamue)
!
      call mfscu(im,ix,km,kmscu,ntcw,ntrac1,dt2,
     &    scuflg,zl,zm,q1,t1,u1,v1,plyr,pix,
     &    thlx,thvx,thlvx,gdx,thetae,radj,
     &    krad,mrad,radmin,buod,xmfd,
     &    tcdo,qcdo,ucdo,vcdo,xlamde)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   compute prandtl number and exchange coefficient varying with height
!
      do k = 1, kmpbl
        do i = 1, im
          if(k < kpbl(i)) then
            tem = phih(i)/phim(i)
            ptem = -3.*(max(zi(i,k+1)-sfcfrac*hpbl(i),0.))**2.
     &               /hpbl(i)**2.
            if(pcnvflg(i)) then
              prn(i,k) =  1. + (tem-1.)*exp(ptem)
            else
              prn(i,k) = tem
            endif
            prn(i,k) = min(prn(i,k),prmax)
            prn(i,k) = max(prn(i,k),prmin)
!
            ckz(i,k) = ck1 + (ck0-ck1)*exp(ptem)
            ckz(i,k) = min(ckz(i,k),ck0)
            ckz(i,k) = max(ckz(i,k),ck1)
            chz(i,k) = ch1 + (ch0-ch1)*exp(ptem)
            chz(i,k) = min(chz(i,k),ch0)
            chz(i,k) = max(chz(i,k),ch1)
          endif
        enddo
      enddo

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute an asymtotic mixing length
!
      do k = 1, km1
        do i = 1, im
          zlup = 0.0
          bsum = 0.0
          mlenflg = .true.
          do n = k, km1
            if(mlenflg) then
              dz = zl(i,n+1) - zl(i,n)
              ptem = gotvx(i,n)*(thvx(i,n+1)-thvx(i,k))*dz
!             ptem = gotvx(i,n)*(thlvx(i,n+1)-thlvx(i,k))*dz
              bsum = bsum + ptem
              zlup = zlup + dz
              if(bsum >= tke(i,k)) then
                if(ptem >= 0.) then
                  tem2 = max(ptem, zfmin)
                else
                  tem2 = min(ptem, -zfmin)
                endif
                ptem1 = (bsum - tke(i,k)) / tem2
                zlup = zlup - ptem1 * dz 
                zlup = max(zlup, 0.)
                mlenflg = .false.
              endif
            endif
          enddo
          zldn = 0.0
          bsum = 0.0
          mlenflg = .true.
          do n = k, 1, -1
            if(mlenflg) then
              if(n == 1) then
                dz = zl(i,1)
                tem1 = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
              else
                dz = zl(i,n) - zl(i,n-1)
                tem1 = thvx(i,n-1)
!               tem1 = thlvx(i,n-1)
              endif
              ptem = gotvx(i,n)*(thvx(i,k)-tem1)*dz
!             ptem = gotvx(i,n)*(thlvx(i,k)-tem1)*dz
              bsum = bsum + ptem
              zldn = zldn + dz
              if(bsum >= tke(i,k)) then
                if(ptem >= 0.) then
                  tem2 = max(ptem, zfmin)
                else
                  tem2 = min(ptem, -zfmin)
                endif
                ptem1 = (bsum - tke(i,k)) / tem2
                zldn = zldn - ptem1 * dz 
                zldn = max(zldn, 0.)
                mlenflg = .false.
              endif
            endif
          enddo
!
          tem = 0.5 * (zi(i,k+1)-zi(i,k))
          tem1 = min(tem, rlmn)
!
          ptem2 = min(zlup,zldn)
          rlam(i,k) = elmfac * ptem2
          rlam(i,k) = max(rlam(i,k), tem1)
          rlam(i,k) = min(rlam(i,k), rlmx)
!
          ptem2 = sqrt(zlup*zldn)
          ele(i,k) = elefac * ptem2
          ele(i,k) = max(ele(i,k), tem1)
          ele(i,k) = min(ele(i,k), elmx)
!
        enddo
      enddo
!
      do k = 1, km1
        do i = 1, im
          tem = vk * zl(i,k)
          if (zol(i) < 0.) then
            ptem = 1. - 100. * zol(i)
            ptem1 = ptem**0.2
            zk = tem * ptem1
          elseif (zol(i) >= 1.) then
            zk = tem / 3.7
          else
            ptem = 1. + 2.7 * zol(i)
            zk = tem / ptem
          endif 
          elm(i,k) = zk*rlam(i,k)/(rlam(i,k)+zk)
!
          dz = zi(i,k+1) - zi(i,k)
          tem = max(gdx(i),dz)
          elm(i,k) = min(elm(i,k), tem)
          ele(i,k) = min(ele(i,k), tem)
!
        enddo
      enddo
      do i = 1, im
        elm(i,km) = elm(i,km1)
        ele(i,km) = ele(i,km1)
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute eddy diffusivities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km1
        do i = 1, im
           tem = 0.5 * (elm(i,k) + elm(i,k+1))
           tem = tem * sqrt(tkeh(i,k))
           if(k < kpbl(i)) then
             if(pblflg(i)) then
               dku(i,k) = ckz(i,k) * tem
               dkt(i,k) = dku(i,k) / prn(i,k)
             else
               dkt(i,k) = chz(i,k) * tem
               dku(i,k) = dkt(i,k) * prn(i,k)
             endif
           else
              ri = max(bf(i,k)/shr2(i,k),rimin)
              if(ri < 0.) then ! unstable regime
                dku(i,k) = ck1 * tem
                dkt(i,k) = rchck * dku(i,k)
              else             ! stable regime
                dkt(i,k) = ch1 * tem
                prnum = 1.0 + 2.1*ri
                prnum = min(prnum,prmax)
                dku(i,k) = dkt(i,k) * prnum
              endif
           endif
!
           if(scuflg(i)) then
             if(k >= mrad(i) .and. k < krad(i)) then
                tem1 = ckz(i,k) * tem
                ptem1 = tem1 / prscu
                dku(i,k) = max(dku(i,k), tem1)
                dkt(i,k) = max(dkt(i,k), ptem1)
             endif
           endif
!
           dkq(i,k) = prtke * dkt(i,k)
!
           dkt(i,k) = min(dkt(i,k),dkmax)
           dkt(i,k) = max(dkt(i,k),xkzo(i,k))
           dkq(i,k) = min(dkq(i,k),dkmax)
           dkq(i,k) = max(dkq(i,k),xkzo(i,k))
           dku(i,k) = min(dku(i,k),dkmax)
           dku(i,k) = max(dku(i,k),xkzmo(i,k))
!
        enddo
      enddo
!
      do i = 1, im
        if(scuflg(i)) then
           k = krad(i)
           tem = bf(i,k) / gotvx(i,k)
           tem1 = max(tem, tdzmin)
           ptem = radj(i) / tem1
           dkt(i,k) = dkt(i,k) + ptem
           dku(i,k) = dku(i,k) + ptem
           dkq(i,k) = dkq(i,k) + ptem
        endif
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute buoyancy and shear productions of tke 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km1
        do i = 1, im
          if (k == 1) then
            tem = -dkt(i,1) * bf(i,1)
!           if(pcnvflg(i)) then
!             ptem1 = xmf(i,1) * buou(i,1)
!           else
              ptem1 = 0.
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem2 = xmfd(i,1) * buod(i,1)
            else
              ptem2 = 0.
            endif
            tem = tem + ptem1 + ptem2
            buop = 0.5 * (gotvx(i,1) * sflux(i) + tem)
!
            tem1 = dku(i,1) * shr2(i,1)
!
            tem = (u1(i,2)-u1(i,1))*rdzt(i,1)
!           if(pcnvflg(i)) then
!             ptem = xmf(i,1) * tem
!             ptem1 = 0.5 * ptem * (u1(i,2)-ucko(i,2))
!           else
              ptem1 = 0.
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem = ucdo(i,1)+ucdo(i,2)-u1(i,1)-u1(i,2)
              ptem = 0.5 * tem * xmfd(i,1) * ptem
            else
              ptem = 0.
            endif
            ptem1 = ptem1 + ptem
!
            tem = (v1(i,2)-v1(i,1))*rdzt(i,1)
!           if(pcnvflg(i)) then
!             ptem = xmf(i,1) * tem
!             ptem2 = 0.5 * ptem * (v1(i,2)-vcko(i,2))
!           else
              ptem2 = 0.
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem = vcdo(i,1)+vcdo(i,2)-v1(i,1)-v1(i,2)
              ptem = 0.5 * tem * xmfd(i,1) * ptem
            else
              ptem = 0.
            endif
            ptem2 = ptem2 + ptem
!
!           tem2 = stress(i)*spd1(i)/zl(i,1)
            tem2 = stress(i)*ustar(i)*phim(i)/(vk*zl(i,1))
            shrp = 0.5 * (tem1 + ptem1 + ptem2 + tem2)
          else
            tem1 = -dkt(i,k-1) * bf(i,k-1)
            tem2 = -dkt(i,k) * bf(i,k)
            tem  = 0.5 * (tem1 + tem2)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = 0.5 * (xmf(i,k-1) + xmf(i,k))
              ptem1 = ptem * buou(i,k)
            else
              ptem1 = 0.
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = 0.5 * (xmfd(i,k-1) + xmfd(i,k))
                ptem2 = ptem0 * buod(i,k)
              else
                ptem2 = 0.
              endif
            else
              ptem2 = 0.
            endif
            buop = tem + ptem1 + ptem2
!
            tem1 = dku(i,k-1) * shr2(i,k-1)
            tem2 = dku(i,k) * shr2(i,k)
            tem  = 0.5 * (tem1 + tem2)
            tem1 = (u1(i,k+1)-u1(i,k))*rdzt(i,k)
            tem2 = (u1(i,k)-u1(i,k-1))*rdzt(i,k-1)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = xmf(i,k) * tem1 + xmf(i,k-1) * tem2
              ptem1 = 0.5 * ptem * (u1(i,k)-ucko(i,k))
            else
              ptem1 = 0.
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = xmfd(i,k) * tem1 + xmfd(i,k-1) * tem2
                ptem2 = 0.5 * ptem0 * (ucdo(i,k)-u1(i,k))
              else
                ptem2 = 0.
              endif
            else
              ptem2 = 0.
            endif
            shrp = tem + ptem1 + ptem2
            tem1 = (v1(i,k+1)-v1(i,k))*rdzt(i,k)
            tem2 = (v1(i,k)-v1(i,k-1))*rdzt(i,k-1)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = xmf(i,k) * tem1 + xmf(i,k-1) * tem2
              ptem1 = 0.5 * ptem * (v1(i,k)-vcko(i,k))
            else
              ptem1 = 0.
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = xmfd(i,k) * tem1 + xmfd(i,k-1) * tem2
                ptem2 = 0.5 * ptem0 * (vcdo(i,k)-v1(i,k))
              else
                ptem2 = 0.
              endif
            else
              ptem2 = 0.
            endif
            shrp = shrp + ptem1 + ptem2
          endif
          prod(i,k) = buop + shrp
        enddo
      enddo
!
!----------------------------------------------------------------------
!     first predict tke due to tke production & dissipation(diss) 
!
      do k = 1,km1
        do i=1,im
           rle(i,k) = ce0 / ele(i,k)
        enddo
      enddo
      kk = max(nint(dt2/cdtn), 1)
      dtn = dt2 / float(kk)
      do n = 1, kk
      do k = 1,km1
        do i=1,im
           tem = sqrt(tke(i,k))
           diss(i,k) = rle(i,k) * tke(i,k) * tem
           tem1 = prod(i,k) + tke(i,k) / dtn
           diss(i,k)=max(min(diss(i,k), tem1), 0.)
           tke(i,k) = tke(i,k) + dtn * (prod(i,k)-diss(i,k))
           tke(i,k) = max(tke(i,k), tkmin)
        enddo
      enddo
      enddo
!
!     compute updraft & downdraft properties for tke
!
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
            qcko(i,k,ntke) = tke(i,k)
          endif
          if(scuflg(i)) then
            qcdo(i,k,ntke) = tke(i,k)
          endif
        enddo
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if (pcnvflg(i) .and. k <= kpbl(i)) then
             dz   = zl(i,k) - zl(i,k-1)
             tem  = 0.5 * xlamue(i,k-1) * dz
             factor = 1. + tem
             qcko(i,k,ntke)=((1.-tem)*qcko(i,k-1,ntke)+tem*
     &                (tke(i,k)+tke(i,k-1)))/factor
          endif
        enddo
      enddo
      do k = kmscu, 1, -1
        do i = 1, im
          if (scuflg(i) .and. k < krad(i)) then
            if(k >= mrad(i)) then
              dz = zl(i,k+1) - zl(i,k)
              tem  = 0.5 * xlamde(i,k) * dz
              factor = 1. + tem
              qcdo(i,k,ntke)=((1.-tem)*qcdo(i,k+1,ntke)+tem*
     &                 (tke(i,k)+tke(i,k+1)))/factor
            endif
          endif
        enddo
      enddo
!
!----------------------------------------------------------------------
!     compute tridiagonal matrix elements for turbulent kinetic energy
!
      do i=1,im
         ad(i,1) = 1.0
         f1(i,1) = tke(i,1)
      enddo
!
      do k = 1,km1
        do i=1,im
          dtodsd  = dt2/del(i,k)
          dtodsu  = dt2/del(i,k+1)
          dsig    = prsl(i,k)-prsl(i,k+1)
          rdz     = rdzt(i,k)
          tem1    = dsig * dkq(i,k) * rdz
          dsdz2   = tem1 * rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5 * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             tem       = tke(i,k) + tke(i,k+1)
             ptem      = qcko(i,k,ntke) + qcko(i,k+1,ntke)
             f1(i,k)   = f1(i,k)-(ptem-tem)*ptem1
             f1(i,k+1) = tke(i,k+1)+(ptem-tem)*ptem2
          else
             f1(i,k+1) = tke(i,k+1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5 * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              tem       = tke(i,k) + tke(i,k+1)
              ptem      = qcdo(i,k,ntke) + qcdo(i,k+1,ntke)
              f1(i,k)   = f1(i,k) + (ptem - tem) * ptem1
              f1(i,k+1) = f1(i,k+1) - (ptem - tem) * ptem2
            endif
          endif
!
        enddo
      enddo
!
!     solve tridiagonal problem for tke
!
      call tridit(im,km,1,al,ad,au,f1,au,f1)
!
!     recover tendency of tke
!
      do k = 1,km
         do i = 1,im
!           f1(i,k) = max(f1(i,k), tkmin)
            qtend = (f1(i,k)-q1(i,k,ntke))*rdt
            rtg(i,k,ntke) = rtg(i,k,ntke)+qtend
         enddo
      enddo
!
!     compute tridiagonal matrix elements for heat and moisture
!
      do i=1,im
         ad(i,1) = 1.
         f1(i,1) = t1(i,1)   + dtdz1(i) * heat(i)
         f2(i,1) = q1(i,1,1) + dtdz1(i) * evap(i)
      enddo
      if(ntrac1 >= 2) then
        do kk = 2, ntrac1
          is = (kk-1) * km
          do i = 1, im
            f2(i,1+is) = q1(i,1,kk)
          enddo
        enddo
      endif
!
      do k = 1,km1
        do i = 1,im
          dtodsd  = dt2/del(i,k)
          dtodsu  = dt2/del(i,k+1)
          dsig    = prsl(i,k)-prsl(i,k+1)
          rdz     = rdzt(i,k)
          tem1    = dsig * dkt(i,k) * rdz
          dsdzt   = tem1 * gocp
          dsdz2   = tem1 * rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5 * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             tem       = t1(i,k) + t1(i,k+1)
             ptem      = tcko(i,k) + tcko(i,k+1)
             f1(i,k)   = f1(i,k)+dtodsd*dsdzt-(ptem-tem)*ptem1
             f1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt+(ptem-tem)*ptem2
             tem       = q1(i,k,1) + q1(i,k+1,1)
             ptem      = qcko(i,k,1) + qcko(i,k+1,1)
             f2(i,k)   = f2(i,k) - (ptem - tem) * ptem1
             f2(i,k+1) = q1(i,k+1,1) + (ptem - tem) * ptem2
          else
             f1(i,k)   = f1(i,k)+dtodsd*dsdzt
             f1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
             f2(i,k+1) = q1(i,k+1,1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5 * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              ptem      = tcdo(i,k) + tcdo(i,k+1)
              tem       = t1(i,k) + t1(i,k+1)
              f1(i,k)   = f1(i,k) + (ptem - tem) * ptem1
              f1(i,k+1) = f1(i,k+1) - (ptem - tem) * ptem2
              tem       = q1(i,k,1) + q1(i,k+1,1)
              ptem      = qcdo(i,k,1) + qcdo(i,k+1,1)
              f2(i,k)   = f2(i,k) + (ptem - tem) * ptem1
              f2(i,k+1) = f2(i,k+1) - (ptem - tem) * ptem2
            endif
          endif
        enddo
      enddo
!
      if(ntrac1 >= 2) then
        do kk = 2, ntrac1
          is = (kk-1) * km
          do k = 1, km1
            do i = 1, im
              if(pcnvflg(i) .and. k < kpbl(i)) then
                dtodsd = dt2/del(i,k)
                dtodsu = dt2/del(i,k+1)
                dsig  = prsl(i,k)-prsl(i,k+1)
                tem   = dsig * rdzt(i,k)
                ptem  = 0.5 * tem * xmf(i,k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem1  = qcko(i,k,kk) + qcko(i,k+1,kk)
                tem2  = q1(i,k,kk) + q1(i,k+1,kk)
                f2(i,k+is) = f2(i,k+is) - (tem1 - tem2) * ptem1
                f2(i,k+1+is)= q1(i,k+1,kk) + (tem1 - tem2) * ptem2
              else
                f2(i,k+1+is) = q1(i,k+1,kk)
              endif
!
              if(scuflg(i)) then
                if(k >= mrad(i) .and. k < krad(i)) then
                  dtodsd = dt2/del(i,k)
                  dtodsu = dt2/del(i,k+1)
                  dsig  = prsl(i,k)-prsl(i,k+1)
                  tem   = dsig * rdzt(i,k)
                  ptem  = 0.5 * tem * xmfd(i,k)
                  ptem1 = dtodsd * ptem
                  ptem2 = dtodsu * ptem
                  tem1  = qcdo(i,k,kk) + qcdo(i,k+1,kk)
                  tem2  = q1(i,k,kk) + q1(i,k+1,kk)
                  f2(i,k+is)  = f2(i,k+is) + (tem1 - tem2) * ptem1
                  f2(i,k+1+is)= f2(i,k+1+is) - (tem1 - tem2) * ptem2
                endif
              endif
!
            enddo
          enddo
        enddo
      endif
!
!     solve tridiagonal problem for heat and moisture
!
      call tridin(im,km,ntrac1,al,ad,au,f1,f2,au,f1,f2)
!
!     recover tendencies of heat and moisture
!
      do  k = 1,km
         do i = 1,im
            ttend      = (f1(i,k)-t1(i,k))*rdt
            qtend      = (f2(i,k)-q1(i,k,1))*rdt
            tdt(i,k)   = tdt(i,k)+ttend
            rtg(i,k,1) = rtg(i,k,1)+qtend
            dtsfc(i)   = dtsfc(i)+cont*del(i,k)*ttend
            dqsfc(i)   = dqsfc(i)+conq*del(i,k)*qtend
         enddo
      enddo
!
      if(ntrac1 >= 2) then
        do kk = 2, ntrac1
          is = (kk-1) * km
          do k = 1, km
            do i = 1, im
              qtend = (f2(i,k+is)-q1(i,k,kk))*rdt
              rtg(i,k,kk) = rtg(i,k,kk)+qtend
            enddo
          enddo
        enddo
      endif
!
!     add tke dissipative heating to temperature tendency
!
      if(dspheat) then
      do k = 1,km1
        do i = 1,im
!         tem = min(diss(i,k), dspmax)
!         ttend = tem / cp
          ttend = diss(i,k) / cp
          tdt(i,k) = tdt(i,k) + dspfac * ttend
        enddo
      enddo
      endif
!
!     compute tridiagonal matrix elements for momentum
!
      do i=1,im
         ad(i,1) = 1.0 + dtdz1(i) * stress(i) / spd1(i)
         f1(i,1) = u1(i,1)
         f2(i,1) = v1(i,1)
      enddo
!
      do k = 1,km1
        do i=1,im
          dtodsd  = dt2/del(i,k)
          dtodsu  = dt2/del(i,k+1)
          dsig    = prsl(i,k)-prsl(i,k+1)
          rdz     = rdzt(i,k)
          tem1    = dsig * dku(i,k) * rdz
          dsdz2   = tem1*rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5 * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             tem       = u1(i,k) + u1(i,k+1)
             ptem      = ucko(i,k) + ucko(i,k+1)
             f1(i,k)   = f1(i,k) - (ptem - tem) * ptem1
             f1(i,k+1) = u1(i,k+1) + (ptem - tem) * ptem2
             tem       = v1(i,k) + v1(i,k+1)
             ptem      = vcko(i,k) + vcko(i,k+1)
             f2(i,k)   = f2(i,k) - (ptem - tem) * ptem1
             f2(i,k+1) = v1(i,k+1) + (ptem - tem) * ptem2
          else
             f1(i,k+1) = u1(i,k+1)
             f2(i,k+1) = v1(i,k+1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5 * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              tem       = u1(i,k) + u1(i,k+1)
              ptem      = ucdo(i,k) + ucdo(i,k+1)
              f1(i,k)   = f1(i,k) + (ptem - tem) *ptem1
              f1(i,k+1) = f1(i,k+1) - (ptem - tem) *ptem2
              tem       = v1(i,k) + v1(i,k+1)
              ptem      = vcdo(i,k) + vcdo(i,k+1)
              f2(i,k)   = f2(i,k) + (ptem - tem) * ptem1
              f2(i,k+1) = f2(i,k+1) - (ptem - tem) * ptem2
            endif
          endif
!
        enddo
      enddo
!
!     solve tridiagonal problem for momentum
!
      call tridi2(im,km,al,ad,au,f1,f2,au,f1,f2)
!
!     recover tendencies of momentum
!
      do k = 1,km
         do i = 1,im
            utend = (f1(i,k)-u1(i,k))*rdt
            vtend = (f2(i,k)-v1(i,k))*rdt
            du(i,k)  = du(i,k)+utend
            dv(i,k)  = dv(i,k)+vtend
            dusfc(i) = dusfc(i)+conw*del(i,k)*utend
            dvsfc(i) = dvsfc(i)+conw*del(i,k)*vtend
         enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  pbl height for diagnostic purpose
!
      do i = 1, im
         hpbl(i) = hpblx(i)
         kpbl(i) = kpblx(i)
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end subroutine satmedmfvdif_run
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine tridit(l,n,nt,cl,cm,cu,rt,au,at)
!-----------------------------------------------------------------------
!!
      use machine     , only : kind_phys
      implicit none
      integer             is,k,kk,n,nt,l,i
      real(kind=kind_phys) fk(l)
!!
      real(kind=kind_phys) cl(l,2:n), cm(l,n), cu(l,n-1),               &
     &                     rt(l,n*nt),                                  &
     &                     au(l,n-1), at(l,n*nt),                       &
     &                     fkk(l,2:n-1)
!-----------------------------------------------------------------------
      do i=1,l
        fk(i)   = 1./cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          at(i,1+is) = fk(i) * rt(i,1+is)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fkk(i,k) = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=2,n-1
          do i=1,l
            at(i,k+is) = fkk(i,k)*(rt(i,k+is)-cl(i,k)*at(i,k+is-1))
          enddo
        enddo
      enddo
      do i=1,l
        fk(i)   = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          at(i,n+is) = fk(i)*(rt(i,n+is)-cl(i,n)*at(i,n+is-1))
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=n-1,1,-1
          do i=1,l
            at(i,k+is) = at(i,k+is) - au(i,k)*at(i,k+is+1)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end subroutine tridit

      end module satmedmfvdif

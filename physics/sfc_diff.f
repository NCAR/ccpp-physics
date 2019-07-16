!>  \file sfc_diff.f
!! This file contains the surface roughness length formulation based on
!! the surface sublayer scheme in
!! Zeng and Dickinson (1998) \cite zeng_and_dickinson_1998.

!> This module contains the CCPP-compliant GFS surface layer scheme.
      module sfc_diff

      use machine , only : kind_phys

      implicit none

      public :: sfc_diff_init, sfc_diff_run, sfc_diff_finalize

      private

      real (kind=kind_phys), parameter :: ca=.4  ! ca - von karman constant

      contains

      subroutine sfc_diff_init
      end subroutine sfc_diff_init

      subroutine sfc_diff_finalize
      end subroutine sfc_diff_finalize

!> \defgroup GFS_diff_main GFS Surface Layer Scheme Module
!> @{
!> \brief This subroutine calculates surface roughness length.
!!
!! This subroutine includes the surface roughness length formulation
!! based on the surface sublayer scheme in
!! Zeng and Dickinson (1998) \cite zeng_and_dickinson_1998.
!> \section arg_table_sfc_diff_run Argument Table
!! | local_name     | standard_name                                                                | long_name                                                        | units      | rank | type      |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                           | count      |    0 | integer   |           | in     | F        |
!! | rvrdm1         | ratio_of_vapor_to_dry_air_gas_constants_minus_one                            | (rv/rd) - 1 (rv = ideal gas constant for water vapor)            | none       |    0 | real      | kind_phys | in     | F        |
!! | eps            | ratio_of_dry_air_to_water_vapor_gas_constants                                | rd/rv                                                            | none       |    0 | real      | kind_phys | in     | F        |
!! | epsm1          | ratio_of_dry_air_to_water_vapor_gas_constants_minus_one                      | (rd/rv) - 1                                                      | none       |    0 | real      | kind_phys | in     | F        |
!! | grav           | gravitational_acceleration                                                   | gravitational acceleration                                       | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | ps             | surface_air_pressure                                                         | surface pressure                                                 | Pa         |    1 | real      | kind_phys | in     | F        |
!! | u1             | x_wind_at_lowest_model_layer                                                 | x component of 1st model layer wind                              | m s-1      |    1 | real      | kind_phys | in     | F        |
!! | v1             | y_wind_at_lowest_model_layer                                                 | y component of 1st model layer wind                              | m s-1      |    1 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature_at_lowest_model_layer                                        | 1st model layer air temperature                                  | K          |    1 | real      | kind_phys | in     | F        |
!! | q1             | water_vapor_specific_humidity_at_lowest_model_layer                          | 1st model layer specific humidity                                | kg kg-1    |    1 | real      | kind_phys | in     | F        |
!! | z1             | height_above_ground_at_lowest_model_layer                                    | height above ground at 1st model layer                           | m          |    1 | real      | kind_phys | in     | F        |
!! | prsl1          | air_pressure_at_lowest_model_layer                                           | Model layer 1 mean pressure                                      | Pa         |    1 | real      | kind_phys | in     | F        |
!! | prslki         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer      | ratio      |    1 | real      | kind_phys | in     | F        |
!! | ddvel          | surface_wind_enhancement_due_to_convection                                   | surface wind enhancement due to convection                       | m s-1      |    1 | real      | kind_phys | in     | F        |
!! | sigmaf         | bounded_vegetation_area_fraction                                             | areal fractional cover of green vegetation bounded on the bottom | frac       |    1 | real      | kind_phys | in     | F        |
!! | vegtype        | vegetation_type_classification                                               | vegetation type at each grid cell                                | index      |    1 | integer   |           | in     | F        |
!! | shdmax         | maximum_vegetation_area_fraction                                             | max fractnl cover of green veg                                   | frac       |    1 | real      | kind_phys | in     | F        |
!! | ivegsrc        | vegetation_type_dataset_choice                                               | land use dataset choice                                          | index      |    0 | integer   |           | in     | F        |
!! | z0pert         | perturbation_of_momentum_roughness_length                                    | perturbation of momentum roughness length                        | frac       |    1 | real      | kind_phys | in     | F        |
!! | ztpert         | perturbation_of_heat_to_momentum_roughness_length_ratio                      | perturbation of heat to momentum roughness length ratio          | frac       |    1 | real      | kind_phys | in     | F        |
!! | flag_iter      | flag_for_iteration                                                           | flag for iteration                                               | flag       |    1 | logical   |           | in     | F        |
!! | redrag         | flag_for_reduced_drag_coefficient_over_sea                                   | flag for reduced drag coefficient over sea                       | flag       |    0 | logical   |           | in     | F        |
!! | wet            | flag_nonzero_wet_surface_fraction                                            | flag indicating presence of some ocean or lake surface area fraction | flag   |    1 | logical   |           | in     | F        |
!! | dry            | flag_nonzero_land_surface_fraction                                           | flag indicating presence of some land surface area fraction      | flag       |    1 | logical   |           | in     | F        |
!! | icy            | flag_nonzero_sea_ice_surface_fraction                                        | flag indicating presence of some sea ice surface area fraction   | flag       |    1 | logical   |           | in     | F        |
!! | tskin_ocn      | surface_skin_temperature_over_ocean_interstitial                             | surface skin temperature over ocean (temporary use as interstitial) | K       |    1 | real      | kind_phys | in     | F        |
!! | tskin_lnd      | surface_skin_temperature_over_land_interstitial                              | surface skin temperature over land  (temporary use as interstitial) | K       |    1 | real      | kind_phys | in     | F        |
!! | tskin_ice      | surface_skin_temperature_over_ice_interstitial                               | surface skin temperature over ice   (temporary use as interstitial) | K       |    1 | real      | kind_phys | in     | F        |
!! | tsurf_ocn      | surface_skin_temperature_after_iteration_over_ocean                          | surface skin temperature after iteration over ocean              | K          |    1 | real      | kind_phys | in     | F        |
!! | tsurf_lnd      | surface_skin_temperature_after_iteration_over_land                           | surface skin temperature after iteration over land               | K          |    1 | real      | kind_phys | in     | F        |
!! | tsurf_ice      | surface_skin_temperature_after_iteration_over_ice                            | surface skin temperature after iteration over ice                | K          |    1 | real      | kind_phys | in     | F        |
!! | snwdph_ocn     | surface_snow_thickness_water_equivalent_over_ocean                           | water equivalent snow depth over ocean                           | mm         |    1 | real      | kind_phys | in     | F        |
!! | snwdph_lnd     | surface_snow_thickness_water_equivalent_over_land                            | water equivalent snow depth over land                            | mm         |    1 | real      | kind_phys | in     | F        |
!! | snwdph_ice     | surface_snow_thickness_water_equivalent_over_ice                             | water equivalent snow depth over ice                             | mm         |    1 | real      | kind_phys | in     | F        |
!! | z0rl_ocn       | surface_roughness_length_over_ocean_interstitial                             | surface roughness length over ocean (temporary use as interstitial) | cm      |    1 | real      | kind_phys | inout  | F        |
!! | z0rl_lnd       | surface_roughness_length_over_land_interstitial                              | surface roughness length over land  (temporary use as interstitial) | cm      |    1 | real      | kind_phys | inout  | F        |
!! | z0rl_ice       | surface_roughness_length_over_ice_interstitial                               | surface roughness length over ice   (temporary use as interstitial) | cm      |    1 | real      | kind_phys | inout  | F        |
!! | ustar_ocn      | surface_friction_velocity_over_ocean                                         | surface friction velocity over ocean                             | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | ustar_lnd      | surface_friction_velocity_over_land                                          | surface friction velocity over land                              | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | ustar_ice      | surface_friction_velocity_over_ice                                           | surface friction velocity over ice                               | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | cm_ocn         | surface_drag_coefficient_for_momentum_in_air_over_ocean                      | surface exchange coeff for momentum over ocean                   | none       |    1 | real      | kind_phys | inout  | F        |
!! | cm_lnd         | surface_drag_coefficient_for_momentum_in_air_over_land                       | surface exchange coeff for momentum over land                    | none       |    1 | real      | kind_phys | inout  | F        |
!! | cm_ice         | surface_drag_coefficient_for_momentum_in_air_over_ice                        | surface exchange coeff for momentum over ice                     | none       |    1 | real      | kind_phys | inout  | F        |
!! | ch_ocn         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ocean             | surface exchange coeff heat & moisture over ocean                | none       |    1 | real      | kind_phys | inout  | F        |
!! | ch_lnd         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_land              | surface exchange coeff heat & moisture over land                 | none       |    1 | real      | kind_phys | inout  | F        |
!! | ch_ice         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ice               | surface exchange coeff heat & moisture over ice                  | none       |    1 | real      | kind_phys | inout  | F        |
!! | rb_ocn         | bulk_richardson_number_at_lowest_model_level_over_ocean                      | bulk Richardson number at the surface over ocean                 | none       |    1 | real      | kind_phys | inout  | F        |
!! | rb_lnd         | bulk_richardson_number_at_lowest_model_level_over_land                       | bulk Richardson number at the surface over land                  | none       |    1 | real      | kind_phys | inout  | F        |
!! | rb_ice         | bulk_richardson_number_at_lowest_model_level_over_ice                        | bulk Richardson number at the surface over ice                   | none       |    1 | real      | kind_phys | inout  | F        |
!! | stress_ocn     | surface_wind_stress_over_ocean                                               | surface wind stress over ocean                                   | m2 s-2     |    1 | real      | kind_phys | inout  | F        |
!! | stress_lnd     | surface_wind_stress_over_land                                                | surface wind stress over land                                    | m2 s-2     |    1 | real      | kind_phys | inout  | F        |
!! | stress_ice     | surface_wind_stress_over_ice                                                 | surface wind stress over ice                                     | m2 s-2     |    1 | real      | kind_phys | inout  | F        |
!! | fm_ocn         | Monin-Obukhov_similarity_function_for_momentum_over_ocean                    | Monin-Obukhov similarity function for momentum over ocean        | none       |    1 | real      | kind_phys | inout  | F        |
!! | fm_lnd         | Monin-Obukhov_similarity_function_for_momentum_over_land                     | Monin-Obukhov similarity function for momentum over land         | none       |    1 | real      | kind_phys | inout  | F        |
!! | fm_ice         | Monin-Obukhov_similarity_function_for_momentum_over_ice                      | Monin-Obukhov similarity function for momentum over ice          | none       |    1 | real      | kind_phys | inout  | F        |
!! | fh_ocn         | Monin-Obukhov_similarity_function_for_heat_over_ocean                        | Monin-Obukhov similarity function for heat over ocean            | none       |    1 | real      | kind_phys | inout  | F        |
!! | fh_lnd         | Monin-Obukhov_similarity_function_for_heat_over_land                         | Monin-Obukhov similarity function for heat over land             | none       |    1 | real      | kind_phys | inout  | F        |
!! | fh_ice         | Monin-Obukhov_similarity_function_for_heat_over_ice                          | Monin-Obukhov similarity function for heat over ice              | none       |    1 | real      | kind_phys | inout  | F        |
!! | fm10_ocn       | Monin-Obukhov_similarity_function_for_momentum_at_10m_over_ocean             | Monin-Obukhov similarity parameter for momentum at 10m over ocean | none      |    1 | real      | kind_phys | inout  | F        |
!! | fm10_lnd       | Monin-Obukhov_similarity_function_for_momentum_at_10m_over_land              | Monin-Obukhov similarity parameter for momentum at 10m over land | none       |    1 | real      | kind_phys | inout  | F        |
!! | fm10_ice       | Monin-Obukhov_similarity_function_for_momentum_at_10m_over_ice               | Monin-Obukhov similarity parameter for momentum at 10m over ice  | none       |    1 | real      | kind_phys | inout  | F        |
!! | fh2_ocn        | Monin-Obukhov_similarity_function_for_heat_at_2m_over_ocean                  | Monin-Obukhov similarity parameter for heat at 2m over ocean     | none       |    1 | real      | kind_phys | inout  | F        |
!! | fh2_lnd        | Monin-Obukhov_similarity_function_for_heat_at_2m_over_land                   | Monin-Obukhov similarity parameter for heat at 2m over land      | none       |    1 | real      | kind_phys | inout  | F        |
!! | fh2_ice        | Monin-Obukhov_similarity_function_for_heat_at_2m_over_ice                    | Monin-Obukhov similarity parameter for heat at 2m over ice       | none       |    1 | real      | kind_phys | inout  | F        |
!! | wind           | wind_speed_at_lowest_model_layer                                             | wind speed at lowest model level                                 | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                           | error message for error handling in CCPP                         | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                              | error flag for error handling in CCPP                            | flag       |    0 | integer   |           | out    | F        |
!!
!>  \section general_diff GFS Surface Layer Scheme General Algorithm
!! - Calculate the thermal roughness length formulation over the ocean (see eq. (25) and (26)
!!  in Zeng et al. (1998) \cite zeng_et_al_1998).
!! - Calculate Zeng's momentum roughness length formulation over land and sea ice.
!! - Calculate the new vegetation-dependent formulation of thermal roughness length
!! (Zheng et al.(2009) \cite zheng_et_al_2009).
!! Zheng et al. (2009) \cite zheng_et_al_2009 proposed a new formulation on
!! \f$ln(Z_{0m}^,/Z_{0t})\f$ as follows:
!! \f[
!!  ln(Z_{0m}^,/Z_{0t})=(1-GVF)^2C_{zil}k(u*Z_{0g}/\nu)^{0.5}
!! \f]
!! where \f$Z_{0m}^,\f$ is the effective momentum roughness length
!! computed in the following equation for each grid, \f$Z_{0t}\f$
!! is the roughness lenghth for heat, \f$C_{zil}\f$ is a coefficient
!! (taken as 0.8), k is the Von Karman constant (0.4),
!! \f$\nu=1.5\times10^{-5}m^{2}s^{-1}\f$ is the molecular viscosity,
!! \f$u*\f$ is the friction velocity, and \f$Z_{0g}\f$ is the bare
!! soil roughness length for momentum (taken as 0.01).
!! \n In order to consider the convergence of \f$Z_{0m}\f$ between
!! fully vegetated and bare soil, the effective \f$Z_{0m}^,\f$ is
!! computed:
!! \f[
!!  \ln(Z_{0m}^,)=(1-GVF)^{2}\ln(Z_{0g})+\left[1-(1-GVF)^{2}\right]\ln(Z_{0m})
!!\f]
!! - Calculate the exchange coefficients:\f$cm\f$, \f$ch\f$, and \f$stress\f$ as inputs of other \a sfc schemes.
!!
      subroutine sfc_diff_run (im,rvrdm1,eps,epsm1,grav,                &  !intent(in)
     &                    ps,u1,v1,t1,q1,z1,                            &  !intent(in)
     &                    prsl1,prslki,ddvel,                           &  !intent(in)
     &                    sigmaf,vegtype,shdmax,ivegsrc,                &  !intent(in)
     &                    z0pert,ztpert,                                &  ! mg, sfc-perts !intent(in)
     &                    flag_iter,redrag,                             &  !intent(in)
     &                    wet,dry,icy,                                  &  !intent(in)
     &                    tskin_ocn, tskin_lnd, tskin_ice,              &  !intent(in)
     &                    tsurf_ocn, tsurf_lnd, tsurf_ice,              &  !intent(in)
     &                   snwdph_ocn,snwdph_lnd,snwdph_ice,              &  !intent(in)
     &                     z0rl_ocn,  z0rl_lnd,  z0rl_ice,              &  !intent(inout)
     &                    ustar_ocn, ustar_lnd, ustar_ice,              &  !intent(inout)
     &                       cm_ocn,    cm_lnd,    cm_ice,              &  !intent(inout)
     &                       ch_ocn,    ch_lnd,    ch_ice,              &  !intent(inout)
     &                       rb_ocn,    rb_lnd,    rb_ice,              &  !intent(inout)
     &                   stress_ocn,stress_lnd,stress_ice,              &  !intent(inout)
     &                       fm_ocn,    fm_lnd,    fm_ice,              &  !intent(inout)
     &                       fh_ocn,    fh_lnd,    fh_ice,              &  !intent(inout)
     &                     fm10_ocn,  fm10_lnd,  fm10_ice,              &  !intent(inout)
     &                      fh2_ocn,   fh2_lnd,   fh2_ice,              &  !intent(inout)
     &                    wind                           ,              &  !intent(inout)
     &                    errmsg, errflg)                                  !intent(out)
!
      use funcphys, only : fpvs

      implicit none
!
      integer, intent(in) :: im, ivegsrc
      integer, dimension(im), intent(in) :: vegtype 

      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
      logical, dimension(im), intent(in) :: flag_iter, wet, dry, icy ! added by s.lu

      real(kind=kind_phys), intent(in) :: rvrdm1, eps, epsm1, grav
      real(kind=kind_phys), dimension(im), intent(in)    ::             &
     &                    ps,u1,v1,t1,q1,z1,prsl1,prslki,ddvel,         &
     &                    sigmaf,shdmax,                                &
     &                    z0pert,ztpert ! mg, sfc-perts
      real(kind=kind_phys), dimension(im), intent(in)    ::             &
     &                    tskin_ocn, tskin_lnd, tskin_ice,              &
     &                    tsurf_ocn, tsurf_lnd, tsurf_ice,              &
     &                   snwdph_ocn,snwdph_lnd,snwdph_ice

      real(kind=kind_phys), dimension(im), intent(inout) ::             &
     &                     z0rl_ocn,  z0rl_lnd,  z0rl_ice,              &
     &                    ustar_ocn, ustar_lnd, ustar_ice,              &
     &                       cm_ocn,    cm_lnd,    cm_ice,              &
     &                       ch_ocn,    ch_lnd,    ch_ice,              &
     &                       rb_ocn,    rb_lnd,    rb_ice,              &
     &                   stress_ocn,stress_lnd,stress_ice,              &
     &                       fm_ocn,    fm_lnd,    fm_ice,              &
     &                       fh_ocn,    fh_lnd,    fh_ice,              &
     &                     fm10_ocn,  fm10_lnd,  fm10_ice,              &
     &                      fh2_ocn,   fh2_lnd,   fh2_ice,              &
     &                      wind
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!     locals
!
      integer   i
!
      real(kind=kind_phys) :: qs1,  rat, thv1, restar,
     &                      czilc, tem1, tem2

      real(kind=kind_phys) :: tvs_ocn,  tvs_lnd,  tvs_ice,              &
     &                         z0_ocn,   z0_lnd,   z0_ice,              &
     &                      z0max_ocn,z0max_lnd,z0max_ice,              &
     &                      ztmax_ocn,ztmax_lnd,ztmax_ice
!
      real(kind=kind_phys), parameter ::
     &              charnock=.014, z0s_max=.317e-2                      &! a limiting value at high winds over sea
     &,             vis=1.4e-5, rnu=1.51e-5, visi=1.0/vis               &
     &,             log01=log(0.01), log05=log(0.05), log07=log(0.07)

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

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals, wind is wind speed,
!  surface roughness length is converted to m from cm
!
      do i=1,im

        ztmax_ocn = 0.; ztmax_lnd = 0.; ztmax_ice = 0.
        if(flag_iter(i)) then
          wind(i) = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))
     &                + max(0.0, min(ddvel(i), 30.0)), 1.0)
          tem1    = 1.0 + rvrdm1 * max(q1(i),1.e-8)
          thv1    = t1(i) * prslki(i) * tem1
          tvs_lnd = 0.5 * (tsurf_lnd(i)+tskin_lnd(i)) * tem1
          tvs_ice = 0.5 * (tsurf_ice(i)+tskin_ice(i)) * tem1
          tvs_ocn = 0.5 * (tsurf_ocn(i)+tskin_ocn(i)) * tem1
          qs1     = fpvs(t1(i))
          qs1     = max(1.0e-8, eps * qs1 / (prsl1(i) + epsm1 * qs1))

          z0_lnd    = 0.01 * z0rl_lnd(i)
          z0max_lnd = max(1.0e-6, min(z0_lnd,z1(i)))
          z0_ice    = 0.01 * z0rl_ice(i)
          z0max_ice = max(1.0e-6, min(z0_ice,z1(i)))
          z0_ocn    = 0.01 * z0rl_ocn(i)
          z0max_ocn = max(1.0e-6, min(z0_ocn,z1(i)))

!  compute stability dependent exchange coefficients
!  this portion of the code is presently suppressed
!

          if (wet(i)) then ! some open ocean
            ustar_ocn(i) = sqrt(grav * z0_ocn / charnock)

!**  test xubin's new z0

!           ztmax  = z0max

            restar = max(ustar_ocn(i)*z0max_ocn*visi, 0.000001)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat    = min(7.0, 2.67 * sqrt(sqrt(restar)) - 2.57)
            ztmax_ocn  = z0max_ocn * exp(-rat)
          endif ! Open ocean
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
     &                     * czilc*ca*sqrt(ustar_lnd(i)*(0.01/1.5e-05)))
            ztmax_ice = z0max_ice*exp( - tem1*tem1
     &                     * czilc*ca*sqrt(ustar_ice(i)*(0.01/1.5e-05)))


! mg, sfc-perts: add surface perturbations to ztmax/z0max ratio over land
            if (dry(i) .and. ztpert(i) /= 0.0) then
              ztmax_lnd = ztmax_lnd * (10.**ztpert(i))
            endif


          endif       ! end of if(sfctype flags) then

          ztmax_ocn  = max(ztmax_ocn,1.0e-6)
          ztmax_lnd  = max(ztmax_lnd,1.0e-6)
          ztmax_ice  = max(ztmax_ice,1.0e-6)

! BWG begin "stability" block, 2019-03-23
      if (wet(i)) then ! Some open ocean
          call stability
!  ---  inputs:
     &     (z1(i),snwdph_ocn(i),thv1,wind(i),
     &      z0max_ocn,ztmax_ocn,tvs_ocn,grav,
!  ---  outputs:
     &      rb_ocn(i),fm_ocn(i),fh_ocn(i),fm10_ocn(i),fh2_ocn(i),
     &      cm_ocn(i),ch_ocn(i),stress_ocn(i),ustar_ocn(i))
      endif ! Open ocean points

      if (dry(i)) then ! Some land
          call stability
!  ---  inputs:
     &     (z1(i),snwdph_lnd(i),thv1,wind(i),
     &      z0max_lnd,ztmax_lnd,tvs_lnd,grav,
!  ---  outputs:
     &      rb_lnd(i),fm_lnd(i),fh_lnd(i),fm10_lnd(i),fh2_lnd(i),
     &      cm_lnd(i),ch_lnd(i),stress_lnd(i),ustar_lnd(i))
      endif ! Dry points

      if (icy(i)) then ! Some ice
          call stability
!  ---  inputs:
     &     (z1(i),snwdph_ice(i),thv1,wind(i),
     &      z0max_ice,ztmax_ice,tvs_ice,grav,
!  ---  outputs:
     &      rb_ice(i),fm_ice(i),fh_ice(i),fm10_ice(i),fh2_ice(i),
     &      cm_ice(i),ch_ice(i),stress_ice(i),ustar_ice(i))
      endif ! Icy points

! BWG: Everything from here to end of subroutine was after
!      the stuff now put into "stability"

!
!  update z0 over ocean
!
          if (wet(i)) then
            z0_ocn = (charnock / grav) * ustar_ocn(i) * ustar_ocn(i)

! mbek -- toga-coare flux algorithm
!           z0 = (charnock / grav) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!           cc = ustar(i) * z0 / rnu
!           pp = cc / (1. + cc)
!           ff = grav * arnu / (charnock * ustar(i) ** 3)
!           z0 = arnu / (ustar(i) * ff ** pp)

            if (redrag) then
              z0rl_ocn(i) = 100.0 * max(min(z0_ocn, z0s_max), 1.e-7)
            else
              z0rl_ocn(i) = 100.0 * max(min(z0_ocn,.1), 1.e-7)
            endif
          endif              ! end of if(open ocean)
        endif                ! end of if(flagiter) loop
      enddo

      return
      end subroutine sfc_diff_run
!> @}

!----------------------------------------
!>\ingroup GFS_diff_main
      subroutine stability                                              &
     &     ( z1, snwdph, thv1, wind, z0max, ztmax, tvs, grav,           & !  ---  inputs:
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar)                !  ---  outputs:

!  ---  inputs:
      real(kind=kind_phys), intent(in) ::                               &
     &       z1, snwdph, thv1, wind, z0max, ztmax, tvs, grav

!  ---  outputs:
      real(kind=kind_phys), intent(out) ::                              &
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar

!  ---  locals:
      real(kind=kind_phys), parameter :: alpha=5., a0=-3.975            &
     &,             a1=12.32, alpha4=4.0*alpha
     &,             b1=-7.755,  b2=6.041,  alpha2=alpha+alpha, beta=1.0
     &,             a0p=-7.941, a1p=24.75, b1p=-8.705, b2p=7.899
     &,             ztmin1=-999.0

      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2,
     &                     z1i,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     hl110,  hlt,    hltinf, olinf,
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
          fm10    = log((z0max+10.)   * tem1)
          fh2     = log((ztmax+2.)    * tem2)
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

!---------------------------------
      end module sfc_diff

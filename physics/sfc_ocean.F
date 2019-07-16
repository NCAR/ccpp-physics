      module sfc_ocean
      implicit none
      private
      public :: sfc_ocean_init, sfc_ocean_run, sfc_ocean_finalize

      contains

!! \section arg_table_sfc_ocean_init  Argument Table
!!
      subroutine sfc_ocean_init()
      end subroutine sfc_ocean_init

!! \section arg_table_sfc_ocean_finalize  Argument Table
!!
      subroutine sfc_ocean_finalize()
      end subroutine sfc_ocean_finalize

#if 0
!! \section arg_table_sfc_ocean_run Argument Table
!! | local_name     | standard_name                                                                | long_name                                                            | units         | rank | type      |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|----------------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                               | count         |    0 | integer   |           | in     | F        |
!! | cp             | specific_heat_of_dry_air_at_constant_pressure                                | specific heat of dry air at constant pressure                        | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | rd             | gas_constant_dry_air                                                         | ideal gas constant for dry air                                       | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | eps            | ratio_of_dry_air_to_water_vapor_gas_constants                                | rd/rv                                                                | none          |    0 | real      | kind_phys | in     | F        |
!! | epsm1          | ratio_of_dry_air_to_water_vapor_gas_constants_minus_one                      | (rd/rv) - 1                                                          | none          |    0 | real      | kind_phys | in     | F        |
!! | hvap           | latent_heat_of_vaporization_of_water_at_0C                                   | latent heat of evaporation/sublimation                               | J kg-1        |    0 | real      | kind_phys | in     | F        |
!! | rvrdm1         | ratio_of_vapor_to_dry_air_gas_constants_minus_one                            | (rv/rd) - 1 (rv = ideal gas constant for water vapor)                | none          |    0 | real      | kind_phys | in     | F        |
!! | ps             | surface_air_pressure                                                         | surface pressure                                                     | Pa            |    1 | real      | kind_phys | in     | F        |
!! | u1             | x_wind_at_lowest_model_layer                                                 | x component of surface layer wind                                    | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1             | y_wind_at_lowest_model_layer                                                 | y component of surface layer wind                                    | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature_at_lowest_model_layer                                        | surface layer mean temperature                                       | K             |    1 | real      | kind_phys | in     | F        |
!! | q1             | water_vapor_specific_humidity_at_lowest_model_layer                          | surface layer mean specific humidity                                 | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | tskin          | surface_skin_temperature_over_ocean_interstitial                             | surface skin temperature over ocean (temporary use as interstitial)  | K             |    1 | real      | kind_phys | in     | F        |
!! | cm             | surface_drag_coefficient_for_momentum_in_air_over_ocean                      | surface exchange coeff for momentum over ocean                       | none          |    1 | real      | kind_phys | in     | F        |
!! | ch             | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ocean             | surface exchange coeff heat & moisture over ocean                    | none          |    1 | real      | kind_phys | in     | F        |
!! | prsl1          | air_pressure_at_lowest_model_layer                                           | surface layer mean pressure                                          | Pa            |    1 | real      | kind_phys | in     | F        |
!! | prslki         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer          | ratio         |    1 | real      | kind_phys | in     | F        |
!! | wet            | flag_nonzero_wet_surface_fraction                                            | flag indicating presence of some ocean or lake surface area fraction | flag          |    1 | logical   |           | in     | F        |
!! | ddvel          | surface_wind_enhancement_due_to_convection                                   | wind enhancement due to convection                                   | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | flag_iter      | flag_for_iteration                                                           | flag for iteration                                                   | flag          |    1 | logical   |           | in     | F        |
!! | qsurf          | surface_specific_humidity_over_ocean                                         | surface air saturation specific humidity over ocean                  | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | cmm            | surface_drag_wind_speed_for_momentum_in_air_over_ocean                       | momentum exchange coefficient over ocean                             | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | chh            | surface_drag_mass_flux_for_heat_and_moisture_in_air_over_ocean               | thermal exchange coefficient over ocean                              | kg m-2 s-1    |    1 | real      | kind_phys | inout  | F        |
!! | gflux          | upward_heat_flux_in_soil_over_ocean                                          | soil heat flux over ocean                                            | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux_over_ocean                         | kinematic surface upward latent heat flux over ocean                 | kg kg-1 m s-1 |    1 | real      | kind_phys | inout  | F        |
!! | hflx           | kinematic_surface_upward_sensible_heat_flux_over_ocean                       | kinematic surface upward sensible heat flux over ocean               | K m s-1       |    1 | real      | kind_phys | inout  | F        |
!! | ep             | surface_upward_potential_latent_heat_flux_over_ocean                         | surface upward potential latent heat flux over ocean                 | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                           | error message for error handling in CCPP                             | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                              | error flag for error handling in CCPP                                | flag          |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine sfc_ocean_run                                          &
!...................................
!  ---  inputs:
     &     ( im, cp, rd, eps, epsm1, hvap, rvrdm1, ps, u1, v1, t1, q1,  &
     &       tskin, cm, ch, prsl1, prslki, wet, ddvel,                  &
     &       flag_iter,                                                 &
!  ---  outputs:
     &       qsurf, cmm, chh, gflux, evap, hflx, ep,                    &
     &       errmsg, errflg                                             &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_ocean                                                     !
!       inputs:                                                         !
!          ( im, ps, u1, v1, t1, q1, tskin, cm, ch,                     !
!            prsl1, prslki, wet, ddvel, flag_iter,                      !
!       outputs:                                                        !
!            qsurf, cmm, chh, gflux, evap, hflx, ep )                   !
!                                                                       !
!                                                                       !
!  subprograms/functions called: fpvs                                   !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!         2005  -- created from the original progtm to account for      !
!                  ocean only                                           !
!    oct  2006  -- h. wei      added cmm and chh to the output          !
!    apr  2009  -- y.-t. hou   modified to match the modified gbphys.f  !
!                  reformatted the code and added program documentation !
!    sep  2009  -- s. moorthi removed rcl and made pa as pressure unit  !
!                  and furthur reformatted the code                     !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horizontal dimension                     1    !
!     ps       - real, surface pressure                            im   !
!     u1, v1   - real, u/v component of surface layer wind         im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     wet      - logical, =T if any ocean/lak, =F otherwise        im   !
!     ddvel    - real, wind enhancement due to convection (m/s)    im   !
!     flag_iter- logical,                                          im   !
!                                                                       !
!  outputs:                                                             !
!     qsurf    - real, specific humidity at sfc                    im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     gflux    - real, ground heat flux (zero for ocean)           im   !
!     evap     - real, evaporation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!     ep       - real, potential evaporation                       im   !
!                                                                       !
! ===================================================================== !
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: im
      real (kind=kind_phys), intent(in) :: cp, rd, eps, epsm1, hvap,    &
     &      rvrdm1

      real (kind=kind_phys), dimension(im), intent(in) :: ps, u1, v1,   &
     &      t1, q1, tskin, cm, ch, prsl1, prslki, ddvel

      logical, dimension(im), intent(in) :: flag_iter, wet

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: qsurf,     &
     &       cmm, chh, gflux, evap, hflx, ep

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:

      real (kind=kind_phys) :: q0, qss, rch, rho, wind, tem, cpinv,     &
     &      hvapi, elocp

      integer :: i

      logical :: flag(im)
!
!===> ...  begin here
      cpinv  = 1.0/cp
      hvapi  = 1.0/hvap
      elocp  = hvap/cp
!
!  -- ...  initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!  --- ...  flag for open water
      do i = 1, im
        flag(i) = (wet(i) .and. flag_iter(i))

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specified
!           ps is in pascals, wind is wind speed,
!           rho is density, qss is sat. hum. at surface

        if ( flag(i) ) then

          wind     = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))                &
     &                 + max( 0.0, min( ddvel(i), 30.0 ) ), 1.0)

          q0       = max( q1(i), 1.0e-8 )
          rho      = prsl1(i) / (rd*t1(i)*(1.0 + rvrdm1*q0))

          qss      = fpvs( tskin(i) )
          qss      = eps*qss / (ps(i) + epsm1*qss)

          evap(i)  = 0.0
          hflx(i)  = 0.0
          ep(i)    = 0.0
          gflux(i) = 0.0

!  --- ...    rcp  = rho cp ch v

          rch      = rho * cp * ch(i) * wind
          cmm(i)   = cm(i) * wind
          chh(i)   = rho * ch(i) * wind

!  --- ...  sensible and latent heat flux over open water

          hflx(i)  = rch * (tskin(i) - t1(i) * prslki(i))

          evap(i)  = elocp*rch * (qss - q0)
          qsurf(i) = qss

          tem      = 1.0 / rho
          hflx(i)  = hflx(i) * tem * cpinv
          evap(i)  = evap(i) * tem * hvapi
        endif
      enddo
!
      return
!...................................
      end subroutine sfc_ocean_run
!-----------------------------------
      end module sfc_ocean

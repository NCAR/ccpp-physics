!>  \file sfc_cice.f
!!  This file contains the sfc_sice for coupling to CICE

!> This module contains the CCPP-compliant GFS sea ice post
!! interstitial codes, which returns updated ice thickness and 
!! concentration to global arrays where there is no ice, and 
!! set temperature to surface skin temperature.

!> This module contains the CCPP-compliant GFS sea ice scheme.
      module sfc_cice

      contains

      subroutine sfc_cice_init
      end subroutine sfc_cice_init
!
      subroutine sfc_cice_finalize
      end subroutine sfc_cice_finalize


!> \defgroup sfc_sice for coupling to CICE
!! @{
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication
!!
!> \brief Brief description of the subroutine
!!
!! \section arg_table_cice_run Arguments
!! | local_name     | standard_name                                                                | long_name                                                       | units         | rank | type      |    kind   | intent | optional | !!
!! |----------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                          | count         |    0 | integer   |           | in     | F        |
!! | cplflx         | flag_for_flux_coupling                                                       | flag controlling cplflx collection (default off)                | flag          |    0 | logical   |           | in     | F
!! | cplchm         | flag_for_chemistry_coupling                                                    | flag controlling cplchm collection (default off)                | flag          |    0 | logical   |           | in     | F        |
!! | hvap           | latent_heat_of_vaporization_of_water_at_0C                                   | latent heat of evaporation/sublimation                          | J kg-1        |    0 | real      | kind_phys | in     | F        | 
!! | cp             | specific_heat_of_dry_air_at_constant_pressure                                | specific heat of dry air at constant pressure                   | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | rvrdm1         | ratio_of_vapor_to_dry_air_gas_constants_minus_one                            | (rv/rd) - 1 (rv = ideal gas constant for water vapor)           | none          |    0 | real      | kind_phys | in     | F        |
!! | rd             | gas_constant_dry_air                                                         | ideal gas constant for dry air                                  | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | u1             | x_wind_at_lowest_model_layer                                                 | u component of surface layer wind                               | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1             | y_wind_at_lowest_model_layer                                                 | v component of surface layer wind                               | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature_at_lowest_model_layer                                        | surface layer mean temperature                                  | K             |    1 | real      | kind_phys | in     | F        |
!! | q1             | water_vapor_specific_humidity_at_lowest_model_layer                          | surface layer mean specific humidity                            | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | cm             | surface_drag_coefficient_for_momentum_in_air_over_ice                        | surface exchange coeff for momentum over ice                    | none          |    1 | real      | kind_phys | in     | F        |
!! | ch             | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ice               | surface exchange coeff heat & moisture over ice                 | none          |    1 | real      | kind_phys | in     | F        
!! | prsl1          | air_pressure_at_lowest_model_layer                                           | surface layer mean pressure                                     | Pa            |    1 | real      | kind_phys | in     | F        |
!! | prslki         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer     | ratio         |    1 | real      | kind_phys | in     | F        |
!! | islimsk        | sea_land_ice_mask                                                            | sea/land/ice mask (=0/1/2)                                      | flag          |    1 | integer   |           | in     | F        |
!! | ddvel          | surface_wind_enhancement_due_to_convection                                   | wind enhancement due to convection                              | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | flag_iter      | flag_for_iteration                                                           | flag for iteration                                              | flag          |    1 | logical   |           | in     | F        |
!! | dqsfc          | dqsfcin                                                                      | aoi_fld%dqsfcin(item,lan)                                       |               |    1 | real      | kind_phys | none   | F        |
!! | dtsfc          | dtsfcin                                                                      | aoi_fld%dtsfcin(item,lan)                                       |               |    1 | real      | kind_phys | none   | F        |
!! | qsurf          | surface_specific_humidity_over_ice                                           | surface air saturation specific humidity over ice               | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | cmm            | surface_drag_wind_speed_for_momentum_in_air_over_ice                         | momentum exchange coefficient over ice                          | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | chh            | surface_drag_mass_flux_for_heat_and_moisture_in_air_over_ice                 | thermal exchange coefficient over ice                           | kg m-2 s-1    |    1 | real      | kind_phys | inout  | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux_over_ice                           | kinematic surface upward latent heat flux over ice              | kg kg-1 m s-1 |    1 | real      | kind_phys | inout  | F        |
!! | hflx           | kinematic_surface_upward_sensible_heat_flux_over_ice                         | kinematic surface upward sensible heat flux over ice            | K m s-1       |    1 | real      | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                           | error message for error handling in CCPP                        | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                              | error flag for error handling in CCPP                           | flag          |    0 | integer   |           | out    | F        |
!!

!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{


!!      use physcons, only : hvap => con_hvap,  cp => con_cp,             &
!!    &                     rvrdm1 => con_fvirt, rd => con_rd
!
!-----------------------------------
      subroutine sfc_cice_run                                           &
     &     ( im, cplflx, cplchm, hvap, cp, rvrdm1, rd,                  & ! ---  inputs:
     &       u1, v1, t1, q1, cm, ch, prsl1, prslki,                     &
     &       islimsk, ddvel, flag_iter, dqsfc, dtsfc,                   &
     &       qsurf, cmm, chh, evap, hflx,                               & ! ---  outputs:
     &       errmsg, errflg
     &     )

! ===================================================================== !
!  description:                                                         !
!  Sep 2015  --  Xingren Wu created from sfc_sice for coupling to CICE  !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_cice                                                      !
!       inputs:                                                         !
!          ( im, u1, v1, t1, q1, cm, ch, prsl1, prslki,                 !
!            islimsk, ddvel, flag_iter, dqsfc, dtsfc,                   !
!       outputs:                                                        !
!            qsurf, cmm, chh, evap, hflx)                               !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:
!     im, - integer, horiz dimension
!     u1, v1   - real, u/v component of surface layer wind
!     t1       - real, surface layer mean temperature ( k )
!     q1       - real, surface layer mean specific humidity
!     cm       - real, surface exchange coeff for momentum (m/s)
!     ch       - real, surface exchange coeff heat & moisture(m/s)
!     prsl1    - real, surface layer mean pressure
!     prslki   - real, ?
!     islimsk  - integer, sea/land/ice mask
!     ddvel    - real, ?
!     flag_iter- logical
!     dqsfc    - real, latent heat flux
!     dtsfc    - real, sensible heat flux
!  outputs:
!     qsurf    - real, specific humidity at sfc
!     cmm      - real, ?
!     chh      - real, ?
!     evap     - real, evaperation from latent heat
!     hflx     - real, sensible heat
!  ====================    end of description    =====================  !
!
!
      use machine , only : kind_phys
      implicit none


      real (kind=kind_phys), intent(in) :: hvap, cp, rvrdm1, rd

!  ---  inputs:
      integer, intent(in) :: im
      logical, intent(in) :: cplflx
      logical, intent(in) :: cplchm

      real (kind=kind_phys), dimension(im), intent(in) :: u1, v1,       &
     &       t1, q1, cm, ch, prsl1, prslki, ddvel, dqsfc, dtsfc

      integer, dimension(im), intent(in) :: islimsk

      logical, intent(in) :: flag_iter(im)

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(out) :: qsurf,       &
     &       cmm, chh, evap, hflx
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: q0, rch, rho, tv1, wind

      real (kind=kind_phys) :: tem

      real(kind=kind_phys) :: cpinv, hvapi, elocp

      integer :: i
 
      logical :: flag(im)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      cpinv = 1.0/cp
      hvapi = 1.0/hvap
      elocp = hvap/cp
!
      if((.not. cplflx).and.(.not.cplchm))then
        return
      endif

      do i = 1, im
         flag(i) = (islimsk(i) == 4) .and. flag_iter(i)
      enddo
!
      do i = 1, im
        if (flag(i)) then

          wind(i)   = sqrt(u1(i)*u1(i) + v1(i)*v1(i))                   &
     &              + max(0.0, min(ddvel(i), 30.0))
          wind(i)   = max(wind(i), 1.0)

          q0(i)     = max(q1(i), 1.0e-8)
          tv1(i)    = t1(i) * (1.0 + rvrdm1*q0(i))
          rho(i)    = prsl1(i) / (rd*tv1(i))

          cmm(i) = cm(i)  * wind(i)
          chh(i) = rho(i) * ch(i) * wind(i)
          rch(i) = chh(i) * cp

          qsurf(i) = q1(i) + dqsfc(i) / (elocp*rch(i))
          tem     = 1.0 / rho(i)
          hflx(i) = dtsfc(i) * tem * cpinv
          evap(i) = dqsfc(i) * tem * hvapi
        endif
      enddo
 
      return
!-----------------------------------
      end subroutine sfc_cice_run
!-----------------------------------

!> @}
      end module sfc_cice

!> \file gmtb_scm_sfc_flux_spec.F90
!!  Contains code to calculate parameters needed by the rest of the GFS physics suite given specified surface fluxes.

module gmtb_scm_sfc_flux_spec

  implicit none

  private

!----------------
! Public entities
!----------------
  public  gmtb_scm_sfc_flux_spec_init, gmtb_scm_sfc_flux_spec_run, gmtb_scm_sfc_flux_spec_finalize

  CONTAINS
!*******************************************************************************************

  subroutine gmtb_scm_sfc_flux_spec_init()
  end subroutine gmtb_scm_sfc_flux_spec_init

  subroutine gmtb_scm_sfc_flux_spec_finalize()
  end subroutine gmtb_scm_sfc_flux_spec_finalize

!> \brief This routine calculates surface-related parameters given specified sensible and latent heat fluxes and a roughness length. Most of the calculation
!! is "backing out" parameters that are calculated in sfc_dff.f from the known surface heat fluxes and roughness length.
!!
!! \section arg_table_gmtb_scm_sfc_flux_spec_run Argument Table
!! | local_name       | standard_name                                                                | long_name                                                       | units         | rank | type      |    kind   | intent | optional |
!! |------------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | u1               | x_wind_at_lowest_model_layer                                                 | x component of 1st model layer wind                             | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1               | y_wind_at_lowest_model_layer                                                 | y component of 1st model layer wind                             | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | z1               | height_above_ground_at_lowest_model_layer                                    | height above ground at 1st model layer                          | m             |    1 | real      | kind_phys | in     | F        |
!! | t1               | air_temperature_at_lowest_model_layer                                        | 1st model layer air temperature                                 | K             |    1 | real      | kind_phys | in     | F        |
!! | q1               | water_vapor_specific_humidity_at_lowest_model_layer                          | 1st model layer specific humidity                               | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | p1               | air_pressure_at_lowest_model_layer                                           | Model layer 1 mean pressure                                     | Pa            |    1 | real      | kind_phys | in     | F        |
!! | roughness_length | surface_roughness_length                                                     | surface roughness length                                        | cm            |    1 | real      | kind_phys | in     | F        |
!! | spec_sh_flux     | specified_kinematic_surface_upward_sensible_heat_flux                        | specified kinematic surface upward sensible heat flux           | K m s-1       |    1 | real      | kind_phys | in     | F        |
!! | spec_lh_flux     | specified_kinematic_surface_upward_latent_heat_flux                          | specified kinematic surface upward latent heat flux             | kg kg-1 m s-1 |    1 | real      | kind_phys | in     | F        |
!! | exner_inverse    | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer     | ratio         |    1 | real      | kind_phys | in     | F        |
!! | T_surf           | surface_skin_temperature                                                     | surface skin temperature                                        | K             |    1 | real      | kind_phys | in     | F        |
!! | cp               | specific_heat_of_dry_air_at_constant_pressure                                | specific heat of dry air at constant pressure                   | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | grav             | gravitational_acceleration                                                   | gravitational acceleration                                      | m s-2         |    0 | real      | kind_phys | in     | F        |
!! | hvap             | latent_heat_of_vaporization_of_water_at_0C                                   | latent heat of vaporization of water at 0C                      | J kg-1        |    0 | real      | kind_phys | in     | F        |
!! | rd               | gas_constant_dry_air                                                         | ideal gas constant for dry air                                  | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | fvirt            | ratio_of_vapor_to_dry_air_gas_constants_minus_one                            | rv/rd - 1 (rv = ideal gas constant for water vapor)             | none          |    0 | real      | kind_phys | in     | F        |
!! | vonKarman        | vonKarman_constant                                                           | vonKarman constant                                              | none          |    0 | real      | kind_phys | in     | F        |
!! | sh_flux          | kinematic_surface_upward_sensible_heat_flux                                  | surface upward sensible heat flux                               | K m s-1       |    1 | real      | kind_phys | out    | F        |
!! | lh_flux          | kinematic_surface_upward_latent_heat_flux                                    | surface upward evaporation flux                                 | kg kg-1 m s-1 |    1 | real      | kind_phys | out    | F        |
!! | u_star           | surface_friction_velocity                                                    | boundary layer parameter                                        | m s-1         |    1 | real      | kind_phys | out    | F        |
!! | sfc_stress       | surface_wind_stress                                                          | surface wind stress                                             | m2 s-2        |    1 | real      | kind_phys | out    | F        |
!! | cm               | surface_drag_coefficient_for_momentum_in_air                                 | surface exchange coeff for momentum                             | none          |    1 | real      | kind_phys | out    | F        |
!! | ch               | surface_drag_coefficient_for_heat_and_moisture_in_air                        | surface exchange coeff heat & moisture                          | none          |    1 | real      | kind_phys | out    | F        |
!! | fm               | Monin-Obukhov_similarity_function_for_momentum                               | Monin-Obukhov similarity function for momentum                  | none          |    1 | real      | kind_phys | out    | F        |
!! | fh               | Monin-Obukhov_similarity_function_for_heat                                   | Monin-Obukhov similarity function for heat                      | none          |    1 | real      | kind_phys | out    | F        |
!! | rb               | bulk_richardson_number_at_lowest_model_level                                 | bulk Richardson number at the surface                           | none          |    1 | real      | kind_phys | out    | F        |
!! | u10m             | x_wind_at_10m                                                                | 10 meter u wind speed                                           | m s-1         |    1 | real      | kind_phys | out    | F        |
!! | v10m             | y_wind_at_10m                                                                | 10 meter v wind speed                                           | m s-1         |    1 | real      | kind_phys | out    | F        |
!! | wind1            | wind_speed_at_lowest_model_layer                                             | wind speed at lowest model level                                | m s-1         |    1 | real      | kind_phys | out    | F        |
!! | qss              | surface_specific_humidity                                                    | surface air saturation specific humidity                        | kg kg-1       |    1 | real      | kind_phys | out    | F        |
!! | t2m              | temperature_at_2m                                                            | 2 meter temperature                                             | K             |    1 | real      | kind_phys | out    | F        |
!! | q2m              | specific_humidity_at_2m                                                      | 2 meter specific humidity                                       | kg kg-1       |    1 | real      | kind_phys | out    | F        |
!! | errmsg           | ccpp_error_message                                                           | error message for error handling in CCPP                        | none          |    0 | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                              | error flag for error handling in CCPP                           | flag          |    0 | integer   |           | out    | F        |
!!
!! \section general_sfc_flux_spec General Algorithm
!!  -# Compute friction velocity from the wind speed at the lowest model layer, the height about the ground, and the roughness length.
!!  -# Compute the surface stress from the friction velocity.
!!  -# Calculate the surface drag coefficient for momentum given the surface stress and wind on the lowest model layer.
!!  -# Calculate the Monin-Obukhov similarity funciton for momentum from the surface drag coefficient.
!!  -# Calculate the Obukhov length from the friction velocity, surface virtual potential temperature, and surface vertical virtual potential temperature flux.
!!  -# Calculate the bulk Richardson number at the lowest model layer.
!!  -# Calculate the Monin-Obukhov similarity function for heat and moisture from the bulk Richardson number and diagnosed similarity function for momentum.
!!  -# Calculate the surface drag coefficient for heat and moisture.
!!  -# Calculate the u and v wind at 10m.
  subroutine gmtb_scm_sfc_flux_spec_run (u1, v1, z1, t1, q1, p1, roughness_length, spec_sh_flux, spec_lh_flux, &
    exner_inverse, T_surf, cp, grav, hvap, rd, fvirt, vonKarman, sh_flux, lh_flux, u_star, sfc_stress, cm, ch, &
    fm, fh, rb, u10m, v10m, wind1, qss, t2m, q2m, errmsg, errflg)

    use machine,             only: kind_phys

    real(kind=kind_phys), intent(in) :: u1(:), v1(:), z1(:), t1(:), q1(:), p1(:), roughness_length(:), &
      spec_sh_flux(:), spec_lh_flux(:), exner_inverse(:), T_surf(:)
    real(kind=kind_phys), intent(in) :: cp, grav, hvap, rd, fvirt, vonKarman
    real(kind=kind_phys), intent(out) :: sh_flux(:), lh_flux(:), u_star(:), sfc_stress(:), &
      cm(:), ch(:), fm(:), fh(:), rb(:), u10m(:), v10m(:), wind1(:), qss(:), t2m(:), q2m(:)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i

    real(kind=kind_phys) :: rho, q1_non_neg, w_thv1, rho_cp_inverse, rho_hvap_inverse, Obukhov_length, thv1, tvs, &
      dtv, adtv, wind10m, u_fraction, roughness_length_m

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

!     !--- set control properties (including namelist read)
  !calculate u_star from wind profiles (need roughness length, and wind and height at lowest model level)
  do i=1, size(z1)
    sh_flux(i) = spec_sh_flux(i)
    lh_flux(i) = spec_lh_flux(i)

    roughness_length_m = 0.01*roughness_length(i)

    wind1(i) = sqrt(u1(i)*u1(i) + v1(i)*v1(i))
    u_star(i) = vonKarman*wind1(i)/(log(z1(i)/roughness_length_m))

    !calculate variables related to u_star
    sfc_stress(i) = u_star(i)*u_star(i)

    if(wind1(i) > 0.0) then
      cm(i) = sfc_stress(i)/(wind1(i)*wind1(i))
    else
      cm(i) = 0.0
    end if
    fm(i) = sqrt((vonKarman*vonKarman)/cm(i))

    !calculate the Obukhov length from the specified surface fluxes
    q1_non_neg       = max( q1(i), 1.0e-8 )
    rho      = p1(i) / (rd*t1(i)*(1.0 + fvirt*q1_non_neg))
    rho_cp_inverse = 1.0/(rho*cp)
    rho_hvap_inverse = 1.0/(rho*hvap)

    thv1    = t1(i) * exner_inverse(i) * (1.0 + fvirt*q1_non_neg)
    ! sh_flux = rho_cp_inverse*sh_flux_wm2(i)
    ! lh_flux = rho_hvap_inverse*lh_flux_wm2(i)
    w_thv1 = sh_flux(i)*exner_inverse(i) + fvirt*t1(i)*lh_flux(i)

    Obukhov_length = -u_star(i)*u_star(i)*u_star(i)*thv1/(vonKarman*grav*w_thv1)

    !calculate the bulk Richardson number and the M-O function for scalars
    tvs     = T_surf(i)*(1.0 + fvirt*q1_non_neg)

    dtv     = thv1 - tvs
    adtv    = max(abs(dtv),0.001)
    dtv     = sign(1._kind_phys,dtv) * adtv
    rb(i)   = max(-5000.0, (grav+grav) * dtv * z1(i) / ((thv1 + tvs) * wind1(i) * wind1(i)))

    fh(i) = rb(i)*fm(i)*fm(i)*Obukhov_length/z1(i)
    ch(i) = vonKarman*vonKarman/(fm(i)*fh(i))

    !calculate sfc_diag variables (bypassing t2m and q2m for now)
    !should calculate qss, but it is not needed in moninedmf (q2m depends on qss - will implement later)

    wind10m = u_star(i)/vonKarman*log(10.0/roughness_length_m)
    u_fraction = sqrt(1.0 - v1(i)**2/wind1(i)**2)
    u10m(i) = u_fraction*wind10m
    v10m(i) = sqrt(wind10m**2 - u10m(i)**2)

    qss(i) = 0.0
    t2m(i) = 0.0
    q2m(i) = 0.0
  end do

  end subroutine gmtb_scm_sfc_flux_spec_run

end module gmtb_scm_sfc_flux_spec

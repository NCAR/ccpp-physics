module gmtb_scm_sfc_flux_spec

  implicit none

  private

!----------------
! Public entities
!----------------
  public  gmtb_scm_sfc_flux_spec_init, gmtb_scm_sfc_flux_spec_run, gmtb_scm_sfc_flux_spec_finalize

  CONTAINS
!*******************************************************************************************

!--------------
! GFS initialze
!--------------

  subroutine gmtb_scm_sfc_flux_spec_init()
  end subroutine gmtb_scm_sfc_flux_spec_init

  subroutine gmtb_scm_sfc_flux_spec_finalize()
  end subroutine gmtb_scm_sfc_flux_spec_finalize

!> \section arg_table_gmtb_scm_sfc_flux_spec_run Argument Table
!! | local_name           | standard_name                                               | long_name                                                               | units         | rank | type                          |    kind   | intent | optional |
!! |----------------------|-------------------------------------------------------------|-------------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | u1             | x_wind_at_lowest_model_layer                                                 | x component of 1st model layer wind                             | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1             | y_wind_at_lowest_model_layer                                                 | y component of 1st model layer wind                             | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | z1             | height_above_mean_sea_level_at_lowest_model_layer                            | height above mean sea level at 1st model layer              | m          |    1 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature_at_lowest_model_layer                                        | 1st model layer air temperature                                 | K             |    1 | real      | kind_phys | in     | F        |
!! | q1             | specific_humidity_at_lowest_model_layer                                      | 1st model layer specific humidity                               | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | p1          | air_pressure_at_lowest_model_layer                                           | Model layer 1 mean pressure                                     | Pa            |    1 | real      | kind_phys | in     | F        |
!! | roughness_length           | surface_roughness_length                                                     | surface roughness length                                        | cm            |    1 | real      | kind_phys | in  | F        |
!! | sh_flux           | kinematic_surface_upward_sensible_heat_flux                                  | surface upward sensible heat flux                               | K m s-1       |    1 | real      | kind_phys | in  | F        |
!! | lh_flux           | kinematic_surface_upward_latent_heat_flux                                    | surface upward evaporation flux                                 | kg kg-1 m s-1 |    1 | real      | kind_phys | in  | F        |
!! | exner_inverse         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer     | ratio         |    1 | real      | kind_phys | in     | F        |
!! | T_surf      | surface_skin_temperature                                               | ocean surface skin temperature                       | K             |    1 | real    | kind_phys | in   | F        |
!! | u_star    | surface_friction_velocity                                              | boundary layer parameter                             | m s-1         |    1 | real    | kind_phys | out   | F        |
!! | errmsg               | error_message                                               | error message for error handling in CCPP                                | none          |    0 | character                     | len=*     | out    | F        |
!! | errflg               | error_flag                                                  | error flag for error handling in CCPP                                   | flag          |    0 | integer                       |           | out    | F        |
!!
  subroutine gmtb_scm_sfc_flux_spec_run (u1, v1, z1, t1, q1, p1, roughness_length, sh_flux, lh_flux, exner_inverse, &
    T_surf, u_star, sfc_stress, cm, ch, fm, fh, rb, u10m, v10m, wind1, qss, t2m, q2m, &
                             errmsg, errflg)

    use machine,             only: kind_phys
    use GFS_typedefs,        only: GFS_init_type,                          &
                                   GFS_statein_type,  GFS_stateout_type,   &
                                   GFS_sfcprop_type,  GFS_coupling_type,   &
                                   GFS_control_type,  GFS_grid_type,       &
                                   GFS_tbd_type,      GFS_cldprop_type,    &
                                   GFS_radtend_type,  GFS_diag_type,       &
                                   GFS_sfccycle_type, GFS_interstitial_type
    use funcphys,            only: gfuncphys
    use module_microphysics, only: gsmconst
    use cldwat2m_micro,      only: ini_micro
    use aer_cloud,           only: aer_cloud_init
    use module_ras,          only: ras_init
    use ozne_def,            only: latsozp, levozp, timeoz, oz_coeff, oz_lat, oz_pres, oz_time, ozplin

    !--- interface variables
    type(GFS_control_type),      intent(inout) :: Model
    type(GFS_statein_type),      intent(inout) :: Statein
    type(GFS_stateout_type),     intent(inout) :: Stateout
    type(GFS_sfcprop_type),      intent(inout) :: Sfcprop
    type(GFS_coupling_type),     intent(inout) :: Coupling
    type(GFS_grid_type),         intent(inout) :: Grid
    type(GFS_tbd_type),          intent(inout) :: Tbd
    type(GFS_cldprop_type),      intent(inout) :: Cldprop
    type(GFS_radtend_type),      intent(inout) :: Radtend
    type(GFS_diag_type),         intent(inout) :: Diag
    type(GFS_sfccycle_type),     intent(inout) :: Sfccycle
    type(GFS_interstitial_type), intent(inout) :: Interstitial
    type(GFS_init_type),         intent(in)    :: Init_parm

    integer, intent(in) :: n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times
    real(kind=kind_phys), intent(in) :: ozone_lat(:), ozone_pres(:), ozone_time(:), ozone_forcing_in(:,:,:,:)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    !
    !     !--- local variables
    real(kind=kind_phys), allocatable :: si(:)
    real(kind=kind_phys), parameter   :: p_ref = 101325.0d0

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    roughness_length_m(:) = 0.01*roughness_length(:)

!     !--- set control properties (including namelist read)
  !calculate u_star from wind profiles (need roughness length, and wind and height at lowest model level)
  do i=1, size(z1)


    wind1(i) = sqrt(u1(i)*u1(i) + v1(i)*v1(i))
    u_star(i) = vonKarman*wind1(i)/(log(z1(i)/roughness_length_m(i)))

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
    rho      = p1(i) / (rd*t1(i)*(1.0 + rvrdm1*q1_non_neg))
    rho_cp_inverse = 1.0/(rho*cp)
    rho_hvap_inverse = 1.0/(rho*hvap)

    thv1    = t1(i) * exner_inverse(i) * (1.0 + rvrdm1*q1_non_neg)
    ! sh_flux = rho_cp_inverse*sh_flux_wm2(i)
    ! lh_flux = rho_hvap_inverse*lh_flux_wm2(i)
    w_thv1 = sh_flux(i)*exner_inverse(i) + rvrdm1*t1(i)*lh_flux(i)

    Obukhov_length = -u_star(i)*u_star(i)*u_star(i)*thv1/(vonKarman*grav*w_thv1)

    !calculate the bulk Richardson number and the M-O function for scalars
    tvs     = T_surf(i)*(1.0 + rvrdm1*q1_non_neg)

    dtv     = thv1 - tvs
    adtv    = max(abs(dtv),0.001)
    dtv     = sign(1._dp,dtv) * adtv
    rb(i)   = max(-5000.0, (grav+grav) * dtv * z1(i) / ((thv1 + tvs) * wind1(i) * wind1(i)))

    fh(i) = rb(i)*fm(i)*fm(i)*Obukhov_length/z1(i)
    ch(i) = vonKarman*vonKarman/(fm(i)*fh(i))

    !calculate sfc_diag variables (bypassing t2m and q2m for now)
    !should calculate qss, but it is not needed in moninedmf (q2m depends on qss - will implement later)

    wind10m = u_star(i)/vonKarman*log(10.0/roughness_length_m(i))
    u_fraction = sqrt(1.0 - v1(i)**2/wind1(i)**2)
    u10m(i) = u_fraction*wind10m
    v10m(i) = sqrt(wind10m**2 - u10m(i)**2)

    qss(i) = 0.0
    t2m(i) = 0.0
    q2m(i) = 0.0
  end do

  end subroutine gmtb_scm_sfc_flux_spec_run

end module gmtb_scm_sfc_flux_spec

!> \file module_myjsfc_wrapper.F90
!!  Contains all of the code related to running the MYJ surface layer scheme 

      MODULE myjsfc_wrapper

      contains

      subroutine myjsfc_wrapper_init ()
      end subroutine myjsfc_wrapper_init

      subroutine myjsfc_wrapper_finalize ()
      end subroutine myjsfc_wrapper_finalize

!!
!> \brief This scheme (1) performs pre-myjsfc work, (20 runs the myj sfc layer scheme, and (3) performs post-myjsfc work
#if 0
!! \section arg_table_myjsfc_wrapper_run Argument Table
!! | local_name     | standard_name                                                               | long_name                                             | units         | rank | type        |    kind   | intent | optional |
!! |----------------|-----------------------------------------------------------------------------|-------------------------------------------------------|---------------|------|-------------|-----------|--------|----------|
!! | restart        | flag_for_restart                                                            | flag for restart (warmstart) or coldstart             | flag          |    0 | logical     |           | in     | F        |
!! | ix             | horizontal_dimension                                                        | horizontal dimension                                  | count         |    0 | integer     |           | in     | F        |
!! | im             | horizontal_loop_extent                                                      | horizontal loop extent                                | count         |    0 | integer     |           | in     | F        |
!! | levs           | vertical_dimension                                                          | vertical layer dimension                              | count         |    0 | integer     |           | in     | F        |
!! | kdt            | index_of_time_step                                                          | current time step index                               | index         |    0 | integer     |           | in     | F        |
!! | ntrac          | number_of_tracers                                                           | number of tracers                                     | count         |    0 | integer     |           | in     | F        |
!! | ntke           | index_for_turbulent_kinetic_energy                                          | tracer index for turbulent kinetic energy             | index         |    0 | integer     |           | in     | F        |
!! | ntcw           | index_for_liquid_cloud_condensate                                           | cloud condensate index in tracer array                | index         |    0 | integer     |           | in     | F        |
!! | ntiw           | index_for_ice_cloud_condensate                                              | tracer index for  ice water                           | index         |    0 | integer     |           | in     | F        |
!! | ntrw           | index_for_rain_water                                                        | tracer index for rain water                           | index         |    0 | integer     |           | in     | F        |
!! | ntsw           | index_for_snow_water                                                        | tracer index for snow water                           | index         |    0 | integer     |           | in     | F        |
!! | ntgl           | index_for_graupel                                                           | tracer index for graupel                              | index         |    0 | integer     |           | in     | F        |
!! | iter           | ccpp_loop_counter                                                           | loop counter for subcycling loops in CCPP             | index         |    0 | integer     |           | in     | F        |
!! | flag_iter      | flag_for_iteration                                                          | flag for iteration                                    | flag          |    1 | logical     |           | in     | F        |
!! | ugrs           | x_wind                                                                      | x component of layer wind                             | m s-1         |    2 | real        | kind_phys | in     | F        |
!! | vgrs           | y_wind                                                                      | y component of layer wind                             | m s-1         |    2 | real        | kind_phys | in     | F        |
!! | tgrs           | air_temperature                                                             | layer mean air temperature                            | K             |    2 | real        | kind_phys | in     | F        |
!! | qgrs           | tracer_concentration                                                        | model layer mean tracer concentration                 | kg kg-1       |    3 | real        | kind_phys | in     | F        |
!! | prsl           | air_pressure                                                                | mean layer pressure                                   | Pa            |    2 | real        | kind_phys | in     | F        |
!! | prsi           | air_pressure_at_interface                                                   | air pressure at model layer interfaces                | Pa            |    2 | real        | kind_phys | in     | F        |
!! | phii           | geopotential_at_interface                                                   | geopotential at model layer interfaces                | m2 s-2        |    2 | real        | kind_phys | in     | F        |
!! | prsik_1        | dimensionless_exner_function_at_lowest_model_interface                      | dimensionless Exner function at lowest model interface| none          |    1 | real        | kind_phys | in     | F        |
!! | prslk_1        | dimensionless_exner_function_at_lowest_model_layer                          | dimensionless Exner function at lowest model layer    | none          |    1 | real        | kind_phys | in     | F        |
!! | tsfc           | surface_skin_temperature                                                    | surface temperature                                   | K             |    1 | real        | kind_phys | in     | F        |
!! | qsfc           | surface_specific_humidity                                                   | surface air saturation specific humidity              | kg kg-1       |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_qsfc   | surface_specific_humidity_for_MYJ_schemes                                   | surface air saturation specific humidity for MYJ schem| kg kg-1       |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_thz0   | potential_temperature_at_viscous_sublayer_top                               | potential temperat at viscous sublayer top over water | K             |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_qz0    | specific_humidity_at_viscous_sublayer_top                                   | specific humidity at_viscous sublayer top over water  | kg kg-1       |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_uz0    | u_wind_component_at_viscous_sublayer_top                                    | u wind component at viscous sublayer top over water   | m s-1         |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_vz0    | v_wind_component_at_viscous_sublayer_top                                    | v wind component at viscous sublayer top over water   | m s-1         |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_z0base | baseline_surface_roughness_length                                           | baseline surface roughness length for momentum in mete| m             |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_akhs   | heat_exchange_coefficient_for_MYJ_schemes                                   | surface heat exchange_coefficient for MYJ schemes     | m s-1         |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_akms   | momentum_exchange_coefficient_for_MYJ_schemes                               | surface momentum exchange_coefficient for MYJ schemes | m s-1         |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_chkqlm | surface_layer_evaporation_switch                                            | surface layer evaporation switch                      | none          |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_elflx  | kinematic_surface_latent_heat_flux                                          | kinematic surface latent heat flux                    | m s-1 kg kg-1 |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_a1u    | weight_for_momentum_at_viscous_sublayer_top                                 | Weight for momentum at viscous layer top              | none          |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_a1t    | weight_for_potental_temperature_at_viscous_sublayer_top                     | Weight for potental temperature at viscous layer top  | none          |    1 | real        | kind_phys | inout  | F        |
!! | phy_myj_a1q    | weight_for_specific_humidity_at_viscous_sublayer_top                        | Weight for Specfic Humidity at viscous layer top      | none          |    1 | real        | kind_phys | inout  | F        |
!! | pblh           | atmosphere_boundary_layer_thickness                                         | PBL thickness                                         | m             |    1 | real        | kind_phys | inout  | F        |
!! | slmsk          | sea_land_ice_mask_real                                                      | landmask: sea/land/ice=0/1/2                          | flag          |    1 | real        | kind_phys | in     | F        |
!! | zorl           | surface_roughness_length                                                    | surface roughness length                              | cm            |    1 | real        | kind_phys | in     | F        |
!! | ustar          | surface_friction_velocity                                                   | boundary layer parameter                              | m s-1         |    1 | real        | kind_phys | inout  | F        |
!! | rib            | bulk_richardson_number_at_lowest_model_level                                | bulk Richardson number at the surface                 | none          |    1 | real        | kind_phys | inout  | F        |
!! | cm             | surface_drag_coefficient_for_momentum_in_air                                | surface exchange coeff for momentum                   | none          |    1 | real        | kind_phys | inout  | F        |
!! | ch             | surface_drag_coefficient_for_heat_and_moisture_in_air                       | surface exchange coeff heat & moisture                | none          |    1 | real        | kind_phys | inout  | F        |
!! | stress         | surface_wind_stress                                                         | surface wind stress                                   | m2 s-2        |    1 | real        | kind_phys | in     | F        |
!! | ffm            | Monin_Obukhov_similarity_function_for_momentum                              | Monin_Obukhov similarity function for momentum        | none          |    1 | real        | kind_phys | inout  | F        |
!! | ffh            | Monin_Obukhov_similarity_function_for_heat                                  | Monin_Obukhov similarity function for heat            | none          |    1 | real        | kind_phys | inout  | F        |
!! | fm10           | Monin_Obukhov_similarity_function_for_momentum_at_10m                       | Monin_Obukhov similarity parameter for momentum at 10m| none          |    1 | real        | kind_phys | inout  | F        |
!! | fh2            | Monin_Obukhov_similarity_function_for_heat_at_2m                            | Monin_Obukhov similarity parameter for heat at 2m     | none          |    1 | real        | kind_phys | inout  | F        |
!! | landfrac       | land_area_fraction                                                          | fraction of horizontal grid area occupied by land     | frac          |    1 | real        | kind_phys | inout  | F        |
!! | lakefrac       | lake_area_fraction                                                          | fraction of horizontal grid area occupied by lake     | frac          |    1 | real        | kind_phys | inout  | F        |
!! | oceanfrac      | sea_area_fraction                                                           | fraction of horizontal grid area occupied by ocean    | frac          |    1 | real        | kind_phys | inout  | F        |
!! | fice           | sea_ice_concentration                                                       | ice fraction over open water                          | frac          |    1 | real        | kind_phys | in     | F        |
!! | z0rl_ocn       | surface_roughness_length_over_ocean_interstitial                            | surface roughness length over ocean (interstitial)    | cm            |    1 | real        | kind_phys | inout  | F        |
!! | z0rl_lnd       | surface_roughness_length_over_land_interstitial                             | surface roughness length over land  (interstitial)    | cm            |    1 | real        | kind_phys | inout  | F        |
!! | z0rl_ice       | surface_roughness_length_over_ice_interstitial                              | surface roughness length over ice   (interstitial)    | cm            |    1 | real        | kind_phys | inout  | F        |
!! | ustar_ocn      | surface_friction_velocity_over_ocean                                        | surface friction velocity over ocean                  | m s-1         |    1 | real        | kind_phys | inout  | F        |
!! | ustar_lnd      | surface_friction_velocity_over_land                                         | surface friction velocity over land                   | m s-1         |    1 | real        | kind_phys | inout  | F        |
!! | ustar_ice      | surface_friction_velocity_over_ice                                          | surface friction velocity over ice                    | m s-1         |    1 | real        | kind_phys | inout  | F        |
!! | cm_ocn         | surface_drag_coefficient_for_momentum_in_air_over_ocean                     | surface exchange coeff for momentum over ocean        | none          |    1 | real        | kind_phys | inout  | F        |
!! | cm_lnd         | surface_drag_coefficient_for_momentum_in_air_over_land                      | surface exchange coeff for momentum over land         | none          |    1 | real        | kind_phys | inout  | F        |
!! | cm_ice         | surface_drag_coefficient_for_momentum_in_air_over_ice                       | surface exchange coeff for momentum over ice          | none          |    1 | real        | kind_phys | inout  | F        |
!! | ch_ocn         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ocean            | surface exchange coeff heat & moisture over ocean     | none          |    1 | real        | kind_phys | inout  | F        |
!! | ch_lnd         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_land             | surface exchange coeff heat & moisture over land      | none          |    1 | real        | kind_phys | inout  | F        |
!! | ch_ice         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ice              | surface exchange coeff heat & moisture over ice       | none          |    1 | real        | kind_phys | inout  | F        |
!! | rb_ocn         | bulk_richardson_number_at_lowest_model_level_over_ocean                     | bulk Richardson number at the surface over ocean      | none          |    1 | real        | kind_phys | inout  | F        |
!! | rb_lnd         | bulk_richardson_number_at_lowest_model_level_over_land                      | bulk Richardson number at the surface over land       | none          |    1 | real        | kind_phys | inout  | F        |
!! | rb_ice         | bulk_richardson_number_at_lowest_model_level_over_ice                       | bulk Richardson number at the surface over ice        | none          |    1 | real        | kind_phys | inout  | F        |
!! | stress_ocn     | surface_wind_stress_over_ocean                                              | surface wind stress over ocean                        | m2 s-2        |    1 | real        | kind_phys | inout  | F        |
!! | stress_lnd     | surface_wind_stress_over_land                                               | surface wind stress over land                         | m2 s-2        |    1 | real        | kind_phys | inout  | F        |
!! | stress_ice     | surface_wind_stress_over_ice                                                | surface wind stress over ice                          | m2 s-2        |    1 | real        | kind_phys | inout  | F        |
!! | fm_ocn         | Monin_Obukhov_similarity_function_for_momentum_over_ocean                   | Monin-Obukhov similarity funct for momentum over ocean| none          |    1 | real        | kind_phys | inout  | F        |
!! | fm_lnd         | Monin_Obukhov_similarity_function_for_momentum_over_land                    | Monin-Obukhov similarity funct for momentum over land | none          |    1 | real        | kind_phys | inout  | F        |
!! | fm_ice         | Monin_Obukhov_similarity_function_for_momentum_over_ice                     | Monin-Obukhov similarity funct for momentum over ice  | none          |    1 | real        | kind_phys | inout  | F        |
!! | fh_ocn         | Monin_Obukhov_similarity_function_for_heat_over_ocean                       | Monin-Obukhov similarity function for heat over ocean | none          |    1 | real        | kind_phys | inout  | F        |
!! | fh_lnd         | Monin_Obukhov_similarity_function_for_heat_over_land                        | Monin-Obukhov similarity function for heat over land  | none          |    1 | real        | kind_phys | inout  | F        |
!! | fh_ice         | Monin_Obukhov_similarity_function_for_heat_over_ice                         | Monin-Obukhov similarity function for heat over ice   | none          |    1 | real        | kind_phys | inout  | F        |
!! | fm10_ocn       | Monin_Obukhov_similarity_function_for_momentum_at_10m_over_ocean            | Monin-Obukhov parameter for momentum at 10m over ocean| none          |    1 | real        | kind_phys | inout  | F        |
!! | fm10_lnd       | Monin_Obukhov_similarity_function_for_momentum_at_10m_over_land             | Monin-Obukhov parameter for momentum at 10m over land | none          |    1 | real        | kind_phys | inout  | F        |
!! | fm10_ice       | Monin_Obukhov_similarity_function_for_momentum_at_10m_over_ice              | Monin-Obukhov parameter for momentum at 10m over ice  | none          |    1 | real        | kind_phys | inout  | F        |
!! | fh2_ocn        | Monin_Obukhov_similarity_function_for_heat_at_2m_over_ocean                 | Monin-Obukhov parameter for heat at 2m over ocean     | none          |    1 | real        | kind_phys | inout  | F        |
!! | fh2_lnd        | Monin_Obukhov_similarity_function_for_heat_at_2m_over_land                  | Monin-Obukhov parameter for heat at 2m over land      | none          |    1 | real        | kind_phys | inout  | F        |
!! | fh2_ice        | Monin_Obukhov_similarity_function_for_heat_at_2m_over_ice                   | Monin-Obukhov parameter for heat at 2m over ice       | none          |    1 | real        | kind_phys | inout  | F        |
!! | wind           | wind_speed_at_lowest_model_layer                                            | wind speed at lowest model level                      | m s-1         |    1 | real        | kind_phys | out    | F        |
!! | con_cp                     | specific_heat_of_dry_air_at_constant_pressure                               | specific heat of dry air at constant pressure                                               | J kg-1 K-1    |    0 | real       | kind_phys | in     | F        |
!! | con_g                      | gravitational_acceleration                                                  | gravitational acceleration                                                                  | m s-2         |    0 | real       | kind_phys | in     | F        |
!! | con_rd                     | gas_constant_dry_air                                                        | ideal gas constant for dry air                                                              | J kg-1 K-1    |    0 | real       | kind_phys | in     | F        |
!! | me             | mpi_rank                                                                    | current MPI-rank                                      | index         |    0 | integer     |           | in     | F        |
!! | lprnt          | flag_print                                                                  | control flag for diagnostic print out                 | flag          |    0 | logical     |           | in     | F        |
!! | errmsg         | ccpp_error_message                                                          | error message for error handling in CCPP              | none          |    0 | character   | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                             | error flag for error handling in CCPP                 | flag          |    0 | integer     |           | out    | F        |
!!
#endif
!###===================================================================
 SUBROUTINE myjsfc_wrapper_run(                    &
     &  restart,                                   &
     &  ix,im,levs,                                &
     &  kdt,ntrac,ntke,                            &
     &  ntcw,ntiw,ntrw,ntsw,ntgl,                  &
     &  iter,flag_iter,                            &
     &  ugrs, vgrs, tgrs, qgrs,                    &
     &  prsl, prsi, phii,                          &
     &  prsik_1, prslk_1, tsfc, qsfc,              &
     &  phy_myj_qsfc, phy_myj_thz0, phy_myj_qz0,   &
     &  phy_myj_uz0, phy_myj_vz0, phy_myj_z0base,  &
     &  phy_myj_akhs, phy_myj_akms,                &
     &  phy_myj_chkqlm, phy_myj_elflx,             &
     &  phy_myj_a1u, phy_myj_a1t, phy_myj_a1q,     &
     &  pblh, slmsk, zorl, ustar, rib,             &
     &  cm,ch,stress,ffm,ffh,fm10,fh2,             &
     &  landfrac,lakefrac,oceanfrac,fice,          &
     &  z0rl_ocn,  z0rl_lnd,  z0rl_ice,            &   ! intent(inout)
     &  ustar_ocn, ustar_lnd, ustar_ice,           &   ! intent(inout)
     &  cm_ocn,    cm_lnd,    cm_ice,              &   ! intent(inout)
     &  ch_ocn,    ch_lnd,    ch_ice,              &   ! intent(inout)
     &  rb_ocn,    rb_lnd,    rb_ice,              &   ! intent(inout)
     &  stress_ocn,stress_lnd,stress_ice,          &   ! intent(inout)
     &  fm_ocn,    fm_lnd,    fm_ice,              &   ! intent(inout)
     &  fh_ocn,    fh_lnd,    fh_ice,              &   ! intent(inout)
     &  fm10_ocn,  fm10_lnd,  fm10_ice,            &   ! intent(inout)
     &  fh2_ocn,   fh2_lnd,   fh2_ice,             &   ! intent(inout)
     &  wind,      con_cp,    con_g,    con_rd,    &
     &  me, lprnt, errmsg, errflg )             ! intent(inout)
! 
      use machine,        only : kind_phys
      use MODULE_SF_JSFC, only: JSFC_INIT,JSFC

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 

      integer,parameter:: &
         klog=4 &                   ! logical variables
        ,kint=4 &                   ! integer variables
        ,kfpt=4 &                   ! floating point variables
        ,kdbl=8                     ! double precision
!
!  ---  constant parameters:
!      real(kind=kind_phys), parameter :: karman  = 0.4

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!For reference
!     real    , parameter :: karman       = 0.4
!     real    , parameter :: g            = 9.81
!     real    , parameter :: r_d          = 287.
!     real    , parameter :: cp           = 7.*r_d/2.
!     real    , parameter :: r_v          = 461.6
!     real    , parameter :: cpv          = 4.*r_v
!     real    , parameter :: rcp          = r_d/cp

!      real, parameter :: g_inv=1/g, cappa=r_d/cp

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

!MYJ-1D
      integer,intent(in) :: im, ix, levs 
      integer,intent(in) :: kdt, iter, me
      integer,intent(in) :: ntrac,ntke,ntcw,ntiw,ntrw,ntsw,ntgl
      logical,intent(in) :: restart, lprnt
      real(kind=kind_phys),intent(in) :: con_cp, con_g, con_rd

!MYJ-2D
      logical,dimension(im),intent(in) :: flag_iter
      real(kind=kind_phys),dimension(im),intent(in)      ::  &
     &        prsik_1, prslk_1, tsfc, qsfc, slmsk
      real(kind=kind_phys),dimension(im),intent(inout)   ::  &
     &        phy_myj_qsfc, phy_myj_thz0, phy_myj_qz0,       &
     &        phy_myj_uz0, phy_myj_vz0, phy_myj_z0base,      &
     &        phy_myj_akhs, phy_myj_akms,                    &
     &        phy_myj_chkqlm, phy_myj_elflx,                 &
     &        phy_myj_a1u, phy_myj_a1t, phy_myj_a1q
      real(kind=kind_phys),dimension(im),intent(inout)   ::  &
     &        pblh, zorl, ustar, rib
      real(kind=kind_phys),dimension(im),intent(out)     ::  &
     &        cm, ch, stress, ffm, ffh, fm10, fh2
      real(kind=kind_phys), dimension(im), intent(inout) ::  &
     &        landfrac, lakefrac, oceanfrac, fice
      real(kind=kind_phys), dimension(im), intent(inout) ::  &
     &                    z0rl_ocn,  z0rl_lnd,  z0rl_ice,    &
     &                   ustar_ocn, ustar_lnd, ustar_ice,    &
     &                      cm_ocn,    cm_lnd,    cm_ice,    &
     &                      ch_ocn,    ch_lnd,    ch_ice,    &
     &                      rb_ocn,    rb_lnd,    rb_ice,    &
     &                  stress_ocn,stress_lnd,stress_ice,    &
     &                      fm_ocn,    fm_lnd,    fm_ice,    &
     &                      fh_ocn,    fh_lnd,    fh_ice,    &
     &                    fm10_ocn,  fm10_lnd,  fm10_ice,    &
     &                     fh2_ocn,   fh2_lnd,   fh2_ice,    &
     &                      wind


!MYJ-3D
      real(kind=kind_phys),dimension(im,levs+1),intent(in) ::  &
              phii, prsi
      real(kind=kind_phys),dimension(im,levs),intent(in)   ::  &
     &        ugrs, vgrs, tgrs, prsl
!MYJ-4D
      real(kind=kind_phys),dimension(im,levs,ntrac),intent(in) :: &
     &       qgrs

!LOCAL
      logical :: lprnt1, lprnt2
      integer :: ntsd, k, k1, i, n, ide, jde, kde

      real(kind=kind_phys) :: g, r_d, g_inv, cappa
      real(kind=kfpt),dimension(levs)       :: epsq2
      real(kind=kfpt),dimension(im)         ::           &
           sfcz,tsk,xland,mavail,rmol,                   &
           ustar1,z0,rib1,sm,pblh_myj
      real(kind=kfpt),dimension(im,13)      ::           &
     &     phy_f2d_myj
      real(kind=kfpt), dimension(im,levs)   ::           &
     &     u_myj, v_myj, t_myj, q_myj, th_myj,           &
     &     cw, dz_myj, pmid, q2, exner
      real(kind=kfpt), dimension(im,levs+1) :: pint
      real(kind=kfpt),dimension(im)         ::           &
     &     cm1,ch1,stress1,ffm1,ffh1,wind1,ffm10,ffh2                    
!      real(kind=kind_phys), dimension(im,levs,ntrac) :: &
!     &     qgrs_myj  

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ntsd = kdt-1

      lprnt1 =.false.
      lprnt2 =.false.

      if (lprnt2) then
         if(me.eq.0)then
           print*,'in myj surface layer wrapper...'
           print*,'ntsd,iter=',ntsd,iter
         end if
      endif

      r_d   = con_rd
      g     = con_g
      g_inv = 1./con_g
      cappa = con_rd/con_cp

      if (ntsd==0.and.iter==1)then
        do i=1,im
           if(flag_iter(i))then
              phy_myj_qsfc(i)   = qgrs(i,1,1)     ! qsfc(:)
              phy_myj_thz0(i)   = tsfc(i)         ! thz0
              phy_myj_qz0(i)    = qgrs(i,1,1)     ! qz0(:)
              phy_myj_uz0(i)    = 0.              ! uz0(:)
              phy_myj_vz0(i)    = 0.              ! vz0(:)
              phy_myj_z0base(i) = zorl(i)*0.01    ! z0base
              phy_myj_akhs(i)   = 0.01            ! akhs(:)
              phy_myj_akms(i)   = 0.01            ! akms(:)
              phy_myj_chkqlm(i) = 0.              ! chkqlm(:)
              phy_myj_elflx(i)  = 0.              ! elflx(:)
              phy_myj_a1u(i)    = 0.              ! a1u
              phy_myj_a1t(i)    = 0.              ! a1t
              phy_myj_a1q(i)    = 0.              ! a1q
           end if
        end do
      end if

!prep MYJ-only variables
      do i=1,im
         sm(i)=1.; if(slmsk(i) > 0.5 ) sm(i)=0.
         xland(i)=sm(i)+1.
         sfcz(i)=phii(i,1)*g_inv
      enddo

      do k=1,levs
         k1=levs+1-k
         do i=1,im
            u_myj(i,k)=ugrs(i,k1)
            v_myj(i,k)=vgrs(i,k1)
            t_myj(i,k)=tgrs(i,k1)
            q_myj(i,k)=qgrs(i,k1,1)
            cw(i,k)   =qgrs(i,k1,ntcw)
!            if(ntrw.gt.0)cw(i,k) = cw(i,k) + qgrs(i,k1,ntrw)
!            if(ntiw.gt.0)cw(i,k) = cw(i,k) + qgrs(i,k1,ntiw)
!            if(ntsw.gt.0)cw(i,k) = cw(i,k) + qgrs(i,k1,ntsw)
!            if(ntgl.gt.0)cw(i,k) = cw(i,k) + qgrs(i,k1,ntgl)
            if(ntke.gt.0)then
              q2(i,k) =qgrs(i,k1,ntke)*2.
            else
              q2(i,k) =0.02
            end if
            pmid(i,k) =prsl(i,k1)
            exner(i,k)=(prsl(i,k1)*1.e-5)**cappa
            th_myj(i,k)=tgrs(i,k1)/exner(i,k)
         end do
      end do
      do k=1,levs+1
         k1=levs+2-k
         do i=1,im
            pint(i,k) =prsi(i,k1)
         end do
      end do

      do k = 1, levs
         k1 = levs-k+1
         do i = 1, im
            dz_myj(i,k) = (phii(i,k1+1)-phii(i,k1)) * g_inv
         enddo
      enddo
         
      if (lprnt1) then
         if(me==0.and.ntsd.lt.2)then
            k=63
            k1=levs+1-k
            print*,'Qingfu starts MYJSFC'
            print*,'ntsd,iter,me,1=',ntsd,iter,me
            print*,'ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,ntke=',   &
                    ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,ntke
            print*,'im,levs,ntsd=',im,levs,ntsd
            do i=10,40,40
            print*,'Qingfu before MYJ surface kdt,i,k1=',kdt,i,k1
            print*,'sfcz,dz_myj,th_myj,tsfc,qsfc=',sfcz(i),dz_myj(i,k),   &
                    th_myj(i,k),tsfc(i),qsfc(i)
             print*,'sm,z0,xland=',                          &
                   sm(i),z0(i),xland(i)
!            print*,'phy_f2d_myj(i,1:13)=',               &
!                   (phy_f2d_myj(i,n),n=1,13)
            print*,'u_myj,v_myj=',                           &
                   u_myj(i,k),v_myj(i,k)
            print*,'t_myj,q_myj,cw,q2=',                     &
                   t_myj(i,k),q_myj(i,k),cw(i,k),q2(i,k)
            print*,'phii,pint,pmid',                         &
                    phii(i,k1),pint(i,k),pmid(i,k)
            print*,'exner,th_myj=',exner(i,k),th_myj(i,k)
            end do
         end if
      endif

!-----------------------------------------------------------------------
      ide=im+1
      jde=2
      kde=levs+1

      do i = 1, im
         epsq2(i)=0.02
         mavail(i)=1.0
         tsk(i)=tsfc(i)
         phy_f2d_myj(i,1)  = phy_myj_qsfc(i)
         phy_f2d_myj(i,2)  = phy_myj_thz0(i)
         phy_f2d_myj(i,3)  = phy_myj_qz0(i)
         phy_f2d_myj(i,4)  = phy_myj_uz0(i)
         phy_f2d_myj(i,5)  = phy_myj_vz0(i)
         phy_f2d_myj(i,6)  = phy_myj_z0base(i)
         phy_f2d_myj(i,7)  = phy_myj_akhs(i)
         phy_f2d_myj(i,8)  = phy_myj_akms(i)
         phy_f2d_myj(i,9)  = phy_myj_chkqlm(i)
         phy_f2d_myj(i,10) = phy_myj_elflx(i)
         phy_f2d_myj(i,11) = phy_myj_a1u(i)
         phy_f2d_myj(i,12) = phy_myj_a1t(i)
         phy_f2d_myj(i,13) = phy_myj_a1q(i)
         z0(i)=zorl(i)*0.01
         rmol(i)=0.
         rib1(I)=rib(i)
         pblh_myj(i)=pblh(i)
         ustar1(i)=ustar(i)
         cm1(i)=0.
         ch1(i)=0.
         stress1(i)=0.
         ffm1(i)=0.
         ffh1(i)=0.
         wind1(i)=0.
         ffm10(i)=0.
         ffh2(i)=0.
      end do

      if((ntsd==0.and.iter.eq.1).or.restart)then
         call JSFC_INIT(ustar1,restart                       &
     &                 ,1,ide,1,jde,1,kde                    &
     &                 ,1,im,1,1,1,levs                      &
     &                 ,1,im,1,1,1,levs)
      end if
               
      call JSFC(flag_iter,iter,me                            &
     &         ,ntsd,epsq2,sfcz,dz_myj                       &
     &         ,pmid,pint,th_myj,t_myj,q_myj,cw              &
     &         ,u_myj,v_myj,q2,tsk                           &
     &         ,phy_f2d_myj(1:im,1),phy_f2d_myj(1:im,2)      &
     &         ,phy_f2d_myj(1:im,3),phy_f2d_myj(1:im,4)      &
     &         ,phy_f2d_myj(1:im,5),xland                    &
     &         ,ustar1,z0,phy_f2d_myj(1:im,6)                &
     &         ,pblh_myj,mavail,rmol                         &
     &         ,phy_f2d_myj(1:im,7),phy_f2d_myj(1:im,8)      &
     &         ,phy_f2d_myj(1:im,9),phy_f2d_myj(1:im,10)     &
     &         ,rib1,cm1,ch1,stress1,ffm1,ffh1,wind1,ffm10,ffh2  &
     &         ,phy_f2d_myj(1:im,11),phy_f2d_myj(1:im,12)    &
     &         ,phy_f2d_myj(1:im,13)                         &
     &         ,1,im,1,1,1,levs                              &
     &         ,1,im,1,1,1,levs                              &
     &         ,1,im,1,1,1,levs)

      do i = 1, im
         if(flag_iter(i))then
            zorl(i) = z0(i)*100.

            phy_myj_qsfc(i)   = phy_f2d_myj(i,1)
            phy_myj_thz0(i)   = phy_f2d_myj(i,2)
            phy_myj_qz0(i)    = phy_f2d_myj(i,3)
            phy_myj_uz0(i)    = phy_f2d_myj(i,4)
            phy_myj_vz0(i)    = phy_f2d_myj(i,5)
            phy_myj_z0base(i) = phy_f2d_myj(i,6)
            phy_myj_akhs(i)   = phy_f2d_myj(i,7)
            phy_myj_akms(i)   = phy_f2d_myj(i,8)
            phy_myj_chkqlm(i) = phy_f2d_myj(i,9)
            phy_myj_elflx(i)  = - phy_f2d_myj(i,10)    ! change flux definition
            phy_myj_a1u(i)    = phy_f2d_myj(i,11)
            phy_myj_a1t(i)    = phy_f2d_myj(i,12)
            phy_myj_a1q(i)    = phy_f2d_myj(i,13)

            rib(I)=rib1(i)
            pblh(I)=pblh_myj(i)
            cm(I)=cm1(i)
            ch(I)=ch1(i)
            stress(I)=stress1(i)
            ffm(I)=ffm1(i)
            ffh(I)=ffh1(i)
            wind(I)=wind1(i)
            fm10(I)=ffm10(i)
            fh2(I)=ffh2(i)
            ustar(i)=ustar1(i)
         end if
      end do

      if (lprnt1) then

        if(me==0.and.ntsd.lt.10)then
           print*,'ntsd,iter,me,2=',ntsd,iter,me
           do i=10,40,40
             if(flag_iter(i))then
               print*,'Qingfu after MYJ surface kdt,i,k1=',kdt,i,k1
               print*,'xland,cm,ch=',xland(i),cm(i),ch(i)
               print*,'ustar,z0,stress=',ustar(i),z0(i),stress(i)
               print*,'ffm,ffh,wind,fm10,fh2=',ffm(i),ffh(i),wind(i),fm10(i),fh2(i)
               print*,'phy_f2d_myj(9,1:13)=',    &
                     (phy_f2d_myj(i,n),n=1,13)
               print*,'u_myj,v_myj=',  &
                     u_myj(i,k),v_myj(i,k)
               print*,'t_myj,q_myj,cw,q2=',  &
                     t_myj(i,k),q_myj(i,k),cw(i,k),q2(i,k)
               print*,'phii,pint,pmid',  &
                     phii(i,k1),pint(i,k),pmid(i,k)
               print*,'exner,th_myj=',exner(i,k),th_myj(i,k)
               print*,'Qingfu finish MYJSFC'
            end if
          end do
        end if

        do k=1,levs
           k1=levs+1-k
           do i=1,im
             if(t_myj(i,k).gt.320..or.t_myj(i,k).lt.150.)then
                print*,'xland,cm,ch=',xland(i),cm(i),ch(i)
                print*,'ustar,z0,stress=',ustar(i),z0(i),stress(i)
                print*,'ffm,ffh,wind,fm10,fh2=',ffm(i),ffh(i),wind(i),fm10(i),fh2(i)
                print*,'phy_f2d_myj(9,1:13)=',    &
                      (phy_f2d_myj(i,n),n=1,13)
                print*,'u_myj,v_myj=',  &
                      u_myj(i,k),v_myj(i,k)
                print*,'t_myj,q_myj,cw,q2=',  &
                      t_myj(i,k),q_myj(i,k),cw(i,k),q2(i,k)
                print*,'phii,pint,pmid',  &
                      phii(i,k1),pint(i,k),pmid(i,k)
                print*,'exner,th_myj=',exner(i,k),th_myj(i,k)
                print*,'Qingfu finish MYJSFC'
             end if
           end do
        end do

      end if

      do i = 1, im
         if(flag_iter(i))then
                z0rl_ocn(i) = zorl(i)
                  cm_ocn(i) = cm(i)
                  ch_ocn(i) = ch(i)
                  rb_ocn(i) = rib(i)
              stress_ocn(i) = stress(i)
                  fm_ocn(i) = ffm(i)
                  fh_ocn(i) = ffh(i)
               ustar_ocn(i) = ustar(i)
                fm10_ocn(i) = fm10(i)
                 fh2_ocn(i) = fh2(i)

                z0rl_lnd(i) = zorl(i)
                  cm_lnd(i) = cm(i)
                  ch_lnd(i) = ch(i)
                  rb_lnd(i) = rib(i)
              stress_lnd(i) = stress(i)
                  fm_lnd(i) = ffm(i)
                  fh_lnd(i) = ffh(i)
               ustar_lnd(i) = ustar(i)
                fm10_lnd(i) = fm10(i)
                 fh2_lnd(i) = fh2(i)

                z0rl_ice(i) = zorl(i)
                  cm_ice(i) = cm(i)
                  ch_ice(i) = ch(i)
                  rb_ice(i) = rib(i)
              stress_ice(i) = stress(i)
                  fm_ice(i) = ffm(i)
                  fh_ice(i) = ffh(i)
               ustar_ice(i) = ustar(i)
                fm10_ice(i) = fm10(i)
                 fh2_ice(i) = fh2(i)
            end if
      end do

      if (lprnt2) then
        if(me==0.and.ntsd.lt.10)then
          print*,'ntsd,iter,me,3=',ntsd,iter,me
          do i=10,40,40
             if(flag_iter(i))then
               print*,'Qingfu after MYJ surface kdt,i,k1,3=',kdt,i,k1
               print*,'Qingfu test after MYJ surface kdt,i=',kdt,i,slmsk(i)
               print*,'a1u,a1t,a1q=',(phy_f2d_myj(i,k),k=11,13)
               print*,'zorl,cm,ch,rb,stress=',z0(i),    &
                       cm(i),ch(i),   &
                       rib(i),stress(i)
               print*,'ffmm,ffhh,ustar,fm10,fh2,wind=', ffm(i), &
                       ffh(i),ustar(i),fm10(i),fh2(i),wind(i)
               print*,'cm(i),ch(i)=',    &
                       (0.4/ffm(i))**2,(0.4/ffm(i)*0.4/ffh(i))
             end if
          end do
        endif
      endif


  END SUBROUTINE myjsfc_wrapper_run

!###=================================================================

END MODULE myjsfc_wrapper

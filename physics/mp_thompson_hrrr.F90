! CCPP license goes here, as well as further documentation
module mp_thompson_hrrr

      use machine, only : kind_phys

      use module_mp_thompson_hrrr, only : thompson_init, mp_gt_driver

      implicit none

   contains

#if 0
!! \section arg_table_mp_thompson_hrrr_init Argument Table
!! | local_name      | standard_name                                                 | long_name                                              | units      | rank | type      |    kind   | intent | optional |
!! |-----------------|---------------------------------------------------------------|--------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ncol            | horizontal_loop_extent                                        | horizontal loop extent                                 | count      |    0 | integer   |           | in     | F        |
!! | nlev            | vertical_dimension                                            | number of vertical levels                              | count      |    0 | integer   |           | in     | F        |
!! | con_g           | gravitational_acceleration                                    | gravitational acceleration                             | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | con_rd          | gas_constant_dry_air                                          | ideal gas constant for dry air                         | J kg-1 K-1 |    0 | real      | kind_phys | in     | F        |
!! | phil            | geopotential                                                  | geopotential at model layer centers                    | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | prsl            | air_pressure                                                  | mean layer pressure                                    | Pa         |    2 | real      | kind_phys | in     | F        |
!! | tgrs            | air_temperature                                               | model layer mean temperature                           | K          |    2 | real      | kind_phys | in     | F        |
!! | is_aerosol_aware| flag_for_aerosol_physics                                      | flag for aerosol-aware physics                         | flag       |    0 | logical   |           | in     | F        |
!! | nwfa2d          | tendency_of_water_friendly_surface_aerosols_at_surface        | instantaneous fake surface aerosol source              | kg-1 s-1   |    1 | real      | kind_phys | inout  | T        |
!! | nwfa            | water_friendly_aerosol_number_concentration                   | number concentration of water-friendly aerosols        | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | nifa            | ice_friendly_aerosol_number_concentration                     | number concentration of ice-friendly aerosols          | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | area            | cell_area                                                     | area of the grid cell                                  | m2         |    1 | real      | kind_phys | in     | F        |
!! | errmsg          | error_message                                                 | error message for error handling in CCPP               | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | error_flag                                                    | error flag for error handling in CCPP                  | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine mp_thompson_hrrr_init(ncol, nlev, con_g, con_rd,  &
                                       phil, prsl, tgrs,           &
                                       is_aerosol_aware,           &
                                       nwfa2d, nwfa, nifa,         &
                                       area, errmsg, errflg)

         implicit none

         ! Interface variables
          ! DH* for all 2-d vars and all 1-d vars running over 1:ncol - third dim 1:1 ok?
         integer,                   intent(in)    :: ncol
         integer,                   intent(in)    :: nlev

         real(kind_phys),           intent(in)    :: con_g
         real(kind_phys),           intent(in)    :: con_rd
         real(kind_phys),           intent(in)    :: phil(1:ncol,1:nlev)
         real(kind_phys),           intent(in)    :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in)    :: tgrs(1:ncol,1:nlev)
         logical,                   intent(in)    :: is_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nwfa2d(1:ncol)
         real(kind_phys), optional, intent(inout) :: nwfa(1:ncol,1:nlev) ! in kg kg-1
         real(kind_phys), optional, intent(inout) :: nifa(1:ncol,1:nlev) ! in kg kg-1
         real(kind_phys),           intent(in)    :: area(1:ncol)
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables
         real (kind=kind_phys) :: hgt(1:ncol,1:nlev)
         real (kind=kind_phys) :: rho(1:ncol,1:nlev)
         real(kind_phys)       :: nwfa_mp(1:ncol,1:nlev) ! in particles per cm3
         real(kind_phys)       :: nifa_mp(1:ncol,1:nlev) ! in particles per cm3
         real(kind_phys)       :: dx(1:ncol)
         logical, parameter    :: is_start = .true.
         ! Dimensions used in thompson_init
         integer               :: ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Geopotential height in m2 s-2 to height in m
         hgt = phil/con_g

         ! Set internal dimensions
         ids = 1
         ims = 1
         its = 1
         ide = ncol
         ime = ncol
         ite = ncol
         jds = 1
         jms = 1
         jts = 1
         jde = 1
         jme = 1
         jte = 1
         kds = 1
         kms = 1
         kts = 1
         kde = nlev
         kme = nlev
         kte = nlev

         if (is_aerosol_aware .and. present(nwfa2d) .and. present(nwfa) .and. present(nifa)) then
            ! Density of air in kg m-3
            rho = prsl/(con_rd*tgrs)
            ! Aerosol number densities [cm-3] = number concentration [kg-1]
            !                                   * air density [kg m-3]
            !                                   * 1E-6 [m3 cm-3]
            nwfa_mp = nwfa*rho*1.0E-6_kind_phys
            nifa_mp = nifa*rho*1.0E-6_kind_phys
            ! Grid cell size
            dx = sqrt(area)
            ! Call init
            call thompson_init(hgt=hgt,                                              &
                               nwfa2d=nwfa2d, nwfa=nwfa_mp, nifa=nifa_mp,            &
                               dx=dx, dy=dx, is_start=is_start,                      &
                               ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                               its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
            if (errflg /= 0) return
            ! Calculate aerosol concentrations as inverse of above
            nwfa = nwfa_mp/rho*1.0E6_kind_phys
            nifa = nwfa_mp/rho*1.0E6_kind_phys
         else if (is_aerosol_aware) then
            write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_hrrr_init:',  &
                                       ' aerosol-aware microphysics require all of the', &
                                       ' following optional arguments: nwfa2d, nwfa, nifa'
            errflg = 1
            return
         else
            call thompson_init(hgt=hgt,                                              &
                               ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                               its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
            if (errflg /= 0) return
         end if

      end subroutine mp_thompson_hrrr_init

#if 0
!! \section arg_table_mp_thompson_hrrr_run Argument Table
!! | local_name      | standard_name                                                 | long_name                                              | units      | rank | type      |    kind   | intent | optional |
!! |-----------------|---------------------------------------------------------------|--------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ncol            | horizontal_loop_extent                                        | horizontal loop extent                                 | count      |    0 | integer   |           | in     | F        |
!! | nlev            | vertical_dimension                                            | number of vertical levels                              | count      |    0 | integer   |           | in     | F        |
!! | con_g           | gravitational_acceleration                                    | gravitational acceleration                             | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | con_rd          | gas_constant_dry_air                                          | ideal gas constant for dry air                         | J kg-1 K-1 |    0 | real      | kind_phys | in     | F        |
!! | spechum         | water_vapor_specific_humidity_updated_by_physics              | water vapor specific humidity                          | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qc              | cloud_condensed_water_mixing_ratio_updated_by_physics         | cloud water mixing ratio wrt dry+vapor (no condensates)| kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qr              | rain_water_mixing_ratio_updated_by_physics                    | rain water mixing ratio wrt dry+vapor (no condensates) | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qi              | ice_water_mixing_ratio_updated_by_physics                     | ice water mixing ratio wrt dry+vapor (no condensates)  | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qs              | snow_water_mixing_ratio_updated_by_physics                    | snow water mixing ratio wrt dry+vapor (no condensates) | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qg              | graupel_mixing_ratio_updated_by_physics                       | graupel mixing ratio wrt dry+vapor (no condensates)    | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | ni              | ice_number_concentration_updated_by_physics                   | ice number concentration                               | kg-1       |    2 | real      | kind_phys | inout  | F        |
!! | nr              | rain_number_concentration_updated_by_physics                  | rain number concentration                              | kg-1       |    2 | real      | kind_phys | inout  | F        |
!! | is_aerosol_aware| flag_for_aerosol_physics                                      | flag for aerosol-aware physics                         | flag       |    0 | logical   |           | in     | F        |
!! | nc              | cloud_droplet_number_concentration_updated_by_physics         | cloud droplet number concentration                     | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | nwfa            | water_friendly_aerosol_number_concentration_updated_by_physics| number concentration of water-friendly aerosols        | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | nifa            | ice_friendly_aerosol_number_concentration_updated_by_physics  | number concentration of ice-friendly aerosols          | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | nwfa2d          | tendency_of_water_friendly_surface_aerosols_at_surface        | instantaneous fake surface aerosol source              | kg-1 s-1   |    1 | real      | kind_phys | in     | T        |
!! | tgrs            | air_temperature_updated_by_physics                            | model layer mean temperature                           | K          |    2 | real      | kind_phys | inout  | F        |
!! | prsl            | air_pressure                                                  | mean layer pressure                                    | Pa         |    2 | real      | kind_phys | in     | F        |
!! | phii            | geopotential_at_interface                                     | geopotential at model layer interfaces                 | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | omega           | omega                                                         | layer mean vertical velocity                           | Pa s-1     |    2 | real      | kind_phys | in     | F        |
!! | dtp             | time_step_for_physics                                         | physics timestep                                       | s          |    0 | real      | kind_phys | in     | F        |
!! | kdt             | index_of_time_step                                            | current forecast iteration                             | index      |    0 | integer   |           | in     | F        |
!! | rain            | lwe_thickness_of_precipitation_amount_on_dynamics_timestep    | total rain at this time step                           | m          |    1 | real      | kind_phys | inout  | F        |
!! | rainst          | lwe_thickness_of_stratiform_precipitation_amount              | stratiform rainfall amount on physics timestep         | m          |    1 | real      | kind_phys | inout  | F        |
!! | snow            | lwe_thickness_of_snow_amount_on_dynamics_timestep             | snow fall at this time step                            | m          |    1 | real      | kind_phys | inout  | F        |
!! | graupel         | lwe_thickness_of_graupel_amount_on_dynamics_timestep          | graupel fall at this time step                         | m          |    1 | real      | kind_phys | inout  | F        |
!! | sr              | ratio_of_snowfall_to_rainfall                                 | ratio of snowfall to large-scale rainfall              | frac       |    1 | real      | kind_phys | out    | F        |
!! | islmsk          | sea_land_ice_mask                                             | sea/land/ice mask (=0/1/2)                             | flag       |    1 | integer   |           | in     | F        |
!! | refl_10cm       | radar_reflectivity_10cm                                       | instantaneous refl_10cm                                | dBZ        |    2 | real      | kind_phys | out    | F        |
!! | do_radar_ref    | flag_for_radar_reflectivity                                   | flag for radar reflectivity                            | flag       |    0 | logical   |           | in     | F        |
!! | re_cloud        | mean_effective_radius_for_liquid_cloud                        | mean effective radius for liquid cloud                 | micron     |    2 | real      | kind_phys | inout  | T        |
!! | re_ice          | mean_effective_radius_for_ice_cloud                           | mean effective radius for ice cloud                    | micron     |    2 | real      | kind_phys | inout  | T        |
!! | re_snow         | mean_effective_radius_for_snow_flake                          | mean effective radius for snow flake                   | micron     |    2 | real      | kind_phys | inout  | T        |
!! | errmsg          | error_message                                                 | error message for error handling in CCPP               | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | error_flag                                                    | error flag for error handling in CCPP                  | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine mp_thompson_hrrr_run(ncol, nlev, con_g, con_rd,         &
                              spechum, qc, qr, qi, qs, qg, ni, nr,       &
                              is_aerosol_aware, nc, nwfa, nifa, nwfa2d,  &
                              tgrs, prsl, phii, omega, dtp, kdt,         &
                              rain, rainst, snow, graupel, sr,           &
                              islmsk,                                    &
                              refl_10cm, do_radar_ref,                   &
                              re_cloud, re_ice, re_snow,                 &
                              errmsg, errflg)

         implicit none

         ! Interface variables
         ! DH* for all 2-d vars and all 1-d vars running over 1:ncol - third dim 1:1 ok?

         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         real(kind_phys),           intent(in   ) :: con_g
         real(kind_phys),           intent(in   ) :: con_rd
         ! Hydrometeors
         real(kind_phys),           intent(inout) :: spechum(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qc(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qr(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qi(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qs(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qg(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: ni(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nr(1:ncol,1:nlev)
         ! Aerosols
         logical,                   intent(in)    :: is_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nc(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nwfa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nifa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(in   ) :: nwfa2d(1:ncol)
         ! State variables and timestep information
         real(kind_phys),           intent(inout) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: phii(1:ncol,1:nlev+1)
         real(kind_phys),           intent(in   ) :: omega(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: dtp
         integer,                   intent(in   ) :: kdt
         ! Rain/snow/graupel fall amounts and fraction of frozen precip
         real(kind_phys),           intent(inout) :: rain(1:ncol)
         real(kind_phys),           intent(inout) :: rainst(1:ncol)
         real(kind_phys),           intent(inout) :: snow(1:ncol)
         real(kind_phys),           intent(inout) :: graupel(1:ncol)
         real(kind_phys),           intent(  out) :: sr(1:ncol)
         ! Sea/land/ice mask (currently not used)
         integer,                   intent(in   ) :: islmsk(1:ncol)
         ! Radar reflectivity
         real(kind_phys),           intent(  out) :: refl_10cm(1:ncol,1:nlev)
         logical,         optional, intent(in   ) :: do_radar_ref
         ! Cloud effective radii
         real(kind_phys), optional, intent(inout) :: re_cloud(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: re_ice(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: re_snow(1:ncol,1:nlev)
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables

         ! Air density
         real(kind_phys) :: rho(1:ncol,1:nlev)              ! kg m-3
         ! Hydrometeors
         real(kind_phys) :: qv_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qc_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qr_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qi_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qs_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qg_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         ! Aerosols
         real(kind_phys) :: nwfa_mp(1:ncol,1:nlev)          ! cm-3
         real(kind_phys) :: nifa_mp(1:ncol,1:nlev)          ! cm-3
         ! Vertical velocity and level width
         real(kind_phys) :: w(1:ncol,1:nlev)                ! m s-1
         real(kind_phys) :: dz(1:ncol,1:nlev)               ! m
         ! Rain/snow/graupel fall amounts
         real(kind_phys) :: rain_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: snow_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: graupel_mp(1:ncol)              ! mm, dummy, not used
         real(kind_phys) :: delta_rain_mp(1:ncol)           ! mm
         real(kind_phys) :: delta_snow_mp(1:ncol)           ! mm
         real(kind_phys) :: delta_graupel_mp(1:ncol)        ! mm
         ! Radar reflectivity
         real(kind_phys) :: vt_dbz_wt(1:ncol,1:nlev)        ! dummy for reflectivity-weighted terminal velocity, which is
                                                            ! in argument list of mp_gt_driver but is not calculated
         logical         :: diagflag                        ! must be true if do_radar_ref is true, not used otherwise
         integer         :: do_radar_ref_mp                 ! integer instead of logical do_radar_ref
         ! Effective cloud radii
         logical         :: do_effective_radii
         real(kind_phys) :: re_cloud_mp(1:ncol,1:nlev)      ! m
         real(kind_phys) :: re_ice_mp(1:ncol,1:nlev)        ! m
         real(kind_phys) :: re_snow_mp(1:ncol,1:nlev)       ! m
         integer         :: has_reqc
         integer         :: has_reqi
         integer         :: has_reqs
         ! Dimensions used in mp_gt_driver
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Convert specific humidity/moist mixing ratios to dry mixing ratios
         qv_mp = spechum/(1.0_kind_phys-spechum)
         qc_mp = qc/(1.0_kind_phys-spechum)
         qr_mp = qr/(1.0_kind_phys-spechum)
         qi_mp = qi/(1.0_kind_phys-spechum)
         qs_mp = qs/(1.0_kind_phys-spechum)
         qg_mp = qg/(1.0_kind_phys-spechum)

         ! Density of air in kg m-3
         rho = prsl/(con_rd*tgrs)

         if (is_aerosol_aware .and. present(nc) .and. present(nwfa) .and. present(nifa) .and. present(nwfa2d)) then
            ! Aerosol number densities [cm-3] = number concentration [kg-1]
            !                                   * air density [kg m-3]
            !                                   * 1E-6 [m3 cm-3]
            nwfa_mp = nwfa*rho*1.0E-6_kind_phys
            nifa_mp = nifa*rho*1.0E-6_kind_phys
         else if (is_aerosol_aware) then
            write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_hrrr_run:',  &
                                       ' aerosol-aware microphysics require all of the', &
                                       ' following optional arguments: nc, nwfa, nifa, nwfa2d'
            errflg = 1
            return
         end if

         ! Convert omega in Pa s-1 to vertical velocity w in m s-1
         w = -omega/(rho*con_g)

         ! Layer width in m from geopotential in m2 s-2
         dz = (phii(:,2:nlev+1) - phii(:,1:nlev)) / con_g

         ! Accumulated values inside Thompson scheme, not used;
         ! only use delta and add to inout variables (different units);
         ! note that the deltas are initialized in mp_gt_driver
         rain_mp    = 0
         snow_mp    = 0
         graupel_mp = 0

         ! Flags for calculating radar reflectivity; diagflag is redundant
         if (do_radar_ref) then
             diagflag = .true.
             do_radar_ref_mp = 1
         else
             diagflag = .false.
             do_radar_ref_mp = 0
         end if

         if (present(re_cloud) .and. present(re_ice) .and. present(re_snow)) then
             do_effective_radii = .true.
             has_reqc = 1
             has_reqi = 1
             has_reqs = 1
             ! Convert micron to m
             re_cloud_mp = re_cloud*1.0E-6_kind_phys
             re_ice_mp   = re_ice*1.0E-6_kind_phys
             re_snow_mp  = re_snow*1.0E-6_kind_phys
         else if (.not.present(re_cloud) .and. .not.present(re_ice) .and. .not.present(re_snow)) then
             do_effective_radii = .false.
             has_reqc = 0
             has_reqi = 0
             has_reqs = 0
             ! Initialize to zero
             re_cloud_mp = 0
             re_ice_mp   = 0
             re_snow_mp  = 0
         else
             write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_hrrr_run:',  &
                                        ' all or none of the following optional', &
                                        ' arguments are required: re_cloud, re_ice, re_snow'
             errflg = 1
             return
         end if

         ! Set internal dimensions
         ids = 1
         ims = 1
         its = 1
         ide = ncol
         ime = ncol
         ite = ncol
         jds = 1
         jms = 1
         jts = 1
         jde = 1
         jme = 1
         jte = 1
         kds = 1
         kms = 1
         kts = 1
         kde = nlev
         kme = nlev
         kte = nlev

         ! Call Thompson MP with or without aerosols
         if (is_aerosol_aware) then
            call mp_gt_driver(qv=qv_mp, qc=qc_mp, qr=qr_mp, qi=qi_mp, qs=qs_mp, qg=qg_mp,    &
                              ni=ni, nr=nr, nc=nc,                                           &
                              nwfa=nwfa_mp, nifa=nifa_mp, nwfa2d=nwfa2d,                     &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtp, itimestep=kdt,         &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm, vt_dbz_wt=vt_dbz_wt,                      &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              re_cloud=re_cloud_mp, re_ice=re_ice_mp, re_snow=re_snow_mp,    &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

         else
            call mp_gt_driver(qv=qv_mp, qc=qc_mp, qr=qr_mp, qi=qi_mp, qs=qs_mp, qg=qg_mp,    &
                              ni=ni, nr=nr, nc=nc,                                           &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtp, itimestep=kdt,         &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm, vt_dbz_wt=vt_dbz_wt,                      &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              re_cloud=re_cloud_mp, re_ice=re_ice_mp, re_snow=re_snow_mp,    &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
         end if

         ! convert dry mixing ratios to specific humidity/moist mixing ratios
         spechum = qv_mp/(1.0_kind_phys+qv_mp)
         qc      = qc_mp/(1.0_kind_phys+qv_mp)
         qr      = qr_mp/(1.0_kind_phys+qv_mp)
         qi      = qi_mp/(1.0_kind_phys+qv_mp)
         qs      = qs_mp/(1.0_kind_phys+qv_mp)
         qg      = qg_mp/(1.0_kind_phys+qv_mp)

         if (is_aerosol_aware) then
            ! Calculate aerosol concentrations as inverse of above
            nwfa = nwfa_mp/rho*1.0E6_kind_phys
            nifa = nwfa_mp/rho*1.0E6_kind_phys
         end if

         ! Convert rainfall from mm to m and add deltas to inout variables
         rain    = rain    + delta_rain_mp/1000.0_kind_phys
         rainst  = rainst  + delta_rain_mp/1000.0_kind_phys
         snow    = snow    + delta_snow_mp/1000.0_kind_phys
         graupel = graupel + delta_graupel_mp/1000.0_kind_phys

         if (do_effective_radii) then
            ! Convert m to micron
            re_cloud = re_cloud_mp*1.0E6_kind_phys
            re_ice   = re_ice_mp*1.0E6_kind_phys
            re_snow  = re_snow_mp*1.0E6_kind_phys
         end if

      end subroutine mp_thompson_hrrr_run

      ! DH* do we need to deallocate stuff? which function to call in module_mp_thompson_hrrr.F90?
      subroutine mp_thompson_hrrr_finalize()
      end subroutine mp_thompson_hrrr_finalize

end module mp_thompson_hrrr

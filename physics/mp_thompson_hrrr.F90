! CCPP license goes here, as well as further documentation

!#define DEBUG_AEROSOLS

module mp_thompson_hrrr

      use machine, only : kind_phys

      use module_mp_thompson_hrrr, only : thompson_init, mp_gt_driver, thompson_finalize

      implicit none

      public :: mp_thompson_hrrr_init, mp_thompson_hrrr_run, mp_thompson_hrrr_finalize

      private

      logical :: is_initialized = .False.

   contains

#if 0
!! \section arg_table_mp_thompson_hrrr_init Argument Table
!! | local_name           | standard_name                                         | long_name                                                | units      | rank | type      |    kind   | intent | optional |
!! |----------------------|-------------------------------------------------------|----------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ncol                 | horizontal_loop_extent                                | horizontal loop extent                                   | count      |    0 | integer   |           | in     | F        |
!! | nlev                 | vertical_dimension                                    | number of vertical levels                                | count      |    0 | integer   |           | in     | F        |
!! | is_aerosol_aware     | flag_for_aerosol_physics                              | flag for aerosol-aware physics                           | flag       |    0 | logical   |           | in     | F        |
!! | nwfa2d               | tendency_of_water_friendly_aerosols_at_surface        | instantaneous fake water-friendly surface aerosol source | kg-1 s-1   |    1 | real      | kind_phys | inout  | T        |
!! | nifa2d               | tendency_of_ice_friendly_aerosols_at_surface          | instantaneous fake ice-friendly surface aerosol source   | kg-1 s-1   |    1 | real      | kind_phys | inout  | T        |
!! | nwfa                 | water_friendly_aerosol_number_concentration           | number concentration of water-friendly aerosols          | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | nifa                 | ice_friendly_aerosol_number_concentration             | number concentration of ice-friendly aerosols            | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | mpicomm              | mpi_comm                                              | MPI communicator                                         | index      |    0 | integer   |           | in     | F        |
!! | mpirank              | mpi_rank                                              | current MPI-rank                                         | index      |    0 | integer   |           | in     | F        |
!! | mpiroot              | mpi_root                                              | master MPI-rank                                          | index      |    0 | integer   |           | in     | F        |
!! | threads              | omp_threads                                           | number of OpenMP threads available to scheme             | count      |    0 | integer   |           | in     | F        |
!! | imp_physics          | flag_for_microphysics_scheme                          | choice of microphysics scheme                            | flag       |    0 | integer   |           | in     | F        |
!! | imp_physics_thompson | flag_for_thompson_microphysics_scheme                 | choice of Thompson microphysics scheme                   | flag       |    0 | integer   |           | in     | F        |
!! | errmsg               | ccpp_error_message                                    | error message for error handling in CCPP                 | none       |    0 | character | len=*     | out    | F        |
!! | errflg               | ccpp_error_flag                                       | error flag for error handling in CCPP                    | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine mp_thompson_hrrr_init(ncol, nlev, is_aerosol_aware, &
                                       nwfa2d, nifa2d, nwfa, nifa,   &
                                       mpicomm, mpirank, mpiroot,    &
                                       imp_physics,                  &
                                       imp_physics_thompson,         &
                                       threads, errmsg, errflg)

         implicit none

         ! Interface variables
         integer,                   intent(in)    :: ncol
         integer,                   intent(in)    :: nlev

         logical,                   intent(in)    :: is_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nwfa2d(1:ncol)
         real(kind_phys), optional, intent(inout) :: nifa2d(1:ncol)
         real(kind_phys), optional, intent(inout) :: nwfa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nifa(1:ncol,1:nlev)
         integer,                   intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         integer,                   intent(in)    :: threads
         integer,                   intent(in)    :: imp_physics
         integer,                   intent(in)    :: imp_physics_thompson
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables: dimensions used in thompson_init
         integer               :: ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return

         ! DH* temporary
         if (mpirank==mpiroot) then
            write(0,*) ' ----------------------------------------------------------------------------------------------------------------'
            write(0,*) ' --- WARNING --- the CCPP Thompson MP scheme is currently under development, use at your own risk --- WARNING ---'
            write(0,*) ' ----------------------------------------------------------------------------------------------------------------'
         end if
         ! *DH temporary

         if (imp_physics/=imp_physics_thompson) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Thompson MP"
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

         if (is_aerosol_aware .and. present(nwfa2d) &
                              .and. present(nifa2d) &
                              .and. present(nwfa)   &
                              .and. present(nifa)   ) then
            ! Call init
            call thompson_init(nwfa2d=nwfa2d, nifa2d=nifa2d, nwfa=nwfa, nifa=nifa,   &
                               ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                               its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, &
                               mpicomm=mpicomm, mpirank=mpirank, mpiroot=mpiroot,    &
                               threads=threads)
            if (errflg /= 0) return
         else if (is_aerosol_aware) then
            write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_hrrr_init:',                    &
                                       ' aerosol-aware microphysics require all of the following', &
                                       ' optional arguments: nifa2d, nwfa2d, nwfa, nifa'
            errflg = 1
            return
         else
            call thompson_init(ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                               its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, &
                               mpicomm=mpicomm, mpirank=mpirank, mpiroot=mpiroot,    &
                               threads=threads)
            if (errflg /= 0) return
         end if

         is_initialized = .true.

      end subroutine mp_thompson_hrrr_init

#if 0
!! \section arg_table_mp_thompson_hrrr_run Argument Table
!! | local_name      | standard_name                                                         | long_name                                                | units      | rank | type      |    kind   | intent | optional |
!! |-----------------|-----------------------------------------------------------------------|----------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ncol            | horizontal_loop_extent                                                | horizontal loop extent                                   | count      |    0 | integer   |           | in     | F        |
!! | nlev            | vertical_dimension                                                    | number of vertical levels                                | count      |    0 | integer   |           | in     | F        |
!! | con_g           | gravitational_acceleration                                            | gravitational acceleration                               | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | con_rd          | gas_constant_dry_air                                                  | ideal gas constant for dry air                           | J kg-1 K-1 |    0 | real      | kind_phys | in     | F        |
!! | spechum         | water_vapor_specific_humidity_updated_by_physics                      | water vapor specific humidity                            | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qc              | cloud_condensed_water_mixing_ratio_updated_by_physics                 | cloud water mixing ratio wrt dry+vapor (no condensates)  | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qr              | rain_water_mixing_ratio_updated_by_physics                            | rain water mixing ratio wrt dry+vapor (no condensates)   | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qi              | ice_water_mixing_ratio_updated_by_physics                             | ice water mixing ratio wrt dry+vapor (no condensates)    | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qs              | snow_water_mixing_ratio_updated_by_physics                            | snow water mixing ratio wrt dry+vapor (no condensates)   | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | qg              | graupel_mixing_ratio_updated_by_physics                               | graupel mixing ratio wrt dry+vapor (no condensates)      | kg kg-1    |    2 | real      | kind_phys | inout  | F        |
!! | ni              | ice_number_concentration_updated_by_physics                           | ice number concentration                                 | kg-1       |    2 | real      | kind_phys | inout  | F        |
!! | nr              | rain_number_concentration_updated_by_physics                          | rain number concentration                                | kg-1       |    2 | real      | kind_phys | inout  | F        |
!! | is_aerosol_aware| flag_for_aerosol_physics                                              | flag for aerosol-aware physics                           | flag       |    0 | logical   |           | in     | F        |
!! | nc              | cloud_droplet_number_concentration_updated_by_physics                 | cloud droplet number concentration                       | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | nwfa            | water_friendly_aerosol_number_concentration_updated_by_physics        | number concentration of water-friendly aerosols          | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | nifa            | ice_friendly_aerosol_number_concentration_updated_by_physics          | number concentration of ice-friendly aerosols            | kg-1       |    2 | real      | kind_phys | inout  | T        |
!! | nwfa2d          | tendency_of_water_friendly_aerosols_at_surface                        | instantaneous fake water-friendly surface aerosol source | kg-1 s-1   |    1 | real      | kind_phys | in     | T        |
!! | nifa2d          | tendency_of_ice_friendly_aerosols_at_surface                          | instantaneous fake ice-friendly surface aerosol source   | kg-1 s-1   |    1 | real      | kind_phys | in     | T        |
!! | tgrs            | air_temperature_updated_by_physics                                    | model layer mean temperature                             | K          |    2 | real      | kind_phys | inout  | F        |
!! | prsl            | air_pressure                                                          | mean layer pressure                                      | Pa         |    2 | real      | kind_phys | in     | F        |
!! | phii            | geopotential_at_interface                                             | geopotential at model layer interfaces                   | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | omega           | omega                                                                 | layer mean vertical velocity                             | Pa s-1     |    2 | real      | kind_phys | in     | F        |
!! | dtp             | time_step_for_physics                                                 | physics timestep                                         | s          |    0 | real      | kind_phys | in     | F        |
!! | rain            | lwe_thickness_of_explicit_precipitation_amount                        | explicit rainfall amount on physics timestep             | m          |    1 | real      | kind_phys | inout  | F        |
!! | graupel         | lwe_thickness_of_graupel_amount                                       | graupel fall on physics timestep                         | m          |    1 | real      | kind_phys | inout  | F        |
!! | ice             | lwe_thickness_of_ice_amount                                           | ice fall on physics timestep                             | m          |    1 | real      | kind_phys | inout  | F        |
!! | snow            | lwe_thickness_of_snow_amount                                          | snow fall on physics timestep                            | m          |    1 | real      | kind_phys | inout  | F        |
!! | sr              | ratio_of_snowfall_to_rainfall                                         | ratio of snowfall to large-scale rainfall                | frac       |    1 | real      | kind_phys | out    | F        |
!! | refl_10cm       | radar_reflectivity_10cm                                               | instantaneous refl_10cm                                  | dBZ        |    2 | real      | kind_phys | out    | F        |
!! | do_radar_ref    | flag_for_radar_reflectivity                                           | flag for radar reflectivity                              | flag       |    0 | logical   |           | in     | F        |
!! | re_cloud        | mean_effective_radius_for_liquid_cloud                                | mean effective radius for liquid cloud                   | micron     |    2 | real      | kind_phys | out    | T        |
!! | re_ice          | mean_effective_radius_for_ice_cloud                                   | mean effective radius for ice cloud                      | micron     |    2 | real      | kind_phys | out    | T        |
!! | re_snow         | mean_effective_radius_for_snow_flake                                  | mean effective radius for snow flake                     | micron     |    2 | real      | kind_phys | out    | T        |
!! | mpicomm         | mpi_comm                                                              | MPI communicator                                         | index      |    0 | integer   |           | in     | F        |
!! | mpirank         | mpi_rank                                                              | current MPI-rank                                         | index      |    0 | integer   |           | in     | F        |
!! | mpiroot         | mpi_root                                                              | master MPI-rank                                          | index      |    0 | integer   |           | in     | F        |
!! | errmsg          | ccpp_error_message                                                    | error message for error handling in CCPP                 | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                                       | error flag for error handling in CCPP                    | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine mp_thompson_hrrr_run(ncol, nlev, con_g, con_rd,         &
                              spechum, qc, qr, qi, qs, qg, ni, nr,       &
                              is_aerosol_aware, nc, nwfa, nifa,          &
                              nwfa2d, nifa2d,                            &
                              tgrs, prsl, phii, omega, dtp,              &
                              rain, graupel, ice, snow, sr,              &
                              refl_10cm, do_radar_ref,                   &
                              re_cloud, re_ice, re_snow,                 &
                              mpicomm, mpirank, mpiroot,                 &
                              errmsg, errflg)

         implicit none

         ! Interface variables

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
         real(kind_phys), optional, intent(in   ) :: nifa2d(1:ncol)
         ! State variables and timestep information
         real(kind_phys),           intent(inout) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: phii(1:ncol,1:nlev+1)
         real(kind_phys),           intent(in   ) :: omega(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: dtp
         ! Rain/snow/graupel fall amounts and fraction of frozen precip
         real(kind_phys),           intent(inout) :: rain(1:ncol)
         real(kind_phys),           intent(inout) :: graupel(1:ncol)
         real(kind_phys),           intent(inout) :: ice(1:ncol)
         real(kind_phys),           intent(inout) :: snow(1:ncol)
         real(kind_phys),           intent(  out) :: sr(1:ncol)
         ! Radar reflectivity
         real(kind_phys),           intent(  out) :: refl_10cm(1:ncol,1:nlev)
         logical,         optional, intent(in   ) :: do_radar_ref
         ! Cloud effective radii
         real(kind_phys), optional, intent(  out) :: re_cloud(1:ncol,1:nlev)
         real(kind_phys), optional, intent(  out) :: re_ice(1:ncol,1:nlev)
         real(kind_phys), optional, intent(  out) :: re_snow(1:ncol,1:nlev)
         ! MPI information
         integer,                   intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
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
         ! Vertical velocity and level width
         real(kind_phys) :: w(1:ncol,1:nlev)                ! m s-1
         real(kind_phys) :: dz(1:ncol,1:nlev)               ! m
         ! Rain/snow/graupel fall amounts
         real(kind_phys) :: rain_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: graupel_mp(1:ncol)              ! mm, dummy, not used
         real(kind_phys) :: ice_mp(1:ncol)                  ! mm, dummy, not used
         real(kind_phys) :: snow_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: delta_rain_mp(1:ncol)           ! mm
         real(kind_phys) :: delta_graupel_mp(1:ncol)        ! mm
         real(kind_phys) :: delta_ice_mp(1:ncol)            ! mm
         real(kind_phys) :: delta_snow_mp(1:ncol)           ! mm
         ! Radar reflectivity
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

         ! Check initialization state
         if (.not.is_initialized) then
            write(errmsg, fmt='((a))') 'mp_thompson_hrrr_run called before mp_thompson_hrrr_init'
            errflg = 1
            return
         end if

         ! Convert specific humidity/moist mixing ratios to dry mixing ratios
         qv_mp = spechum/(1.0_kind_phys-spechum)
         qc_mp = qc/(1.0_kind_phys-spechum)
         qr_mp = qr/(1.0_kind_phys-spechum)
         qi_mp = qi/(1.0_kind_phys-spechum)
         qs_mp = qs/(1.0_kind_phys-spechum)
         qg_mp = qg/(1.0_kind_phys-spechum)

         if (is_aerosol_aware .and. .not. (present(nc)     .and. &
                                           present(nwfa)   .and. &
                                           present(nifa)   .and. &
                                           present(nwfa2d) .and. &
                                           present(nifa2d)       )) then
            write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_hrrr_run:',  &
                                       ' aerosol-aware microphysics require all of the', &
                                       ' following optional arguments:', &
                                       ' nc, nwfa, nifa, nwfa2d, nifa2d'
            errflg = 1
            return
         end if

         ! Density of air in kg m-3
         rho = prsl/(con_rd*tgrs)

         ! Convert omega in Pa s-1 to vertical velocity w in m s-1
         w = -omega/(rho*con_g)

         ! Layer width in m from geopotential in m2 s-2
         dz = (phii(:,2:nlev+1) - phii(:,1:nlev)) / con_g

         ! Accumulated values inside Thompson scheme, not used;
         ! only use delta and add to inout variables (different units)
         rain_mp          = 0
         graupel_mp       = 0
         ice_mp           = 0
         snow_mp          = 0
         delta_rain_mp    = 0
         delta_graupel_mp = 0
         delta_ice_mp     = 0
         delta_snow_mp    = 0

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
         else if (.not.present(re_cloud) .and. .not.present(re_ice) .and. .not.present(re_snow)) then
             do_effective_radii = .false.
             has_reqc = 0
             has_reqi = 0
             has_reqs = 0
         else
             write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_hrrr_run:',  &
                                        ' all or none of the following optional', &
                                        ' arguments are required: re_cloud, re_ice, re_snow'
             errflg = 1
             return
         end if
         ! Initialize to zero, intent(out) variables
         re_cloud_mp = 0
         re_ice_mp   = 0
         re_snow_mp  = 0

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

#ifdef DEBUG_AEROSOLS
         if (mpirank==mpiroot) then
             write(0,*) "AEROSOL DEBUG: called mp_thompson_hrrr_run, is_aerosol_aware=",  is_aerosol_aware, &
                      & ", do_effective_radii=", do_effective_radii, ", do_radar_ref=", do_radar_ref, &
                      & ", diagflag=", diagflag, ", do_radar_ref_mp=", do_radar_ref_mp
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: prsl min/mean/max =", &
                              & minval(prsl), sum(prsl)/real(size(prsl)), maxval(prsl)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: tgrs min/mean/max =", &
                              & minval(tgrs), sum(tgrs)/real(size(tgrs)), maxval(tgrs)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: rho min/mean/max =", &
                              & minval(rho), sum(rho)/real(size(rho)), maxval(rho)
             if (is_aerosol_aware) then
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: nwfa2d min/mean/max =", &
                                    & minval(nwfa2d), sum(nwfa2d)/real(size(nwfa2d)), maxval(nwfa2d)
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: nifa2d min/mean/max =", &
                                    & minval(nifa2d), sum(nifa2d)/real(size(nifa2d)), maxval(nifa2d)
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: nwfa min/mean/max =", &
                                    & minval(nwfa), sum(nwfa)/real(size(nwfa)), maxval(nwfa)
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: nifa min/mean/max =", &
                                    & minval(nifa), sum(nifa)/real(size(nifa)), maxval(nifa)
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: nc min/mean/max =", &
                                    & minval(nc), sum(nc)/real(size(nc)), maxval(nc)
             end if
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: ni min/mean/max =", &
                                 & minval(ni), sum(ni)/real(size(ni)), maxval(ni)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: nr min/mean/max =", &
                                 & minval(nr), sum(nr)/real(size(nr)), maxval(nr)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: omega min/mean/max =", &
                                 & minval(omega), sum(omega)/real(size(omega)), maxval(omega)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: w min/mean/max =", &
                                 & minval(w), sum(w)/real(size(w)), maxval(w)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: re_cloud_mp min/mean/max =", &
                                 & minval(re_cloud_mp), sum(re_cloud_mp)/real(size(re_cloud_mp)), maxval(re_cloud_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: re_ice_mp min/mean/max =", &
                                 & minval(re_ice_mp), sum(re_ice_mp)/real(size(re_ice_mp)), maxval(re_ice_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: re_snow_mp min/mean/max =", &
                                 & minval(re_snow_mp), sum(re_snow_mp)/real(size(re_snow_mp)), maxval(re_snow_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: delta_rain_mp min/mean/max =", &
                                 & minval(delta_rain_mp), sum(delta_rain_mp)/real(size(delta_rain_mp)), maxval(delta_rain_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: delta_snow_mp min/mean/max =", &
                                 & minval(delta_snow_mp), sum(delta_snow_mp)/real(size(delta_snow_mp)), maxval(delta_snow_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: delta_ice_mp min/mean/max =", &
                                 & minval(delta_ice_mp), sum(delta_ice_mp)/real(size(delta_ice_mp)), maxval(delta_ice_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run before: delta_graupel_mp min/mean/max =", &
                                 & minval(delta_graupel_mp), sum(delta_graupel_mp)/real(size(delta_graupel_mp)), maxval(delta_graupel_mp)
         end if
#endif

         ! Call Thompson MP with or without aerosols
         if (is_aerosol_aware) then
            call mp_gt_driver(qv=qv_mp, qc=qc_mp, qr=qr_mp, qi=qi_mp, qs=qs_mp, qg=qg_mp,    &
                              ni=ni, nr=nr, nc=nc,                                           &
                              nwfa=nwfa, nifa=nifa, nwfa2d=nwfa2d, nifa2d=nifa2d,            &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtp,                        &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              icenc=ice_mp, icencv=delta_ice_mp,                             &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm,                                           &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              re_cloud=re_cloud_mp, re_ice=re_ice_mp, re_snow=re_snow_mp,    &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,          &
                              errmsg=errmsg, errflg=errflg)

         else
            call mp_gt_driver(qv=qv_mp, qc=qc_mp, qr=qr_mp, qi=qi_mp, qs=qs_mp, qg=qg_mp,    &
                              ni=ni, nr=nr, nc=nc,                                           &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtp,                        &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              icenc=ice_mp, icencv=delta_ice_mp,                             &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm,                                           &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              re_cloud=re_cloud_mp, re_ice=re_ice_mp, re_snow=re_snow_mp,    &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,          &
                              errmsg=errmsg, errflg=errflg)
         end if
         if (errflg/=0) return

         ! convert dry mixing ratios to specific humidity/moist mixing ratios
         spechum = qv_mp/(1.0_kind_phys+qv_mp)
         qc      = qc_mp/(1.0_kind_phys+qv_mp)
         qr      = qr_mp/(1.0_kind_phys+qv_mp)
         qi      = qi_mp/(1.0_kind_phys+qv_mp)
         qs      = qs_mp/(1.0_kind_phys+qv_mp)
         qg      = qg_mp/(1.0_kind_phys+qv_mp)

#ifdef DEBUG_AEROSOLS
         if (mpirank==mpiroot) then
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: prsl min/mean/max =", &
                              & minval(prsl), sum(prsl)/real(size(prsl)), maxval(prsl)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: tgrs min/mean/max =", &
                              & minval(tgrs), sum(tgrs)/real(size(tgrs)), maxval(tgrs)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: rho min/mean/max =", &
                              & minval(rho), sum(rho)/real(size(rho)), maxval(rho)
             if (is_aerosol_aware) then
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: nwfa2d min/mean/max =", &
                                    & minval(nwfa2d), sum(nwfa2d)/real(size(nwfa2d)), maxval(nwfa2d)
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: nifa2d min/mean/max =", &
                                    & minval(nifa2d), sum(nifa2d)/real(size(nifa2d)), maxval(nifa2d)
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: nwfa min/mean/max =", &
                                    & minval(nwfa), sum(nwfa)/real(size(nwfa)), maxval(nwfa)
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: nifa min/mean/max =", &
                                    & minval(nifa), sum(nifa)/real(size(nifa)), maxval(nifa)
                write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: nc min/mean/max =", &
                                    & minval(nc), sum(nc)/real(size(nc)), maxval(nc)
             end if
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: ni min/mean/max =", &
                                 & minval(ni), sum(ni)/real(size(ni)), maxval(ni)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: nr min/mean/max =", &
                                 & minval(nr), sum(nr)/real(size(nr)), maxval(nr)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: omega min/mean/max =", &
                                 & minval(omega), sum(omega)/real(size(omega)), maxval(omega)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: w min/mean/max =", &
                                 & minval(w), sum(w)/real(size(w)), maxval(w)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: re_cloud_mp min/mean/max =", &
                                 & minval(re_cloud_mp), sum(re_cloud_mp)/real(size(re_cloud_mp)), maxval(re_cloud_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: re_ice_mp min/mean/max =", &
                                 & minval(re_ice_mp), sum(re_ice_mp)/real(size(re_ice_mp)), maxval(re_ice_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: re_snow_mp min/mean/max =", &
                                 & minval(re_snow_mp), sum(re_snow_mp)/real(size(re_snow_mp)), maxval(re_snow_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: delta_rain_mp min/mean/max =", &
                                 & minval(delta_rain_mp), sum(delta_rain_mp)/real(size(delta_rain_mp)), maxval(delta_rain_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: delta_snow_mp min/mean/max =", &
                                 & minval(delta_snow_mp), sum(delta_snow_mp)/real(size(delta_snow_mp)), maxval(delta_snow_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: delta_ice_mp min/mean/max =", &
                                 & minval(delta_ice_mp), sum(delta_ice_mp)/real(size(delta_ice_mp)), maxval(delta_ice_mp)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp thompson run after: delta_graupel_mp min/mean/max =", &
                                 & minval(delta_graupel_mp), sum(delta_graupel_mp)/real(size(delta_graupel_mp)), maxval(delta_graupel_mp)
         end if
#endif

         ! Convert rainfall deltas from mm to m (on physics timestep); add to inout variables
         rain    = rain    + delta_rain_mp/1000.0_kind_phys
         graupel = graupel + delta_graupel_mp/1000.0_kind_phys
         ice     = ice     + delta_ice_mp/1000.0_kind_phys
         snow    = snow    + delta_snow_mp/1000.0_kind_phys

         if (do_effective_radii) then
            ! Convert m to micron
            re_cloud = re_cloud_mp*1.0E6_kind_phys
            re_ice   = re_ice_mp*1.0E6_kind_phys
            re_snow  = re_snow_mp*1.0E6_kind_phys
         end if

      end subroutine mp_thompson_hrrr_run

#if 0
!! \section arg_table_mp_thompson_hrrr_finalize Argument Table
!! | local_name      | standard_name                                                 | long_name                                              | units      | rank | type      |    kind   | intent | optional |
!! |-----------------|---------------------------------------------------------------|--------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | errmsg          | ccpp_error_message                                            | error message for error handling in CCPP               | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                               | error flag for error handling in CCPP                  | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine mp_thompson_hrrr_finalize(errmsg, errflg)

         implicit none

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (.not.is_initialized) return

         call thompson_finalize()

         is_initialized = .false.

      end subroutine mp_thompson_hrrr_finalize

end module mp_thompson_hrrr

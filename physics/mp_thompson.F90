!>\file mp_thompson.F90
!! This file contains aerosol-aware Thompson MP scheme.


!>\defgroup aathompson Aerosol-Aware Thompson MP Module
!! This module contains the aerosol-aware Thompson microphysics scheme.
module mp_thompson

      use machine, only : kind_phys

      use module_mp_thompson, only : thompson_init, mp_gt_driver, thompson_finalize, calc_effectRad
      use module_mp_thompson, only : naIN0, naIN1, naCCN0, naCCN1, eps, Nt_c
      use module_mp_thompson, only : re_qc_min, re_qc_max, re_qi_min, re_qi_max, re_qs_min, re_qs_max

      use module_mp_thompson_make_number_concentrations, only: make_IceNumber, make_DropletNumber, make_RainNumber

      implicit none

      public :: mp_thompson_init, mp_thompson_run, mp_thompson_finalize

      private

      logical :: is_initialized = .False.

   contains

!> This subroutine is a wrapper around the actual thompson_init().
!! \section arg_table_mp_thompson_init Argument Table
!! \htmlinclude mp_thompson_init.html
!!
      subroutine mp_thompson_init(ncol, nlev, con_g, con_rd, con_eps,   &
                                  restart, imp_physics,                 &
                                  imp_physics_thompson, convert_dry_rho,&
                                  spechum, qc, qr, qi, qs, qg, ni, nr,  &
                                  is_aerosol_aware, nc, nwfa2d, nifa2d, &
                                  nwfa, nifa, tgrs, prsl, phil, area,   &
                                  re_cloud, re_ice, re_snow,            &
                                  mpicomm, mpirank, mpiroot,            &
                                  threads, errmsg, errflg)

         implicit none

         ! Interface variables
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         real(kind_phys),           intent(in   ) :: con_g, con_rd, con_eps
         logical,                   intent(in   ) :: restart
         integer,                   intent(in   ) :: imp_physics
         integer,                   intent(in   ) :: imp_physics_thompson
         ! Hydrometeors
         logical,                   intent(in   ) :: convert_dry_rho
         real(kind_phys),           intent(inout) :: spechum(:,:)
         real(kind_phys),           intent(inout) :: qc(:,:)
         real(kind_phys),           intent(inout) :: qr(:,:)
         real(kind_phys),           intent(inout) :: qi(:,:)
         real(kind_phys),           intent(inout) :: qs(:,:)
         real(kind_phys),           intent(inout) :: qg(:,:)
         real(kind_phys),           intent(inout) :: ni(:,:)
         real(kind_phys),           intent(inout) :: nr(:,:)
         ! Aerosols
         logical,                   intent(in   ) :: is_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nc(:,:)
         real(kind_phys), optional, intent(inout) :: nwfa(:,:)
         real(kind_phys), optional, intent(inout) :: nifa(:,:)
         real(kind_phys), optional, intent(inout) :: nwfa2d(:)
         real(kind_phys), optional, intent(inout) :: nifa2d(:)
         ! State variables
         real(kind_phys),           intent(in   ) :: tgrs(:,:)
         real(kind_phys),           intent(in   ) :: prsl(:,:)
         real(kind_phys),           intent(in   ) :: phil(:,:)
         real(kind_phys),           intent(in   ) :: area(:)
         ! Cloud effective radii
         real(kind_phys), optional, intent(  out) :: re_cloud(:,:)
         real(kind_phys), optional, intent(  out) :: re_ice(:,:)
         real(kind_phys), optional, intent(  out) :: re_snow(:,:)
         ! MPI information
         integer,                   intent(in   ) :: mpicomm
         integer,                   intent(in   ) :: mpirank
         integer,                   intent(in   ) :: mpiroot
         ! Threading/blocking information
         integer,                   intent(in   ) :: threads
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         !
         real(kind_phys) :: qv(1:ncol,1:nlev)       ! kg kg-1 (water vapor mixing ratio)
         real(kind_phys) :: hgt(1:ncol,1:nlev)      ! m
         real(kind_phys) :: rho(1:ncol,1:nlev)      ! kg m-3
         real(kind_phys) :: orho(1:ncol,1:nlev)     ! m3 kg-1
         real(kind_phys) :: nc_local(1:ncol,1:nlev) ! needed because nc is only allocated if is_aerosol_aware is true
         !
         real (kind=kind_phys) :: h_01, airmass, niIN3, niCCN3
         integer :: i, k

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

         ! Consistency checks
         if (imp_physics/=imp_physics_thompson) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Thompson MP"
            errflg = 1
            return
         end if

         ! Call Thompson init
         call thompson_init(is_aerosol_aware_in=is_aerosol_aware, mpicomm=mpicomm, &
                            mpirank=mpirank, mpiroot=mpiroot, threads=threads,     &
                            errmsg=errmsg, errflg=errflg)
         if (errflg /= 0) return

         ! For restart runs, the init is done here
         if (restart) then
           is_initialized = .true.
           return
         end if

         ! Geopotential height in m2 s-2 to height in m
         hgt = phil/con_g

         ! Ensure non-negative mass mixing ratios of all water variables
         where(spechum<0) spechum = 1.0E-10     ! COMMENT, gthompsn, spechum should *never* be identically zero.
         where(qc<0)      qc = 0.0
         where(qr<0)      qr = 0.0
         where(qi<0)      qi = 0.0
         where(qs<0)      qs = 0.0
         where(qg<0)      qg = 0.0

         !> - Convert specific humidity to water vapor mixing ratio.
         !> - Also, hydrometeor variables are mass or number mixing ratio
         !> - either kg of species per kg of dry air, or per kg of (dry + vapor).

         qv = spechum/(1.0_kind_phys-spechum)

         if (convert_dry_rho) then
           qc = qc/(1.0_kind_phys-spechum)
           qr = qr/(1.0_kind_phys-spechum)
           qi = qi/(1.0_kind_phys-spechum)
           qs = qs/(1.0_kind_phys-spechum)
           qg = qg/(1.0_kind_phys-spechum)

           ni = ni/(1.0_kind_phys-spechum)
           nr = nr/(1.0_kind_phys-spechum)
           if (is_aerosol_aware) then
              nc = nc/(1.0_kind_phys-spechum)
              nwfa = nwfa/(1.0_kind_phys-spechum)
              nifa = nifa/(1.0_kind_phys-spechum)
           end if
         end if

         ! Density of moist air in kg m-3 and inverse density of air
         rho = con_eps*prsl/(con_rd*tgrs*(qv+con_eps))
         orho = 1.0/rho

         ! Ensure we have 1st guess ice number where mass non-zero but no number.
         where(qi .LE. 0.0) ni=0.0
         where(qi .GT. 0 .and. ni .LE. 0.0) ni = make_IceNumber(qi*rho, tgrs) * orho
         where(qi .EQ. 0.0 .and. ni .GT. 0.0) ni=0.0

         ! Ensure we have 1st guess rain number where mass non-zero but no number.
         where(qr .LE. 0.0) nr=0.0
         where(qr .GT. 0 .and. nr .LE. 0.0) nr = make_RainNumber(qr*rho, tgrs) * orho
         where(qr .EQ. 0.0 .and. nr .GT. 0.0) nr=0.0


         !..Check for existing aerosol data, both CCN and IN aerosols.  If missing
         !.. fill in just a basic vertical profile, somewhat boundary-layer following.
         if (is_aerosol_aware) then

           ! Potential cloud condensation nuclei (CCN)
           if (MAXVAL(nwfa) .lt. eps) then
             if (mpirank==mpiroot) write(*,*) ' Apparently there are no initial CCN aerosols.'
             do i = 1, ncol
               if (hgt(i,1).le.1000.0) then
                 h_01 = 0.8
               elseif (hgt(i,1).ge.2500.0) then
                 h_01 = 0.01
               else
                 h_01 = 0.8*cos(hgt(i,1)*0.001 - 1.0)
               endif
               niCCN3 = -1.0*ALOG(naCCN1/naCCN0)/h_01
               nwfa(i,1) = naCCN1+naCCN0*exp(-((hgt(i,2)-hgt(i,1))/1000.)*niCCN3)
               airmass = 1./orho(i,1) * (hgt(i,2)-hgt(i,1))*area(i) ! kg
               nwfa2d(i) = nwfa(i,1) * 0.000196 * (airmass*5.E-11)
               do k = 2, nlev
                 nwfa(i,k) = naCCN1+naCCN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niCCN3)
               enddo
             enddo
           else
             if (mpirank==mpiroot) write(*,*) ' Apparently initial CCN aerosols are present.'
             if (MAXVAL(nwfa2d) .lt. eps) then
               !+---+-----------------------------------------------------------------+
               !..Scale the lowest level aerosol data into an emissions rate.  This is
               !.. very far from ideal, but need higher emissions where larger amount
               !.. of (climo) existing and lesser emissions where there exists fewer to
               !.. begin as a first-order simplistic approach.  Later, proper connection to
               !.. emission inventory would be better, but, for now, scale like this:
               !.. where: Nwfa=50 per cc, emit 0.875E4 aerosols per second per grid box unit
               !..        that was tested as ~(20kmx20kmx50m = 2.E10 m**-3)
               !+---+-----------------------------------------------------------------+
               if (mpirank==mpiroot) write(*,*) ' Apparently there are no initial CCN aerosol surface emission rates.'
               do i = 1, ncol
                  airmass = 1./orho(i,1) * (hgt(i,2)-hgt(i,1))*area(i) ! kg
                  nwfa2d(i) = nwfa(i,1) * 0.000196 * (airmass*5.E-11)
               enddo
             else
                if (mpirank==mpiroot) write(*,*) ' Apparently initial CCN aerosol surface emission rates are present.'
             endif
           endif

           ! Potential ice nuclei (IN)
           if (MAXVAL(nifa) .lt. eps) then
             if (mpirank==mpiroot) write(*,*) ' Apparently there are no initial IN aerosols.'
             do i = 1, ncol
               if (hgt(i,1).le.1000.0) then
                  h_01 = 0.8
               elseif (hgt(i,1).ge.2500.0) then
                  h_01 = 0.01
               else
                  h_01 = 0.8*cos(hgt(i,1)*0.001 - 1.0)
               endif
               niIN3 = -1.0*ALOG(naIN1/naIN0)/h_01
               nifa(i,1) = naIN1+naIN0*exp(-((hgt(i,2)-hgt(i,1))/1000.)*niIN3)
               nifa2d(i) = 0.
               do k = 2, nlev
                  nifa(i,k) = naIN1+naIN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niIN3)
               enddo
             enddo
           else
             if (mpirank==mpiroot) write(*,*) ' Apparently initial IN aerosols are present.'
             if (MAXVAL(nifa2d) .lt. eps) then
               if (mpirank==mpiroot) write(*,*) ' Apparently there are no initial IN aerosol surface emission rates, set to zero.'
               ! calculate IN surface flux here, right now just set to zero
               nifa2d = 0.
             else
               if (mpirank==mpiroot) write(*,*) ' Apparently initial IN aerosol surface emission rates are present.'
             endif
           endif

           ! Ensure we have 1st guess cloud droplet number where mass non-zero but no number.
           where(qc .LE. 0.0) nc=0.0
           where(qc .GT. 0 .and. nc .LE. 0.0) nc = make_DropletNumber(qc*rho, nwfa*rho) * orho
           where(qc .EQ. 0.0 .and. nc .GT. 0.0) nc = 0.0

           ! Ensure non-negative aerosol number concentrations.
           where(nwfa .LE. 0.0) nwfa = 1.1E6
           where(nifa .LE. 0.0) nifa = naIN1*0.01

           ! Copy to local array for calculating cloud effective radii below
           nc_local = nc

         else

           ! Constant droplet concentration for single moment cloud water as in
           ! module_mp_thompson.F90, only needed for effective radii calculation
           nc_local = Nt_c/rho

         end if

         ! Calculate initial cloud effective radii if requested
         if (present(re_cloud) .and. present(re_ice) .and. present(re_snow)) then
           ! Effective radii [m] are now intent(out), bounds applied in calc_effectRad
           do i = 1, ncol
             call calc_effectRad (tgrs(i,:), prsl(i,:), qv(i,:), qc(i,:),          &
                                  nc_local(i,:), qi(i,:), ni(i,:), qs(i,:),        &
                                  re_cloud(i,:), re_ice(i,:), re_snow(i,:), 1, nlev)
             do k = 1, nlev
               re_cloud(i,k) = MAX(re_qc_min, MIN(re_cloud(i,k), re_qc_max))
               re_ice(i,k)   = MAX(re_qi_min, MIN(re_ice(i,k),   re_qi_max))
               re_snow(i,k)  = MAX(re_qs_min, MIN(re_snow(i,k),  re_qs_max))
             end do
           end do
           !! Convert to micron: required for bit-for-bit identical restarts;
           !! otherwise entering mp_thompson_init and converting mu to m and
           !! back (without updating re_*) introduces b4b differences.
           !! If this code is used, change units in metadata from m to um!
           !re_cloud = 1.0E6*re_cloud
           !re_ice   = 1.0E6*re_ice
           !re_snow  = 1.0E6*re_snow
         else if (present(re_cloud) .or. present(re_ice) .or. present(re_snow)) then
           write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_init:',  &
                                      ' all or none of the following optional', &
                                      ' arguments are required: re_cloud, re_ice, re_snow'
           errflg = 1
           return
         end if

         if (convert_dry_rho) then
           !qc = qc/(1.0_kind_phys+qv)
           !qr = qr/(1.0_kind_phys+qv)
           !qi = qi/(1.0_kind_phys+qv)
           !qs = qs/(1.0_kind_phys+qv)
           !qg = qg/(1.0_kind_phys+qv)

           ni = ni/(1.0_kind_phys+qv)
           nr = nr/(1.0_kind_phys+qv)
           if (is_aerosol_aware) then
              nc = nc/(1.0_kind_phys+qv)
              nwfa = nwfa/(1.0_kind_phys+qv)
              nifa = nifa/(1.0_kind_phys+qv)
           end if
         end if

         is_initialized = .true.

      end subroutine mp_thompson_init


!> \section arg_table_mp_thompson_run Argument Table
!! \htmlinclude mp_thompson_run.html
!!
!>\ingroup aathompson
!>\section gen_thompson_hrrr Thompson MP General Algorithm
!>@{
      subroutine mp_thompson_run(ncol, nlev, con_g, con_rd,        &
                              con_eps, convert_dry_rho,            &
                              spechum, qc, qr, qi, qs, qg, ni, nr, &
                              is_aerosol_aware, nc, nwfa, nifa,    &
                              nwfa2d, nifa2d,                      &
                              tgrs, prsl, phii, omega, dtp,        &
                              prcp, rain, graupel, ice, snow, sr,  &
                              refl_10cm, reset, do_radar_ref,      &
                              re_cloud, re_ice, re_snow,           &
                              mpicomm, mpirank, mpiroot,           &
                              errmsg, errflg)

         implicit none

         ! Interface variables

         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         real(kind_phys),           intent(in   ) :: con_g
         real(kind_phys),           intent(in   ) :: con_rd
         real(kind_phys),           intent(in   ) :: con_eps
         ! Hydrometeors
         logical,                   intent(in   ) :: convert_dry_rho
         real(kind_phys),           intent(inout) :: spechum(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qc(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qr(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qi(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qs(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qg(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: ni(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nr(1:ncol,1:nlev)
         ! Aerosols
         logical,                   intent(in)    :: is_aerosol_aware, reset
         ! The following arrays are not allocated if is_aerosol_aware is false
         real(kind_phys), optional, intent(inout) :: nc(:,:)
         real(kind_phys), optional, intent(inout) :: nwfa(:,:)
         real(kind_phys), optional, intent(inout) :: nifa(:,:)
         real(kind_phys), optional, intent(in   ) :: nwfa2d(:)
         real(kind_phys), optional, intent(in   ) :: nifa2d(:)
         ! State variables and timestep information
         real(kind_phys),           intent(inout) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: phii(1:ncol,1:nlev+1)
         real(kind_phys),           intent(in   ) :: omega(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: dtp
         ! Precip/rain/snow/graupel fall amounts and fraction of frozen precip
         real(kind_phys),           intent(  out) :: prcp(1:ncol)
         real(kind_phys),           intent(  out) :: rain(1:ncol)
         real(kind_phys),           intent(  out) :: graupel(1:ncol)
         real(kind_phys),           intent(  out) :: ice(1:ncol)
         real(kind_phys),           intent(  out) :: snow(1:ncol)
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
         real(kind_phys) :: rho(1:ncol,1:nlev)              !< kg m-3
         ! Water vapor mixing ratio (instead of specific humidity)
         real(kind_phys) :: qv(1:ncol,1:nlev)               !< kg kg-1
         ! Vertical velocity and level width
         real(kind_phys) :: w(1:ncol,1:nlev)                !< m s-1
         real(kind_phys) :: dz(1:ncol,1:nlev)               !< m
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
         integer         :: has_reqc
         integer         :: has_reqi
         integer         :: has_reqs
         ! DH* 2020-06-05 hardcode these values for not using random perturbations,
         ! hasn't been tested yet with this version of module_mp_thompson.F90
         integer, parameter :: rand_perturb_on = 0
         integer, parameter :: kme_stoch = 1
         !real(kind_phys) :: rand_pert(1:ncol,1:kme_stoch)
         ! *DH 2020-06-05
         ! Dimensions used in mp_gt_driver
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Check initialization state
         if (.not.is_initialized) then
            write(errmsg, fmt='((a))') 'mp_thompson_run called before mp_thompson_init'
            errflg = 1
            return
         end if

         if (is_aerosol_aware .and. .not. (present(nc)     .and. &
                                           present(nwfa)   .and. &
                                           present(nifa)   .and. &
                                           present(nwfa2d) .and. &
                                           present(nifa2d)       )) then
            write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_run:',  &
                                       ' aerosol-aware microphysics require all of the', &
                                       ' following optional arguments:', &
                                       ' nc, nwfa, nifa, nwfa2d, nifa2d'
            errflg = 1
            return
         end if

         !> - Convert specific humidity to water vapor mixing ratio.
         !> - Also, hydrometeor variables are mass or number mixing ratio
         !> - either kg of species per kg of dry air, or per kg of (dry + vapor).

         qv = spechum/(1.0_kind_phys-spechum)

         if (convert_dry_rho) then
           qc = qc/(1.0_kind_phys-spechum)
           qr = qr/(1.0_kind_phys-spechum)
           qi = qi/(1.0_kind_phys-spechum)
           qs = qs/(1.0_kind_phys-spechum)
           qg = qg/(1.0_kind_phys-spechum)

           ni = ni/(1.0_kind_phys-spechum)
           nr = nr/(1.0_kind_phys-spechum)
           if (is_aerosol_aware) then
              nc = nc/(1.0_kind_phys-spechum)
              nwfa = nwfa/(1.0_kind_phys-spechum)
              nifa = nifa/(1.0_kind_phys-spechum)
           end if
         end if

         !> - Density of air in kg m-3
         rho = con_eps*prsl/(con_rd*tgrs*(qv+con_eps))

         !> - Convert omega in Pa s-1 to vertical velocity w in m s-1
         w = -omega/(rho*con_g)

         !> - Layer width in m from geopotential in m2 s-2
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
             ! Initialize to zero, intent(out) variables
             re_cloud = 0
             re_ice   = 0
             re_snow  = 0
         else if (.not.present(re_cloud) .and. .not.present(re_ice) .and. .not.present(re_snow)) then
             do_effective_radii = .false.
             has_reqc = 0
             has_reqi = 0
             has_reqs = 0
         else
             write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_run:',  &
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

         !> - Call mp_gt_driver() with or without aerosols
         if (is_aerosol_aware) then
            call mp_gt_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr,        &
                              nc=nc, nwfa=nwfa, nifa=nifa, nwfa2d=nwfa2d, nifa2d=nifa2d,     &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtp,                        &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              icenc=ice_mp, icencv=delta_ice_mp,                             &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm,                                           &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow,             &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              rand_perturb_on=rand_perturb_on, kme_stoch=kme_stoch,          &
                              ! DH* 2020-06-05 not passing this optional argument, see
                              !       comment in module_mp_thompson.F90 / mp_gt_driver
                              !rand_pert=rand_pert,                                          &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,          &
                              errmsg=errmsg, errflg=errflg, reset=reset)

         else
            call mp_gt_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr,        &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtp,                        &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              icenc=ice_mp, icencv=delta_ice_mp,                             &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm,                                           &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow,             &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              rand_perturb_on=rand_perturb_on, kme_stoch=kme_stoch,          &
                              ! DH* 2020-06-05 not passing this optional argument, see
                              !       comment in module_mp_thompson.F90 / mp_gt_driver
                              !rand_pert=rand_pert,                                          &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,          &
                              errmsg=errmsg, errflg=errflg, reset=reset)
         end if
         if (errflg/=0) return

         !> - Convert water vapor mixing ratio back to specific humidity
         spechum = qv/(1.0_kind_phys+qv)

         if (convert_dry_rho) then
           qc = qc/(1.0_kind_phys+qv)
           qr = qr/(1.0_kind_phys+qv)
           qi = qi/(1.0_kind_phys+qv)
           qs = qs/(1.0_kind_phys+qv)
           qg = qg/(1.0_kind_phys+qv)

           ni = ni/(1.0_kind_phys+qv)
           nr = nr/(1.0_kind_phys+qv)
           if (is_aerosol_aware) then
              nc = nc/(1.0_kind_phys+qv)
              nwfa = nwfa/(1.0_kind_phys+qv)
              nifa = nifa/(1.0_kind_phys+qv)
           end if
         end if

         !> - Convert rainfall deltas from mm to m (on physics timestep); add to inout variables
         ! "rain" in Thompson MP refers to precipitation (total of liquid rainfall+snow+graupel+ice)
         prcp    = max(0.0, delta_rain_mp/1000.0_kind_phys)
         graupel = max(0.0, delta_graupel_mp/1000.0_kind_phys)
         ice     = max(0.0, delta_ice_mp/1000.0_kind_phys)
         snow    = max(0.0, delta_snow_mp/1000.0_kind_phys)
         rain    = max(0.0, (delta_rain_mp - (delta_graupel_mp + delta_ice_mp + delta_snow_mp))/1000.0_kind_phys)

      end subroutine mp_thompson_run
!>@}

!! \section arg_table_mp_thompson_finalize Argument Table
!! \htmlinclude mp_thompson_finalize.html
!!
      subroutine mp_thompson_finalize(errmsg, errflg)

         implicit none

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (.not.is_initialized) return

         call thompson_finalize()

         is_initialized = .false.

      end subroutine mp_thompson_finalize

end module mp_thompson

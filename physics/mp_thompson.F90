!>\file mp_thompson.F90
!! This file contains aerosol-aware Thompson MP scheme.


!>\defgroup aathompson Aerosol-Aware Thompson MP Module
!! This module contains the aerosol-aware Thompson microphysics scheme.
module mp_thompson

      use machine, only : kind_phys

      use module_mp_thompson, only : thompson_init, mp_gt_driver, thompson_finalize, calc_effectRad
      use module_mp_thompson, only : naIN0, naIN1, naCCN0, naCCN1, eps, Nt_c_l, Nt_c_o
      use module_mp_thompson, only : re_qc_min, re_qc_max, re_qi_min, re_qi_max, re_qs_min, re_qs_max

      use module_mp_thompson_make_number_concentrations, only: make_IceNumber, make_DropletNumber, make_RainNumber

      implicit none

      public :: mp_thompson_init, mp_thompson_run, mp_thompson_finalize

      private

      logical :: is_initialized = .False.

      integer, parameter :: ext_ndiag3d = 37

   contains

!> This subroutine is a wrapper around the actual thompson_init().
!! \section arg_table_mp_thompson_init Argument Table
!! \htmlinclude mp_thompson_init.html
!!
      subroutine mp_thompson_init(ncol, nlev, con_g, con_rd, con_eps,      &
                                  restart, imp_physics,                    &
                                  imp_physics_thompson, convert_dry_rho,   &
                                  spechum, qc, qr, qi, qs, qg, ni, nr,     &
                                  is_aerosol_aware,  merra2_aerosol_aware, &
                                  nc, nwfa2d, nifa2d,                      &
                                  nwfa, nifa, tgrs, prsl, phil, area,      &
                                  aerfld, mpicomm, mpirank, mpiroot,       &
                                  threads, ext_diag, diag3d,               &
                                  errmsg, errflg)

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
         logical,                   intent(in   ) :: merra2_aerosol_aware
         real(kind_phys),           intent(inout) :: nc(:,:)
         real(kind_phys),           intent(inout) :: nwfa(:,:)
         real(kind_phys),           intent(inout) :: nifa(:,:)
         real(kind_phys),           intent(inout) :: nwfa2d(:)
         real(kind_phys),           intent(inout) :: nifa2d(:)
         real(kind_phys),           intent(in)    :: aerfld(:,:,:)
         ! State variables
         real(kind_phys),           intent(in   ) :: tgrs(:,:)
         real(kind_phys),           intent(in   ) :: prsl(:,:)
         real(kind_phys),           intent(in   ) :: phil(:,:)
         real(kind_phys),           intent(in   ) :: area(:)
         ! MPI information
         integer,                   intent(in   ) :: mpicomm
         integer,                   intent(in   ) :: mpirank
         integer,                   intent(in   ) :: mpiroot
         ! Threading/blocking information
         integer,                   intent(in   ) :: threads
         ! Extended diagnostics
         logical,                   intent(in   ) :: ext_diag
         real(kind_phys),           intent(in   ) :: diag3d(:,:,:)
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
         real (kind=kind_phys) :: h_01, z1, niIN3, niCCN3
         integer :: i, k

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return

         ! Consistency checks
         if (imp_physics/=imp_physics_thompson) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Thompson MP"
            errflg = 1
            return
         end if

         if (ext_diag) then
            if (size(diag3d,dim=3) /= ext_ndiag3d) then
               write(errmsg,'(*(a))') "Logic error: number of diagnostic 3d arrays from model does not match requirements"
               errflg = 1
               return
            end if
         end if

         if (is_aerosol_aware .and. merra2_aerosol_aware) then
            write(errmsg,'(*(a))') "Logic error: Only one Thompson aerosol option can be true, either is_aerosol_aware or merra2_aerosol_aware)"
            errflg = 1
            return
         end if

         ! Call Thompson init
         call thompson_init(is_aerosol_aware_in=is_aerosol_aware,              &
                            merra2_aerosol_aware_in=merra2_aerosol_aware,      &
                            mpicomm=mpicomm, mpirank=mpirank, mpiroot=mpiroot, &
                            threads=threads, errmsg=errmsg, errflg=errflg)
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
               z1 = hgt(i,2)-hgt(i,1)
               nwfa2d(i) = nwfa(i,1) * 0.000196 * (50./z1)
               do k = 2, nlev
                 nwfa(i,k) = naCCN1+naCCN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niCCN3)
               enddo
             enddo
           else if (merra2_aerosol_aware) then
             call get_niwfa(aerfld, nifa, nwfa, ncol, nlev)
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
                  z1 = hgt(i,2)-hgt(i,1)
                  nwfa2d(i) = nwfa(i,1) * 0.000196 * (50./z1)
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
 
        else if (merra2_aerosol_aware) then

           ! Ensure we have 1st guess cloud droplet number where mass non-zero but no number.
           where(qc .LE. 0.0) nc=0.0
           where(qc .GT. 0 .and. nc .LE. 0.0) nc = make_DropletNumber(qc*rho, nwfa*rho) * orho
           where(qc .EQ. 0.0 .and. nc .GT. 0.0) nc = 0.0

         else

           ! Constant droplet concentration for single moment cloud water as in
           ! module_mp_thompson.F90, only needed for effective radii calculation
           nc_local = Nt_c_l/rho

         end if

         if (convert_dry_rho) then
           !qc = qc/(1.0_kind_phys+qv)
           !qr = qr/(1.0_kind_phys+qv)
           !qi = qi/(1.0_kind_phys+qv)
           !qs = qs/(1.0_kind_phys+qv)
           !qg = qg/(1.0_kind_phys+qv)

           ni = ni/(1.0_kind_phys+qv)
           nr = nr/(1.0_kind_phys+qv)
           if (is_aerosol_aware .or. merra2_aerosol_aware) then
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
                              is_aerosol_aware,                    &
                              merra2_aerosol_aware, nc, nwfa, nifa,&
                              nwfa2d, nifa2d, aero_ind_fdb,        &
                              tgrs, prsl, phii, omega,             &
                              sedi_semi, decfl, islmsk, dtp,       &
                              dt_inner,                            &
                              first_time_step, istep, nsteps,      &
                              prcp, rain, graupel, ice, snow, sr,  &
                              refl_10cm, fullradar_diag,           &
                              do_radar_ref, aerfld,                &
                              mpicomm, mpirank, mpiroot, blkno,    &
                              ext_diag, diag3d, reset_diag3d,      &
                              spp_wts_mp, spp_mp, n_var_spp,       &
                              spp_prt_list, spp_var_list,          &
                              spp_stddev_cutoff,                   &
                              cplchm, pfi_lsan, pfl_lsan,          &
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
         real(kind_phys),           intent(inout) :: spechum(:,:)
         real(kind_phys),           intent(inout) :: qc(:,:)
         real(kind_phys),           intent(inout) :: qr(:,:)
         real(kind_phys),           intent(inout) :: qi(:,:)
         real(kind_phys),           intent(inout) :: qs(:,:)
         real(kind_phys),           intent(inout) :: qg(:,:)
         real(kind_phys),           intent(inout) :: ni(:,:)
         real(kind_phys),           intent(inout) :: nr(:,:)
         ! Aerosols
         logical,                   intent(in)    :: is_aerosol_aware, fullradar_diag 
         logical,                   intent(in)    :: merra2_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nc(:,:)
         real(kind_phys), optional, intent(inout) :: nwfa(:,:)
         real(kind_phys), optional, intent(inout) :: nifa(:,:)
         real(kind_phys), optional, intent(in   ) :: nwfa2d(:)
         real(kind_phys), optional, intent(in   ) :: nifa2d(:)
         real(kind_phys),           intent(in)    :: aerfld(:,:,:)
         logical,         optional, intent(in   ) :: aero_ind_fdb
         ! State variables and timestep information
         real(kind_phys),           intent(inout) :: tgrs(:,:)
         real(kind_phys),           intent(in   ) :: prsl(:,:)
         real(kind_phys),           intent(in   ) :: phii(:,:)
         real(kind_phys),           intent(in   ) :: omega(:,:)
         integer,                   intent(in   ) :: islmsk(:)
         real(kind_phys),           intent(in   ) :: dtp
         logical,                   intent(in   ) :: first_time_step
         integer,                   intent(in   ) :: istep, nsteps
         real,                      intent(in   ) :: dt_inner
         ! Precip/rain/snow/graupel fall amounts and fraction of frozen precip
         real(kind_phys),           intent(inout) :: prcp(:)
         real(kind_phys),           intent(inout) :: rain(:)
         real(kind_phys),           intent(inout) :: graupel(:)
         real(kind_phys),           intent(inout) :: ice(:)
         real(kind_phys),           intent(inout) :: snow(:)
         real(kind_phys),           intent(  out) :: sr(:)
         ! Radar reflectivity
         real(kind_phys),           intent(inout) :: refl_10cm(:,:)
         logical,                   intent(in   ) :: do_radar_ref
         logical,                   intent(in)    :: sedi_semi
         integer,                   intent(in)    :: decfl
         ! MPI and block information
         integer,                   intent(in)    :: blkno
         integer,                   intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         ! Extended diagnostic output
         logical,                   intent(in)    :: ext_diag
         real(kind_phys), target,   intent(inout) :: diag3d(:,:,:)
         logical,                   intent(in)    :: reset_diag3d

         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg
         
         ! SPP
         integer,                   intent(in) :: spp_mp
         integer,                   intent(in) :: n_var_spp
         real(kind_phys),           intent(in) :: spp_wts_mp(:,:)
         real(kind_phys),           intent(in) :: spp_prt_list(:)
         character(len=3),          intent(in) :: spp_var_list(:)
         real(kind_phys),           intent(in) :: spp_stddev_cutoff(:)

         logical, intent (in) :: cplchm
         ! ice and liquid water 3d precipitation fluxes - only allocated if cplchm is .true.
         real(kind=kind_phys), intent(inout), dimension(:,:) :: pfi_lsan
         real(kind=kind_phys), intent(inout), dimension(:,:) :: pfl_lsan

         ! Local variables

         ! Reduced time step if subcycling is used
         real(kind_phys) :: dtstep
         integer         :: ndt
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

         real(kind_phys) :: pfils(1:ncol,1:nlev,1)
         real(kind_phys) :: pflls(1:ncol,1:nlev,1)
         ! Radar reflectivity
         logical         :: diagflag                        ! must be true if do_radar_ref is true, not used otherwise
         integer         :: do_radar_ref_mp                 ! integer instead of logical do_radar_ref
         ! Effective cloud radii - turned off in CCPP (taken care off in radiation)
         logical, parameter :: do_effective_radii = .false.
         integer, parameter :: has_reqc = 0
         integer, parameter :: has_reqi = 0
         integer, parameter :: has_reqs = 0
         integer, parameter :: kme_stoch = 1
         integer         :: spp_mp_opt 
         ! Dimensions used in mp_gt_driver
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte
         ! Pointer arrays for extended diagnostics
         !real(kind_phys), dimension(:,:,:), pointer :: vts1       => null()
         !real(kind_phys), dimension(:,:,:), pointer :: txri       => null()
         !real(kind_phys), dimension(:,:,:), pointer :: txrc       => null()
         real(kind_phys), dimension(:,:,:), pointer :: prw_vcdc   => null()
         real(kind_phys), dimension(:,:,:), pointer :: prw_vcde   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tpri_inu   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tpri_ide_d => null()
         real(kind_phys), dimension(:,:,:), pointer :: tpri_ide_s => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprs_ide   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprs_sde_d => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprs_sde_s => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprg_gde_d => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprg_gde_s => null()
         real(kind_phys), dimension(:,:,:), pointer :: tpri_iha   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tpri_wfz   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tpri_rfz   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprg_rfz   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprs_scw   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprg_scw   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprg_rcs   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprs_rcs   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprr_rci   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprg_rcg   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprw_vcd_c => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprw_vcd_e => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprr_sml   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprr_gml   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprr_rcg   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprr_rcs   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tprv_rev   => null()
         real(kind_phys), dimension(:,:,:), pointer :: tten3      => null()
         real(kind_phys), dimension(:,:,:), pointer :: qvten3     => null()
         real(kind_phys), dimension(:,:,:), pointer :: qrten3     => null()
         real(kind_phys), dimension(:,:,:), pointer :: qsten3     => null()
         real(kind_phys), dimension(:,:,:), pointer :: qgten3     => null()
         real(kind_phys), dimension(:,:,:), pointer :: qiten3     => null()
         real(kind_phys), dimension(:,:,:), pointer :: niten3     => null()
         real(kind_phys), dimension(:,:,:), pointer :: nrten3     => null()
         real(kind_phys), dimension(:,:,:), pointer :: ncten3     => null()
         real(kind_phys), dimension(:,:,:), pointer :: qcten3     => null()

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (first_time_step .and. istep==1 .and. blkno==1) then
            ! Check initialization state
            if (.not.is_initialized) then
               write(errmsg, fmt='((a))') 'mp_thompson_run called before mp_thompson_init'
               errflg = 1
               return
            end if
            ! Check forr optional arguments of aerosol-aware microphysics
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
            else if (merra2_aerosol_aware .and. .not. (present(nc)     .and. &
                                                       present(nwfa)   .and. &
                                                       present(nifa)         )) then
              write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_run:', &
                                         ' merra2 aerosol-aware microphysics require the', &
                                         ' following optional arguments: nc, nwfa, nifa'
              errflg = 1
              return
            end if
            ! Consistency cheecks - subcycling and inner loop at the same time are not supported
            if (nsteps>1 .and. dt_inner < dtp) then
               write(errmsg,'(*(a))') "Logic error: Subcycling and inner loop cannot be used at the same time"
               errflg = 1
               return
            else if (mpirank==mpiroot .and. nsteps>1) then
               write(*,'(a,i0,a,a,f6.2,a)') 'Thompson MP is using ', nsteps, ' substep(s) per time step with an ', &
                                            'effective time step of ', dtp/real(nsteps, kind=kind_phys), ' seconds'
            else if (mpirank==mpiroot .and. dt_inner < dtp) then
               ndt = max(nint(dtp/dt_inner),1)
               write(*,'(a,i0,a,a,f6.2,a)') 'Thompson MP is using ', ndt, ' inner loops per time step with an ', &
                                            'effective time step of ', dtp/real(ndt, kind=kind_phys), ' seconds'
            end if
         end if

         ! Set stochastic physics selection to apply all perturbations
         if ( spp_mp==7 ) then
            spp_mp_opt=7
         else
            spp_mp_opt=0
         endif

         ! Set reduced time step if subcycling is used
         if (nsteps>1) then
            dtstep = dtp/real(nsteps, kind=kind_phys)
         else
            dtstep = dtp
         end if

         !> - Convert specific humidity to water vapor mixing ratio.
         !> - Also, hydrometeor variables are mass or number mixing ratio
         !> - either kg of species per kg of dry air, or per kg of (dry + vapor).

         ! DH* - do this only if istep == 1? Would be ok if it was
         ! guaranteed that nothing else in the same subcycle group
         ! was using these arrays, but it is somewhat dangerous.
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
         ! *DH

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
         if(cplchm) then
           pfi_lsan = 0.0
           pfl_lsan = 0.0
         end if

         ! Set pointers for extended diagnostics
         set_extended_diagnostic_pointers: if (ext_diag) then
            if (reset_diag3d) then
               diag3d = 0.0
            end if
            !vts1       => diag3d(:,:,X:X)
            !txri       => diag3d(:,:,X:X)
            !txrc       => diag3d(:,:,X:X)
            prw_vcdc   => diag3d(:,:,1:1)
            prw_vcde   => diag3d(:,:,2:2)
            tpri_inu   => diag3d(:,:,3:3)
            tpri_ide_d => diag3d(:,:,4:4)
            tpri_ide_s => diag3d(:,:,5:5)
            tprs_ide   => diag3d(:,:,6:6)
            tprs_sde_d => diag3d(:,:,7:7)
            tprs_sde_s => diag3d(:,:,8:8)
            tprg_gde_d => diag3d(:,:,9:9)
            tprg_gde_s => diag3d(:,:,10:10)
            tpri_iha   => diag3d(:,:,11:11)
            tpri_wfz   => diag3d(:,:,12:12)
            tpri_rfz   => diag3d(:,:,13:13)
            tprg_rfz   => diag3d(:,:,14:14)
            tprs_scw   => diag3d(:,:,15:15)
            tprg_scw   => diag3d(:,:,16:16)
            tprg_rcs   => diag3d(:,:,17:17)
            tprs_rcs   => diag3d(:,:,18:18)
            tprr_rci   => diag3d(:,:,19:19)
            tprg_rcg   => diag3d(:,:,20:20)
            tprw_vcd_c => diag3d(:,:,21:21)
            tprw_vcd_e => diag3d(:,:,22:22)
            tprr_sml   => diag3d(:,:,23:23)
            tprr_gml   => diag3d(:,:,24:24)
            tprr_rcg   => diag3d(:,:,25:25)
            tprr_rcs   => diag3d(:,:,26:26)
            tprv_rev   => diag3d(:,:,27:27)
            tten3      => diag3d(:,:,28:28)
            qvten3     => diag3d(:,:,29:29)
            qrten3     => diag3d(:,:,30:30)
            qsten3     => diag3d(:,:,31:31)
            qgten3     => diag3d(:,:,32:32)
            qiten3     => diag3d(:,:,33:33)
            niten3     => diag3d(:,:,34:34)
            nrten3     => diag3d(:,:,35:35)
            ncten3     => diag3d(:,:,36:36)
            qcten3     => diag3d(:,:,37:37)
         end if set_extended_diagnostic_pointers
         if (merra2_aerosol_aware) then
           call get_niwfa(aerfld, nifa, nwfa, ncol, nlev)
         end if
         !> - Call mp_gt_driver() with or without aerosols, with or without effective radii, ...
         if (is_aerosol_aware .or. merra2_aerosol_aware) then
            call mp_gt_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr,        &
                              nc=nc, nwfa=nwfa, nifa=nifa, nwfa2d=nwfa2d, nifa2d=nifa2d,     &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtstep, dt_inner=dt_inner,  &
                              sedi_semi=sedi_semi, decfl=decfl, lsm=islmsk,                  &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              icenc=ice_mp, icencv=delta_ice_mp,                             &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm,                                           &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              aero_ind_fdb=aero_ind_fdb, rand_perturb_on=spp_mp_opt,         &
                              kme_stoch=kme_stoch,                                           &
                              rand_pert=spp_wts_mp, spp_var_list=spp_var_list,               &
                              spp_prt_list=spp_prt_list, n_var_spp=n_var_spp,                &
                              spp_stddev_cutoff=spp_stddev_cutoff,                           &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,          &
                              fullradar_diag=fullradar_diag, istep=istep, nsteps=nsteps,     &
                              first_time_step=first_time_step, errmsg=errmsg, errflg=errflg, &
                              ! Extended diagnostics
                              ext_diag=ext_diag,                                             &
                              ! vts1=vts1, txri=txri, txrc=txrc,                             &
                              prw_vcdc=prw_vcdc,                                             &
                              prw_vcde=prw_vcde, tpri_inu=tpri_inu, tpri_ide_d=tpri_ide_d,   &
                              tpri_ide_s=tpri_ide_s, tprs_ide=tprs_ide,                      &
                              tprs_sde_d=tprs_sde_d,                                         &
                              tprs_sde_s=tprs_sde_s, tprg_gde_d=tprg_gde_d,                  &
                              tprg_gde_s=tprg_gde_s, tpri_iha=tpri_iha,                      &
                              tpri_wfz=tpri_wfz, tpri_rfz=tpri_rfz, tprg_rfz=tprg_rfz,       &
                              tprs_scw=tprs_scw, tprg_scw=tprg_scw, tprg_rcs=tprg_rcs,       &
                              tprs_rcs=tprs_rcs,                                             &
                              tprr_rci=tprr_rci, tprg_rcg=tprg_rcg, tprw_vcd_c=tprw_vcd_c,   &
                              tprw_vcd_e=tprw_vcd_e, tprr_sml=tprr_sml, tprr_gml=tprr_gml,   &
                              tprr_rcg=tprr_rcg, tprr_rcs=tprr_rcs,                          &
                              tprv_rev=tprv_rev, tten3=tten3,                                &
                              qvten3=qvten3, qrten3=qrten3, qsten3=qsten3, qgten3=qgten3,    &
                              qiten3=qiten3, niten3=niten3, nrten3=nrten3, ncten3=ncten3,    &
                              qcten3=qcten3, pfils=pfils, pflls=pflls)
         else
            call mp_gt_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr,        &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtstep, dt_inner=dt_inner,  &
                              sedi_semi=sedi_semi, decfl=decfl, lsm=islmsk,                  &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              icenc=ice_mp, icencv=delta_ice_mp,                             &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm,                                           &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              rand_perturb_on=spp_mp_opt, kme_stoch=kme_stoch,               &
                              rand_pert=spp_wts_mp, spp_var_list=spp_var_list,               &
                              spp_prt_list=spp_prt_list, n_var_spp=n_var_spp,                &
                              spp_stddev_cutoff=spp_stddev_cutoff,                           &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,          &
                              fullradar_diag=fullradar_diag, istep=istep, nsteps=nsteps,     &
                              first_time_step=first_time_step, errmsg=errmsg, errflg=errflg, &
                              ! Extended diagnostics
                              ext_diag=ext_diag,                                             &
                              ! vts1=vts1, txri=txri, txrc=txrc,                             &
                              prw_vcdc=prw_vcdc,                                             &
                              prw_vcde=prw_vcde, tpri_inu=tpri_inu, tpri_ide_d=tpri_ide_d,   &
                              tpri_ide_s=tpri_ide_s, tprs_ide=tprs_ide,                      &
                              tprs_sde_d=tprs_sde_d,                                         &
                              tprs_sde_s=tprs_sde_s, tprg_gde_d=tprg_gde_d,                  &
                              tprg_gde_s=tprg_gde_s, tpri_iha=tpri_iha,                      &
                              tpri_wfz=tpri_wfz, tpri_rfz=tpri_rfz, tprg_rfz=tprg_rfz,       &
                              tprs_scw=tprs_scw, tprg_scw=tprg_scw, tprg_rcs=tprg_rcs,       &
                              tprs_rcs=tprs_rcs,                                             &
                              tprr_rci=tprr_rci, tprg_rcg=tprg_rcg, tprw_vcd_c=tprw_vcd_c,   &
                              tprw_vcd_e=tprw_vcd_e, tprr_sml=tprr_sml, tprr_gml=tprr_gml,   &
                              tprr_rcg=tprr_rcg, tprr_rcs=tprr_rcs,                          &
                              tprv_rev=tprv_rev, tten3=tten3,                                &
                              qvten3=qvten3, qrten3=qrten3, qsten3=qsten3, qgten3=qgten3,    &
                              qiten3=qiten3, niten3=niten3, nrten3=nrten3, ncten3=ncten3,    &
                              qcten3=qcten3, pfils=pfils, pflls=pflls)
         end if
         if (errflg/=0) return

         ! DH* - do this only if istep == nsteps? Would be ok if it was
         ! guaranteed that nothing else in the same subcycle group
         ! was using these arrays, but it is somewhat dangerous.

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
           if (is_aerosol_aware .or. merra2_aerosol_aware) then
              nc = nc/(1.0_kind_phys+qv)
              nwfa = nwfa/(1.0_kind_phys+qv)
              nifa = nifa/(1.0_kind_phys+qv)
           end if
         end if
         ! *DH

         !> - Convert rainfall deltas from mm to m (on physics timestep); add to inout variables
         ! "rain" in Thompson MP refers to precipitation (total of liquid rainfall+snow+graupel+ice)
         prcp    = prcp    + max(0.0, delta_rain_mp/1000.0_kind_phys)
         graupel = graupel + max(0.0, delta_graupel_mp/1000.0_kind_phys)
         ice     = ice     + max(0.0, delta_ice_mp/1000.0_kind_phys)
         snow    = snow    + max(0.0, delta_snow_mp/1000.0_kind_phys)
         rain    = rain    + max(0.0, (delta_rain_mp - (delta_graupel_mp + delta_ice_mp + delta_snow_mp))/1000.0_kind_phys)

         ! Recompute sr at last subcycling step
         if (nsteps>1 .and. istep == nsteps) then
           ! Unlike inside mp_gt_driver, rain does not contain frozen precip
           sr = (snow + graupel + ice)/(rain + snow + graupel + ice +1.e-12)
         end if

         ! output instantaneous ice/snow and rain water 3d precipitation fluxes
         if(cplchm) then
           pfi_lsan(:,:) = pfils(:,:,1)
           pfl_lsan(:,:) = pflls(:,:,1)
         end if

         unset_extended_diagnostic_pointers: if (ext_diag) then
           !vts1       => null()
           !txri       => null()
           !txrc       => null()
           prw_vcdc   => null()
           prw_vcde   => null()
           tpri_inu   => null()
           tpri_ide_d => null()
           tpri_ide_s => null()
           tprs_ide   => null()
           tprs_sde_d => null()
           tprs_sde_s => null()
           tprg_gde_d => null()
           tprg_gde_s => null()
           tpri_iha   => null()
           tpri_wfz   => null()
           tpri_rfz   => null()
           tprg_rfz   => null()
           tprs_scw   => null()
           tprg_scw   => null()
           tprg_rcs   => null()
           tprs_rcs   => null()
           tprr_rci   => null()
           tprg_rcg   => null()
           tprw_vcd_c => null()
           tprw_vcd_e => null()
           tprr_sml   => null()
           tprr_gml   => null()
           tprr_rcg   => null()
           tprr_rcs   => null()
           tprv_rev   => null()
           tten3      => null()
           qvten3     => null()
           qrten3     => null()
           qsten3     => null()
           qgten3     => null()
           qiten3     => null()
           niten3     => null()
           nrten3     => null()
           ncten3     => null()
           qcten3     => null()
         end if unset_extended_diagnostic_pointers

      end subroutine mp_thompson_run
!>@}

!> \section arg_table_mp_thompson_finalize Argument Table
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

      subroutine get_niwfa(aerfld, nifa, nwfa, ncol, nlev)
         ! To calculate nifa and nwfa from bins of aerosols.
         ! In GOCART and MERRA2, aerosols are given as mixing ratio (kg/kg). To
         ! convert from kg/kg to #/kg, the "unit mass" (mass of one particle)
         ! within the mass bins is calculated. A lognormal size distribution
         ! within aerosol bins is used to find the size based upon the median
         ! mass. NIFA is mainly summarized over five dust bins and NWFA over the
         ! other 10 bins. The parameters besides each bins are carefully tuned
         ! for a good performance of the scheme.
         !
         ! The fields for the last index of the aerfld array
         ! are specified as below.
         ! 1: dust bin 1,                     0.1 to 1.0  micrometers
         ! 2: dust bin 2,                     1.0 to 1.8  micrometers
         ! 3: dust bin 3,                     1.8 to 3.0  micrometers
         ! 4: dust bin 4,                     3.0 to 6.0  micrometers
         ! 5: dust bin 5,                     6.0 to 10.0 micrometers
         ! 6: sea salt bin 1,                 0.03 to 0.1 micrometers
         ! 7: sea salt bin 2,                 0.1 to 0.5  micrometers
         ! 8: sea salt bin 3,                 0.5 to 1.5  micrometers 
         ! 9: sea salt bin 4,                 1.5 to 5.0  micrometers
         ! 10: sea salt bin 5,                5.0 to 10.0 micrometers
         ! 11: Sulfate,                       0.35 (mean) micrometers
         ! 15: water-friendly organic carbon, 0.35 (mean) micrometers
         !
         ! Bin densities are as follows:
         ! 1:    dust bin 1:         2500 kg/m2
         ! 2-5:  dust bin 2-5:       2650 kg/m2
         ! 6-10: sea salt bins 6-10: 2200 kg/m2
         ! 11:   sulfate:            1700 kg/m2
         ! 15:   organic carbon:     1800 kg/m2
         
         implicit none
         integer, intent(in)::ncol, nlev
         real (kind=kind_phys), dimension(:,:,:), intent(in)  :: aerfld
         real (kind=kind_phys), dimension(:,:),   intent(out ):: nifa, nwfa

         nifa=(aerfld(:,:,1)/4.0737762+aerfld(:,:,2)/30.459203+aerfld(:,:,3)/153.45048+ &
              aerfld(:,:,4)/1011.5142+ aerfld(:,:,5)/5683.3501)*1.e15

         nwfa=((aerfld(:,:,6)/0.0045435214+aerfld(:,:,7)/0.2907854+aerfld(:,:,8)/12.91224+ &
              aerfld(:,:,9)/206.2216+ aerfld(:,:,10)/4326.23)*1.+aerfld(:,:,11)/0.3053104*5+ &
              aerfld(:,:,15)/0.3232698*1)*1.e15
      end subroutine get_niwfa

end module mp_thompson

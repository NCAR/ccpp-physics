!>\file mp_tempo.F90
!! This file contains aerosol-aware TEMPO MP scheme.


!>\defgroup aatempo Aerosol-Aware TEMPO MP Module
!! This module contains the aerosol-aware TEMPO microphysics scheme.
module mp_tempo

      use mpi_f08
      use machine, only : kind_phys

      !physical constants that are set from the host
      use module_mp_tempo_params, only : pi, lvap0, lfus, lsub, rv, rdry, cp, t0, r_uni, rho_w
      use module_mp_tempo_params, only : roverrv, eps, naccn0, naccn1, nain0, nain1
      use module_mp_tempo_params, only : initialize_parameters
      use module_mp_tempo_cfgs, only : ty_tempo_cfgs
      use module_mp_tempo_driver, only : tempo_init, tempo_run, ty_tempo_driver_diags, tempo_aerosol_surface_emissions

      implicit none

      public :: mp_tempo_init, mp_tempo_run, mp_tempo_final

      private

   contains

!> This subroutine is a wrapper around the actual tempo_init().
!! \section arg_table_mp_tempo_init Argument Table
!! \htmlinclude mp_tempo_init.html
!!
      subroutine mp_tempo_init(ncol, nlev, &
           imp_physics, imp_physics_tempo, &
           mpirank, mpiroot, &
           tgrs, prsl, phil, con_pi, con_hvap, con_hfus, &
           con_rv, con_g, con_rd, con_cp, &
           con_t0c, con_rgas, rhowater, &
           restart, convert_dry_rho, is_aerosol_aware, &
           is_hail_aware, do_sat_adj, semi_sedi, &
           spechum, nwfa, nifa, nwfa2d, nifa2d, &
           tempo_cfgs, is_initialized, errmsg, errflg)
         
         ! Interface variables
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         logical,                   intent(in   ) :: restart
         logical,                   intent(inout) :: is_initialized
         integer,                   intent(in   ) :: imp_physics
         integer,                   intent(in   ) :: imp_physics_tempo
         logical,                   intent(in   ) :: do_sat_adj
         logical,                   intent(in   ) :: semi_sedi
         logical,                   intent(in   ) :: convert_dry_rho
         logical,                   intent(in   ) :: is_aerosol_aware
         logical,                   intent(in   ) :: is_hail_aware
         real(kind_phys),           intent(in   ) :: con_pi, con_hvap, con_hfus, &
                                                     con_rv, con_g, con_rd, con_cp, &
                                                     con_t0c, con_rgas, rhowater
         ! Hydrometeors
         real(kind_phys),           intent(inout) :: spechum(:,:)
         ! Aerosols
         real(kind_phys),           intent(inout), optional :: nwfa(:,:)
         real(kind_phys),           intent(inout), optional :: nifa(:,:)
         real(kind_phys),           intent(inout), optional :: nwfa2d(:)
         real(kind_phys),           intent(inout), optional :: nifa2d(:)

         ! State variables
         real(kind_phys),           intent(in   ) :: tgrs(:,:)
         real(kind_phys),           intent(in   ) :: prsl(:,:)
         real(kind_phys),           intent(in   ) :: phil(:,:)
         ! MPI information
         integer,                   intent(in   ) :: mpirank
         integer,                   intent(in   ) :: mpiroot
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg
         type(ty_tempo_cfgs),       intent(inout) :: tempo_cfgs
         
         real(kind_phys) :: qv(1:ncol,1:nlev)       ! kg kg-1 (water vapor mixing ratio)
         real(kind_phys) :: hgt(1:ncol,1:nlev)      ! m
         real(kind_phys) :: rho(1:ncol,1:nlev)      ! kg m-3
         real(kind_phys) :: orho(1:ncol,1:nlev)     ! m3 kg-1
         
         real (kind=kind_phys) :: h_01, z1, niIN3, niCCN3
         integer :: i, k
         
         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (do_sat_adj) then
            if ((is_aerosol_aware) .or. (is_hail_aware)) then
               write(errmsg, fmt='((a))') 'do_sat_adj should be run with is_aerosol_aware=F and is_hail_aware=F'
               errflg = 1
               return
            endif
         end if

         if (is_initialized) return
         
         ! Consistency checks
         if (imp_physics/=imp_physics_tempo) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Tempo MP"
            errflg = 1
            return
         end if

         ! Call tempo init (also sets initial default values of physical constants)
         if (mpirank==mpiroot) write(*,*) 'Calling tempo_init() with ltaerosol= ', is_aerosol_aware, &
              ' lthailaware= ', is_hail_aware, ' sedi_semi= ', semi_sedi, ' do_sat_adj= ', do_sat_adj

         ! Main call to tempo_init()
         call tempo_init(aerosolaware_flag=is_aerosol_aware, hailaware_flag=is_hail_aware, &
              semi_sedi_flag=semi_sedi, cloud_condensation_flag=(.not. do_sat_adj), &
              tempo_cfgs=tempo_cfgs)

         ! Set local TEMPO MP module constants from host model and overwrite derived constants calculated in module_mp_tempo_params/initialize_parameters()
         pi = con_pi
         lvap0 = con_hvap
         lfus = con_hfus
         lsub = lvap0 + lfus
         rv   = con_rv
         rdry = con_rd
         cp = con_cp
         t0 = con_t0c
         r_uni = con_rgas
         rho_w = rhowater

         ! Although initialize_parameters() is already called during the call to tempo_init() above, it needs to be called again with the host-set constants to recalculate dependent parameters
         call initialize_parameters()

         if (errflg /= 0) return

         ! For restart runs, the init is done here
         if (restart) then
           is_initialized = .true.
           return
         end if

         where(spechum<0) spechum = 1.0e-10
         qv = spechum/(1.0_kind_phys-spechum)         
         if (convert_dry_rho) then
           if (is_aerosol_aware) then
              nwfa = nwfa/(1.0_kind_phys-spechum)
              nifa = nifa/(1.0_kind_phys-spechum)
           end if
         end if

         ! Geopotential height in m2 s-2 to height in m
         hgt = phil/con_g
         
         ! Density of moist air in kg m-3 and inverse density of air
         rho = roverrv*prsl/(rdry*tgrs*(qv+roverrv))
         orho = 1.0/rho

         ! Check for existing aerosol data, both CCN and IN aerosols.  If missing
         ! fill in just a basic vertical profile, somewhat boundary-layer following.
         if (is_aerosol_aware) then

           ! Potential cloud condensation nuclei (CCN)
           if (MAXVAL(nwfa) .lt. eps) then
             if (mpirank==mpiroot) write(*,*) ' There are no initial CCN aerosols. A basic vertical profile will be created.'
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
           else
             if (mpirank==mpiroot) write(*,*) ' Initial CCN aerosols are present.'
             if (MAXVAL(nwfa2d) .lt. eps) then
               !+---+-----------------------------------------------------------------+
               !..Scale the lowest level aerosol data into an emissions rate.  This is
               !.. very far from ideal, but need higher emissions where larger amount
               !.. of (climo) existing and lesser emissions where there exists fewer to
               !.. begin as a first-order simplistic approach.  Later, proper connection to
               !.. emission inventory would be better.
               !+---+-----------------------------------------------------------------+
               if (mpirank==mpiroot) write(*,*) ' There are no initial CCN aerosol surface emission rates. Rates will be created from surface values.'
               do i = 1, ncol
                  z1 = hgt(i,2)-hgt(i,1)
                  nwfa2d(i) = nwfa(i,1) * 0.000196 * (5./z1)
               enddo
             else
                if (mpirank==mpiroot) write(*,*) ' Initial CCN aerosol surface emission rates are present.'
             endif
           endif

           ! Potential ice nuclei (IN)
           if (MAXVAL(nifa) .lt. eps) then
             if (mpirank==mpiroot) write(*,*) ' There are no initial IN aerosols. A basic vertical profile will be created.'
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
             if (mpirank==mpiroot) write(*,*) ' Initial IN aerosols are present.'
             if (MAXVAL(nifa2d) .lt. eps) then
               if (mpirank==mpiroot) write(*,*) ' There are no initial IN aerosol surface emission rates. Rates will be set to zero.'
               ! calculate IN surface flux here, right now just set to zero
               nifa2d = 0.
             else
               if (mpirank==mpiroot) write(*,*) ' Initial IN aerosol surface emission rates are present.'
             endif
           endif

           ! Ensure non-negative aerosol number concentrations.
           where(nwfa .LE. 0.0) nwfa = 1.1E6
           where(nifa .LE. 0.0) nifa = naIN1*0.01
         end if

         if (convert_dry_rho) then
           if (is_aerosol_aware) then
              nwfa = nwfa/(1.0_kind_phys+qv)
              nifa = nifa/(1.0_kind_phys+qv)
           end if
         end if

         is_initialized = .true.

      end subroutine mp_tempo_init


!> \section arg_table_mp_tempo_run Argument Table
!! \htmlinclude mp_tempo_run.html
!!
!>\ingroup aatempo
!>\section gen_tempo TEMPO MP General Algorithm
      subroutine mp_tempo_run(ncol, nlev, &
        convert_dry_rho, dtp, dt_inner, &
        spechum, qc, qr, qi, qs, qg, ni, nr, &
        nc, nwfa, nifa, nwfa2d, nifa2d, ng, volg, &
        con_g, first_time_step, &
        tgrs, prsl, phii, omega, &
        is_aerosol_aware, is_hail_aware, &
        prcp, rain, graupel, ice, snow, sr, refl_10cm, &
        do_radar_ref, &
        is_initialized, tempo_cfgs, ten_q, ten_t, ten_u, ten_v, &
        dspechum, dqc, dqr, dqi, dqs, dqg, dni, dnr, dnc, dnwfa, &
        dnifa, dng, dvolg, errmsg, errflg)


         ! Interface variables
         logical,                   intent(in   ) :: is_initialized
         logical,                   intent(in   ) :: convert_dry_rho
         logical,                   intent(in   ) :: do_radar_ref
         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         real(kind_phys),           intent(in   ) :: con_g
         ! Hydrometeors
         real(kind_phys),           intent(in) :: spechum(:,:)
         real(kind_phys),           intent(in) :: qc(:,:)
         real(kind_phys),           intent(in) :: qr(:,:)
         real(kind_phys),           intent(in) :: qi(:,:)
         real(kind_phys),           intent(in) :: qs(:,:)
         real(kind_phys),           intent(in) :: qg(:,:)
         real(kind_phys),           intent(in) :: ni(:,:)
         real(kind_phys),           intent(in) :: nr(:,:)
         real(kind_phys), optional, intent(in) :: nc(:,:)
         real(kind_phys), optional, intent(in) :: nwfa(:,:)
         real(kind_phys), optional, intent(in) :: nifa(:,:)
         real(kind_phys), optional, intent(in   ) :: nwfa2d(:)
         real(kind_phys), optional, intent(in   ) :: nifa2d(:)
         real(kind_phys), optional, intent(in) :: ng(:,:)
         real(kind_phys), optional, intent(in) :: volg(:,:)
         logical,                   intent(in)    :: is_aerosol_aware
         logical,                   intent(in)    :: is_hail_aware
         ! Precip/rain/snow/graupel fall amounts and fraction of frozen precip
         real(kind_phys),           intent(inout) :: prcp(:)
         real(kind_phys),           intent(inout) :: rain(:)
         real(kind_phys),           intent(inout) :: graupel(:)
         real(kind_phys),           intent(inout) :: ice(:)
         real(kind_phys),           intent(inout) :: snow(:)
         real(kind_phys),           intent(  out) :: sr(:)
         ! Radar reflectivity
         real(kind_phys),           intent(inout) :: refl_10cm(:,:)         
         ! State variables and timestep information
         real(kind_phys),           intent(in   ) :: tgrs(:,:)
         real(kind_phys),           intent(in   ) :: prsl(:,:)
         real(kind_phys),           intent(in   ) :: phii(:,:)
         real(kind_phys),           intent(in   ) :: omega(:,:)
         real(kind_phys),           intent(in   ) :: dtp
         real(kind=kind_phys),      intent(in   ) :: dt_inner
         logical,                   intent(in   ) :: first_time_step

         real(kind_phys),           intent(  out) :: ten_q(:,:,:)
         real(kind_phys),           intent(  out) :: ten_t(:,:)
         real(kind_phys),           intent(  out) :: ten_u(:,:)
         real(kind_phys),           intent(  out) :: ten_v(:,:)
         real(kind_phys),           intent(  out) :: dspechum(:,:)
         real(kind_phys),           intent(  out) :: dqc(:,:)
         real(kind_phys),           intent(  out) :: dqr(:,:)
         real(kind_phys),           intent(  out) :: dqi(:,:)
         real(kind_phys),           intent(  out) :: dqs(:,:)
         real(kind_phys),           intent(  out) :: dqg(:,:)
         real(kind_phys),           intent(  out) :: dni(:,:)
         real(kind_phys),           intent(  out) :: dnr(:,:)
         real(kind_phys), optional, intent(  out) :: dnc(:,:)
         real(kind_phys), optional, intent(  out) :: dnwfa(:,:)
         real(kind_phys), optional, intent(  out) :: dnifa(:,:)
         real(kind_phys), optional, intent(  out) :: dng(:,:)
         real(kind_phys), optional, intent(  out) :: dvolg(:,:)
         
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
         type(ty_tempo_driver_diags) :: tempo_driver_diags
         ! Local variables

         ! Reduced time step if dt_inner
         real(kind_phys) :: dt
         ! Air density
         real(kind_phys) :: rho(1:ncol,1:nlev)              !< kg m-3
         ! Water vapor mixing ratio (instead of specific humidity)
         real(kind_phys) :: qv(1:ncol,1:nlev)               !< kg kg-1
         ! Vertical velocity and level width
         real(kind_phys) :: w(1:ncol,1:nlev)                !< m s-1
         real(kind_phys) :: dz(1:ncol,1:nlev)               !< m
         real(kind_phys) :: xnwfa(1:ncol,1:nlev,1)
         real(kind_phys) :: xnwfa2d(1:ncol,1)
         
         !temporary new states used to calculate tendencies
         real(kind_phys) :: new_spechum(1:ncol,1:nlev)
         real(kind_phys) :: new_qc(1:ncol,1:nlev)
         real(kind_phys) :: new_qr(1:ncol,1:nlev)
         real(kind_phys) :: new_qi(1:ncol,1:nlev)
         real(kind_phys) :: new_qs(1:ncol,1:nlev)
         real(kind_phys) :: new_qg(1:ncol,1:nlev)
         real(kind_phys) :: new_ni(1:ncol,1:nlev)
         real(kind_phys) :: new_nr(1:ncol,1:nlev)
         real(kind_phys), allocatable :: new_nc(:,:)
         real(kind_phys), allocatable :: new_nwfa(:,:)
         real(kind_phys), allocatable :: new_nifa(:,:)
         real(kind_phys), allocatable :: new_ng(:,:)
         real(kind_phys), allocatable :: new_volg(:,:)
         real(kind_phys) :: new_tgrs(1:ncol,1:nlev)

         ! Dimensions
         integer :: ndt, i, k, it
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte
         integer :: itimestep = 1
         
         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ten_q    = 0.0 ! Since this scheme is outputting tracer tendencies individually,
                        ! we also need to initialize the entire array to 0, so that when
                        ! tendencies are applied, all tracer tendencies other than those
                        ! set in this scheme are 0.
         ten_t    = 0.0
         ten_u    = 0.0
         ten_v    = 0.0
         
         dspechum = 0.0
         dqc      = 0.0
         dqr      = 0.0
         dqi      = 0.0
         dqs      = 0.0
         dqg      = 0.0
         dni      = 0.0
         dnr      = 0.0
         
         new_spechum = spechum
         new_qc = qc
         new_qr = qr
         new_qi = qi
         new_qs = qs
         new_qg = qg
         new_ni = ni
         new_nr = nr
         new_tgrs = tgrs

         if (is_aerosol_aware) then
           dnc      = 0.0
           dnwfa    = 0.0
           dnifa    = 0.0
           
           allocate(new_nc(ncol,nlev))
           allocate(new_nwfa(ncol,nlev))
           allocate(new_nifa(ncol,nlev))
           new_nc   = nc
           new_nwfa = nwfa
           new_nifa = nifa
         endif

         if (is_hail_aware) then
           dng = 0.0
           dvolg  = 0.0
           
           allocate(new_ng(ncol,nlev))
           allocate(new_volg(ncol,nlev))
           new_ng = ng
           new_volg  = volg
         endif

         if (first_time_step) then
            ! Check initialization state
            if (.not.is_initialized) then
               write(errmsg, fmt='((a))') 'mp_tempo_run called before mp_tempo_init'
               errflg = 1
               return
            end if
         endif
         
         ndt = max(nint(dtp/dt_inner), 1)
         dt = dtp/ndt
         if (dt <= dt_inner) dt = dt_inner

         !> - Convert specific humidity to water vapor mixing ratio.
         !> - Also, hydrometeor variables are mass or number mixing ratio
         !> - either kg of species per kg of dry air, or per kg of (dry + vapor).
         qv = new_spechum/(1.0_kind_phys-new_spechum)

         if (convert_dry_rho) then
           new_qc = new_qc/(1.0_kind_phys-new_spechum)
           new_qr = new_qr/(1.0_kind_phys-new_spechum)
           new_qi = new_qi/(1.0_kind_phys-new_spechum)
           new_qs = new_qs/(1.0_kind_phys-new_spechum)
           new_qg = new_qg/(1.0_kind_phys-new_spechum)
           new_ni = new_ni/(1.0_kind_phys-new_spechum)
           new_nr = new_nr/(1.0_kind_phys-new_spechum)
           if (is_hail_aware) then
              new_ng = new_ng/(1.0_kind_phys-new_spechum)
              new_volg = new_volg/(1.0_kind_phys-new_spechum)
           endif
           if (is_aerosol_aware) then
              new_nc = new_nc/(1.0_kind_phys-new_spechum)
              new_nwfa = new_nwfa/(1.0_kind_phys-new_spechum)
              new_nifa = new_nifa/(1.0_kind_phys-new_spechum)
           end if
         end if

         !> - Layer width in m from geopotential in m2 s-2
         dz = (phii(:,2:nlev+1) - phii(:,1:nlev)) / con_g

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

         ! handle dt_inner < dtp
         do it = 1, ndt

            !> - Density of air in kg m-3
            rho = roverrv*prsl/(rdry*new_tgrs*(qv+roverrv))

            !> - Convert omega in Pa s-1 to vertical velocity w in m s-1
            w = -omega/(rho*con_g)

            if (present(nwfa) .and. present(nwfa2d)) then
               xnwfa(:,:,1) = nwfa(:,:)
               xnwfa2d(:,1) = nwfa2d(:)
               call tempo_aerosol_surface_emissions(dt=dt, nwfa=xnwfa, nwfa2d=xnwfa2d, ims=ims, ime=ime, &
                    jms=jms, jme=jme, kms=kms, kme=kme, kts=kts)
               new_nwfa(:,:) = xnwfa(:,:,1)
            endif
            
            call tempo_run(tempo_cfgs=tempo_cfgs, &
                 dt=dt, itimestep=itimestep , &
                 qv=qv, qc=new_qc, qr=new_qr, qi=new_qi, qs=new_qs, qg=new_qg, ni=new_ni, nr=new_nr, &
                 nc=new_nc, nwfa=new_nwfa, nifa=new_nifa, &
                 ng=new_ng, qb=new_volg, &
                 w=w, t=new_tgrs, p=prsl, dz=dz, &
                 ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
                 ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme , &
                 its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte , &
                 tempo_diags=tempo_driver_diags)
            
            ice = ice + max(0.0, tempo_driver_diags%ice_liquid_equiv_precip(:,1)/1000.0_kind_phys)
            snow = snow + (max(0.0, tempo_driver_diags%ice_liquid_equiv_precip(:,1)) + &
                 max(0.0, tempo_driver_diags%snow_liquid_equiv_precip(:,1)))/1000.0_kind_phys
            graupel = graupel + max(0.0, tempo_driver_diags%graupel_liquid_equiv_precip(:,1)/1000.0_kind_phys)
            rain = rain + max(0.0, tempo_driver_diags%rain_precip(:,1)/1000.0_kind_phys)
            prcp = prcp + (max(0.0, tempo_driver_diags%ice_liquid_equiv_precip(:,1)) + &
                 max(0.0, tempo_driver_diags%snow_liquid_equiv_precip(:,1)) + &
                 max(0.0, tempo_driver_diags%graupel_liquid_equiv_precip(:,1)) + &
                 max(0.0, tempo_driver_diags%rain_precip(:,1)))/1000._kind_phys
         enddo

         ! diagnostics that are not precipitation don't need to be in the inner time loop
         sr = tempo_driver_diags%frozen_fraction(:,1)

         if (do_radar_ref) then
            refl_10cm = tempo_driver_diags%refl10cm(:,:,1)
         endif

         itimestep = itimestep + 1
         
         if (errflg/=0) return

         !> - Convert water vapor mixing ratio back to specific humidity
         new_spechum = qv/(1.0_kind_phys+qv)

         if (convert_dry_rho) then
           new_qc = new_qc/(1.0_kind_phys+qv)
           new_qr = new_qr/(1.0_kind_phys+qv)
           new_qi = new_qi/(1.0_kind_phys+qv)
           new_qs = new_qs/(1.0_kind_phys+qv)
           new_qg = new_qg/(1.0_kind_phys+qv)
           new_ni = new_ni/(1.0_kind_phys+qv)
           new_nr = new_nr/(1.0_kind_phys+qv)
           if (is_hail_aware) then
              new_ng = new_ng/(1.0_kind_phys+qv)
              new_volg = new_volg/(1.0_kind_phys+qv)
           endif
           if (is_aerosol_aware) then
              new_nc = new_nc/(1.0_kind_phys+qv)
              new_nwfa = new_nwfa/(1.0_kind_phys+qv)
              new_nifa = new_nifa/(1.0_kind_phys+qv)
           end if
         end if

         dspechum = (new_spechum - spechum)/dtp
         dqc = (new_qc - qc)/dtp
         dqr = (new_qr - qr)/dtp
         dqi = (new_qi - qi)/dtp
         dqs = (new_qs - qs)/dtp
         dqg = (new_qg - qg)/dtp
         dni = (new_ni - ni)/dtp
         dnr = (new_nr - nr)/dtp
         ten_t = (new_tgrs - tgrs)/dtp
         if (is_hail_aware) then
           dng = (new_ng - ng)/dtp
           dvolg  = (new_volg - volg)/dtp
           
           deallocate (new_ng, new_volg)
         end if
         if (is_aerosol_aware) then
           dnc = (new_nc - nc)/dtp
           dnwfa = (new_nwfa - nwfa)/dtp
           dnifa = (new_nifa - nifa)/dtp
           
           deallocate(new_nc, new_nwfa, new_nifa)
         end if

      end subroutine mp_tempo_run

!> \section arg_table_mp_tempo_final Argument Table
!! \htmlinclude mp_tempo_final.html
!!
      subroutine mp_tempo_final(is_initialized, errmsg, errflg)
        
         logical,                   intent(inout) :: is_initialized
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (.not.is_initialized) return

         is_initialized = .false.

      end subroutine mp_tempo_final

end module mp_tempo

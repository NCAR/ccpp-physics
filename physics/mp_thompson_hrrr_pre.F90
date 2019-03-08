! CCPP license goes here, as well as further documentation

!#define DEBUG_AEROSOLS

module mp_thompson_hrrr_pre

      use machine, only : kind_phys

      use module_mp_thompson_hrrr, only : naIN0, naIN1, naCCN0, naCCN1, eps

      implicit none

      public :: mp_thompson_hrrr_pre_init, mp_thompson_hrrr_pre_run, mp_thompson_hrrr_pre_finalize

      private

   contains

      subroutine mp_thompson_hrrr_pre_init()
      end subroutine mp_thompson_hrrr_pre_init

#if 0
!! \section arg_table_mp_thompson_hrrr_pre_run Argument Table
!! | local_name      | standard_name                                                         | long_name                                                | units      | rank | type      |    kind   | intent | optional |
!! |-----------------|-----------------------------------------------------------------------|----------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ncol            | horizontal_loop_extent                                                | horizontal loop extent                                   | count      |    0 | integer   |           | in     | F        |
!! | nlev            | vertical_dimension                                                    | number of vertical levels                                | count      |    0 | integer   |           | in     | F        |
!! | kdt             | index_of_time_step                                                    | current forecast iteration                               | index      |    0 | integer   |           | in     | F        |
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
!! | nwfa2d          | tendency_of_water_friendly_aerosols_at_surface                        | instantaneous fake water-friendly surface aerosol source | kg-1 s-1   |    1 | real      | kind_phys | inout  | T        |
!! | nifa2d          | tendency_of_ice_friendly_aerosols_at_surface                          | instantaneous fake ice-friendly surface aerosol source   | kg-1 s-1   |    1 | real      | kind_phys | inout  | T        |
!! | tgrs            | air_temperature_updated_by_physics                                    | model layer mean temperature                             | K          |    2 | real      | kind_phys | in     | F        |
!! | tgrs_save       | air_temperature_save                                                  | air temperature before entering a physics scheme         | K          |    2 | real      | kind_phys | out    | F        |
!! | prsl            | air_pressure                                                          | mean layer pressure                                      | Pa         |    2 | real      | kind_phys | in     | F        |
!! | phil            | geopotential                                                          | geopotential at model layer centers                      | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | area            | cell_area                                                             | area of the grid cell                                    | m2         |    1 | real      | kind_phys | in     | F        |
!! | mpicomm         | mpi_comm                                                              | MPI communicator                                         | index      |    0 | integer   |           | in     | F        |
!! | mpirank         | mpi_rank                                                              | current MPI-rank                                         | index      |    0 | integer   |           | in     | F        |
!! | mpiroot         | mpi_root                                                              | master MPI-rank                                          | index      |    0 | integer   |           | in     | F        |
!! | blkno           | ccpp_block_number                                                     | for explicit data blocking: block number of this block   | index      |    0 | integer   |           | in     | F        |
!! | errmsg          | ccpp_error_message                                                    | error message for error handling in CCPP                 | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                                       | error flag for error handling in CCPP                    | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine mp_thompson_hrrr_pre_run(ncol, nlev, kdt, con_g, con_rd,    &
                                  spechum, qc, qr, qi, qs, qg, ni, nr,       &
                                  is_aerosol_aware, nc, nwfa, nifa, nwfa2d,  &
                                  nifa2d, tgrs, tgrs_save, prsl, phil, area, &
                                  mpicomm, mpirank, mpiroot, blkno,          &
                                  errmsg, errflg)

         implicit none

         ! Interface variables
         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         integer,                   intent(in   ) :: kdt
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
         logical,                   intent(in   ) :: is_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nc(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nwfa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nifa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nwfa2d(1:ncol)
         real(kind_phys), optional, intent(inout) :: nifa2d(1:ncol)
         ! State variables and timestep information
         real(kind_phys),           intent(in   ) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(  out) :: tgrs_save(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: phil(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: area(1:ncol)
         ! MPI information
         integer,                   intent(in   ) :: mpicomm
         integer,                   intent(in   ) :: mpirank
         integer,                   intent(in   ) :: mpiroot
         ! Blocking information
         integer,                   intent(in   ) :: blkno
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables
         real (kind=kind_phys) :: hgt(1:ncol,1:nlev)  ! m
         real (kind=kind_phys) :: rho(1:ncol,1:nlev)  ! kg m-3
         real (kind=kind_phys) :: orho(1:ncol,1:nlev) ! m3 kg-1
         real (kind=kind_phys) :: h_01, airmass, niIN3, niCCN3
         integer :: i, k

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Return if not first timestep
         if (kdt > 1) return

         ! Fix initial values of hydrometeors
         where(spechum<0) spechum = 0.0
         where(qc<0)      qc = 0.0
         where(qr<0)      qr = 0.0
         where(qi<0)      qi = 0.0
         where(qs<0)      qs = 0.0
         where(qg<0)      qg = 0.0
         where(ni<0)      ni = 0.0
         where(nr<0)      nr = 0.0
         ! If qi is in boundary conditions but ni is not, reset qi to zero (and vice versa)
         if (maxval(qi)>0.0 .and. maxval(ni)==0.0) qi = 0.0
         if (maxval(ni)>0.0 .and. maxval(qi)==0.0) ni = 0.0
         ! If qr is in boundary conditions but nr is not, reset qr to zero (and vice versa)
         if (maxval(qr)>0.0 .and. maxval(nr)==0.0) qr = 0.0
         if (maxval(nr)>0.0 .and. maxval(qr)==0.0) nr = 0.0

         ! Return if aerosol-aware option is not used
         if (.not. is_aerosol_aware) return

         if (.not.present(nc)     .or. &
             .not.present(nwfa2d) .or. &
             .not.present(nifa2d) .or. &
             .not.present(nwfa)   .or. &
             .not.present(nifa)        ) then
             write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_hrrr_pre_run:',                 &
                                        ' aerosol-aware microphysics require all of the following', &
                                        ' optional arguments: nc, nifa2d, nwfa2d, nwfa, nifa'
             errflg = 1
             return
          end if

          ! Fix initial values of aerosols
          where(nc<0)     nc = 0.0
          where(nwfa<0)   nwfa = 0.0
          where(nifa<0)   nifa = 0.0
          where(nwfa2d<0) nwfa2d = 0.0
          where(nifa2d<0) nifa2d = 0.0

#ifdef DEBUG_AEROSOLS
         if (mpirank==mpiroot) then
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp_thompson_hrrr_pre_run before: nwfa2d min/mean/max =", &
                                 & minval(nwfa2d), sum(nwfa2d)/real(size(nwfa2d)), maxval(nwfa2d)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp_thompson_hrrr_pre_run before: nifa2d min/mean/max =", &
                                 & minval(nifa2d), sum(nifa2d)/real(size(nifa2d)), maxval(nifa2d)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp_thompson_hrrr_pre_run before: nwfa min/mean/max =",   &
                                 & minval(nwfa), sum(nwfa)/real(size(nwfa)), maxval(nwfa)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp_thompson_hrrr_pre_run before: nifa min/mean/max =",   &
                                 & minval(nifa), sum(nifa)/real(size(nifa)), maxval(nifa)
         end if
#endif

         ! Save current air temperature for tendency limiters in mp_thompson_hrrr_post
         tgrs_save = tgrs

         ! Geopotential height in m2 s-2 to height in m
         hgt = phil/con_g

         ! Density of air in kg m-3 and inverse density of air
         rho = prsl/(con_rd*tgrs)
         orho = 1.0/rho

!..Check for existing aerosol data, both CCN and IN aerosols.  If missing
!.. fill in just a basic vertical profile, somewhat boundary-layer following.

!.. CCN
         if (MAXVAL(nwfa) .lt. eps) then
            if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial CCN aerosols.'
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
               nwfa2d(i) = nwfa(i,1) * 0.000196 * (airmass*2.E-10)
               do k = 2, nlev
                  nwfa(i,k) = naCCN1+naCCN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niCCN3)
               enddo
            enddo
         else
            if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently initial CCN aerosols are present.'
            if (MAXVAL(nwfa2d) .lt. eps) then
! Hard-coded switch between new (from WRFv4.0, top) and old (until WRFv3.9.1.1, bottom) surface emission rate calculations
#if 0
               !+---+-----------------------------------------------------------------+
               !..Scale the lowest level aerosol data into an emissions rate.  This is
               !.. very far from ideal, but need higher emissions where larger amount
               !.. of (climo) existing and lesser emissions where there exists fewer to
               !.. begin as a first-order simplistic approach.  Later, proper connection to
               !.. emission inventory would be better, but, for now, scale like this:
               !.. where: Nwfa=50 per cc, emit 0.875E4 aerosols per second per grid box unit
               !..        that was tested as ~(20kmx20kmx50m = 2.E10 m**-3)
               !+---+-----------------------------------------------------------------+
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial CCN aerosol surface emission rates.'
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Use new (WRFv4+) formula to calculate CCN surface emission rates.'
               do i = 1, ncol
                  airmass = 1./orho(i,1) * (hgt(i,2)-hgt(i,1))*area(i) ! kg
                  nwfa2d(i) = nwfa(i,1) * 0.000196 * (airmass*2.E-10)
               enddo
#else
#if 1
               !+---+-----------------------------------------------------------------+
               !..Scale the lowest level aerosol data into an emissions rate.  This is
               !.. very far from ideal, but need higher emissions where larger amount
               !.. of existing and lesser emissions where not already lots of aerosols
               !.. for first-order simplistic approach.  Later, proper connection to
               !.. emission inventory would be better, but, for now, scale like this:
               !.. where: Nwfa=50 per cc, emit 0.875E4 aerosols per kg per second
               !..        Nwfa=500 per cc, emit 0.875E5 aerosols per kg per second
               !..        Nwfa=5000 per cc, emit 0.875E6 aerosols per kg per second
               !.. for a grid with 20km spacing and scale accordingly for other spacings.
               !+---+-----------------------------------------------------------------+
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial CCN aerosol surface emission rates.'
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Use old (pre WRFv4) formula to calculate CCN surface emission rates.'
               do i = 1, ncol
                  if (SQRT(area(i))/20000.0 .ge. 1.0) then
                     h_01 = 0.875
                  else
                     h_01 = (0.875 + 0.125*((20000.-SQRT(area(i)))/16000.)) * SQRT(area(i))/20000.
                  endif
                  nwfa2d(i) = 10.0**(LOG10(nwfa(i,1)*1.E-6)-3.69897)
                  nwfa2d(i) = nwfa2d(i)*h_01 * 1.E6
               enddo
#else
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial CCN aerosol surface emission rates, set to zero.'
               ! calculate CCN surface flux here, right now just set to zero
               nwfa2d = 0.
#endif
#endif
            else
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently initial CCN aerosol surface emission rates are present.'
            endif
         endif

!.. IN
         if (MAXVAL(nifa) .lt. eps) then
            if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial IN aerosols.'
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
            if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently initial IN aerosols are present.' 
            if (MAXVAL(nifa2d) .lt. eps) then
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial IN aerosol surface emission rates, set to zero.'
               ! calculate IN surface flux here, right now just set to zero
               nifa2d = 0.
            else
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently initial IN aerosol surface emission rates are present.'
            endif
         endif

#ifdef DEBUG_AEROSOLS
         if (mpirank==mpiroot) then
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp_thompson_hrrr_pre_run after: nwfa2d min/mean/max =", &
                                 & minval(nwfa2d), sum(nwfa2d)/real(size(nwfa2d)), maxval(nwfa2d)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp_thompson_hrrr_pre_run after: nifa2d min/mean/max =", &
                                 & minval(nifa2d), sum(nifa2d)/real(size(nifa2d)), maxval(nifa2d)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp_thompson_hrrr_pre_run after: nwfa min/mean/max =",   &
                                 & minval(nwfa), sum(nwfa)/real(size(nwfa)), maxval(nwfa)
             write(0,'(a,3e16.7)') "AEROSOL DEBUG mp_thompson_hrrr_pre_run after: nifa min/mean/max =",   &
                                 & minval(nifa), sum(nifa)/real(size(nifa)), maxval(nifa)
         end if
#endif

      end subroutine mp_thompson_hrrr_pre_run

      subroutine mp_thompson_hrrr_pre_finalize()
      end subroutine mp_thompson_hrrr_pre_finalize

end module mp_thompson_hrrr_pre

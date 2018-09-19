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
!! | is_aerosol_aware| flag_for_aerosol_physics                                              | flag for aerosol-aware physics                           | flag       |    0 | logical   |           | in     | F        |
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
!! | errmsg          | ccpp_error_message                                                    | error message for error handling in CCPP                 | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                                       | error flag for error handling in CCPP                    | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine mp_thompson_hrrr_pre_run(ncol, nlev, kdt, con_g, con_rd,    &
                                  is_aerosol_aware, nwfa, nifa, nwfa2d,      &
                                  nifa2d, tgrs, tgrs_save, prsl, phil, area, &
                                  mpicomm, mpirank, mpiroot, errmsg, errflg)

         implicit none

         ! Interface variables
         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         integer,                   intent(in   ) :: kdt
         real(kind_phys),           intent(in   ) :: con_g
         real(kind_phys),           intent(in   ) :: con_rd
         ! Aerosols
         logical,                   intent(in   ) :: is_aerosol_aware
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

         ! Return if aerosol-aware option is not used
         if (.not. is_aerosol_aware) return

         if (.not.present(nwfa2d) .or. &
             .not.present(nifa2d) .or. &
             .not.present(nwfa)   .or. &
             .not.present(nifa)        ) then
             write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_hrrr_pre_run:',                 &
                                        ' aerosol-aware microphysics require all of the following', &
                                        ' optional arguments: nifa2d, nwfa2d, nwfa, nifa'
             errflg = 1
             return
          end if

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
            write(*,*) ' Apparently there are no initial CCN aerosols.'
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
            write(*,*) ' Apparently initial CCN aerosols are present.' 
            if (MAXVAL(nwfa2d) .lt. eps) then
               write(*,*) ' Apparently there are no initial CCN aerosol surface emission rates.'
               do i = 1, ncol
                  airmass = 1./orho(i,1) * (hgt(i,2)-hgt(i,1))*area(i) ! kg
                  nwfa2d(i) = nwfa(i,1) * 0.000196 * (airmass*2.E-10)
               enddo
            else
               write(*,*) ' Apparently initial CCN aerosol surface emission rates are present.'
            endif
         endif

!.. IN
         if (MAXVAL(nifa) .lt. eps) then
            write(*,*) ' Apparently there are no initial IN aerosols.'
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
            write(*,*) ' Apparently initial IN aerosols are present.' 
            if (MAXVAL(nifa2d) .lt. eps) then
               write(*,*) ' Apparently there are no initial IN aerosol surface emission rates.'
               ! calculate IN surface flux here, right now just set to zero
               nifa2d = 0.
            else
               write(*,*) ' Apparently initial IN aerosol surface emission rates are present.'
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

!>\file  rrfs_smoke_data.F90
!! This file contains data for the RRFS-Smoke modules.

module rrfs_smoke_data
  use machine ,        only : kind_phys
  implicit none
  INTEGER, PARAMETER :: dep_seasons = 5
  INTEGER, PARAMETER :: nlu = 25

  type wesely_pft
    integer          :: npft
    integer          :: months
    INTEGER, pointer :: seasonal_wes(:,:,:,:) => NULL()
  contains
    final :: wesely_pft_destructor
  end type wesely_pft
  
  interface wesely_pft
    procedure :: wesely_pft_constructor
  end interface wesely_pft

!--------------------------------------------------
! many of these parameters will depend on the RADM mechanism!
! if you change it, lets talk about it and get it done!!!
!--------------------------------------------------

  REAL(kind_phys), parameter    :: small_value = 1.e-36
  REAL(kind_phys), parameter    :: large_value = 1.e36

!--------------------------------------------------
! following currently hardwired to USGS
!--------------------------------------------------
  integer, parameter :: isice_temp   = 24
  integer, parameter :: iswater_temp = 16
  integer, parameter :: wrf2mz_lt_map(nlu) = (/ 1, 2, 2, 2, 2, &
                                                4, 3, 3, 3, 3, &
                                                4, 5, 4, 5, 6, &
                                                7, 9, 6, 8, 9, &
                                                6, 6, 8, 0, 0 /)
  real(kind_phys), parameter    :: wh2o = 18.0153
  real(kind_phys), parameter    :: wpan = 121.04793
  real(kind_phys), PARAMETER ::  KARMAN=0.4
  INTEGER,  parameter :: luse2usgs(21) = (/14,13,12,11,15,8,9,10,10,7, &
                                           17,4,1,5,24,19,16,21,22,23,16 /)
  character(len=4), parameter :: mminlu = 'USGS'

  ! integer, parameter :: pan_seasons = 5
  ! integer, parameter :: pan_lands = 11

  type smoke_data
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Taken from dep_simple_mod
    INTEGER :: ixxxlu(nlu)
    REAL(KIND_PHYS) :: kpart(nlu)
    REAL(KIND_PHYS) :: rac(nlu,dep_seasons), rclo(nlu,dep_seasons), rcls(nlu,dep_seasons)
    REAL(KIND_PHYS) :: rgso(nlu,dep_seasons), rgss(nlu,dep_seasons)
    REAL(KIND_PHYS) :: ri(nlu,dep_seasons), rlu(nlu,dep_seasons)
    ! REAL(KIND_PHYS) :: ri_pan(pan_seasons,pan_lands)
    ! never used: real(kind_phys) :: c0_pan(pan_lands)
    ! never used: real(kind_phys) :: k_pan (pan_lands)

    ! never used: integer :: month
    REAL(KIND_PHYS) :: dratio(1000), hstar(1000), hstar4(1000)
    REAL(KIND_PHYS) :: f0(1000), dhr(1000), scpr23(1000)

    ! Note: scpr23 is only read, never written

    ! never used: type(wesely_pft) :: seasonal_pft

    ! never used: logical, pointer :: is_aerosol(:) => NULL()

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Taken from dep_wet_ls_mod
    real(kind_phys), dimension(:), pointer :: alpha => NULL()
  contains
    final :: smoke_data_destructor
    procedure :: dep_init
  end type smoke_data

  interface smoke_data
    procedure :: smoke_data_constructor
  end interface smoke_data

  type(smoke_data), target, private :: private_thread_data
  logical, private :: rrfs_smoke_data_initialized = .false.

  !$OMP THREADPRIVATE(private_thread_data)
  !$OMP THREADPRIVATE(rrfs_smoke_data_initialized)

contains

  function get_thread_smoke_data() result(data)
    implicit none
    class(smoke_data), pointer :: data
    if(.not. rrfs_smoke_data_initialized) then
      private_thread_data = smoke_data()
      rrfs_smoke_data_initialized = .true.
    endif
    data => private_thread_data
  end function get_thread_smoke_data

  subroutine wesely_pft_destructor(this)
    implicit none
    type(wesely_pft) :: this
    if(associated(this%seasonal_wes)) then
      deallocate(this%seasonal_wes)
      nullify(this%seasonal_wes)
    endif
  end subroutine wesely_pft_destructor

  function wesely_pft_constructor() result(this)
    implicit none
    class(wesely_pft), pointer :: this
    nullify(this%seasonal_wes)
  end function wesely_pft_constructor

  function smoke_data_constructor() result(this)
    implicit none
    type(smoke_data) :: this
    ! These are never used:
    ! this%c0_pan = (/ 0.000, 0.006, 0.002, 0.009, 0.015, &
    !               0.006, 0.000, 0.000, 0.000, 0.002, 0.002 /)
    ! this%k_pan = (/ 0.000, 0.010, 0.005, 0.004, 0.003, &
    !              0.005, 0.000, 0.000, 0.000, 0.075, 0.002 /)
    ! this%month = 0
    ! this%seasonal_pft = wesely_pft()
    ! nullify(this%is_aerosol)
    nullify(this%alpha)
    ! This is not called in the original non-thread-safe code:
    ! call this%dep_init()
  end function smoke_data_constructor

  subroutine smoke_data_destructor(this)
    implicit none
    type(smoke_data) :: this
    if(associated(this%alpha)) then
      deallocate(this%alpha)
      nullify(this%alpha)
    endif
    ! Never used:
    ! if(associated(this%is_aerosol)) then
    !   deallocate(this%is_aerosol)
    !   nullify(this%is_aerosolo)
    ! endif
  end subroutine smoke_data_destructor


!      SUBROUTINE dep_init( id, numgas, mminlu_loc, &
!                           ips, ipe, jps, jpe, ide, jde )
      SUBROUTINE dep_init(this,errmsg,errflg)
        ! Lifted out of dep_simple_mod, this initializes
        ! member variables that were module variables in
        ! that module.
!--
        implicit none
        class(smoke_data) :: this
        character(*), intent(inout) :: errmsg
        integer, intent(inout) :: errflg

!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        ! Unused:
        ! integer, intent(in) ::  numgas
        ! integer, intent(in) :: ips, ipe, jps, jpe
        ! integer, intent(in) :: ide, jde
        ! mmin_lu_loc had no definition, but is also unused

!--------------------------------------------------
! .. Local Scalars
!--------------------------------------------------
        INTEGER :: iland, iseason, l
        integer :: iprt
        integer :: astat
        integer :: ncid
        integer :: dimid
        integer :: varid
        integer :: cpos, slen
        integer :: lon_e, lat_e
        integer :: iend, jend
        integer :: chem_opt
        integer, allocatable :: input_wes_seasonal(:,:,:,:)
        REAL(KIND_PHYS) :: sc
        character(len=128) :: err_msg
        character(len=128) :: filename
        character(len=3)   :: id_num
!--------------------------------------------------
! .. Local Arrays
!--------------------------------------------------
        REAL(KIND_PHYS) :: dat1(nlu,dep_seasons), dat2(nlu,dep_seasons),         &
                dat3(nlu,dep_seasons), dat4(nlu,dep_seasons),         &
                dat5(nlu,dep_seasons), dat6(nlu,dep_seasons),         &
                dat7(nlu,dep_seasons)
        ! REAL(KIND_PHYS) :: dat8(pan_seasons,pan_lands)
      chem_opt = chem_opt

!--------------------------------------------------
! .. Data Statements ..
!     THIS%RI for stomatal resistance
!      data ((this%ri(ILAND,ISEASON),ILAND=1,nlu),ISEASON=1,dep_seasons)/0.10E+11, &
        DATA ((dat1(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
          0.60E+02, 0.60E+02, 0.60E+02, 0.60E+02, 0.70E+02, 0.12E+03, &
          0.12E+03, 0.12E+03, 0.12E+03, 0.70E+02, 0.13E+03, 0.70E+02, &
          0.13E+03, 0.10E+03, 0.10E+11, 0.80E+02, 0.10E+03, 0.10E+11, &
          0.80E+02, 0.10E+03, 0.10E+03, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.10E+11, 0.10E+11, &
          0.70E+02, 0.25E+03, 0.50E+03, 0.10E+11, 0.10E+11, 0.50E+03, &
          0.10E+11, 0.10E+11, 0.50E+03, 0.50E+03, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.10E+11, &
          0.10E+11, 0.70E+02, 0.25E+03, 0.50E+03, 0.10E+11, 0.10E+11, &
          0.50E+03, 0.10E+11, 0.10E+11, 0.50E+03, 0.50E+03, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.70E+02, 0.40E+03, 0.80E+03, 0.10E+11, &
          0.10E+11, 0.80E+03, 0.10E+11, 0.10E+11, 0.80E+03, 0.80E+03, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.12E+03, &
          0.12E+03, 0.12E+03, 0.14E+03, 0.24E+03, 0.24E+03, 0.24E+03, &
          0.12E+03, 0.14E+03, 0.25E+03, 0.70E+02, 0.25E+03, 0.19E+03, &
          0.10E+11, 0.16E+03, 0.19E+03, 0.10E+11, 0.16E+03, 0.19E+03, &
          0.19E+03, 0.10E+11, 0.10E+11, 0.10E+11/
! ..
        IF (nlu/=25) THEN
          errmsg='number of land use classifications not correct '
          errflg=1
          return
        END IF
        IF (dep_seasons/=5) THEN
          errmsg='number of dep_seasons not correct '
          errflg=1
          return
        END IF

!     SURFACE RESISTANCE DATA FOR DEPOSITION MODEL OF
!     M. L. WESELY, ATMOSPHERIC ENVIRONMENT 23 (1989) 1293-1304

!     Seasonal categories:
!     1: midsummer with lush vegetation
!     2: autumn with unharvested cropland
!     3: late autumn with frost, no snow
!     4: winter, snow on ground and subfreezing
!     5: transitional spring with partially green short annuals

!     Land use types:
!     USGS type                                Wesely type
!      1: Urban and built-up land              1
!      2: Dryland cropland and pasture         2
!      3: Irrigated cropland and pasture       2
!      4: Mix. dry/irrg. cropland and pasture  2
!      5: Cropland/grassland mosaic            2
!      6: Cropland/woodland mosaic             4
!      7: Grassland                            3
!      8: Shrubland                            3
!      9: Mixed shrubland/grassland            3
!     10: Savanna                              3, always summer
!     11: Deciduous broadleaf forest           4
!     12: Deciduous needleleaf forest          5, autumn and winter modi
!     13: Evergreen broadleaf forest           4, always summer
!     14: Evergreen needleleaf forest          5
!     15: Mixed Forest                         6
!     16: Water Bodies                         7
!     17: Herbaceous wetland                   9
!     18: Wooded wetland                       6
!     19: Barren or sparsely vegetated         8
!     20: Herbaceous Tundra                    9
!     21: Wooded Tundra                        6
!     22: Mixed Tundra                         6
!     23: Bare Ground Tundra                   8
!     24: Snow or Ice                          -, always winter
!     25: No data                              8


!     Order of data:
!      |
!      |   seasonal category
!     \|/
!     ---> landuse type
!     1       2       3       4       5       6       7       8       9
!     THIS%RLU for outer surfaces in the upper canopy
        DO iseason = 1, dep_seasons
           this%ri(1:nlu,iseason) = dat1(1:nlu,iseason)
        END DO
!      data ((this%rlu(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
        DATA ((dat2(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.10E+11, 0.25E+04, 0.20E+04, 0.10E+11, &
          0.25E+04, 0.20E+04, 0.20E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
          0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, 0.90E+04, &
          0.20E+04, 0.40E+04, 0.80E+04, 0.10E+11, 0.90E+04, 0.80E+04, &
          0.10E+11, 0.90E+04, 0.80E+04, 0.80E+04, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
          0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, &
          0.90E+04, 0.20E+04, 0.40E+04, 0.80E+04, 0.10E+11, 0.90E+04, &
          0.80E+04, 0.10E+11, 0.90E+04, 0.80E+04, 0.80E+04, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.20E+04, 0.60E+04, 0.90E+04, 0.10E+11, &
          0.90E+04, 0.90E+04, 0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.40E+04, 0.40E+04, &
          0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, &
          0.20E+04, 0.40E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.30E+04, &
          0.10E+11, 0.40E+04, 0.30E+04, 0.10E+11, 0.40E+04, 0.30E+04, &
          0.30E+04, 0.10E+11, 0.10E+11, 0.10E+11/
        DO iseason = 1, dep_seasons
           this%rlu(1:nlu,iseason) = dat2(1:nlu,iseason)
        END DO
!     THIS%RAC for transfer that depends on canopy height and density
!      data ((this%rac(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+03, &
        DATA ((dat3(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+04, 0.10E+03, &
          0.10E+03, 0.10E+03, 0.10E+03, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.00E+00, 0.30E+03, 0.20E+04, 0.00E+00, &
          0.30E+03, 0.20E+04, 0.20E+04, 0.00E+00, 0.00E+00, 0.00E+00, &
          0.10E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+04, &
          0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.15E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.17E+04, 0.00E+00, 0.20E+03, 0.17E+04, &
          0.00E+00, 0.20E+03, 0.17E+04, 0.17E+04, 0.00E+00, 0.00E+00, &
          0.00E+00, 0.10E+03, 0.10E+02, 0.10E+02, 0.10E+02, 0.10E+02, &
          0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+04, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, 0.00E+00, 0.10E+03, &
          0.15E+04, 0.00E+00, 0.10E+03, 0.15E+04, 0.15E+04, 0.00E+00, &
          0.00E+00, 0.00E+00, 0.10E+03, 0.10E+02, 0.10E+02, 0.10E+02, &
          0.10E+02, 0.10E+04, 0.10E+02, 0.10E+02, 0.10E+02, 0.10E+02, &
          0.10E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, 0.00E+00, &
          0.50E+02, 0.15E+04, 0.00E+00, 0.50E+02, 0.15E+04, 0.15E+04, &
          0.00E+00, 0.00E+00, 0.00E+00, 0.10E+03, 0.50E+02, 0.50E+02, &
          0.50E+02, 0.50E+02, 0.12E+04, 0.80E+02, 0.80E+02, 0.80E+02, &
          0.10E+03, 0.12E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, &
          0.00E+00, 0.20E+03, 0.15E+04, 0.00E+00, 0.20E+03, 0.15E+04, &
          0.15E+04, 0.00E+00, 0.00E+00, 0.00E+00/
        DO iseason = 1, dep_seasons
           this%rac(1:nlu,iseason) = dat3(1:nlu,iseason)
        END DO
!     THIS%RGSS for ground surface  SO2
!      data ((this%rgss(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.40E+03, &
        DATA ((dat4(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.40E+03, &
          0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.50E+03, 0.35E+03, &
          0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, 0.50E+03, 0.50E+03, &
          0.50E+03, 0.10E+03, 0.10E+01, 0.10E+01, 0.10E+03, 0.10E+04, &
          0.10E+01, 0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, 0.10E+04, &
          0.40E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.50E+03, &
          0.35E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, 0.50E+03, &
          0.50E+03, 0.50E+03, 0.10E+03, 0.10E+01, 0.10E+01, 0.10E+03, &
          0.10E+04, 0.10E+01, 0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, &
          0.10E+04, 0.40E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, &
          0.50E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, &
          0.50E+03, 0.50E+03, 0.50E+03, 0.20E+03, 0.10E+01, 0.10E+01, &
          0.20E+03, 0.10E+04, 0.10E+01, 0.20E+03, 0.20E+03, 0.10E+04, &
          0.10E+03, 0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, &
          0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, &
          0.10E+03, 0.10E+03, 0.50E+03, 0.10E+03, 0.10E+03, 0.10E+01, &
          0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, &
          0.10E+04, 0.10E+03, 0.10E+04, 0.50E+03, 0.15E+03, 0.15E+03, &
          0.15E+03, 0.15E+03, 0.50E+03, 0.35E+03, 0.35E+03, 0.35E+03, &
          0.35E+03, 0.50E+03, 0.50E+03, 0.50E+03, 0.50E+03, 0.20E+03, &
          0.10E+01, 0.10E+01, 0.20E+03, 0.10E+04, 0.10E+01, 0.20E+03, &
          0.20E+03, 0.10E+04, 0.10E+03, 0.10E+04/
        DO iseason = 1, dep_seasons
           this%rgss(1:nlu,iseason) = dat4(1:nlu,iseason)
        END DO
!     THIS%RGSO for ground surface  O3
!      data ((this%rgso(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.30E+03, &
        DATA ((dat5(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.30E+03, &
          0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.30E+03, 0.20E+04, 0.10E+04, 0.30E+03, 0.40E+03, &
          0.10E+04, 0.30E+03, 0.30E+03, 0.40E+03, 0.35E+04, 0.40E+03, &
          0.30E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.30E+03, 0.20E+04, 0.80E+03, 0.30E+03, &
          0.40E+03, 0.80E+03, 0.30E+03, 0.30E+03, 0.40E+03, 0.35E+04, &
          0.40E+03, 0.30E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.30E+03, 0.20E+04, 0.10E+04, &
          0.30E+03, 0.40E+03, 0.10E+04, 0.30E+03, 0.30E+03, 0.40E+03, &
          0.35E+04, 0.40E+03, 0.60E+03, 0.35E+04, 0.35E+04, 0.35E+04, &
          0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, &
          0.35E+04, 0.35E+04, 0.20E+03, 0.35E+04, 0.35E+04, 0.20E+04, &
          0.35E+04, 0.35E+04, 0.40E+03, 0.35E+04, 0.35E+04, 0.35E+04, &
          0.40E+03, 0.35E+04, 0.40E+03, 0.30E+03, 0.15E+03, 0.15E+03, &
          0.15E+03, 0.15E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.30E+03, &
          0.20E+04, 0.10E+04, 0.30E+03, 0.40E+03, 0.10E+04, 0.30E+03, &
          0.30E+03, 0.40E+03, 0.35E+04, 0.40E+03/
        DO iseason = 1, dep_seasons
           this%rgso(1:nlu,iseason) = dat5(1:nlu,iseason)
        END DO
!     THIS%RCLS for exposed surfaces in the lower canopy  SO2
!      data ((this%rcls(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
        DATA ((dat6(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.10E+11, 0.25E+04, 0.20E+04, 0.10E+11, &
          0.25E+04, 0.20E+04, 0.20E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
          0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, 0.90E+04, &
          0.20E+04, 0.20E+04, 0.40E+04, 0.10E+11, 0.90E+04, 0.40E+04, &
          0.10E+11, 0.90E+04, 0.40E+04, 0.40E+04, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, &
          0.90E+04, 0.20E+04, 0.30E+04, 0.60E+04, 0.10E+11, 0.90E+04, &
          0.60E+04, 0.10E+11, 0.90E+04, 0.60E+04, 0.60E+04, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.90E+04, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.90E+04, 0.90E+04, 0.20E+04, 0.20E+03, 0.40E+03, 0.10E+11, &
          0.90E+04, 0.40E+03, 0.10E+11, 0.90E+04, 0.40E+03, 0.40E+03, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.40E+04, 0.40E+04, &
          0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, &
          0.20E+04, 0.40E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.30E+04, &
          0.10E+11, 0.40E+04, 0.30E+04, 0.10E+11, 0.40E+04, 0.30E+04, &
          0.30E+04, 0.10E+11, 0.10E+11, 0.10E+11/
        DO iseason = 1, dep_seasons
           this%rcls(1:nlu,iseason) = dat6(1:nlu,iseason)
        END DO
!     THIS%RCLO for exposed surfaces in the lower canopy  O3
!      data ((this%rclo(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
        DATA ((dat7(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
          0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.10E+04, 0.10E+04, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+11, &
          0.10E+04, 0.10E+04, 0.10E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, &
          0.40E+03, 0.40E+03, 0.40E+03, 0.10E+04, 0.40E+03, 0.40E+03, &
          0.10E+04, 0.10E+04, 0.60E+03, 0.10E+11, 0.40E+03, 0.60E+03, &
          0.10E+11, 0.40E+03, 0.60E+03, 0.60E+03, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, 0.10E+04, 0.40E+03, &
          0.40E+03, 0.10E+04, 0.10E+04, 0.60E+03, 0.10E+11, 0.80E+03, &
          0.60E+03, 0.10E+11, 0.80E+03, 0.60E+03, 0.60E+03, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.10E+04, 0.40E+03, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.40E+03, 0.40E+03, 0.10E+04, 0.15E+04, 0.60E+03, 0.10E+11, &
          0.80E+03, 0.60E+03, 0.10E+11, 0.80E+03, 0.60E+03, 0.60E+03, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, &
          0.10E+04, 0.10E+04, 0.50E+03, 0.50E+03, 0.50E+03, 0.50E+03, &
          0.10E+04, 0.50E+03, 0.15E+04, 0.10E+04, 0.15E+04, 0.70E+03, &
          0.10E+11, 0.60E+03, 0.70E+03, 0.10E+11, 0.60E+03, 0.70E+03, &
          0.70E+03, 0.10E+11, 0.10E+11, 0.10E+11/

        DO iseason = 1, dep_seasons
           this%rclo(1:nlu,iseason) = dat7(1:nlu,iseason)
        END DO

    ! data ((dat8(iseason,iland),iseason=1,pan_seasons),iland=1,pan_lands) / &
    !                      1.e36,  60., 120.,  70., 130., 100.,1.e36,1.e36,  80., 100., 150., &
    !                      1.e36,1.e36,1.e36,1.e36, 250., 500.,1.e36,1.e36,1.e36,1.e36,1.e36, &
    !                      1.e36,1.e36,1.e36,1.e36, 250., 500.,1.e36,1.e36,1.e36,1.e36,1.e36, &
    !                      1.e36,1.e36,1.e36,1.e36, 400., 800.,1.e36,1.e36,1.e36,1.e36,1.e36, &
    !                      1.e36, 120., 240., 140., 250., 190.,1.e36,1.e36, 160., 200., 300. /
    !   this%ri_pan(:,:) = dat8(:,:)

!--------------------------------------------------
!	Initialize parameters
!--------------------------------------------------
        this%hstar  = 0.
        this%hstar4 = 0.
        this%dhr    = 0.
        this%f0     = 0.
        this%dratio = 1.0 ! FIXME: IS THIS RIGHT?
        this%scpr23 = 1.0 ! FIXME: IS THIS RIGHT?

!--------------------------------------------------
!     HENRY''S LAW COEFFICIENTS
!     Effective Henry''s law coefficient at pH 7
!     [KH298]=mole/(l atm)
!--------------------------------------------------

!     DATA FOR AEROSOL PARTICLE DEPOSITION FOR THE MODEL OF
!     J. W. ERISMAN, A. VAN PUL AND P. WYERS
!     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607

!     vd = (u* / k) * CORRECTION FACTORS

!     CONSTANT K FOR LANDUSE TYPES:
! urban and built-up land                  
        this%kpart(1) = 500.
! dryland cropland and pasture             
        this%kpart(2) = 500.
! irrigated cropland and pasture           
        this%kpart(3) = 500.
! mixed dryland/irrigated cropland and past
        this%kpart(4) = 500.
! cropland/grassland mosaic                
        this%kpart(5) = 500.
! cropland/woodland mosaic                 
        this%kpart(6) = 100.
! grassland                                
        this%kpart(7) = 500.
! shrubland                                
        this%kpart(8) = 500.
! mixed shrubland/grassland                
        this%kpart(9) = 500.
! savanna                                  
        this%kpart(10) = 500.
! deciduous broadleaf forest               
        this%kpart(11) = 100.
! deciduous needleleaf forest              
        this%kpart(12) = 100.
! evergreen broadleaf forest               
        this%kpart(13) = 100.
! evergreen needleleaf forest              
        this%kpart(14) = 100.
! mixed forest                             
        this%kpart(15) = 100.
! water bodies                             
        this%kpart(16) = 500.
! herbaceous wetland                       
        this%kpart(17) = 500.
! wooded wetland                           
        this%kpart(18) = 500.
! barren or sparsely vegetated             
        this%kpart(19) = 500.
! herbaceous tundra                        
        this%kpart(20) = 500.
! wooded tundra                            
        this%kpart(21) = 100.
! mixed tundra                             
        this%kpart(22) = 500.
! bare ground tundra                       
        this%kpart(23) = 500.
! snow or ice                              
        this%kpart(24) = 500.
!     Comments:
        this%kpart(25) = 500.
!     Erisman et al. (1994) give
!     k = 500 for low vegetation and k = 100 for forests.

!     For desert k = 500 is taken according to measurements
!     on bare soil by
!     J. Fontan, A. Lopez, E. Lamaud and A. Druilhet (1997)
!     Vertical Flux Measurements of the Submicronic Aerosol Particles
!     and Parametrisation of the Dry Deposition Velocity
!     in: Biosphere-Atmosphere Exchange of Pollutants
!     and Trace Substances
!     Editor: S. Slanina. Springer-Verlag Berlin, Heidelberg, 1997
!     pp. 381-390

!     For coniferous forest the Erisman value of  k = 100 is taken.
!     Measurements of Erisman et al. (1997) in a coniferous forest
!     in the Netherlands, lead to values of k between 20 and 38
!     (Atmospheric Environment 31 (1997), 321-332).
!     However, these high values of vd may be reached during
!     instable cases. The eddy correlation measurements
!     of Gallagher et al. (1997) made during the same experiment
!     show for stable cases (L>0) values of k between 200 and 250
!     at minimum (Atmospheric Environment 31 (1997), 359-373).
!     Fontan et al. (1997) found k = 250 in a forest
!     of maritime pine in southwestern France.

!     For gras, model calculations of Davidson et al. support
!     the value of 500.
!     C. I. Davidson, J. M. Miller and M. A. Pleskov
!     The Influence of Surface Structure on Predicted Particles
!     Dry Deposition to Natural Gras Canopies
!     Water, Air, and Soil Pollution 18 (1982) 25-43

!     Snow covered surface: The experiment of Ibrahim et al. (1983)
!     gives k = 436 for 0.7 um diameter particles.
!     The deposition velocity of Milford and Davidson (1987)
!     gives k = 154 for continental sulfate aerosol.
!     M. Ibrahim, L. A. Barrie and F. Fanaki
!     Atmospheric Environment 17 (1983), 781-788

!     J. B. Milford and C. I. Davidson
!     The Sizes of Particulate Sulfate and Nitrate in the Atmosphere
!     - A Review
!     JAPCA 37 (1987), 125-134
! no data                                  
!       WRITE (0,*) ' return from rcread '
!     *********************************************************

!     Simplified landuse scheme for deposition and biogenic emission
!     subroutines
!     (ISWATER and ISICE are already defined elsewhere,
!     therefore water and ice are not considered here)

!     1 urban or bare soil
!     2 agricultural
!     3 grassland
!     4 deciduous forest
!     5 coniferous and mixed forest
!     6 other natural landuse categories


        IF (mminlu=='OLD ') THEN
          this%ixxxlu(1) = 1
          this%ixxxlu(2) = 2
          this%ixxxlu(3) = 3
          this%ixxxlu(4) = 4
          this%ixxxlu(5) = 5
          this%ixxxlu(6) = 5
          this%ixxxlu(7) = 0
          this%ixxxlu(8) = 6
          this%ixxxlu(9) = 1
          this%ixxxlu(10) = 6
          this%ixxxlu(11) = 0
          this%ixxxlu(12) = 4
          this%ixxxlu(13) = 6
        END IF
        IF (mminlu=='USGS') THEN
          this%ixxxlu(1) = 1
          this%ixxxlu(2) = 2
          this%ixxxlu(3) = 2
          this%ixxxlu(4) = 2
          this%ixxxlu(5) = 2
          this%ixxxlu(6) = 4
          this%ixxxlu(7) = 3
          this%ixxxlu(8) = 6
          this%ixxxlu(9) = 3
          this%ixxxlu(10) = 6
          this%ixxxlu(11) = 4
          this%ixxxlu(12) = 5
          this%ixxxlu(13) = 4
          this%ixxxlu(14) = 5
          this%ixxxlu(15) = 5
          this%ixxxlu(16) = 0
          this%ixxxlu(17) = 6
          this%ixxxlu(18) = 4
          this%ixxxlu(19) = 1
          this%ixxxlu(20) = 6
          this%ixxxlu(21) = 4
          this%ixxxlu(22) = 6
          this%ixxxlu(23) = 1
          this%ixxxlu(24) = 0
          this%ixxxlu(25) = 1
        END IF
        IF (mminlu=='SiB ') THEN
          this%ixxxlu(1) = 4
          this%ixxxlu(2) = 4
          this%ixxxlu(3) = 4
          this%ixxxlu(4) = 5
          this%ixxxlu(5) = 5
          this%ixxxlu(6) = 6
          this%ixxxlu(7) = 3
          this%ixxxlu(8) = 6
          this%ixxxlu(9) = 6
          this%ixxxlu(10) = 6
          this%ixxxlu(11) = 1
          this%ixxxlu(12) = 2
          this%ixxxlu(13) = 6
          this%ixxxlu(14) = 1
          this%ixxxlu(15) = 0
          this%ixxxlu(16) = 0
          this%ixxxlu(17) = 1
        END IF

      END SUBROUTINE dep_init
end module rrfs_smoke_data

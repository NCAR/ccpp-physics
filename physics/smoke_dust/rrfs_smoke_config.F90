!>\file  rrfs_smoke_config.F90
!! This file contains the configuration for RRFS-Smoke.
!
! Haiqin.Li@noaa.gov  
! 06/2021
! constant parameters and chemistry configurations and tracers
! (This will be splited into three subroutines for configuration, constant and tracers later)
! 06/2021 move configuration into chem nml
!
module rrfs_smoke_config

  use machine ,       only : kind_phys

  implicit none

  !-- constant paramters
  real(kind=kind_phys), parameter :: epsilc     = 1.e-12
  !-- aerosol module configurations
  integer :: chem_opt = 1   
  integer :: kemit = 1
  integer :: dust_opt = 5
  integer :: seas_opt = 0   ! turn off by default
  logical :: do_plumerise  = .true.
  integer :: addsmoke_flag = 1
  integer :: smoke_forecast = 1
  integer :: plumerisefire_frq=60
  integer :: n_dbg_lines = 3
  integer :: wetdep_ls_opt = 1
  integer :: drydep_opt  = 1
  integer :: pm_settling = 1
  integer :: nfire_types = 5
  integer :: ebb_dcycle  = 2 ! 1: read in ebb_smoke(i,24), 2: daily
  logical :: dbg_opt     = .true.
  logical :: aero_ind_fdb = .false.
  logical :: add_fire_heat_flux= .false.
  logical :: do_rrfs_sd = .true.
!  integer :: wind_eff_opt = 1
  logical :: extended_sd_diags = .false.
  real(kind_phys) :: wetdep_ls_alpha = .5 ! scavenging factor

  ! --
  integer, parameter :: CHEM_OPT_GOCART= 1
  integer, parameter :: num_moist=2, num_chem=20, num_emis_seas=5, num_emis_dust=5

  ! -- hydrometeors
  integer, parameter :: p_qv=1
  integer, parameter :: p_qc=2
  integer, parameter :: p_qi=3
  ! -- set pointers to predefined atmospheric tracers
  ! -- FV3 GFDL microphysics
  integer, parameter :: p_atm_shum = 1
  integer, parameter :: p_atm_cldq = 2

  integer :: numgas = 0

  !-- tracers
   integer, parameter :: p_smoke=5
   integer, parameter :: p_dust_1=10
   integer, parameter :: p_dust_2=11
   integer, parameter :: p_dust_3=12
   integer, parameter :: p_dust_4=13
   integer, parameter :: p_dust_5=14, p_coarse_pm=14
   integer, parameter :: p_seas_1=15
   integer, parameter :: p_seas_2=16
   integer, parameter :: p_seas_3=17
   integer, parameter :: p_seas_4=18
   integer, parameter :: p_seas_5=19

   integer, parameter :: p_edust1=1,p_edust2=2,p_edust3=3,p_edust4=4,p_edust5=5

  ! -- fire options
!  integer, parameter :: num_plume_data = 1


end module

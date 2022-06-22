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

  !-- chemistyr module configurations
  integer :: chem_opt = 1
  integer :: kemit = 1
  integer :: dust_opt = 5
  integer :: dmsemis_opt = 1
  integer :: seas_opt = 2
  integer :: biomass_burn_opt=1
  logical :: do_plumerise  = .true.
  integer :: addsmoke_flag = 1
  integer :: plumerisefire_frq=60  ! Let's add to the namelist
  integer :: chem_conv_tr  = 0
  integer :: aer_ra_feedback=1 !0
  integer :: aer_ra_frq  = 60
  integer :: wetdep_ls_opt = 1
  integer :: drydep_opt  = 1
  logical :: bb_dcycle   = .false.
  logical :: smoke_forecast = .false.
  logical :: aero_ind_fdb = .false.
  logical :: dbg_opt     = .true.

  real(kind=kind_phys), parameter :: depo_fact=0.
  integer, parameter :: CHEM_OPT_GOCART= 1
  INTEGER, PARAMETER :: gocartracm_kpp = 301
  integer, parameter :: chem_tune_tracers = 20
  integer, parameter :: DUST_OPT_NONE = 0
  integer, parameter :: SEAS_OPT_NONE = 0
  ! -- DMS emissions
  integer, parameter :: DMSE_OPT_NONE   = 0
  integer, parameter :: DMSE_OPT_ENABLE = 1
  ! -- subgrid convective transport
  integer, parameter :: CTRA_OPT_NONE  = 0
  integer, parameter :: CTRA_OPT_GRELL = 2
  ! -- large scale wet deposition
  integer, parameter :: WDLS_OPT_NONE  = 0
  integer, parameter :: WDLS_OPT_GSD   = 1
  integer, parameter :: WDLS_OPT_NGAC  = 2

  ! --
  integer, parameter :: call_chemistry     = 1
  integer, parameter :: num_moist=3, num_chem=20, num_emis_seas=5, num_emis_dust=5
  integer, parameter :: num_emis_ant = 7

  integer, parameter :: SEAS_OPT_DEFAULT = 1

  integer, parameter :: DUST_OPT_GOCART  = 1
  integer, parameter :: DUST_OPT_AFWA    = 3
  integer, parameter :: DUST_OPT_FENGSHA = 5

  ! -- biomass burning emissions
  integer, parameter :: BURN_OPT_ENABLE = 1
  integer, parameter :: FIRE_OPT_MODIS  = 1
  integer, parameter :: FIRE_OPT_GBBEPx = 2

  ! -- hydrometeors
  integer, parameter :: p_qv=1
  integer, parameter :: p_qc=2
  integer, parameter :: p_qi=3
  ! -- set pointers to predefined atmospheric tracers
  ! -- FV3 GFDL microphysics
  integer, parameter :: p_atm_shum = 1
  integer, parameter :: p_atm_cldq = 2
  integer, parameter :: p_atm_o3mr = 7

  integer :: numgas = 0

  real(kind=kind_phys) :: wetdep_ls_alpha(chem_tune_tracers)=-999.

  !-- tracers
  integer, parameter :: p_so2=1
  integer, parameter :: p_sulf=2
  integer, parameter :: p_dms=3
  integer, parameter :: p_msa=4
  integer, parameter :: p_p25=5, p_smoke=5
  integer, parameter :: p_bc1=6
  integer, parameter :: p_bc2=7
  integer, parameter :: p_oc1=8
  integer, parameter :: p_oc2=9
  integer, parameter :: p_dust_1=10
  integer, parameter :: p_dust_2=11
  integer, parameter :: p_dust_3=12
  integer, parameter :: p_dust_4=13
  integer, parameter :: p_dust_5=14
  integer, parameter :: p_seas_1=15
  integer, parameter :: p_seas_2=16
  integer, parameter :: p_seas_3=17
  integer, parameter :: p_seas_4=18
  integer, parameter :: p_seas_5=19
  integer, parameter :: p_p10   =20

  integer, parameter :: p_edust1=1,p_edust2=2,p_edust3=3,p_edust4=4,p_edust5=5
  integer, parameter :: p_eseas1=1,p_eseas2=2,p_eseas3=3,p_eseas4=4,p_eseas5=5
 
  integer :: p_ho=0,p_h2o2=0,p_no3=0

  ! constants
  real(kind=kind_phys), PARAMETER :: airmw      = 28.97
  real(kind=kind_phys), PARAMETER :: mw_so2_aer = 64.066
  real(kind=kind_phys), PARAMETER :: mw_so4_aer = 96.066
  real(kind=kind_phys), parameter :: smw        = 32.00
  real(kind=kind_phys), parameter :: mwdry      = 28.
!  <mw>d is the molecular weight of dry air (28.966), <mw>w/<mw>d = 0.62197, and
!  (<mw>d - <mw>w)/<mw>d = 0.37803
!  http://atmos.nmsu.edu/education_and_outreach/encyclopedia/humidity.htm

  ! -- fire options
!  integer, parameter :: num_plume_data = 1


end module

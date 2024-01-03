!>\file bl_mynn_common.f90
!! Define Model-specific constants/parameters.
!! This module will be used at the initialization stage
!! where all model-specific constants are read and saved into
!! memory. This module is then used again in the MYNN-EDMF. All
!! MYNN-specific constants are declared globally in the main 
!! module (module_bl_mynn) further below:

!>\ingroup gp_mynnedmf
!! Define Model-specific constants/parameters
 module bl_mynn_common

!------------------------------------------
!
!------------------------------------------

! The following 5-6 lines are the only lines in this file that are not 
! universal for all dycores... Any ideas how to universalize it?
! For MPAS:
! use mpas_kind_types,only: kind_phys => RKIND
! For CCPP:
  use machine,  only : kind_phys

 implicit none
 save

! To be specified from dycore
 real(kind=kind_phys):: cp           !<= 7.*r_d/2. (J/kg/K)
 real(kind=kind_phys):: cpv          !<= 4.*r_v    (J/kg/K) Spec heat H2O gas
 real(kind=kind_phys):: cice         !<= 2106.     (J/kg/K) Spec heat H2O ice
 real(kind=kind_phys):: cliq         !<= 4190.     (J/kg/K) Spec heat H2O liq
 real(kind=kind_phys):: p608         !<= R_v/R_d-1.
 real(kind=kind_phys):: ep_2         !<= R_d/R_v
 real(kind=kind_phys):: grav         !<= accel due to gravity
 real(kind=kind_phys):: karman       !<= von Karman constant
 real(kind=kind_phys):: t0c          !<= temperature of water at freezing, 273.15 K
 real(kind=kind_phys):: rcp          !<= r_d/cp
 real(kind=kind_phys):: r_d          !<= 287.  (J/kg/K) gas const dry air
 real(kind=kind_phys):: r_v          !<= 461.6 (J/kg/K) gas const water
 real(kind=kind_phys):: xlf          !<= 0.35E6 (J/kg) fusion at 0 C
 real(kind=kind_phys):: xlv          !<= 2.50E6 (J/kg) vaporization at 0 C
 real(kind=kind_phys):: xls          !<= 2.85E6 (J/kg) sublimation
 real(kind=kind_phys):: rvovrd       !<= r_v/r_d != 1.608

! Specified locally
 real(kind=kind_phys),parameter:: zero   = 0.0
 real(kind=kind_phys),parameter:: half   = 0.5
 real(kind=kind_phys),parameter:: one    = 1.0
 real(kind=kind_phys),parameter:: two    = 2.0
 real(kind=kind_phys),parameter:: onethird  = 1./3.
 real(kind=kind_phys),parameter:: twothirds = 2./3.
 real(kind=kind_phys),parameter:: tref  = 300.0   !< reference temperature (K)
 real(kind=kind_phys),parameter:: TKmin = 253.0   !< for total water conversion, Tripoli and Cotton (1981)
 real(kind=kind_phys),parameter:: p1000mb=100000.0
 real(kind=kind_phys),parameter:: svp1  = 0.6112 !<(kPa)
 real(kind=kind_phys),parameter:: svp2  = 17.67  !<(dimensionless)
 real(kind=kind_phys),parameter:: svp3  = 29.65  !<(K)
 real(kind=kind_phys),parameter:: tice  = 240.0  !<-33 (C), temp at saturation w.r.t. ice

! To be derived in the init routine
 real(kind=kind_phys):: ep_3         !<= 1.-ep_2 != 0.378
 real(kind=kind_phys):: gtr          !<= grav/tref
 real(kind=kind_phys):: rk           !<= cp/r_d
 real(kind=kind_phys):: tv0          !<= p608*tref
 real(kind=kind_phys):: tv1          !<= (1.+p608)*tref
 real(kind=kind_phys):: xlscp        !<= (xlv+xlf)/cp
 real(kind=kind_phys):: xlvcp        !<= xlv/cp
 real(kind=kind_phys):: g_inv        !<= 1./grav

 end module bl_mynn_common

!>\file  dep_data_mod.F90
!! This file contains data for the dry deposition modules.
module dep_data_mod

  use machine , only : kind_phys

  integer, parameter :: nvegtype = 25
  real(kind_phys), dimension(nvegtype), parameter :: &
        kpart = (/500., 500., 500., 500., 500., &
                  500., 500., 500., 500., 500., &
                  500., 500., 500., 500., 500., &
                  500., 500., 500., 500., 500., &
                  500., 500., 500., 500., 500.    /)
  real(kind_phys), parameter :: max_dep_vel = 0.005                   ! m/s (may need to set per species)
  real(kind_phys), parameter :: dep_ref_hgt = 2.0                     ! Meters 
  real(kind_phys), parameter :: pi = 3.1415926536
!  3*PI
  REAL(kind_phys), PARAMETER :: threepi=3.0*pi
  real(kind_phys), parameter :: gravity =  9.81
! mean gravitational acceleration [ m/sec**2 ]
  REAL(kind_phys), PARAMETER :: grav=9.80622
  real(kind_phys), parameter :: boltzmann = 1.3807e-16
! universal gas constant [ J/mol-K ]
  REAL(kind_phys), PARAMETER :: rgasuniv=8.314510
! Avogadro's Constant [ 1/mol ]
  REAL, PARAMETER :: avo=6.0221367E23
  ! Boltzmann's Constant [ J / K ]i\
  REAL(kind_phys), PARAMETER :: boltz=rgasuniv/avo
  real(kind_phys), parameter :: Cb = 2., Cim = 0.4, alpha = 0.8, Cin = 2.5, vv = 0.8
  real(kind_phys), parameter :: A_for = 0.1 ! forest
  real(kind_phys), parameter :: A_grs = 0.2 ! grass
  real(kind_phys), parameter :: A_wat = 100. ! water
  real(kind_phys), parameter :: eps0_for = 0.8*0.01 ! forest
  real(kind_phys), parameter :: eps0_grs = 0.4*0.01 ! grass
  real(kind_phys), parameter :: eps0_wat = 0.6*0.01 ! water

  REAL(kind_phys), PARAMETER :: one3=1.0/3.0
  REAL(kind_phys), PARAMETER :: two3=2.0/3.0
!  SQRT( 2 )
  REAL(kind_phys), PARAMETER :: sqrt2=1.4142135623731
!  SQRT( PI )
  REAL(kind_phys), PARAMETER :: sqrtpi=1.7724539
  REAL(kind_phys) :: karman = 0.4                             ! von Karman constant
  REAL(kind_phys), PARAMETER :: conmin= 1.E-16
  REAL(kind_phys), PARAMETER :: pirs=3.14159265358979324
  REAL(kind_phys), PARAMETER :: f6dpi=6.0/pirs
  REAL(kind_phys), PARAMETER :: f6dpim9=1.0E-9*f6dpi
  REAL(kind_phys), PARAMETER :: rhosmoke = 1.4E3
  REAL(kind_phys), PARAMETER :: rhodust  = 2.6E3
  REAL(kind_phys), PARAMETER :: smokefac=f6dpim9/rhosmoke
  REAL(kind_phys), PARAMETER :: dustfac=f6dpim9/rhodust
!  starting standard surface temperature [ K ]
  REAL(kind_phys), PARAMETER :: tss0=288.15
  REAL(kind_phys), PARAMETER :: sigma1 = 1.8
  REAL(kind_phys), PARAMETER :: mean_diameter1 = 4.e-8
  REAL(kind_phys), PARAMETER :: fact_wfa = 1.e-9*6.0/pirs*exp(4.5*log(sigma1)**2)/mean_diameter1**3
  REAL(kind_phys), PARAMETER :: sginia=2.00
!  initial sigma-G for nucleimode                 
  REAL(kind_phys), PARAMETER :: sginin=1.70
! initial sigma-G for coarse mode               
  REAL(kind_phys), PARAMETER :: sginic=2.5
!  starting standard surface pressure [ Pa ]  
  REAL(kind_phys), PARAMETER :: pss0=101325.0
! lowest particle diameter ( m )   
  REAL(kind_phys), PARAMETER :: dgmin=1.0E-09
! lowest particle density ( Kg/m**3 )
  REAL(kind_phys), PARAMETER :: densmin=1.0E03
! index for Aitken mode number                  
  INTEGER, PARAMETER :: vdnnuc=1
! index for accumulation mode number            
  INTEGER, PARAMETER :: vdnacc=2
! index for coarse mode number                  
  INTEGER, PARAMETER :: vdncor=3
! index for Aitken mode mass                    
  INTEGER, PARAMETER :: vdmnuc=4
! index for accumulation mode                   
  INTEGER, PARAMETER :: vdmacc=5
! index for fine mode mass (Aitken + accumulation)
  INTEGER, PARAMETER :: vdmfine=6
! index for coarse mode mass                    
  INTEGER, PARAMETER :: vdmcor=7
! index for Aitken mode number                  
  INTEGER, PARAMETER :: vsnnuc=1
! index for Accumulation mode number            
  INTEGER, PARAMETER :: vsnacc=2
! index for coarse mode number                  
  INTEGER, PARAMETER :: vsncor=3
! index for Aitken mode mass                     
  INTEGER, PARAMETER :: vsmnuc=4
! index for accumulation mode mass              
  INTEGER, PARAMETER :: vsmacc=5
! index for coarse mass                         
  INTEGER, PARAMETER :: vsmcor=6
! coarse mode exp( log^2( sigmag )/8 )  
! nuclei        **4                    
      REAL(kind_phys) :: esn04
! accumulation                         
      REAL(kind_phys) :: esa04
      REAL(kind_phys) :: esc04
! coarse                               
! nuclei        **5                    
      REAL(kind_phys) :: esn05
      REAL(kind_phys) :: esa05
! accumulation                         
! nuclei        **8                    
      REAL(kind_phys) :: esn08
! accumulation                         
      REAL(kind_phys) :: esa08
      REAL(kind_phys) :: esc08
! coarse                               
! nuclei        **9                    
      REAL(kind_phys) :: esn09
      REAL(kind_phys) :: esa09
! accumulation                         
! nuclei        **12                   
      REAL(kind_phys) :: esn12
! accumulation                         
      REAL(kind_phys) :: esa12
      REAL(kind_phys) :: esc12
! coarse mode                          
! nuclei        **16                   
      REAL(kind_phys) :: esn16
! accumulation                         
      REAL(kind_phys) :: esa16
      REAL(kind_phys) :: esc16
! coarse                               
! nuclei        **20                   
      REAL(kind_phys) :: esn20
! accumulation                         
      REAL(kind_phys) :: esa20
      REAL(kind_phys) :: esc20
! coarse                               
! nuclei        **25                   
      REAL(kind_phys) :: esn25
      REAL(kind_phys) :: esa25
! accumulation                         
! nuclei        **24                   
      REAL(kind_phys) :: esn24
! accumulation                         
      REAL(kind_phys) :: esa24
      REAL(kind_phys) :: esc24
! coarse                               
! nuclei        **28                   
      REAL(kind_phys) :: esn28
! accumulation                         
      REAL(kind_phys) :: esa28
      REAL(kind_phys) :: esc28
! coarse                               
! nuclei        **32                   
      REAL(kind_phys) :: esn32
! accumulation                         
      REAL(kind_phys) :: esa32
      REAL(kind_phys) :: esc32
! coarese                              
! nuclei        **36                   
      REAL(kind_phys) :: esn36
! accumulation                         
      REAL(kind_phys) :: esa36
      REAL(kind_phys) :: esc36
! coarse                               
! nuclei        **49                   
      REAL(kind_phys) :: esn49
      REAL(kind_phys) :: esa49
! accumulation                         
! nuclei        **52                   
      REAL(kind_phys) :: esn52
      REAL(kind_phys) :: esa52
! accumulation                         
! nuclei        **64                   
      REAL(kind_phys) :: esn64
! accumulation                         
      REAL(kind_phys) :: esa64
      REAL(kind_phys) :: esc64
! coarse                               
      REAL(kind_phys) :: esn100
! nuclei        **100                  
! nuclei        **(-20)                
      REAL(kind_phys) :: esnm20
! accumulation                         
      REAL(kind_phys) :: esam20
      REAL(kind_phys) :: escm20
! coarse                               
! nuclei        **(-32)                
      REAL(kind_phys) :: esnm32
! accumulation                         
      REAL(kind_phys) :: esam32
      REAL(kind_phys) :: escm32
!SAM 10/08 Gaussian quadrature constants for SOA_VBS deposition numerical
!integration
  INTEGER, PARAMETER ::  NGAUSdv= 7   ! Number of Gaussian Quadrature Points
  REAL(kind_phys) :: xxlsgn, xxlsga, xxlsgc
  REAL(kind_phys) :: Y_GQ(NGAUSdv), WGAUS(NGAUSdv)
end module dep_data_mod

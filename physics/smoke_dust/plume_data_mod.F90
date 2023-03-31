!>\file  plume_data_mod.F90
!! This file contains data for the fire plume rise module.

module plume_data_mod

  use machine , only : kind_phys

  implicit none

  ! -- FRP parameters
  integer, dimension(0:20), parameter :: &
    catb = (/                         &
      0,                              &
      2, 1, 2, 1,                     & !floresta tropical 2 and 4 / extra trop fores 1,3,5
      2, 3, 3, 3, 3,                  & !cerrado/woody savanna :6 a 9
      4, 4, 4, 4, 4, 0, 4, 0, 0, 0, 0 & !pastagem/lavouras: 10 ...
    /)

  real(kind=kind_phys), dimension(0:4), parameter :: &
    flaming = (/ &
      0.00, & ! 
      0.45, & ! % biomass burned at flaming phase : tropical forest igbp 2 and 4
      0.45, & ! % biomass burned at flaming phase : extratropical forest igbp 1 , 3 and 5
      0.75, & ! % biomass burned at flaming phase : cerrado/woody savanna igbp 6 to 9
      0.00  & ! % biomass burned at flaming phase : pastagem/lavoura: igbp 10 a 17
    /)

  real(kind=kind_phys), dimension(0:20), parameter :: &
    msize= (/  &
      0.00021, & !0near water,1Evergreen needleleaf,2EvergreenBroadleaf,!3Deciduous Needleleaf,4Deciduous Broadleaf
      0.00021, 0.00021, 0.00021, 0.00021, & !5Mixed forest,6Closed shrublands,7Open shrublands,8Woody savannas,9Savannas,
      0.00023, 0.00022, 0.00022, 0.00022, 0.00029, &!  10Grassland,11Permanent wetlands,12cropland,13'Urban and Built-Up'
      0.00029, 0.00021, 0.00026, 0.00021, 0.00026, &!14cropland/natural vegetation mosaic,15Snow and ice,16Barren or sparsely vegetated
      0.00021, 0.00021, 0.00021, 0.00021, 0.00021, 0.00021 & !17Water,18Wooded Tundra,19Mixed Tundra,20Bare Ground Tundra
   /)

  ! -- FRP buffer indices
  integer, parameter :: p_frp_hr      = 1
  integer, parameter :: p_frp_std     = 2
  integer, parameter :: num_frp_plume = 2

  ! -- plumerise parameters
  integer, parameter :: tropical_forest = 1
  integer, parameter :: boreal_forest   = 2
  integer, parameter :: savannah        = 3
  integer, parameter :: grassland       = 4
  integer, parameter :: nveg_agreg      = 4
  integer, parameter :: wind_eff        = 1

  public

end module plume_data_mod

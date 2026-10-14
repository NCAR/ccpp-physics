!> This module provides the Fortran KIND parameters for REAL variables 
!! in the SAMF scheme, synchronized with the model's machine definitions.
module mo_conv_kind
  use machine, only : kind_phys, kind_sngl_prec, kind_dbl_prec
  implicit none
  public

  ! Define standard single and double precision kinds from machine module
  integer, parameter :: sp = kind_sngl_prec  ! 4
  integer, parameter :: dp = kind_dbl_prec  ! 8

  ! Floating point working precision (conv_wp)
  ! 1. If SAMFCNV_USE_SP is defined, force 4-byte precision for SAMF.
  ! 2. Otherwise, match kind_phys to ensure B4B with the legacy model.
#ifdef SAMFCNV_USE_SP
  integer, parameter :: conv_wp = sp
#else
  integer, parameter :: conv_wp = kind_phys
#endif

end module mo_conv_kind

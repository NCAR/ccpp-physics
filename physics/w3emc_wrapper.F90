!> \file w3emc_wrapper.f90
!! Wrapper with generic interfaces for w3emc library to reduce warnings
!  Modules wrap w3emc routines:
!    - w3difdat: https://www.nco.ncep.noaa.gov/pmb/docs/libs/w3lib/w3difdat.html
!    - w3movdat: https://www.nco.ncep.noaa.gov/pmb/docs/libs/w3lib/w3movdat.html
!
!  Example of w3difdat:
!   1. w3emc module has w3difdat interface which calls w3difdat_wrapper
!   2. w3emc_wrapper has w3difdat_wrapper interface which calls true w3difdat
module w3emc_wrapper
  use iso_fortran_env, only: real32, real64

  interface w3difdat_wrapper
     module procedure :: w3difdat32
     module procedure :: w3difdat64
  end interface w3difdat_wrapper

  interface w3movdat_wrapper
     module procedure :: w3movdat32
     module procedure :: w3movdat64
  end interface w3movdat_wrapper

contains
  subroutine w3difdat32(jdat, idat, it, rinc)
    integer, intent(in) :: jdat(8), idat(8), it
    real(real32), intent(out) :: rinc(5)
    call w3difdat(jdat, idat, it, rinc)
  end subroutine w3difdat32

  subroutine w3difdat64(jdat, idat, it, rinc)
    integer, intent(in) :: jdat(8), idat(8), it
    real(real64), intent(out) :: rinc(5)
    call w3difdat(jdat, idat, it, rinc)
  end subroutine w3difdat64

  subroutine w3movdat32(rinc, idat, jdat)
    real(real32), intent(in) :: rinc(5)
    integer, intent(in) :: idat(8)
    integer, intent(out) :: jdat(8)
    call w3movdat(rinc, idat, jdat)
  end subroutine w3movdat32

  subroutine w3movdat64(rinc, idat, jdat)
    real(real64), intent(in) :: rinc(5)
    integer, intent(in) :: idat(8)
    integer, intent(out) :: jdat(8)
    call w3movdat(rinc, idat, jdat)
  end subroutine w3movdat64
end module w3emc_wrapper

! Module to be loaded 
module w3emc
  use iso_fortran_env, only: real32, real64
  use w3emc_wrapper, only: w3difdat_wrapper, w3movdat_wrapper
  implicit none

  interface w3difdat
     module procedure :: w3difdat32
     module procedure :: w3difdat64
  end interface w3difdat

  interface w3movdat
     module procedure :: w3movdat32
     module procedure :: w3movdat64
  end interface w3movdat

contains
  subroutine w3difdat32(jdat, idat, it, rinc)
    integer, intent(in) :: jdat(8), idat(8), it
    real(real32), intent(out) :: rinc(5)
    call w3difdat_wrapper(jdat, idat, it, rinc)
  end subroutine w3difdat32

  subroutine w3difdat64(jdat, idat, it, rinc)
    integer, intent(in) :: jdat(8), idat(8), it
    real(real64), intent(out) :: rinc(5)
    call w3difdat_wrapper(jdat, idat, it, rinc)
  end subroutine w3difdat64

  subroutine w3movdat32(rinc, idat, jdat)
    real(real32), intent(in) :: rinc(5)
    integer, intent(in) :: idat(8)
    integer, intent(out) :: jdat(8)
    call w3movdat_wrapper(rinc, idat, jdat)
  end subroutine w3movdat32

  subroutine w3movdat64(rinc, idat, jdat)
    real(real64), intent(in) :: rinc(5)
    integer, intent(in) :: idat(8)
    integer, intent(out) :: jdat(8)
    call w3movdat_wrapper(rinc, idat, jdat)
  end subroutine w3movdat64
end module w3emc

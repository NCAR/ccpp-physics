!>\file mersenne_twister.F90
!! This file contains the module that calculates random numbers using the
!! Mersenne twister

! Module: mersenne_twister   Modern random number generator
!\author Iredell          Org: W/NX23     date: 2005-06-14
!> This module calculates random numbers using the Mersenne twister.
!!
!!   (It has been adapted to a Fortran 90 module from open source software.
!!   The comments from the original software are given below in the remarks.)
!!   The Mersenne twister (aka MT19937) is a state-of-the-art random number
!!   generator based on Mersenne primes and originally developed in 1997 by
!!   Matsumoto and Nishimura. It has a period before repeating of 2^19937-1,
!!   which certainly should be good enough for geophysical purposes. :-)
!!   Considering the algorithm's robustness, it runs fairly speedily.
!!   (Some timing statistics are given below in the remarks.)
!!   This adaptation uses the standard Fortran 90 random number interface,
!!   which can generate an arbitrary number of random numbers at one time.
!!   The random numbers generated are uniformly distributed between 0 and 1.
!!   The module also can generate random numbers from a Gaussian distribution
!!   with mean 0 and standard deviation 1, using a Numerical Recipes algorithm.
!!   The module also can generate uniformly random integer indices.
!!   There are also thread-safe versions of the generators in this adaptation,
!!   necessitating the passing of generator states which must be kept private.
!
! Program History Log:
!   2005-06-14  Mark Iredell
!   2026-10-08  Converted from fixed form to free-form Fortran 2008
!
! Usage:
!   The module can be compiled with 4-byte reals or with 8-byte reals, but
!   4-byte integers are required. The module should be endian-independent.
!   The Fortran 90 interfaces random_seed and random_number are overloaded
!   and can be used as in the standard by adding the appropriate use statement
!     use mersenne_twister
!   In the below use cases, harvest is a real array of arbitrary size,
!   and iharvest is an integer array of arbitrary size.
!   To generate uniformly distributed random numbers between 0 and 1,
!     call random_number(harvest)
!   To generate Gaussian distributed random numbers with 0 mean and 1 sigma,
!     call random_gauss(harvest)
!   To generate uniformly distributed random integer indices between 0 and n,
!     call random_index(n,iharvest)
!   In standard "saved" mode, the random number generator can be used without
!   setting a seed. But to set a seed, only 1 non-zero integer is required, e.g.
!     call random_setseed(4357) ! set default seed
!   The full generator state can be set via the standard interface random_seed,
!   but it is recommended to use this method only to restore saved states, e.g.
!     call random_seed(size=lsave)  ! get size of generator state seed array
!     allocate isave(lsave)         ! allocate seed array
!     call random_seed(get=isave)   ! fill seed array (then maybe save to disk)
!     call random_seed(put=isave)   ! restore state (after read from disk maybe)
!   Locally kept generator states can also be saved in a seed array, e.g.
!     type(random_stat):: stat
!     call random_seed(get=isave,stat=stat)  ! fill seed array
!     call random_seed(put=isave,stat=stat)  ! restore state
!   To generate random numbers in a threaded region, the "thread-safe" mode
!   must be used where generator states of type random_state are passed, e.g.
!     type(random_stat):: stat(8)
!     do i=1,8                               ! threadable loop
!       call random_setseed(7171*i,stat(i))  ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_number(harvest,stat(i))  ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_gauss(harvest,stat(i))   ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_index(n,iharvest,stat(i))! thread-safe call
!     enddo
!   There is also a relatively inefficient "interactive" mode available, where
!   setting seeds and generating random numbers are done in the same call.
!   There is also a functional mode available, returning one value at a time.
!
! Public Defined Types:
!   random_stat       Generator state (private contents)
!
! Public Subprograms:
!   random_seed         determine size or put or get state
!     size              optional integer output size of seed array
!     put               optional integer(:) input seed array
!     get               optional integer(:) output seed array
!     stat              optional type(random_stat) (thread-safe mode)
!   random_setseed      set seed (thread-safe mode)
!     inseed            integer seed input
!     stat              type(random_stat) output
!   random_setseed      set seed (saved mode)
!     inseed            integer seed input
!   random_number       get mersenne twister random numbers (thread-safe mode)
!     harvest           real(:) numbers output
!     stat              type(random_stat) input
!   random_number       get mersenne twister random numbers (saved mode)
!     harvest           real(:) numbers output
!   random_number       get mersenne twister random numbers (interactive mode)
!     harvest           real(:) numbers output
!     inseed            integer seed input
!   random_number_f     get mersenne twister random number (functional mode)
!     harvest           real number output
!   random_gauss        get gaussian random numbers (thread-safe mode)
!     harvest           real(:) numbers output
!     stat              type(random_stat) input
!   random_gauss        get gaussian random numbers (saved mode)
!     harvest           real(:) numbers output
!   random_gauss        get gaussian random numbers (interactive mode)
!     harvest           real(:) numbers output
!     inseed            integer seed input
!   random_gauss_f      get gaussian random number (functional mode)
!     harvest           real number output
!   random_index        get random indices (thread-safe mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!     stat              type(random_stat) input
!   random_index        get random indices (saved mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!   random_index        get random indices (interactive mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!     inseed            integer seed input
!   random_index_f      get random index (functional mode)
!     imax              integer maximum index input
!     iharvest          integer number output
!
! Remarks:
!   (1) Here are the comments in the original open source code:
!     A C-program for MT19937: Real number version
!     genrand() generates one pseudorandom real number (double)
!     which is uniformly distributed on [0,1]-interval, for each
!     call. sgenrand(seed) set initial values to the working area
!     of 624 words. Before genrand(), sgenrand(seed) must be
!     called once. (seed is any 32-bit integer except for 0).
!     Integer generator is obtained by modifying two lines.
!     Coded by Takuji Nishimura, considering the suggestions by
!     Topher Cooper and Marc Rieffel in July-Aug. 1997.
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Library General Public
!     License as published by the Free Software Foundation; either
!     version 2 of the License, or (at your option) any later
!     version.
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!     See the GNU Library General Public License for more details.
!     You should have received a copy of the GNU Library General
!     Public License along with this library; if not, write to the
!     Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
!     02111-1307  USA
!     Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
!     When you use this, send an email to: matumoto@math.keio.ac.jp
!     with an appropriate reference to your work.
!     Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   (2) On a single IBM Power4 processor on the NCEP operational cluster (2005)
!     each Mersenne twister random number takes less than 30 ns, about 3 times
!     slower than the default random number generator, and each random number
!     from a Gaussian distribution takes less than 150 ns.
!
! Attributes:
!   Language: Fortran 2008 (free form)
!
!$$$
module mersenne_twister
  use machine, only: kind_dbl_prec
  implicit none
  private

  !  Public declarations
  public :: random_stat
  public :: random_seed
  public :: random_setseed
  public :: random_number
  public :: random_number_f
  public :: random_gauss
  public :: random_gauss_f
  public :: random_index
  public :: random_index_f

  !  Parameters
  integer, parameter :: n = 624
  integer, parameter :: m = 397
  integer, parameter :: mata   = -1727483681     !< constant vector a
  ! Written as -huge-1 because the literal 2147483648 is out of range for a
  ! default (4-byte) integer; the value is unchanged (-2147483648).
  integer, parameter :: umask  = -huge(1) - 1    !< most significant w-r bits
  integer, parameter :: lmask  = 2147483647      !< least significant r bits
  integer, parameter :: tmaskb = -1658038656     !< tempering parameter
  integer, parameter :: tmaskc = -272236544      !< tempering parameter
  integer, parameter :: mag01(0:1) = [0, mata]
  integer, parameter :: iseed = 4357
  integer, parameter :: nrest = n + 4
  real(kind_dbl_prec), parameter :: two32 = 2.0_kind_dbl_prec**32

  !  Defined types
  type :: random_stat                            !< Generator state
    private
    integer :: mti = n + 1
    integer :: mt(0:n-1)
    integer :: iset
    real(kind_dbl_prec) :: gset
  end type random_stat

  !  Saved data
  type(random_stat), save :: sstat

  !  Overloaded interfaces
  interface random_setseed
    module procedure random_setseed_s
    module procedure random_setseed_t
  end interface random_setseed

  interface random_number
    module procedure random_number_i
    module procedure random_number_s
    module procedure random_number_t
  end interface random_number

  interface random_gauss
    module procedure random_gauss_i
    module procedure random_gauss_s
    module procedure random_gauss_t
  end interface random_gauss

  interface random_index
    module procedure random_index_i
    module procedure random_index_s
    module procedure random_index_t
  end interface random_index

!  All the subprograms
contains

  !  Subprogram random_seed
  !> This subroutine sets and gets state; overloads Fortran 90 standard.
  subroutine random_seed(size, put, get, stat)
    integer, intent(out), optional :: size
    integer, intent(in),  optional :: put(nrest)
    integer, intent(out), optional :: get(nrest)
    type(random_stat), intent(inout), optional :: stat

    if (present(size)) then          ! return size of seed array
      ! if (present(put) .or. present(get)) &
      !   call errmsg('RANDOM_SEED: more than one option set - some ignored')
      size = nrest
    else if (present(put)) then      ! restore from seed array
      ! if (present(get)) &
      !   call errmsg('RANDOM_SEED: more than one option set - some ignored')
      if (present(stat)) then
        call put_state(stat)
      else
        call put_state(sstat)
      end if
    else if (present(get)) then      ! save to seed array
      if (present(stat)) then
        call get_state(stat)
      else
        call get_state(sstat)
      end if
    else                             ! reset default seed
      if (present(stat)) then
        call random_setseed_t(iseed, stat)
      else
        call random_setseed_t(iseed, sstat)
      end if
    end if

  contains

    !> Restore a generator state from the seed array put.
    subroutine put_state(s)
      type(random_stat), intent(inout) :: s
      s%mti  = put(1)
      s%mt   = put(2:n+1)
      s%iset = put(n+2)
      s%gset = transfer(put(n+3:nrest), s%gset)
      if (s%mti < 0 .or. s%mti > n .or. any(s%mt == 0) .or. &
          s%iset < 0 .or. s%iset > 1) then
        call random_setseed_t(iseed, s)
        ! call errmsg('RANDOM_SEED: invalid seeds put - default seeds used')
      end if
    end subroutine put_state

    !> Save a generator state to the seed array get.
    subroutine get_state(s)
      type(random_stat), intent(inout) :: s
      if (s%mti == n + 1) call random_setseed_t(iseed, s)
      get(1)     = s%mti
      get(2:n+1) = s%mt
      get(n+2)   = s%iset
      get(n+3:nrest) = transfer(s%gset, get, nrest - (n + 3) + 1)
    end subroutine get_state

  end subroutine random_seed

  !  Subprogram random_setseed_s
  !> This subroutine sets seed in saved mode.
  subroutine random_setseed_s(inseed)
    integer, intent(in) :: inseed
    call random_setseed_t(inseed, sstat)
  end subroutine random_setseed_s

  !  Subprogram random_setseed_t
  !> This subroutine sets seed in thread-safe mode.
  subroutine random_setseed_t(inseed, stat)
    integer, intent(in) :: inseed
    type(random_stat), intent(out) :: stat
    integer :: ii, mti

    ii = inseed
    if (ii == 0) ii = iseed
    stat%mti = n
    stat%mt(0) = iand(ii, -1)
    ! NOTE: this recurrence relies on 32-bit integer multiplication wrapping
    ! around on overflow, as in the original code. That is not guaranteed by
    ! the standard; do not compile this file with overflow trapping (-ftrapv).
    do mti = 1, n - 1
      stat%mt(mti) = iand(69069*stat%mt(mti-1), -1)
    end do
    stat%iset = 0
    stat%gset = 0.0_kind_dbl_prec
  end subroutine random_setseed_t

  !  Subprogram random_number_f
  !> This function generates random numbers in functional mode.
  function random_number_f() result(harvest)
    real(kind_dbl_prec) :: harvest
    real(kind_dbl_prec) :: h(1)

    if (sstat%mti == n + 1) call random_setseed_t(iseed, sstat)
    call random_number_t(h, sstat)
    harvest = h(1)
  end function random_number_f

  !  Subprogram random_number_i
  !> This subroutine generates random numbers in interactive mode.
  subroutine random_number_i(harvest, inseed)
    real(kind_dbl_prec), intent(out) :: harvest(:)
    integer, intent(in) :: inseed
    type(random_stat) :: stat

    call random_setseed_t(inseed, stat)
    call random_number_t(harvest, stat)
  end subroutine random_number_i

  !  Subprogram random_number_s
  !> This subroutine generates random numbers in saved mode; overloads Fortran 90 standard.
  subroutine random_number_s(harvest)
    real(kind_dbl_prec), intent(out) :: harvest(:)

    if (sstat%mti == n + 1) call random_setseed_t(iseed, sstat)
    call random_number_t(harvest, sstat)
  end subroutine random_number_s

  !  Subprogram random_number_t
  !> This subroutine generates random numbers in thread-safe mode.
  subroutine random_number_t(harvest, stat)
    real(kind_dbl_prec), intent(out) :: harvest(:)
    type(random_stat), intent(inout) :: stat
    integer :: j, kk, y

    do j = 1, size(harvest)
      if (stat%mti >= n) then        ! regenerate the full state vector
        do kk = 0, n - m - 1
          y = ior(iand(stat%mt(kk), umask), iand(stat%mt(kk+1), lmask))
          stat%mt(kk) = ieor(ieor(stat%mt(kk+m), shiftr(y, 1)), &
                             mag01(iand(y, 1)))
        end do
        do kk = n - m, n - 2
          y = ior(iand(stat%mt(kk), umask), iand(stat%mt(kk+1), lmask))
          stat%mt(kk) = ieor(ieor(stat%mt(kk+(m-n)), shiftr(y, 1)), &
                             mag01(iand(y, 1)))
        end do
        y = ior(iand(stat%mt(n-1), umask), iand(stat%mt(0), lmask))
        stat%mt(n-1) = ieor(ieor(stat%mt(m-1), shiftr(y, 1)), &
                            mag01(iand(y, 1)))
        stat%mti = 0
      end if

      ! Tempering (formerly the statement functions tshftu/s/t/l)
      y = stat%mt(stat%mti)
      y = ieor(y, shiftr(y, 11))
      y = ieor(y, iand(shiftl(y, 7), tmaskb))
      y = ieor(y, iand(shiftl(y, 15), tmaskc))
      y = ieor(y, shiftr(y, 18))

      if (y < 0) then
        harvest(j) = (real(y, kind=kind_dbl_prec) + two32) / &
                     (two32 - 1.0_kind_dbl_prec)
      else
        ! NOTE: real(y) converts via DEFAULT real, exactly as in the original
        ! code, so this branch is rounded to single precision unless default
        ! reals are promoted (e.g. -fdefault-real-8 / -r8). Kept as is so
        ! results stay bit-for-bit identical; real(y, kind=kind_dbl_prec)
        ! would be the consistent form but changes the random sequence.
        harvest(j) = real(y) / (two32 - 1.0_kind_dbl_prec)
      end if
      stat%mti = stat%mti + 1
    end do
  end subroutine random_number_t

  !  Subprogram random_gauss_f
  !> This function generates Gaussian random numbers in functional mode.
  function random_gauss_f() result(harvest)
    real(kind_dbl_prec) :: harvest
    real(kind_dbl_prec) :: h(1)

    if (sstat%mti == n + 1) call random_setseed_t(iseed, sstat)
    call random_gauss_t(h, sstat)
    harvest = h(1)
  end function random_gauss_f

  !  Subprogram random_gauss_i
  !> This subroutine generates Gaussian random numbers in interactive mode.
  subroutine random_gauss_i(harvest, inseed)
    real(kind_dbl_prec), intent(out) :: harvest(:)
    integer, intent(in) :: inseed
    type(random_stat) :: stat

    call random_setseed_t(inseed, stat)
    call random_gauss_t(harvest, stat)
  end subroutine random_gauss_i

  !  Subprogram random_gauss_s
  !> This subroutine generates Gaussian random numbers in saved mode.
  subroutine random_gauss_s(harvest)
    real(kind_dbl_prec), intent(out) :: harvest(:)

    if (sstat%mti == n + 1) call random_setseed_t(iseed, sstat)
    call random_gauss_t(harvest, sstat)
  end subroutine random_gauss_s

  !  Subprogram random_gauss_t
  !> This subroutine generates Gaussian random numbers in thread-safe mode.
  subroutine random_gauss_t(harvest, stat)
    real(kind_dbl_prec), intent(out) :: harvest(:)
    type(random_stat), intent(inout) :: stat
    integer :: mx, my, mz, j
    real(kind_dbl_prec) :: r2(2), r, g1, g2

    mz = size(harvest)
    if (mz <= 0) return
    mx = 0
    if (stat%iset == 1) then
      mx = 1
      harvest(1) = stat%gset
      stat%iset = 0
    end if
    my = (mz - mx)/2*2 + mx
    do
      call random_number_t(harvest(mx+1:my), stat)
      do j = mx, my - 2, 2
        call rgauss(harvest(j+1), harvest(j+2), r, g1, g2)
        if (r < 1.0_kind_dbl_prec) then
          harvest(mx+1) = g1
          harvest(mx+2) = g2
          mx = mx + 2
        end if
      end do
      if (mx == my) exit
    end do
    if (my < mz) then
      do
        call random_number_t(r2, stat)
        call rgauss(r2(1), r2(2), r, g1, g2)
        if (r < 1.0_kind_dbl_prec) exit
      end do
      harvest(mz) = g1
      stat%gset = g2
      stat%iset = 1
    end if

  contains

    !> This subroutine contains numerical Recipes algorithm to generate Gaussian random numbers.
    !! g1 and g2 are only defined on return when r < 1.
    subroutine rgauss(r1, r2, r, g1, g2)
      real(kind_dbl_prec), intent(in)  :: r1, r2
      real(kind_dbl_prec), intent(out) :: r, g1, g2
      real(kind_dbl_prec) :: v1, v2, fac

      v1 = 2._kind_dbl_prec*r1 - 1._kind_dbl_prec
      v2 = 2._kind_dbl_prec*r2 - 1._kind_dbl_prec
      r = v1**2 + v2**2
      if (r < 1.0_kind_dbl_prec) then
        fac = sqrt(-2._kind_dbl_prec*log(r)/r)
        g1 = v1*fac
        g2 = v2*fac
      end if
    end subroutine rgauss

  end subroutine random_gauss_t

  !  Subprogram random_index_f
  !> This function generates random indices in functional mode.
  function random_index_f(imax) result(iharvest)
    integer, intent(in) :: imax
    integer :: iharvest
    integer :: ih(1)

    if (sstat%mti == n + 1) call random_setseed_t(iseed, sstat)
    call random_index_t(imax, ih, sstat)
    iharvest = ih(1)
  end function random_index_f

  !  Subprogram random_index_i
  !> This subroutine generates random indices in interactive mode.
  subroutine random_index_i(imax, iharvest, inseed)
    integer, intent(in)  :: imax
    integer, intent(out) :: iharvest(:)
    integer, intent(in)  :: inseed
    type(random_stat) :: stat

    call random_setseed_t(inseed, stat)
    call random_index_t(imax, iharvest, stat)
  end subroutine random_index_i

  !  Subprogram random_index_s
  !> This subroutine generates random indices in saved mode.
  subroutine random_index_s(imax, iharvest)
    integer, intent(in)  :: imax
    integer, intent(out) :: iharvest(:)

    if (sstat%mti == n + 1) call random_setseed_t(iseed, sstat)
    call random_index_t(imax, iharvest, sstat)
  end subroutine random_index_s

  !  Subprogram random_index_t
  !> This subroutine generates random indices in thread-safe mode.
  subroutine random_index_t(imax, iharvest, stat)
    integer, intent(in)  :: imax
    integer, intent(out) :: iharvest(:)
    type(random_stat), intent(inout) :: stat
    integer, parameter :: mh = n
    integer :: i1, i2, mz
    real(kind_dbl_prec) :: h(mh)

    mz = size(iharvest)
    do i1 = 1, mz, mh
      i2 = min((i1 - 1) + mh, mz)
      call random_number_t(h(:i2-(i1-1)), stat)
      iharvest(i1:i2) = max(ceiling(h(:i2-(i1-1))*imax), 1)
    end do
  end subroutine random_index_t

end module mersenne_twister

module mpiutil

  use iso_fortran_env, only : real32, real64
  use iso_fortran_env, only : error_unit, output_unit
  use mpi_f08

  implicit none

  private
  public ccpp_bcast

  interface ccpp_bcast
     procedure :: bcast_i32d0
     procedure :: bcast_i32d1
     procedure :: bcast_i32d2
     procedure :: bcast_i32d3
     procedure :: bcast_r32d0
     procedure :: bcast_r64d0
     procedure :: bcast_r32d1
     procedure :: bcast_r64d1
     procedure :: bcast_r32d2
     procedure :: bcast_r64d2
     procedure :: bcast_r32d3
     procedure :: bcast_r64d3
     procedure :: bcast_r32d4
     procedure :: bcast_r64d4
     procedure :: bcast_r32d5
     procedure :: bcast_r64d5
     procedure :: bcast_ld0
  end interface ccpp_bcast

contains

! Helper routines for MPI broadcasting

   subroutine bcast_i32d0(arr, root, comm, ierr)
      integer, intent(inout) :: arr
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, 1, MPI_INTEGER, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_i32d0")
      end if
   end subroutine bcast_i32d0

   subroutine bcast_i32d1(arr, root, comm, ierr)
      integer, intent(inout) :: arr(:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_INTEGER, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_i32d1")
      end if
   end subroutine bcast_i32d1

   subroutine bcast_i32d2(arr, root, comm, ierr)
      integer, intent(inout) :: arr(:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_INTEGER, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_i32d2")
      end if
   end subroutine bcast_i32d2

   subroutine bcast_i32d3(arr, root, comm, ierr)
      integer, intent(inout) :: arr(:,:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_INTEGER, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_i32d3")
      end if
   end subroutine bcast_i32d3

   subroutine bcast_r32d0(arr, root, comm, ierr)
      real(kind=real32), intent(inout) :: arr
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, 1, MPI_REAL, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r32d0")
      end if
   end subroutine bcast_r32d0

   subroutine bcast_r64d0(arr, root, comm, ierr)
      real(kind=real64), intent(inout) :: arr
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, 1, MPI_DOUBLE_PRECISION, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r64d0")
      end if
   end subroutine bcast_r64d0

   subroutine bcast_r32d1(arr, root, comm, ierr)
      real(kind=real32), intent(inout) :: arr(:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_REAL, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r32d1")
      end if
   end subroutine bcast_r32d1

   subroutine bcast_r64d1(arr, root, comm, ierr)
      real(kind=real64), intent(inout) :: arr(:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_DOUBLE_PRECISION, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r64d1")
      end if
   end subroutine bcast_r64d1

   subroutine bcast_r32d2(arr, root, comm, ierr)
      real(kind=real32), intent(inout) :: arr(:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_REAL, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r32d2")
      end if
   end subroutine bcast_r32d2

   subroutine bcast_r64d2(arr, root, comm, ierr)
      real(kind=real64), intent(inout) :: arr(:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_DOUBLE_PRECISION, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r64d2")
      end if
   end subroutine bcast_r64d2

   subroutine bcast_r32d3(arr, root, comm, ierr)
      real(kind=real32), intent(inout) :: arr(:,:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_REAL, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r32d3")
      end if
   end subroutine bcast_r32d3

   subroutine bcast_r64d3(arr, root, comm, ierr)
      real(kind=real64), intent(inout) :: arr(:,:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_DOUBLE_PRECISION, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r64d3")
      end if
   end subroutine bcast_r64d3

   subroutine bcast_r32d4(arr, root, comm, ierr)
      real(kind=real32), intent(inout) :: arr(:,:,:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_REAL, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r32d4")
      end if
   end subroutine bcast_r32d4

   subroutine bcast_r64d4(arr, root, comm, ierr)
      real(kind=real64), intent(inout) :: arr(:,:,:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_DOUBLE_PRECISION, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r64d4")
      end if
   end subroutine bcast_r64d4

   subroutine bcast_r32d5(arr, root, comm, ierr)
      real(kind=real32), intent(inout) :: arr(:,:,:,:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_REAL, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r32d5")
      end if
   end subroutine bcast_r32d5

   subroutine bcast_r64d5(arr, root, comm, ierr)
      real(kind=real64), intent(inout) :: arr(:,:,:,:,:)
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, size(arr), MPI_DOUBLE_PRECISION, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_r64d5")
      end if
   end subroutine bcast_r64d5

   subroutine bcast_ld0(arr, root, comm, ierr)
      logical, intent(inout) :: arr
      integer, intent(in) :: root
      type(MPI_Comm), intent(in) :: comm
      integer, intent(out) :: ierr
      call MPI_BCAST(arr, 1, MPI_LOGICAL, root, comm, ierr)
      if (ierr/=MPI_SUCCESS) then
         call ccpp_abort(comm, "mpiutil.F90:bcast_ld0")
      end if
   end subroutine bcast_ld0

! Temporary: helper routine to abort code until
! discussion about potential redesign of how to
! abort model runs for CCPP errors is settled.

   subroutine ccpp_abort(comm, str)
#ifdef __INTEL_COMPILER
      use ifcore
#endif
      implicit none
      type(MPI_Comm), intent(in) :: comm
      character(len=*), intent(in) :: str
      integer :: ierr
      write(output_unit,'(a)') "ccpp_abort: " // trim(str)
      write(error_unit,'(a)') "ccpp_abort: " // trim(str)
#if defined(__INTEL_COMPILER)
      call tracebackqq("ccpp_abort" // trim(str), user_exit_code=-1)
#else if defined(__GFORTRAN__)
      call backtrace()
#endif
      call MPI_ABORT(comm, ierr)
   end subroutine ccpp_abort

end module mpiutil

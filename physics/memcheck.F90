!> \file memcheck.F90
!!  Contains code to check memory usage with/without CCPP.

    module memcheck

      use machine, only: kind_phys

      implicit none

      private
 
      public memcheck_init, memcheck_run, memcheck_finalize

      ! Can use larger time frame to track memory leaks
      real(kind_phys), parameter :: SECONDS_ELAPSED_MIN = 3500.0
      real(kind_phys), parameter :: SECONDS_ELAPSED_MAX = 3700.0

      contains

      subroutine memcheck_init ()
      end subroutine memcheck_init

      subroutine memcheck_finalize ()
      end subroutine memcheck_finalize

!> \section arg_table_memcheck_run Argument Table
!! | local_name      | standard_name                                          | long_name                                               | units         | rank | type      |    kind   | intent | optional |
!! |-----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | seconds_elapsed | seconds_elapsed_since_model_initialization             | seconds elapsed since model initialization              | s             |    0 | real      | kind_phys | in     | F        |
!! | block_number    | block_number                                           | for explicit data blocking: block number of this block  | index         |    0 | integer   |           | in     | F        |
!! | mpicomm         | mpi_comm                                               | MPI communicator                                        | index         |    0 | integer   |           | in     | F        |
!! | errmsg          | error_message                                          | error message for error handling in CCPP                | none          |    0 | character | len=*     | out    | F        |
!! | errflg          | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer   |           | out    | F        |
!!
      subroutine memcheck_run (seconds_elapsed, block_number, mpicomm, errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use ccpp_api, only: ccpp_memory_usage

         implicit none

         !--- interface variables
         real(kind=kind_phys),       intent(in)  :: seconds_elapsed
         integer,                    intent(in)  :: block_number
         integer,                    intent(in)  :: mpicomm
         character(len=*),           intent(out) :: errmsg
         integer,                    intent(out) :: errflg

         !--- local variables
         integer :: impi, ierr
         integer :: mpirank, mpisize
         integer :: ompthread
         character(len=1024) :: memory_usage

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (seconds_elapsed < SECONDS_ELAPSED_MIN .or. &
            seconds_elapsed > SECONDS_ELAPSED_MAX) return

         if (block_number>1) return

#ifdef MPI
         call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr)
#else
         mpirank = 0
         mpisize = 1
#endif

#ifdef OPENMP
         ompthread = OMP_GET_THREAD_NUM()
#else
         ompthread = 0
#endif
 
         ierr = ccpp_memory_usage(mpicomm, memory_usage)

         ! Output ordered by MPI rank
         do impi=0,mpisize-1
            if (mpirank==impi .and. ompthread==0) then
                write(0,'(a)') trim(memory_usage)
            end if
#ifdef MPI
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
         end do

      end subroutine memcheck_run

    end module memcheck

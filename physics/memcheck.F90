!> \file memcheck.F90
!!  Contains code to check memory usage with/without CCPP.

    module memcheck

      implicit none

      private
 
      public memcheck_init, memcheck_run, memcheck_finalize

   contains

      subroutine memcheck_init ()
      end subroutine memcheck_init

      subroutine memcheck_finalize ()
      end subroutine memcheck_finalize

!> \section arg_table_memcheck_run Argument Table
!! | local_name      | standard_name                                          | long_name                                               | units | rank | type      |    kind   | intent | optional |
!! |-----------------|--------------------------------------------------------|---------------------------------------------------------|-------|------|-----------|-----------|--------|----------|
!! | mpicomm         | mpi_comm                                               | MPI communicator                                        | index |    0 | integer   |           | in     | F        |
!! | mpirank         | mpi_rank                                               | current MPI-rank                                        | index |    0 | integer   |           | in     | F        |
!! | mpisize         | mpi_size                                               | number of MPI tasks in communicator                     | count |    0 | integer   |           | in     | F        |
!! | mpiroot         | mpi_root                                               | master MPI-rank                                         | index |    0 | integer   |           | in     | T        |
!! | ompthreads      | omp_threads                                            | number of OpenMP threads available for physics schemes  | count |    0 | integer   |           | in     | F        |
!! | errmsg          | ccpp_error_message                                     | error message for error handling in CCPP                | none  |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag  |    0 | integer   |           | out    | F        |
!!
      subroutine memcheck_run (mpicomm, mpirank, mpisize, mpiroot, ompthreads, errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use ccpp_api, only: ccpp_memory_usage

         implicit none

         !--- interface variables
         integer,                    intent(in)  :: mpicomm
         integer,                    intent(in)  :: mpirank
         integer,                    intent(in)  :: mpisize
         integer, optional,          intent(in)  :: mpiroot
         integer,                    intent(in)  :: ompthreads
         character(len=*),           intent(out) :: errmsg
         integer,                    intent(out) :: errflg

         !--- local variables
         integer :: impi, ierr
         integer :: ompthread
         character(len=1024) :: memory_usage

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef OPENMP
         ompthread = OMP_GET_THREAD_NUM()
#else
         ompthread = 0
#endif

         ierr = ccpp_memory_usage(mpicomm, memory_usage)
         if (present(mpiroot) .and. mpirank==mpiroot) then
            write(0,'(a)') trim(memory_usage)
         else if (.not.present(mpiroot)) then
            ! Output ordered by MPI rank
            do impi=0,mpisize-1
               if (mpirank==impi) then
                   write(0,'(a)') trim(memory_usage)
               end if
#ifdef MPI
               call MPI_BARRIER(mpicomm,ierr)
#endif
            end do
         end if

#ifdef MPI
         call MPI_BARRIER(mpicomm,ierr)
#endif

      end subroutine memcheck_run

    end module memcheck

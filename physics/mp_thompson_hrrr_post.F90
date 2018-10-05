module mp_thompson_hrrr_post

   use machine, only : kind_phys

   implicit none

   public :: mp_thompson_hrrr_post_init, mp_thompson_hrrr_post_run, mp_thompson_hrrr_post_finalize

   private

   logical :: is_initialized = .False.

   real(kind_phys), dimension(:), allocatable :: mp_tend_lim

contains

#if 0
!! \section arg_table_mp_thompson_hrrr_post_init Argument Table
!! | local_name      | standard_name                                         | long_name                                                | units      | rank | type      |    kind   | intent | optional |
!! |-----------------|-------------------------------------------------------|----------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ncol            | horizontal_loop_extent                                | horizontal loop extent                                   | count      |    0 | integer   |           | in     | F        |
!! | area            | cell_area                                             | area of the grid cell                                    | m2         |    1 | real      | kind_phys | in     | F        |
!! | errmsg          | ccpp_error_message                                    | error message for error handling in CCPP                 | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                       | error flag for error handling in CCPP                    | flag       |    0 | integer   |           | out    | F        |
!!
#endif
   subroutine mp_thompson_hrrr_post_init(ncol, area, errmsg, errflg)

      implicit none

      ! Interface variables
      integer,                            intent(in) :: ncol
      real(kind_phys), dimension(1:ncol), intent(in) :: area

      ! CCPP error handling
      character(len=*), intent(  out) :: errmsg
      integer,          intent(  out) :: errflg

      ! Local variables
      real(kind_phys), dimension(1:ncol) :: dx
      integer :: i

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (is_initialized) return

      allocate(mp_tend_lim(1:ncol))

      ! Cell size in m as square root of cell area
      dx = sqrt(area)

      do i=1,ncol
         ! DH* testing/development: limiters on temperature tendency
         ! depending on the grid spacing - no limit for >50km *DH
         if (dx(i)<=5000) then
            mp_tend_lim(i) = 0.07    ! [K/s], 3-km HRRR value
         else if (dx(i)<=50000) then
            mp_tend_lim(i) = 0.002   ! [K/s], 13-km RAP value
         else
            ! no limit for grid spacings >50km
            !mp_tend_lim(i) = 0.00006 ! [K/s], guess for >50km
            mp_tend_lim(i) = huge(mp_tend_lim(i))
         end if
      end do

      is_initialized = .True.

   end subroutine mp_thompson_hrrr_post_init

#if 0
!! \section arg_table_mp_thompson_hrrr_post_run Argument Table
!! | local_name      | standard_name                                         | long_name                                                | units      | rank | type      |    kind   | intent | optional |
!! |-----------------|-------------------------------------------------------|----------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ncol            | horizontal_loop_extent                                | horizontal loop extent                                   | count      |    0 | integer   |           | in     | F        |
!! | nlev            | vertical_dimension                                    | number of vertical levels                                | count      |    0 | integer   |           | in     | F        |
!! | tgrs_save       | air_temperature_save                                  | air temperature before entering a physics scheme         | K          |    2 | real      | kind_phys | in     | F        |
!! | tgrs            | air_temperature_updated_by_physics                    | model layer mean temperature                             | K          |    2 | real      | kind_phys | inout  | F        |
!! | prslk           | dimensionless_exner_function_at_model_layers          | dimensionless Exner function at model layer centers      | none       |    2 | real      | kind_phys | in     | F        |
!! | dtp             | time_step_for_physics                                 | physics timestep                                         | s          |    0 | real      | kind_phys | in     | F        |
!! | mpicomm         | mpi_comm                                              | MPI communicator                                         | index      |    0 | integer   |           | in     | F        |
!! | mpirank         | mpi_rank                                              | current MPI-rank                                         | index      |    0 | integer   |           | in     | F        |
!! | mpiroot         | mpi_root                                              | master MPI-rank                                          | index      |    0 | integer   |           | in     | F        |
!! | errmsg          | ccpp_error_message                                    | error message for error handling in CCPP                 | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                       | error flag for error handling in CCPP                    | flag       |    0 | integer   |           | out    | F        |
!!
#endif
   subroutine mp_thompson_hrrr_post_run(ncol, nlev, tgrs_save, tgrs, prslk, dtp, &
                                        mpicomm, mpirank, mpiroot, errmsg, errflg)

      implicit none

      ! Interface variables
      integer,                                   intent(in)    :: ncol
      integer,                                   intent(in)    :: nlev
      real(kind_phys), dimension(1:ncol,1:nlev), intent(in)    :: tgrs_save
      real(kind_phys), dimension(1:ncol,1:nlev), intent(inout) :: tgrs
      real(kind_phys), dimension(1:ncol,1:nlev), intent(in)    :: prslk
      real(kind_phys),                           intent(in)    :: dtp
      ! MPI information
      integer,          intent(in   ) :: mpicomm
      integer,          intent(in   ) :: mpirank
      integer,          intent(in   ) :: mpiroot
      ! CCPP error handling
      character(len=*), intent(  out) :: errmsg
      integer,          intent(  out) :: errflg

      ! Local variables
      real(kind_phys), dimension(1:ncol,1:nlev) :: mp_tend
      integer :: i, k

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (.not.is_initialized) then
         write(errmsg, fmt='((a))') 'mp_thompson_hrrr_post_run called before mp_thompson_hrrr_post_init'
         errflg = 1
         return
      end if

      ! mp_tend and mp_tend_lim are expressed in potential temperature
      mp_tend = (tgrs - tgrs_save)/prslk

      do k=1,nlev
         do i=1,ncol
            mp_tend(i,k) = max( -mp_tend_lim(i)*dtp, min( mp_tend_lim(i)*dtp, mp_tend(i,k) ) )
            ! DH*
            if ( mpirank==mpiroot .and. (tgrs_save(i,k) + mp_tend(i,k)*prslk(i,k) .ne. tgrs(i,k)) ) then
               write(0,*) "DH DEBUG mp_thompson_hrrr_post_run mp_tend limiter: i, k, t_old, t_new, t_lim:", &
                                  & i, k, tgrs_save(i,k), tgrs(i,k), tgrs_save(i,k) + mp_tend(i,k)*prslk(i,k)
            end if
            ! *DH
            tgrs(i,k) = tgrs_save(i,k) + mp_tend(i,k)*prslk(i,k)
         end do
      end do

   end subroutine mp_thompson_hrrr_post_run

#if 0
!! \section arg_table_mp_thompson_hrrr_post_finalize Argument Table
!! | local_name      | standard_name                                         | long_name                                                | units      | rank | type      |    kind   | intent | optional |
!! |-----------------|-------------------------------------------------------|----------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | errmsg          | ccpp_error_message                                    | error message for error handling in CCPP                 | none       |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                       | error flag for error handling in CCPP                    | flag       |    0 | integer   |           | out    | F        |
!!
#endif
   subroutine mp_thompson_hrrr_post_finalize(errmsg, errflg)

      implicit none

      ! CCPP error handling
      character(len=*),          intent(  out) :: errmsg
      integer,                   intent(  out) :: errflg

      ! Check initialization state
      if (.not. is_initialized) return

      if (allocated(mp_tend_lim)) deallocate(mp_tend_lim)

      is_initialized = .False.

   end subroutine mp_thompson_hrrr_post_finalize

end module mp_thompson_hrrr_post

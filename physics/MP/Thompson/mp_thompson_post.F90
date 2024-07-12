module mp_thompson_post

   use mpi_f08
   use machine, only : kind_phys

   implicit none

   public :: mp_thompson_post_init, mp_thompson_post_run, mp_thompson_post_finalize

   private

   logical :: is_initialized = .false.

   logical :: apply_limiter

contains

!! \section arg_table_mp_thompson_post_init Argument Table
!! \htmlinclude mp_thompson_post_init.html
!!
   subroutine mp_thompson_post_init(ttendlim, errmsg, errflg)

      implicit none

      ! Interface variables
      real(kind_phys), intent(in) :: ttendlim

      ! CCPP error handling
      character(len=*), intent(  out) :: errmsg
      integer,          intent(  out) :: errflg

      ! Local variables
      integer :: i

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (is_initialized) return

      if (ttendlim < 0) then
          apply_limiter = .false.
      else
          apply_limiter = .true.
      end if

      is_initialized = .true.

   end subroutine mp_thompson_post_init

!> \section arg_table_mp_thompson_post_run Argument Table
!! \htmlinclude mp_thompson_post_run.html
!!
   subroutine mp_thompson_post_run(ncol, nlev, tgrs_save, tgrs, prslk, dtp, ttendlim, &
                                   kdt, mpicomm, mpirank, mpiroot, errmsg, errflg)

      implicit none

      ! Interface variables
      integer,                         intent(in)    :: ncol
      integer,                         intent(in)    :: nlev
      real(kind_phys), dimension(:,:), intent(in)    :: tgrs_save
      real(kind_phys), dimension(:,:), intent(inout) :: tgrs
      real(kind_phys), dimension(:,:), intent(in)    :: prslk
      real(kind_phys),                 intent(in)    :: dtp
      real(kind_phys),                 intent(in)    :: ttendlim
      integer,                         intent(in)    :: kdt
      ! MPI information
      type(MPI_Comm),   intent(in   ) :: mpicomm
      integer,          intent(in   ) :: mpirank
      integer,          intent(in   ) :: mpiroot
      ! CCPP error handling
      character(len=*), intent(  out) :: errmsg
      integer,          intent(  out) :: errflg

      ! Local variables
      real(kind_phys), dimension(1:ncol,1:nlev) :: mp_tend
      integer :: i, k
#ifdef DEBUG
      integer :: events
#endif

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (.not.is_initialized) then
         write(errmsg, fmt='((a))') 'mp_thompson_post_run called before mp_thompson_post_init'
         errflg = 1
         return
      end if

      ! If limiter is deactivated, return immediately
      if (.not.apply_limiter) return

      ! mp_tend and ttendlim are expressed in potential temperature
      mp_tend = (tgrs - tgrs_save)/prslk

#ifdef DEBUG
      events = 0
#endif
      do k=1,nlev
         do i=1,ncol
            mp_tend(i,k) = max( -ttendlim*dtp, min( ttendlim*dtp, mp_tend(i,k) ) )

#ifdef DEBUG
            if (tgrs_save(i,k) + mp_tend(i,k)*prslk(i,k) .ne. tgrs(i,k)) then
              write(0,'(a,3i6,3e16.7)') "mp_thompson_post_run mp_tend limiter: kdt, i, k, t_old, t_new, t_lim:", &
                                      & kdt, i, k, tgrs_save(i,k), tgrs(i,k), tgrs_save(i,k) + mp_tend(i,k)*prslk(i,k)
              events = events + 1
            end if
#endif
            tgrs(i,k) = tgrs_save(i,k) + mp_tend(i,k)*prslk(i,k)
         end do
      end do

#ifdef DEBUG
      if (events > 0) then
        write(0,'(a,i0,a,i0,a,i0)') "mp_thompson_post_run: ttendlim applied ", events, "/", nlev*ncol, &
                                  & " times at timestep ", kdt
      end if
#endif

   end subroutine mp_thompson_post_run

!! \section arg_table_mp_thompson_post_finalize Argument Table
!! \htmlinclude mp_thompson_post_finalize.html
!!
   subroutine mp_thompson_post_finalize(errmsg, errflg)

      implicit none

      ! CCPP error handling
      character(len=*),          intent(  out) :: errmsg
      integer,                   intent(  out) :: errflg
      
      ! initialize ccpp error handling variables
      errmsg = ''
      errflg = 0

      ! Check initialization state
      if (.not. is_initialized) return

      is_initialized = .false.

   end subroutine mp_thompson_post_finalize

end module mp_thompson_post

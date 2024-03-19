module CCPP_driver

  use ccpp_types,         only: ccpp_t

  use ccpp_static_api,    only: ccpp_physics_init,                   &
                                ccpp_physics_timestep_init,          &
                                ccpp_physics_run,                    &
                                ccpp_physics_timestep_finalize,      &
                                ccpp_physics_finalize

  use CCPP_data,          only: cdata_tile,                          &
                                cdata_domain,                        &
                                cdata_block,                         &
                                ccpp_suite,                          &
                                GFS_control,                         &
                                GFS_data

  implicit none

!--------------------------------------------------------!
!  Pointer to CCPP containers defined in CCPP_data       !
!--------------------------------------------------------!
  type(ccpp_t), pointer :: cdata => null()

!--------------------------------------------------------!
!  Flag for non-uniform block sizes (last block smaller) !
!  and number of OpenMP threads (with special thread     !
!  number nthrdsX in case of non-uniform block sizes)    !
!--------------------------------------------------------!
  logical :: non_uniform_blocks
  integer :: nthrds, nthrdsX

!----------------
! Public Entities
!----------------
! functions
  public CCPP_step
! module variables
  public non_uniform_blocks

  CONTAINS
!*******************************************************************************************

!-------------------------------
!  CCPP step
!-------------------------------
  subroutine CCPP_step (step, nblks, ierr)

#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    character(len=*),         intent(in)  :: step
    integer,                  intent(in)  :: nblks
    integer,                  intent(out) :: ierr
    ! Local variables
    integer :: nb, nt, ntX
    integer :: ierr2
    ! DH* 20210104 - remove kdt_rad when code to clear diagnostic buckets is removed
    integer :: kdt_rad

    ierr = 0

    if (trim(step)=="init") then

      ! Get and set number of OpenMP threads (module
      ! variable) that are available to run physics
#ifdef _OPENMP
      nthrds = omp_get_max_threads()
#else
      nthrds = 1
#endif

      ! For non-uniform blocksizes, we use index nthrds+1
      ! for the interstitial data type with different length
      if (non_uniform_blocks) then
        nthrdsX = nthrds+1
      else
        nthrdsX = nthrds
      end if

      ! For physics running over the entire domain, block and thread
      ! number are not used; set to safe values
      cdata_domain%blk_no = 1
      cdata_domain%thrd_no = 1

      ! Allocate cdata structures for blocks and threads
      if (.not.allocated(cdata_block)) allocate(cdata_block(1:nblks,1:nthrdsX))

      ! Loop over all blocks and threads
      do nt=1,nthrdsX
        do nb=1,nblks
          ! Assign the correct block and thread numbers
          cdata_block(nb,nt)%blk_no = nb
          cdata_block(nb,nt)%thrd_no = nt
        end do
      end do

    else if (trim(step)=="physics_init") then

      ! Since the physics init step is independent of the blocking structure,
      ! we can use cdata_domain. And since we don't use threading on the host
      ! model side, we can allow threading inside the physics init routines.
      GFS_control%nthreads = nthrds

      call ccpp_physics_init(cdata_domain, suite_name=trim(ccpp_suite), ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_init"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

    ! Timestep init = time_vary
    else if (trim(step)=="timestep_init") then

      ! Since the physics timestep init step is independent of the blocking structure,
      ! we can use cdata_domain. And since we don't use threading on the host
      ! model side, we can allow threading inside the timestep init (time_vary) routines.
      GFS_control%nthreads = nthrds

      call ccpp_physics_timestep_init(cdata_domain, suite_name=trim(ccpp_suite), group_name="time_vary", ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_timestep_init for group time_vary"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

      ! call timestep_init for "physics"
      call ccpp_physics_timestep_init(cdata_domain, suite_name=trim(ccpp_suite),group_name="physics", ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_timestep_init for group physics"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DH* 20210104 - this block of code will be removed once the CCPP framework    !
      ! fully supports handling diagnostics through its metadata, work in progress   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !--- determine if radiation diagnostics buckets need to be cleared
      if (nint(GFS_control%fhzero*3600) >= nint(max(GFS_control%fhswr,GFS_control%fhlwr))) then
        if (mod(GFS_control%kdt,GFS_control%nszero) == 1) then
          do nb = 1,nblks
            call GFS_data(nb)%Intdiag%rad_zero(GFS_control)
          end do
        endif
      else
        kdt_rad = nint(min(GFS_control%fhswr,GFS_control%fhlwr)/GFS_control%dtp)
        if (mod(GFS_control%kdt,kdt_rad) == 1) then
          do nb = 1,nblks
            call GFS_data(nb)%Intdiag%rad_zero(GFS_control)
          enddo
        endif
      endif

      !--- determine if physics diagnostics buckets need to be cleared
      if ((mod(GFS_control%kdt-1,GFS_control%nszero)) == 0) then
        do nb = 1,nblks
          call GFS_data(nb)%Intdiag%phys_zero(GFS_control)
        end do
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! *DH 20210104                                                                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Radiation, physics and and stochastic physics - threaded regions using blocked data structures
    else if (trim(step)=="radiation" .or. trim(step)=="physics" .or. trim(step)=="stochastics") then

      ! Set number of threads available to physics schemes to one,
      ! because threads are used on the host model side for blocking
      GFS_control%nthreads = 1

!$OMP parallel num_threads (nthrds)      &
!$OMP          default (shared)          &
!$OMP          private (nb,nt,ntX,ierr2) &
!$OMP          reduction (+:ierr)
#ifdef _OPENMP
      nt = omp_get_thread_num()+1
#else
      nt = 1
#endif
!$OMP do schedule (dynamic,1)
      do nb = 1,nblks
        ! For non-uniform blocks, the last block has a different (shorter)
        ! length than the other blocks; use special CCPP_Interstitial(nthrdsX)
        if (non_uniform_blocks .and. nb==nblks) then
            ntX = nthrdsX
        else
            ntX = nt
        end if
        !--- Call CCPP radiation/physics/stochastics group
        call ccpp_physics_run(cdata_block(nb,ntX), suite_name=trim(ccpp_suite), group_name=trim(step), ierr=ierr2)
        if (ierr2/=0) then
           write(0,'(2a,3(a,i4),a)') "An error occurred in ccpp_physics_run for group ", trim(step), &
                                     ", block ", nb, " and thread ", nt, " (ntX=", ntX, "):"
           write(0,'(a)') trim(cdata_block(nb,ntX)%errmsg)
           ierr = ierr + ierr2
        end if
      end do
!$OMP end do

!$OMP end parallel
      if (ierr/=0) return

    ! Timestep finalize = time_vary
    else if (trim(step)=="timestep_finalize") then

      ! Since the physics timestep finalize step is independent of the blocking structure,
      ! we can use cdata_domain. And since we don't use threading on the host model side,
      ! we can allow threading inside the timestep finalize (time_vary) routines.
      GFS_control%nthreads = nthrds

      call ccpp_physics_timestep_finalize(cdata_domain, suite_name=trim(ccpp_suite), group_name="time_vary", ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_timestep_finalize for group time_vary"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

    ! Physics finalize
    else if (trim(step)=="physics_finalize") then

      ! Since the physics finalize step is independent of the blocking structure,
      ! we can use cdata_domain. And since we don't use threading on the host
      ! model side, we can allow threading inside the physics finalize routines.
      GFS_control%nthreads = nthrds

      call ccpp_physics_finalize(cdata_domain, suite_name=trim(ccpp_suite), ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_finalize"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

    ! Finalize
    else if (trim(step)=="finalize") then
      ! Deallocate cdata structure for blocks and threads
      if (allocated(cdata_block)) deallocate(cdata_block)

    else

      write(0,'(2a)') 'Error, undefined CCPP step ', trim(step)
      ierr = 1
      return

    end if

  end subroutine CCPP_step

end module CCPP_driver

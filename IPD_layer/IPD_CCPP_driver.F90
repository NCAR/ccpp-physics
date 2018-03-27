module IPD_CCPP_driver

  use IPD_typedefs,       only: IPD_init_type,                       &
                                IPD_control_type,  IPD_data_type,    &
                                IPD_diag_type,     IPD_restart_type, &
                                IPD_interstitial_type

  use ccpp_types,         only: ccpp_t
  use ccpp_errors,        only: ccpp_error, ccpp_debug
  use ccpp,               only: ccpp_init, ccpp_finalize
  use ccpp_fcall,         only: ccpp_run
  use ccpp_fields,        only: ccpp_field_add

! Begin include auto-generated list of modules for ccpp
#include "ccpp_modules.inc"
! End include auto-generated list of modules for ccpp

  use iso_c_binding,      only: c_loc

  implicit none

!------------------------------------------------------!
!  CCPP container                                      !
!------------------------------------------------------!
  type(ccpp_t), save, target :: cdata
  type(ccpp_t), dimension(:,:), allocatable, save, target :: cdata_block

!----------------
! Public Entities
!----------------
! functions
  public IPD_step

  CONTAINS
!*******************************************************************************************

  !-------------------------------
  !  IPD step generalized for CCPP
  !-------------------------------
  subroutine IPD_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, IPD_Interstitial, &
                       nBlocks, Init_parm, l_salp_data, l_snupx, ccpp_suite, step, ierr)

    use namelist_soilveg,  only: salp_data, snupx, max_vegtyp
    use block_control_mod, only: block_control_type
    use IPD_typedefs,      only: kind_phys
#ifdef OPENMP
    use omp_lib
#endif

    implicit none

    type(IPD_control_type),      target, intent(inout)           :: IPD_Control
    type(IPD_data_type),         target, intent(inout)           :: IPD_Data(:)
    type(IPD_diag_type),         target, intent(inout)           :: IPD_Diag(:)
    type(IPD_restart_type),      target, intent(inout)           :: IPD_Restart
    type(IPD_interstitial_type), target, intent(inout)           :: IPD_Interstitial(:)
    integer,                     target, intent(in)              :: nBlocks
    type(IPD_init_type),         target, intent(in)   , optional :: Init_parm
    real(kind=kind_phys),                intent(inout), optional :: l_salp_data
    real(kind=kind_phys),                intent(inout), optional :: l_snupx(max_vegtyp)
    character(len=256),                  intent(in),    optional :: ccpp_suite
    integer,                             intent(in)              :: step
    integer,                             intent(out)             :: ierr
    ! Local variables
    integer :: nb
    integer :: nThreads, nt

    ierr = 0

#ifdef OPENMP
    nThreads = omp_get_max_threads()
#else
    nThreads = 1
#endif

    if (step==0) then

      if (.not. present(Init_parm)) then
        call ccpp_error('Error, IPD init step called without mandatory Init_parm argument')
        ierr = 1
        return
      else if (.not. present(l_salp_data)) then
        call ccpp_error('Error, IPD init step called without mandatory l_salp_data argument')
        ierr = 1
        return
      else if (.not. present(l_snupx)) then
        call ccpp_error('Error, IPD init step called without mandatory l_snupx argument')
        ierr = 1
        return
      else if (.not. present(ccpp_suite)) then
        call ccpp_error('Error, IPD init step called without mandatory ccpp_suite argument')
        ierr = 1
        return
      end if

      call ccpp_init(ccpp_suite, cdata, ierr)
      if (ierr/=0) return

      !--- Add the DDTs to the CCPP data structure for IPD initialization
      call ccpp_field_add(cdata, 'IPD_Control',      '', c_loc(IPD_Control),      ierr=ierr)
      if (ierr/=0) return
      call ccpp_field_add(cdata, 'IPD_Data',         '', c_loc(IPD_Data),         rank=size(shape(IPD_Data)),         dims=shape(IPD_Data),         ierr=ierr)
      if (ierr/=0) return
      call ccpp_field_add(cdata, 'IPD_Diag',         '', c_loc(IPD_Diag),         rank=size(shape(IPD_Diag)),         dims=shape(IPD_Diag),         ierr=ierr)
      if (ierr/=0) return
      call ccpp_field_add(cdata, 'IPD_Restart',      '', c_loc(IPD_Restart),      ierr=ierr)
      if (ierr/=0) return
      call ccpp_field_add(cdata, 'IPD_Interstitial', '', c_loc(IPD_Interstitial), rank=size(shape(IPD_Interstitial)), dims=shape(IPD_Interstitial), ierr=ierr)
      if (ierr/=0) return
      call ccpp_field_add(cdata, 'Init_parm',        '', c_loc(Init_parm),        ierr=ierr)
      if (ierr/=0) return
      call ccpp_field_add(cdata, 'salp_data',            l_salp_data,             ierr=ierr)
      if (ierr/=0) return
      call ccpp_field_add(cdata, 'snupx',                l_snupx,                 ierr=ierr)
      if (ierr/=0) return

      call ccpp_run(cdata%suite%init, cdata, ierr)
      if (ierr/=0) return

      ! Allocate cdata structures
      allocate(cdata_block(1:nBlocks,1:nThreads))

#ifndef __PGI
      ! Loop over blocks for each of the threads
!$OMP parallel num_threads (nThreads) &
!$OMP          default (shared) &
!$OMP          private (nb,nt) &
!$OMP          reduction (+:ierr)
#ifdef OPENMP
      nt = omp_get_thread_num()+1
#else
      nt = 1
#endif
#else
      do nt=1,nThreads
#endif
      do nb = 1,nBlocks
#ifndef __PGI
        !--- Initialize CCPP, use suite from scalar cdata to avoid reading the SDF multiple times
        call ccpp_init(ccpp_suite, cdata_block(nb,nt), ierr, suite=cdata%suite)
#else
        !--- Initialize CCPP, cannot use suite from scalar cdata with PGI (crashes)
        call ccpp_init(ccpp_suite, cdata_block(nb,nt), ierr)
#endif
        if (ierr/=0) then
          write(0,'(2(a,i4))') "An error occurred in IPD_step 0 for block ", nb, " and thread ", nt
          exit
        end if
! Begin include auto-generated list of calls to ccpp_field_add
#include "ccpp_fields.inc"
! End include auto-generated list of calls to ccpp_field_add
      end do
#ifndef __PGI
!$OMP end parallel
#else
      end do
#endif
      if (ierr/=0) return

    ! Time vary steps
    else if (step==1) then

      ! Loop over blocks; cannot use OpenMP for this step (sfcsub.F!);
      ! however, threading may be implemented inside the IPD_setup_step
      do nb = 1,nBlocks
        nt = 1
        call ccpp_run(cdata_block(nb,nt)%suite%ipds(step), cdata_block(nb,nt), ierr)
        if (ierr/=0) then
            write(0,'(2(a,i4),a)') "An error occurred in IPD_step 1 for block ", nb, " and thread ", nt, &
                                 & "; error message: '" // trim(IPD_Interstitial(nt)%errmsg) // "'"
            return
        end if
      end do

    ! Radiation, physics and stochastics
    else if (step==2 .or. step==3 .or. step==4) then

!$OMP parallel do num_threads (nThreads) &
!$OMP            default (none) &
!$OMP            schedule (dynamic,1) &
!$OMP            shared   (nBlocks, cdata_block, step, IPD_Interstitial) &
!$OMP            private  (nb, nt) &
!$OMP            reduction (+:ierr)
      do nb = 1,nBlocks
#ifdef OPENMP
        nt = omp_get_thread_num()+1
#else
        nt = 1
#endif
        call ccpp_run(cdata_block(nb,nt)%suite%ipds(step), cdata_block(nb,nt), ierr)
        if (ierr/=0) then
          write(0,'(3(a,i4),a)') "An error occurred in IPD_step ", step, " for block ", nb, " and thread ", nt, &
                               & "; error message: '" // trim(IPD_Interstitial(nt)%errmsg) // "'"
        end if
      end do
!$OMP end parallel do
      if (ierr/=0) return

    ! Finalize
    else if (step==5) then

!$OMP parallel num_threads (nThreads) &
!$OMP          default (shared) &
!$OMP          private (nb,nt) &
!$OMP          reduction (+:ierr)
#ifdef OPENMP
      nt = omp_get_thread_num()+1
#else
      nt = 1
#endif
      do nb = 1,nBlocks
        !--- Initialize CCPP
        call ccpp_finalize(cdata_block(nb,nt), ierr)
        if (ierr/=0) then
           write(0,'(a,i4,a,i4)') "An error occurred in IPD_step 5 for block ", nb, " and thread ", nt
           exit
        end if
      end do
!$OMP end parallel
      if (ierr/=0) return

      ! Deallocate cdata structure for blocks and threads
      deallocate(cdata_block)

      call ccpp_finalize(cdata, ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in IPD_step 5"
      end if

    else

      call ccpp_error('Error, undefined step for ccpp_run')
      ierr = 1
      return

    end if

    ! DH* TODO CLEAN UP STDIO (USE FV3 MESSAGING? WRITE STATEMENTS? BE CONSISTENT!) *DH
  end subroutine IPD_step

end module IPD_CCPP_driver

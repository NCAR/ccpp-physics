module IPD_CCPP_driver

  use IPD_typedefs,       only: IPD_init_type,                       &
                                IPD_control_type,  IPD_data_type,    &
                                IPD_diag_type,     IPD_restart_type
  use ccpp_types,         only: ccpp_t
  use ccpp_errors,        only: ccpp_error, ccpp_debug
  use ccpp,               only: ccpp_init
  use ccpp_fcall,         only: ccpp_run
  use ccpp_fields,        only: ccpp_fields_add

#ifdef CCXX
! Begin include auto-generated list of modules for ccpp
#include "ccpp_modules.inc"
! End include auto-generated list of modules for ccpp
#endif

  use iso_c_binding,      only: c_loc

  implicit none

!------------------------------------------------------!
!  CCPP container                                      !
!------------------------------------------------------!
  type(ccpp_t), save, target :: cdata
  type(ccpp_t), dimension(:), allocatable, save, target :: cdata_block

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
  subroutine IPD_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, nBlocks, Atm_block, &
                       Init_parm, l_salp_data, l_snupx, ccpp_suite, step, ierr)

    use namelist_soilveg,  only: salp_data, snupx, max_vegtyp
    use block_control_mod, only: block_control_type
    use IPD_typedefs,      only: kind_phys

    implicit none

    type(IPD_control_type),    target, intent(inout)           :: IPD_Control
    type(IPD_data_type),       target, intent(inout)           :: IPD_Data(:)
    type(IPD_diag_type),       target, intent(inout)           :: IPD_Diag(:)
    type(IPD_restart_type),    target, intent(inout)           :: IPD_Restart
    integer,                   target, intent(in)              :: nBlocks
    type (block_control_type), target, intent(in)   , optional :: Atm_block
    type(IPD_init_type),       target, intent(in)   , optional :: Init_parm
    real(kind=kind_phys),              intent(inout), optional :: l_salp_data
    real(kind=kind_phys),              intent(inout), optional :: l_snupx(max_vegtyp)
    character(len=256),                intent(in),    optional :: ccpp_suite
    integer,                           intent(in)              :: step
    integer,                           intent(out)             :: ierr
    ! Local variables
    integer :: nb

    ierr = 0

    if (step==0) then

      if (.not. present(Atm_block)) then
        call ccpp_error('Error, IPD init step called without mandatory Atm_block argument')
        ierr = 1
        return
      else if (.not. present(Init_parm)) then
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

      !--- Add the DDTs to the CCPP data structure
      call ccpp_fields_add(cdata, 'IPD_Control', '', c_loc(IPD_Control), ierr=ierr)
      call ccpp_fields_add(cdata, 'IPD_Data',    '', c_loc(IPD_Data), rank=size(shape(IPD_Data)), dims=shape(IPD_Data), ierr=ierr)
      call ccpp_fields_add(cdata, 'IPD_Diag',    '', c_loc(IPD_Diag), rank=size(shape(IPD_Diag)), dims=shape(IPD_Diag), ierr=ierr)
      call ccpp_fields_add(cdata, 'IPD_Restart', '', c_loc(IPD_Restart), ierr=ierr)
      call ccpp_fields_add(cdata, 'Atm_block',   '', c_loc(Atm_block),   ierr=ierr)
      call ccpp_fields_add(cdata, 'Init_parm',   '', c_loc(Init_parm),   ierr=ierr)
      call ccpp_fields_add(cdata, 'salp_data',       l_salp_data,        ierr=ierr)
      call ccpp_fields_add(cdata, 'snupx',           l_snupx,            ierr=ierr)

      call ccpp_run(cdata%suite%init, cdata, ierr)

      ! Allocate cdata structures
      allocate(cdata_block(1:nBlocks))

      ! Loop over blocks - in general, cannot use OpenMP for this step;
      ! however, threading may be implemented inside the ccpp_init,
      ! suite_init and scheme_init routines.
      do nb = 1,nBlocks

         !--- Initialize CCPP
         call ccpp_init(ccpp_suite, cdata_block(nb), ierr)

#ifdef CCXX
! Begin include auto-generated list of calls to ccpp_fields_add
#include "ccpp_fields.inc"
! End include auto-generated list of calls to ccpp_fields_add
#endif
#ifdef CCPP
         !--- Add the DDTs to the CCPP data structure for this block
         call ccpp_fields_add(cdata_block(nb), 'IPD_Control', '', c_loc(IPD_Control),  ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'IPD_Data',    '', c_loc(IPD_Data(nb)), ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'IPD_Diag',    '', c_loc(IPD_Diag), rank=size(shape(IPD_Diag)), dims=shape(IPD_Diag), ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'IPD_Restart', '', c_loc(IPD_Restart),  ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'Atm_block',   '', c_loc(Atm_block),    ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'Init_parm',   '', c_loc(Init_parm),    ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'salp_data',       l_salp_data,         ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'snupx',           l_snupx,             ierr=ierr)
#endif
      end do

    else if (step==1) then

      ! Loop over blocks - in general, cannot use OpenMP for this step;
      ! however, threading may be implemented inside the IPD_setup_step
      ! DH* TODO - figure out in how to determine inside physics code
      ! whether OpenMP is supposed to be used at this level or whether
      ! threading is already implemented outside (i.e. here) *DH
      do nb = 1,nBlocks
#ifdef CCXX
        call ccpp_run(cdata_block(nb)%suite%ipds(step), cdata_block(nb), ierr)
#else
        call ccpp_run(cdata_block(nb)%suite%ipds(1)%subcycles(1)%schemes(step), cdata_block(nb), ierr)
#endif
      end do

    else if (step==2 .or. step==3 .or. step==4) then
      ! Radiation, physics and stochastics

!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (nBlocks, cdata_block, step) &
!$OMP            private  (nb, ierr)
      do nb = 1,nBlocks
#ifdef CCXX
        call ccpp_run(cdata_block(nb)%suite%ipds(step), cdata_block(nb), ierr)
#else
        call ccpp_run(cdata_block(nb)%suite%ipds(1)%subcycles(1)%schemes(step), cdata_block(nb), ierr)
#endif
      end do
!$OMP end parallel do

    else if (step==5) then
      ! DH* ccpp_run(cdata%suite%finalize, ...) not yet implemented
      deallocate(cdata_block)
    else
      call ccpp_error('Error, undefined step for ccpp_run')
      ierr = 1
      return
    end if

  end subroutine IPD_step

end module IPD_CCPP_driver

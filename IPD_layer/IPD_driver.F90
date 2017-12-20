module IPD_driver

  use IPD_typedefs,               only: IPD_init_type,                       &
                                        IPD_control_type,  IPD_data_type,    &
                                        IPD_diag_type,     IPD_restart_type

  use physics_abstraction_layer,  only: initialize,        time_vary_step,   &
                                        radiation_step1,   physics_step1,    &
                                        physics_step2

  use physics_diag_layer,         only: diag_populate

  use physics_restart_layer,      only: restart_populate

#ifdef CCPP
  use fms_mod,            only: error_mesg, FATAL
  use ccpp_types,         only: ccpp_t
  use ccpp,               only: ccpp_init
  use ccpp_fcall,         only: ccpp_run
  use ccpp_fields,        only: ccpp_fields_add
! Begin include auto-generated list of modules for ccpp
#include "ccpp_modules.inc"
! End include auto-generated list of modules for ccpp
  use iso_c_binding,      only: c_loc
#endif

       implicit none

#ifdef CCPP
!------------------------------------------------------!
!  CCPP container                                      !
!------------------------------------------------------!
type(ccpp_t), save, target :: cdata
type(ccpp_t), dimension(:), allocatable, save, target :: cdata_block
#endif

!------------------------------------------------------!
!  IPD containers                                      !
!------------------------------------------------------!
!  type(GFS_control_type)              :: IPD_Control  !
!  type(IPD_data_type)     allocatable :: IPD_Data(:)  !
!  type(IPD_diag_type),                :: IPD_Diag(:)  !
!  type(IPD_restart_type),             :: IPD_Restart  !
!------------------------------------------------------!

!----------------
! Public Entities
!----------------
! functions
  public IPD_initialize
  public IPD_setup_step 
  public IPD_radiation_step
  public IPD_physics_step1
  public IPD_physics_step2
#ifdef CCPP
  public IPD_step
#endif

  CONTAINS
!*******************************************************************************************


!----------------
!  IPD Initialize 
!----------------
  subroutine IPD_initialize (IPD_control, IPD_Data, IPD_Diag, IPD_Restart, IPD_init_parm)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data(:)
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart
    type(IPD_init_type),    intent(in)    :: IPD_init_parm

    !--- initialize the physics suite
    call initialize (IPD_Control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout,      &
                     IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid, &
                     IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,   &
                     IPD_Data(:)%Intdiag, IPD_Data(:)%Sfccycle, IPD_init_parm)


    !--- populate/associate the Diag container elements
    call diag_populate (IPD_Diag, IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                                  IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                                  IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                                  IPD_Data%Intdiag, IPD_Data%Sfccycle, IPD_init_parm)


    !--- allocate and populate/associate the Restart container elements
    call restart_populate (IPD_Restart, IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                                        IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                                        IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                                        IPD_Data%Intdiag, IPD_Data%Sfccycle, IPD_init_parm)

  end subroutine IPD_initialize


!---------------------------------------------
!  IPD setup step
!    surface data cycling, random streams, etc
!---------------------------------------------
  subroutine IPD_setup_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call time_vary_step (IPD_Control, IPD_Data%Statein, IPD_Data%Stateout,   &
                         IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                         IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                         IPD_Data%Intdiag, IPD_Data%Sfccycle)

  end subroutine IPD_setup_step


!--------------------
!  IPD radiation step
!--------------------
  subroutine IPD_radiation_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call radiation_step1 (IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                          IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                          IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                          IPD_Data%Intdiag)

  end subroutine IPD_radiation_step


!-------------------
!  IPD physics step1
!-------------------
  subroutine IPD_physics_step1 (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call physics_step1 (IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                        IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                        IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                        IPD_Data%Intdiag)

  end subroutine IPD_physics_step1


!-------------------
!  IPD physics step2
!-------------------
  subroutine IPD_physics_step2 (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call physics_step2 (IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                        IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                        IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                        IPD_Data%Intdiag)

  end subroutine IPD_physics_step2


#ifdef CCPP
  !-------------------------------
  !  IPD step generalized for CCPP
  !-------------------------------
  subroutine IPD_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, nBlocks, Atm_block, Init_parm, l_salp_data, l_snupx, ccpp_suite, step)

    use namelist_soilveg,  only: salp_data, snupx, max_vegtyp
    use block_control_mod, only: block_control_type
    use IPD_typedefs,      only: kind_phys

    implicit none

    type(IPD_control_type),    intent(inout)           :: IPD_Control
    type(IPD_data_type),       intent(inout)           :: IPD_Data(:)
    type(IPD_diag_type),       intent(inout)           :: IPD_Diag(:)
    type(IPD_restart_type),    intent(inout)           :: IPD_Restart
    integer,                   intent(in)              :: nBlocks
    type (block_control_type), intent(in)   , optional :: Atm_block
    type(IPD_init_type),       intent(in)   , optional :: Init_parm
    real(kind=kind_phys),      intent(inout), optional :: l_salp_data
    real(kind=kind_phys),      intent(inout), optional :: l_snupx(max_vegtyp)
    character(len=256),        intent(in),    optional :: ccpp_suite
    integer,                   intent(in)              :: step
    ! Local variables
    integer :: nb
    integer :: ierr

    if (step==0) then

      if (.not. present(Atm_block)) then
        call error_mesg('ccpp-ipd', 'IPD init step called without mandatory Atm_block argument', FATAL)
      else if (.not. present(Init_parm)) then
        call error_mesg('ccpp-ipd', 'IPD init step called without mandatory Init_parm argument', FATAL)
      else if (.not. present(l_salp_data)) then
        call error_mesg('ccpp-ipd', 'IPD init step called without mandatory l_salp_data argument', FATAL)
      else if (.not. present(l_snupx)) then
        call error_mesg('ccpp-ipd', 'IPD init step called without mandatory l_snupx argument', FATAL)
      else if (.not. present(ccpp_suite)) then
        call error_mesg('ccpp-ipd', 'IPD init step called without mandatory ccpp_suite argument', FATAL)
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

! Begin include auto-generated list of calls to ccpp_fields_add
#include "ccpp_fields.inc"
! End include auto-generated list of calls to ccpp_fields_add

         !--- Add the DDTs to the CCPP data structure for this block
         call ccpp_fields_add(cdata_block(nb), 'IPD_Control', '', c_loc(IPD_Control),  ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'IPD_Data',    '', c_loc(IPD_Data(nb)), ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'IPD_Diag',    '', c_loc(IPD_Diag), rank=size(shape(IPD_Diag)), dims=shape(IPD_Diag), ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'IPD_Restart', '', c_loc(IPD_Restart),  ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'Atm_block',   '', c_loc(Atm_block),    ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'Init_parm',   '', c_loc(Init_parm),    ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'salp_data',       l_salp_data,         ierr=ierr)
         call ccpp_fields_add(cdata_block(nb), 'snupx',           l_snupx,             ierr=ierr)

      end do

    else if (step==1) then

      ! Loop over blocks - in general, cannot use OpenMP for this step;
      ! however, threading may be implemented inside the IPD_setup_step
      ! DH* TODO - figure out in how to determine inside physics code
      ! whether OpenMP is supposed to be used at this level or whether
      ! threading is already implemented outside (i.e. here) *DH
      do nb = 1,nBlocks
        call ccpp_run(cdata_block(nb)%suite%ipds(1)%subcycles(1)%schemes(step), cdata_block(nb), ierr)
      end do

    ! DH* TODO: is the number of steps available from CCPP? then do step>1 and step < N-1 here
    else if (step==2 .or. step==3 .or. step==4) then

!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (nBlocks, cdata_block, step) &
!$OMP            private  (nb, ierr)
      do nb = 1,nBlocks
        call ccpp_run(cdata_block(nb)%suite%ipds(1)%subcycles(1)%schemes(step), cdata_block(nb), ierr)
      end do
!$OMP end parallel do

    ! DH* TODO: is the number of steps available from CCPP? then do step=N here
    else if (step==5) then
      ! DH* ccpp_run(cdata%suite%finalize, ...) not yet implemented
      deallocate(cdata_block)
    else
      call error_mesg('ccpp-ipd', 'IPD init step called without mandatory ccpp_suite argument', FATAL)
    end if
  end subroutine IPD_step
#endif

end module IPD_driver

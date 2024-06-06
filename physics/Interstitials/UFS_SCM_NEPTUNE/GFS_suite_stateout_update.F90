! #########################################################################################
!> \file GFS_suite_stateout_update.f90
!!  Update the state variables due to process-split physics from accumulated tendencies 
!!  during that phase.
!!  Update gas concentrations, if using prognostic photolysis schemes.
!!  Also, set bounds on the mass-weighted rime factor when using Ferrier-Aligo microphysics.
! #########################################################################################
module GFS_suite_stateout_update
  use machine,       only: kind_phys
  use module_ozphys, only: ty_ozphys
  implicit none
contains
! #########################################################################################
!> \section arg_table_GFS_suite_stateout_update_run Argument Table
!! \htmlinclude GFS_suite_stateout_update_run.html
!!
! #########################################################################################
  subroutine GFS_suite_stateout_update_run (im, levs, ntrac, dtp, tgrs, ugrs, vgrs, qgrs, &
       dudt, dvdt, dtdt, dqdt, gt0, gu0, gv0, gq0, oz0, ntiw, nqrimef, imp_physics,       &
       imp_physics_fer_hires, epsq, ozphys, oz_phys_2015, oz_phys_2006, con_1ovg, prsl,   &
       dp, ozpl, qdiag3d, do3_dt_prd, do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz, errmsg, errflg)

    ! Inputs
    integer,              intent(in )                   :: im
    integer,              intent(in )                   :: levs
    integer,              intent(in )                   :: ntrac
    integer,              intent(in )                   :: imp_physics,imp_physics_fer_hires
    integer,              intent(in )                   :: ntiw, nqrimef
    real(kind=kind_phys), intent(in )                   :: dtp, epsq, con_1ovg
    real(kind=kind_phys), intent(in ), dimension(:,:)   :: tgrs, ugrs, vgrs, prsl, dp
    real(kind=kind_phys), intent(in ), dimension(:,:,:) :: qgrs, ozpl
    real(kind=kind_phys), intent(in ), dimension(:,:)   :: dudt, dvdt, dtdt
    real(kind=kind_phys), intent(in ), dimension(:,:,:) :: dqdt
    logical,              intent(in)                    :: qdiag3d
    logical,              intent(in)                    :: oz_phys_2015
    logical,              intent(in)                    :: oz_phys_2006
    type(ty_ozphys),      intent(in)                    :: ozphys

    ! Outputs (optional)
    real(kind=kind_phys), intent(inout), dimension(:,:), optional :: &
         do3_dt_prd,  & ! Physics tendency: production and loss effect
         do3_dt_ozmx, & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp, & ! Physics tendency: temperature effect
         do3_dt_ohoz    ! Physics tendency: overhead ozone effect

    ! Outputs
    real(kind=kind_phys), intent(out), dimension(:,:)   :: gt0, gu0, gv0, oz0
    real(kind=kind_phys), intent(out), dimension(:,:,:) :: gq0
    character(len=*),     intent(out)                   :: errmsg
    integer,              intent(out)                   :: errflg

    ! Locals
    integer :: i, k

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Update prognostic state varaibles using accumulated tendencies from "process-split"
    ! section of GFS suite.
    gt0(:,:)   = tgrs(:,:)   + dtdt(:,:)   * dtp
    gu0(:,:)   = ugrs(:,:)   + dudt(:,:)   * dtp
    gv0(:,:)   = vgrs(:,:)   + dvdt(:,:)   * dtp
    gq0(:,:,:) = qgrs(:,:,:) + dqdt(:,:,:) * dtp

    ! If using photolysis physics schemes, update (prognostic) gas concentrations using 
    ! updated state.
    if (oz_phys_2015) then
       call ozphys%run_o3prog_2015(con_1ovg, dtp, prsl, gt0, dp, ozpl, oz0, qdiag3d, &
            do3_dt_prd, do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz)
    endif
    if (oz_phys_2006) then
       call ozphys%run_o3prog_2006(con_1ovg, dtp, prsl, gt0, dp, ozpl, oz0, qdiag3d, &
            do3_dt_prd, do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz)
    endif

    ! If using Ferrier-Aligo microphysics, set bounds on the mass-weighted rime factor.
    if (imp_physics == imp_physics_fer_hires) then
       do k=1,levs
          do i=1,im
             if(gq0(i,k,ntiw) > epsq) then
                gq0(i,k,nqrimef) = max(1., gq0(i,k,nqrimef)/gq0(i,k,ntiw))
             else
                gq0(i,k,nqrimef) = 1.
             end if
          end do
       end do
    end if

  end subroutine GFS_suite_stateout_update_run

end module GFS_suite_stateout_update

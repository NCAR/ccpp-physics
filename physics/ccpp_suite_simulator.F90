! ########################################################################################
!
! Description: This suite simulates the evolution of the internal physics state 
!              represented by a CCPP Suite Definition File (SDF).
!
! To activate this suite it must be a) embedded within the SDF and b) activated through
! the physics namelist.
! The derived-data type "base_physics_process" contains the metadata needed to reconstruct 
! the temporal evolution of the state. An array of base_physics_process, physics_process,
! is populated by the host during initialization and passed to the physics. Additionally, 
! this type holds any data, or type-bound procedures, required by the suite simulator(s).
!
! For this initial demonstration we are using 2-dimensional (height, time) forcing data,
! which is on the same native vertical grid as the SCM. The dataset has a temporal 
! resolution of 1-hour, created by averaging all local times from a Tropical Warm Pool 
! International Cloud Experiment (TWPICE) case. This was to create a dataset with a 
! (constant) diurnal cycle.
!
! ########################################################################################
module ccpp_suite_simulator
  use machine, only: kind_phys
  use module_ccpp_suite_simulator, only: base_physics_process, sim_LWRAD, sim_SWRAD,    &
       sim_PBL, sim_GWD, sim_DCNV, sim_SCNV, sim_cldMP
  implicit none
  public ccpp_suite_simulator_run
contains

  ! ######################################################################################
  !
  ! SUBROUTINE ccpp_suite_simulator_run
  !
  ! ######################################################################################
!! \section arg_table_ccpp_suite_simulator_run
!! \htmlinclude ccpp_suite_simulator_run.html
!!
  subroutine ccpp_suite_simulator_run(do_ccpp_suite_sim, kdt, nCol, nLay, dtp, jdat,     &
       iactive_T, iactive_u, iactive_v, iactive_q, proc_start, proc_end, physics_process,&
       in_pre_active, in_post_active, tgrs, ugrs, vgrs, qgrs, active_phys_tend, gt0, gu0,&
       gv0, gq0, errmsg, errflg)

    ! Inputs
    logical,           intent(in)  :: do_ccpp_suite_sim
    integer,           intent(in)  :: kdt, nCol, nLay, jdat(8), iactive_T, iactive_u,    &
                                      iactive_v, iactive_q
    real(kind_phys),   intent(in)  :: dtp, tgrs(:,:), ugrs(:,:), vgrs(:,:), qgrs(:,:,:), &
                                      active_phys_tend(:,:,:)
    ! Outputs
    type(base_physics_process),intent(inout) :: physics_process(:)
    real(kind_phys), intent(inout) :: gt0(:,:), gu0(:,:), gv0(:,:), gq0(:,:)
    character(len=*),intent(out)   :: errmsg
    integer,         intent(out)   :: errflg
    integer,         intent(inout) :: proc_start, proc_end
    logical,         intent(inout) :: in_pre_active, in_post_active

    ! Locals
    integer :: iCol, year, month, day, hour, min, sec, iprc
    real(kind_phys), dimension(nCol,nLay) :: gt1, gu1, gv1, dTdt, dudt, dvdt, gq1, dqdt

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ccpp_suite_sim) return

    ! Current forecast time (Data-format specific)
    year  = jdat(1)
    month = jdat(2)
    day   = jdat(3)
    hour  = jdat(5)
    min   = jdat(6)
    sec   = jdat(7)

    ! Set state at beginning of the physics timestep.
    gt1(:,:)  = tgrs(:,:)
    gu1(:,:)  = ugrs(:,:)
    gv1(:,:)  = vgrs(:,:)
    gq1(:,:)  = qgrs(:,:,1)
    dTdt(:,:) = 0.
    dudt(:,:) = 0.
    dvdt(:,:) = 0.
    dqdt(:,:) = 0.

    !
    ! Set bookeeping indices
    !
    if (in_pre_active) then
       proc_start = 1
       proc_end   = max(1,physics_process(1)%iactive_scheme-1)
    endif
    if (in_post_active) then
       proc_start = physics_process(1)%iactive_scheme
       proc_end   = size(physics_process)
    endif

    !
    ! Simulate internal physics timestep evolution.
    !
    do iprc = proc_start,proc_end
       do iCol = 1,nCol

          ! Reset locals
          physics_process(iprc)%tend1d%T(:) = 0.
          physics_process(iprc)%tend1d%u(:) = 0.
          physics_process(iprc)%tend1d%v(:) = 0.
          physics_process(iprc)%tend1d%q(:) = 0.

          ! Using scheme simulator
          ! Very simple...
          ! Interpolate 2D data (time,level) tendency to local time. 
          ! Here the data is already on the SCM vertical coordinate. 
          !
          ! In theory the data can be of any dimensionality and the onus falls on the 
          ! developer to extend the type "base_physics_process" to work with for their
          ! application. 
          !
          if (physics_process(iprc)%use_sim) then
             if (physics_process(iprc)%name == "LWRAD") then
                call sim_LWRAD(year, month, day, hour, min, sec, physics_process(iprc))
             endif
             if (physics_process(iprc)%name == "SWRAD")then
                call sim_SWRAD(year, month, day, hour, min, sec, physics_process(iprc))
             endif
             if (physics_process(iprc)%name == "GWD")then
                call sim_GWD(year, month, day, hour, min, sec, physics_process(iprc))
             endif
             if (physics_process(iprc)%name == "PBL")then
                call sim_PBL(year, month, day, hour, min, sec, physics_process(iprc))
             endif
             if (physics_process(iprc)%name == "SCNV")then
                call sim_SCNV(year, month, day, hour, min, sec, physics_process(iprc))
             endif
             if (physics_process(iprc)%name == "DCNV")then
                call sim_DCNV(year, month, day, hour, min, sec, physics_process(iprc))
             endif
             if (physics_process(iprc)%name == "cldMP")then
                call sim_cldMP(year, month, day, hour, min, sec, physics_process(iprc))
             endif

          ! Using data tendency from "active" scheme(s).
          else
             if (iactive_T > 0) physics_process(iprc)%tend1d%T = active_phys_tend(iCol,:,iactive_T)
             if (iactive_u > 0) physics_process(iprc)%tend1d%u = active_phys_tend(iCol,:,iactive_u)
             if (iactive_v > 0) physics_process(iprc)%tend1d%v = active_phys_tend(iCol,:,iactive_v)
             if (iactive_q > 0) physics_process(iprc)%tend1d%q = active_phys_tend(iCol,:,iactive_q)
          endif

          ! Update state now? (time-split scheme)
          if (physics_process(iprc)%time_split) then
             gt1(iCol,:)  = gt1(iCol,:) + (dTdt(iCol,:) + physics_process(iprc)%tend1d%T)*dtp
             gu1(iCol,:)  = gu1(iCol,:) + (dudt(iCol,:) + physics_process(iprc)%tend1d%u)*dtp
             gv1(iCol,:)  = gv1(iCol,:) + (dvdt(iCol,:) + physics_process(iprc)%tend1d%v)*dtp
             gq1(iCol,:)  = gq1(iCol,:) + (dqdt(iCol,:) + physics_process(iprc)%tend1d%q)*dtp
             dTdt(iCol,:) = 0.
             dudt(iCol,:) = 0.
             dvdt(iCol,:) = 0.
             dqdt(iCol,:) = 0.
          ! Accumulate tendencies, update later? (process-split scheme)
          else
             dTdt(iCol,:) = dTdt(iCol,:) + physics_process(iprc)%tend1d%T
             dudt(iCol,:) = dudt(iCol,:) + physics_process(iprc)%tend1d%u
             dvdt(iCol,:) = dvdt(iCol,:) + physics_process(iprc)%tend1d%v
             dqdt(iCol,:) = dqdt(iCol,:) + physics_process(iprc)%tend1d%q
          endif
       enddo ! END: Loop over columns

       ! Print diagnostics
       if (physics_process(iprc)%use_sim) then
          if (physics_process(iprc)%time_split) then
             write(*,'(a25,i2,a4,i2,a5,a10,a35)') 'CCPP suite simulator: ',iprc,' of ',proc_end,' ',physics_process(iprc)%name,'time split scheme     (simulated)'
          else
             write(*,'(a25,i2,a4,i2,a5,a10,a35)') 'CCPP suite simulator: ',iprc,' of ',proc_end,' ',physics_process(iprc)%name,'process split scheme  (simulated)'
          endif
       else
          if (physics_process(iprc)%time_split) then
             write(*,'(a25,i2,a4,i2,a5,a10,a35)') 'CCPP suite simulator: ',iprc,' of ',proc_end,' ',physics_process(iprc)%name,'time split scheme        (active)'
          else
             write(*,'(a25,i2,a4,i2,a5,a10,a35)') 'CCPP suite simulator: ',iprc,' of ',proc_end,' ',physics_process(iprc)%name,'process split scheme     (active)'
          endif
          write(*,'(a25,i2)')                     '       # prog. vars.: ',physics_process(1)%nprg_active
       endif
    enddo    ! END: Loop over physics processes

    !
    ! Update state with accumulated tendencies (process-split only)
    ! (Suites where active scheme is last physical process)
    !
    iprc = minval([iprc,proc_end])
    if (.not. physics_process(iprc)%time_split) then
       do iCol = 1,nCol
          gt0(iCol,:) = gt1(iCol,:) + dTdt(iCol,:)*dtp
          gu0(iCol,:) = gu1(iCol,:) + dudt(iCol,:)*dtp
          gv0(iCol,:) = gv1(iCol,:) + dvdt(iCol,:)*dtp
          gq0(iCol,:) = gq1(iCol,:) + dqdt(iCol,:)*dtp
       enddo
    endif

    !
    ! Update bookeeping indices 
    !
    if (in_pre_active) then
       in_pre_active  = .false.
       in_post_active = .true.
    endif

    if (size(physics_process) == proc_end) then
       in_pre_active  = .true.
       in_post_active = .false.
    endif

  end subroutine ccpp_suite_simulator_run

end module ccpp_suite_simulator

! ########################################################################################
! 
! CCPP scheme to replace physics schemes with simulated data tendencies.
!
! Description:
!
! ########################################################################################
module ccpp_scheme_simulator
  use machine, only: kind_phys

  implicit none

  ! ########################################################################################
  ! Types used by the scheme simulator
  ! ########################################################################################
  ! Type containing 1D (time) physics tendencies.
  type phys_tend_1d
     real(kind_phys), dimension(:), pointer :: T
     real(kind_phys), dimension(:), pointer :: u
     real(kind_phys), dimension(:), pointer :: v
     real(kind_phys), dimension(:), pointer :: q
  end type phys_tend_1d

  ! Type containing 2D (lev,time) physics tendencies.
  type phys_tend_2d
     real(kind_phys), dimension(:),   pointer :: time
     real(kind_phys), dimension(:,:), pointer :: T
     real(kind_phys), dimension(:,:), pointer :: u
     real(kind_phys), dimension(:,:), pointer :: v
     real(kind_phys), dimension(:,:), pointer :: q
  end type phys_tend_2d

  ! Type containing 3D (loc,lev,time) physics tendencies.
  type phys_tend_3d
     real(kind_phys), dimension(:),     pointer :: time
     real(kind_phys), dimension(:),     pointer :: lon
     real(kind_phys), dimension(:),     pointer :: lat
     real(kind_phys), dimension(:,:,:), pointer :: T
     real(kind_phys), dimension(:,:,:), pointer :: u
     real(kind_phys), dimension(:,:,:), pointer :: v
     real(kind_phys), dimension(:,:,:), pointer :: q
  end type phys_tend_3d

  ! Type containing 4D (lon, lat,lev,time) physics tendencies.
  type phys_tend_4d
     real(kind_phys), dimension(:),       pointer :: time
     real(kind_phys), dimension(:,:),     pointer :: lon
     real(kind_phys), dimension(:,:),     pointer :: lat
     real(kind_phys), dimension(:,:,:,:), pointer :: T
     real(kind_phys), dimension(:,:,:,:), pointer :: u
     real(kind_phys), dimension(:,:,:,:), pointer :: v
     real(kind_phys), dimension(:,:,:,:), pointer :: q
  end type phys_tend_4d

  ! This type contains the meta information and data for each physics process.
  type base_physics_process
     character(len=16)  :: name
     logical            :: time_split = .false.
     logical            :: use_sim    = .false.
     integer            :: order
     type(phys_tend_1d) :: tend1d
     type(phys_tend_2d) :: tend2d
     type(phys_tend_3d) :: tend3d
     type(phys_tend_4d) :: tend4d
   contains
     generic,   public  :: linterp => linterp_1D, linterp_2D
     procedure, private :: linterp_1D
     procedure, private :: linterp_2D
     procedure, public  :: find_nearest_loc_2d_1d
     procedure, public  :: cmp_time_wts
  end type base_physics_process

  ! This array contains the governing information on how to advance the physics timestep.
  type(base_physics_process), dimension(:), allocatable :: &
       physics_process

  ! For time-split physics process we need to call this scheme twice in the SDF, once
  ! before the "active" scheme is called, and once after. This is because the active
  ! scheme uses an internal physics state that has been advanced forward by a subsequent
  ! physics process(es).
  character(len=16) :: active_name
  integer :: iactive_scheme
  integer :: proc_start, proc_end
  logical :: active_time_split_process=.false.
  logical :: in_pre_active = .true.
  logical :: in_post_active = .false.

  ! Set to true in data was loaded into "physics_process"
  logical :: do_ccpp_scheme_simulator=.false.

  public ccpp_scheme_simulator_run

contains

  ! ######################################################################################
  !
  ! SUBROUTINE ccpp_scheme_simulator_run
  !
  ! ######################################################################################
!! \section arg_table_ccpp_scheme_simulator_run
!! \htmlinclude ccpp_scheme_simulator_run.html
!!
  subroutine ccpp_scheme_simulator_run(kdt, dtp, jdat, tgrs, ugrs, vgrs, qgrs, dtidx,    &
       dtend, index_of_process_dcnv, index_of_process_longwave,                          &
       index_of_process_shortwave, index_of_process_scnv,                                &
       index_of_process_orographic_gwd, index_of_process_pbl, index_of_process_mp,       &
       index_of_temperature, index_of_x_wind, index_of_y_wind, ntqv, gt0, gu0, gv0, gq0, &
       dtdq_pbl, dtdq_mp, errmsg, errflg)

    ! Inputs
    integer,         intent(in) :: kdt, ntqv, index_of_process_dcnv,                     &
         index_of_process_longwave, index_of_process_shortwave, index_of_process_scnv,   &
         index_of_process_orographic_gwd, index_of_process_pbl, index_of_process_mp,     &
         index_of_temperature, index_of_x_wind, index_of_y_wind
    integer,         intent(in), dimension(8) :: jdat
    integer,         intent(in), dimension(:,:) :: dtidx
    real(kind_phys), intent(in) :: dtp
    real(kind_phys), intent(in), dimension(:,:) :: tgrs, ugrs, vgrs
    real(kind_phys), intent(in), dimension(:,:,:) :: qgrs, dtend

    ! Outputs
    real(kind_phys), intent(inout), dimension(:,:) :: gt0, gu0, gv0
    real(kind_phys), intent(inout), dimension(:,:) :: gq0, dtdq_pbl, dtdq_mp
    character(len=*),intent(out) :: errmsg
    integer,         intent(out) :: errflg

    ! Locals
    integer :: iCol, iLay, nCol, nLay, idtend, fcst_year, fcst_month, fcst_day,          &
         fcst_hour, fcst_min, fcst_sec, iprc, index_of_active_process
    real(kind_phys) :: w1, w2,hrofday
    real(kind_phys), dimension(:,:), allocatable :: gt1, gu1, gv1, dTdt, dudt, dvdt
    real(kind_phys), dimension(:,:), allocatable :: gq1, dqdt

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ccpp_scheme_simulator) return

    ! Current forecast time (Data-format specific)
    fcst_year  = jdat(1)
    fcst_month = jdat(2)
    fcst_day   = jdat(3)
    fcst_hour  = jdat(5)
    fcst_min   = jdat(6)
    fcst_sec   = jdat(7)

    ! Dimensions
    nCol = size(gq0(:,1))
    nLay = size(gq0(1,:))

    ! Allocate temporaries
    allocate(gt1(nCol,nLay), gu1(nCol,nLay), gv1(nCol,nLay), gq1(nCol,nLay))
    allocate(dTdt(nCol,nLay), dudt(nCol,nLay), dvdt(nCol,nLay), dqdt(nCol,nLay))

    ! Get tendency for "active" process.
    ! DJS2023: For the UFS and SCM, the physics tendencies are stored in a multi-dimensional
    ! array, CCPP standard_name = cumulative_change_of_state_variables.
    ! These are not the instantaneous physics tendencies that are applied to the state by the 
    ! physics schemes. Not all schemes output physics tendencies...
    ! Rather these are intended for diagnostic puposes and are accumulated over some interval.
    ! In the UFS/SCM this is controlled by the diagnostic bucket interval, namelist option "fhzero".
    ! For this to work, you need to clear the diagnostic buckets after each physics timestep when
    ! running in the UFS/SCM.
    ! In the SCM this is done by adding the following runtime options:
    ! --n_itt_out 1 --n_itt_diag 1
    !
    if (active_name == "LWRAD") index_of_active_process = index_of_process_longwave
    if (active_name == "SWRAD") index_of_active_process = index_of_process_shortwave
    if (active_name == "PBL")   index_of_active_process = index_of_process_pbl
    if (active_name == "GWD")   index_of_active_process = index_of_process_orographic_gwd
    if (active_name == "SCNV")  index_of_active_process = index_of_process_scnv
    if (active_name == "DCNV")  index_of_active_process = index_of_process_dcnv
    if (active_name == "cldMP") index_of_active_process = index_of_process_mp

    ! Set state at beginning of the physics timestep.
    gt1(:,:)  = tgrs(:,:)
    gu1(:,:)  = ugrs(:,:)
    gv1(:,:)  = vgrs(:,:)
    gq1(:,:)  = qgrs(:,:,1)
    dTdt(:,:) = 0.
    dudt(:,:) = 0.
    dvdt(:,:) = 0.
    dqdt(:,:) = 0.

    if (in_pre_active) then
       proc_start = 1
       proc_end   = iactive_scheme-1
    endif
    if (in_post_active) then
       proc_start = iactive_scheme
       proc_end   = size(physics_process)
    endif

    ! Internal physics timestep evolution.
    do iprc = proc_start,proc_end
       if (iprc == iactive_scheme .and. active_time_split_process) then
          print*,'Reached active process. ', iprc
       else
          print*,'Simulating ',iprc,' of ',proc_end
       endif

       do iCol = 1,nCol
          ! Reset locals
          physics_process(iprc)%tend1d%T(:) = 0.
          physics_process(iprc)%tend1d%u(:) = 0.
          physics_process(iprc)%tend1d%v(:) = 0.
          physics_process(iprc)%tend1d%q(:) = 0.

          ! Using scheme simulator (very simple, interpolate data tendency to local time)
          if (physics_process(iprc)%use_sim) then
             if (associated(physics_process(iprc)%tend2d%T)) then
                errmsg = physics_process(iprc)%linterp("T", fcst_year, fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec)
             endif
             if (associated(physics_process(iprc)%tend2d%u)) then
                errmsg = physics_process(iprc)%linterp("u", fcst_year, fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec)
             endif
             if (associated(physics_process(iprc)%tend2d%v)) then
                errmsg = physics_process(iprc)%linterp("v", fcst_year, fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec)
             endif
             if (associated(physics_process(iprc)%tend2d%q)) then
                errmsg = physics_process(iprc)%linterp("q", fcst_year, fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec)
             endif

          ! Using data tendency from "active" scheme(s).
          ! DJS2023: This block is very ufs specific. See Note Above.
          else
             idtend = dtidx(index_of_temperature,index_of_active_process)
             if (idtend >= 1) physics_process(iprc)%tend1d%T = dtend(iCol,:,idtend)/dtp
             !
             idtend = dtidx(index_of_x_wind,index_of_active_process)
             if (idtend >= 1) physics_process(iprc)%tend1d%u = dtend(iCol,:,idtend)/dtp
             !
             idtend = dtidx(index_of_y_wind,index_of_active_process)
             if (idtend >= 1) physics_process(iprc)%tend1d%v = dtend(iCol,:,idtend)/dtp
             !
             idtend = dtidx(100+ntqv,index_of_active_process)
             if (idtend >= 1) physics_process(iprc)%tend1d%q = dtend(iCol,:,idtend)/dtp
          endif

          ! Update state now?
          if (physics_process(iprc)%time_split) then
             gt1(iCol,:) = gt1(iCol,:) + (dTdt(iCol,:) + physics_process(iprc)%tend1d%T)*dtp
             gu1(iCol,:) = gu1(iCol,:) + (dudt(iCol,:) + physics_process(iprc)%tend1d%u)*dtp
             gv1(iCol,:) = gv1(iCol,:) + (dvdt(iCol,:) + physics_process(iprc)%tend1d%v)*dtp
             gq1(iCol,:) = gq1(iCol,:) + (dqdt(iCol,:) + physics_process(iprc)%tend1d%q)*dtp
             !dTdt(iCol,:) = 0.
             !dudt(iCol,:) = 0.
             !dvdt(iCol,:) = 0.
             !dqdt(iCol,:) = 0.
          ! Accumulate tendencies, update later?
          else
             dTdt(iCol,:) = dTdt(iCol,:) + physics_process(iprc)%tend1d%T
             dudt(iCol,:) = dudt(iCol,:) + physics_process(iprc)%tend1d%u
             dvdt(iCol,:) = dvdt(iCol,:) + physics_process(iprc)%tend1d%v
             dqdt(iCol,:) = dqdt(iCol,:) + physics_process(iprc)%tend1d%q
          endif
          ! These are needed by samfshalcnv
          if (trim(physics_process(iprc)%name) == "PBL") then
             dtdq_pbl(iCol,:) = physics_process(iprc)%tend1d%q
          endif
          if (trim(physics_process(iprc)%name) == "cldMP") then
             dtdq_mp(iCol,:) = physics_process(iprc)%tend1d%q
          endif
       enddo
       !
       do iLay=1,nLay
          !write(*,'(i3,4f13.6)') ilay, gq0(iCol,iLay) , gq1(iCol,iLay) , dqdt(iCol,iLay)*dtp, physics_process(iprc)%tend1d%q(iLay)*dtp
       enddo
       gt0(iCol,:) = gt1(iCol,:) + dTdt(iCol,:)*dtp
       gu0(iCol,:) = gu1(iCol,:) + dudt(iCol,:)*dtp
       gv0(iCol,:) = gv1(iCol,:) + dvdt(iCol,:)*dtp
       gq0(iCol,:) = gq1(iCol,:) + dqdt(iCol,:)*dtp
    enddo

    if (in_pre_active) then
       in_pre_active  = .false.
       in_post_active = .true.
    endif

    if (size(physics_process)+1 == iprc) then
       in_pre_active  = .true.
       in_post_active = .false.
    endif

    !
  end subroutine ccpp_scheme_simulator_run

  ! ####################################################################################
  ! Type-bound procedure to compute tendency profile for time-of-day. 
  !
  ! For use with 1D data (level, time) tendencies with diurnal (24-hr) forcing.
  ! ####################################################################################
  function linterp_1D(this, var_name, year, month, day, hour, minute, second) result(err_message)
    class(base_physics_process), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer, intent(in) :: year, month, day, hour, minute, second
    character(len=128) :: err_message
    integer :: ti(1), tf(1)
    real(kind_phys) :: w1, w2

    ! Interpolation weights
    call this%cmp_time_wts(year, month, day, hour, minute, second, w1, w2, ti, tf)

    select case(var_name)
    case("T")
       this%tend1d%T = w1*this%tend2d%T(:,ti(1)) + w2*this%tend2d%T(:,tf(1))
    case("u")
       this%tend1d%u = w1*this%tend2d%u(:,ti(1)) + w2*this%tend2d%u(:,tf(1))
    case("v")
       this%tend1d%v = w1*this%tend2d%v(:,ti(1)) + w2*this%tend2d%v(:,tf(1))
    case("q")
       this%tend1d%q = w1*this%tend2d%q(:,ti(1)) + w2*this%tend2d%q(:,tf(1))
    end select

  end function linterp_1D

  ! ####################################################################################
  ! Type-bound procedure to compute tendency profile for time-of-day.
  !
  ! For use with 2D data (location, level, time) tendencies with diurnal (24-hr) forcing.
  ! This assumes that the location dimension has a [longitude, latitude] associated with
  ! each location.
  ! ####################################################################################
  function linterp_2D(this, var_name, lon, lat, year, month, day, hour, minute, second) result(err_message)
    class(base_physics_process), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer, intent(in) :: year, month, day, hour, minute, second
    real(kind_phys), intent(in) :: lon, lat
    character(len=128) :: err_message
    integer :: ti(1), tf(1), iNearest
    real(kind_phys) :: w1, w2

    ! Interpolation weights (temporal)
    call this%cmp_time_wts(year, month, day, hour, minute, second, w1, w2, ti, tf)

    ! Grab data tendency closest to column [lon,lat]
    iNearest = this%find_nearest_loc_2d_1d(lon,lat)
    
    select case(var_name)
    case("T")
       this%tend1d%T = w1*this%tend3d%T(iNearest,:,ti(1)) + w2*this%tend3d%T(iNearest,:,tf(1))
    case("u")
       this%tend1d%u = w1*this%tend3d%u(iNearest,:,ti(1)) + w2*this%tend3d%u(iNearest,:,tf(1))
    case("v")
       this%tend1d%v = w1*this%tend3d%v(iNearest,:,ti(1)) + w2*this%tend3d%v(iNearest,:,tf(1))
    case("q")
       this%tend1d%q = w1*this%tend3d%q(iNearest,:,ti(1)) + w2*this%tend3d%q(iNearest,:,tf(1))
    end select
  end function linterp_2D

  ! ####################################################################################
  ! Type-bound procedure to find nearest location.
  !
  ! For use with linterp_2D, NOT YET IMPLEMENTED.
  ! ####################################################################################
  pure function find_nearest_loc_2d_1d(this, lon, lat)
    class(base_physics_process), intent(in) :: this
    real(kind_phys), intent(in) :: lon, lat
    integer :: find_nearest_loc_2d_1d

    find_nearest_loc_2d_1d = 1 
  end function find_nearest_loc_2d_1d

  ! ####################################################################################
  ! Type-bound procedure to compute linear interpolation weights for a diurnal (24-hour)
  ! forcing.
  ! ####################################################################################
  subroutine cmp_time_wts(this, year, month, day, hour, minute, second, w1, w2, ti, tf)
    ! Inputs
    class(base_physics_process), intent(in) :: this
    integer, intent(in) :: year, month, day, hour, minute, second
    ! Outputs
    integer,intent(out) :: ti(1), tf(1)
    real(kind_phys),intent(out) :: w1, w2
    ! Locals
    real(kind_phys) :: hrofday

    hrofday = hour*3600. + minute*60. + second
    ti = max(hour,1)
    tf = min(ti + 1,24)
    w1 = ((hour+1)*3600 - hrofday)/3600
    w2 = 1 - w1

  end subroutine cmp_time_wts

end module ccpp_scheme_simulator

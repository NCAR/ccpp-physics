! ########################################################################################
!
! This module contains the type, base_physics_process, and supporting subroutines needed
! by the ccpp suite simulator.
!
! ########################################################################################
module module_ccpp_suite_simulator
!> \section arg_table_module_ccpp_suite_simulator Argument table
!! \htmlinclude module_ccpp_suite_simulator.html
!!
  use machine, only : kind_phys
  implicit none
  
  public base_physics_process

  ! Type containing 1D (time) physics tendencies. 
  type phys_tend_1d
     real(kind_phys), dimension(:), allocatable :: T
     real(kind_phys), dimension(:), allocatable :: u
     real(kind_phys), dimension(:), allocatable :: v
     real(kind_phys), dimension(:), allocatable :: q
     real(kind_phys), dimension(:), allocatable :: p
     real(kind_phys), dimension(:), allocatable :: z
  end type phys_tend_1d

  ! Type containing 2D (lev,time) physics tendencies.
  type phys_tend_2d
     real(kind_phys), dimension(:),   allocatable :: time
     real(kind_phys), dimension(:,:), allocatable :: T
     real(kind_phys), dimension(:,:), allocatable :: u
     real(kind_phys), dimension(:,:), allocatable :: v
     real(kind_phys), dimension(:,:), allocatable :: q
     real(kind_phys), dimension(:,:), allocatable :: p
     real(kind_phys), dimension(:,:), allocatable :: z
  end type phys_tend_2d

  ! Type containing 3D (loc,lev,time) physics tendencies.
  type phys_tend_3d
     real(kind_phys), dimension(:),     allocatable :: time
     real(kind_phys), dimension(:),     allocatable :: lon
     real(kind_phys), dimension(:),     allocatable :: lat
     real(kind_phys), dimension(:,:,:), allocatable :: T
     real(kind_phys), dimension(:,:,:), allocatable :: u
     real(kind_phys), dimension(:,:,:), allocatable :: v
     real(kind_phys), dimension(:,:,:), allocatable :: q
  end type phys_tend_3d

  ! Type containing 4D (lon,lat,lev,time) physics tendencies.
  type phys_tend_4d
     real(kind_phys), dimension(:),       allocatable :: time
     real(kind_phys), dimension(:,:),     allocatable :: lon
     real(kind_phys), dimension(:,:),     allocatable :: lat
     real(kind_phys), dimension(:,:,:,:), allocatable :: T
     real(kind_phys), dimension(:,:,:,:), allocatable :: u
     real(kind_phys), dimension(:,:,:,:), allocatable :: v
     real(kind_phys), dimension(:,:,:,:), allocatable :: q
  end type phys_tend_4d

! This type contains the meta information and data for each physics process.

!> \section arg_table_base_physics_process Argument Table
!! \htmlinclude base_physics_process.html
!!
  type base_physics_process
     character(len=16)  :: name                 ! Physics process name
     logical            :: time_split = .false. ! Is process time-split?
     logical            :: use_sim    = .false. ! Is process "active"?
     integer            :: order                ! Order of process in process-loop 
     type(phys_tend_1d) :: tend1d               ! Instantaneous data
     type(phys_tend_2d) :: tend2d               ! 2-dimensional data
     type(phys_tend_3d) :: tend3d               ! Not used. Placeholder for 3-dimensional spatial data.
     type(phys_tend_4d) :: tend4d               ! Not used. Placeholder for 4-dimensional spatio-tempo data.
     character(len=16)  :: active_name          ! "Active" scheme: Physics process name
     integer            :: iactive_scheme       ! "Active" scheme: Order of process in process-loop
     logical            :: active_tsp           ! "Active" scheme: Is process time-split?
     integer            :: nprg_active          ! "Active" scheme: Number of prognostic variables
   contains
     generic,   public  :: linterp => linterp_1D, linterp_2D
     procedure, private :: linterp_1D
     procedure, private :: linterp_2D
     procedure, public  :: find_nearest_loc_2d_1d
     procedure, public  :: cmp_time_wts
  end type base_physics_process

contains

  ! ####################################################################################
  ! Type-bound procedure to compute tendency profile for time-of-day.
  !
  ! For use with 1D data (level, time) tendencies with diurnal (24-hr) forcing.
  ! ####################################################################################
  function linterp_1D(this, var_name, year, month, day, hour, min, sec) result(err_message)
    class(base_physics_process), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer, intent(in) :: year, month, day, hour, min, sec
    character(len=128) :: err_message
    integer :: ti(1), tf(1), ntime
    real(kind_phys) :: w1, w2

    ! Interpolation weights
    call this%cmp_time_wts(year, month, day, hour, min, sec, w1, w2, ti, tf)

    ntime = size(this%tend2d%T(1,:))

    select case(var_name)
    case("T")
       if (tf(1) .le. ntime) then
          this%tend1d%T = w1*this%tend2d%T(:,ti(1)) + w2*this%tend2d%T(:,tf(1))
       else
          this%tend1d%T = this%tend2d%T(:,1)
       endif
    case("u")
       if (tf(1) .le. ntime) then
          this%tend1d%u = w1*this%tend2d%u(:,ti(1)) + w2*this%tend2d%u(:,tf(1))
       else
          this%tend1d%u = this%tend2d%u(:,1)
       endif
    case("v")
       if (tf(1) .le. ntime) then
          this%tend1d%v = w1*this%tend2d%v(:,ti(1)) + w2*this%tend2d%v(:,tf(1))
       else
          this%tend1d%v = this%tend2d%v(:,1)
       endif
    case("q")
       if (tf(1) .le. ntime) then
          this%tend1d%q = w1*this%tend2d%q(:,ti(1)) + w2*this%tend2d%q(:,tf(1))
       else
          this%tend1d%q = this%tend2d%q(:,1)
       endif
    end select

  end function linterp_1D

  ! ####################################################################################
  ! Type-bound procedure to compute tendency profile for time-of-day.
  !
  ! For use with 2D data (location, level, time) tendencies with diurnal (24-hr) forcing.
  ! This assumes that the location dimension has a [longitude, latitude] allocated with
  ! each location.
  ! ####################################################################################
  function linterp_2D(this, var_name, lon, lat, year, month, day, hour, min, sec) result(err_message)
    class(base_physics_process), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer, intent(in) :: year, month, day, hour, min, sec
    real(kind_phys), intent(in) :: lon, lat
    character(len=128) :: err_message
    integer :: ti(1), tf(1), iNearest
    real(kind_phys) :: w1, w2

    ! Interpolation weights (temporal)
    call this%cmp_time_wts(year, month, day, hour, min, sec, w1, w2, ti, tf)

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
  subroutine cmp_time_wts(this, year, month, day, hour, minute, sec, w1, w2, ti, tf)
    ! Inputs
    class(base_physics_process), intent(in) :: this
    integer, intent(in) :: year, month, day, hour, minute, sec
    ! Outputs
    integer,intent(out) :: ti(1), tf(1)
    real(kind_phys),intent(out) :: w1, w2
    ! Locals
    real(kind_phys) :: hrofday

    hrofday = hour*3600. + minute*60. + sec
    ti = max(hour,1)
    tf = min(ti + 1,24)
    w1 = ((hour+1)*3600 - hrofday)/3600
    w2 = 1 - w1

  end subroutine cmp_time_wts

  ! ####################################################################################
  ! ####################################################################################
  subroutine sim_LWRAD( year, month, day, hour, min, sec, process)
    type(base_physics_process), intent(inout) :: process
    integer, intent(in) :: year, month, day, hour, min, sec
    character(len=128) :: errmsg

    if (allocated(process%tend2d%T)) then
       errmsg = process%linterp("T", year,month,day,hour,min,sec)
    endif

  end subroutine sim_LWRAD

  ! ####################################################################################
  ! ####################################################################################
  subroutine sim_SWRAD( year, month, day, hour, min, sec, process)
    type(base_physics_process), intent(inout) :: process
    integer, intent(in) :: year, month, day, hour, min, sec
    character(len=128) :: errmsg

    if (allocated(process%tend2d%T)) then
       errmsg = process%linterp("T", year,month,day,hour,min,sec)
    endif

  end subroutine sim_SWRAD

  ! ####################################################################################
  ! ####################################################################################
  subroutine sim_GWD( year, month, day, hour, min, sec, process)
    type(base_physics_process), intent(inout) :: process
    integer, intent(in) :: year, month, day, hour, min, sec
    character(len=128) :: errmsg

    if (allocated(process%tend2d%T)) then
       errmsg = process%linterp("T", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%u)) then
       errmsg = process%linterp("u", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%v)) then
       errmsg = process%linterp("v", year,month,day,hour,min,sec)
    endif

  end subroutine sim_GWD

  ! ####################################################################################
  ! ####################################################################################
  subroutine sim_PBL( year, month, day, hour, min, sec, process)
    type(base_physics_process), intent(inout) :: process
    integer, intent(in) :: year, month, day, hour, min, sec
    character(len=128) :: errmsg

    if (allocated(process%tend2d%T)) then
       errmsg = process%linterp("T", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%u)) then
       errmsg = process%linterp("u", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%v)) then
       errmsg = process%linterp("v", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%q)) then
       errmsg = process%linterp("q", year,month,day,hour,min,sec)
    endif

  end subroutine sim_PBL

  ! ####################################################################################
  ! ####################################################################################
  subroutine sim_DCNV( year, month, day, hour, min, sec, process)
    type(base_physics_process), intent(inout) :: process
    integer, intent(in) :: year, month, day, hour, min, sec
    character(len=128) :: errmsg

    if (allocated(process%tend2d%T)) then
       errmsg = process%linterp("T", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%u)) then
       errmsg = process%linterp("u", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%v)) then
       errmsg = process%linterp("v", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%q)) then
       errmsg = process%linterp("q", year,month,day,hour,min,sec)
    endif

  end subroutine sim_DCNV

  ! ####################################################################################
  ! ####################################################################################
  subroutine sim_SCNV( year, month, day, hour, min, sec, process)
    type(base_physics_process), intent(inout) :: process
    integer, intent(in) :: year, month, day, hour, min, sec
    character(len=128) :: errmsg

    if (allocated(process%tend2d%T)) then
       errmsg = process%linterp("T", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%u)) then
       errmsg = process%linterp("u", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%v)) then
       errmsg = process%linterp("v", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%q)) then
       errmsg = process%linterp("q", year,month,day,hour,min,sec)
    endif

  end subroutine sim_SCNV

  ! ####################################################################################
  ! ####################################################################################
  subroutine sim_cldMP( year, month, day, hour, min, sec, process)
    type(base_physics_process), intent(inout) :: process
    integer, intent(in) :: year, month, day, hour, min, sec
    character(len=128) :: errmsg

    if (allocated(process%tend2d%T)) then
       errmsg = process%linterp("T", year,month,day,hour,min,sec)
    endif
    if (allocated(process%tend2d%q)) then
       errmsg = process%linterp("q", year,month,day,hour,min,sec)
    endif
  end subroutine sim_cldMP

end module module_ccpp_suite_simulator

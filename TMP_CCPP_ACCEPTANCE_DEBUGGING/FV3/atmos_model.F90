!#define DHDEBUG
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module atmos_model_mod
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
use mpp_mod,            only: mpp_clock_end, CLOCK_COMPONENT, MPP_CLOCK_SYNC
use mpp_mod,            only: FATAL, mpp_min, mpp_max, mpp_error, mpp_chksum
use mpp_domains_mod,    only: domain2d
use mpp_mod,            only: mpp_get_current_pelist_name
#ifdef INTERNAL_FILE_NML
use mpp_mod,            only: input_nml_file
#else
use fms_mod,            only: open_namelist_file
#endif
use fms_mod,            only: file_exist, error_mesg
use fms_mod,            only: close_file, write_version_number, stdlog, stdout
use fms_mod,            only: clock_flag_default
use fms_mod,            only: check_nml_error
use diag_manager_mod,   only: diag_send_complete_instant
use time_manager_mod,   only: time_type, get_time, get_date, &
                              operator(+), operator(-),real_to_time_type
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers, get_tracer_names, &
                              get_tracer_index
use xgrid_mod,          only: grid_box_type
use atmosphere_mod,     only: atmosphere_init
use atmosphere_mod,     only: atmosphere_restart
use atmosphere_mod,     only: atmosphere_end
use atmosphere_mod,     only: atmosphere_state_update
use atmosphere_mod,     only: atmos_phys_driver_statein
use atmosphere_mod,     only: atmosphere_control_data
use atmosphere_mod,     only: atmosphere_resolution, atmosphere_domain
use atmosphere_mod,     only: atmosphere_grid_bdry, atmosphere_grid_ctr
use atmosphere_mod,     only: atmosphere_dynamics, atmosphere_diag_axes
use atmosphere_mod,     only: atmosphere_etalvls, atmosphere_hgt
!rab use atmosphere_mod,     only: atmosphere_tracer_postinit
use atmosphere_mod,     only: atmosphere_diss_est, atmosphere_nggps_diag
use atmosphere_mod,     only: atmosphere_scalar_field_halo
use atmosphere_mod,     only: atmosphere_get_bottom_layer
use atmosphere_mod,     only: set_atmosphere_pelist
use atmosphere_mod,     only: Atm, mytile
use block_control_mod,  only: block_control_type, define_blocks_packed
use DYCORE_typedefs,    only: DYCORE_data_type, DYCORE_diag_type
#ifdef CCPP
use IPD_typedefs,       only: IPD_init_type, IPD_diag_type,    &
                              IPD_restart_type, IPD_kind_phys, &
                              IPD_func0d_proc, IPD_func1d_proc
#else
use IPD_typedefs,       only: IPD_init_type, IPD_control_type, &
                              IPD_data_type, IPD_diag_type,    &
                              IPD_restart_type, IPD_kind_phys, &
                              IPD_func0d_proc, IPD_func1d_proc
#endif

#ifdef CCPP
use CCPP_data,          only: ccpp_suite,                      &
                              IPD_control => GFS_control,      &
                              IPD_data => GFS_data,            &
                              IPD_interstitial => GFS_interstitial
use IPD_driver,         only: IPD_initialize, IPD_step, IPD_finalize
use CCPP_driver,        only: CCPP_step, non_uniform_blocks
#ifdef HYBRID
use physics_abstraction_layer, only: physics_step1
#endif
#else
use IPD_driver,         only: IPD_initialize, IPD_step
use physics_abstraction_layer, only: time_vary_step, radiation_step1, physics_step1, physics_step2
#endif
use FV3GFS_io_mod,      only: FV3GFS_restart_read, FV3GFS_restart_write, &
                              FV3GFS_IPD_checksum,                       &
                              FV3GFS_diag_register, FV3GFS_diag_output,  &
                              DIAG_SIZE
use fv_iau_mod, only: iau_external_data_type,getiauforcing,iau_initialize
use module_fv3_config, only:  output_1st_tstep_rst, first_kdt

!-----------------------------------------------------------------------

implicit none

! DH*
#ifdef DHDEBUG

#ifdef __GFORTRAN__
#define PRINT_SUM
#else
#define PRINT_CHKSUM
#endif

interface print_var
  module procedure print_logic_0d
  module procedure print_int_0d
  module procedure print_int_1d
  module procedure print_real_0d
  module procedure print_real_1d
  module procedure print_real_2d
  module procedure print_real_3d
end interface

#endif
! *DH

private

public update_atmos_radiation_physics
public update_atmos_model_state
public update_atmos_model_dynamics
public atmos_model_init, atmos_model_end, atmos_data_type
public atmos_model_exchange_phase_1, atmos_model_exchange_phase_2
public atmos_model_restart
public get_atmos_model_ungridded_dim
public addLsmask2grid
!-----------------------------------------------------------------------

!<PUBLICTYPE >
 type atmos_data_type
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     integer                       :: layout(2)          ! computer task laytout
     logical                       :: regional           ! true if domain is regional
     logical                       :: nested             ! true if there is a nest
     integer                       :: mlon, mlat
     logical                       :: pe                 ! current pe.
     real(kind=8),             pointer, dimension(:)     :: ak, bk
     real,                     pointer, dimension(:,:)   :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real,                     pointer, dimension(:,:)   :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: lon      => null() ! local longitude axis grid box centers in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: lat      => null() ! local latitude axis grid box centers in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: dx, dy
     real(kind=8),             pointer, dimension(:,:)   :: area
     real(kind=8),             pointer, dimension(:,:,:) :: layer_hgt, level_hgt
     type(domain2d)                :: domain             ! domain decomposition
     type(time_type)               :: Time               ! current time
     type(time_type)               :: Time_step          ! atmospheric time step.
     type(time_type)               :: Time_init          ! reference time.
     type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange 
     type(IPD_diag_type), pointer, dimension(:) :: Diag
 end type atmos_data_type
                                                         ! to calculate gradient on cubic sphere grid.
!</PUBLICTYPE >

integer :: fv3Clock, getClock, updClock, setupClock, radClock, physClock

!-----------------------------------------------------------------------
integer :: blocksize    = 1
logical :: chksum_debug = .false.
logical :: dycore_only  = .false.
logical :: debug        = .false.
logical :: sync         = .false.
integer, parameter     :: maxhr = 4096
real, dimension(maxhr) :: fdiag = 0.
real                   :: fhmax=384.0, fhmaxhf=120.0, fhout=3.0, fhouthf=1.0, avg_max_length=3600.
#ifdef CCPP
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, fdiag, fhmax, fhmaxhf, fhout, fhouthf, ccpp_suite, avg_max_length
#else
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, fdiag, fhmax, fhmaxhf, fhout, fhouthf, avg_max_length
#endif

type (time_type) :: diag_time

!--- concurrent and decoupled radiation and physics variables
!-------------------
!  DYCORE containers
!-------------------
type(DYCORE_data_type),    allocatable :: DYCORE_Data(:)  ! number of blocks
type(DYCORE_diag_type)                 :: DYCORE_Diag(25)

!----------------
!  IPD containers
!----------------
#ifndef CCPP
type(IPD_control_type)              :: IPD_Control
type(IPD_data_type),    allocatable :: IPD_Data(:)  ! number of blocks
type(IPD_diag_type),    target      :: IPD_Diag(DIAG_SIZE)
type(IPD_restart_type)              :: IPD_Restart
#else
! IPD_Control and IPD_Data are coming from CCPP_data
type(IPD_diag_type),    target      :: IPD_Diag(DIAG_SIZE)
type(IPD_restart_type)              :: IPD_Restart
#endif

!--------------
! IAU container
!--------------
type(iau_external_data_type)        :: IAU_Data ! number of blocks

!-----------------
!  Block container
!-----------------
type (block_control_type), target   :: Atm_block

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

#ifdef NAM_phys
  logical,parameter :: flip_vc = .false.
#else
  logical,parameter :: flip_vc = .true.
#endif

contains

!#######################################################################
! <SUBROUTINE NAME="update_radiation_physics">
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down". 
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_radiation_physics (Atmos)
!   </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

subroutine update_atmos_radiation_physics (Atmos)
!-----------------------------------------------------------------------
  type (atmos_data_type), intent(in) :: Atmos
!--- local variables---
    integer :: nb, jdat(8), rc
    procedure(IPD_func0d_proc), pointer :: Func0d => NULL()
    procedure(IPD_func1d_proc), pointer :: Func1d => NULL()
#ifdef CCPP
    integer :: ierr
#endif
    if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "statein driver"
!--- get atmospheric state from the dynamic core
    call set_atmosphere_pelist()
    call mpp_clock_begin(getClock)
    if (IPD_control%do_skeb) call atmosphere_diss_est (IPD_control%skeb_npass) !  do smoothing for SKEB
    call atmos_phys_driver_statein (IPD_data, Atm_block, flip_vc)
    call mpp_clock_end(getClock)

!--- if dycore only run, set up the dummy physics output state as the input state
    if (dycore_only) then
      do nb = 1,Atm_block%nblks
        IPD_Data(nb)%Stateout%gu0 = IPD_Data(nb)%Statein%ugrs
        IPD_Data(nb)%Stateout%gv0 = IPD_Data(nb)%Statein%vgrs
        IPD_Data(nb)%Stateout%gt0 = IPD_Data(nb)%Statein%tgrs
        IPD_Data(nb)%Stateout%gq0 = IPD_Data(nb)%Statein%qgrs
      enddo
    else
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "setup step"

!--- update IPD_Control%jdat(8)
      jdat(:) = 0
      call get_date (Atmos%Time, jdat(1), jdat(2), jdat(3),  &
                                 jdat(5), jdat(6), jdat(7))
      IPD_Control%jdat(:) = jdat(:)

      ! DH*
#ifdef DHDEBUG
      write(0,*) "Calling MY_DIAGTOSCREEN before time_vary step"
      call MY_DIAGTOSCREEN()
#endif
      ! *DH

!--- execute the IPD atmospheric setup step
      call mpp_clock_begin(setupClock)
#ifdef CCPP
      call CCPP_step (step="time_vary", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to IPD-CCPP time_vary step failed')
#else
      Func1d => time_vary_step
      call IPD_step (IPD_Control, IPD_Data(:), IPD_Diag, IPD_Restart, IPD_func1d=Func1d)
#endif
!--- if coupled, assign coupled fields
      if( IPD_Control%cplflx .or. IPD_Control%cplwav ) then
!        print *,'in atmos_model,nblks=',Atm_block%nblks
!        print *,'in atmos_model,IPD_Data size=',size(IPD_Data)
!        print *,'in atmos_model,tsfc(1)=',IPD_Data(1)%sfcprop%tsfc(1)
!        print *,'in atmos_model, tsfc size=',size(IPD_Data(1)%sfcprop%tsfc)
        call assign_importdata(rc)
!        print *,'in atmos_model, after assign_importdata, rc=',rc
      endif

      call mpp_clock_end(setupClock)

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "radiation driver"

      ! DH*
#ifdef DHDEBUG
      write(0,*) "Calling MY_DIAGTOSCREEN before radiation step"
      call MY_DIAGTOSCREEN()
#endif
      ! *DH

!--- execute the IPD atmospheric radiation subcomponent (RRTM)

      call mpp_clock_begin(radClock)
#ifdef CCPP
      ! Performance improvement. Only enter if it is time to call the radiation physics.
      if (IPD_Control%lsswr .or. IPD_Control%lslwr) then
        call CCPP_step (step="radiation", nblks=Atm_block%nblks, ierr=ierr)
        if (ierr/=0)  call mpp_error(FATAL, 'Call to IPD-CCPP radiation step failed')
      endif
#else
      Func0d => radiation_step1
!$OMP parallel do default (none)       &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
      enddo
#endif
      call mpp_clock_end(radClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RADIATION STEP  ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "physics driver"

      ! DH*
#ifdef DHDEBUG
      write(0,*) "Calling MY_DIAGTOSCREEN before physics step"
      call MY_DIAGTOSCREEN()
#endif
      ! *DH

!--- execute the IPD atmospheric physics step1 subcomponent (main physics driver)

      call mpp_clock_begin(physClock)
#if defined(CCPP) && !defined(HYBRID)
      call CCPP_step (step="physics", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to IPD-CCPP physics step failed')
#else
      Func0d => physics_step1
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
#ifdef CCPP
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, IPD_Interstitial, Func0d) &
#else
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
#endif
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
#ifdef CCPP
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_Interstitial, IPD_func0d=Func0d)
#else
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
#endif
      enddo
#endif
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP1   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "stochastic physics driver"

      ! DH*
#ifdef DHDEBUG
      write(0,*) "Calling MY_DIAGTOSCREEN before stochastics step"
      call MY_DIAGTOSCREEN()
#endif
      ! *DH

!--- execute the IPD atmospheric physics step2 subcomponent (stochastic physics driver)

      call mpp_clock_begin(physClock)
#ifdef CCPP
      call CCPP_step (step="stochastics", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to IPD-CCPP stochastics step failed')
#else
      Func0d => physics_step2
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
      enddo
#endif
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP2   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif
      call getiauforcing(IPD_Control,IAU_data)
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "end of radiation and physics step"
    endif

    ! DH*
#ifdef DHDEBUG
    write(0,*) "Calling MY_DIAGTOSCREEN after stochastics step"
    call MY_DIAGTOSCREEN()
#endif
    ! *DH

#ifdef CCPP
    ! Update flag for first time step of time integration
    IPD_Control%first_time_step = .false.
#endif
!-----------------------------------------------------------------------
 end subroutine update_atmos_radiation_physics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

#ifdef CCPP
#ifdef OPENMP
  use omp_lib
#endif
  use fv_mp_mod, only: commglobal
  use mpp_mod, only: mpp_npes
#endif

  type (atmos_data_type), intent(inout) :: Atmos
  type (time_type), intent(in) :: Time_init, Time, Time_step
!--- local variables ---
  integer :: unit, ntdiag, ntfamily, i, j, k
  integer :: mlon, mlat, nlon, nlat, nlev, sec, dt
  integer :: ierr, io, logunit
  integer :: idx, tile_num
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: blk, ibs, ibe, jbs, jbe
  real(kind=IPD_kind_phys) :: dt_phys
  real, allocatable    :: q(:,:,:,:), p_half(:,:,:)
  character(len=80)    :: control
  character(len=64)    :: filename, filename2, pelist_name
  character(len=132)   :: text
  logical              :: p_hydro, hydro, fexist
  logical, save        :: block_message = .true.
  type(IPD_init_type)  :: Init_parm
  integer              :: bdat(8), cdat(8)
  integer              :: ntracers, maxhf, maxh
  character(len=32), allocatable, target :: tracer_names(:)
#ifdef CCPP
  integer :: nthrds
#endif

!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
   call get_time (Atmos % Time_step, sec)
   dt_phys = real(sec)      ! integer seconds

   logunit = stdlog()

!-----------------------------------------------------------------------
! initialize atmospheric model -----

#ifndef CCPP
!---------- initialize atmospheric dynamics -------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area)
#endif

   IF ( file_exist('input.nml')) THEN
#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=atmos_model_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_model_nml')
#else
      unit = open_namelist_file ( )
      ierr=1
      do while (ierr /= 0)
         read  (unit, nml=atmos_model_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'atmos_model_nml')
      enddo
 10     call close_file (unit)
#endif
   endif

#ifdef CCPP
!---------- initialize atmospheric dynamics after reading the namelist -------
!---------- (need name of CCPP suite definition file from input.nml) ---------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area)
#endif

!-----------------------------------------------------------------------
   call atmosphere_resolution (nlon, nlat, global=.false.)
   call atmosphere_resolution (mlon, mlat, global=.true.)
   call alloc_atmos_data_type (nlon, nlat, Atmos)
   call atmosphere_domain (Atmos%domain, Atmos%layout, Atmos%regional, Atmos%nested, Atmos%pelist)
   call atmosphere_diag_axes (Atmos%axes)
   call atmosphere_etalvls (Atmos%ak, Atmos%bk, flip=flip_vc)
   call atmosphere_grid_bdry (Atmos%lon_bnd, Atmos%lat_bnd, global=.false.)
   call atmosphere_grid_ctr (Atmos%lon, Atmos%lat)
   call atmosphere_hgt (Atmos%layer_hgt, 'layer', relative=.false., flip=flip_vc)
   call atmosphere_hgt (Atmos%level_hgt, 'level', relative=.false., flip=flip_vc)

   Atmos%mlon = mlon
   Atmos%mlat = mlat
!-----------------------------------------------------------------------
!--- before going any further check definitions for 'blocks'
!-----------------------------------------------------------------------
   call atmosphere_control_data (isc, iec, jsc, jec, nlev, p_hydro, hydro, tile_num)
   call define_blocks_packed ('atmos_model', Atm_block, isc, iec, jsc, jec, nlev, &
                              blocksize, block_message)
   
   allocate(DYCORE_Data(Atm_block%nblks))
   allocate(IPD_Data(Atm_block%nblks))
#ifdef CCPP
#ifdef OPENMP
   nthrds = omp_get_max_threads()
#else
   nthrds = 1
#endif

   ! This logic deals with non-uniform block sizes for CCPP.
   ! When non-uniform block sizes are used, it is required
   ! that only the last block has a different (smaller)
   ! size than all other blocks. This is the standard in
   ! FV3. If this is the case, set non_uniform_blocks (a
   ! variable imported from CCPP_driver) to .true. and
   ! allocate nthreads+1 elements of the interstitial array.
   ! The extra element will be used by the thread that
   ! runs over the last, smaller block.
   if (minval(Atm_block%blksz)==maxval(Atm_block%blksz)) then
      non_uniform_blocks = .false.
      allocate(IPD_Interstitial(nthrds))
   else if (all(minloc(Atm_block%blksz)==(/size(Atm_block%blksz)/))) then
      non_uniform_blocks = .true.
      allocate(IPD_Interstitial(nthrds+1))
   else
      call mpp_error(FATAL, 'For non-uniform blocksizes, only the last element ' // &
                            'in Atm_block%blksz can be different from the others')
   end if

#endif

!--- update IPD_Control%jdat(8)
   bdat(:) = 0
   call get_date (Time_init, bdat(1), bdat(2), bdat(3),  &
                             bdat(5), bdat(6), bdat(7))
   cdat(:) = 0
   call get_date (Time,      cdat(1), cdat(2), cdat(3),  &
                             cdat(5), cdat(6), cdat(7))
   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   allocate (tracer_names(ntracers))
   do i = 1, ntracers
     call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
   enddo
!--- setup IPD Init_parm
   Init_parm%me              =  mpp_pe()
   Init_parm%master          =  mpp_root_pe()
   Init_parm%tile_num        =  tile_num
   Init_parm%isc             =  isc
   Init_parm%jsc             =  jsc
   Init_parm%nx              =  nlon
   Init_parm%ny              =  nlat
   Init_parm%levs            =  nlev
   Init_parm%cnx             =  mlon
   Init_parm%cny             =  mlat
   Init_parm%gnx             =  Init_parm%cnx*4
   Init_parm%gny             =  Init_parm%cny*2
   Init_parm%nlunit          =  9999
   Init_parm%logunit         =  logunit
   Init_parm%bdat(:)         =  bdat(:)
   Init_parm%cdat(:)         =  cdat(:)
   Init_parm%dt_dycore       =  dt_phys
   Init_parm%dt_phys         =  dt_phys
   Init_parm%blksz           => Atm_block%blksz
   Init_parm%ak              => Atmos%ak
   Init_parm%bk              => Atmos%bk
   Init_parm%xlon            => Atmos%lon
   Init_parm%xlat            => Atmos%lat
   Init_parm%area            => Atmos%area
   Init_parm%tracer_names    => tracer_names
#ifdef CCPP
   Init_parm%restart         = Atm(mytile)%flagstruct%warm_start
   Init_parm%hydrostatic     = Atm(mytile)%flagstruct%hydrostatic
#endif

#ifdef INTERNAL_FILE_NML
   Init_parm%input_nml_file  => input_nml_file
   Init_parm%fn_nml='using internal file'
#else
   pelist_name=mpp_get_current_pelist_name()
   Init_parm%fn_nml='input_'//trim(pelist_name)//'.nml'
   inquire(FILE=Init_parm%fn_nml, EXIST=fexist)
   if (.not. fexist ) then
      Init_parm%fn_nml='input.nml'
   endif
#endif

#ifdef CCPP
   call IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, &
                        IPD_Interstitial, commglobal, mpp_npes(), Init_parm)
#else
   call IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Init_parm)
#endif

#ifdef CCPP
   ! Initialize the CCPP framework
   call CCPP_step (step="init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to IPD-CCPP init step failed')
   ! Doing the init here requires logic in thompson aerosol init if no aerosol
   ! profiles are specified and internal profiles are calculated, because these
   ! require temperature/geopotential etc which are not yet set. Sim. for RUC LSM.
   call CCPP_step (step="physics_init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to IPD-CCPP physics_init step failed')
#endif

   Atmos%Diag => IPD_Diag

   Atm(mytile)%flagstruct%do_skeb = IPD_Control%do_skeb

!  initialize the IAU module
   call iau_initialize (IPD_Control,IAU_data,Init_parm)

   Init_parm%blksz           => null()
   Init_parm%ak              => null()
   Init_parm%bk              => null()
   Init_parm%xlon            => null()
   Init_parm%xlat            => null()
   Init_parm%area            => null()
   Init_parm%tracer_names    => null()
   deallocate (tracer_names)

   !--- update tracers in FV3 with any initialized during the physics/radiation init phase
!rab   call atmosphere_tracer_postinit (IPD_Data, Atm_block)

   call atmosphere_nggps_diag (Time, init=.true.)
   call FV3GFS_diag_register (IPD_Diag, Time, Atm_block, IPD_Control, Atmos%lon, Atmos%lat, Atmos%axes)
#ifdef CCPP
   call FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, IPD_Control, Atmos%domain, Atm(mytile)%flagstruct%warm_start)
#else
   call FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, IPD_Control, Atmos%domain)
#endif

   !--- set the initial diagnostic timestamp
   diag_time = Time 
   if (output_1st_tstep_rst) then
     diag_time = Time - real_to_time_type(mod(int((first_kdt - 1)*dt_phys/3600.),6)*3600.0)
   endif

   !---- print version number to logfile ----

   call write_version_number ( version, tagname )
   !--- write the namelist to a log file
   if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog( )
      write (unit, nml=atmos_model_nml)
      call close_file (unit)
   endif

   !--- get fdiag
#ifdef GFS_PHYS
!--- check fdiag to see if it is an interval or a list
   if (nint(fdiag(2)) == 0) then
     if (fhmaxhf > 0) then
       maxhf = fhmaxhf / fhouthf
       maxh  = maxhf + (fhmax-fhmaxhf) / fhout
       fdiag(1) = fhouthf
       do i=2,maxhf
        fdiag(i) = fdiag(i-1) + fhouthf
       enddo
       do i=maxhf+1,maxh
         fdiag(i) = fdiag(i-1) + fhout
       enddo
     else
       maxh  = fhmax / fhout
       do i = 2, maxh
         fdiag(i) = fdiag(i-1) + fhout
       enddo
     endif
   endif
   if (mpp_pe() == mpp_root_pe()) write(6,*) "---fdiag",fdiag(1:40)
#endif

   setupClock = mpp_clock_id( 'GFS Step Setup        ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   radClock   = mpp_clock_id( 'GFS Radiation         ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   physClock  = mpp_clock_id( 'GFS Physics           ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   getClock   = mpp_clock_id( 'Dynamics get state    ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   updClock   = mpp_clock_id( 'Dynamics update state ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   if (sync) then
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default+MPP_CLOCK_SYNC, grain=CLOCK_COMPONENT )
   else
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   endif

#ifdef CCPP
   ! Set flag for first time step of time integration
   IPD_Control%first_time_step = .true.
#endif
!-----------------------------------------------------------------------
end subroutine atmos_model_init
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_dynamics"
!
! <OVERVIEW>
subroutine update_atmos_model_dynamics (Atmos)
! run the atmospheric dynamics to advect the properties
  type (atmos_data_type), intent(in) :: Atmos

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call atmosphere_dynamics (Atmos%Time)
    call mpp_clock_end(fv3Clock)

end subroutine update_atmos_model_dynamics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_exchange_phase_1"
!
! <OVERVIEW>
!   Perform data exchange with coupled components in run phase 1
! </OVERVIEW>
!
! <DESCRIPTION>
!  This subroutine currently exports atmospheric fields and tracers
!  to the chemistry component during the model's run phase 1, i.e.
!  before chemistry is run.
! </DESCRIPTION>

subroutine atmos_model_exchange_phase_1 (Atmos, rc)

  use ESMF

  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    !--- if coupled, exchange coupled fields
    if( IPD_Control%cplchm ) then
      ! -- export fields to chemistry
      call update_atmos_chemistry('export', rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
    endif

 end subroutine atmos_model_exchange_phase_1
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_exchange_phase_2"
!
! <OVERVIEW>
!   Perform data exchange with coupled components in run phase 2
! </OVERVIEW>
!
! <DESCRIPTION>
!  This subroutine currently imports fields updated by the coupled
!  chemistry component back into the atmospheric model during run
!  phase 2.
! </DESCRIPTION>

subroutine atmos_model_exchange_phase_2 (Atmos, rc)

  use ESMF

  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    !--- if coupled, exchange coupled fields
    if( IPD_Control%cplchm ) then
      ! -- import fields from chemistry
      call update_atmos_chemistry('import', rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
    endif

 end subroutine atmos_model_exchange_phase_2
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_state"
!
! <OVERVIEW>
subroutine update_atmos_model_state (Atmos)
! to update the model state after all concurrency is completed
  type (atmos_data_type), intent(inout) :: Atmos
!--- local variables
  integer :: isec,seconds
  integer :: rc
  real(kind=IPD_kind_phys) :: time_int, time_intfull
!
    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call mpp_clock_begin(updClock)
    call atmosphere_state_update (Atmos%Time, IPD_Data, IAU_Data, Atm_block, flip_vc)
    call mpp_clock_end(updClock)
    call mpp_clock_end(fv3Clock)

    if (chksum_debug) then
      if (mpp_pe() == mpp_root_pe()) print *,'UPDATE STATE    ', IPD_Control%kdt, IPD_Control%fhour
      if (mpp_pe() == mpp_root_pe()) print *,'in UPDATE STATE    ', size(IPD_Data(1)%SfcProp%tsfc),'nblks=',Atm_block%nblks
      call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
    endif

    !--- advance time ---
    Atmos % Time = Atmos % Time + Atmos % Time_step

    call get_time (Atmos%Time - diag_time, isec)
    call get_time (Atmos%Time - Atmos%Time_init, seconds)
    call atmosphere_nggps_diag(Atmos%Time,ltavg=.true.,avg_max_length=avg_max_length)
    if (ANY(nint(fdiag(:)*3600.0) == seconds) .or. (IPD_Control%kdt == first_kdt) ) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---isec,seconds",isec,seconds
      time_int = real(isec)
      time_intfull = real(seconds)
      if (mpp_pe() == mpp_root_pe()) write(6,*) ' gfs diags time since last bucket empty: ',time_int/3600.,'hrs'
      call atmosphere_nggps_diag(Atmos%Time)
      call FV3GFS_diag_output(Atmos%Time, IPD_DIag, Atm_block, IPD_Control%nx, IPD_Control%ny, &
                            IPD_Control%levs, 1, 1, 1.d0, time_int, time_intfull,              &
                            IPD_Control%fhswr, IPD_Control%fhlwr)
      if (mod(isec,3600*nint(IPD_Control%fhzero)) == 0) diag_time = Atmos%Time
      call diag_send_complete_instant (Atmos%Time)
    endif

    !--- this may not be necessary once write_component is fully implemented
    !!!call diag_send_complete_extra (Atmos%Time)

    !--- get bottom layer data from dynamical core for coupling
    call atmosphere_get_bottom_layer (Atm_block, DYCORE_Data) 

    !if in coupled mode, set up coupled fields
    if (IPD_Control%cplflx .or. IPD_Control%cplwav) then
      if (mpp_pe() == mpp_root_pe()) print *,'COUPLING: IPD layer'
!jw       call setup_exportdata(IPD_Control, IPD_Data, Atm_block)
      call setup_exportdata(rc)
    endif

 end subroutine update_atmos_model_state
! </SUBROUTINE>



!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_end (Atmos)
  type (atmos_data_type), intent(inout) :: Atmos
!---local variables
  integer :: idx
#ifdef CCPP
  integer :: ierr
#endif

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----
    call atmosphere_end (Atmos % Time, Atmos%grid)
    call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
                               IPD_Control, Atmos%domain)

#ifdef CCPP
!   Fast physics (from dynamics) are finalized in atmosphere_end above;
!   standard/slow physics (from IPD) are finalized in CCPP_step 'finalize'.
!   The CCPP framework for all cdata structures is finalized in CCPP_step 'finalize'.
    call CCPP_step (step="finalize", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to IPD-CCPP finalize step failed')
#endif

end subroutine atmos_model_end

! </SUBROUTINE>
!#######################################################################
! <SUBROUTINE NAME="atmos_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine atmos_model_restart(Atmos, timestamp)
  type (atmos_data_type),   intent(inout) :: Atmos
  character(len=*),  intent(in)           :: timestamp

    call atmosphere_restart(timestamp)
    call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
                               IPD_Control, Atmos%domain, timestamp)

end subroutine atmos_model_restart
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="get_atmos_model_ungridded_dim">
!
! <DESCRIPTION>
!  Retrieve ungridded dimensions of atmospheric model arrays
! </DESCRIPTION>

subroutine get_atmos_model_ungridded_dim(nlev, ntracers, nsoillev)

  integer, optional, intent(out) :: nlev, ntracers, nsoillev

  if (present(nlev))     nlev = Atm_block%npz
  if (present(ntracers)) call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
  if (present(nsoillev)) then
    nsoillev = 0
    if (allocated(IPD_Data)) then
      if (associated(IPD_Data(1)%Sfcprop%slc)) &
        nsoillev = size(IPD_Data(1)%Sfcprop%slc, 2)
    end if
  end if

end subroutine get_atmos_model_ungridded_dim
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_chemistry">
! <DESCRIPTION>
!  Populate exported chemistry fields with current atmospheric state
!  data (state='export'). Update tracer concentrations for atmospheric
!  chemistry with values from chemistry component (state='import').
!  Fields should be exported/imported from/to the atmospheric state
!  after physics calculations.
!
!  NOTE: It is assumed that all the chemical tracers follow the standard
!  atmospheric tracers, which end with ozone. The order of the chemical
!  tracers must match their order in the chemistry component.
!
!  Requires:
!         IPD_Data
!         Atm_block
! </DESCRIPTION>
subroutine update_atmos_chemistry(state, rc)

  use ESMF
  use module_cplfields,   only: cplFieldGet

  character(len=*),  intent(in)  :: state
  integer, optional, intent(out) :: rc

  !--- local variables
  integer :: localrc
  integer :: ni, nj, nk, nt, ntoz
  integer :: nb, ix, i, j, k, it
  integer :: ib, jb

  real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: prsl, phil, &
                                                     prsi, phii, &
                                                     temp, &
                                                     ua, va, vvl, &
                                                     dkt, slc
  real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: q

  real(ESMF_KIND_R8), dimension(:,:), pointer :: hpbl, area, stype, rainc, &
    uustar, rain, sfcdsw, slmsk, tsfc, shfsfc, snowd, vtype, vfrac, zorl

  logical, parameter :: diag = .true.

  ! -- begin
  if (present(rc)) rc = ESMF_SUCCESS

  ni  = Atm_block%iec - Atm_block%isc + 1
  nj  = Atm_block%jec - Atm_block%jsc + 1
  nk  = Atm_block%npz
  call get_number_tracers(MODEL_ATMOS, num_tracers=nt)

  select case (trim(state))
    case ('import')
      !--- retrieve references to allocated memory for each field
      call cplFieldGet(state,'inst_tracer_mass_frac', &
        farrayPtr4d=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      !--- tracers quantities
      !--- locate the end location of standard atmospheric tracers, marked by ozone
      ntoz = get_tracer_index(MODEL_ATMOS, 'o3mr')

      do it = ntoz + 1, nt
!$OMP parallel do default (none) &
!$OMP             shared  (it, nk, nj, ni, Atm_block, IPD_Data, q)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
        do k = 1, nk
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              IPD_Data(nb)%Stateout%gq0(ix,k,it) = q(i,j,k,it)
            enddo
          enddo
        enddo
      enddo

      if (diag) then
        write(6,'("update_atmos: ",a,": qgrs - min/max/avg",3g16.6)') &
          trim(state), minval(q), maxval(q), sum(q)/size(q)
      end if

    case ('export')
      !--- retrieve references to allocated memory for each field
      call cplFieldGet(state,'inst_pres_interface', farrayPtr3d=prsi, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_pres_levels', &
        farrayPtr3d=prsl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_geop_interface', farrayPtr3d=phii, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_geop_levels', &
        farrayPtr3d=phil, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_temp_levels', farrayPtr3d=temp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_zonal_wind_levels', farrayPtr3d=ua, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_merid_wind_levels', farrayPtr3d=va, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_omega_levels', farrayPtr3d=vvl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_tracer_mass_frac', &
        farrayPtr4d=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_soil_moisture_content', farrayPtr3d=slc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'soil_type', farrayPtr2d=stype, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_pbl_height', &
        farrayPtr2d=hpbl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'surface_cell_area', farrayPtr2d=area, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_convective_rainfall_amount', &
        farrayPtr2d=rainc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_exchange_coefficient_heat_levels', &
        farrayPtr3d=dkt, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_friction_velocity', farrayPtr2d=uustar, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_rainfall_amount', farrayPtr2d=rain, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_down_sw_flx', &
        farrayPtr2d=sfcdsw, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_land_sea_mask', farrayPtr2d=slmsk, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_temp_height_surface', farrayPtr2d=tsfc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_up_sensi_heat_flx', &
        farrayPtr2d=shfsfc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_lwe_snow_thickness', &
        farrayPtr2d=snowd, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'vegetation_type', farrayPtr2d=vtype, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_vegetation_area_frac', &
        farrayPtr2d=vfrac, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_surface_roughness', farrayPtr2d=zorl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      !--- handle all three-dimensional variables
!$OMP parallel do default (none) &
!$OMP             shared  (nk, nj, ni, Atm_block, IPD_Data, prsi, phii, prsl, phil, temp, ua, va, vvl, dkt)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
      do k = 1, nk
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            !--- interface values
            prsi(i,j,k) = IPD_Data(nb)%Statein%prsi(ix,k)
            phii(i,j,k) = IPD_Data(nb)%Statein%phii(ix,k)
            !--- layer values
            prsl(i,j,k) = IPD_Data(nb)%Statein%prsl(ix,k)
            phil(i,j,k) = IPD_Data(nb)%Statein%phil(ix,k)
            temp(i,j,k) = IPD_Data(nb)%Stateout%gt0(ix,k)
            ua  (i,j,k) = IPD_Data(nb)%Stateout%gu0(ix,k)
            va  (i,j,k) = IPD_Data(nb)%Stateout%gv0(ix,k)
            vvl (i,j,k) = IPD_Data(nb)%Statein%vvl (ix,k)
            dkt (i,j,k) = IPD_Data(nb)%Coupling%dkt(ix,k)
          enddo
        enddo
      enddo

      !--- top interface values
      k = nk+1
      do j = 1, nj
        jb = j + Atm_block%jsc - 1
        do i = 1, ni
          ib = i + Atm_block%isc - 1
          nb = Atm_block%blkno(ib,jb)
          ix = Atm_block%ixp(ib,jb)
          prsi(i,j,k) = IPD_Data(nb)%Statein%prsi(ix,k)
          phii(i,j,k) = IPD_Data(nb)%Statein%phii(ix,k)
        enddo
      enddo

      !--- tracers quantities
      do it = 1, nt
!$OMP parallel do default (none) &
!$OMP             shared  (it, nk, nj, ni, Atm_block, IPD_Data, q)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
        do k = 1, nk
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              q(i,j,k,it) = IPD_Data(nb)%Stateout%gq0(ix,k,it)
            enddo
          enddo
        enddo
      enddo

!$OMP parallel do default (none) &
!$OMP             shared  (nj, ni, Atm_block, IPD_Data, &
!$OMP                      hpbl, area, stype, rainc, rain, uustar, sfcdsw, &
!$OMP                      slmsk, snowd, tsfc, shfsfc, vtype, vfrac, zorl, slc) &
!$OMP             private (j, jb, i, ib, nb, ix)
      do j = 1, nj
        jb = j + Atm_block%jsc - 1
        do i = 1, ni
          ib = i + Atm_block%isc - 1
          nb = Atm_block%blkno(ib,jb)
          ix = Atm_block%ixp(ib,jb)
          hpbl(i,j)    = IPD_Data(nb)%IntDiag%hpbl(ix)
          area(i,j)    = IPD_Data(nb)%Grid%area(ix)
          stype(i,j)   = IPD_Data(nb)%Sfcprop%stype(ix)
          rainc(i,j)   = IPD_Data(nb)%Coupling%rainc_cpl(ix)
          rain(i,j)    = IPD_Data(nb)%Coupling%rain_cpl(ix)
          uustar(i,j)  = IPD_Data(nb)%Sfcprop%uustar(ix)
          sfcdsw(i,j)  = IPD_Data(nb)%Coupling%sfcdsw(ix)
          slmsk(i,j)   = IPD_Data(nb)%Sfcprop%slmsk(ix)
          snowd(i,j)   = IPD_Data(nb)%Sfcprop%snowd(ix)
          tsfc(i,j)    = IPD_Data(nb)%Sfcprop%tsfc(ix)
          shfsfc(i,j)  = IPD_Data(nb)%Coupling%ushfsfci(ix)
          vtype(i,j)   = IPD_Data(nb)%Sfcprop%vtype(ix)
          vfrac(i,j)   = IPD_Data(nb)%Sfcprop%vfrac(ix)
          zorl(i,j)    = IPD_Data(nb)%Sfcprop%zorl(ix)
          slc(i,j,:)   = IPD_Data(nb)%Sfcprop%slc(ix,:)
        enddo
      enddo

      if (diag) then
        ! -- diagnostics
        write(6,'("update_atmos: prsi - min/max/avg",3g16.6)') minval(prsi), maxval(prsi), sum(prsi)/size(prsi)
        write(6,'("update_atmos: phii - min/max/avg",3g16.6)') minval(phii), maxval(phii), sum(phii)/size(phii)
        write(6,'("update_atmos: prsl - min/max/avg",3g16.6)') minval(prsl), maxval(prsl), sum(prsl)/size(prsl)
        write(6,'("update_atmos: phil - min/max/avg",3g16.6)') minval(phil), maxval(phil), sum(phil)/size(phil)
        write(6,'("update_atmos: tgrs - min/max/avg",3g16.6)') minval(temp), maxval(temp), sum(temp)/size(temp)
        write(6,'("update_atmos: ugrs - min/max/avg",3g16.6)') minval(ua), maxval(ua), sum(ua)/size(ua)
        write(6,'("update_atmos: vgrs - min/max/avg",3g16.6)') minval(va), maxval(va), sum(va)/size(va)
        write(6,'("update_atmos: vvl  - min/max/avg",3g16.6)') minval(vvl), maxval(vvl), sum(vvl)/size(vvl)
        write(6,'("update_atmos: qgrs - min/max/avg",3g16.6)') minval(q), maxval(q), sum(q)/size(q)

        write(6,'("update_atmos: hpbl - min/max/avg",3g16.6)') minval(hpbl), maxval(hpbl), sum(hpbl)/size(hpbl)
        write(6,'("update_atmos: rainc - min/max/avg",3g16.6)') minval(rainc), maxval(rainc), sum(rainc)/size(rainc)
        write(6,'("update_atmos: rain - min/max/avg",3g16.6)') minval(rain), maxval(rain), sum(rain)/size(rain)
        write(6,'("update_atmos: shfsfc - min/max/avg",3g16.6)') minval(shfsfc), maxval(shfsfc), sum(shfsfc)/size(shfsfc)
        write(6,'("update_atmos: sfcdsw - min/max/avg",3g16.6)') minval(sfcdsw), maxval(sfcdsw), sum(sfcdsw)/size(sfcdsw)
        write(6,'("update_atmos: slmsk - min/max/avg",3g16.6)') minval(slmsk), maxval(slmsk), sum(slmsk)/size(slmsk)
        write(6,'("update_atmos: snowd - min/max/avg",3g16.6)') minval(snowd), maxval(snowd), sum(snowd)/size(snowd)
        write(6,'("update_atmos: tsfc - min/max/avg",3g16.6)') minval(tsfc), maxval(tsfc), sum(tsfc)/size(tsfc)
        write(6,'("update_atmos: vtype - min/max/avg",3g16.6)') minval(vtype), maxval(vtype), sum(vtype)/size(vtype)
        write(6,'("update_atmos: vfrac - min/max/avg",3g16.6)') minval(vfrac), maxval(vfrac), sum(vfrac)/size(vfrac)
        write(6,'("update_atmos: area - min/max/avg",3g16.6)') minval(area), maxval(area), sum(area)/size(area)
        write(6,'("update_atmos: stype - min/max/avg",3g16.6)') minval(stype), maxval(stype), sum(stype)/size(stype)
        write(6,'("update_atmos: zorl - min/max/avg",3g16.6)') minval(zorl), maxval(zorl), sum(zorl)/size(zorl)
        write(6,'("update_atmos: slc - min/max/avg",3g16.6)') minval(slc), maxval(slc), sum(slc)/size(slc)
      end if

    case default
      ! -- do nothing
  end select

end subroutine update_atmos_chemistry
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_data_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the atmos_data_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the atmos_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, atm)
! </TEMPLATE>

! <IN NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields in the atmos_data_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine atmos_data_type_chksum(id, timestep, atm)
type(atmos_data_type), intent(in) :: atm 
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    integer :: n, outunit

100 format("CHECKSUM::",A32," = ",Z20)
101 format("CHECKSUM::",A16,a,'%',a," = ",Z20)

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(Atmos_data_type):: ', id, timestep
  write(outunit,100) ' atm%lon_bnd                ', mpp_chksum(atm%lon_bnd               )
  write(outunit,100) ' atm%lat_bnd                ', mpp_chksum(atm%lat_bnd               )
  write(outunit,100) ' atm%lon                    ', mpp_chksum(atm%lon                   )
  write(outunit,100) ' atm%lat                    ', mpp_chksum(atm%lat                   )

end subroutine atmos_data_type_chksum

! </SUBROUTINE>

  subroutine alloc_atmos_data_type (nlon, nlat, Atmos)
   integer, intent(in) :: nlon, nlat
   type(atmos_data_type), intent(inout) :: Atmos
    allocate ( Atmos % lon_bnd  (nlon+1,nlat+1), &
               Atmos % lat_bnd  (nlon+1,nlat+1), &
               Atmos % lon      (nlon,nlat),     &
               Atmos % lat      (nlon,nlat)      )

  end subroutine alloc_atmos_data_type

  subroutine dealloc_atmos_data_type (Atmos)
   type(atmos_data_type), intent(inout) :: Atmos
    deallocate (Atmos%lon_bnd, &
                Atmos%lat_bnd, &
                Atmos%lon,     &
                Atmos%lat      )
  end subroutine dealloc_atmos_data_type

  subroutine assign_importdata(rc)

    use module_cplfields,  only: importFields, nImportFields, QueryFieldList, &
                                 ImportFieldsList, importFieldsValid
    use ESMF
!
    implicit none
    integer, intent(out) :: rc

    !--- local variables
    integer :: n, j, i, ix, nb, isc, iec, jsc, jec, dimCount, findex
    character(len=128) :: impfield_name, fldname
    type(ESMF_TypeKind_Flag)                           :: datatype
    real(kind=ESMF_KIND_R4), dimension(:,:), pointer   :: datar42d
    real(kind=ESMF_KIND_R8), dimension(:,:), pointer   :: datar82d
    real(kind=IPD_kind_phys), dimension(:,:), pointer  :: datar8
    logical found, lcpl_fice
!
!------------------------------------------------------------------------------
!
    ! set up local dimension
    rc=-999
    isc = IPD_control%isc
    iec = IPD_control%isc+IPD_control%nx-1
    jsc = IPD_control%jsc
    jec = IPD_control%jsc+IPD_control%ny-1
    lcpl_fice = .false.

    allocate(datar8(isc:iec,jsc:jec))
    if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,dim=',isc,iec,jsc,jec
    if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,IPD_Data, size', size(IPD_Data)
    if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,tsfc, size', size(IPD_Data(1)%sfcprop%tsfc)

    do n=1,nImportFields

      ! Each import field is only available if it was connected in the
      ! import state.
      found = .false.
      if (ESMF_FieldIsCreated(importFields(n))) then

        ! put the data from local cubed sphere grid to column grid for phys
        datar8 = -99999.0
        call ESMF_FieldGet(importFields(n), dimCount=dimCount ,typekind=datatype, &
          name=impfield_name, rc=rc)

        if ( dimCount == 2) then
          if ( datatype == ESMF_TYPEKIND_R8) then
            call ESMF_FieldGet(importFields(n),farrayPtr=datar82d,localDE=0, rc=rc)
            datar8=datar82d
            if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplIMP,atmos gets ',trim(impfield_name),' datar8=',maxval(datar8),minval(datar8), &
               datar8(isc,jsc)
            found = .true.
! gfs physics runs with r8
!          else
!            call ESMF_FieldGet(importFields(n),farrayPtr=datar42d,localDE=0,
!            rc=rc)
!            datar8=datar42d
          endif
        endif
!
        ! get sea land mask: in order to update the coupling fields over the ocean/ice
!        fldname = 'land_mask'
!        findex = QueryFieldList(ImportFieldsList,fldname)
!        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
!          if (trim(impfield_name) == trim(fldname) .and. found) then
!            do j=jsc,jec
!            do i=isc,iec
!              nb = Atm_block%blkno(i,j)
!              ix = Atm_block%ixp(i,j)
!              IPD_Data(nb)%Coupling%slimskin_cpl(ix) = datar8(i,j)
!            enddo
!            enddo
!            if( mpp_pe()==mpp_root_pe()) print *,'get land mask from mediator'
!          endif
!        endif

        ! get sea ice surface temperature 
        fldname = 'sea_ice_surface_temperature'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              IPD_Data(nb)%Coupling%tisfcin_cpl(ix) = datar8(i,j)
            enddo
            enddo
          endif
        endif

        ! get sst:  sst needs to be adjusted by land sea mask before passing to
        ! fv3
        fldname = 'sea_surface_temperature'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
!
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
!if it is ocean or ice get sst from mediator
              if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
                IPD_Data(nb)%Coupling%tseain_cpl(ix) = datar8(i,j)
                IPD_Data(nb)%Sfcprop%tsfc(ix) = datar8(i,j)
              endif
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get sst from mediator'
          endif
        endif

        ! get sea ice fraction:  fice or sea ice concentration from the mediator
        fldname = 'ice_fraction'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            lcpl_fice = .true.
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              IPD_Data(nb)%Coupling%ficein_cpl(ix) = 0.
              IPD_Data(nb)%Coupling%slimskin_cpl(ix) = 0.
!if it is ocean or ice get sst from mediator
              if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
                if( datar8(i,j) > 0.15 .and. IPD_Data(nb)%Sfcprop%lakemsk(ix) /= 1 ) then
                  IPD_Data(nb)%Coupling%ficein_cpl(ix) = datar8(i,j)
                  IPD_Data(nb)%Sfcprop%slmsk(ix) = 2.0
                  IPD_Data(nb)%Coupling%slimskin_cpl(ix) = 4.
                endif
              endif
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get fice from mediator'
          endif
        endif

        ! get upward LW flux:  for sea ice covered area
        fldname = 'mean_up_lw_flx'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
                IPD_Data(nb)%Coupling%ulwsfcin_cpl(ix) = -datar8(i,j)
              endif
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get lwflx from mediator'
          endif
        endif

        ! get latent heat flux:  for sea ice covered area
        fldname = 'mean_laten_heat_flx'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
                IPD_Data(nb)%Coupling%dqsfcin_cpl(ix) = -datar8(i,j)
              endif
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get laten_heat from mediator'
          endif
        endif

        ! get sensible heat flux:  for sea ice covered area
        fldname = 'mean_sensi_heat_flx'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
                IPD_Data(nb)%Coupling%dtsfcin_cpl(ix) = -datar8(i,j)
              endif
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensi_heat from mediator'
          endif
        endif

        ! get zonal compt of momentum flux:  for sea ice covered area
        fldname = 'mean_zonal_moment_flx'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
                IPD_Data(nb)%Coupling%dusfcin_cpl(ix) = -datar8(i,j)
              endif
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get zonal_moment_flx from mediator'
          endif
        endif

        ! get meridional compt of momentum flux:  for sea ice covered area
        fldname = 'mean_merid_moment_flx'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
                IPD_Data(nb)%Coupling%dvsfcin_cpl(ix) = -datar8(i,j)
              endif
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get merid_moment_flx from mediator'
          endif
        endif

        ! get sea ice volume:  for sea ice covered area
        fldname = 'mean_ice_volume'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              IPD_Data(nb)%Coupling%hicein_cpl(ix) = datar8(i,j)
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug) print *,'fv3 assign_import: get ice_volume  from mediator'
          endif
        endif

        ! get snow volume:  for sea ice covered area
        fldname = 'mean_snow_volume'
        findex = QueryFieldList(ImportFieldsList,fldname)
        if (importFieldsValid(findex) .and. datar8(isc,jsc) > -99999.0) then
          if (trim(impfield_name) == trim(fldname) .and. found) then
            do j=jsc,jec
            do i=isc,iec
              nb = Atm_block%blkno(i,j)
              ix = Atm_block%ixp(i,j)
              IPD_Data(nb)%Coupling%hsnoin_cpl(ix) = datar8(i,j)
            enddo
            enddo
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get snow_volume  from mediator'
          endif
        endif

      endif
    enddo
!
    deallocate(datar8)

! update sea ice related fields:
    if( lcpl_fice ) then
      do j=jsc,jec
      do i=isc,iec
        nb = Atm_block%blkno(i,j)
        ix = Atm_block%ixp(i,j)
!if it is ocean or ice get sst from mediator
        if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or.  IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
           IPD_Data(nb)%Sfcprop%tisfc(ix) = IPD_Data(nb)%Coupling%tisfcin_cpl(ix)
           if( IPD_Data(nb)%Sfcprop%lakemsk(ix) /= 1 ) then
             if( IPD_Data(nb)%Coupling%ficein_cpl(ix) > 0.15 ) then
               IPD_Data(nb)%Sfcprop%fice(ix)  = IPD_Data(nb)%Coupling%ficein_cpl(ix)
               IPD_Data(nb)%Sfcprop%hice(ix)  = IPD_Data(nb)%Coupling%hicein_cpl(ix)
               IPD_Data(nb)%Sfcprop%snowd(ix) = IPD_Data(nb)%Coupling%hsnoin_cpl(ix)
             else
               IPD_Data(nb)%Sfcprop%fice(ix)  = 0.
               IPD_Data(nb)%Sfcprop%hice(ix)  = 0.
               IPD_Data(nb)%Sfcprop%snowd(ix) = 0.
               IPD_Data(nb)%Sfcprop%slmsk(ix) = 0.
             endif
           endif
        endif
      enddo
      enddo
    endif

    rc=0
!
    if (mpp_pe() == mpp_root_pe()) print *,'end of assign_importdata'
  end subroutine assign_importdata

!
  subroutine setup_exportdata (rc)

    use module_cplfields,  only: exportData, nExportFields, exportFieldsList, &
                                 queryFieldList, fillExportFields

    implicit none

!------------------------------------------------------------------------------

    !--- interface variables
    integer, intent(out) :: rc

    !--- local variables
    integer            :: j, i, ix, nb, isc, iec, jsc, jec, idx
    real(IPD_kind_phys)    :: rtime
!
    if (mpp_pe() == mpp_root_pe()) print *,'enter setup_exportdata'

    isc = IPD_control%isc
    iec = IPD_control%isc+IPD_control%nx-1
    jsc = IPD_control%jsc
    jec = IPD_control%jsc+IPD_control%ny-1

    rtime  = 1./IPD_control%dtp
!    print *,'in cplExp,dim=',isc,iec,jsc,jec,'nExportFields=',nExportFields
!    print *,'in cplExp,IPD_Data, size', size(IPD_Data)
!    print *,'in cplExp,u10micpl, size', size(IPD_Data(1)%coupling%u10mi_cpl)

    if(.not.allocated(exportData)) then
      allocate(exportData(isc:iec,jsc:jec,nExportFields))
    endif

    ! set cpl fields to export Data

    if (IPD_Control%cplflx .or. IPD_Control%cplwav) then 
    ! Instantaneous u wind (m/s) 10 m above ground
    idx = queryfieldlist(exportFieldsList,'inst_zonal_wind_height10m')
    if (idx > 0 ) then
      if (mpp_pe() == mpp_root_pe() .and. debug) print *,'cpl, in get u10mi_cpl'
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%u10mi_cpl(ix)
        enddo
      enddo
    endif

    ! Instantaneous v wind (m/s) 10 m above ground
    idx = queryfieldlist(exportFieldsList,'inst_merid_wind_height10m')
    if (idx > 0 ) then
      if (mpp_pe() == mpp_root_pe() .and. debug) print *,'cpl, in get v10mi_cpl'
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%v10mi_cpl(ix)
        enddo
      enddo
      if (mpp_pe() == mpp_root_pe() .and. debug) print *,'cpl, get v10mi_cpl, exportData=',exportData(isc,jsc,idx),'idx=',idx
    endif

    endif !if cplflx or cplwav 

    if (IPD_Control%cplflx) then
    ! MEAN Zonal compt of momentum flux (N/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_zonal_moment_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dusfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Merid compt of momentum flux (N/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_merid_moment_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Sensible heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_sensi_heat_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dtsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Latent heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_laten_heat_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dqsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Downward LW heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_lw_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dlwsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Downward SW heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dswsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN precipitation rate (kg/m2) ?????? checking unit ??????
    idx = queryfieldlist(exportFieldsList,'mean_prec_rate')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%rain_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! Instataneous Zonal compt of momentum flux (N/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_zonal_moment_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dusfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Merid compt of momentum flux (N/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_merid_moment_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Sensible heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_sensi_heat_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dtsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Latent heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_laten_heat_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dqsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Downward long wave radiation flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_lw_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dlwsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Downward solar radiation flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dswsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Temperature (K) 2 m above ground
    idx = queryfieldlist(exportFieldsList,'inst_temp_height2m')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%t2mi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Specific humidity (kg/kg) 2 m above ground
    idx = queryfieldlist(exportFieldsList,'inst_spec_humid_height2m')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%q2mi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Temperature (K) at surface
    idx = queryfieldlist(exportFieldsList,'inst_temp_height_surface')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%tsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Pressure (Pa) land and sea surface
    idx = queryfieldlist(exportFieldsList,'inst_pres_height_surface')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%psurfi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Surface height (m)
    idx = queryfieldlist(exportFieldsList,'inst_surface_height')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%oro_cpl(ix)
        enddo
      enddo
    endif

    ! MEAN NET long wave radiation flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_lw_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nlwsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN NET solar radiation flux over the ocean (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nswsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! Instataneous NET long wave radiation flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_lw_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nlwsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous NET solar radiation flux over the ocean (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nswsfci_cpl(ix)
        enddo
      enddo
    endif

    ! MEAN sfc downward nir direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_ir_dir_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dnirbm_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN sfc downward nir diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_ir_dif_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dnirdf_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN sfc downward uv+vis direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_vis_dir_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvisbm_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN sfc downward uv+vis diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_vis_dif_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvisdf_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! Instataneous sfc downward nir direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_ir_dir_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dnirbmi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous sfc downward nir diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_ir_dif_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dnirdfi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous sfc downward uv+vis direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_vis_dir_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvisbmi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous sfc downward uv+vis diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_vis_dif_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvisdfi_cpl(ix)
        enddo
      enddo
    endif

    ! MEAN NET sfc nir direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_ir_dir_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nnirbm_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN NET sfc nir diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_ir_dif_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nnirdf_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN NET sfc uv+vis direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_vis_dir_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nvisbm_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN NET sfc uv+vis diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_vis_dif_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nvisdf_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! Instataneous net sfc nir direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_ir_dir_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nnirbmi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous net sfc nir diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_ir_dif_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nnirdfi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous net sfc uv+vis direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_vis_dir_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nvisbmi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous net sfc uv+vis diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_vis_dif_flx')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nvisdfi_cpl(ix)
        enddo
      enddo
    endif

    ! Land/Sea mask (sea:0,land:1)
    idx = queryfieldlist(exportFieldsList,'inst_land_sea_mask')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%slmsk_cpl(ix)
        enddo
      enddo
    endif

! Data from DYCORE:

    ! bottom layer temperature (t)
    idx = queryfieldlist(exportFieldsList,'inst_temp_height_lowest')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%t_bot)) then 
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%t_bot(ix)
          else 
            exportData(i,j,idx) = 0.0
          endif 
        enddo
      enddo
    endif

    ! bottom layer specific humidity (q)
    !!! CHECK if tracer 1 is for specific humidity !!!
    idx = queryfieldlist(exportFieldsList,'inst_spec_humid_height_lowest')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%tr_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%tr_bot(ix,1)
          else 
            exportData(i,j,idx)=0.0
          endif 
        enddo
      enddo
    endif

    ! bottom layer zonal wind (u)
    idx = queryfieldlist(exportFieldsList,'inst_zonal_wind_height_lowest')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%u_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%u_bot(ix)
          else
            exportData(i,j,idx) = 0.0
          endif 
        enddo
      enddo
    endif

    ! bottom layer meridionalw wind (v)
    idx = queryfieldlist(exportFieldsList,'inst_merid_wind_height_lowest')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%v_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%v_bot(ix)
          else 
            exportData(i,j,idx) = 0.0 
          endif 
        enddo
      enddo
    endif

    ! bottom layer pressure (p)
    idx = queryfieldlist(exportFieldsList,'inst_pres_height_lowest')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%p_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%p_bot(ix)
          else 
            exportData(i,j,idx) = 0.0
          endif 
        enddo
      enddo
    endif

    ! bottom layer height (z)
    idx = queryfieldlist(exportFieldsList,'inst_height_lowest')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%z_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%z_bot(ix)
          else 
            exportData(i,j,idx) = 0.0 
          endif 
        enddo
      enddo
    endif

! END Data from DYCORE.

    ! MEAN snow precipitation rate (kg/m2) ?????? checking unit ??????
    idx = queryfieldlist(exportFieldsList,'mean_fprec_rate')
    if (idx > 0 ) then
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%snow_cpl(ix) * rtime
        enddo
      enddo
    endif
    endif !cplflx 

!---
    ! Fill the export Fields for ESMF/NUOPC style coupling
    call fillExportFields(exportData)

!---
    if (IPD_Control%cplflx) then 
    ! zero out accumulated fields
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          IPD_Data(nb)%coupling%dusfc_cpl(ix)  = 0.
          IPD_Data(nb)%coupling%dvsfc_cpl(ix)  = 0.
          IPD_Data(nb)%coupling%dtsfc_cpl(ix)  = 0.
          IPD_Data(nb)%coupling%dqsfc_cpl(ix)  = 0.
          IPD_Data(nb)%coupling%dlwsfc_cpl(ix) = 0.
          IPD_Data(nb)%coupling%dswsfc_cpl(ix) = 0.
          IPD_Data(nb)%coupling%rain_cpl(ix)   = 0.
          IPD_Data(nb)%coupling%nlwsfc_cpl(ix) = 0.
          IPD_Data(nb)%coupling%nswsfc_cpl(ix) = 0.
          IPD_Data(nb)%coupling%dnirbm_cpl(ix) = 0.
          IPD_Data(nb)%coupling%dnirdf_cpl(ix) = 0.
          IPD_Data(nb)%coupling%dvisbm_cpl(ix) = 0.
          IPD_Data(nb)%coupling%dvisdf_cpl(ix) = 0.
          IPD_Data(nb)%coupling%nnirbm_cpl(ix) = 0.
          IPD_Data(nb)%coupling%nnirdf_cpl(ix) = 0.
          IPD_Data(nb)%coupling%nvisbm_cpl(ix) = 0.
          IPD_Data(nb)%coupling%nvisdf_cpl(ix) = 0.
          IPD_Data(nb)%coupling%snow_cpl(ix)   = 0.
        enddo
      enddo
    endif !cplflx
    if (mpp_pe() == mpp_root_pe()) print *,'end of setup_exportdata'

  end subroutine setup_exportdata

  subroutine addLsmask2grid(fcstgrid, rc)

    use ESMF
!
    implicit none
    type(ESMF_Grid)      :: fcstgrid
    integer, optional, intent(out) :: rc
!
!  local vars
    integer isc, iec, jsc, jec
    integer i, j, nb, ix
!    integer CLbnd(2), CUbnd(2), CCount(2), TLbnd(2), TUbnd(2), TCount(2)
    type(ESMF_StaggerLoc) :: staggerloc
    integer, allocatable :: lsmask(:,:)
    integer(kind=ESMF_KIND_I4), pointer  :: maskPtr(:,:)
!
    isc = IPD_control%isc
    iec = IPD_control%isc+IPD_control%nx-1
    jsc = IPD_control%jsc
    jec = IPD_control%jsc+IPD_control%ny-1
    allocate(lsmask(isc:iec,jsc:jec))
!
    do j=jsc,jec
      do i=isc,iec
        nb = Atm_block%blkno(i,j)
        ix = Atm_block%ixp(i,j)
! use land sea mask: land:1, ocean:0
        lsmask(i,j) = IPD_Data(nb)%SfcProp%slmsk(ix)
      enddo
    enddo
!
! Get mask
    call ESMF_GridAddItem(fcstgrid, itemflag=ESMF_GRIDITEM_MASK,   &
         staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

!    call ESMF_GridGetItemBounds(fcstgrid, itemflag=ESMF_GRIDITEM_MASK,   &
!         staggerloc=ESMF_STAGGERLOC_CENTER, computationalLBound=ClBnd,  &
!         computationalUBound=CUbnd, computationalCount=Ccount,  &
!         totalLBound=TLbnd, totalUBound=TUbnd, totalCount=Tcount, rc=rc)
!    print *,'in set up grid, aft add esmfgridadd item mask, rc=',rc, &
!     'ClBnd=',ClBnd,'CUbnd=',CUbnd,'Ccount=',Ccount, &
!     'TlBnd=',TlBnd,'TUbnd=',TUbnd,'Tcount=',Tcount
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridGetItem(fcstgrid, itemflag=ESMF_GRIDITEM_MASK,   &
         staggerloc=ESMF_STAGGERLOC_CENTER,farrayPtr=maskPtr, rc=rc)
!    print *,'in set up grid, aft get maskptr, rc=',rc, 'size=',size(maskPtr,1),size(maskPtr,2), &
!      'bound(maskPtr)=', LBOUND(maskPtr,1),LBOUND(maskPtr,2),UBOUND(maskPtr,1),UBOUND(maskPtr,2)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
!    
    do j=jsc,jec
      do i=isc,iec
        maskPtr(i-isc+1,j-jsc+1) = lsmask(i,j)
      enddo
    enddo
!      print *,'in set set lsmask, maskPtr=', maxval(maskPtr), minval(maskPtr)
!
    deallocate(lsmask)  

  end subroutine addLsmask2grid
!------------------------------------------------------------------------------

! DH*
#ifdef DHDEBUG
  subroutine MY_DIAGTOSCREEN()

     implicit none

     integer, parameter :: ISTART = 1
     integer, parameter :: IEND = 9999999

     integer, parameter :: KSTART = 1
     integer, parameter :: KEND = 9999999

     integer :: nblks, nb

     do nb=1,Atm_block%nblks
         call MY_DIAGTOSCREEN_RUN(IPD_Control,           &
                                  IPD_Data(nb)%Statein,  &
                                  IPD_Data(nb)%Stateout, &
                                  IPD_Data(nb)%Sfcprop,  &
                                  IPD_Data(nb)%Coupling, &
                                  IPD_Data(nb)%Grid,     &
                                  IPD_Data(nb)%Tbd,      &
                                  IPD_Data(nb)%Cldprop,  &
                                  IPD_Data(nb)%Radtend,  &
                                  IPD_Data(nb)%Intdiag,  &
                                  nb)
     end do

     contains

      subroutine MY_DIAGTOSCREEN_RUN (Model, Statein, Stateout, Sfcprop, Coupling,     &
                                      Grid, Tbd, Cldprop, Radtend, Diag, blkno)

         use mpi
         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type, GFS_statein_type,  &
                                          GFS_stateout_type, GFS_sfcprop_type, &
                                          GFS_coupling_type, GFS_grid_type,    &
                                          GFS_tbd_type, GFS_cldprop_type,      &
                                          GFS_radtend_type, GFS_diag_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in   ) :: Model
         type(GFS_statein_type),   intent(in   ) :: Statein
         type(GFS_stateout_type),  intent(in   ) :: Stateout
         type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
         type(GFS_coupling_type),  intent(in   ) :: Coupling
         type(GFS_grid_type),      intent(in   ) :: Grid
         type(GFS_tbd_type),       intent(in   ) :: Tbd
         type(GFS_cldprop_type),   intent(in   ) :: Cldprop
         type(GFS_radtend_type),   intent(in   ) :: Radtend
         type(GFS_diag_type),      intent(in   ) :: Diag
         integer,                  intent(in   ) :: blkno

         !--- local variables
         integer :: impi, iomp, ierr, n
         integer :: mpirank, mpisize, mpicomm
         integer, parameter :: omprank = 0

         mpicomm = MPI_COMM_WORLD
         mpirank = Model%me
         call MPI_COMM_SIZE(mpicomm, mpisize, ierr)

         do impi=0,mpisize-1
            if (mpirank==impi .and. mpirank<10) then
                ! Sfcprop
                call print_var(mpirank,omprank, blkno, 'Sfcprop%slmsk'    , Sfcprop%slmsk)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%lakemsk'  , Sfcprop%lakemsk)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%tsfc'     , Sfcprop%tsfc)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%tisfc'    , Sfcprop%tisfc)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%snowd'    , Sfcprop%snowd)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%zorl'     , Sfcprop%zorl)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%fice'     , Sfcprop%fice)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%hprim'    , Sfcprop%hprim)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%hprime'   , Sfcprop%hprime)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%sncovr'   , Sfcprop%sncovr)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%snoalb'   , Sfcprop%snoalb)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%alvsf'    , Sfcprop%alvsf)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%alnsf'    , Sfcprop%alnsf)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%alvwf'    , Sfcprop%alvwf)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%alnwf'    , Sfcprop%alnwf)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%facsf'    , Sfcprop%facsf)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%facwf'    , Sfcprop%facwf)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%slope'    , Sfcprop%slope)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%shdmin'   , Sfcprop%shdmin)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%shdmax'   , Sfcprop%shdmax)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%tg3'      , Sfcprop%tg3)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%vfrac'    , Sfcprop%vfrac)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%vtype'    , Sfcprop%vtype)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%stype'    , Sfcprop%stype)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%uustar'   , Sfcprop%uustar)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%oro'      , Sfcprop%oro)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%oro_uf'   , Sfcprop%oro_uf)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%hice'     , Sfcprop%hice)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%weasd'    , Sfcprop%weasd)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%canopy'   , Sfcprop%canopy)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%ffmm'     , Sfcprop%ffmm)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%ffhh'     , Sfcprop%ffhh)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%f10m'     , Sfcprop%f10m)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%tprcp'    , Sfcprop%tprcp)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%srflag'   , Sfcprop%srflag)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%sr'       , Sfcprop%sr)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%slc'      , Sfcprop%slc)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%smc'      , Sfcprop%smc)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%stc'      , Sfcprop%stc)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%t2m'      , Sfcprop%t2m)
                call print_var(mpirank,omprank, blkno, 'Sfcprop%q2m'      , Sfcprop%q2m)
                if (Model%nstf_name(1)>0) then
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%tref    ', Sfcprop%tref)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%z_c     ', Sfcprop%z_c)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%c_0     ', Sfcprop%c_0)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%c_d     ', Sfcprop%c_d)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%w_0     ', Sfcprop%w_0)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%w_d     ', Sfcprop%w_d)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%xt      ', Sfcprop%xt)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%xs      ', Sfcprop%xs)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%xu      ', Sfcprop%xu)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%xv      ', Sfcprop%xv)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%xz      ', Sfcprop%xz)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%zm      ', Sfcprop%zm)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%xtts    ', Sfcprop%xtts)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%xzts    ', Sfcprop%xzts)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%d_conv  ', Sfcprop%d_conv)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%ifd     ', Sfcprop%ifd)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%dt_cool ', Sfcprop%dt_cool)
                   call print_var(mpirank,omprank, blkno, 'Sfcprop%qrain   ', Sfcprop%qrain)
                end if
                ! CCPP only
                !if (Model%lsm == Model%lsm_ruc) then
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%sh2o',        Sfcprop%sh2o)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%smois',       Sfcprop%smois)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%tslb',        Sfcprop%tslb)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%zs',          Sfcprop%zs)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%clw_surf',    Sfcprop%clw_surf)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%qwv_surf',    Sfcprop%qwv_surf)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%cndm_surf',   Sfcprop%cndm_surf)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%flag_frsoil', Sfcprop%flag_frsoil)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%rhofr',       Sfcprop%rhofr)
                !   call print_var(mpirank,omprank, blkno, 'Sfcprop%tsnow',       Sfcprop%tsnow)
                !end if
                ! Radtend
                call print_var(mpirank,omprank, blkno, 'Radtend%sfcfsw%upfxc', Radtend%sfcfsw(:)%upfxc)
                call print_var(mpirank,omprank, blkno, 'Radtend%sfcfsw%dnfxc', Radtend%sfcfsw(:)%dnfxc)
                call print_var(mpirank,omprank, blkno, 'Radtend%sfcfsw%upfx0', Radtend%sfcfsw(:)%upfx0)
                call print_var(mpirank,omprank, blkno, 'Radtend%sfcfsw%dnfx0', Radtend%sfcfsw(:)%dnfx0)
                call print_var(mpirank,omprank, blkno, 'Radtend%sfcflw%upfxc', Radtend%sfcflw(:)%upfxc)
                call print_var(mpirank,omprank, blkno, 'Radtend%sfcflw%upfx0', Radtend%sfcflw(:)%upfx0)
                call print_var(mpirank,omprank, blkno, 'Radtend%sfcflw%dnfxc', Radtend%sfcflw(:)%dnfxc)
                call print_var(mpirank,omprank, blkno, 'Radtend%sfcflw%dnfx0', Radtend%sfcflw(:)%dnfx0)
                call print_var(mpirank,omprank, blkno, 'Radtend%htrsw',        Radtend%htrsw)
                call print_var(mpirank,omprank, blkno, 'Radtend%htrlw',        Radtend%htrlw)
                call print_var(mpirank,omprank, blkno, 'Radtend%sfalb',        Radtend%sfalb)
                call print_var(mpirank,omprank, blkno, 'Radtend%coszen',       Radtend%coszen)
                call print_var(mpirank,omprank, blkno, 'Radtend%tsflw',        Radtend%tsflw)
                call print_var(mpirank,omprank, blkno, 'Radtend%semis',        Radtend%semis)
                call print_var(mpirank,omprank, blkno, 'Radtend%coszdg',       Radtend%coszdg)
                call print_var(mpirank,omprank, blkno, 'Radtend%swhc',         Radtend%swhc)
                call print_var(mpirank,omprank, blkno, 'Radtend%lwhc',         Radtend%lwhc)
                call print_var(mpirank,omprank, blkno, 'Radtend%lwhd',         Radtend%lwhd)
                ! Tbd
                call print_var(mpirank,omprank, blkno, 'Tbd%icsdsw'          , Tbd%icsdsw)
                call print_var(mpirank,omprank, blkno, 'Tbd%icsdlw'          , Tbd%icsdlw)
                call print_var(mpirank,omprank, blkno, 'Tbd%ozpl'            , Tbd%ozpl)
                call print_var(mpirank,omprank, blkno, 'Tbd%h2opl'           , Tbd%h2opl)
                call print_var(mpirank,omprank, blkno, 'Tbd%rann'            , Tbd%rann)
                call print_var(mpirank,omprank, blkno, 'Tbd%acv'             , Tbd%acv)
                call print_var(mpirank,omprank, blkno, 'Tbd%acvb'            , Tbd%acvb)
                call print_var(mpirank,omprank, blkno, 'Tbd%acvt'            , Tbd%acvt)
                if (Model%do_sppt) then
                  call print_var(mpirank,omprank, blkno, 'Tbd%dtdtr'         , Tbd%dtdtr)
                  call print_var(mpirank,omprank, blkno, 'Tbd%dtotprcp'      , Tbd%dtotprcp)
                  call print_var(mpirank,omprank, blkno, 'Tbd%dcnvprcp'      , Tbd%dcnvprcp)
                  call print_var(mpirank,omprank, blkno, 'Tbd%drain_cpl'     , Tbd%drain_cpl)
                  call print_var(mpirank,omprank, blkno, 'Tbd%dsnow_cpl'     , Tbd%dsnow_cpl)
                end if
                call print_var(mpirank,omprank, blkno, 'Tbd%phy_fctd'        , Tbd%phy_fctd)
                call print_var(mpirank,omprank, blkno, 'Tbd%phy_f2d'         , Tbd%phy_f2d)
                call print_var(mpirank,omprank, blkno, 'Tbd%phy_f3d'         , Tbd%phy_f3d)
                do n=1,size(Tbd%phy_f3d(1,1,:))
                    call print_var(mpirank,omprank, blkno, 'Tbd%phy_f3d_n'   , Tbd%phy_f3d(:,:,n))
                end do
                !call print_var(mpirank,omprank, blkno, 'Tbd%in_nm'           , Tbd%in_nm)
                !call print_var(mpirank,omprank, blkno, 'Tbd%ccn_nm'          , Tbd%ccn_nm)
                call print_var(mpirank,omprank, blkno, 'Tbd%aer_nm'          , Tbd%aer_nm)
                ! Diag
                !call print_var(mpirank,omprank, blkno, 'Diag%fluxr       ',    Diag%fluxr)
                !do n=1,size(Diag%fluxr(1,:))
                !    call print_var(mpirank,omprank, blkno, 'Diag%fluxr_n ',    Diag%fluxr(:,n))
                !end do
                !call print_var(mpirank,omprank, blkno, 'Diag%srunoff     ',    Diag%srunoff)
                call print_var(mpirank,omprank, blkno, 'Diag%evbsa       ',    Diag%evbsa)
                call print_var(mpirank,omprank, blkno, 'Diag%evcwa       ',    Diag%evcwa)
                call print_var(mpirank,omprank, blkno, 'Diag%snohfa      ',    Diag%snohfa)
                call print_var(mpirank,omprank, blkno, 'Diag%transa      ',    Diag%transa)
                call print_var(mpirank,omprank, blkno, 'Diag%sbsnoa      ',    Diag%sbsnoa)
                call print_var(mpirank,omprank, blkno, 'Diag%snowca      ',    Diag%snowca)
                call print_var(mpirank,omprank, blkno, 'Diag%soilm       ',    Diag%soilm)
                call print_var(mpirank,omprank, blkno, 'Diag%tmpmin      ',    Diag%tmpmin)
                call print_var(mpirank,omprank, blkno, 'Diag%tmpmax      ',    Diag%tmpmax)
                call print_var(mpirank,omprank, blkno, 'Diag%dusfc       ',    Diag%dusfc)
                call print_var(mpirank,omprank, blkno, 'Diag%dvsfc       ',    Diag%dvsfc)
                call print_var(mpirank,omprank, blkno, 'Diag%dtsfc       ',    Diag%dtsfc)
                call print_var(mpirank,omprank, blkno, 'Diag%dqsfc       ',    Diag%dqsfc)
                call print_var(mpirank,omprank, blkno, 'Diag%totprcp     ',    Diag%totprcp)
                call print_var(mpirank,omprank, blkno, 'Diag%totice      ',    Diag%totice)
                call print_var(mpirank,omprank, blkno, 'Diag%totsnw      ',    Diag%totsnw)
                call print_var(mpirank,omprank, blkno, 'Diag%totgrp      ',    Diag%totgrp)
                call print_var(mpirank,omprank, blkno, 'Diag%totprcpb    ',    Diag%totprcpb)
                call print_var(mpirank,omprank, blkno, 'Diag%toticeb     ',    Diag%toticeb)
                call print_var(mpirank,omprank, blkno, 'Diag%totsnwb     ',    Diag%totsnwb)
                call print_var(mpirank,omprank, blkno, 'Diag%totgrpb     ',    Diag%totgrpb)
                call print_var(mpirank,omprank, blkno, 'Diag%suntim      ',    Diag%suntim)
                !call print_var(mpirank,omprank, blkno, 'Diag%runoff      ',    Diag%runoff)
                call print_var(mpirank,omprank, blkno, 'Diag%ep          ',    Diag%ep)
                call print_var(mpirank,omprank, blkno, 'Diag%cldwrk      ',    Diag%cldwrk)
                call print_var(mpirank,omprank, blkno, 'Diag%dugwd       ',    Diag%dugwd)
                call print_var(mpirank,omprank, blkno, 'Diag%dvgwd       ',    Diag%dvgwd)
                call print_var(mpirank,omprank, blkno, 'Diag%psmean      ',    Diag%psmean)
                call print_var(mpirank,omprank, blkno, 'Diag%cnvprcp     ',    Diag%cnvprcp)
                call print_var(mpirank,omprank, blkno, 'Diag%cnvprcpb    ',    Diag%cnvprcpb)
                call print_var(mpirank,omprank, blkno, 'Diag%spfhmin     ',    Diag%spfhmin)
                call print_var(mpirank,omprank, blkno, 'Diag%spfhmax     ',    Diag%spfhmax)
                call print_var(mpirank,omprank, blkno, 'Diag%u10mmax     ',    Diag%u10mmax)
                call print_var(mpirank,omprank, blkno, 'Diag%v10mmax     ',    Diag%v10mmax)
                call print_var(mpirank,omprank, blkno, 'Diag%wind10mmax  ',    Diag%wind10mmax)
                call print_var(mpirank,omprank, blkno, 'Diag%rain        ',    Diag%rain)
                call print_var(mpirank,omprank, blkno, 'Diag%rainc       ',    Diag%rainc)
                call print_var(mpirank,omprank, blkno, 'Diag%ice         ',    Diag%ice)
                call print_var(mpirank,omprank, blkno, 'Diag%snow        ',    Diag%snow)
                call print_var(mpirank,omprank, blkno, 'Diag%graupel     ',    Diag%graupel)
                call print_var(mpirank,omprank, blkno, 'Diag%u10m        ',    Diag%u10m)
                call print_var(mpirank,omprank, blkno, 'Diag%v10m        ',    Diag%v10m)
                call print_var(mpirank,omprank, blkno, 'Diag%dpt2m       ',    Diag%dpt2m)
                call print_var(mpirank,omprank, blkno, 'Diag%zlvl        ',    Diag%zlvl)
                call print_var(mpirank,omprank, blkno, 'Diag%psurf       ',    Diag%psurf)
                call print_var(mpirank,omprank, blkno, 'Diag%hpbl        ',    Diag%hpbl)
                call print_var(mpirank,omprank, blkno, 'Diag%pwat        ',    Diag%pwat)
                call print_var(mpirank,omprank, blkno, 'Diag%t1          ',    Diag%t1)
                call print_var(mpirank,omprank, blkno, 'Diag%q1          ',    Diag%q1)
                call print_var(mpirank,omprank, blkno, 'Diag%u1          ',    Diag%u1)
                call print_var(mpirank,omprank, blkno, 'Diag%v1          ',    Diag%v1)
                call print_var(mpirank,omprank, blkno, 'Diag%chh         ',    Diag%chh)
                call print_var(mpirank,omprank, blkno, 'Diag%cmm         ',    Diag%cmm)
                call print_var(mpirank,omprank, blkno, 'Diag%epi         ',    Diag%epi)
                call print_var(mpirank,omprank, blkno, 'Diag%smcwlt2     ',    Diag%smcwlt2)
                call print_var(mpirank,omprank, blkno, 'Diag%smcref2     ',    Diag%smcref2)
                call print_var(mpirank,omprank, blkno, 'Diag%tdomr       ',    Diag%tdomr)
                call print_var(mpirank,omprank, blkno, 'Diag%tdomzr      ',    Diag%tdomzr)
                call print_var(mpirank,omprank, blkno, 'Diag%tdomip      ',    Diag%tdomip)
                call print_var(mpirank,omprank, blkno, 'Diag%tdoms       ',    Diag%tdoms)
                ! CCPP only
                !call print_var(mpirank,omprank, blkno, 'Diag%snowfallac  ',    Diag%snowfallac)
                !call print_var(mpirank,omprank, blkno, 'Diag%acsnow      ',    Diag%acsnow)
                call print_var(mpirank,omprank, blkno, 'Diag%skebu_wts   ',    Diag%skebu_wts)
                call print_var(mpirank,omprank, blkno, 'Diag%skebv_wts   ',    Diag%skebv_wts)
                call print_var(mpirank,omprank, blkno, 'Diag%sppt_wts    ',    Diag%sppt_wts)
                call print_var(mpirank,omprank, blkno, 'Diag%shum_wts    ',    Diag%shum_wts)
                call print_var(mpirank,omprank, blkno, 'Diag%zmtnblck    ',    Diag%zmtnblck)
                if (Model%ldiag3d) then
                  call print_var(mpirank,omprank, blkno, 'Diag%du3dt       ',    Diag%du3dt)
                  do n=1,size(Diag%du3dt(1,1,:))
                    call print_var(mpirank,omprank, blkno, 'Diag%du3dt_n     ',  Diag%du3dt(:,:,n))
                  end do
                  call print_var(mpirank,omprank, blkno, 'Diag%dv3dt       ',    Diag%dv3dt)
                  do n=1,size(Diag%dv3dt(1,1,:))
                    call print_var(mpirank,omprank, blkno, 'Diag%dv3dt_n     ',  Diag%dv3dt(:,:,n))
                  end do
                  call print_var(mpirank,omprank, blkno, 'Diag%dt3dt       ',    Diag%dt3dt)
                  do n=1,size(Diag%dt3dt(1,1,:))
                    call print_var(mpirank,omprank, blkno, 'Diag%dt3dt_n     ',  Diag%dt3dt(:,:,n))
                  end do
                  call print_var(mpirank,omprank, blkno, 'Diag%dq3dt       ',    Diag%dq3dt)
                  do n=1,size(Diag%dq3dt(1,1,:))
                    call print_var(mpirank,omprank, blkno, 'Diag%dq3dt_n     ',  Diag%dq3dt(:,:,n))
                  end do
                  call print_var(mpirank,omprank, blkno, 'Diag%upd_mf      ',    Diag%upd_mf)
                  call print_var(mpirank,omprank, blkno, 'Diag%dwn_mf      ',    Diag%dwn_mf)
                  call print_var(mpirank,omprank, blkno, 'Diag%det_mf      ',    Diag%det_mf)
                  call print_var(mpirank,omprank, blkno, 'Diag%cldcov      ',    Diag%cldcov)
                end if
                if(Model%lradar) then
                  call print_var(mpirank,omprank, blkno, 'Diag%refl_10cm   ',  Diag%refl_10cm)
                end if
                ! CCPP only
                !if (Model%do_mynnedmf) then
                !  call print_var(mpirank,omprank, blkno, 'Diag%edmf_a      ',  Diag%edmf_a)
                !  call print_var(mpirank,omprank, blkno, 'Diag%edmf_w      ',  Diag%edmf_w)
                !  call print_var(mpirank,omprank, blkno, 'Diag%edmf_qt     ',  Diag%edmf_qt)
                !  call print_var(mpirank,omprank, blkno, 'Diag%edmf_thl    ',  Diag%edmf_thl)
                !  call print_var(mpirank,omprank, blkno, 'Diag%edmf_ent    ',  Diag%edmf_ent)
                !  call print_var(mpirank,omprank, blkno, 'Diag%edmf_qc     ',  Diag%edmf_qc)
                !  call print_var(mpirank,omprank, blkno, 'Diag%nupdraft    ',  Diag%nupdraft)
                !  call print_var(mpirank,omprank, blkno, 'Diag%maxMF       ',  Diag%maxMF)
                !  call print_var(mpirank,omprank, blkno, 'Diag%ktop_shallow',  Diag%ktop_shallow)
                !  call print_var(mpirank,omprank, blkno, 'Diag%exch_h      ',  Diag%exch_h)
                !  call print_var(mpirank,omprank, blkno, 'Diag%exch_m      ',  Diag%exch_m)
                !end if
                ! Statein
                call print_var(mpirank,omprank, blkno, 'Statein%phii'    ,     Statein%phii)
                call print_var(mpirank,omprank, blkno, 'Statein%prsi'    ,     Statein%prsi)
                call print_var(mpirank,omprank, blkno, 'Statein%prsik'   ,     Statein%prsik)
                call print_var(mpirank,omprank, blkno, 'Statein%phil'    ,     Statein%phil)
                call print_var(mpirank,omprank, blkno, 'Statein%prsl'    ,     Statein%prsl)
                call print_var(mpirank,omprank, blkno, 'Statein%prslk'   ,     Statein%prslk)
                call print_var(mpirank,omprank, blkno, 'Statein%pgr'     ,     Statein%pgr)
                call print_var(mpirank,omprank, blkno, 'Statein%ugrs'    ,     Statein%ugrs)
                call print_var(mpirank,omprank, blkno, 'Statein%vgrs'    ,     Statein%vgrs)
                call print_var(mpirank,omprank, blkno, 'Statein%vvl'     ,     Statein%vvl)
                call print_var(mpirank,omprank, blkno, 'Statein%tgrs'    ,     Statein%tgrs)
                call print_var(mpirank,omprank, blkno, 'Statein%qgrs'    ,     Statein%qgrs)
                do n=1,size(Statein%qgrs(1,1,:))
                   call print_var(mpirank,omprank, blkno, 'Statein%qgrs_n',    Statein%qgrs(:,:,n))
                end do
                call print_var(mpirank,omprank, blkno, 'Statein%diss_est',     Statein%diss_est)
                call print_var(mpirank,omprank, blkno, 'Statein%smc'     ,     Statein%smc)
                call print_var(mpirank,omprank, blkno, 'Statein%stc'     ,     Statein%stc)
                call print_var(mpirank,omprank, blkno, 'Statein%slc'     ,     Statein%slc)
                ! Stateout
                call print_var(mpirank,omprank, blkno, 'Stateout%gu0',         Stateout%gu0)
                call print_var(mpirank,omprank, blkno, 'Stateout%gv0',         Stateout%gv0)
                call print_var(mpirank,omprank, blkno, 'Stateout%gt0',         Stateout%gt0)
                call print_var(mpirank,omprank, blkno, 'Stateout%gq0',         Stateout%gq0)
                do n=1,size(Stateout%gq0(1,1,:))
                   call print_var(mpirank,omprank, blkno, 'Stateout%gq0_n',    Stateout%gq0(:,:,n))
                end do
                ! Coupling
                call print_var(mpirank,omprank, blkno, 'Coupling%nirbmdi', Coupling%nirbmdi)
                call print_var(mpirank,omprank, blkno, 'Coupling%nirdfdi', Coupling%nirdfdi)
                call print_var(mpirank,omprank, blkno, 'Coupling%visbmdi', Coupling%visbmdi)
                call print_var(mpirank,omprank, blkno, 'Coupling%visdfdi', Coupling%visdfdi)
                call print_var(mpirank,omprank, blkno, 'Coupling%nirbmui', Coupling%nirbmui)
                call print_var(mpirank,omprank, blkno, 'Coupling%nirdfui', Coupling%nirdfui)
                call print_var(mpirank,omprank, blkno, 'Coupling%visbmui', Coupling%visbmui)
                call print_var(mpirank,omprank, blkno, 'Coupling%visdfui', Coupling%visdfui)
                call print_var(mpirank,omprank, blkno, 'Coupling%sfcdsw ', Coupling%sfcdsw )
                call print_var(mpirank,omprank, blkno, 'Coupling%sfcnsw ', Coupling%sfcnsw )
                call print_var(mpirank,omprank, blkno, 'Coupling%sfcdlw ', Coupling%sfcdlw )
                if (Model%cplflx .or. Model%do_sppt) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%rain_cpl', Coupling%rain_cpl)
                   call print_var(mpirank,omprank, blkno, 'Coupling%snow_cpl', Coupling%snow_cpl)
                end if
                if (Model%cplflx) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%slimskin_cpl', Coupling%slimskin_cpl )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dusfcin_cpl ', Coupling%dusfcin_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dvsfcin_cpl ', Coupling%dvsfcin_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dtsfcin_cpl ', Coupling%dtsfcin_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dqsfcin_cpl ', Coupling%dqsfcin_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%ulwsfcin_cpl', Coupling%ulwsfcin_cpl )
                   call print_var(mpirank,omprank, blkno, 'Coupling%tseain_cpl  ', Coupling%tseain_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%tisfcin_cpl ', Coupling%tisfcin_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%ficein_cpl  ', Coupling%ficein_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%hicein_cpl  ', Coupling%hicein_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%hsnoin_cpl  ', Coupling%hsnoin_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dusfc_cpl   ', Coupling%dusfc_cpl    )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dvsfc_cpl   ', Coupling%dvsfc_cpl    )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dtsfc_cpl   ', Coupling%dtsfc_cpl    )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dqsfc_cpl   ', Coupling%dqsfc_cpl    )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dlwsfc_cpl  ', Coupling%dlwsfc_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dswsfc_cpl  ', Coupling%dswsfc_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dnirbm_cpl  ', Coupling%dnirbm_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dnirdf_cpl  ', Coupling%dnirdf_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dvisbm_cpl  ', Coupling%dvisbm_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dvisdf_cpl  ', Coupling%dvisdf_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nlwsfc_cpl  ', Coupling%nlwsfc_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nswsfc_cpl  ', Coupling%nswsfc_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nnirbm_cpl  ', Coupling%nnirbm_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nnirdf_cpl  ', Coupling%nnirdf_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nvisbm_cpl  ', Coupling%nvisbm_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nvisdf_cpl  ', Coupling%nvisdf_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dusfci_cpl  ', Coupling%dusfci_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dvsfci_cpl  ', Coupling%dvsfci_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dtsfci_cpl  ', Coupling%dtsfci_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dqsfci_cpl  ', Coupling%dqsfci_cpl   )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dlwsfci_cpl ', Coupling%dlwsfci_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dswsfci_cpl ', Coupling%dswsfci_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dnirbmi_cpl ', Coupling%dnirbmi_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dnirdfi_cpl ', Coupling%dnirdfi_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dvisbmi_cpl ', Coupling%dvisbmi_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dvisdfi_cpl ', Coupling%dvisdfi_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nlwsfci_cpl ', Coupling%nlwsfci_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nswsfci_cpl ', Coupling%nswsfci_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nnirbmi_cpl ', Coupling%nnirbmi_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nnirdfi_cpl ', Coupling%nnirdfi_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nvisbmi_cpl ', Coupling%nvisbmi_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%nvisdfi_cpl ', Coupling%nvisdfi_cpl  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%t2mi_cpl    ', Coupling%t2mi_cpl     )
                   call print_var(mpirank,omprank, blkno, 'Coupling%q2mi_cpl    ', Coupling%q2mi_cpl     )
                   call print_var(mpirank,omprank, blkno, 'Coupling%u10mi_cpl   ', Coupling%u10mi_cpl    )
                   call print_var(mpirank,omprank, blkno, 'Coupling%v10mi_cpl   ', Coupling%v10mi_cpl    )
                   call print_var(mpirank,omprank, blkno, 'Coupling%tsfci_cpl   ', Coupling%tsfci_cpl    )
                   call print_var(mpirank,omprank, blkno, 'Coupling%psurfi_cpl  ', Coupling%psurfi_cpl   )
                end if
                if (Model%cplchm) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%rain_cpl ', Coupling%rain_cpl )
                   call print_var(mpirank,omprank, blkno, 'Coupling%rainc_cpl', Coupling%rainc_cpl)
                   call print_var(mpirank,omprank, blkno, 'Coupling%ushfsfci ', Coupling%ushfsfci )
                   call print_var(mpirank,omprank, blkno, 'Coupling%dkt      ', Coupling%dkt      )
                end if
                if (Model%do_sppt) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%sppt_wts', Coupling%sppt_wts)
                end if
                if (Model%do_shum) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%shum_wts', Coupling%shum_wts)
                end if
                if (Model%do_skeb) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%skebu_wts', Coupling%skebu_wts)
                   call print_var(mpirank,omprank, blkno, 'Coupling%skebv_wts', Coupling%skebv_wts)
                end if
                if (Model%do_sfcperts) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%sfc_wts', Coupling%sfc_wts)
                end if
                if (Model%lgocart .or. Model%ldiag3d) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%dqdti  ', Coupling%dqdti  )
                   call print_var(mpirank,omprank, blkno, 'Coupling%cnvqci ', Coupling%cnvqci )
                   call print_var(mpirank,omprank, blkno, 'Coupling%upd_mfi', Coupling%upd_mfi)
                   call print_var(mpirank,omprank, blkno, 'Coupling%dwn_mfi', Coupling%dwn_mfi)
                   call print_var(mpirank,omprank, blkno, 'Coupling%det_mfi', Coupling%det_mfi)
                   call print_var(mpirank,omprank, blkno, 'Coupling%cldcovi', Coupling%cldcovi)
                end if
                if(Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
                   call print_var(mpirank,omprank, blkno, 'Coupling%nwfa2d', Coupling%nwfa2d)
                   call print_var(mpirank,omprank, blkno, 'Coupling%nifa2d', Coupling%nifa2d)
                end if
                ! Grid
                call print_var(mpirank,omprank, blkno, 'Grid%xlon  ', Grid%xlon  )
                call print_var(mpirank,omprank, blkno, 'Grid%xlat  ', Grid%xlat  )
                call print_var(mpirank,omprank, blkno, 'Grid%xlat_d', Grid%xlat_d)
                call print_var(mpirank,omprank, blkno, 'Grid%sinlat', Grid%sinlat)
                call print_var(mpirank,omprank, blkno, 'Grid%coslat', Grid%coslat)
                call print_var(mpirank,omprank, blkno, 'Grid%area  ', Grid%area  )
                call print_var(mpirank,omprank, blkno, 'Grid%dx    ', Grid%dx    )
                if (Model%ntoz > 0) then
                   call print_var(mpirank,omprank, blkno, 'Grid%ddy_o3   ', Grid%ddy_o3   )
                   call print_var(mpirank,omprank, blkno, 'Grid%jindx1_o3', Grid%jindx1_o3)
                   call print_var(mpirank,omprank, blkno, 'Grid%jindx2_o3', Grid%jindx2_o3)
                endif
                if (Model%h2o_phys) then
                   call print_var(mpirank,omprank, blkno, 'Grid%ddy_h   ', Grid%ddy_h   )
                   call print_var(mpirank,omprank, blkno, 'Grid%jindx1_h', Grid%jindx1_h)
                   call print_var(mpirank,omprank, blkno, 'Grid%jindx2_h', Grid%jindx2_h)
                endif
                ! Model/Control
                ! not yet
            end if
         end do

      end subroutine MY_DIAGTOSCREEN_RUN

  end subroutine MY_DIAGTOSCREEN

      subroutine print_logic_0d(mpirank,omprank,blkno,name,var)

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          logical, intent(in) :: var

          write(0,'(2a,3i6,1x,l)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_logic_0d

      subroutine print_int_0d(mpirank,omprank,blkno,name,var)

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          integer, intent(in) :: var

          write(0,'(2a,3i6,i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_int_0d

      subroutine print_int_1d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          integer, intent(in) :: var(:)

          integer :: i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_int(size(var),var), minval(var), maxval(var)
#else
          do i=ISTART,min(IEND,size(var(:)))
              write(0,'(2a,3i6,i6,i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, var(i)
          end do
#endif

      end subroutine print_int_1d

      subroutine print_real_0d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var

          write(0,'(2a,3i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_real_0d

      subroutine print_real_1d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:)

          integer :: i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),var), minval(var), maxval(var)
#else
          do i=ISTART,min(IEND,size(var(:)))
              write(0,'(2a,3i6,i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, var(i)
          end do
#endif

      end subroutine print_real_1d

      subroutine print_real_2d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:)
       
          integer :: k, i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),reshape(var,(/size(var)/))), minval(var), maxval(var)
#else
          do i=ISTART,min(IEND,size(var(:,1)))
              do k=KSTART,min(KEND,size(var(1,:)))
                  write(0,'(2a,3i6,2i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, var(i,k)
              end do
          end do
#endif

      end subroutine print_real_2d

      subroutine print_real_3d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:,:)

          integer :: k, i, l

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),reshape(var,(/size(var)/))), minval(var), maxval(var)
#else
          do i=ISTART,min(IEND,size(var(:,1,1)))
              do k=KSTART,min(KEND,size(var(1,:,1)))
                  do l=1,size(var(1,1,:))
                      write(0,'(2a,3i6,3i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, l, var(i,k,l)
                  end do
              end do
          end do
#endif

      end subroutine print_real_3d

      function chksum_int(N, var) result(hash)
          implicit none
          integer, intent(in) :: N
          integer, dimension(1:N), intent(in) :: var
          integer*8, dimension(1:N) :: int_var
          integer*8 :: a, b, i, hash
          integer*8, parameter :: mod_adler=65521

          a=1
          b=0
          i=1
          hash = 0
          int_var = TRANSFER(var, a, N)

          do i= 1, N
              a = MOD(a + int_var(i), mod_adler)
              b = MOD(b+a, mod_adler)
          end do

          hash = ior(b * 65536, a)

      end function chksum_int

      function chksum_real(N, var) result(hash)
          use machine,               only: kind_phys
          implicit none
          integer, intent(in) :: N
          real(kind_phys), dimension(1:N), intent(in) :: var
          integer*8, dimension(1:N) :: int_var
          integer*8 :: a, b, i, hash
          integer*8, parameter :: mod_adler=65521

          a=1
          b=0
          i=1
          hash = 0
          int_var = TRANSFER(var, a, N)

          do i= 1, N
              a = MOD(a + int_var(i), mod_adler)
              b = MOD(b+a, mod_adler)
          end do

          hash = ior(b * 65536, a)

     end function chksum_real
#endif
! *DH

end module atmos_model_mod

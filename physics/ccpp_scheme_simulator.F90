! ########################################################################################
! 
! CCPP scheme to replace physics schemes with simulated data tendencies.
!
! Description:
!
! ########################################################################################
module ccpp_scheme_simulator
  use machine, only: kind_phys
  use netcdf
#ifdef MPI
  use mpi
#endif
  implicit none

  ! Type containing physics tendencies for a physics process.
  type phys_tend
     real(kind_phys), dimension(:,:),   pointer :: T
     real(kind_phys), dimension(:,:),   pointer :: u
     real(kind_phys), dimension(:,:),   pointer :: v
     real(kind_phys), dimension(:,:,:), pointer :: q
  end type phys_tend

  ! This type contains the meta information and data for each physics process.
  type base_physics_process
     character(len=16) :: name
     logical           :: time_split = .false.
     logical           :: use_sim    = .false.
     integer           :: order
     type(phys_tend)   :: tend
  end type base_physics_process

  ! This array contains the governing information on how to advance the physics timestep.
  type(base_physics_process),dimension(:), allocatable :: &
       physics_process

  integer :: nPhysProcess

  ! ########################################################################################
  !
  ! Configuration for CCPP scheme simulator. Set in namelist. Used during initialization to 
  ! populate "physics_processes" type array.
  !
  ! ########################################################################################

  ! For each process there is a corresponding namelist entry, which is constructed as follows:
  ! {use_scheme_sim[0(no)/1(yes)], time_split[0(no)/1(yes)], order[1:nPhysProcess]}
  integer, dimension(3) ::    &
       proc_LWRAD_config = (/0,0,0/), &
       proc_SWRAD_config = (/0,0,0/), &
       proc_PBL_config   = (/0,0,0/), &
       proc_GWD_config   = (/0,0,0/), &
       proc_SCNV_config  = (/0,0,0/), &
       proc_DCNV_config  = (/0,0,0/), &
       proc_cldMP_config = (/0,0,0/)

  ! Activation flag for scheme.
  logical :: do_ccpp_scheme_simulator = .false.

  ! Data driven physics tendencies
  integer :: nlev_data, ntime_data
  real(kind_phys), allocatable, dimension(:)   :: time_data
  real(kind_phys), allocatable, dimension(:,:), target :: dTdt_LWRAD_data,               &
       dTdt_SWRAD_data, dTdt_PBL_data, dudt_PBL_data, dvdt_PBL_data, dTdt_GWD_data,      &
       dudt_GWD_data, dvdt_GWD_data, dTdt_SCNV_data, dudt_SCNV_data, dvdt_SCNV_data,     &
       dTdt_DCNV_data, dudt_DCNV_data, dvdt_DCNV_data, dTdt_cldMP_data
  real(kind_phys), allocatable, dimension(:,:,:), target :: dqdt_PBL_data,               &
       dqdt_SCNV_data, dqdt_DCNV_data, dqdt_cldMP_data

  ! Scheme initialization flag.
  logical :: module_initialized = .false.

  ! Order in process loop for "active" physics process.
  integer :: iactive_scheme

  public ccpp_scheme_simulator_init, ccpp_scheme_simulator_run
contains

  ! ######################################################################################
  !
  ! SUBROUTINE ccpp_scheme_simulator_init
  !
  ! ######################################################################################
!! \section arg_table_ccpp_scheme_simulator_init
!! \htmlinclude ccpp_scheme_simulator_init.html
!!
  subroutine ccpp_scheme_simulator_init(mpirank, mpiroot, mpicomm, nlunit, nml_file,     &
       errmsg, errflg)

    ! Inputs
    integer,          intent (in) :: mpirank, mpiroot, mpicomm, nlunit
    character(len=*), intent (in) :: nml_file

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    integer :: ncid, dimID, varID, status, nlon, nlat, ios, iprc
    character(len=256) :: fileIN
    logical :: exists
    integer,parameter :: nTrc = 1 ! Only specific humodty for now, but preserve 3 dimensionality

    ! Switches for input data
    logical :: have_dTdt_LWRAD_data     = .false., &
               have_dTdt_SWRAD_data     = .false., &
               have_dTdt_PBL_data       = .false., &
               have_dqdt_PBL_data       = .false., &
               have_dudt_PBL_data       = .false., &
               have_dvdt_PBL_data       = .false., &
               have_dTdt_GWD_data       = .false., &
               have_dudt_GWD_data       = .false., &
               have_dvdt_GWD_data       = .false., &
               have_dTdt_SCNV_data      = .false., &
               have_dudt_SCNV_data      = .false., &
               have_dvdt_SCNV_data      = .false., &
               have_dqdt_SCNV_data      = .false., &
               have_dTdt_DCNV_data      = .false., &
               have_dudt_DCNV_data      = .false., &
               have_dvdt_DCNV_data      = .false., &
               have_dqdt_DCNV_data      = .false., &
               have_dTdt_cldMP_data     = .false., &
               have_dqdt_cldMP_data     = .false.

    ! Namelist
    namelist / scm_data_nml / fileIN, nPhysProcess, proc_LWRAD_config, proc_SWRAD_config,  &
         proc_PBL_config, proc_GWD_config, proc_SCNV_config, proc_DCNV_config,             &
         proc_cldMP_config

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (module_initialized) return
    module_initialized = .true.

    ! ######################################################################################
    !
    ! Read in namelist
    !
    ! ######################################################################################
    inquire (file = trim (nml_file), exist = exists)
    if (.not. exists) then
        errmsg = 'SCM data tendency :: namelist file: '//trim(nml_file)//' does not exist'
        errflg = 1
        return
    else
        open (unit = nlunit, file = nml_file, action = 'read', status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read (nlunit, nml = scm_data_nml)
    close (nlunit)

    ! ######################################################################################
    ! 
    ! Error checking
    !
    ! ######################################################################################
    ! Only proceed if scheme simulator requested.
    if (proc_SWRAD_config(1) .or. proc_LWRAD_config(1) .or. proc_PBL_config(1)  .or.       &
         proc_GWD_config(1)  .or. proc_SCNV_config(1)  .or. proc_DCNV_config(1) .or.       &
         proc_cldMP_config(1)) then
       do_ccpp_scheme_simulator = .true.
    else
       return
    endif
    
    ! Check that input data file exists
    inquire (file = trim (fileIN), exist = exists)
    if (.not. exists) then
       errmsg = 'SCM data tendency file: '//trim(fileIN)//' does not exist'
       errflg = 1
       return
    endif

    ! #######################################################################################
    !
    ! Read mandatory information from data file...
    ! (ONLY master processor(0), if MPI enabled)
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif

       ! Open file (required)
       status = nf90_open(trim(fileIN), NF90_NOWRITE, ncid)
       if (status /= nf90_noerr) then
          errmsg = 'Error reading in SCM data tendency file: '//trim(fileIN)
          errflg = 1
          return
       endif
       
       ! Get dimensions (required)
       status = nf90_inq_dimid(ncid, 'time', dimid)
       if (status == nf90_noerr) then
          status = nf90_inquire_dimension(ncid, dimid, len = ntime_data)
       else
          errmsg = 'SCM data tendency file: '//trim(fileIN)//' does not contain [time] dimension'
          errflg = 1
          return
       endif
       !
       status = nf90_inq_dimid(ncid, 'lev', dimid)
       if (status == nf90_noerr) then
          status = nf90_inquire_dimension(ncid, dimid, len = nlev_data)
       else
          errmsg = 'SCM data tendency file: '//trim(fileIN)//' does not contain [lev] dimension'
          errflg = 1
          return
       endif
#ifdef MPI
    endif ! On master processor

    ! Other processors waiting...
    call mpi_barrier(mpicomm, mpierr)

    ! #######################################################################################
    !
    ! Broadcast dimensions...
    ! (ALL processors)
    !
    ! #######################################################################################
    call mpi_bcast(ntime_data, 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nlev_data,  1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_barrier(mpicomm, mpierr)

    if (mpirank .eq. mpiroot) then
#endif

       ! ####################################################################################
       !
       ! What data fields do we have?
       !
       ! ####################################################################################

       !
       status = nf90_inq_varid(ncid, 'dT_dt_lwrad', varID)
       if (status == nf90_noerr) have_dTdt_LWRAD_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dT_dt_swrad', varID)
       if (status == nf90_noerr) have_dTdt_SWRAD_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dT_dt_pbl', varID)
       if (status == nf90_noerr) have_dTdt_PBL_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dq_dt_pbl', varID)
       if (status == nf90_noerr) have_dqdt_PBL_data = .true.
       !
       status = nf90_inq_varid(ncid, 'du_dt_pbl', varID)
       if (status == nf90_noerr) have_dudt_PBL_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dv_dt_pbl', varID)
       if (status == nf90_noerr) have_dvdt_PBL_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dT_dt_cgwd', varID)
       if (status == nf90_noerr) have_dTdt_GWD_data = .true.
       !
       status = nf90_inq_varid(ncid, 'du_dt_cgwd', varID)
       if (status == nf90_noerr) have_dudt_GWD_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dv_dt_cgwd', varID)
       if (status == nf90_noerr) have_dvdt_GWD_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dT_dt_shalconv', varID)
       if (status == nf90_noerr) have_dTdt_SCNV_data = .true.
       !
       status = nf90_inq_varid(ncid, 'du_dt_shalconv', varID)
       if (status == nf90_noerr) have_dudt_SCNV_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dv_dt_shalconv', varID)
       if (status == nf90_noerr) have_dvdt_SCNV_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dq_dt_shalconv', varID)
       if (status == nf90_noerr) have_dqdt_SCNV_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dT_dt_deepconv', varID)
       if (status == nf90_noerr) have_dTdt_DCNV_data = .true.
       !
       status = nf90_inq_varid(ncid, 'du_dt_deepconv', varID)
       if (status == nf90_noerr) have_dudt_DCNV_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dv_dt_deepconv', varID)
       if (status == nf90_noerr) have_dvdt_DCNV_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dq_dt_deepconv', varID)
       if (status == nf90_noerr) have_dqdt_DCNV_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dT_dt_micro', varID)
       if (status == nf90_noerr) have_dTdt_cldMP_data = .true.
       !
       status = nf90_inq_varid(ncid, 'dq_dt_micro', varID)
       if (status == nf90_noerr) have_dqdt_cldMP_data = .true.

#ifdef MPI
    endif ! Master process
#endif

    ! Allocate space for data
    allocate(time_data(ntime_data))
    if (have_dTdt_LWRAD_data) allocate(dTdt_LWRAD_data(nlev_data, ntime_data))
    if (have_dTdt_SWRAD_data) allocate(dTdt_SWRAD_data(nlev_data, ntime_data))
    if (have_dTdt_PBL_data)   allocate(dTdt_PBL_data(  nlev_data, ntime_data))
    if (have_dqdt_PBL_data)   allocate(dqdt_PBL_data(  nlev_data, ntime_data, nTrc))
    if (have_dudt_PBL_data)   allocate(dudt_PBL_data(  nlev_data, ntime_data))
    if (have_dvdt_PBL_data)   allocate(dvdt_PBL_data(  nlev_data, ntime_data))
    if (have_dTdt_GWD_data)   allocate(dTdt_GWD_data(  nlev_data, ntime_data))
    if (have_dudt_GWD_data)   allocate(dudt_GWD_data(  nlev_data, ntime_data))
    if (have_dvdt_GWD_data)   allocate(dvdt_GWD_data(  nlev_data, ntime_data))
    if (have_dTdt_SCNV_data)  allocate(dTdt_SCNV_data( nlev_data, ntime_data))
    if (have_dudt_SCNV_data)  allocate(dudt_SCNV_data( nlev_data, ntime_data))
    if (have_dvdt_SCNV_data)  allocate(dvdt_SCNV_data( nlev_data, ntime_data))
    if (have_dqdt_SCNV_data)  allocate(dqdt_SCNV_data( nlev_data, ntime_data, nTrc))
    if (have_dTdt_DCNV_data)  allocate(dTdt_DCNV_data( nlev_data, ntime_data))
    if (have_dudt_DCNV_data)  allocate(dudt_DCNV_data( nlev_data, ntime_data))
    if (have_dvdt_DCNV_data)  allocate(dvdt_DCNV_data( nlev_data, ntime_data))
    if (have_dqdt_DCNV_data)  allocate(dqdt_DCNV_data( nlev_data, ntime_data, nTrc))
    if (have_dTdt_cldMP_data) allocate(dTdt_cldMP_data(nlev_data, ntime_data))
    if (have_dqdt_cldMP_data) allocate(dqdt_cldMP_data(nlev_data, ntime_data, nTrc))

    ! #######################################################################################
    !
    ! Read in data ...
    ! (ONLY master processor(0), if MPI enabled) 
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif

       ! Temporal info (required)
       status = nf90_inq_varid(ncid, 'times', varID)
       if (status == nf90_noerr) then
          status = nf90_get_var(  ncid, varID, time_data)
       else
          errmsg = 'SCM data tendency file: '//trim(fileIN)//' does not contain times variable'
          errflg = 1
          return
       endif
       
       ! Read in physics data tendencies (optional)
       status = nf90_inq_varid(ncid, 'dT_dt_lwrad', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dTdt_LWRAD_data)
       !
       status = nf90_inq_varid(ncid, 'dT_dt_swrad', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dTdt_SWRAD_data)
       !
       status = nf90_inq_varid(ncid, 'dT_dt_pbl', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dTdt_PBL_data)
       !
       status = nf90_inq_varid(ncid, 'dq_dt_pbl', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dqdt_PBL_data)
       !
       status = nf90_inq_varid(ncid, 'du_dt_pbl', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dudt_PBL_data)
       !
       status = nf90_inq_varid(ncid, 'dv_dt_pbl', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dvdt_PBL_data)
       !
       status = nf90_inq_varid(ncid, 'dT_dt_cgwd', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dTdt_GWD_data)
       !
       status = nf90_inq_varid(ncid, 'du_dt_cgwd', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dudt_GWD_data)
       !
       status = nf90_inq_varid(ncid, 'dv_dt_cgwd', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dvdt_GWD_data)
       !
       status = nf90_inq_varid(ncid, 'dT_dt_shalconv', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dTdt_SCNV_data)
       !
       status = nf90_inq_varid(ncid, 'du_dt_shalconv', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dudt_SCNV_data)
       !
       status = nf90_inq_varid(ncid, 'dv_dt_shalconv', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dvdt_SCNV_data)
       !
       status = nf90_inq_varid(ncid, 'dq_dt_shalconv', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dqdt_SCNV_data)
       !
       status = nf90_inq_varid(ncid, 'dT_dt_deepconv', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dTdt_DCNV_data)
       !
       status = nf90_inq_varid(ncid, 'du_dt_deepconv', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dudt_DCNV_data)
       !
       status = nf90_inq_varid(ncid, 'dv_dt_deepconv', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dvdt_DCNV_data)
       !
       status = nf90_inq_varid(ncid, 'dq_dt_deepconv', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dqdt_DCNV_data)
       !
       status = nf90_inq_varid(ncid, 'dT_dt_micro', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dTdt_cldMP_data)
       !
       status = nf90_inq_varid(ncid, 'dq_dt_micro', varID)
       if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, dqdt_cldMP_data)
       !
       status = nf90_close(ncid)

#ifdef MPI
    endif ! Master process

    ! Other processors waiting...
    call mpi_barrier(mpicomm, mpierr)
    ! #######################################################################################
    !
    ! Broadcast data... 
    ! (ALL processors)
    !
    ! #######################################################################################

    if (have_dTdt_LWRAD_data) then
       call mpi_bcast(dTdt_LWRAD_data, size(dTdt_LWRAD_data), MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dTdt_SWRAD_data) then
       call mpi_bcast(dTdt_SWRAD_data, size(dTdt_SWRAD_data), MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dTdt_PBL_data) then
       call mpi_bcast(dTdt_PBL_data,   size(dTdt_PBL_data),   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dqdt_PBL_data) then
       call mpi_bcast(dqdt_PBL_data,   size(dqdt_PBL_data),   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dudt_PBL_data) then
       call mpi_bcast(dudt_PBL_data,   size(dudt_PBL_data),   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dvdt_PBL_data) then
       call mpi_bcast(dvdt_PBL_data,   size(dvdt_PBL_data),   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dTdt_GWD_data) then
       call mpi_bcast(dTdt_GWD_data,   size(dTdt_GWD_data),   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dudt_GWD_data) then
       call mpi_bcast(dudt_GWD_data,   size(dudt_GWD_data),   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dvdt_GWD_data) then
       call mpi_bcast(dvdt_GWD_data,   size(dvdt_GWD_data),   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dTdt_SCNV_data) then
       call mpi_bcast(dTdt_SCNV_data,  size(dTdt_SCNV_data),  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dudt_SCNV_data) then
       call mpi_bcast(dudt_SCNV_data,  size(dudt_SCNV_data),  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dvdt_SCNV_data) then
       call mpi_bcast(dvdt_SCNV_data,  size(dvdt_SCNV_data),  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dqdt_SCNV_data) then
       call mpi_bcast(dqdt_SCNV_data,  size(dqdt_SCNV_data),  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dTdt_DCNV_data) then
       call mpi_bcast(dTdt_DCNV_data,  size(dTdt_DCNV_data),  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dudt_DCNV_data) then
       call mpi_bcast(dudt_DCNV_data,  size(dudt_DCNV_data),  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dvdt_DCNV_data) then
       call mpi_bcast(dvdt_DCNV_data,  size(dvdt_DCNV_data),  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dqdt_DCNV_data) then
       call mpi_bcast(dqdt_DCNV_data,  size(dqdt_DCNV_data),  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dTdt_cldMP_data) then
       call mpi_bcast(dTdt_cldMP_data, size(dTdt_cldMP_data), MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (have_dqdt_cldMP_data) then
       call mpi_bcast(dqdt_cldMP_data, size(dqdt_cldMP_data), MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    !
    call mpi_barrier(mpicomm, mpierr)
#endif

    ! #######################################################################################
    !
    ! Populate physics_process type.
    !
    ! #######################################################################################

    ! Allocate
    allocate(physics_process(nPhysProcess))

    ! Metadata
    do iprc = 1,nPhysProcess
       if (iprc == proc_SWRAD_config(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "SWRAD"
          if (proc_SWRAD_config(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (proc_SWRAD_config(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == proc_LWRAD_config(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "LWRAD"
          if (proc_LWRAD_config(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (proc_LWRAD_config(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == proc_GWD_config(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "GWD"
          if (proc_GWD_config(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (proc_GWD_config(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == proc_PBL_config(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "PBL"
          if (proc_PBL_config(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (proc_PBL_config(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == proc_SCNV_config(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "SCNV"
          if (proc_SCNV_config(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (proc_SCNV_config(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == proc_DCNV_config(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "DCNV"
          if (proc_DCNV_config(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (proc_DCNV_config(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == proc_cldMP_config(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "cldMP"
          if (proc_cldMP_config(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (proc_cldMP_config(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
    enddo

    ! Load data
    if (have_dTdt_LWRAD_data) physics_process(proc_SWRAD_config(3))%tend%T => dTdt_LWRAD_data
    if (have_dTdt_SWRAD_data) physics_process(proc_LWRAD_config(3))%tend%T => dTdt_SWRAD_data
    if (have_dTdt_PBL_data)   physics_process(proc_PBL_config(3))%tend%T   => dTdt_PBL_data
    if (have_dudt_PBL_data)   physics_process(proc_PBL_config(3))%tend%u   => dudt_PBL_data
    if (have_dvdt_PBL_data)   physics_process(proc_PBL_config(3))%tend%v   => dvdt_PBL_data
    if (have_dqdt_PBL_data)   physics_process(proc_PBL_config(3))%tend%q   => dqdt_PBL_data
    if (have_dTdt_GWD_data)   physics_process(proc_GWD_config(3))%tend%T   => dTdt_GWD_data
    if (have_dudt_GWD_data)   physics_process(proc_GWD_config(3))%tend%u   => dudt_GWD_data
    if (have_dvdt_GWD_data)   physics_process(proc_GWD_config(3))%tend%v   => dvdt_GWD_data
    if (have_dTdt_SCNV_data)  physics_process(proc_SCNV_config(3))%tend%T  => dTdt_SCNV_data
    if (have_dudt_SCNV_data)  physics_process(proc_SCNV_config(3))%tend%u  => dudt_SCNV_data
    if (have_dvdt_SCNV_data)  physics_process(proc_SCNV_config(3))%tend%v  => dvdt_SCNV_data
    if (have_dqdt_SCNV_data)  physics_process(proc_SCNV_config(3))%tend%q  => dqdt_SCNV_data
    if (have_dTdt_DCNV_data)  physics_process(proc_DCNV_config(3))%tend%T  => dTdt_DCNV_data
    if (have_dudt_DCNV_data)  physics_process(proc_DCNV_config(3))%tend%u  => dudt_DCNV_data
    if (have_dvdt_DCNV_data)  physics_process(proc_DCNV_config(3))%tend%v  => dvdt_DCNV_data
    if (have_dqdt_DCNV_data)  physics_process(proc_DCNV_config(3))%tend%q  => dqdt_DCNV_data
    if (have_dTdt_cldMP_data) physics_process(proc_cldMP_config(3))%tend%T => dTdt_cldMP_data
    if (have_dqdt_cldMP_data) physics_process(proc_cldMP_config(3))%tend%q => dqdt_cldMP_data

    ! Which process-scheme is "Active"?
    do iprc = 1,nPhysProcess
       if (.not. physics_process(iprc)%use_sim) then
          iactive_scheme = iprc
       endif
    enddo

    !
    if (mpirank .eq. mpiroot) then
       print*, "----------------------------------"
       print*, "--- Using CCPP data tendencies ---"
       print*, "----------------------------------" 
       do iprc = 1,nPhysProcess
          if (physics_process(iprc)%use_sim) then
             print*,"  simulate_scheme: ", trim(physics_process(iprc)%name)
             print*,"      order:       ", physics_process(iprc)%order
             print*,"      time_split:  ", physics_process(iprc)%time_split
          endif
       enddo
       print*, "  active_scheme:   ", trim(physics_process(iactive_scheme)%name)
       print*, "      order:       ", physics_process(iactive_scheme)%order
       print*, "      time_split : ", physics_process(iactive_scheme)%time_split
       print*, "----------------------------------"
       print*, "----------------------------------"
    endif

  end subroutine ccpp_scheme_simulator_init

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
       errmsg, errflg)

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
    real(kind_phys), intent(inout), dimension(:,:,:) :: gq0
    character(len=*),intent(out) :: errmsg
    integer,         intent(out) :: errflg

    ! Locals
    integer :: iCol, iLay, iTrc, nCol, nLay,  nTrc, ti(1), tf(1), idtend, fcst_year,     &
         fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec, iprc, index_of_process
    real(kind_phys) :: w1, w2,hrofday
    real(kind_phys), dimension(:),     allocatable :: dT, du, dv, dq
    real(kind_phys), dimension(:,:),   allocatable :: gt1, gu1, gv1, dTdt, dudt, dvdt
    real(kind_phys), dimension(:,:,:), allocatable :: gq1, dqdt

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ccpp_scheme_simulator) return

    ! Current forecast time
    fcst_year  = jdat(1)
    fcst_month = jdat(2)
    fcst_day   = jdat(3)
    fcst_hour  = jdat(5)
    fcst_min   = jdat(6)
    fcst_sec   = jdat(7)

    ! Dimensions
    nCol = size(gq0(:,1,1))
    nLay = size(gq0(1,:,1))
    nTrc = size(gq0(1,1,:))

    ! Allocate temporaries
    allocate(gt1(nCol,nLay), gu1(nCol,nLay), gv1(nCol,nLay), gq1(nCol,nLay,1))
    allocate(dTdt(nCol,nLay), dudt(nCol,nLay), dvdt(nCol,nLay), dqdt(nCol,nLay,1))
    allocate(dT(nLay), du(nLay), dv(nLay), dq(nLay))

    ! Set state
    gt1(:,:)   = tgrs(:,:)
    gu1(:,:)   = ugrs(:,:)
    gv1(:,:)   = vgrs(:,:)
    gq1(:,:,1) = qgrs(:,:,1)
    dTdt(:,:)  = 0.
    dudt(:,:)  = 0.
    dvdt(:,:)  = 0.
    dqdt(:,:,1)= 0.

    ! Model internal physics timestep evolution of "state".
    do iprc = 1,nPhysProcess
       do iCol = 1,nCol
          !
          dT = 0.
          du = 0.
          dv = 0.
          dq = 0.

          ! Using scheme simulator (very simple, interpolate data tendency to local time)
          if (physics_process(iprc)%use_sim) then
             if (associated(physics_process(iprc)%tend%T)) then
                call linterp_data_tend("T", physics_process(iprc)%name, iprc, fcst_year, fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec, dT)
             endif
             if (associated(physics_process(iprc)%tend%u)) then
                call linterp_data_tend("u", physics_process(iprc)%name, iprc, fcst_year, fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec, du)
             endif
             if (associated(physics_process(iprc)%tend%v)) then
                call linterp_data_tend("v", physics_process(iprc)%name, iprc, fcst_year, fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec, dv)
             endif
             if (associated(physics_process(iprc)%tend%q)) then
                call linterp_data_tend("q", physics_process(iprc)%name, iprc, fcst_year, fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec, dq)
             endif

          ! Using data tendency from "active" scheme(s).
          ! DJS2023: This block is very ufs specific. Need to tidy this up.
          else
             if (physics_process(iprc)%name == "LWRAD") index_of_process = index_of_process_longwave
             if (physics_process(iprc)%name == "SWRAD") index_of_process = index_of_process_shortwave
             if (physics_process(iprc)%name == "PBL")   index_of_process = index_of_process_pbl
             if (physics_process(iprc)%name == "GWD")   index_of_process = index_of_process_orographic_gwd
             if (physics_process(iprc)%name == "SCNV")  index_of_process = index_of_process_scnv
             if (physics_process(iprc)%name == "DCNV")  index_of_process = index_of_process_dcnv
             if (physics_process(iprc)%name == "cldMP") index_of_process = index_of_process_mp
             !
             idtend = dtidx(index_of_temperature,index_of_process)
             if (idtend >= 1) dT = dtend(iCol,:,idtend)/dtp
             !
             idtend = dtidx(index_of_x_wind,index_of_process)
             if (idtend >= 1) du = dtend(iCol,:,idtend)/dtp
             !
             idtend = dtidx(index_of_y_wind,index_of_process)
             if (idtend >= 1) dv = dtend(iCol,:,idtend)/dtp
             !
             idtend = dtidx(100+ntqv,index_of_process)
             if (idtend >= 1) dq = dtend(iCol,:,idtend)/dtp
          endif

          ! Update state now?
          if (physics_process(iprc)%time_split) then
             gt1(iCol,:)    = gt1(iCol,:)   + (dTdt(iCol,:)   + dT)*dtp
             gu1(iCol,:)    = gu1(iCol,:)   + (dudt(iCol,:)   + du)*dtp
             gv1(iCol,:)    = gv1(iCol,:)   + (dvdt(iCol,:)   + dv)*dtp
             gq1(iCol,:,1)  = gq1(iCol,:,1) + (dqdt(iCol,:,1) + dq)*dtp
             dTdt(iCol,:)   = 0.
             dudt(iCol,:)   = 0.
             dvdt(iCol,:)   = 0.
             dqdt(iCol,:,1) = 0.
          ! Accumulate tendencies, update later?
          else
             dTdt(iCol,:)   = dTdt(iCol,:)   + dT
             dudt(iCol,:)   = dudt(iCol,:)   + du
             dvdt(iCol,:)   = dvdt(iCol,:)   + dv
             dqdt(iCol,:,1) = dqdt(iCol,:,1) + dq
          endif
       enddo
       !
       gt0(iCol,:)    = gt1(iCol,:)   + dTdt(iCol,:)*dtp
       gu0(iCol,:)    = gu1(iCol,:)   + dudt(iCol,:)*dtp
       gv0(iCol,:)    = gv1(iCol,:)   + dvdt(iCol,:)*dtp
       gq0(iCol,:,1)  = gq1(iCol,:,1) + dqdt(iCol,:,1)*dtp

    enddo
    !
  end subroutine ccpp_scheme_simulator_run

  ! ####################################################################################
  ! Utility functions/routines
  ! ####################################################################################
  ! The routine interpolates the data-tendencies
  subroutine linterp_data_tend(var_name, process_name, iprc, year, month, day, hour,   &
       minute, second, var_out)
    ! Inputs
    character(len=*), intent(in) :: var_name, process_name
    integer, intent(in) :: year, month, day, hour, minute, second, iprc
    
    ! Outputs
    real(kind_phys),dimension(:),intent(out) :: var_out

    ! Locals
    integer :: ti(1), tf(1)
    real(kind_phys) :: w1, w2, hrofday

    ! Linear interpolation weights
    hrofday = hour*3600. + minute*60. + second
    ti = findloc(abs(time_data-hrofday),minval(abs(time_data-hrofday)))
    if (hrofday - time_data(ti(1)) .le. 0) ti = ti-1
    tf = ti + 1
    w1 = (time_data(tf(1))-hrofday) / (time_data(tf(1)) - time_data(ti(1)))
    w2 = 1 - w1

    !
    select case(var_name)
       case("T")
          var_out = w1*physics_process(iprc)%tend%T(:,ti(1)) + w2*physics_process(iprc)%tend%T(:,tf(1))
       case("u")
          var_out = w1*physics_process(iprc)%tend%u(:,ti(1)) + w2*physics_process(iprc)%tend%u(:,tf(1))
       case("v")
          var_out = w1*physics_process(iprc)%tend%v(:,ti(1)) + w2*physics_process(iprc)%tend%v(:,tf(1))
       case("q")
          var_out = w1*physics_process(iprc)%tend%q(:,ti(1),1) + w2*physics_process(iprc)%tend%q(:,tf(1),1)
    end select

  end subroutine linterp_data_tend

end module ccpp_scheme_simulator

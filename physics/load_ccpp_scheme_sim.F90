! ########################################################################################
! 
! CCPP scheme to read and load data for ccpp_scheme_simulator
!
! ########################################################################################
module load_ccpp_scheme_sim
  use machine, only: kind_phys
  use netcdf
  use ccpp_scheme_simulator, only: do_ccpp_scheme_simulator, physics_process, active_name,&
       iactive_scheme, proc_start, proc_end, active_time_split_process
#ifdef MPI
  use mpi
#endif
  implicit none

  ! ########################################################################################
  !
  ! Configuration for CCPP scheme simulator. Set in namelist. Used during initialization to 
  ! populate "physics_process" type array, defined in ccpp_scheme_simulator.F90
  !
  ! ########################################################################################

  ! Number of physics process (set in namelist)
  integer :: nPhysProcess

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
  logical :: do_load_ccpp_scheme = .false.

  ! Data driven physics tendencies
  integer :: nlev_data, ntime_data
  real(kind_phys), allocatable, dimension(:), target   :: time_data
  real(kind_phys), allocatable, dimension(:,:), target :: dTdt_LWRAD_data,               &
       dTdt_SWRAD_data, dTdt_PBL_data, dudt_PBL_data, dvdt_PBL_data, dTdt_GWD_data,      &
       dudt_GWD_data, dvdt_GWD_data, dTdt_SCNV_data, dudt_SCNV_data, dvdt_SCNV_data,     &
       dTdt_DCNV_data, dudt_DCNV_data, dvdt_DCNV_data, dTdt_cldMP_data
  real(kind_phys), allocatable, dimension(:,:,:), target :: dqdt_PBL_data,               &
       dqdt_SCNV_data, dqdt_DCNV_data, dqdt_cldMP_data

  ! Scheme initialization flag.
  logical :: module_initialized = .false.

  public load_ccpp_scheme_sim_init
contains

  ! ######################################################################################
  !
  ! SUBROUTINE load_ccpp_scheme_sim_init
  !
  ! ######################################################################################
!! \section arg_table_load_ccpp_scheme_sim_init
!! \htmlinclude load_ccpp_scheme_sim_init.html
!!
  subroutine load_ccpp_scheme_sim_init(mpirank, mpiroot, mpicomm, nlunit, nml_file,          &
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
    ! Part A) Read in namelist and data.
    !
    ! ######################################################################################

    ! Read in namelist
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

    ! Read mandatory information from data file...
    ! (ONLY master processor(0), if MPI enabled)
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

    ! Broadcast dimensions...
    ! (ALL processors)
    call mpi_bcast(ntime_data, 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nlev_data,  1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_barrier(mpicomm, mpierr)

    if (mpirank .eq. mpiroot) then
#endif

       !
       ! What data fields do we have?
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

    ! Read in data ...
    ! (ONLY master processor(0), if MPI enabled) 
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

    ! Broadcast data... 
    ! (ALL processors)
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
    ! Part B) Populate physics_process type.
    !
    ! #######################################################################################
    ! Default process extent (no time-split physics processes)
    proc_start = 1
    proc_end   = nPhysProcess

    ! Allocate
    allocate(physics_process(nPhysProcess))

    ! Metadata
    do iprc = 1,nPhysProcess
       allocate(physics_process(iprc)%tend1d%T(nlev_data))
       allocate(physics_process(iprc)%tend1d%u(nlev_data))
       allocate(physics_process(iprc)%tend1d%v(nlev_data))
       allocate(physics_process(iprc)%tend1d%q(nlev_data,1))
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
    physics_process(proc_LWRAD_config(3))%tend2d%time => time_data
    physics_process(proc_SWRAD_config(3))%tend2d%time => time_data
    physics_process(proc_PBL_config(3))%tend2d%time   => time_data
    physics_process(proc_GWD_config(3))%tend2d%time   => time_data
    physics_process(proc_DCNV_config(3))%tend2d%time  => time_data
    physics_process(proc_SCNV_config(3))%tend2d%time  => time_data
    physics_process(proc_cldMP_config(3))%tend2d%time => time_data
    if (have_dTdt_LWRAD_data) physics_process(proc_SWRAD_config(3))%tend2d%T => dTdt_LWRAD_data
    if (have_dTdt_SWRAD_data) physics_process(proc_LWRAD_config(3))%tend2d%T => dTdt_SWRAD_data
    if (have_dTdt_PBL_data)   physics_process(proc_PBL_config(3))%tend2d%T   => dTdt_PBL_data
    if (have_dudt_PBL_data)   physics_process(proc_PBL_config(3))%tend2d%u   => dudt_PBL_data
    if (have_dvdt_PBL_data)   physics_process(proc_PBL_config(3))%tend2d%v   => dvdt_PBL_data
    if (have_dqdt_PBL_data)   physics_process(proc_PBL_config(3))%tend2d%q   => dqdt_PBL_data
    if (have_dTdt_GWD_data)   physics_process(proc_GWD_config(3))%tend2d%T   => dTdt_GWD_data
    if (have_dudt_GWD_data)   physics_process(proc_GWD_config(3))%tend2d%u   => dudt_GWD_data
    if (have_dvdt_GWD_data)   physics_process(proc_GWD_config(3))%tend2d%v   => dvdt_GWD_data
    if (have_dTdt_SCNV_data)  physics_process(proc_SCNV_config(3))%tend2d%T  => dTdt_SCNV_data
    if (have_dudt_SCNV_data)  physics_process(proc_SCNV_config(3))%tend2d%u  => dudt_SCNV_data
    if (have_dvdt_SCNV_data)  physics_process(proc_SCNV_config(3))%tend2d%v  => dvdt_SCNV_data
    if (have_dqdt_SCNV_data)  physics_process(proc_SCNV_config(3))%tend2d%q  => dqdt_SCNV_data
    if (have_dTdt_DCNV_data)  physics_process(proc_DCNV_config(3))%tend2d%T  => dTdt_DCNV_data
    if (have_dudt_DCNV_data)  physics_process(proc_DCNV_config(3))%tend2d%u  => dudt_DCNV_data
    if (have_dvdt_DCNV_data)  physics_process(proc_DCNV_config(3))%tend2d%v  => dvdt_DCNV_data
    if (have_dqdt_DCNV_data)  physics_process(proc_DCNV_config(3))%tend2d%q  => dqdt_DCNV_data
    if (have_dTdt_cldMP_data) physics_process(proc_cldMP_config(3))%tend2d%T => dTdt_cldMP_data
    if (have_dqdt_cldMP_data) physics_process(proc_cldMP_config(3))%tend2d%q => dqdt_cldMP_data

    ! Which process-scheme is "Active"? Is it a time-split process?
    do iprc = 1,nPhysProcess
       if (.not. physics_process(iprc)%use_sim) then
          iactive_scheme = iprc
          active_name    = physics_process(iprc)%name
          if (physics_process(iprc)%time_split) then
             active_time_split_process = .true.
          endif
       endif
    enddo

    !
    if (mpirank .eq. mpiroot) then
       print*, "-----------------------------------"
       print*, "--- Using CCPP scheme simulator ---"
       print*, "-----------------------------------" 
       do iprc = 1,nPhysProcess
          if (physics_process(iprc)%use_sim) then
             print*,"  simulate_scheme: ", trim(physics_process(iprc)%name)
             print*,"      order:       ", physics_process(iprc)%order
             print*,"      time_split:  ", physics_process(iprc)%time_split
          endif
       enddo
       print*, "  active_scheme:   ", trim(active_name)
       print*, "      order:       ", physics_process(iactive_scheme)%order
       print*, "      time_split : ", active_time_split_process
       print*, "-----------------------------------"
       print*, "-----------------------------------"
    endif

  end subroutine load_ccpp_scheme_sim_init

end module load_ccpp_scheme_sim

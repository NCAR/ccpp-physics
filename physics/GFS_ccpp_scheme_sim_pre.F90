! ########################################################################################
! 
! Interstitial CCPP scheme to couple UFS physics to ccpp_scheme_simulator.
! ) _init: read and load data into type used by ccpp_scheme_simulator
! ) _run:  prepare GFS diagnostic physics tendencies for ccpp_scheme_simulator 
!
! ########################################################################################
module GFS_ccpp_scheme_sim_pre
  use machine, only: kind_phys
  use netcdf
  use module_ccpp_scheme_simulator, only: base_physics_process
#ifdef MPI
  use mpi
#endif
  implicit none

  public GFS_ccpp_scheme_sim_pre_init, GFS_ccpp_scheme_sim_pre_run
contains

  ! ######################################################################################
  !
  ! SUBROUTINE GFS_ccpp_scheme_sim_pre_init
  !
  ! ######################################################################################
!! \section arg_table_GFS_ccpp_scheme_sim_pre_init
!! \htmlinclude GFS_ccpp_scheme_sim_pre_init.html
!!
  subroutine GFS_ccpp_scheme_sim_pre_init(mpirank, mpiroot, mpicomm, do_ccpp_scheme_sim, &
       scheme_sim_data, nprg_active, nprc_sim, prc_LWRAD_cfg, prc_SWRAD_cfg, prc_PBL_cfg,&
       prc_GWD_cfg, prc_SCNV_cfg, prc_DCNV_cfg, prc_cldMP_cfg, active_name,              &
       iactive_scheme, active_time_split_process, physics_process, errmsg, errflg)

    ! Inputs
    integer,               intent (in) :: mpirank, mpiroot, mpicomm, nprg_active, nprc_sim
    logical,               intent (in) :: do_ccpp_scheme_sim
    character(len=256),    intent (in) :: scheme_sim_data
    integer, dimension(3), intent (in) :: prc_LWRAD_cfg, prc_SWRAD_cfg, prc_PBL_cfg,     &
         prc_GWD_cfg, prc_SCNV_cfg, prc_DCNV_cfg, prc_cldMP_cfg

    ! Outputs
    type(base_physics_process),intent(inout) :: physics_process(:)
    character(len=16),intent(inout) :: active_name(:)
    integer, intent(inout)          :: iactive_scheme(:)
    logical, intent(inout)          :: active_time_split_process(:)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! Local variables
    integer :: ncid, dimID, varID, status, nlon, nlat, ios, iprc, iactive
    logical :: exists

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

    ! Data driven physics tendencies
    integer :: nlev_data, ntime_data
    real(kind_phys), allocatable, dimension(:), target   :: time_data
    real(kind_phys), allocatable, dimension(:,:), target :: dTdt_LWRAD_data,             &
         dTdt_SWRAD_data, dTdt_PBL_data, dudt_PBL_data, dvdt_PBL_data, dTdt_GWD_data,    &
         dudt_GWD_data, dvdt_GWD_data, dTdt_SCNV_data, dudt_SCNV_data, dvdt_SCNV_data,   &
         dTdt_DCNV_data, dudt_DCNV_data, dvdt_DCNV_data, dTdt_cldMP_data
    real(kind_phys), allocatable, dimension(:,:), target :: dqdt_PBL_data,               &
         dqdt_SCNV_data, dqdt_DCNV_data, dqdt_cldMP_data

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ccpp_scheme_sim) return
    
    ! ######################################################################################
    !
    ! Part A) Read in data.
    !
    ! ######################################################################################
    
    ! Check that input data file exists
    inquire (file = trim (scheme_sim_data), exist = exists)
    if (.not. exists) then
       errmsg = 'SCM data tendency file: '//trim(scheme_sim_data)//' does not exist'
       errflg = 1
       return
    endif

    ! Read mandatory information from data file...
    ! (ONLY master processor(0), if MPI enabled)
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif

       ! Open file (required)
       status = nf90_open(trim(scheme_sim_data), NF90_NOWRITE, ncid)
       if (status /= nf90_noerr) then
          errmsg = 'Error reading in SCM data tendency file: '//trim(scheme_sim_data)
          errflg = 1
          return
       endif
       
       ! Get dimensions (required)
       status = nf90_inq_dimid(ncid, 'time', dimid)
       if (status == nf90_noerr) then
          status = nf90_inquire_dimension(ncid, dimid, len = ntime_data)
       else
          errmsg = 'SCM data tendency file: '//trim(scheme_sim_data)//' does not contain [time] dimension'
          errflg = 1
          return
       endif
       !
       status = nf90_inq_dimid(ncid, 'lev', dimid)
       if (status == nf90_noerr) then
          status = nf90_inquire_dimension(ncid, dimid, len = nlev_data)
       else
          errmsg = 'SCM data tendency file: '//trim(scheme_sim_data)//' does not contain [lev] dimension'
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
    if (have_dqdt_PBL_data)   allocate(dqdt_PBL_data(  nlev_data, ntime_data))
    if (have_dudt_PBL_data)   allocate(dudt_PBL_data(  nlev_data, ntime_data))
    if (have_dvdt_PBL_data)   allocate(dvdt_PBL_data(  nlev_data, ntime_data))
    if (have_dTdt_GWD_data)   allocate(dTdt_GWD_data(  nlev_data, ntime_data))
    if (have_dudt_GWD_data)   allocate(dudt_GWD_data(  nlev_data, ntime_data))
    if (have_dvdt_GWD_data)   allocate(dvdt_GWD_data(  nlev_data, ntime_data))
    if (have_dTdt_SCNV_data)  allocate(dTdt_SCNV_data( nlev_data, ntime_data))
    if (have_dudt_SCNV_data)  allocate(dudt_SCNV_data( nlev_data, ntime_data))
    if (have_dvdt_SCNV_data)  allocate(dvdt_SCNV_data( nlev_data, ntime_data))
    if (have_dqdt_SCNV_data)  allocate(dqdt_SCNV_data( nlev_data, ntime_data))
    if (have_dTdt_DCNV_data)  allocate(dTdt_DCNV_data( nlev_data, ntime_data))
    if (have_dudt_DCNV_data)  allocate(dudt_DCNV_data( nlev_data, ntime_data))
    if (have_dvdt_DCNV_data)  allocate(dvdt_DCNV_data( nlev_data, ntime_data))
    if (have_dqdt_DCNV_data)  allocate(dqdt_DCNV_data( nlev_data, ntime_data))
    if (have_dTdt_cldMP_data) allocate(dTdt_cldMP_data(nlev_data, ntime_data))
    if (have_dqdt_cldMP_data) allocate(dqdt_cldMP_data(nlev_data, ntime_data))

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
          errmsg = 'SCM data tendency file: '//trim(scheme_sim_data)//' does not contain times variable'
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

    ! Metadata
    do iprc = 1,nprc_sim
       if (iprc == prc_SWRAD_cfg(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "SWRAD"
          if (prc_SWRAD_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (prc_SWRAD_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == prc_LWRAD_cfg(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "LWRAD"
          if (prc_LWRAD_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (prc_LWRAD_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == prc_GWD_cfg(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "GWD"
          if (prc_GWD_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (prc_GWD_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == prc_PBL_cfg(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "PBL"
          if (prc_PBL_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (prc_PBL_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == prc_SCNV_cfg(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "SCNV"
          if (prc_SCNV_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (prc_SCNV_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == prc_DCNV_cfg(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "DCNV"
          if (prc_DCNV_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (prc_DCNV_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
       if (iprc == prc_cldMP_cfg(3)) then
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "cldMP"
          if (prc_cldMP_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          endif
          if (prc_cldMP_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
       endif
    enddo

    ! Load data
    physics_process(prc_LWRAD_cfg(3))%tend2d%time => time_data
    physics_process(prc_SWRAD_cfg(3))%tend2d%time => time_data
    physics_process(prc_PBL_cfg(3))%tend2d%time   => time_data
    physics_process(prc_GWD_cfg(3))%tend2d%time   => time_data
    physics_process(prc_DCNV_cfg(3))%tend2d%time  => time_data
    physics_process(prc_SCNV_cfg(3))%tend2d%time  => time_data
    physics_process(prc_cldMP_cfg(3))%tend2d%time => time_data
    if (have_dTdt_LWRAD_data) physics_process(prc_SWRAD_cfg(3))%tend2d%T => dTdt_LWRAD_data
    if (have_dTdt_SWRAD_data) physics_process(prc_LWRAD_cfg(3))%tend2d%T => dTdt_SWRAD_data
    if (have_dTdt_PBL_data)   physics_process(prc_PBL_cfg(3))%tend2d%T   => dTdt_PBL_data
    if (have_dudt_PBL_data)   physics_process(prc_PBL_cfg(3))%tend2d%u   => dudt_PBL_data
    if (have_dvdt_PBL_data)   physics_process(prc_PBL_cfg(3))%tend2d%v   => dvdt_PBL_data
    if (have_dqdt_PBL_data)   physics_process(prc_PBL_cfg(3))%tend2d%q   => dqdt_PBL_data
    if (have_dTdt_GWD_data)   physics_process(prc_GWD_cfg(3))%tend2d%T   => dTdt_GWD_data
    if (have_dudt_GWD_data)   physics_process(prc_GWD_cfg(3))%tend2d%u   => dudt_GWD_data
    if (have_dvdt_GWD_data)   physics_process(prc_GWD_cfg(3))%tend2d%v   => dvdt_GWD_data
    if (have_dTdt_SCNV_data)  physics_process(prc_SCNV_cfg(3))%tend2d%T  => dTdt_SCNV_data
    if (have_dudt_SCNV_data)  physics_process(prc_SCNV_cfg(3))%tend2d%u  => dudt_SCNV_data
    if (have_dvdt_SCNV_data)  physics_process(prc_SCNV_cfg(3))%tend2d%v  => dvdt_SCNV_data
    if (have_dqdt_SCNV_data)  physics_process(prc_SCNV_cfg(3))%tend2d%q  => dqdt_SCNV_data
    if (have_dTdt_DCNV_data)  physics_process(prc_DCNV_cfg(3))%tend2d%T  => dTdt_DCNV_data
    if (have_dudt_DCNV_data)  physics_process(prc_DCNV_cfg(3))%tend2d%u  => dudt_DCNV_data
    if (have_dvdt_DCNV_data)  physics_process(prc_DCNV_cfg(3))%tend2d%v  => dvdt_DCNV_data
    if (have_dqdt_DCNV_data)  physics_process(prc_DCNV_cfg(3))%tend2d%q  => dqdt_DCNV_data
    if (have_dTdt_cldMP_data) physics_process(prc_cldMP_cfg(3))%tend2d%T => dTdt_cldMP_data
    if (have_dqdt_cldMP_data) physics_process(prc_cldMP_cfg(3))%tend2d%q => dqdt_cldMP_data

    ! Which process-scheme(s) is(are) "Active"? Are they time-split process?
    iactive = 0
    active_time_split_process(:) = .false.
    do iprc = 1,nprc_sim
       if (.not. physics_process(iprc)%use_sim) then
          iactive = iactive + 1
          iactive_scheme(iactive) = iprc
          active_name(iactive)    = physics_process(iprc)%name
          if (physics_process(iprc)%time_split) then
             active_time_split_process(iactive) = .true.
          endif
       endif
    enddo

    !
    if (mpirank .eq. mpiroot) then
       print*, "-----------------------------------"
       print*, "--- Using CCPP scheme simulator ---"
       print*, "-----------------------------------" 
       iactive = 1
       do iprc = 1,nprc_sim
          if (physics_process(iprc)%use_sim) then
             print*,"  simulate_scheme: ", trim(physics_process(iprc)%name)
             print*,"      order:       ", physics_process(iprc)%order
             print*,"      time_split:  ", physics_process(iprc)%time_split
          else
             print*, "  active_scheme:   ", trim(active_name(iactive))
             print*, "      order:       ", physics_process(iactive_scheme(iactive))%order
             print*, "      time_split : ", active_time_split_process(iactive)
             iactive = iactive + 1
          endif
       enddo
       print*, "-----------------------------------"
       print*, "-----------------------------------"
    endif

  end subroutine GFS_ccpp_scheme_sim_pre_init

  ! ######################################################################################
  !
  ! SUBROUTINE GFS_ccpp_scheme_sim_pre_run
  !
  ! ######################################################################################
!! \section arg_table_GFS_ccpp_scheme_sim_pre_run
!! \htmlinclude GFS_ccpp_scheme_sim_pre_run.html
!! 
  subroutine GFS_ccpp_scheme_sim_pre_run(dtend, ntqv, dtidx, dtp, index_of_process_dcnv, &
       index_of_process_longwave, index_of_process_shortwave, index_of_process_scnv,     &
       index_of_process_orographic_gwd, index_of_process_pbl, index_of_process_mp,       &
       index_of_temperature, index_of_x_wind, index_of_y_wind,                           &
       active_name, iactive_scheme_inloop, active_phys_tend, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: ntqv, index_of_process_dcnv, index_of_process_longwave,       &
         index_of_process_shortwave, index_of_process_scnv,                              &
         index_of_process_orographic_gwd, index_of_process_pbl, index_of_process_mp,     &
         index_of_temperature, index_of_x_wind, index_of_y_wind, iactive_scheme_inloop
    integer, intent(in), dimension(:,:) :: dtidx
    real(kind_phys), intent(in) :: dtp
    real(kind_phys), intent(in), dimension(:,:,:) :: dtend
    character(len=16),intent(in), dimension(:) :: active_name

    ! Outputs
    real(kind_phys), intent(out) :: active_phys_tend(:,:,:)
    character(len=*),intent(out) :: errmsg
    integer,         intent(out) :: errflg

    ! Locals
    integer :: idtend, iactive

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Get tendency for "active" process.

    ! ######################################################################################
    ! DJS2023: For the UFS and SCM, the physics tendencies are stored in a multi-dimensional
    ! array, CCPP standard_name = cumulative_change_of_state_variables.
    ! These are not the instantaneous physics tendencies that are applied to the state by 
    ! the physics schemes. Not all schemes output physics tendencies...
    ! Rather these are intended for diagnostic puposes and are accumulated over some 
    ! interval.
    ! In the UFS/SCM this is controlled by the diagnostic bucket interval, namelist option 
    ! "fhzero". For this to work, you need to clear the diagnostic buckets after each
    ! physics timestep when running in the UFS/SCM.
    !
    ! In the SCM this is done by adding the following runtime options:
    ! --n_itt_out 1 --n_itt_diag 1
    !
    ! ######################################################################################
    if (active_name(iactive_scheme_inloop) == "LWRAD") iactive = index_of_process_longwave
    if (active_name(iactive_scheme_inloop) == "SWRAD") iactive = index_of_process_shortwave
    if (active_name(iactive_scheme_inloop) == "PBL")   iactive = index_of_process_pbl
    if (active_name(iactive_scheme_inloop) == "GWD")   iactive = index_of_process_orographic_gwd
    if (active_name(iactive_scheme_inloop) == "SCNV")  iactive = index_of_process_scnv
    if (active_name(iactive_scheme_inloop) == "DCNV")  iactive = index_of_process_dcnv
    if (active_name(iactive_scheme_inloop) == "cldMP") iactive = index_of_process_mp

    ! Heat
    idtend = dtidx(index_of_temperature,iactive)
    if (idtend >= 1) active_phys_tend(:,:,1) = dtend(:,:,idtend)/dtp
    ! u-wind
    idtend = dtidx(index_of_x_wind,iactive)
    if (idtend >= 1) active_phys_tend(:,:,2) = dtend(:,:,idtend)/dtp
    ! v-wind
    idtend = dtidx(index_of_y_wind,iactive)
    if (idtend >= 1) active_phys_tend(:,:,3) = dtend(:,:,idtend)/dtp
    ! Moisture
    idtend = dtidx(100+ntqv,iactive)
    if (idtend >= 1) active_phys_tend(:,:,4) = dtend(:,:,idtend)/dtp


  end subroutine GFS_ccpp_scheme_sim_pre_run

end module GFS_ccpp_scheme_sim_pre

! ########################################################################################
! 
! Description: Interstitial CCPP suite to couple UFS physics to ccpp_suite_simulator.
!
! Contains:
! - load_ccpp_suite_sim(): read and load data into type used by ccpp_suite_simulator.
!      called once during model initialization
! - GFS_ccpp_suite_sim_pre_run(): prepare GFS diagnostic physics tendencies for 
!      ccpp_suite_simulator. 
!
! ########################################################################################
module GFS_ccpp_suite_sim_pre
  use machine, only: kind_phys
  use module_ccpp_suite_simulator, only: base_physics_process
  use netcdf
  implicit none
  public GFS_ccpp_suite_sim_pre_run, load_ccpp_suite_sim
contains

  ! ######################################################################################
  !
  ! SUBROUTINE GFS_ccpp_suite_sim_pre_run
  !
  ! ######################################################################################
!! \section arg_table_GFS_ccpp_suite_sim_pre_run
!! \htmlinclude GFS_ccpp_suite_sim_pre_run.html
!! 
  subroutine GFS_ccpp_suite_sim_pre_run(do_ccpp_suite_sim, dtend, ntqv, dtidx, dtp,      &
       index_of_process_dcnv, index_of_process_longwave, index_of_process_shortwave,     &
       index_of_process_scnv, index_of_process_orographic_gwd, index_of_process_pbl,     &
       index_of_process_mp, index_of_temperature, index_of_x_wind, index_of_y_wind,      &
       physics_process, iactive_T, iactive_u, iactive_v, iactive_q, active_phys_tend,    &
       errmsg, errflg)

    ! Inputs
    logical, intent(in) :: do_ccpp_suite_sim
    integer, intent(in) :: ntqv, index_of_process_dcnv, index_of_process_longwave,       &
         index_of_process_shortwave, index_of_process_scnv,                              &
         index_of_process_orographic_gwd, index_of_process_pbl, index_of_process_mp,     &
         index_of_temperature, index_of_x_wind, index_of_y_wind
    integer, intent(in), dimension(:,:) :: dtidx
    real(kind_phys), intent(in) :: dtp
    real(kind_phys), intent(in), dimension(:,:,:), optional :: dtend
    type(base_physics_process),intent(in) :: physics_process(:)
    integer,         intent(in) :: iactive_T, iactive_u, iactive_v, iactive_q

    ! Outputs
    real(kind_phys), intent(out) :: active_phys_tend(:,:,:)
    character(len=*),intent(out) :: errmsg
    integer,         intent(out) :: errflg

    ! Locals
    integer :: idtend, iactive

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ccpp_suite_sim) return

    ! Get tendency for "active" process.

    ! ######################################################################################
    ! DJS2023: For the UFS and SCM, the physics tendencies are stored in a multi-dimensional
    ! array, CCPP standard_name = cumulative_change_of_state_variables.
    ! These are not the instantaneous physics tendencies that are applied to the state by 
    ! the physics suites. Not all suites output physics tendencies...
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
    if (physics_process(1)%active_name == "LWRAD") iactive = index_of_process_longwave
    if (physics_process(1)%active_name == "SWRAD") iactive = index_of_process_shortwave
    if (physics_process(1)%active_name == "PBL")   iactive = index_of_process_pbl
    if (physics_process(1)%active_name == "GWD")   iactive = index_of_process_orographic_gwd
    if (physics_process(1)%active_name == "SCNV")  iactive = index_of_process_scnv
    if (physics_process(1)%active_name == "DCNV")  iactive = index_of_process_dcnv
    if (physics_process(1)%active_name == "cldMP") iactive = index_of_process_mp

    ! Heat
    idtend = dtidx(index_of_temperature,iactive)
    if (idtend >= 1) then
       active_phys_tend(:,:,iactive_T) = dtend(:,:,idtend)/dtp
    endif

    ! u-wind
    idtend = dtidx(index_of_x_wind,iactive)
    if (idtend >= 1) then
       active_phys_tend(:,:,iactive_u) = dtend(:,:,idtend)/dtp
    endif

    ! v-wind
    idtend = dtidx(index_of_y_wind,iactive)
    if (idtend >= 1) then
       active_phys_tend(:,:,iactive_v) = dtend(:,:,idtend)/dtp
    endif

    ! Moisture
    idtend = dtidx(100+ntqv,iactive)
    if (idtend >= 1) then
       active_phys_tend(:,:,iactive_q) = dtend(:,:,idtend)/dtp
    endif

  end subroutine GFS_ccpp_suite_sim_pre_run

  ! ######################################################################################
  subroutine load_ccpp_suite_sim(nlunit, nml_file, physics_process, iactive_T,           &
       iactive_u, iactive_v, iactive_q, errmsg, errflg)

    ! Inputs
    integer,          intent (in) :: nlunit
    character(len=*), intent (in) :: nml_file

    ! Outputs
    type(base_physics_process),intent(inout),allocatable :: physics_process(:)
    integer, intent(inout)        :: iactive_T, iactive_u, iactive_v, iactive_q
    integer, intent(out)          :: errflg
    character(len=256), intent(out) :: errmsg

    ! Local variables
    integer :: ncid, dimID, varID, status, ios, iprc, nlev_data, ntime_data
    character(len=256) :: suite_sim_file
    logical :: exists, do_ccpp_suite_sim
    integer :: nprc_sim

    ! For each process there is a corresponding namelist entry, which is constructed as 
    ! follows:
    ! {use_suite_sim[0(no)/1(yes)], time_split[0(no)/1(yes)], order[1:nPhysProcess]}
    integer, dimension(3) ::    &
         prc_LWRAD_cfg = (/0,0,0/), &
         prc_SWRAD_cfg = (/0,0,0/), &
         prc_PBL_cfg   = (/0,0,0/), &
         prc_GWD_cfg   = (/0,0,0/), &
         prc_SCNV_cfg  = (/0,0,0/), &
         prc_DCNV_cfg  = (/0,0,0/), &
         prc_cldMP_cfg = (/0,0,0/)

    ! Namelist
    namelist / ccpp_suite_sim_nml / do_ccpp_suite_sim, suite_sim_file, nprc_sim,         &
         prc_LWRAD_cfg, prc_SWRAD_cfg, prc_PBL_cfg, prc_GWD_cfg, prc_SCNV_cfg,           &
         prc_DCNV_cfg, prc_cldMP_cfg

    errmsg = ''
    errflg = 0

    ! Read in namelist
    inquire (file = trim (nml_file), exist = exists)
    if (.not. exists) then
        errmsg = 'CCPP suite simulator namelist file: '//trim(nml_file)//' does not exist.'
        errflg = 1
        return
    else
        open (unit = nlunit, file = nml_file, action = 'read', status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read (nlunit, nml = ccpp_suite_sim_nml, iostat=status)
    close (nlunit)

    ! Only proceed if suite simulator requested.
    if (prc_SWRAD_cfg(1)  == 1 .or. prc_LWRAD_cfg(1) == 1 .or. prc_PBL_cfg(1)  == 1 .or. &
         prc_GWD_cfg(1)   == 1 .or. prc_SCNV_cfg(1)  == 1 .or. prc_DCNV_cfg(1) == 1 .or. &
         prc_cldMP_cfg(1) == 1 ) then
    else
       return
    endif

    ! Check that input data file exists.
    inquire (file = trim (suite_sim_file), exist = exists)
    if (.not. exists) then
       errmsg = 'CCPP suite simulator file: '//trim(suite_sim_file)//' does not exist'
       errflg = 1
       return
    endif

    !
    ! Read data file...
    !

    ! Open file
    status = nf90_open(trim(suite_sim_file), NF90_NOWRITE, ncid)
    if (status /= nf90_noerr) then
       errmsg = 'Error reading in CCPP suite simulator file: '//trim(suite_sim_file)
       errflg = 1
       return
    endif

    ! Metadata (dimensions)
    status = nf90_inq_dimid(ncid, 'time', dimid)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dimid, len = ntime_data)
    else
       errmsg = 'CCPP suite simulator file: '//trim(suite_sim_file)//' does not contain [time] dimension'
       errflg = 1
       return
    endif
 
    status = nf90_inq_dimid(ncid, 'lev', dimid)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dimid, len = nlev_data)
    else
       errmsg = 'CCPP suite simulator file: '//trim(suite_sim_file)//' does not contain [lev] dimension'
       errflg = 1
       return
    endif

    ! Allocate space and read in data
    allocate(physics_process(nprc_sim))
    physics_process(1)%active_name    = ''
    physics_process(1)%iactive_scheme = 0
    physics_process(1)%active_tsp     = .false.
    do iprc = 1,nprc_sim
       allocate(physics_process(iprc)%tend1d%T(   nlev_data            ))
       allocate(physics_process(iprc)%tend1d%u(   nlev_data            ))
       allocate(physics_process(iprc)%tend1d%v(   nlev_data            ))
       allocate(physics_process(iprc)%tend1d%q(   nlev_data            ))
       allocate(physics_process(iprc)%tend2d%time(           ntime_data))
       allocate(physics_process(iprc)%tend2d%T(   nlev_data, ntime_data))
       allocate(physics_process(iprc)%tend2d%u(   nlev_data, ntime_data))
       allocate(physics_process(iprc)%tend2d%v(   nlev_data, ntime_data))
       allocate(physics_process(iprc)%tend2d%q(   nlev_data, ntime_data))

       ! Temporal info
       status = nf90_inq_varid(ncid, 'times', varID)
       if (status == nf90_noerr) then
          status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%time)
       else
          errmsg = 'SCM data tendency file: '//trim(suite_sim_file)//' does not contain times variable'
          errflg = 1
          return
       endif

       if (iprc == prc_SWRAD_cfg(3)) then
          ! Metadata
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "SWRAD"
          if (prc_SWRAD_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          else
             physics_process(1)%nprg_active = 1
             iactive_T   = 1
          endif
          if (prc_SWRAD_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif

          ! Data
          status = nf90_inq_varid(ncid, 'dT_dt_swrad', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%T)
       endif

       if (iprc == prc_LWRAD_cfg(3)) then
          ! Metadata
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "LWRAD"
          if (prc_LWRAD_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          else
             physics_process(1)%nprg_active = 1
             iactive_T   = 1
          endif
          if (prc_LWRAD_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif

          ! Data
          status = nf90_inq_varid(ncid, 'dT_dt_lwrad', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%T)
       endif

       if (iprc == prc_GWD_cfg(3)) then
          ! Metadata
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "GWD"
          if (prc_GWD_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          else
             physics_process(1)%nprg_active = 3
             iactive_T   = 1
             iactive_u   = 2
             iactive_v   = 3
          endif
          if (prc_GWD_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif

          ! Data
          status = nf90_inq_varid(ncid, 'dT_dt_cgwd', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%T)
          status = nf90_inq_varid(ncid, 'du_dt_cgwd', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%u)
          status = nf90_inq_varid(ncid, 'dv_dt_cgwd', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%v)
       endif

       if (iprc == prc_PBL_cfg(3)) then
          ! Metadata
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "PBL"
          if (prc_PBL_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          else
             physics_process(1)%nprg_active = 4
             iactive_T   = 1
             iactive_u   = 2
             iactive_v   = 3
             iactive_q   = 4
          endif
          if (prc_PBL_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif

          ! Data
          status = nf90_inq_varid(ncid, 'dT_dt_pbl', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%T)
          status = nf90_inq_varid(ncid, 'dq_dt_pbl', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%q)
          status = nf90_inq_varid(ncid, 'du_dt_pbl', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%u)
          status = nf90_inq_varid(ncid, 'dv_dt_pbl', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%v)
       endif

       if (iprc == prc_SCNV_cfg(3)) then
          ! Metadata
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "SCNV"
          if (prc_SCNV_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          else
             physics_process(1)%nprg_active = 4
             iactive_T   = 1
             iactive_u   = 2
             iactive_v   = 3
             iactive_q   = 4
          endif
          if (prc_SCNV_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif

          ! Data
          status = nf90_inq_varid(ncid, 'dT_dt_shalconv', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%T)
          status = nf90_inq_varid(ncid, 'du_dt_shalconv', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%u)
          status = nf90_inq_varid(ncid, 'dv_dt_shalconv', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%v)
          status = nf90_inq_varid(ncid, 'dq_dt_shalconv', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%q)
       endif

       if (iprc == prc_DCNV_cfg(3)) then
          ! Metadata
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "DCNV"
          if (prc_DCNV_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          else
             physics_process(1)%nprg_active = 4
             iactive_T   = 1
             iactive_u   = 2
             iactive_v   = 3
             iactive_q   = 4
          endif
          if (prc_DCNV_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif
          ! Data
          status = nf90_inq_varid(ncid, 'dT_dt_deepconv', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%T)
          status = nf90_inq_varid(ncid, 'du_dt_deepconv', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%u)
          status = nf90_inq_varid(ncid, 'dv_dt_deepconv', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%v)
          status = nf90_inq_varid(ncid, 'dq_dt_deepconv', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%q)
       endif

       if (iprc == prc_cldMP_cfg(3)) then
          ! Metadata
          physics_process(iprc)%order      = iprc
          physics_process(iprc)%name       = "cldMP"
          if (prc_cldMP_cfg(1) == 1) then
             physics_process(iprc)%use_sim = .true.
          else
             physics_process(1)%nprg_active = 2
             iactive_T   = 1
             iactive_q   = 2
          endif
          if (prc_cldMP_cfg(2) == 1) then
             physics_process(iprc)%time_split = .true.
          endif

          ! Data
          status = nf90_inq_varid(ncid, 'dT_dt_micro', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%T)
          status = nf90_inq_varid(ncid, 'dq_dt_micro', varID)
          if (status == nf90_noerr) status = nf90_get_var(  ncid, varID, physics_process(iprc)%tend2d%q)
       endif

       ! Which process-suite is "active"? Is process time-split?
       if (.not. physics_process(iprc)%use_sim) then
          physics_process(1)%iactive_scheme = iprc
          physics_process(1)%active_name    = physics_process(iprc)%name
          if (physics_process(iprc)%time_split) then
             physics_process(1)%active_tsp = .true.
          endif
       endif

    enddo

    if (physics_process(1)%iactive_scheme == 0) then
       errflg = 1
       errmsg = "ERROR: No active suite set for CCPP suite simulator"
       return
    endif

    print*, "-----------------------------------"
    print*, "--- Using CCPP suite simulator ---"
    print*, "-----------------------------------"
    do iprc = 1,nprc_sim
       if (physics_process(iprc)%use_sim) then
          print*,"  simulate_suite: ", trim(physics_process(iprc)%name)
          print*,"      order:       ", physics_process(iprc)%order
          print*,"      time_split:  ", physics_process(iprc)%time_split
       else
          print*, "  active_suite:   ", trim(physics_process(1)%active_name)
          print*, "      order:       ", physics_process(physics_process(1)%iactive_scheme)%order
          print*, "      time_split : ", physics_process(1)%active_tsp
       endif
    enddo
    print*, "-----------------------------------"
    print*, "-----------------------------------"

  end subroutine load_ccpp_suite_sim

end module GFS_ccpp_suite_sim_pre

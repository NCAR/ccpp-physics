! ########################################################################################
! 
! CCPP scheme to replace physics schemes with simulated data tendencies.
!
! ########################################################################################
module ccpp_scheme_ccpp_scheme_simulator
  use machine, only: kind_phys
  use netcdf
  implicit none

  !
  ! Data driven phsyics tendencies
  !
  real(kind_phys), allocatable, dimension(:)   :: time_data
  real(kind_phys), allocatable, dimension(:,:) :: dTdt_LWRAD_data, dTdt_SWRAD_data,      &
       dTdt_PBL_data, dudt_PBL_data, dvdt_PBL_data, dTdt_GWD_data, dudt_GWD_data,        &
       dvdt_GWD_data, dTdt_SCNV_data, dudt_SCNV_data, dvdt_SCNV_data, dTdt_DCNV_data,    &
       dudt_DCNV_data, dvdt_DCNV_data, dTdt_cldMP_data
  real(kind_phys), allocatable, dimension(:,:,:) :: dqdt_PBL_data, dqdt_SCNV_data,       &
       dqdt_DCNV_data, dqdt_cldMP_data

  !
  ! Logical switches for CCPP scheme simulator(s)
  !
  logical :: use_RAD_scheme_sim       = .false., &
             use_PBL_scheme_sim       = .false., &
             use_GWD_scheme_sim       = .false., &
             use_SCNV_scheme_sim      = .false., &
             use_DCNV_scheme_sim      = .false., &
             use_cldMP_scheme_sim     = .false.
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
  logical :: do_ccpp_scheme_simulator = .false.

  ! Host-model initial time information
  integer :: init_year, init_month, init_day, init_hour, init_min, init_sec

  public ccpp_scheme_ccpp_scheme_simulator_init, ccpp_scheme_ccpp_scheme_simulator_run
contains

  ! ######################################################################################
  !
  ! SUBROUTINE ccpp_scheme_ccpp_scheme_simulator_init
  !
  ! ######################################################################################
!! \section arg_table_ccpp_scheme_ccpp_scheme_simulator_init
!! \htmlinclude ccpp_scheme_ccpp_scheme_simulator_init.html
!!
  subroutine ccpp_scheme_ccpp_scheme_simulator_init(me, master, nlunit, nml_file, idat, errmsg, errflg)

    ! Inputs
    integer,          intent (in) :: me, master, nlunit
    character(len=*), intent (in) :: nml_file
    integer,          intent (in), dimension(8) :: idat

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    integer :: ncid, dimID, varID, status, nlon, nlat, nlev, ntime, ios
    character(len=256) :: fileIN
    logical :: exists
    integer,parameter :: nTrc = 1 ! Only specific humodty for now, but preserve 3 dimensionality

    ! Namelist
    namelist / scm_data_nml / &
         fileIN, use_RAD_scheme_sim, use_PBL_scheme_sim, use_GWD_scheme_sim, use_SCNV_scheme_sim, &
         use_DCNV_scheme_sim, use_cldMP_scheme_sim

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Store model initialization time.
    init_year  = idat(1)
    init_month = idat(2)
    init_day   = idat(3)
    init_hour  = idat(5)
    init_min   = idat(6)
    init_sec   = idat(7)

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
    if (use_RAD_scheme_sim   .or. use_PBL_scheme_sim   .or. use_GWD_scheme_sim .or. &
         use_SCNV_scheme_sim  .or. use_DCNV_scheme_sim  .or. use_cldMP_scheme_sim) then
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
       status = nf90_inquire_dimension(ncid, dimid, len = ntime)
    else
       errmsg = 'SCM data tendency file: '//trim(fileIN)//' does not contain time dimension'
       errflg = 1
       return
    endif
    !
    status = nf90_inq_dimid(ncid, 'lev', dimid)
    if (status == nf90_noerr) then
       status = nf90_inquire_dimension(ncid, dimid, len = nlev)
    else
       errmsg = 'SCM data tendency file: '//trim(fileIN)//' does not contain lev dimension'
       errflg = 1
       return
    endif

    ! Temporal info (required)
    status = nf90_inq_varid(ncid, 'times', varID)
    if (status == nf90_noerr) then
       allocate(time_data(ntime))
       status = nf90_get_var(  ncid, varID, time_data)
    else
       errmsg = 'SCM data tendency file: '//trim(fileIN)//' does not contain times variable'
       errflg = 1
       return
    endif

    ! Read in physics data tendencies (optional)
    status = nf90_inq_varid(ncid, 'dT_dt_lwrad', varID)
    if (status == nf90_noerr) then
       allocate(dTdt_LWRAD_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_LWRAD_data)
       have_dTdt_LWRAD_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dT_dt_swrad', varID)
    if (status == nf90_noerr) then
       allocate(dTdt_SWRAD_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_SWRAD_data)
       have_dTdt_SWRAD_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dT_dt_pbl', varID)
    if (status == nf90_noerr) then
       allocate(dTdt_PBL_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_PBL_data)
       have_dTdt_PBL_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dq_dt_pbl', varID)
    if (status == nf90_noerr) then
       allocate(dqdt_PBL_data(nlev, ntime, nTrc))
       status = nf90_get_var(  ncid, varID, dqdt_PBL_data)
       have_dqdt_PBL_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'du_dt_pbl', varID)
    if (status == nf90_noerr) then
       allocate(dudt_PBL_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dudt_PBL_data)
       have_dudt_PBL_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dv_dt_pbl', varID)
    if (status == nf90_noerr) then
       allocate(dvdt_PBL_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dvdt_PBL_data)
       have_dvdt_PBL_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dT_dt_cgwd', varID)
    if (status == nf90_noerr) then
       allocate(dTdt_GWD_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_GWD_data)
       have_dTdt_GWD_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'du_dt_cgwd', varID)
    if (status == nf90_noerr) then
       allocate(dudt_GWD_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dudt_GWD_data)
       have_dudt_GWD_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dv_dt_cgwd', varID)
    if (status == nf90_noerr) then
       allocate(dvdt_GWD_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dvdt_GWD_data)
       have_dvdt_GWD_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dT_dt_shalconv', varID)
    if (status == nf90_noerr) then
       allocate(dTdt_SCNV_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_SCNV_data)
       have_dTdt_SCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'du_dt_shalconv', varID)
    if (status == nf90_noerr) then
       allocate(dudt_SCNV_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dudt_SCNV_data)
       have_dudt_SCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dv_dt_shalconv', varID)
    if (status == nf90_noerr) then
       allocate(dvdt_SCNV_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dvdt_SCNV_data)
       have_dvdt_SCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dq_dt_shalconv', varID)
    if (status == nf90_noerr) then
       allocate(dqdt_SCNV_data(nlev, ntime, nTrc))
       status = nf90_get_var(  ncid, varID, dqdt_SCNV_data)
       have_dqdt_SCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dT_dt_deepconv', varID)
    if (status == nf90_noerr) then
       allocate(dTdt_DCNV_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_DCNV_data)
       have_dTdt_DCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'du_dt_deepconv', varID)
    if (status == nf90_noerr) then
       allocate(dudt_DCNV_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dudt_DCNV_data)
       have_dudt_DCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dv_dt_deepconv', varID)
    if (status == nf90_noerr) then
       allocate(dvdt_DCNV_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dvdt_DCNV_data)
       have_dvdt_DCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dq_dt_deepconv', varID)
    if (status == nf90_noerr) then
       allocate(dqdt_DCNV_data(nlev, ntime, nTrc))
       status = nf90_get_var(  ncid, varID, dqdt_DCNV_data)
       have_dqdt_DCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dT_dt_micro', varID)
    if (status == nf90_noerr) then
       allocate(dTdt_cldMP_data(nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_cldMP_data)
       have_dTdt_cldMP_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dq_dt_micro', varID)
    if (status == nf90_noerr) then
       allocate(dqdt_cldMP_data(nlev, ntime, nTrc))
       status = nf90_get_var(  ncid, varID, dqdt_cldMP_data)
       have_dqdt_cldMP_data = .true.
    endif

    !
    if (me == 0) then
       print*, "--- Using SCM data tendencies ---"
       print*, "---------------------------------"
       print*, "                 "
       print*, "use_RAD_scheme_sim:   ", use_RAD_scheme_sim
       print*, "   dTdt_LWRAD_data:   ", have_dTdt_LWRAD_data
       print*, "   dTdt_SWRAD_data:   ", have_dTdt_SWRAD_data
       print*, "use_PBL_scheme_sim:   ", use_PBL_scheme_sim
       print*, "   dTdt_PBL_data:     ", have_dTdt_PBL_data
       print*, "   dqdt_PBL_data:     ", have_dqdt_PBL_data
       print*, "   dudt_PBL_data:     ", have_dudt_PBL_data
       print*, "   dvdt_PBL_data:     ", have_dvdt_PBL_data
       print*, "use_GWD_scheme_sim:   ", use_GWD_scheme_sim
       print*, "   dTdt_gwd_data:     ", have_dTdt_GWD_data
       print*, "   dudt_gwd_data:     ", have_dudt_GWD_data
       print*, "   dvdt_gwd_data:     ", have_dvdt_GWD_data
       print*, "use_SCNV_scheme_sim:  ", use_SCNV_scheme_sim
       print*, "   dTdt_SCNV_data:    ", have_dTdt_SCNV_data
       print*, "   dudt_SCNV_data:    ", have_dudt_SCNV_data
       print*, "   dvdt_SCNV_data:    ", have_dvdt_SCNV_data
       print*, "   dqdt_SCNV_data:    ", have_dqdt_SCNV_data
       print*, "use_DCNV_scheme_sim:  ", use_DCNV_scheme_sim
       print*, "   dTdt_DCNV_data:    ", have_dTdt_DCNV_data
       print*, "   dudt_DCNV_data:    ", have_dudt_DCNV_data
       print*, "   dvdt_DCNV_data:    ", have_dvdt_DCNV_data
       print*, "   dqdt_DCNV_data:    ", have_dqdt_DCNV_data
       print*, "use_cldMP_scheme_sim: ", use_cldMP_scheme_sim
       print*, "   dTdt_cldMP_data:   ", have_dTdt_cldMP_data
       print*, "   dqdt_cldMP_data:   ", have_dqdt_cldMP_data
       print*, "---------------------------------"
    endif

  end subroutine ccpp_scheme_ccpp_scheme_simulator_init

  ! ######################################################################################
  !
  ! SUBROUTINE ccpp_scheme_ccpp_scheme_simulator_run
  !
  ! ######################################################################################
!! \section arg_table_ccpp_scheme_ccpp_scheme_simulator_run
!! \htmlinclude ccpp_scheme_ccpp_scheme_simulator_run.html
!!
  subroutine ccpp_scheme_ccpp_scheme_simulator_run(solhr, kdt, dtp, jdat, tgrs, ugrs, vgrs, qgrs, dtidx, &
       dtend, index_of_process_dcnv, index_of_process_longwave,                          &
       index_of_process_shortwave, index_of_process_scnv,                                &
       index_of_process_orographic_gwd, index_of_process_pbl, index_of_process_mp,       &
       index_of_temperature, index_of_x_wind, index_of_y_wind, ntqv, gt0, gu0, gv0, gq0, &
       errmsg, errflg)

    ! Inputs
    integer,         intent(in   ) :: kdt
    integer,         intent (in), dimension(8) :: jdat
    real(kind_phys), intent(in   ) :: dtp, solhr
    real(kind_phys), intent(in   ), dimension(:,:) :: tgrs, ugrs, vgrs
    real(kind=kind_phys), dimension(:,:,:), intent(in) :: qgrs, dtend
    integer, intent(in) :: dtidx(:,:), index_of_process_dcnv, index_of_process_longwave, &
         index_of_process_shortwave, index_of_process_scnv,                              &
         index_of_process_orographic_gwd, index_of_process_pbl, index_of_process_mp,     &
         index_of_temperature, index_of_x_wind, index_of_y_wind, ntqv

    ! Outputs
    real(kind_phys), intent(inout), dimension(:,:) :: gt0, gu0, gv0
    real(kind_phys), intent(inout), dimension(:,:,:) :: gq0
    character(len=*),intent(out  ) :: errmsg
    integer,         intent(out  ) :: errflg

    ! Locals
    integer :: iCol, iLay, iTrc, nCol, nLay,  nTrc, ti(1), tf(1), idtend, fcst_year,     &
         fcst_month, fcst_day, fcst_hour, fcst_min, fcst_sec
    real(kind_phys) :: w1, w2,hrofday
    real(kind_phys), dimension(:,:),   allocatable :: gt1, gu1, gv1
    real(kind_phys), dimension(:,:,:), allocatable :: gq1

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
    allocate(gt1(nCol,nLay), gu1(nCol,nLay), gv1(nCol,nLay), gq1(nCol,nLay,1)) ! *only specific humidity to start (ntrc=1).

    ! Determine temporal interpolation weights for data-tendecies.
    ! DJS: The data tendencies have a temporal dimension, to capture the diurnal cycle, 
    ! which is needed for reasonable solar forcing.
    hrofday = fcst_hour*3600. + fcst_min*60. + fcst_sec
    ti = findloc(abs(time_data-hrofday),minval(abs(time_data-hrofday)))
    if (hrofday - time_data(ti(1)) .le. 0) ti = ti-1
    tf = ti + 1
    w1 = (time_data(tf(1))-hrofday) / (time_data(tf(1)) - time_data(ti(1)))
    w2 = 1 - w1

    do iCol = 1,nCol
       ! Set state
       gt1(iCol,:)   = tgrs(iCol,:)
       gu1(iCol,:)   = ugrs(iCol,:)
       gv1(iCol,:)   = vgrs(iCol,:)
       gq1(iCol,:,1) = qgrs(iCol,:,1)

       ! ###############################################################################
       ! Radiation
       ! ###############################################################################
       if (use_RAD_scheme_sim) then
          if (have_dTdt_LWRAD_data) then
             gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_LWRAD_data(:,ti(1)) + w2*dTdt_LWRAD_data(:,tf(1))) * dtp
          endif
          if (have_dTdt_SWRAD_data) then
             gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_SWRAD_data(:,ti(1)) + w2*dTdt_SWRAD_data(:,tf(1))) * dtp
          endif
       else
          idtend = dtidx(index_of_temperature,index_of_process_longwave)
          if (idtend >= 1) then
             gt1(iCol,:) = gt1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_temperature,index_of_process_shortwave)
          if (idtend >=1) then
             gt1(iCol,:) = gt1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
       endif

       ! ###############################################################################
       ! PBL
       ! ###############################################################################
       if (use_PBL_scheme_sim) then
          if (have_dTdt_PBL_data) then
             gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_PBL_data(:,ti(1)) + w2*(dTdt_PBL_data(:,tf(1)))) * dtp
          endif
          if (have_dudt_PBL_data) then
             gu1(iCol,:) = gu1(iCol,:) + (w1*dudt_PBL_data(:,ti(1)) + w2*(dudt_PBL_data(:,tf(1)))) * dtp
          endif
          if (have_dvdt_PBL_data) then
             gv1(iCol,:) = gv1(iCol,:) + (w1*dvdt_PBL_data(:,ti(1)) + w2*(dvdt_PBL_data(:,tf(1)))) * dtp
          endif
          if (have_dqdt_PBL_data) then
             gq1(iCol,:,1) = gq1(iCol,:,1) + (w1*dqdt_PBL_data(:,ti(1),1) + w2*(dqdt_PBL_data(:,tf(1),1))) * dtp
          endif
       else
          idtend = dtidx(index_of_temperature,index_of_process_pbl)
          if (idtend >= 1) then
             gt1(iCol,:) = gt1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_x_wind,index_of_process_pbl)
          if (idtend >= 1) then
             gu1(iCol,:) = gu1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_y_wind,index_of_process_pbl)
          if (idtend >= 1) then
             gv1(iCol,:) = gv1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(100+ntqv, index_of_process_pbl)
          if (idtend >= 1) then
             gq1(iCol,:,1) = gq1(iCol,:,1) + dtend(iCol,:,idtend)! * dtp
          endif
       endif

       ! ###############################################################################
       ! Gravity wave drag
       ! ###############################################################################
       if (use_GWD_scheme_sim) then
          if (have_dTdt_GWD_data) then
             gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_GWD_data(:,ti(1)) + w2*(dTdt_GWD_data(:,tf(1)))) * dtp
          endif
          if (have_dudt_GWD_data) then
             gu1(iCol,:) = gu1(iCol,:) + (w1*dudt_GWD_data(:,ti(1)) + w2*(dudt_GWD_data(:,tf(1)))) * dtp
          endif
          if (have_dvdt_GWD_data) then
             gv1(iCol,:) = gv1(iCol,:) + (w1*dvdt_GWD_data(:,ti(1)) + w2*(dvdt_GWD_data(:,tf(1)))) * dtp
          endif
       else
          idtend = dtidx(index_of_temperature,index_of_process_orographic_gwd)
          if (idtend >= 1) then
             gt1(iCol,:) = gt1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_x_wind,index_of_process_orographic_gwd)
          if (idtend >= 1) then
             gu1(iCol,:) = gu1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_y_wind,index_of_process_orographic_gwd)
          if (idtend >= 1) then
             gv1(iCol,:) = gv1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
       endif

       ! ###############################################################################
       ! Shallow convection
       ! ###############################################################################
       if (use_SCNV_scheme_sim) then
          if (have_dTdt_SCNV_data) then
             gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_SCNV_data(:,ti(1)) + w2*(dTdt_SCNV_data(:,tf(1)))) * dtp
          endif
          if (have_dudt_SCNV_data) then
             gu1(iCol,:) = gu1(iCol,:) + (w1*dudt_SCNV_data(:,ti(1)) + w2*(dudt_SCNV_data(:,tf(1)))) * dtp
          endif
          if (have_dvdt_SCNV_data) then
             gv1(iCol,:) = gv1(iCol,:) + (w1*dvdt_SCNV_data(:,ti(1)) + w2*(dvdt_SCNV_data(:,tf(1)))) * dtp
          endif
          if (have_dqdt_SCNV_data) then
             gq1(iCol,:,1) = gq1(iCol,:,1) + (w1*dqdt_SCNV_data(:,ti(1),1) + w2*(dqdt_SCNV_data(:,tf(1),1))) * dtp
          endif
       else
          idtend = dtidx(index_of_temperature,index_of_process_scnv)
          if (idtend >= 1) then
             gt1(iCol,:) = gt1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_x_wind,index_of_process_scnv)
          if (idtend >= 1) then
             gu1(iCol,:) = gu1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_y_wind,index_of_process_scnv)
          if (idtend >= 1) then
             gv1(iCol,:) = gv1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(100+ntqv,index_of_process_scnv)
          if (idtend >= 1) then
             gq1(iCol,:,1) = gq1(iCol,:,1) + dtend(iCol,:,idtend)! * dtp
          endif
       endif

       ! ###############################################################################
       ! Deep convection
       ! ###############################################################################
       if (use_DCNV_scheme_sim) then
          if (have_dTdt_DCNV_data) then
             gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_DCNV_data(:,ti(1)) + w2*(dTdt_DCNV_data(:,tf(1)) )) * dtp
          endif
          if (have_dudt_DCNV_data) then
             gu1(iCol,:) = gu1(iCol,:) + (w1*dudt_DCNV_data(:,ti(1)) + w2*(dudt_DCNV_data(:,tf(1)))) * dtp
          endif
          if (have_dvdt_DCNV_data) then
             gv1(iCol,:) = gv1(iCol,:) + (w1*dvdt_DCNV_data(:,ti(1)) + w2*(dvdt_DCNV_data(:,tf(1)) )) * dtp
          endif
          if (have_dqdt_DCNV_data) then
             gq1(iCol,:,1) = gq1(iCol,:,1) + (w1*dqdt_DCNV_data(:,ti(1),1) + w2*(dqdt_DCNV_data(:,tf(1),1))) * dtp
          endif
       else
          idtend = dtidx(index_of_temperature,index_of_process_dcnv)
          if (idtend >= 1) then
             gt1(iCol,:) = gt1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_x_wind,index_of_process_dcnv)
          if (idtend >= 1) then
             gu1(iCol,:) = gu1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(index_of_y_wind,index_of_process_dcnv)
          if (idtend >= 1) then
             gv1(iCol,:) = gv1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(100+ntqv,index_of_process_dcnv)
          if (idtend >= 1) then
             gq1(iCol,:,1) = gq1(iCol,:,1) + dtend(iCol,:,idtend)! * dtp
          endif
       endif

       ! ###############################################################################
       ! Cloud microphysics
       ! ###############################################################################
       if (use_cldMP_scheme_sim) then
          if (have_dTdt_cldMP_data) then
             gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_cldMP_data(:,ti(1)) + w2*(dTdt_cldMP_data(:,tf(1)))) * dtp
          endif
          if (have_dqdt_cldMP_data) then
             gq1(iCol,:,1) = gq1(iCol,:,1) + (w1*dqdt_cldMP_data(:,ti(1),1) + w2*(dqdt_cldMP_data(:,tf(1),1))) * dtp
          endif
       else
          idtend = dtidx(index_of_temperature,index_of_process_mp)
          if (idtend >= 1) then
             gt1(iCol,:) = gt1(iCol,:) + dtend(iCol,:,idtend)! * dtp
          endif
          idtend = dtidx(100+ntqv,index_of_process_mp)
          if (idtend >= 1) then
             gq1(iCol,:,1) = gq1(iCol,:,1) + dtend(iCol,:,idtend)! * dtp
          endif
       endif

    enddo ! columns
    !
  end subroutine ccpp_scheme_ccpp_scheme_simulator_run

end module ccpp_scheme_ccpp_scheme_simulator

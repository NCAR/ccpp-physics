module GFS_driver

  use machine,                  only: kind_phys
  use GFS_typedefs,             only: GFS_init_type,                       &
                                      GFS_statein_type, GFS_stateout_type, &
                                      GFS_sfcprop_type, GFS_coupling_type, &
                                      GFS_control_type, GFS_grid_type,     &
                                      GFS_tbd_type,     GFS_cldprop_type,  &
                                      GFS_radtend_type, GFS_diag_type
  use module_radiation_driver,  only: GFS_radiation_driver, radupdate
  use module_physics_driver,    only: GFS_physics_driver
  use module_radsw_parameters,  only: topfsw_type, sfcfsw_type
  use module_radlw_parameters,  only: topflw_type, sfcflw_type
  use funcphys,                 only: gfuncphys

  implicit none

  private

!--------------------------------------------------------------------------------
! GFS_init_type
!--------------------------------------------------------------------------------
!   This container is the minimum set of data required from the dycore/atmosphere
!   component to allow proper initialization of the GFS physics
!  
!   Type is defined in GFS_typedefs.F90
!--------------------------------------------------------------------------------
! type GFS_init_type
!   public
!   integer :: me                                !< my MPI-rank
!   integer :: master                            !< master MPI-rank
!   integer :: isc                               !< starting i-index for this MPI-domain
!   integer :: jsc                               !< starting j-index for this MPI-domain
!   integer :: nx                                !< number of points in i-dir for this MPI rank
!   integer :: ny                                !< number of points in j-dir for this MPI rank
!   integer :: levs                              !< number of vertical levels
!   integer :: cnx                               !< number of points in i-dir for this cubed-sphere face
!                                                !< equal to gnx for lat-lon grids
!   integer :: cny                               !< number of points in j-dir for this cubed-sphere face
!                                                !< equal to gny for lat-lon grids
!   integer :: gnx                               !< number of global points in x-dir (i) along the equator
!   integer :: gny                               !< number of global points in y-dir (j) along any meridian
!   integer :: nlunit                            !< fortran unit number for file opens
!   integer :: logunit                           !< fortran unit number for writing logfile
!   integer :: dt_dycore                         !< dynamics time step in seconds
!   integer :: dt_phys                           !< physics  time step in seconds
!   integer :: bdat(8)                           !< model begin date in GFS format   (same as idat)
!   integer :: cdat(8)                           !< model current date in GFS format (same as jdat)
!   !--- blocking data
!   integer, pointer :: blksz(:)                 !< for explicit data blocking
!                                                !< default blksz(1)=[nx*ny]
!   !--- ak/bk for pressure level calculations
!   integer, pointer :: ak(:)                    !< from surface (k=1) to TOA (k=levs)
!   integer, pointer :: bk(:)                    !< from surface (k=1) to TOA (k=levs)
!   !--- grid metrics
!   real(kind=kind_phys), pointer :: xlon(:,:)   !< column longitude for MPI rank
!   real(kind=kind_phys), pointer :: xlat(:,:)   !< column latitude  for MPI rank
!   real(kind=kind_phys), pointer :: area(:,:)   !< column area for length scale calculations
!
!   character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
!                                                 !< based on name location in array
!   character(len=65) :: fn_nml                  !< namelist filename
! end type GFS_init_type
!--------------------------------------------------------------------------------
    
!------------------
! Module parameters
!------------------

!----------------------------
! Module variable definitions
!----------------------------
  real(kind=kind_phys), parameter :: con_24  =   24.0_kind_phys
  real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
  real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
  real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys
  real(kind=kind_phys), parameter :: qmin    =    1.0e-10

  integer, allocatable :: blksz(:)

!----------------
! Public entities
!----------------
  public  GFS_initialize              !< GFS initialization routine
  public  GFS_time_vary_step          !< perform operations needed prior radiation or physics
  public  GFS_radiation_driver        !< radiation_driver (was grrad)
  public  GFS_physics_driver          !< physics_driver (was gbphys)
  public  GFS_stochastic_driver       !< stochastic physics


  CONTAINS
!*******************************************************************************************


!--------------
! GFS initialze
!--------------
  subroutine GFS_initialize (Model, Statein, Stateout, Sfcprop,    &
                             Coupling, Grid, Tbd, Cldprop, Radtend, & 
                             Diag, Init_parm)

    use module_microphysics, only: gsmconst
    use cldwat2m_micro,      only: ini_micro
    use aer_cloud,           only: aer_cloud_init
    use module_ras,          only: ras_init


    !--- interface variables
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(inout) :: Statein(:)
    type(GFS_stateout_type),  intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),  intent(inout) :: Coupling(:)
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),   intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),   intent(inout) :: Radtend(:)
    type(GFS_diag_type),      intent(inout) :: Diag(:)
    type(GFS_init_type),      intent(in)    :: Init_parm

    !--- local variables
    integer :: nb
    integer :: nblks
    integer :: ntrac
    real(kind=kind_phys), allocatable :: si(:)
    real(kind=kind_phys), parameter   :: p_ref = 101325.0d0


    nblks = size(Init_parm%blksz)
    ntrac = size(Init_parm%tracer_names)
    allocate (blksz(nblks))
    blksz(:) = Init_parm%blksz(:)

    !--- set control properties (including namelist read)
    call Model%init (Init_parm%nlunit, Init_parm%fn_nml,           &
                     Init_parm%me, Init_parm%master,               &
                     Init_parm%logunit, Init_parm%isc,             &
                     Init_parm%jsc, Init_parm%nx, Init_parm%ny,    &
                     Init_parm%levs, Init_parm%cnx, Init_parm%cny, &
                     Init_parm%gnx, Init_parm%gny,                 &
                     Init_parm%dt_dycore, Init_parm%dt_phys,       &
                     Init_parm%bdat, Init_parm%cdat,               &
                     Init_parm%tracer_names)

    call read_o3data  (Model%ntoz, Model%me, Model%master)
    call read_h2odata (Model%h2o_phys, Model%me, Model%master)

    do nb = 1,nblks
      call Statein  (nb)%create (Init_parm%blksz(nb), Model)
      call Stateout (nb)%create (Init_parm%blksz(nb), Model)
      call Sfcprop  (nb)%create (Init_parm%blksz(nb), Model)
      call Coupling (nb)%create (Init_parm%blksz(nb), Model)
      call Grid     (nb)%create (Init_parm%blksz(nb), Model)
      call Tbd      (nb)%create (Init_parm%blksz(nb), Model)
      call Cldprop  (nb)%create (Init_parm%blksz(nb), Model)
      call Radtend  (nb)%create (Init_parm%blksz(nb), Model)
      !--- internal representation of diagnostics
      call Diag     (nb)%create (Init_parm%blksz(nb), Model)
    enddo

    !--- populate the grid components
    call GFS_grid_populate (Grid, Init_parm%xlon, Init_parm%xlat, Init_parm%area)
     
    !--- read in and initialize ozone and water
    if (Model%ntoz > 0) then
      do nb = 1, nblks
        call setindxoz (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_o3, &
                        Grid(nb)%jindx2_o3, Grid(nb)%ddy_o3)
      enddo
    endif

    if (Model%h2o_phys) then
      do nb = 1, nblks
        call setindxh2o (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_h, &
                         Grid(nb)%jindx2_h, Grid(nb)%ddy_h)
      enddo
    endif

    !--- Call gfuncphys (funcphys.f) to compute all physics function tables.
    call gfuncphys ()

    call gsmconst (Model%dtp, Model%me, .TRUE.)

    !--- define sigma level for radiation initialization 
    !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
    !--- ps is replaced with p0. The value of p0 uses that in http://www.emc.ncep.noaa.gov/officenotes/newernotes/on461.pdf
    !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
    allocate(si(Model%levr+1))
    si = (Init_parm%ak + Init_parm%bk * p_ref - Init_parm%ak(Model%levr+1)) &
             / (p_ref - Init_parm%ak(Model%levr+1))
    call rad_initialize (si, Model%levr, Model%ictm, Model%isol, &
           Model%ico2, Model%iaer, Model%ialb, Model%iems,       &
           Model%ntcw, Model%num_p3d, Model%npdf3d, Model%ntoz,  &
           Model%iovr_sw, Model%iovr_lw, Model%isubc_sw,         &
           Model%isubc_lw, Model%crick_proof, Model%ccnorm,      &
           Model%norad_precip, Model%idate,Model%iflip, Model%me)
    deallocate (si)

    !--- initialize Morrison-Gettleman microphysics
    if (Model%ncld == 2) then
      call ini_micro (Model%mg_dcs, Model%mg_qcvar, Model%mg_ts_auto_ice)
      call aer_cloud_init ()
    endif

    !--- initialize ras
    if (Model%ras) call ras_init (Model%levs, Model%me)

    !--- initialize soil vegetation
    call set_soilveg(Model%me, Model%isot, Model%ivegsrc, Model%nlunit)

    !--- lsidea initialization
    if (Model%lsidea) then
      print *,' LSIDEA is active but needs to be reworked for FV3 - shutting down'
      stop
      !--- NEED TO get the logic from the old phys/gloopb.f initialization area
    endif

    !--- sncovr may not exist in ICs from chgres.
    !--- FV3GFS handles this as part of the IC ingest
    !--- this not is placed here to alert users to the need to study
    !--- the FV3GFS_io.F90 module

  end subroutine GFS_initialize


!-------------------------------------------------------------------------
! time_vary_step
!-------------------------------------------------------------------------
!    routine called prior to radiation and physics steps to handle:
!      1) sets up various time/date variables
!      2) sets up various triggers
!      3) defines random seed indices for radiation (in a reproducible way)
!      5) interpolates coefficients for prognostic ozone calculation
!      6) performs surface data cycling via the GFS gcycle routine
!-------------------------------------------------------------------------
  subroutine GFS_time_vary_step (Model, Statein, Stateout, Sfcprop, Coupling, & 
                                 Grid, Tbd, Cldprop, Radtend, Diag)

    implicit none

    !--- interface variables
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(inout) :: Statein(:)
    type(GFS_stateout_type),  intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),  intent(inout) :: Coupling(:)
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),   intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),   intent(inout) :: Radtend(:)
    type(GFS_diag_type),      intent(inout) :: Diag(:)
    !--- local variables
    integer :: nb, nblks
    real(kind=kind_phys) :: rinc(5)
    real(kind=kind_phys) :: sec

    nblks = size(blksz)
    !--- Model%jdat is being updated directly inside of FV3GFS_cap.F90
    !--- update calendars and triggers
    rinc(1:5)   = 0
    call w3difdat(Model%jdat,Model%idat,4,rinc)
    sec = rinc(4)
    Model%phour = sec/con_hr
    !--- set current bucket hour
    Model%zhour = Model%phour
    Model%fhour = (sec + Model%dtp)/con_hr
    Model%kdt   = nint((sec + Model%dtp)/Model%dtp)

    Model%ipt    = 1
    Model%lprnt  = .false.
    Model%lssav  = .true.

    !--- radiation triggers
    Model%lsswr  = (mod(Model%kdt, Model%nsswr) == 1)
    Model%lslwr  = (mod(Model%kdt, Model%nslwr) == 1)

    !--- set the solar hour based on a combination of phour and time initial hour
    Model%solhr  = mod(Model%phour+Model%idate(1),con_24)

    if ((Model%debug) .and. (Model%me == Model%master)) then
      print *,'   sec ', sec
      print *,'   kdt ', Model%kdt
      print *,' nsswr ', Model%nsswr
      print *,' nslwr ', Model%nslwr
      print *,' nscyc ', Model%nscyc
      print *,' lsswr ', Model%lsswr
      print *,' lslwr ', Model%lslwr
      print *,' fhour ', Model%fhour
      print *,' phour ', Model%phour
      print *,' solhr ', Model%solhr
    endif

    !--- radiation time varying routine
    if (Model%lsswr .or. Model%lslwr) then
      call GFS_rad_time_vary (Model, Statein, Tbd, sec)
    endif

    !--- physics time varying routine
    call GFS_phys_time_vary (Model, Grid, Tbd)

    !--- repopulate specific time-varying sfc properties for AMIP/forecast runs
    if (Model%nscyc >  0) then
      if (mod(Model%kdt,Model%nscyc) == 1) THEN
        call gcycle (nblks, Model, Grid(:), Sfcprop(:), Cldprop(:))
      endif
    endif

    !--- determine if diagnostics buckets need to be cleared
    if (mod(Model%kdt,Model%nszero) == 1) then
      do nb = 1,nblks
        call Diag(nb)%rad_zero  (Model)
        call Diag(nb)%phys_zero (Model)
    !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
      enddo
    endif

  end subroutine GFS_time_vary_step


!-------------------------------------------------------------------------
! GFS stochastic_driver
!-------------------------------------------------------------------------
!    routine called prior to radiation and physics steps to handle:
!      1) sets up various time/date variables
!      2) sets up various triggers
!      3) defines random seed indices for radiation (in a reproducible way)
!      5) interpolates coefficients for prognostic ozone calculation
!      6) performs surface data cycling via the GFS gcycle routine
!-------------------------------------------------------------------------
  subroutine GFS_stochastic_driver (Model, Statein, Stateout, Sfcprop, Coupling, &
                                    Grid, Tbd, Cldprop, Radtend, Diag)

    implicit none

    !--- interface variables
    type(GFS_control_type),   intent(in   ) :: Model
    type(GFS_statein_type),   intent(in   ) :: Statein
    type(GFS_stateout_type),  intent(in   ) :: Stateout
    type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
    type(GFS_coupling_type),  intent(inout) :: Coupling
    type(GFS_grid_type),      intent(in   ) :: Grid
    type(GFS_tbd_type),       intent(in   ) :: Tbd
    type(GFS_cldprop_type),   intent(in   ) :: Cldprop
    type(GFS_radtend_type),   intent(in   ) :: Radtend
    type(GFS_diag_type),      intent(inout) :: Diag
    !--- local variables
    integer :: k, i
    real(kind=kind_phys) :: upert, vpert, tpert, qpert, qnew

     if (Model%do_sppt) then
       do k = 1,size(Statein%tgrs,2)
         do i = 1,size(Statein%tgrs,1)
      
           upert = (Stateout%gu0(i,k)   - Statein%ugrs(i,k))   * Coupling%sppt_wts(i,k)
           vpert = (Stateout%gv0(i,k)   - Statein%vgrs(i,k))   * Coupling%sppt_wts(i,k)
           tpert = (Stateout%gt0(i,k)   - Statein%tgrs(i,k))   * Coupling%sppt_wts(i,k) - Tbd%dtdtr(i,k)
           qpert = (Stateout%gq0(i,k,1) - Statein%qgrs(i,k,1)) * Coupling%sppt_wts(i,k)
 
           Stateout%gu0(i,k)  = Statein%ugrs(i,k)+upert
           Stateout%gv0(i,k)  = Statein%vgrs(i,k)+vpert
 
           !negative humidity check
           qnew = Statein%qgrs(i,k,1)+qpert
           if (qnew .GE. 1.0e-10) then
              Stateout%gq0(i,k,1) = qnew
              Stateout%gt0(i,k)   = Statein%tgrs(i,k) + tpert + Tbd%dtdtr(i,k)
           endif
         enddo
       enddo
 
       Diag%totprcp(:)      = Diag%totprcp(:)      + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%dtotprcp(:)
       Diag%cnvprcp(:)      = Diag%cnvprcp(:)      + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%dcnvprcp(:)
       Coupling%rain_cpl(:) = Coupling%rain_cpl(:) + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%drain_cpl(:)
       Coupling%snow_cpl(:) = Coupling%snow_cpl(:) + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%dsnow_cpl(:)
     endif

     if (Model%do_shum) then
       Stateout%gq0(:,:,1) = Stateout%gq0(:,:,1)*(1.0 + Coupling%shum_wts(:,:))
     endif

  end subroutine GFS_stochastic_driver



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------
! GFS_rad_time_vary
!-----------------------------------------------------------------------
!
!  Routine containing all of the setup logic originally in phys/gloopr.f
!
!-----------------------------------------------------------------------
  subroutine GFS_rad_time_vary (Model, Statein, Tbd, sec)

    use physparam,        only: ipsd0, ipsdlim, iaerflg
    use mersenne_twister, only: random_setseed, random_index, random_stat

    implicit none

    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(in)    :: Statein(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    real(kind=kind_phys),     intent(in)    :: sec
    !--- local variables
    type (random_stat) :: stat
    integer :: ix, nb, j, i, nblks, ipseed
    integer :: numrdm(Model%cnx*Model%cny*2)

    nblks = size(blksz,1)

    call radupdate (Model%idat, Model%jdat, Model%fhswr, Model%dtf, Model%lsswr, &
                    Model%me, Model%slag, Model%sdec, Model%cdec, Model%solcon )

    !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
    if ((Model%isubc_lw==2) .or. (Model%isubc_sw==2)) then
      ipseed = mod(nint(con_100*sqrt(sec)), ipsdlim) + 1 + ipsd0
      call random_setseed (ipseed, stat)
      call random_index (ipsdlim, numrdm, stat)

      !--- set the random seeds for each column in a reproducible way
      ix = 0
      nb = 1
      do j = 1,Model%ny
        do i = 1,Model%nx
          ix = ix + 1
          if (ix .gt. blksz(nb)) then
            ix = 1
            nb = nb + 1
          endif
          !--- for testing purposes, replace numrdm with '100'
          Tbd(nb)%icsdsw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx)
          Tbd(nb)%icsdlw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + Model%cnx*Model%cny)
        enddo
      enddo
    endif  ! isubc_lw and isubc_sw

    if (Model%num_p3d == 4) then
      if (Model%kdt == 1) then
        do nb = 1,nblks
          Tbd(nb)%phy_f3d(:,:,1) = Statein(nb)%tgrs
          Tbd(nb)%phy_f3d(:,:,2) = max(qmin,Statein(nb)%qgrs(:,:,1))
          Tbd(nb)%phy_f3d(:,:,3) = Statein(nb)%tgrs
          Tbd(nb)%phy_f3d(:,:,4) = max(qmin,Statein(nb)%qgrs(:,:,1))
          Tbd(nb)%phy_f2d(:,1)   = Statein(nb)%prsi(:,1)
          Tbd(nb)%phy_f2d(:,2)   = Statein(nb)%prsi(:,1)
        enddo
      endif
    endif

  end subroutine GFS_rad_time_vary


!-----------------------------------------------------------------------
! GFS_phys_time_vary
!-----------------------------------------------------------------------
!
!  Routine containing all of the setup logic originally in phys/gloopb.f
!
!-----------------------------------------------------------------------
  subroutine GFS_phys_time_vary (Model, Grid, Tbd)
    use mersenne_twister, only: random_setseed, random_number

    implicit none
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    !--- local variables
    integer :: nb, ix, k, j, i, nblks, iseed, iskip
    real(kind=kind_phys) :: wrk(1)
    real(kind=kind_phys) :: rannie(Model%cny)
    real(kind=kind_phys) :: rndval(Model%cnx*Model%cny*Model%nrcm)

    nblks = size(blksz,1)

    !--- switch for saving convective clouds - cnvc90.f 
    !--- aka Ken Campana/Yu-Tai Hou legacy
    if ((mod(Model%kdt,Model%nsswr) == 0) .and. (Model%lsswr)) then
      !--- initialize,accumulate,convert
      Model%clstp = 1100 + min(Model%fhswr/con_hr,Model%fhour,con_99)
    elseif (mod(Model%kdt,Model%nsswr) == 0) then
      !--- accumulate,convert
      Model%clstp = 0100 + min(Model%fhswr/con_hr,Model%fhour,con_99)
    elseif (Model%lsswr) then
      !--- initialize,accumulate
      Model%clstp = 1100
    else
      !--- accumulate
      Model%clstp = 0100
    endif

    !--- random number needed for RAS and old SAS and when cal_pre=.true.
    if ( ((Model%imfdeepcnv <= 0) .or. (Model%cal_pre)) .and. (Model%random_clds) ) then
      iseed = mod(con_100*sqrt(Model%fhour*con_hr),1.0d9) + Model%seed0
      call random_setseed(iseed)
      call random_number(wrk)
      do i = 1,Model%cnx*Model%nrcm
        iseed = iseed + nint(wrk(1)) * i
        call random_setseed(iseed)
        call random_number(rannie)
        rndval(1+(i-1)*Model%cny:i*Model%cny) = rannie(1:Model%cny)
      enddo

      do k = 1,Model%nrcm
        iskip = (k-1)*Model%cnx*Model%cny
        ix = 0
        nb = 1
        do j = 1,Model%ny
          do i = 1,Model%nx
            ix = ix + 1
            if (ix .gt. blksz(nb)) then
              ix = 1
              nb = nb + 1
            endif
            Tbd(nb)%rann(ix,k) = rndval(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + iskip)
          enddo
        enddo
      enddo
    endif  ! imfdeepcnv, cal_re, random_clds

    !--- o3 interpolation
    if (Model%ntoz > 0) then
      do nb = 1, nblks
        call ozinterpol (Model%me, blksz(nb), Model%idate, Model%fhour, &
                         Grid(nb)%jindx1_o3, Grid(nb)%jindx2_o3,            &
                         Tbd(nb)%ozpl, Grid(nb)%ddy_o3)
      enddo
    endif

    !--- h2o interpolation
     if (Model%h2o_phys) then
       do nb = 1, nblks
         call h2ointerpol (Model%me, blksz(nb), Model%idate, Model%fhour, &
                           Grid(nb)%jindx1_h, Grid(nb)%jindx2_h,          &
                           Tbd(nb)%h2opl, Grid(nb)%ddy_h)
       enddo
     endif

  end subroutine GFS_phys_time_vary


!------------------
! GFS_grid_populate
!------------------
  subroutine GFS_grid_populate (Grid, xlon, xlat, area)
    use physcons,                 only: pi => con_pi

    implicit none

    type(GFS_grid_type)              :: Grid(:)
    real(kind=kind_phys), intent(in) :: xlon(:,:)
    real(kind=kind_phys), intent(in) :: xlat(:,:)
    real(kind=kind_phys), intent(in) :: area(:,:)

    !--- local variables
    integer :: nb, ix, blksz, i, j

    blksz = size(Grid(1)%xlon)

    nb = 1
    ix = 0
    do j = 1,size(xlon,2)
      do i = 1,size(xlon,1)
        ix=ix+1
        if (ix .gt. blksz) then
          nb = nb + 1
          ix = 1
        endif
        Grid(nb)%xlon(ix)   = xlon(i,j)
        Grid(nb)%xlat(ix)   = xlat(i,j)
        Grid(nb)%xlat_d(ix) = xlat(i,j) * 180.0_kind_phys/pi
        Grid(nb)%sinlat(ix) = sin(Grid(nb)%xlat(ix))
        Grid(nb)%coslat(ix) = sqrt(1.0_kind_phys - Grid(nb)%sinlat(ix)*Grid(nb)%sinlat(ix))
        Grid(nb)%area(ix)   = area(i,j)
        Grid(nb)%dx(ix)     = sqrt(area(i,j))
      enddo
    enddo

  end subroutine GFS_grid_populate

end module GFS_driver


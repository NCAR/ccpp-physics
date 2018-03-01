module GFS_driver_gmtb_scm
  use machine,                  only: kind_phys
  use GFS_typedefs,             only: GFS_init_type,                          &
                                      GFS_statein_type,  GFS_stateout_type,   &
                                      GFS_sfcprop_type,  GFS_coupling_type,   &
                                      GFS_control_type,  GFS_grid_type,       &
                                      GFS_tbd_type,      GFS_cldprop_type,    &
                                      GFS_radtend_type,  GFS_diag_type,       &
                                      GFS_sfccycle_type, GFS_interstitial_type

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
  public  GFS_ozone_and_h20_setup, GFS_initialize              !< GFS initialization routine

  CONTAINS
!*******************************************************************************************

  subroutine GFS_ozone_and_h20_setup(n_columns, n_levels, lats, pres, &
    n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times, &
    ozone_lat, ozone_pres, ozone_time, ozone_forcing_in)
    use ozne_def, only: latsozp, levozp, timeoz, oz_coeff, oz_lat, oz_pres, oz_time, ozplin
    integer, intent(in) :: n_columns, n_levels
    real(kind_phys), intent(in) :: lats(:), pres(:)

    integer, intent(inout) :: n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times
    real(kind_phys), allocatable, intent(inout) :: ozone_lat(:), ozone_pres(:), ozone_time(:), ozone_forcing_in(:,:,:,:)

    !set ozone forcing array dimensions
    n_ozone_lats = n_columns
    n_ozone_layers = n_levels
    n_ozone_coefficients = 4
    n_ozone_times = 2

    !reallocate ozone forcing latitutes, pressure levels, and times using the array dimensions set here
    if(allocated(ozone_lat)) then
      deallocate(ozone_lat)
      allocate(ozone_lat(n_ozone_lats))
    endif
    ozone_lat = lats

    if(allocated(ozone_pres)) then
      deallocate(ozone_pres)
      allocate(ozone_pres(n_ozone_layers))
    endif
    ozone_pres = pres

    if(allocated(ozone_time)) then
      deallocate(ozone_time)
      allocate(ozone_time(n_ozone_times+1))
    endif
    ozone_time = (/12.0, 13.0, 14.0/)

    if(allocated(ozone_forcing_in))then
      deallocate(ozone_forcing_in)
      allocate(ozone_forcing_in(n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times))
    endif
    ozone_forcing_in = 0.0

    !allocate memory for the variables stored in ozne_def and set them
    allocate(oz_lat(n_ozone_lats), oz_pres(n_ozone_layers), oz_time(n_ozone_times+1))
    allocate(ozplin(n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times))
    latsozp = n_ozone_lats
    levozp = n_ozone_layers
    timeoz = n_ozone_times
    oz_coeff = n_ozone_coefficients
    oz_lat = ozone_lat
    oz_pres = ozone_pres
    oz_time = ozone_time
    ozplin = ozone_forcing_in

  end subroutine GFS_ozone_and_h20_setup


!--------------
! GFS initialze
!--------------
  subroutine GFS_initialize (Model, Statein, Stateout, Sfcprop,     &
                             Coupling, Grid, Tbd, Cldprop, Radtend, &
                             Diag, Sfccycle, Interstitial, Init_parm, n_columns)

    use module_microphysics, only: gsmconst
    use cldwat2m_micro,      only: ini_micro
    use aer_cloud,           only: aer_cloud_init
    use module_ras,          only: ras_init
#ifdef OPENMP
    use omp_lib
#endif


    !--- interface variables
    type(GFS_control_type),      intent(inout) :: Model
    type(GFS_statein_type),      intent(inout) :: Statein
    type(GFS_stateout_type),     intent(inout) :: Stateout
    type(GFS_sfcprop_type),      intent(inout) :: Sfcprop
    type(GFS_coupling_type),     intent(inout) :: Coupling
    type(GFS_grid_type),         intent(inout) :: Grid
    type(GFS_tbd_type),          intent(inout) :: Tbd
    type(GFS_cldprop_type),      intent(inout) :: Cldprop
    type(GFS_radtend_type),      intent(inout) :: Radtend
    type(GFS_diag_type),         intent(inout) :: Diag
    type(GFS_sfccycle_type),     intent(inout) :: Sfccycle
    type(GFS_interstitial_type), intent(inout) :: Interstitial
    type(GFS_init_type),         intent(in)    :: Init_parm
    integer,                     intent(in)    :: n_columns
!
!     !--- local variables
!     integer :: nb
!     integer :: nblks
!     integer :: blkszmax
!     integer :: nt
!     integer :: nthreads
     integer :: ntrac
     real(kind=kind_phys), allocatable :: si(:)
     real(kind=kind_phys), parameter   :: p_ref = 101325.0d0
!
!     nblks = size(Init_parm%blksz)
     ntrac = size(Init_parm%tracer_names)
!     allocate (blksz(nblks))
!     blksz(:) = Init_parm%blksz(:)
!
!     !--- set control properties (including namelist read)
    call Model%init (Init_parm%nlunit, Init_parm%fn_nml,           &
                     Init_parm%me, Init_parm%master,               &
                     Init_parm%logunit, Init_parm%isc,             &
                     Init_parm%jsc, Init_parm%nx, Init_parm%ny,    &
                     Init_parm%levs, Init_parm%cnx, Init_parm%cny, &
                     Init_parm%gnx, Init_parm%gny,                 &
                     Init_parm%dt_dycore, Init_parm%dt_phys,       &
                     Init_parm%bdat, Init_parm%cdat,               &
                     Init_parm%tracer_names, Init_parm%blksz)



!
!     call read_o3data  (Model%ntoz, Model%me, Model%master)
!     call read_h2odata (Model%h2o_phys, Model%me, Model%master)
!
!     blkszmax = 0
!     do nb = 1,nblks
      call Statein%create(1, Model)
      call Stateout%create(1, Model)
      call Sfcprop%create(1, Model)
      call Coupling%create(1, Model)
      call Grid%create(1, Model)
      call Tbd%create(1, 1, Model)
      call Cldprop%create(1, Model)
      call Radtend%create(1, Model)
      !--- internal representation of diagnostics
      call Diag%create(1, Model)
      !--- internal representation of sfccycle
      call Sfccycle%create(1, Model)
!       !--- maximum blocksize
!       blkszmax = max(blkszmax, Init_parm%blksz(nb))
!     enddo
!
! #ifdef CCPP
! #ifdef OPENMP
!     nthreads = omp_get_max_threads()
! #else
!     nthreads = 1
! #endif
!     do nt=1,nthreads


  call Interstitial%create(1, Model)
!     enddo
! #endif
!
!     !--- populate the grid components
    call GFS_grid_populate (Grid, Init_parm%xlon, Init_parm%xlat, Init_parm%area)

!
!     !--- read in and initialize ozone and water
!     if (Model%ntoz > 0) then
!       do nb = 1, nblks
!         call setindxoz (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_o3, &
!                         Grid(nb)%jindx2_o3, Grid(nb)%ddy_o3)
!       enddo
!     endif
!
!     if (Model%h2o_phys) then
!       do nb = 1, nblks
!         call setindxh2o (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_h, &
!                          Grid(nb)%jindx2_h, Grid(nb)%ddy_h)
!       enddo
!     endif
!
!     !--- Call gfuncphys (funcphys.f) to compute all physics function tables.
     call gfuncphys ()
!
     call gsmconst (Model%dtp, Model%me, .TRUE.)
!
!     !--- define sigma level for radiation initialization
!     !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
!     !--- ps is replaced with p0. The value of p0 uses that in http://www.emc.ncep.noaa.gov/officenotes/newernotes/on461.pdf
!     !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
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
!
!     !--- initialize Morrison-Gettleman microphysics
    if (Model%ncld == 2) then
      call ini_micro (Model%mg_dcs, Model%mg_qcvar, Model%mg_ts_auto_ice)
      call aer_cloud_init ()
    endif
!
!     !--- initialize ras
     if (Model%ras) call ras_init (Model%levs, Model%me)
!
!     !--- initialize soil vegetation
     call set_soilveg(Model%me, Model%isot, Model%ivegsrc, Model%nlunit)
!
!     !--- lsidea initialization
    if (Model%lsidea) then
      print *,' LSIDEA is active but needs to be reworked for FV3 - shutting down'
      stop
      !--- NEED TO get the logic from the old phys/gloopb.f initialization area
    endif
!
!     !--- sncovr may not exist in ICs from chgres.
!     !--- FV3GFS handles this as part of the IC ingest
!     !--- this not is placed here to alert users to the need to study
!     !--- the FV3GFS_io.F90 module
!
   end subroutine GFS_initialize
!
!
! #ifndef CCPP
! !-------------------------------------------------------------------------
! ! time_vary_step
! !-------------------------------------------------------------------------
! !    routine called prior to radiation and physics steps to handle:
! !      1) sets up various time/date variables
! !      2) sets up various triggers
! !      3) defines random seed indices for radiation (in a reproducible way)
! !      5) interpolates coefficients for prognostic ozone calculation
! !      6) performs surface data cycling via the GFS gcycle routine
! !-------------------------------------------------------------------------
!   subroutine GFS_time_vary_step (Model, Statein, Stateout, Sfcprop, Coupling, &
!                                  Grid, Tbd, Cldprop, Radtend, Diag, Sfccycle)
!
!     use GFS_phys_time_vary_1,  only: GFS_phys_time_vary_1_run
!     use GFS_phys_time_vary_2,  only: GFS_phys_time_vary_2_run
!     use GFS_rad_time_vary,     only: GFS_rad_time_vary_run
!     implicit none
!
!     !--- interface variables
!     type(GFS_control_type),   intent(inout) :: Model
!     type(GFS_statein_type),   intent(inout) :: Statein
!     type(GFS_stateout_type),  intent(inout) :: Stateout
!     type(GFS_sfcprop_type),   intent(inout) :: Sfcprop
!     type(GFS_coupling_type),  intent(inout) :: Coupling
!     type(GFS_grid_type),      intent(inout) :: Grid
!     type(GFS_tbd_type),       intent(inout) :: Tbd
!     type(GFS_cldprop_type),   intent(inout) :: Cldprop
!     type(GFS_radtend_type),   intent(inout) :: Radtend
!     type(GFS_diag_type),      intent(inout) :: Diag
!     type(GFS_sfccycle_type),  intent(inout) :: Sfccycle
!
!     call GFS_phys_time_vary_1_run (Model, Tbd)
!
!     call GFS_rad_time_vary_run (Model, Statein, Tbd)
!
!     call GFS_phys_time_vary_2_run (Grid, Model, Tbd, Sfcprop, Cldprop, Diag, Sfccycle)
!
!   end subroutine GFS_time_vary_step
!
!
! !-------------------------------------------------------------------------
! ! GFS stochastic_driver
! !-------------------------------------------------------------------------
! !    routine called prior to radiation and physics steps to handle:
! !      1) sets up various time/date variables
! !      2) sets up various triggers
! !      3) defines random seed indices for radiation (in a reproducible way)
! !      5) interpolates coefficients for prognostic ozone calculation
! !      6) performs surface data cycling via the GFS gcycle routine
! !-------------------------------------------------------------------------
!   subroutine GFS_stochastic_driver (Model, Statein, Stateout, Sfcprop, Coupling, &
!                                     Grid, Tbd, Cldprop, Radtend, Diag)
!
!     use GFS_stochastics, only: GFS_stochastics_run
!
!     implicit none
!
!     !--- interface variables
!     type(GFS_control_type),   intent(in   ) :: Model
!     type(GFS_statein_type),   intent(in   ) :: Statein
!     type(GFS_stateout_type),  intent(in   ) :: Stateout
!     type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
!     type(GFS_coupling_type),  intent(inout) :: Coupling
!     type(GFS_grid_type),      intent(in   ) :: Grid
!     type(GFS_tbd_type),       intent(in   ) :: Tbd
!     type(GFS_cldprop_type),   intent(in   ) :: Cldprop
!     type(GFS_radtend_type),   intent(in   ) :: Radtend
!     type(GFS_diag_type),      intent(inout) :: Diag
!
!     call GFS_stochastics_run(Model, Statein, Stateout, Sfcprop, Coupling, &
!                              Grid, Tbd, Cldprop, Radtend, Diag)
!
!   end subroutine GFS_stochastic_driver
! #endif
!
!
!------------------
! GFS_grid_populate
!------------------
  subroutine GFS_grid_populate (Grid, xlon, xlat, area)
    use physcons,                 only: pi => con_pi

    implicit none

    type(GFS_grid_type)              :: Grid
    real(kind=kind_phys), intent(in) :: xlon(:,:)
    real(kind=kind_phys), intent(in) :: xlat(:,:)
    real(kind=kind_phys), intent(in) :: area(:,:)

    !--- local variables
    !integer :: nb, ix, blksz, i, j
    integer :: n_columns, i

    n_columns = size(Grid%xlon)

    do i=1, n_columns
      Grid%xlon = xlon(i,1)
      Grid%xlat = xlat(i,1)
      Grid%xlat_d(i) = xlat(i,1) * 180.0_kind_phys/pi
      Grid%sinlat(i) = sin(Grid%xlat(i))
      Grid%coslat(i) = sqrt(1.0_kind_phys - Grid%sinlat(i)*Grid%sinlat(i))
      Grid%area(i)   = area(i,1)
      Grid%dx(i)     = sqrt(area(i,1))
    end do

    ! blksz = size(Grid(1)%xlon)

    ! nb = 1
    ! ix = 0
    ! do j = 1,size(xlon,2)
    !   do i = 1,size(xlon,1)
    !     ix=ix+1
    !     if (ix .gt. blksz) then
    !       nb = nb + 1
    !       ix = 1
    !     endif
    !     Grid(nb)%xlon(ix)   = xlon(i,j)
    !     Grid(nb)%xlat(ix)   = xlat(i,j)
    !     Grid(nb)%xlat_d(ix) = xlat(i,j) * 180.0_kind_phys/pi
    !     Grid(nb)%sinlat(ix) = sin(Grid(nb)%xlat(ix))
    !     Grid(nb)%coslat(ix) = sqrt(1.0_kind_phys - Grid(nb)%sinlat(ix)*Grid(nb)%sinlat(ix))
    !     Grid(nb)%area(ix)   = area(i,j)
    !     Grid(nb)%dx(ix)     = sqrt(area(i,j))
    !   enddo
    ! enddo



  end subroutine GFS_grid_populate

end module GFS_driver_gmtb_scm

!***********************************************************************
!> TODO: replace with appropriate licence for CCPP
!*   GNU Lesser General Public License
!*   <http://www.gnu.org/licenses/>.
!***********************************************************************

!> @brief Land IAU (Incremental Analysis Update) module, 
!> for the NoahMP soil/snow temperature (can be extended to include soil moisture) 

!! \section land_iau_mod
!> - reads settings from namelist file (which indicates if IAU increments are available or not)
!> - reads in DA increments from GSI/JEDI DA at the start of (the DA) cycle 
!> - interpolates increments to FV3 grid points (if increments are in Gaussian grid)
!> - interpolates temporally (with filter, weights if required by configuration)
!> - updates states with the interpolated increments

!> March, 2024: Tseganeh Z. Gichamo, (EMC) based on the FV3 IAU mod 
!> by Xi.Chen <xi.chen@noaa.gov> and Philip Pegion, PSL <philip.pegion@noaa.gov>
!-------------------------------------------------------------------------------

module land_iau_mod
  
  use machine,                  only: kind_phys, kind_dyn
  use physcons,                 only: pi => con_pi
  use netcdf

  implicit none

  private

  real(kind=kind_phys),allocatable::s2c(:,:,:)

  integer,allocatable,dimension(:,:) :: id1,id2,jdc

  real(kind=kind_phys) :: deg2rad,dt,rdt
  integer :: im,jm,km,nfiles,ncid
  integer:: jbeg, jend

  integer :: n_soill, n_snowl              !soil and snow layers
  logical :: do_land_iau 

  integer :: is,  ie,  js,  je
  integer :: npz     

  real(kind=kind_phys), allocatable:: wk3_stc(:, :, :, :), wk3_slc(:, :, :, :)

  type land_iau_internal_data_type
      real(kind=kind_phys),allocatable :: stc_inc(:,:,:)
      real(kind=kind_phys),allocatable :: slc_inc(:,:,:) 
  end type land_iau_internal_data_type

  type land_iau_external_data_type
      real(kind=kind_phys),allocatable :: stc_inc(:,:,:)   
      real(kind=kind_phys),allocatable :: slc_inc(:,:,:)   
      logical          :: in_interval = .false.
  end type land_iau_external_data_type

  type land_iau_state_type
      type(land_iau_internal_data_type):: inc1
      type(land_iau_internal_data_type):: inc2
      real(kind=kind_phys)        :: hr1
      real(kind=kind_phys)        :: hr2
      real(kind=kind_phys)        :: wt
      real(kind=kind_phys)        :: wt_normfact
  end type land_iau_state_type

  type land_iau_control_type      
      integer :: isc
      integer :: jsc
      integer :: nx
      integer :: ny
      integer :: tile_num
      integer :: nblks
      integer, allocatable :: blksz(:)    ! this could vary for the last block
      integer, allocatable :: blk_strt_indx(:)

      integer :: lsoil  !< number of soil layers
      ! this is the max dim (TBC: check it is consitent for noahmpdrv)
      integer              :: lsnow_lsm       !< maximum number of snow layers internal to land surface model
      logical              :: do_land_iau
      real(kind=kind_phys) :: iau_delthrs     ! iau time interval (to scale increments) in hours
      character(len=240)   :: iau_inc_files(7)! list of increment files
      real(kind=kind_phys) :: iaufhrs(7)      ! forecast hours associated with increment files
      logical              :: iau_filter_increments   
      integer              :: lsoil_incr    ! soil layers (from top) updated by DA   
      !, iau_drymassfixer
      integer              :: me              !< MPI rank designator
      integer              :: mpi_root          !< MPI rank of master atmosphere processor
      character(len=64)    :: fn_nml          !< namelist filename for surface data cycling
      real(kind=kind_phys) :: dtp             !< physics timestep in seconds
      real(kind=kind_phys) :: fhour           !< current forecast hour
      character(len=:), pointer, dimension(:) :: input_nml_file => null() !<character string containing full namelist
                                                                          !< for use with internal file reads
      integer              :: input_nml_file_length    !<length (number of lines) in namelist for internal reads
      logical              :: gaussian_inc_file   ! if true increment in gaussian grids
                              !  else fv3

  end type land_iau_control_type

  type(land_iau_state_type) :: Land_IAU_state
  public land_iau_control_type, land_iau_external_data_type, land_iau_mod_set_control, &
         land_iau_mod_init, land_iau_mod_getiauforcing, land_iau_mod_finalize

contains

subroutine land_iau_mod_set_control(Land_IAU_Control,fn_nml,input_nml_file_i, me, mpi_root, &
                                   isc, jsc, nx, ny, tile_num, nblks, blksz, &
                                   lsoil, lsnow_lsm, dtp, fhour, errmsg, errflg)          !nlunit

   type (land_iau_control_type), intent(inout) :: Land_IAU_Control
   character(*), intent(in)                   :: fn_nml          !< namelist filename for surface data cycling
   character(len=:), intent(in), dimension(:), pointer :: input_nml_file_i
   integer, intent(in)                        :: me, mpi_root          !< MPI rank of master atmosphere processor   
   integer, intent(in)                        :: isc, jsc, nx, ny, tile_num, nblks, lsoil, lsnow_lsm
   integer, dimension(:),          intent(in) :: blksz   !(one:) !GFS_Control%blksz
   real(kind=kind_phys), intent(in)           :: dtp             !< physics timestep in seconds
   real(kind=kind_phys), intent(in)           :: fhour           !< current forecast hour
   character(len=*),              intent(out) :: errmsg
   integer,                       intent(out) :: errflg
   
   integer                                    :: nb, ix
   integer                                    :: nlunit = 360          ! unit for namelist  !, intent(in)
   integer                                    :: ios
   logical                                    :: exists
   character(len=512)                         :: ioerrmsg
   !character(len=32)                          :: fn_nml = "input.nml"
   character(len=:), pointer, dimension(:)    :: input_nml_file => null()
   integer                    :: input_nml_file_length    !< length(number of lines) in namelist for internal reads
   

   !> these are not available through the CCPP interface so need to read them from namelist file
   !> vars to read from namelist
   logical               :: do_land_iau               = .false.
   real(kind=kind_phys)  :: land_iau_delthrs              = 0           !< iau time interval (to scale increments)
   character(len=240)    :: land_iau_inc_files(7)         = ''          !< list of increment files
   real(kind=kind_phys)  :: land_iaufhrs(7)               = -1          !< forecast hours associated with increment files
   logical               :: land_iau_filter_increments    = .false.     !< filter IAU increments
   logical               :: gaussian_inc_file             = .false.
   integer               :: lsoil_incr = 4
  
   NAMELIST /lnd_iau_nml/ do_land_iau, land_iau_delthrs, land_iau_inc_files, land_iaufhrs, &
                        land_iau_filter_increments, &  !, lnd_iau_drymassfixer   
                        lsoil_incr, gaussian_inc_file                                      
   
   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0

!3.11.24: copied from GFS_typedefs.F90 
#ifdef INTERNAL_FILE_NML
    ! allocate required to work around GNU compiler bug 100886
    ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100886
    allocate(input_nml_file, mold=input_nml_file_i)
    input_nml_file => input_nml_file_i
    read(input_nml_file, nml=lnd_iau_nml)
    ! Set length (number of lines) in namelist for internal reads
    input_nml_file_length = size(input_nml_file)
#else
   ! if (file_exist(fn_nml)) then 
   inquire (file=trim(fn_nml), exist=exists)    ! TBCL: this maybe be replaced by nlunit passed from ccpp
   if (.not. exists) then
      ! call mpp_error(FATAL, 'lnd_iau_mod_set_control: namelist file ',trim(fn_nml),' does not exist')
      write(6,*) 'lnd_iau_mod_set_control: namelist file ',trim(fn_nml),' does not exist'      
      errmsg = 'lnd_iau_mod_set_control: namelist file '//trim(fn_nml)//' does not exist'
      errflg = 1
      return
   else
      Land_IAU_Control%fn_nml = trim(fn_nml)   ! maynot need this
      open (unit=nlunit, file=trim(fn_nml), action='READ', status='OLD', iostat=ios, iomsg=ioerrmsg)
      rewind(nlunit)
      read (nlunit, nml=lnd_iau_nml)
      close (nlunit)
      if (ios /= 0) then
         ! call mpp_error(FATAL, 'lnd_iau_mod_set_control: error reading namelist file ',trim(fn_nml))
         ! write(6,*) 'lnd_iau_mod_set_control: error reading namelist file ',trim(fn_nml)
         write(6,*) trim(ioerrmsg)         
         errmsg = 'lnd_iau_mod_set_control: error reading namelist file '//trim(fn_nml)  &
                  // 'the error message from file handler:' //trim(ioerrmsg) 
         errflg = 1
         return
       end if
   endif 
#endif

   if (me == mpi_root) then
      write(6,*) "lnd_iau_nml"
      write(6, lnd_iau_nml)
   endif
   
   Land_IAU_Control%do_land_iau = do_land_iau
   Land_IAU_Control%iau_delthrs = land_iau_delthrs
   Land_IAU_Control%iau_inc_files = land_iau_inc_files
   Land_IAU_Control%iaufhrs = land_iaufhrs   
   Land_IAU_Control%iau_filter_increments = land_iau_filter_increments
   ! Land_IAU_Control%iau_drymassfixer = lnd_iau_drymassfixer
   Land_IAU_Control%lsoil_incr = lsoil_incr

   Land_IAU_Control%me = me
   Land_IAU_Control%mpi_root = mpi_root
   Land_IAU_Control%isc = isc
   Land_IAU_Control%jsc = jsc
   Land_IAU_Control%nx = nx
   Land_IAU_Control%ny = ny
   Land_IAU_Control%tile_num = tile_num
   Land_IAU_Control%nblks = nblks
   Land_IAU_Control%lsoil = lsoil
   Land_IAU_Control%lsnow_lsm = lsnow_lsm
   Land_IAU_Control%dtp = dtp
   Land_IAU_Control%fhour = fhour

   Land_IAU_Control%input_nml_file = input_nml_file
   Land_IAU_Control%input_nml_file_length = input_nml_file_length
   Land_IAU_Control%gaussian_inc_file = gaussian_inc_file

   allocate(Land_IAU_Control%blksz(nblks))
   allocate(Land_IAU_Control%blk_strt_indx(nblks))

   ! Land_IAU_Control%blk_strt_indx: start index of each block, for flattened (ncol=nx*ny) arrays 
   ! required in noahmpdriv_run to get subsection of the stc array for each
   ! proces/thread
   ix = 1
   do nb=1, nblks
      Land_IAU_Control%blksz(nb) = blksz(nb)
      Land_IAU_Control%blk_strt_indx(nb) = ix
      ix = ix + blksz(nb)
   enddo

end subroutine land_iau_mod_set_control

subroutine land_iau_mod_init (Land_IAU_Control, Land_IAU_Data, xlon, xlat, errmsg, errflg)     !nlunit, ncols, IPD_Data,,Init_parm)
   ! integer,                              intent(in) :: me, mpi_root
   type (land_iau_control_type),          intent(in) :: Land_IAU_Control
   type (land_iau_external_data_type), intent(inout) :: Land_IAU_Data  
   real(kind=kind_phys), dimension(:), intent(in)  :: xlon    ! longitude  !GFS_Data(cdata%blk_no)%Grid%xlon
   real(kind=kind_phys), dimension(:), intent(in)  :: xlat    ! latitude
   character(len=*),              intent(out) :: errmsg
   integer,                       intent(out) :: errflg

   ! local
   character(len=128) :: fname
   real(kind=kind_phys), allocatable:: lat(:), lon(:),agrid(:,:,:)
   real(kind=kind_phys) sx,wx,wt,normfact,dtp
   integer:: ib, i, j, k, nstep, kstep
   integer:: i1, i2, j1
   logical:: found
   integer nfilesall
   integer, allocatable :: idt(:)
   real (kind=kind_phys), allocatable :: Init_parm_xlon (:, :)   
   real (kind=kind_phys), allocatable :: Init_parm_xlat (:, :)   
   integer :: nlon, nlat
   ! integer :: nb, ix, nblks, blksz  
   logical                            :: exists
   ! necdf vars
   integer           :: ncid, dimid, varid, status, IDIM

   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0

   do_land_iau = Land_IAU_Control%do_land_iau
   n_soill = Land_IAU_Control%lsoil     !4  for sfc updates
!  n_snowl = Land_IAU_Control%lsnowl 
   npz = Land_IAU_Control%lsoil
   
   is  = Land_IAU_Control%isc
   ie  = is + Land_IAU_Control%nx-1
   js  = Land_IAU_Control%jsc
   je  = js + Land_IAU_Control%ny-1
   nlon = Land_IAU_Control%nx
   nlat = Land_IAU_Control%ny
   !nblks = Land_IAU_Control%nblks
   !blksz = Land_IAU_Control%blksz(1)

   allocate(Init_parm_xlon(nlon,nlat), Init_parm_xlat(nlon,nlat))
   ib = 1
   do j = 1, nlat  !ny
      ! do i = 1, nx      
         Init_parm_xlon (:,j) = xlon(ib:ib+nlon-1)  
         Init_parm_xlat (:,j) = xlat(ib:ib+nlon-1) 
         ib = ib+nlon
      ! enddo
   enddo
 
   allocate(s2c(is:ie,js:je,4))
   allocate(id1(is:ie,js:je))
   allocate(id2(is:ie,js:je))
   allocate(jdc(is:ie,js:je))
   allocate(agrid(is:ie,js:je,2))

! determine number of increment files to read, and the valid forecast hours
   nfilesall = size(Land_IAU_Control%iau_inc_files)
   nfiles = 0
   if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print*,'in land_iau_init incfile1 iaufhr1 ', &
                                 trim(Land_IAU_Control%iau_inc_files(1)),Land_IAU_Control%iaufhrs(1)
   do k=1,nfilesall
      if (trim(Land_IAU_Control%iau_inc_files(k)) .eq. '' .or. Land_IAU_Control%iaufhrs(k) .lt. 0) exit   
      if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then
         print *,k, " ", trim(adjustl(Land_IAU_Control%iau_inc_files(k)))
      endif
      nfiles = nfiles + 1
   enddo
   if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'nfiles = ',nfiles
   if (nfiles < 1) then
      return
   endif
   if (nfiles > 1) then
      allocate(idt(nfiles-1))
      idt = Land_IAU_Control%iaufhrs(2:nfiles)-Land_IAU_Control%iaufhrs(1:nfiles-1)
      do k=1,nfiles-1
         if (idt(k) .ne. Land_IAU_Control%iaufhrs(2)-Land_IAU_Control%iaufhrs(1)) then
           print *,'in land_iau_init: forecast intervals in iaufhrs must be constant'
         !   call mpp_error (FATAL,' forecast intervals in iaufhrs must be constant')
           errmsg = 'Fatal error in land_iau_init. forecast intervals in iaufhrs must be constant'
           errflg = 1
           return
         endif
      enddo
      deallocate(idt)
   endif
   if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'land_iau interval = ',Land_IAU_Control%iau_delthrs,' hours'
   dt = (Land_IAU_Control%iau_delthrs*3600.)
   rdt = 1.0/dt

   allocate(Land_IAU_Data%stc_inc(is:ie, js:je, km))
   allocate(Land_IAU_Data%slc_inc(is:ie, js:je, km))
! allocate arrays that will hold iau state
   allocate (Land_IAU_state%inc1%stc_inc(is:ie, js:je, km))
   allocate (Land_IAU_state%inc1%slc_inc(is:ie, js:je, km))
   Land_IAU_state%hr1=Land_IAU_Control%iaufhrs(1)
   Land_IAU_state%wt = 1.0 ! IAU increment filter weights (default 1.0)
   Land_IAU_state%wt_normfact = 1.0
   if (Land_IAU_Control%iau_filter_increments) then
      ! compute increment filter weights, sum to obtain normalization factor
      dtp=Land_IAU_Control%dtp
      nstep = 0.5*Land_IAU_Control%iau_delthrs*3600/dtp
      ! compute normalization factor for filter weights
      normfact = 0.
      do k=1,2*nstep+1
         kstep = k-1-nstep
         sx     = acos(-1.)*kstep/nstep
         wx     = acos(-1.)*kstep/(nstep+1)
         if (kstep .ne. 0) then
            wt = sin(wx)/wx*sin(sx)/sx
         else
            wt = 1.0
         endif
         normfact = normfact + wt
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'filter wts',k,kstep,wt
      enddo
      Land_IAU_state%wt_normfact = (2*nstep+1)/normfact
   endif

   ! Find bounding latitudes:
   jbeg = jm-1
   jend = 2
   do j=js,je
      do i=is,ie
            j1 = jdc(i,j)
         jbeg = min(jbeg, j1)
         jend = max(jend, j1+1)
      enddo
   enddo
   
   if (Land_IAU_Control%gaussian_inc_file) then 
      !set up interpolation weights to go from GSI's gaussian grid to cubed sphere
      deg2rad = pi/180.
      !  npz = Land_IAU_Control%levs
      fname = 'INPUT/'//trim(Land_IAU_Control%iau_inc_files(1))    
      inquire (file=trim(fname), exist=exists)    
      if (exists) then   !  if( file_exist(fname) ) then
         ! call open_ncfile( fname, ncid )        
         status = nf90_open(trim(fname), NF90_NOWRITE, ncid)  ! open the file
         call netcdf_err(status, ' opening file '//trim(fname), errflg, errmsg) 
         if (errflg .ne. 0) return     
         call get_nc_dimlen(ncid, "longitude", im, errflg, errmsg) 
         if (errflg .ne. 0) return     
         call get_nc_dimlen(ncid, "latitude", jm, errflg, errmsg) 
         if (errflg .ne. 0) return     
         km = n_soill
         if(Land_IAU_Control%me == Land_IAU_Control%mpi_root)  write(*,*) fname, ' DA increment dimensions:', im,jm,km
         allocate (lon(im) )
         allocate (lat(jm) )
         call get_var1d(ncid, im, "longitude", lon, errflg, errmsg) 
         if (errflg .ne. 0) return     
         call get_var1d(ncid, jm, "latitude", lat, errflg, errmsg) 
         if (errflg .ne. 0) return     
         status = nf90_close(ncid)
         CALL netcdf_err(status, 'closing file: '//trim(fname) , errflg, errmsg) 
         if (errflg .ne. 0) return     
         ! Convert to radians
         do i=1,im
         lon(i) = lon(i) * deg2rad
         enddo
         do j=1,jm
         lat(j) = lat(j) * deg2rad
         enddo
      else
         errmsg = 'FATAL Error in Land_IAU_initialize: Expected file '// trim(fname)//' for DA increment does not exist'
         errflg = 1
         return
      endif

      ! Initialize lat-lon to Cubed bi-linear interpolation coeff:
      do j = 1,size(Init_parm_xlon,2)
         do i = 1,size(Init_parm_xlon,1)
      !         print*,i,j,is-1+j,js-1+j
            agrid(is-1+i,js-1+j,1)=Init_parm_xlon(i,j)
            agrid(is-1+i,js-1+j,2)=Init_parm_xlat(i,j)
         enddo
      enddo
      call remap_coef( is, ie, js, je, is, ie, js, je, &
         im, jm, lon, lat, id1, id2, jdc, s2c, &
         agrid)

      if (allocated(lon)) deallocate (lon)
      if (allocated(lat)) deallocate (lat)
      if (allocated(agrid)) deallocate (agrid)
      if (allocated(Init_parm_xlon)) deallocate(Init_parm_xlon)
      if (allocated(Init_parm_xlat)) deallocate(Init_parm_xlat)
      
      ! Read all increment files at iau init time (at beginning of cycle) and interpolate to target grid
      allocate (wk3_stc(nfiles, 1:im,jbeg:jend, 1:km))
      allocate (wk3_slc(nfiles, 1:im,jbeg:jend, 1:km))   
      do k=1, nfiles
         call read_iau_forcing_all_timesteps(Land_IAU_Control, &
         'INPUT/'//trim(Land_IAU_Control%iau_inc_files(k)), errmsg, errflg, &
            wk3_stc(k, :, :, :), wk3_slc(k, :, :, :)) 
      enddo   
      call interp_inc_at_timestep(Land_IAU_Control, km, wk3_stc(1, :, :, :), Land_IAU_state%inc1%stc_inc, errmsg, errflg)
      call interp_inc_at_timestep(Land_IAU_Control, km, wk3_slc(1, :, :, :), Land_IAU_state%inc1%slc_inc, errmsg, errflg)
      if (nfiles.EQ.1) then  ! only need to get incrments once since constant forcing over window
         call setiauforcing(Land_IAU_Control, Land_IAU_Data, Land_IAU_state%wt)
      endif
      if (nfiles.GT.1) then  !have multiple files, but only need 2 at a time and interpoalte for timesteps between them
         allocate (Land_IAU_state%inc2%stc_inc(is:ie, js:je, km))
         allocate (Land_IAU_state%inc2%slc_inc(is:ie, js:je, km))      
         Land_IAU_state%hr2=Land_IAU_Control%iaufhrs(2)      
         call interp_inc_at_timestep(Land_IAU_Control, km, wk3_stc(2, :, :, :), Land_IAU_state%inc2%stc_inc, errmsg, errflg)
         call interp_inc_at_timestep(Land_IAU_Control, km, wk3_slc(2, :, :, :), Land_IAU_state%inc2%slc_inc, errmsg, errflg)
      endif
   else   ! increment files in fv3 tiles 
      ! Read all increment files at iau init time (at beginning of cycle) and interpolate to target grid
      ! increments already in the fv3 modele grid--no need for interpolation 
      ! Read all increment files at iau init time (at beginning of cycle) and interpolate to target grid
      allocate (wk3_stc(nfiles, 1:im,jbeg:jend, 1:km))
      allocate (wk3_slc(nfiles, 1:im,jbeg:jend, 1:km))   
      do k=1, nfiles
         call read_iau_forcing_fv3(Land_IAU_Control, &
         'INPUT/'//trim(Land_IAU_Control%iau_inc_files(k)), errmsg, errflg, &
             
            Land_IAU_state%inc1%stc_inc, Land_IAU_state%inc1%slc_inc, errmsg, errflg)
      enddo
      call interp_inc_at_timestep(Land_IAU_Control, km, wk3_stc(1, :, :, :), Land_IAU_state%inc1%stc_inc, errmsg, errflg)
      call interp_inc_at_timestep(Land_IAU_Control, km, wk3_slc(1, :, :, :), Land_IAU_state%inc1%slc_inc, errmsg, errflg)
      if (nfiles.EQ.1) then  ! only need to get incrments once since constant forcing over window
         call setiauforcing(Land_IAU_Control, Land_IAU_Data, Land_IAU_state%wt)
      endif
      if (nfiles.GT.1) then  !have multiple files, but only need 2 at a time and interpoalte for timesteps between them
         allocate (Land_IAU_state%inc2%stc_inc(is:ie, js:je, km))
         allocate (Land_IAU_state%inc2%slc_inc(is:ie, js:je, km))      
         Land_IAU_state%hr2=Land_IAU_Control%iaufhrs(2)      
         call interp_inc_at_timestep(Land_IAU_Control, km, wk3_stc(2, :, :, :), Land_IAU_state%inc2%stc_inc, errmsg, errflg)
         call interp_inc_at_timestep(Land_IAU_Control, km, wk3_slc(2, :, :, :), Land_IAU_state%inc2%slc_inc, errmsg, errflg)
      endif

   endif


!   print*,'end of IAU init',dt,rdt

end subroutine land_iau_mod_init

subroutine land_iau_mod_finalize(Land_IAU_Control, Land_IAU_Data, errmsg, errflg)

   implicit none

   type (land_iau_control_type),          intent(in) :: Land_IAU_Control
   type(land_iau_external_data_type),  intent(inout) :: Land_IAU_Data
   character(len=*),                    intent(out) :: errmsg
   integer,                             intent(out) :: errflg

   if (allocated (wk3_stc)) deallocate (wk3_stc)
   if (allocated (wk3_slc)) deallocate (wk3_slc)

   if (allocated(Land_IAU_Data%stc_inc)) deallocate (Land_IAU_Data%stc_inc)
   if (allocated(Land_IAU_Data%slc_inc)) deallocate (Land_IAU_Data%slc_inc)

   if (allocated(Land_IAU_state%inc1%stc_inc)) deallocate(Land_IAU_state%inc1%stc_inc)
   if (allocated(Land_IAU_state%inc1%slc_inc)) deallocate(Land_IAU_state%inc1%slc_inc)

   if (allocated(Land_IAU_state%inc2%stc_inc)) deallocate(Land_IAU_state%inc2%stc_inc)
   if (allocated(Land_IAU_state%inc2%slc_inc)) deallocate(Land_IAU_state%inc2%slc_inc)

end subroutine land_iau_mod_finalize

 subroutine land_iau_mod_getiauforcing(Land_IAU_Control, Land_IAU_Data, errmsg, errflg)

   implicit none
   type (land_iau_control_type),          intent(in) :: Land_IAU_Control
   type(land_iau_external_data_type),  intent(inout) :: Land_IAU_Data
   character(len=*),              intent(out) :: errmsg
   integer,                       intent(out) :: errflg
   real(kind=kind_phys) t1,t2,sx,wx,wt,dtp
   integer n,i,j,k,kstep,nstep,itnext

   Land_IAU_Data%in_interval=.false.
   if (nfiles.LE.0) then
       return
   endif

   if (nfiles .eq. 1) then 
       t1 = Land_IAU_Control%iaufhrs(1)-0.5*Land_IAU_Control%iau_delthrs
       t2 = Land_IAU_Control%iaufhrs(1)+0.5*Land_IAU_Control%iau_delthrs
   else
       t1 = Land_IAU_Control%iaufhrs(1)
       t2 = Land_IAU_Control%iaufhrs(nfiles)
   endif
   if (Land_IAU_Control%iau_filter_increments) then
      ! compute increment filter weight
      ! t1 is beginning of window, t2 end of window
      ! Land_IAU_Control%fhour current time
      ! in window kstep=-nstep,nstep (2*nstep+1 total)
      ! time step Land_IAU_Control%dtp
      dtp=Land_IAU_Control%dtp
      nstep = 0.5*Land_IAU_Control%iau_delthrs*3600/dtp
      ! compute normalized filter weight
      kstep = ((Land_IAU_Control%fhour-t1) - 0.5*Land_IAU_Control%iau_delthrs)*3600./dtp
      if (Land_IAU_Control%fhour >= t1 .and. Land_IAU_Control%fhour < t2) then
         sx     = acos(-1.)*kstep/nstep
         wx     = acos(-1.)*kstep/(nstep+1)
         if (kstep .ne. 0) then
            wt = (sin(wx)/wx*sin(sx)/sx)
         else
            wt = 1.
         endif
         Land_IAU_state%wt = Land_IAU_state%wt_normfact*wt
         !if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'kstep,t1,t,t2,filter wt=',kstep,t1,Land_IAU_Control%fhour,t2,Land_IAU_state%wt/Land_IAU_state%wt_normfact
      else
         Land_IAU_state%wt = 0.
      endif
   endif

   if (nfiles.EQ.1) then
      !  check to see if we are in the IAU window,  
      ! no need to update the states since they are fixed over the window
      if ( Land_IAU_Control%fhour < t1 .or. Land_IAU_Control%fhour >= t2 ) then
!         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'no iau forcing',t1,Land_IAU_Control%fhour,t2
         Land_IAU_Data%in_interval=.false.
      else
         if (Land_IAU_Control%iau_filter_increments) call setiauforcing(Land_IAU_Control,Land_IAU_Data,Land_IAU_state%wt)
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'apply lnd iau forcing t1,t,t2,filter wt= ',t1,Land_IAU_Control%fhour,t2,Land_IAU_state%wt/Land_IAU_state%wt_normfact
         Land_IAU_Data%in_interval=.true.
      endif
      return
   endif

   if (nfiles > 1) then
      itnext=2
      if (Land_IAU_Control%fhour < t1 .or. Land_IAU_Control%fhour >= t2) then
!         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'no iau forcing',Land_IAU_Control%iaufhrs(1),Land_IAU_Control%fhour,Land_IAU_Control%iaufhrs(nfiles)
         Land_IAU_Data%in_interval=.false.
      else
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'apply lnd iau forcing t1,t,t2,filter wt= ',t1,Land_IAU_Control%fhour,t2,Land_IAU_state%wt/Land_IAU_state%wt_normfact
         Land_IAU_Data%in_interval=.true.
         do k=nfiles, 1, -1
            if (Land_IAU_Control%iaufhrs(k) > Land_IAU_Control%fhour) then
               itnext=k
            endif
         enddo
!         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'itnext=',itnext
         if (Land_IAU_Control%fhour >= Land_IAU_state%hr2) then ! need to read in next increment file
            Land_IAU_state%hr1=Land_IAU_state%hr2
            Land_IAU_state%hr2=Land_IAU_Control%iaufhrs(itnext)
            Land_IAU_state%inc1=Land_IAU_state%inc2
     
            ! if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'reading next lnd iau increment file',trim(Land_IAU_Control%iau_inc_files(itnext))
            if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'interpolating next lnd iau increment ', itnext  !trim(Land_IAU_Control%iau_inc_files(itnext))
            call interp_inc_at_timestep(Land_IAU_Control, km, wk3_stc(itnext, :, :, :), Land_IAU_state%inc2%stc_inc, errmsg, errflg)
            call interp_inc_at_timestep(Land_IAU_Control, km, wk3_slc(itnext, :, :, :), Land_IAU_state%inc2%slc_inc, errmsg, errflg)
            call interp_inc_at_timestep(Land_IAU_Control, 1, wk3_t2m(itnext, :, :, :), Land_IAU_state%inc2%tmp2m_inc, errmsg, errflg)
            call interp_inc_at_timestep(Land_IAU_Control, 1, wk3_q2m(itnext, :, :, :), Land_IAU_state%inc2%spfh2m_inc, errmsg, errflg)
         endif
         call updateiauforcing(Land_IAU_Control,Land_IAU_Data,Land_IAU_state%wt)
      endif
   endif

 end subroutine land_iau_mod_getiauforcing

subroutine updateiauforcing(Land_IAU_Control, Land_IAU_Data, wt)

   implicit none
   type (land_iau_control_type),        intent(in) :: Land_IAU_Control
   type(land_iau_external_data_type),  intent(inout) :: Land_IAU_Data
   real(kind=kind_phys) delt, wt
   integer i,j,k,l

!   if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'in updateiauforcing',nfiles,Land_IAU_Control%iaufhrs(1:nfiles)
   delt = (Land_IAU_state%hr2-(Land_IAU_Control%fhour))/(Land_IAU_state%hr2-Land_IAU_state%hr1)
   do j = js,je
      do i = is,ie
         do k = 1,npz  ! do k = 1,n_soill    !         
            Land_IAU_Data%stc_inc(i,j,k)  =(delt*Land_IAU_state%inc1%stc_inc(i,j,k)  + (1.-delt)* Land_IAU_state%inc2%stc_inc(i,j,k))*rdt*wt
            Land_IAU_Data%slc_inc(i,j,k)  =(delt*Land_IAU_state%inc1%slc_inc(i,j,k)  + (1.-delt)* Land_IAU_state%inc2%slc_inc(i,j,k))*rdt*wt
         end do
         Land_IAU_Data%tmp2m_inc(i,j,1)  =(delt*Land_IAU_state%inc1%tmp2m_inc(i,j,1)  + (1.-delt)* Land_IAU_state%inc2%tmp2m_inc(i,j,1))*rdt*wt
         Land_IAU_Data%spfh2m_inc(i,j,1)  =(delt*Land_IAU_state%inc1%spfh2m_inc(i,j,1)  + (1.-delt)* Land_IAU_state%inc2%spfh2m_inc(i,j,1))*rdt*wt
       enddo
   enddo
 end subroutine updateiauforcing

 subroutine setiauforcing(Land_IAU_Control, Land_IAU_Data, wt)

   implicit none
   type (land_iau_control_type),        intent(in)   :: Land_IAU_Control
   type(land_iau_external_data_type),  intent(inout) :: Land_IAU_Data
   real(kind=kind_phys) delt, dt,wt
   integer i,j,k,l
   !  this is only called if using 1 increment file
   if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,'in land_iau setiauforcing rdt = ',rdt
   do j = js,je
      do i = is,ie
         do k = 1,npz   !  do k = 1,n_soill    !         
            Land_IAU_Data%stc_inc(i,j,k) = wt*Land_IAU_state%inc1%stc_inc(i,j,k)*rdt
            Land_IAU_Data%slc_inc(i,j,k) = wt*Land_IAU_state%inc1%slc_inc(i,j,k)*rdt
         end do
         Land_IAU_Data%tmp2m_inc(i,j,1) = wt*Land_IAU_state%inc1%tmp2m_inc(i,j,1)*rdt
         Land_IAU_Data%spfh2m_inc(i,j,1) = wt*Land_IAU_state%inc1%spfh2m_inc(i,j,1)*rdt
      enddo
   enddo

 end subroutine setiauforcing

subroutine read_iau_forcing_all_timesteps(Land_IAU_Control, fname, errmsg, errflg, &
                                          wk3_out_stc, wk3_out_slc)   !, fname_sfc) is, ie, js, je, ks,ke, 
   type (land_iau_control_type),   intent(in) :: Land_IAU_Control
   character(len=*),              intent(in) :: fname
   character(len=*),             intent(inout) :: errmsg
   integer,                      intent(inout) :: errflg
   real(kind=kind_phys),                 intent(out) :: wk3_out_stc(1:im, jbeg:jend, 1:km)
   real(kind=kind_phys),                 intent(out) :: wk3_out_slc(1:im, jbeg:jend, 1:km)
   
   integer  :: i, j, k, l, npz
   integer  :: i1, i2, j1
   logical  :: exists
   integer  :: ncid, status, varid
   integer  :: ierr

   character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
   character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']

   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0
   
   inquire (file=trim(fname), exist=exists)    
   if (exists) then
      status = nf90_open(trim(fname), NF90_NOWRITE, ncid)  ! open the file
      call netcdf_err(status, ' opening file '//trim(fname), errflg, errmsg) 
      if (errflg .ne. 0) return
   else
      errmsg = 'FATAL Error in land iau read_iau_forcing_all_timesteps: Expected file '//trim(fname)//' for DA increment does not exist'
      errflg = 1
      return
   endif

   do i = 1, size(stc_vars)
      print *, trim(stc_vars(i))
      ! call check_var_exists(ncid, trim(stc_vars(i)), ierr)
      status = nf90_inq_varid(ncid, trim(stc_vars(i)), varid)
      if (status == nf90_noerr) then   !if (ierr == 0) then
         call get_var3d_values(ncid, varid, 1,im, jbeg,jend, 1,1, wk3_out_stc(:, :, i), status)
         call netcdf_err(status, 'reading var: '//trim(stc_vars(i)), errflg, errmsg_out)
         if (errflg .ne. 0) return 
      else
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *, &
         'warning: no increment for ',trim(stc_vars(i)),' found, assuming zero'
         wk3_out_stc(:, :, i) = 0.
      endif
   enddo
   do i = 1, size(slc_vars)
      print *, trim(slc_vars(i))
      status = nf90_inq_varid(ncid, trim(slc_vars(i)), varid)
      if (status == nf90_noerr) then   !if (ierr == 0) then
         ! call get_var3_r4( ncid, trim(slc_vars(i)), 1,im, jbeg,jend, 1,1, wk3_out_slc(:, :, i) )
         call get_var3d_values(ncid, varid, 1,im, jbeg,jend, 1,1, wk3_out_slc(:, :, i), status)
         call netcdf_err(status, 'reading var: '//trim(slc_vars(i)), errflg, errmsg_out)
         if (errflg .ne. 0) return          
      else
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,&
         'warning: no increment for ',trim(slc_vars(i)),' found, assuming zero'
         wk3_out_slc(:, :, i) = 0.
      endif
   enddo

   call close_ncfile(ncid)
   
end subroutine read_iau_forcing_all_timesteps

subroutine read_iau_forcing_fv3(Land_IAU_Control, stc_inc_out, slc_inc_out, errmsg, errflg)

   type (land_iau_control_type), intent(in) :: Land_IAU_Control
   ! character(len=*),             intent(in) :: fname
   character(len=*),          intent(inout) :: errmsg
   integer,                   intent(inout) :: errflg
   real(kind=kind_phys),        intent(out) :: stc_inc_out(1:im, jbeg:jend, 1:km)
   real(kind=kind_phys),        intent(out) :: slc_inc_out(1:im, jbeg:jend, 1:km)

   integer  :: i, j, k, l, npz
   integer  :: i1, i2, j1
   logical  :: exists
   integer  :: ncid, status, varid
   integer  :: ierr
   character(len=500)  :: fname
   character(len=2)    :: tile_str

   character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
   character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']

   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0

   write(tile_str, '(I0)') Land_IAU_Control%tile_num

   fname = 'INPUT/'//trim(Land_IAU_Control%iau_inc_files(1))//".tile"//tile_str//".nc"
   
   inquire (file=trim(fname), exist=exists)    
   if (exists) then
      status = nf90_open(trim(fname), NF90_NOWRITE, ncid)  ! open the file
      call netcdf_err(status, ' opening file '//trim(fname), errflg, errmsg) 
      if (errflg .ne. 0) return
   else
      errmsg = 'FATAL Error in land iau read_iau_forcing_fv3: Expected file '//trim(fname)//' for DA increment does not exist'
      errflg = 1
      return
   endif

   do i = 1, size(stc_vars)
      print *, trim(stc_vars(i))
      ! call check_var_exists(ncid, trim(stc_vars(i)), ierr)
      status = nf90_inq_varid(ncid, trim(stc_vars(i)), varid)
      if (status == nf90_noerr) then   !if (ierr == 0) then
         ! var stored as soilt1_inc(Time, yaxis_1, xaxis_1)
         call get_var3d_values(ncid, varid, 1,im, jbeg,jend, 1,1, wk3_out_stc(:, :, i), status)
         call netcdf_err(status, 'reading var: '//trim(stc_vars(i)), errflg, errmsg_out)
         if (errflg .ne. 0) return 
      else
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *, &
         'warning: no increment for ',trim(stc_vars(i)),' found, assuming zero'
         wk3_out_stc(:, :, i) = 0.
      endif
   enddo
   do i = 1, size(slc_vars)
      print *, trim(slc_vars(i))
      status = nf90_inq_varid(ncid, trim(slc_vars(i)), varid)
      if (status == nf90_noerr) then   !if (ierr == 0) then
         ! call get_var3_r4( ncid, trim(slc_vars(i)), 1,im, jbeg,jend, 1,1, wk3_out_slc(:, :, i) )
         call get_var3d_values(ncid, varid, 1,im, jbeg,jend, 1,1, wk3_out_slc(:, :, i), status)
         call netcdf_err(status, 'reading var: '//trim(slc_vars(i)), errflg, errmsg_out)
         if (errflg .ne. 0) return          
      else
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *,&
         'warning: no increment for ',trim(slc_vars(i)),' found, assuming zero'
         wk3_out_slc(:, :, i) = 0.
      endif
   enddo

   call close_ncfile(ncid)


end subroutine read_iau_forcing_fv3

subroutine interp_inc_at_timestep(Land_IAU_Control, km_in, wk3_in, var, errmsg, errflg)   !field_name, , jbeg, jend)
   ! interpolate increment from GSI gaussian grid to cubed sphere
   ! everying is on the A-grid, earth relative
   type (land_iau_control_type), intent(in) :: Land_IAU_Control
   ! character(len=*), intent(in) :: field_name
   integer,                                intent(in) :: km_in        !jbeg,jend
   real(kind=kind_phys),                           intent(in) :: wk3_in(1:im,jbeg:jend, 1:km_in) 
   real(kind=kind_phys), dimension(is:ie, js:je, 1:km), intent(inout) :: var
   
   character(len=*),              intent(inout) :: errmsg
   integer,                       intent(inout) :: errflg
   integer:: i1, i2, j1, k, j, i
   
   do k=1,km_in
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            var(i,j,k) = s2c(i,j,1)*wk3_in(i1,j1  ,k) + s2c(i,j,2)*wk3_in(i2,j1  ,k)+&
                        s2c(i,j,3)*wk3_in(i2,j1+1,k) + s2c(i,j,4)*wk3_in(i1,j1+1,k)
         enddo
      enddo
   enddo
end subroutine interp_inc_at_timestep

!> This subroutine is copied from 'fv_treat_da_inc.F90 by Xi.Chen <xi.chen@noaa.gov>
! copying it here, due to inability to 'include' from the original module when the land iau mod is called through CCPP frameowrk
!> @author Xi.Chen <xi.chen@noaa.gov> !> @date 02/12/2016
  !=============================================================================
  !>@brief The subroutine 'remap_coef' calculates the coefficients for horizonal regridding.
  subroutine remap_coef( is, ie, js, je, isd, ied, jsd, jed, &
      im, jm, lon, lat, id1, id2, jdc, s2c, agrid )

    integer, intent(in):: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(in):: im, jm
    real(kind=kind_phys),    intent(in):: lon(im), lat(jm)
    real(kind=kind_phys),    intent(out):: s2c(is:ie,js:je,4)
    integer, intent(out), dimension(is:ie,js:je):: id1, id2, jdc
    real(kind=kind_phys),    intent(in):: agrid(isd:ied,jsd:jed,2)
    ! local:
    real(kind=kind_phys) :: rdlon(im)
    real(kind=kind_phys) :: rdlat(jm)
    real(kind=kind_phys):: a1, b1
    integer i,j, i1, i2, jc, i0, j0
    do i=1,im-1
      rdlon(i) = 1. / (lon(i+1) - lon(i))
    enddo
    rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

    do j=1,jm-1
      rdlat(j) = 1. / (lat(j+1) - lat(j))
    enddo

    ! * Interpolate to cubed sphere cell center
    do 5000 j=js,je

      do i=is,ie

        if ( agrid(i,j,1)>lon(im) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
        elseif ( agrid(i,j,1)<lon(1) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
        else
          do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
              i1 = i0;  i2 = i0+1
              a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
              go to 111
            endif
          enddo
        endif
111     continue

        if ( agrid(i,j,2)<lat(1) ) then
          jc = 1
          b1 = 0.
        elseif ( agrid(i,j,2)>lat(jm) ) then
          jc = jm-1
          b1 = 1.
        else
          do j0=1,jm-1
            if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
              jc = j0
              b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
              go to 222
            endif
          enddo
        endif
222     continue

        if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
!TODO uncomment and fix mpp_pe             write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
        endif

        s2c(i,j,1) = (1.-a1) * (1.-b1)
        s2c(i,j,2) =     a1  * (1.-b1)
        s2c(i,j,3) =     a1  *     b1
        s2c(i,j,4) = (1.-a1) *     b1
        id1(i,j) = i1
        id2(i,j) = i2
        jdc(i,j) = jc
      enddo   !i-loop
5000 continue   ! j-loop

  end subroutine remap_coef

  !> Calculate soil mask for land on model grid.
!! Output is 1  - soil, 2 - snow-covered, 0 - land ice, -1  not land.
!!
!! @param[in] lensfc  Number of land points for this tile 
!! @param[in] veg_type_landice Value of vegetion class that indicates land-ice
!! @param[in] stype Soil type
!! @param[in] swe Model snow water equivalent
!! @param[in] vtype Model vegetation type
!! @param[out] mask Land mask for increments
!! @author Clara Draper @date March 2021
!! @author Yuan Xue: introduce stype to make the mask calculation more generic
subroutine calculate_landinc_mask(swe,vtype,stype,lensfc,veg_type_landice,mask)
 
   implicit none

   integer, intent(in)           :: lensfc, veg_type_landice
   real, intent(in)              :: swe(lensfc)
   integer, intent(in)           :: vtype(lensfc),stype(lensfc)
   integer, intent(out)          :: mask(lensfc)

   integer :: i

   mask = -1 ! not land

   ! land (but not land-ice)
   do i=1,lensfc
       if (stype(i) .GT. 0) then
         if (swe(i) .GT. 0.001) then ! snow covered land
               mask(i) = 2
         else                        ! non-snow covered land
               mask(i) = 1
         endif
       end if ! else should work here too
       if ( vtype(i) ==  veg_type_landice  ) then ! land-ice
               mask(i) = 0
       endif
   end do

end subroutine calculate_landinc_mask

  SUBROUTINE NETCDF_ERR(ERR, STRING, errflg, errmsg_out)

   !--------------------------------------------------------------
   ! IF AT NETCDF CALL RETURNS AN ERROR, PRINT OUT A MESSAGE
   ! AND STOP PROCESSING.
   !--------------------------------------------------------------
      IMPLICIT NONE

      include 'mpif.h'

      INTEGER, INTENT(IN) :: ERR
      CHARACTER(LEN=*), INTENT(IN) :: STRING
      CHARACTER(LEN=80) :: ERRMSG
      integer :: errflg
      character(len=*) :: errmsg_out

      !Errors messages handled through CCPP error handling variables
      errmsg_out = ''
      errflg = 0

      IF (ERR == NF90_NOERR) RETURN
      ERRMSG = NF90_STRERROR(ERR)
      PRINT*,'FATAL ERROR in Land IAU ', TRIM(STRING), ': ', TRIM(ERRMSG)
      errmsg_out = 'FATAL ERROR in Land IAU ', TRIM(STRING), ': ', TRIM(ERRMSG)
   !  CALL MPI_ABORT(MPI_COMM_WORLD, 999)
      errflg = 1
      return

   END SUBROUTINE NETCDF_ERR

   subroutine get_nc_dimlen(ncid, dim_name, dim_len, errflg, errmsg_out )
      integer, intent(in):: ncid
      character(len=*), intent(in)::  dim_name
      integer, intent(out):: dim_len
      integer :: dimid
      integer :: errflg
      character(len=*) :: errmsg_out

      !Errors messages handled through CCPP error handling variables
      errmsg_out = ''
      errflg = 0

      status = nf90_inq_dimid(ncid, dim_name, dimid)
      CALL netcdf_err(status, 'reading dim id '//trim(dim_name), errflg, errmsg_out)
      if (errflg .ne. 0) return
      status = nf90_inquire_dimension(ncid, dimid, len = dim_len)
      CALL netcdf_err(status, 'reading dim length '//trim(dim_name), errflg, errmsg_out)

   end subroutine get_nc_dimlen
      ! status = nf90_inq_dimid(ncid, "longitude", dimid)
      ! CALL netcdf_err(status, 'reading longitude dim id')
      ! status = nf90_inquire_dimension(ncid, dimid, len = im)
      ! CALL netcdf_err(status, 'reading dim longitude')
      ! status = nf90_inq_dimid(ncid, "latitude", dimid)
      ! CALL netcdf_err(status, 'reading latitude dim id')
      ! status = nf90_inquire_dimension(ncid, dimid, len = jm)
      ! CALL netcdf_err(status, 'reading dim latitude')
   subroutine get_var1d(ncid, dim_len, var_name, var_arr, errflg, errmsg_out)
      integer, intent(in):: ncid, dim_len
      character(len=*), intent(in)::  var_name
      real(kind=kind_phys), intent(out):: var_arr(dim_len)
      integer :: errflg
      character(len=*) :: errmsg_out
      integer :: var_id

      !Errors messages handled through CCPP error handling variables
      errmsg_out = ''
      errflg = 0

      status = nf90_inq_varid(ncid, trim(var_name), varid)
      CALL NETCDF_ERR(status, 'getting varid: '//trim(var_name), errflg, errmsg_out)
      if (errflg .ne. 0) return
      status = nf90_get_var(ncid, varid, var_arr)
      CALL NETCDF_ERR(status, 'reading var: '//trim(var_name), errflg, errmsg_out)

   end subroutine get_var1d

   subroutine get_var3d_values(ncid, varid, is,ie, js,je, ks,ke, var3d, status)
      integer, intent(in):: ncid, varid
      integer, intent(in):: is, ie, js, je, ks,ke
      real(kind=kind_phys), intent(out):: var3d(is:ie,js:je,ks:ke)
      integer, intent(out):: status 
      ! integer, dimension(3):: start, nreco
      ! start(1) = is; start(2) = js; start(3) = ks
      ! nreco(1) = ie - is + 1
      ! nreco(2) = je - js + 1
      ! nreco(3) = ke - ks + 1

      status = nf90_get_var(ncid, varid, var3d, &  !start = start, count = nreco)
               start = (/is, js, ks/), count = (/ie - is + 1, je - js + 1, ke - ks + 1/))

   end subroutine get_var3d_values
  
end module land_iau_mod



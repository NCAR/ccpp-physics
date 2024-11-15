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
!> - maps increments to FV3 grid points belonging to mpi process
!> - interpolates temporally (with filter-weights if required by configuration)
!> - updates states with the interpolated increments

!> March, 2024: Tseganeh Z. Gichamo, (EMC) based on the FV3 IAU mod 
!> by Xi.Chen <xi.chen@noaa.gov> and Philip Pegion, PSL <philip.pegion@noaa.gov>
!-------------------------------------------------------------------------------

!> \section arg_table_land_iau_mod Argument table                               
!! \htmlinclude land_iau_mod.html                                               
!!
module land_iau_mod
  
  use machine,                  only: kind_phys, kind_dyn
  use netcdf

  implicit none

  private

!> \section arg_table_land_iau_external_data_type Argument Table 
!! \htmlinclude land_iau_external_data_type.html
!!
  type land_iau_external_data_type
      real(kind=kind_phys),allocatable :: stc_inc(:,:,:)   
      real(kind=kind_phys),allocatable :: slc_inc(:,:,:)   
      logical                          :: in_interval = .false.
      real(kind=kind_phys)              :: hr1        
      real(kind=kind_phys)              :: hr2      
      real(kind=kind_phys)              :: wt
      real(kind=kind_phys)              :: wt_normfact
      real(kind=kind_phys)              :: rdt      
      integer                           :: itnext ! track the increment steps here
  end type land_iau_external_data_type

!!> \section arg_table_land_iau_state_type Argument Table
!! \htmlinclude land_iau_state_type.html
!!
  ! land_iau_state_type holds 'raw' (not interpolated) inrements, 
  ! read during land_iau_mod_init
  type land_iau_state_type
      real(kind=kind_phys),allocatable :: stc_inc(:,:,:,:)
      real(kind=kind_phys),allocatable :: slc_inc(:,:,:,:) 
  end type land_iau_state_type


!!!> \section arg_table_land_iau_control_type Argument Table
!! \htmlinclude land_iau_control_type.html
!!
  type land_iau_control_type      
      integer :: isc
      integer :: jsc
      integer :: nx
      integer :: ny
      integer :: tile_num
      integer :: nblks
      integer, allocatable :: blksz(:)        ! this could vary for the last block
      integer, allocatable :: blk_strt_indx(:)

      integer              :: lsoil           !< number of soil layers
      integer              :: lsnow_lsm       !< maximum number of snow layers internal to land surface model
      logical              :: do_land_iau
      real(kind=kind_phys) :: iau_delthrs     ! iau time interval (to scale increments) in hours
      character(len=240)   :: iau_inc_files(7)     ! list of increment files
      real(kind=kind_phys) :: iaufhrs(7)      ! forecast hours associated with increment files
      logical              :: iau_filter_increments   
      integer              :: lsoil_incr      ! soil layers (from top) updated by DA   
      logical              :: upd_stc
      logical              :: upd_slc
      logical              :: do_stcsmc_adjustment  !do moisture/temperature adjustment for consistency after increment add
      real(kind=kind_phys) :: min_T_increment
 
      integer              :: me              !< MPI rank designator
      integer              :: mpi_root        !< MPI rank of master atmosphere processor
      character(len=64)    :: fn_nml          !< namelist filename for surface data cycling
      real(kind=kind_phys) :: dtp             !< physics timestep in seconds
      real(kind=kind_phys) :: fhour           !< current forecast hour

      integer              :: ntimes

  end type land_iau_control_type

  public land_iau_control_type, land_iau_external_data_type, land_iau_state_type, land_iau_mod_set_control, &
         land_iau_mod_init, land_iau_mod_getiauforcing, land_iau_mod_finalize, calculate_landinc_mask

contains

subroutine land_iau_mod_set_control(Land_IAU_Control,fn_nml,input_nml_file_i, me, mpi_root, &
                                   isc, jsc, nx, ny, tile_num, nblks, blksz, &
                                   lsoil, lsnow_lsm, dtp, fhour, errmsg, errflg)          

   type (land_iau_control_type), intent(inout) :: Land_IAU_Control
   character(*), intent(in)                    :: fn_nml               !< namelist filename for surface data cycling
   character(len=:), intent(in), dimension(:), pointer :: input_nml_file_i
   integer, intent(in)                        :: me, mpi_root          !< MPI rank of master atmosphere processor   
   integer, intent(in)                        :: isc, jsc, nx, ny, tile_num, nblks, lsoil, lsnow_lsm
   integer, dimension(:),          intent(in) :: blksz                 !(one:) !GFS_Control%blksz
   real(kind=kind_phys), intent(in)           :: dtp                   !< physics timestep in seconds
   real(kind=kind_phys), intent(in)           :: fhour                 !< current forecast hour
   character(len=*),              intent(out) :: errmsg
   integer,                       intent(out) :: errflg
   
   integer                                    :: nb, ix
   integer                                    :: nlunit = 360          ! unit for namelist  !, intent(in)
   integer                                    :: ios
   logical                                    :: exists
   character(len=512)                         :: ioerrmsg

   character(len=:), pointer, dimension(:)    :: input_nml_file => null()
   character(len=4)                           :: iosstr

   !> land iau setting read from namelist
   logical               :: do_land_iau               = .false.
   real(kind=kind_phys)  :: land_iau_delthrs              = 0           !< iau time interval (to scale increments)
   character(len=240)    :: land_iau_inc_files(7)         = ''          !< list of increment files
   real(kind=kind_phys)  :: land_iau_fhrs(7)               = -1         !< forecast hours associated with increment files
   logical               :: land_iau_filter_increments    = .false.     !< filter IAU increments

   integer               :: lsoil_incr = 4
   logical               :: land_iau_upd_stc = .false.
   logical               :: land_iau_upd_slc = .false.
   logical               :: land_iau_do_stcsmc_adjustment = .false.
   real(kind=kind_phys)  :: land_iau_min_T_increment = 0.0001
  
   NAMELIST /land_iau_nml/ do_land_iau, land_iau_delthrs, land_iau_inc_files, land_iau_fhrs,   &  
                        land_iau_filter_increments, &  
                        lsoil_incr, land_iau_upd_stc, land_iau_upd_slc, land_iau_do_stcsmc_adjustment, land_iau_min_T_increment                                    
   
   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0

!3.11.24: copied from GFS_typedefs.F90 
#ifdef INTERNAL_FILE_NML
    ! allocate required to work around GNU compiler bug 100886
    ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100886
    allocate(input_nml_file, mold=input_nml_file_i)
    input_nml_file => input_nml_file_i
    read(input_nml_file, nml=land_iau_nml, ERR=888, END=999, iostat=ios)
#else
   inquire (file=trim(fn_nml), exist=exists)    ! TODO: this maybe be replaced by nlunit passed from ccpp
   if (.not. exists) then    
      errmsg = 'lnd_iau_mod_set_control: namelist file '//trim(fn_nml)//' does not exist'
      errflg = 1
      return
   else
      Land_IAU_Control%fn_nml = trim(fn_nml)  
      open (unit=nlunit, file=trim(fn_nml), action='READ', status='OLD', iostat=ios, iomsg=ioerrmsg)
      rewind(nlunit)
      read (nlunit, nml=land_iau_nml, ERR=888, END=999, iostat=ios)
      close (nlunit)
      if (ios /= 0) then      
         errmsg = 'lnd_iau_mod_set_control: error reading namelist file '//trim(fn_nml)  &
                  // 'the error message from file handler:' //trim(ioerrmsg) 
         errflg = 1
         return
       end if
   endif 
#endif
 
888 if (ios /= 0) then  
         write(iosstr, '(I0)') ios         
         errmsg = 'lnd_iau_mod_set_control: I/O error code '//trim(iosstr)//' at land_iau namelist read'  
         errflg = 1
         return
    end if 
       
999 if (ios /= 0) then  
        write(iosstr, '(I0)') ios
        if (me == mpi_root) then
          WRITE(6, * ) 'lnd_iau_mod_set_control: Warning! EoF ('//trim(iosstr)//') while reading land_iau namelist,' &
                  // ' likely because land_iau_nml was not found in input.nml. It will be set to default.' 
        endif
   endif

   if (me == mpi_root) then
      write(6,*) "land_iau_nml"
      write(6, land_iau_nml)
   endif
   
   Land_IAU_Control%do_land_iau = do_land_iau
   Land_IAU_Control%iau_delthrs = land_iau_delthrs
   Land_IAU_Control%iau_inc_files = land_iau_inc_files
   Land_IAU_Control%iaufhrs = land_iau_fhrs   
   Land_IAU_Control%iau_filter_increments = land_iau_filter_increments
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

   Land_IAU_Control%upd_stc = land_iau_upd_stc
   Land_IAU_Control%upd_slc = land_iau_upd_slc
   Land_IAU_Control%do_stcsmc_adjustment = land_iau_do_stcsmc_adjustment
   Land_IAU_Control%min_T_increment = land_iau_min_T_increment

   allocate(Land_IAU_Control%blksz(nblks))
   allocate(Land_IAU_Control%blk_strt_indx(nblks))

   ! Land_IAU_Control%blk_strt_indx = start index of each block, for flattened (ncol=nx*ny) arrays 
   ! It's required in noahmpdriv_run to get subsection of the stc array for each proces/thread
   ix = 1
   do nb=1, nblks
      Land_IAU_Control%blksz(nb) = blksz(nb)
      Land_IAU_Control%blk_strt_indx(nb) = ix
      ix = ix + blksz(nb)
   enddo

end subroutine land_iau_mod_set_control

subroutine land_iau_mod_init (Land_IAU_Control, Land_IAU_Data, Land_IAU_State, errmsg, errflg)    
   type (land_iau_control_type),       intent(inout) :: Land_IAU_Control
   type (land_iau_external_data_type), intent(inout) :: Land_IAU_Data
   type(land_iau_state_type),          intent(inout) :: Land_IAU_state  
   character(len=*),                   intent(  out) :: errmsg
   integer,                            intent(  out) :: errflg

   ! local
   character(len=128)   :: fname
   real(kind=kind_phys) :: sx, wx, wt, normfact, dtp
   integer              :: k, nstep, kstep
   integer              :: nfilesall, ntimesall
   integer, allocatable :: idt(:) 
   integer              :: nlon, nlat
   logical              :: exists
   integer              :: ncid, dimid, varid, status, IDIM
   
   real(kind=kind_phys) :: dt   !, rdt
   integer              :: im, jm, km, nfiles, ntimes

   integer :: is,  ie,  js,  je
   integer :: npz
   integer :: i, j  

   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0

   npz = Land_IAU_Control%lsoil
   km = Land_IAU_Control%lsoil
   
   is  = Land_IAU_Control%isc
   ie  = is + Land_IAU_Control%nx-1
   js  = Land_IAU_Control%jsc
   je  = js + Land_IAU_Control%ny-1
   nlon = Land_IAU_Control%nx
   nlat = Land_IAU_Control%ny

   ! allocate arrays that will hold iau state
   allocate(Land_IAU_Data%stc_inc(nlon, nlat, km))
   allocate(Land_IAU_Data%slc_inc(nlon, nlat, km))

   Land_IAU_Data%hr1=Land_IAU_Control%iaufhrs(1)
   Land_IAU_Data%wt = 1.0 ! IAU increment filter weights (default 1.0)
   Land_IAU_Data%wt_normfact = 1.0
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
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then 
            print *,'Land IAU init: IAU filter weights params k, kstep, wt ',k, kstep, wt
         endif
      enddo
      Land_IAU_Data%wt_normfact = (2*nstep+1)/normfact
   endif

   ! increment files are in fv3 tiles 
   if (trim(Land_IAU_Control%iau_inc_files(1)) .eq. '' .or. Land_IAU_Control%iaufhrs(1) .lt. 0) then ! only 1 file expected
      errmsg = "Error! in Land IAU init: increment file name is empty or iaufhrs(1) is negative"
      errflg = 1
      return   
   endif    
   if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then
      print*,"Land_iau_init: Increment file name: ", trim(adjustl(Land_IAU_Control%iau_inc_files(1)))
   endif

   ! determine number of valid forecast hours; read from the increment file ("Time" dim)
   if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then
      print *, "Land_iau_init: timesetps and forecast times (in hours) with valid increment values"
   endif
   ntimesall = size(Land_IAU_Control%iaufhrs)
   ntimes = 0
   do k=1,ntimesall
      if (Land_IAU_Control%iaufhrs(k) .lt. 0) exit
      if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then
         print *,k, " fhour ", Land_IAU_Control%iaufhrs(k)
      endif
      ntimes = ntimes + 1
   enddo

   Land_IAU_Control%ntimes = ntimes
   if (ntimes < 1) then
      errmsg = "Error! in Land IAU init: ntimes < 1 (no valid hour with increments); do_land_iau should not be .true."
      errflg = 1
      return
   endif
   if (ntimes > 1) then
      allocate(idt(ntimes-1))
      idt = Land_IAU_Control%iaufhrs(2:ntimes)-Land_IAU_Control%iaufhrs(1:ntimes-1)
      do k=1,ntimes-1
         if (idt(k) .ne. Land_IAU_Control%iaufhrs(2)-Land_IAU_Control%iaufhrs(1)) then
            errmsg = 'Fatal error in land_iau_init. forecast intervals in iaufhrs must be constant'
            errflg = 1
            return
         endif
      enddo
      deallocate(idt)
   endif
   dt = (Land_IAU_Control%iau_delthrs*3600.)
   Land_IAU_Data%rdt = 1.0/dt   !rdt
   
   ! Read all increment files at iau init time (at beginning of cycle) 
   ! increments are already in the fv3 grid--no need for interpolation     
   call read_iau_forcing_fv3(Land_IAU_Control, Land_IAU_state%stc_inc, Land_IAU_state%slc_inc, errmsg, errflg)  
   if (errflg .ne. 0) return

   if (ntimes.EQ.1) then  ! only need to get incrments once since constant forcing over window
      call setiauforcing(Land_IAU_Control, Land_IAU_Data, Land_IAU_state)
      Land_IAU_Data%itnext = 0
   endif
   if (ntimes.GT.1) then  !have increments at multiple forecast hours, 
      ! but only need 2 at a time and interpoalte for timesteps between them  
      ! interpolation is done in land_iau_mod_getiauforcing 
      Land_IAU_Data%hr2=Land_IAU_Control%iaufhrs(2)   
      Land_IAU_Data%itnext = 2
   endif

end subroutine land_iau_mod_init

subroutine land_iau_mod_finalize(Land_IAU_Control, Land_IAU_Data, Land_IAU_state, errmsg, errflg)

   implicit none

   type(land_iau_control_type),           intent(in) :: Land_IAU_Control
   type(land_iau_external_data_type),  intent(inout) :: Land_IAU_Data
   type(land_iau_state_type),          intent(inout) :: Land_IAU_state
   character(len=*),                     intent(out) :: errmsg
   integer,                              intent(out) :: errflg

   ! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   if (allocated(Land_IAU_Data%stc_inc)) deallocate (Land_IAU_Data%stc_inc)
   if (allocated(Land_IAU_Data%slc_inc)) deallocate (Land_IAU_Data%slc_inc)

   if (allocated(Land_IAU_state%stc_inc)) deallocate(Land_IAU_state%stc_inc)
   if (allocated(Land_IAU_state%slc_inc)) deallocate(Land_IAU_state%slc_inc)

end subroutine land_iau_mod_finalize

 subroutine land_iau_mod_getiauforcing(Land_IAU_Control, Land_IAU_Data, Land_IAU_State, errmsg, errflg)

   implicit none
   type(land_iau_control_type),        intent(inout) :: Land_IAU_Control
   type(land_iau_external_data_type),  intent(inout) :: Land_IAU_Data
   type(land_iau_state_type),             intent(in) :: Land_IAU_State
   character(len=*),                     intent(out) :: errmsg
   integer,                              intent(out) :: errflg
   real(kind=kind_phys) t1,t2,sx,wx,wt,dtp
   integer n,i,j,k,kstep,nstep     
   integer :: ntimes

    ! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   ntimes = Land_IAU_Control%ntimes

   Land_IAU_Data%in_interval=.false.
   if (ntimes.LE.0) then
       errmsg = 'called land_iau_mod_getiauforcing, but ntimes <=0, probably there is no increment file. Exiting.'
       errflg = 1
       return
   endif

   if (ntimes .eq. 1) then 
       t1 = Land_IAU_Control%iaufhrs(1)-0.5*Land_IAU_Control%iau_delthrs
       t2 = Land_IAU_Control%iaufhrs(1)+0.5*Land_IAU_Control%iau_delthrs
   else
       t1 = Land_IAU_Control%iaufhrs(1)
       t2 = Land_IAU_Control%iaufhrs(ntimes)
   endif
   if (Land_IAU_Control%iau_filter_increments) then
      ! compute increment filter weight
      ! t1 is beginning of window, t2 end of window, and Land_IAU_Control%fhour is current time
      ! in window kstep=-nstep,nstep (2*nstep+1 total) with time step of Land_IAU_Control%dtp
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
         Land_IAU_Data%wt = Land_IAU_Data%wt_normfact*wt
      else
         Land_IAU_Data%wt = 0.
      endif
   endif

   if (ntimes.EQ.1) then
   !  check to see if we are in the IAU window, no need to update the states since they are fixed over the window
      if ( Land_IAU_Control%fhour <= t1 .or. Land_IAU_Control%fhour > t2 ) then
         Land_IAU_Data%in_interval=.false.
      else
         Land_IAU_Data%in_interval=.true.
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then 
            print *,'land_iau_mod_getiauforcing: applying forcing at t for t1,t,t2,filter wt rdt ', &
            t1,Land_IAU_Control%fhour,t2,Land_IAU_Data%wt/Land_IAU_Data%wt_normfact,Land_IAU_Data%rdt
         endif         
         if (Land_IAU_Control%iau_filter_increments) call setiauforcing(Land_IAU_Control, Land_IAU_Data, Land_IAU_state)       
      endif
      return
   endif

   if (ntimes > 1) then
      if ( Land_IAU_Control%fhour <= t1 .or. Land_IAU_Control%fhour > t2 ) then
         Land_IAU_Data%in_interval=.false.
      else
         Land_IAU_Data%in_interval=.true.
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then 
            print *,'land_iau_mod_getiauforcing: applying forcing at t for t1,t,t2,filter wt rdt ', &
            t1,Land_IAU_Control%fhour,t2,Land_IAU_Data%wt/Land_IAU_Data%wt_normfact,Land_IAU_Data%rdt
         endif                 
         if (Land_IAU_Control%fhour > Land_IAU_Data%hr2) then ! need to read in next increment file
            Land_IAU_Data%itnext = Land_IAU_Data%itnext + 1
            Land_IAU_Data%hr1=Land_IAU_Data%hr2
            Land_IAU_Data%hr2=Land_IAU_Control%iaufhrs(Land_IAU_Data%itnext)
         endif
         
         call updateiauforcing(Land_IAU_Control, Land_IAU_Data, Land_IAU_State)
      endif
   endif

 end subroutine land_iau_mod_getiauforcing

subroutine updateiauforcing(Land_IAU_Control, Land_IAU_Data, Land_IAU_State)

   implicit none
   
   type (land_iau_control_type),          intent(in) :: Land_IAU_Control
   type(land_iau_external_data_type),  intent(inout) :: Land_IAU_Data
   type(land_iau_state_type),             intent(in) :: Land_IAU_State
   real(kind=kind_phys) delt_t  
   integer i,j,k
   integer :: is,  ie,  js,  je, npz, t1, t2
  
   t2 = Land_IAU_Data%itnext
   t1  = t2 - 1
   is  = 1  ! Land_IAU_Control%isc
   ie  = is + Land_IAU_Control%nx-1
   js  = 1  ! Land_IAU_Control%jsc
   je  = js + Land_IAU_Control%ny-1
   npz = Land_IAU_Control%lsoil

   delt_t = (Land_IAU_Data%hr2-(Land_IAU_Control%fhour))/(Land_IAU_Data%hr2-Land_IAU_Data%hr1)
   
   do j = js,je
      do i = is,ie
         do k = 1,npz  ! do k = 1,n_soill    !         
            Land_IAU_Data%stc_inc(i,j,k)  =(delt_t*Land_IAU_State%stc_inc(t1,i,j,k)  + (1.-delt_t)* Land_IAU_State%stc_inc(t2,i,j,k))*Land_IAU_Data%rdt*Land_IAU_Data%wt
            Land_IAU_Data%slc_inc(i,j,k)  =(delt_t*Land_IAU_State%slc_inc(t1,i,j,k)  + (1.-delt_t)* Land_IAU_State%slc_inc(t2,i,j,k))*Land_IAU_Data%rdt*Land_IAU_Data%wt
         end do
       enddo
   enddo
 end subroutine updateiauforcing

 subroutine setiauforcing(Land_IAU_Control, Land_IAU_Data, Land_IAU_State)

   implicit none
   type(land_iau_control_type),       intent(in   ) :: Land_IAU_Control
   type(land_iau_external_data_type), intent(inout) :: Land_IAU_Data
   type(land_iau_state_type),         intent(in   ) :: Land_IAU_State
   real(kind=kind_phys) delt
   integer i, j, k
   integer :: is,  ie,  js,  je, npz
   
   is  = 1 
   ie  = is + Land_IAU_Control%nx-1
   js  = 1 
   je  = js + Land_IAU_Control%ny-1
   npz = Land_IAU_Control%lsoil

   do j = js, je
      do i = is, ie
         do k = 1, npz          
            Land_IAU_Data%stc_inc(i,j,k) = Land_IAU_Data%wt*Land_IAU_State%stc_inc(1,i,j,k)*Land_IAU_Data%rdt
            Land_IAU_Data%slc_inc(i,j,k) = Land_IAU_Data%wt*Land_IAU_State%slc_inc(1,i,j,k)*Land_IAU_Data%rdt
         end do
      enddo
   enddo

 end subroutine setiauforcing

subroutine read_iau_forcing_fv3(Land_IAU_Control, wk3_stc, wk3_slc, errmsg, errflg)  

   type (land_iau_control_type),       intent(in) :: Land_IAU_Control
   real(kind=kind_phys), allocatable, intent(out) :: wk3_stc(:, :, :, :), wk3_slc(:, :, :, :)
   character(len=*),          intent(inout) :: errmsg
   integer,                   intent(inout) :: errflg

   integer  :: i, it, km 
   logical  :: exists
   integer  :: ncid, status, varid
   integer  :: ierr
   character(len=500)  :: fname
   character(len=2)    :: tile_str
   integer             :: n_t, n_y, n_x

   character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
   character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']
   character(len=32) :: slsn_mask = "soilsnow_mask"

   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0

   km = Land_IAU_Control%lsoil

   write(tile_str, '(I0)') Land_IAU_Control%tile_num

   fname = 'INPUT/'//trim(Land_IAU_Control%iau_inc_files(1))//".tile"//trim(tile_str)//".nc"
   
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
   ! var stored as soilt1_inc(Time, yaxis_1, xaxis_1)
   call get_nc_dimlen(ncid, "Time", n_t, errflg, errmsg) 
   if (errflg .ne. 0) return
   call get_nc_dimlen(ncid, "yaxis_1", n_y, errflg, errmsg) 
   if (errflg .ne. 0) return
   call get_nc_dimlen(ncid, "xaxis_1", n_x, errflg, errmsg) 
   if (errflg .ne. 0) return

   if (n_x .lt. Land_IAU_Control%nx) then 
      errmsg = 'Error in land iau read_iau_forcing_fv3: Land_IAU_Control%nx bigger than dim xaxis_1 in file '//trim(fname)
      errflg = 1
      return
   endif
   if (n_y .lt. Land_IAU_Control%ny) then 
      errmsg = 'Error in land iau read_iau_forcing_fv3: Land_IAU_Control%ny bigger than dim yaxis_1 in file '//trim(fname)
      errflg = 1
      return
   endif

   allocate(wk3_stc(n_t, Land_IAU_Control%nx, Land_IAU_Control%ny, km))
   allocate(wk3_slc(n_t, Land_IAU_Control%nx, Land_IAU_Control%ny, km))
   
   do i = 1, size(stc_vars)
      if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *, trim(stc_vars(i))
      status = nf90_inq_varid(ncid, trim(stc_vars(i)), varid)
      if (status == nf90_noerr) then   
         do it = 1, n_t
            ! var stored as soilt1_inc(Time, yaxis_1, xaxis_1)
            call get_var3d_values(ncid, varid, trim(stc_vars(i)), Land_IAU_Control%isc, Land_IAU_Control%nx, Land_IAU_Control%jsc, Land_IAU_Control%ny, &
                 it, 1, wk3_stc(it,:, :, i), status, errflg, errmsg)
            if (errflg .ne. 0) return 
         enddo
      else
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then 
            print *, 'warning! No increment for ',trim(stc_vars(i)),' found, assuming zero'
         endif
         wk3_stc(:, :, :, i) = 0.
      endif
   enddo
   do i = 1, size(slc_vars)
      if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) print *, trim(slc_vars(i))
      status = nf90_inq_varid(ncid, trim(slc_vars(i)), varid)
      if (status == nf90_noerr) then   
         do it = 1, n_t
            call get_var3d_values(ncid, varid, trim(slc_vars(i)), Land_IAU_Control%isc, Land_IAU_Control%nx, Land_IAU_Control%jsc, Land_IAU_Control%ny, &
                 it, 1, wk3_slc(it, :, :, i), status, errflg, errmsg)
            if (errflg .ne. 0) return   
         end do       
      else
         if (Land_IAU_Control%me == Land_IAU_Control%mpi_root) then 
            print *, 'warning! No increment for ',trim(slc_vars(i)),' found, assuming zero'
         endif
         wk3_slc(:, :, :, i) = 0.
      endif
   enddo
   
   !set too small increments to zero
   where(abs(wk3_stc) < Land_IAU_Control%min_T_increment) wk3_stc = 0.0

   status =nf90_close(ncid) 
   call netcdf_err(status, 'closing file '//trim(fname), errflg, errmsg) 

end subroutine read_iau_forcing_fv3

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
   subroutine calculate_landinc_mask(swe,vtype,stype,lensfc,veg_type_landice, mask)
   
      implicit none

      integer, intent(in)           :: lensfc, veg_type_landice
      real(kind=kind_phys), intent(in)              :: swe(lensfc)
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

   subroutine netcdf_err(ERR, STRING, errflg, errmsg_out)

   !--------------------------------------------------------------
   ! Process the error flag from a NETCDF call and return it as (human readable) MESSAGE
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
      errmsg_out = 'FATAL ERROR in Land IAU '//TRIM(STRING)//': '//TRIM(ERRMSG)
      errflg = 1
      return

   end subroutine netcdf_err

   subroutine get_nc_dimlen(ncid, dim_name, dim_len, errflg, errmsg_out )
      integer, intent(in):: ncid
      character(len=*), intent(in)::  dim_name
      integer, intent(out):: dim_len
      integer :: dimid
      integer :: errflg
      character(len=*) :: errmsg_out
      integer :: status

      !Errors messages handled through CCPP error handling variables
      errmsg_out = ''
      errflg = 0

      status = nf90_inq_dimid(ncid, dim_name, dimid)
      CALL netcdf_err(status, 'reading dim id '//trim(dim_name), errflg, errmsg_out)
      if (errflg .ne. 0) return
      status = nf90_inquire_dimension(ncid, dimid, len = dim_len)
      CALL netcdf_err(status, 'reading dim length '//trim(dim_name), errflg, errmsg_out)

   end subroutine get_nc_dimlen

   subroutine get_var1d(ncid, dim_len, var_name, var_arr, errflg, errmsg_out)
      integer, intent(in):: ncid, dim_len
      character(len=*), intent(in)::  var_name
      real(kind=kind_phys), intent(out):: var_arr(dim_len)
      integer :: errflg
      character(len=*) :: errmsg_out
      integer :: varid, status

      !Errors messages handled through CCPP error handling variables
      errmsg_out = ''
      errflg = 0

      status = nf90_inq_varid(ncid, trim(var_name), varid)
      call netcdf_err(status, 'getting varid: '//trim(var_name), errflg, errmsg_out)
      if (errflg .ne. 0) return

      status = nf90_get_var(ncid, varid, var_arr)
      call netcdf_err(status, 'reading var: '//trim(var_name), errflg, errmsg_out)

   end subroutine get_var1d

   subroutine get_var3d_values(ncid, varid, var_name, is,ix, js,jy, ks,kz, var3d, status, errflg, errmsg_out)
      integer, intent(in):: ncid, varid
      integer, intent(in):: is, ix, js, jy, ks,kz
      character(len=*), intent(in)::  var_name
      real(kind=kind_phys), intent(out):: var3d(ix, jy, kz)   
      integer, intent(out):: status 
      integer :: errflg
      character(len=*) :: errmsg_out

      !Errors messages handled through CCPP error handling variables
      errmsg_out = ''
      errflg = 0

      status = nf90_get_var(ncid, varid, var3d, &  
               start = (/is, js, ks/), count = (/ix, jy, kz/))

      call netcdf_err(status, 'get_var3d_values '//trim(var_name), errflg, errmsg_out)
      

   end subroutine get_var3d_values

   subroutine get_var3d_values_int(ncid, varid, var_name, is,ix, js,jy, ks,kz, var3d, status, errflg, errmsg_out)
      integer, intent(in):: ncid, varid
      integer, intent(in):: is, ix, js, jy, ks,kz
      character(len=*), intent(in)::  var_name
      integer, intent(out):: var3d(ix, jy, kz)   
      integer, intent(out):: status 
      integer :: errflg
      character(len=*) :: errmsg_out

      !Errors messages handled through CCPP error handling variables
      errmsg_out = ''
      errflg = 0

      status = nf90_get_var(ncid, varid, var3d, &  !start = start, count = nreco)
               start = (/is, js, ks/), count = (/ix, jy, kz/))
      
      call netcdf_err(status, 'get_var3d_values_int '//trim(var_name), errflg, errmsg_out)

   end subroutine get_var3d_values_int
  
end module land_iau_mod



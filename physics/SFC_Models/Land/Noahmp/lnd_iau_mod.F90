!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!> The routine 'remapcoeff is copied from 'fv_treat_da_inc.F90 by Xi.Chen <xi.chen@noaa.gov> 
! and put at the end of this module because, due to the compile order in CCPP framework it wasn't possible to 'include' 
! the original module when the land iau mod is called through CCPP frameowrk
!


!-------------------------------------------------------------------------------
!> @brief incremental analysis update module
!> @author Xi.Chen - author of fv_treat_da_inc.F90
!> @author Philip Pegion <philip.pegion@noaa.gov>
!> @date 09/13/2017
!
!>  REVISION HISTORY:
!>  09/13/2017 - Initial Version based on fv_treat_da_inc.F90
!-------------------------------------------------------------------------------

#ifdef OVERLOAD_R4
#define _GET_VAR1 get_var1_real
#else
#define _GET_VAR1 get_var1_double
#endif

module lnd_iau_mod

!   use fms_mod,             only: file_exist
!   use mpp_mod,             only: mpp_error, FATAL, NOTE, mpp_pe
!   use mpp_domains_mod,     only: domain2d
!   use constants_mod,       only: pi=>pi_8
!   use fv_arrays_mod,       only: R_GRID         !, &
                                 ! fv_atmos_type,       &
                                 ! fv_grid_type,        &
                                 ! fv_grid_bounds_type, &                                 
!   use fv_mp_mod,           only: is_master
  use sim_nc_mod_lnd,          only: open_ncfile,         &
                                 close_ncfile,        &
                                 get_ncdim1,          &
                                 get_var1_double,     &
                                 get_var3_r4,         &
                                 get_var1_real, check_var_exists
! #ifdef GFS_TYPES
!   use GFS_typedefs,        only: IPD_init_type => GFS_init_type, &
!                                  LND_IAU_Control_type => GFS_control_type, &
!                                  kind_phys, &
!                                  IPD_Data_type => GFS_data_type
! #else
!   use IPD_typedefs,        only: IPD_init_type, LND_IAU_Control_type, &
!                                  kind_phys => IPD_kind_phys
! #endif

!   use block_control_mod,   only: block_control_type
!   use fv_treat_da_inc_mod, only: remap_coef
!   use tracer_manager_mod,  only: get_tracer_names,get_tracer_index, get_number_tracers
!   use field_manager_mod,   only: MODEL_ATMOS
  
  use machine,                  only: kind_phys, kind_dyn
  use physcons,                 only: pi => con_pi

  implicit none

  private

  real,allocatable::s2c(:,:,:)
!  real:: s2c(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je,4)
!  integer, dimension(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je):: &
!      id1, id2, jdc
  integer,allocatable,dimension(:,:) :: id1,id2,jdc

  real :: deg2rad,dt,rdt
  integer :: im,jm,km,nfiles,ncid
  integer:: jbeg, jend

  integer :: n_soill, n_snowl              !1.27.24 soil and snow layers
  logical :: do_lnd_iau_inc   !do_lnd_iau_inc

  integer :: is,  ie,  js,  je
  integer :: npz     !, ntracers
!   character(len=32), allocatable :: tracer_names(:)
!   integer, allocatable :: tracer_indicies(:)

!   real(kind=4), allocatable:: wk3(:, :,:,:)
  real(kind=4), allocatable:: wk3_stc(:, :, :, :), wk3_slc(:, :, :, :), wk3_t2m(:, :, :, :), wk3_q2m(:, :, :, :)

  type iau_internal_data_type
      ! real,allocatable :: ua_inc(:,:,:)
      ! real,allocatable :: va_inc(:,:,:)
      ! real,allocatable :: temp_inc(:,:,:)
      ! real,allocatable :: delp_inc(:,:,:)
      ! real,allocatable :: delz_inc(:,:,:)
      ! real,allocatable :: tracer_inc(:,:,:,:)
      real,allocatable :: stc_inc(:,:,:)   
      real,allocatable :: slc_inc(:,:,:) 
      real,allocatable :: tmp2m_inc(:,:, :) 
      real,allocatable :: spfh2m_inc(:,:, :) 
  end type iau_internal_data_type

  type lnd_iau_external_data_type
      real,allocatable :: stc_inc(:,:,:)   
      real,allocatable :: slc_inc(:,:,:) 
      real,allocatable :: tmp2m_inc(:,:,:) 
      real,allocatable :: spfh2m_inc(:,:,:)    
      logical          :: in_interval = .false.
      ! logical          :: drymassfixer = .false.
  end type lnd_iau_external_data_type

  type iau_state_type
      type(iau_internal_data_type):: inc1
      type(iau_internal_data_type):: inc2
      real(kind=kind_phys)        :: hr1
      real(kind=kind_phys)        :: hr2
      real(kind=kind_phys)        :: wt
      real(kind=kind_phys)        :: wt_normfact
  end type iau_state_type

  type lnd_iau_control_type      
      integer :: isc
      integer :: jsc
      integer :: nx
      integer :: ny
      integer :: nblks
      ! integer :: blksz   ! this could vary for the last block
      integer, allocatable :: blksz(:)
      integer, allocatable :: blk_strt_indx(:)

      integer :: lsoil  !< number of soil layers
      ! this is the max dim (TBC: check it is consitent for noahmpdrv)
      integer              :: lsnow_lsm       !< maximum number of snow layers internal to land surface model
      logical              :: do_lnd_iau_inc
      real(kind=kind_phys) :: iau_delthrs     ! iau time interval (to scale increments) in hours
      character(len=240)   :: iau_inc_files(7)! list of increment files
      real(kind=kind_phys) :: iaufhrs(7)      ! forecast hours associated with increment files
      logical              :: iau_filter_increments          
      !, iau_drymassfixer
      integer              :: me              !< MPI rank designator
      integer              :: mpi_root          !< MPI rank of master atmosphere processor
      character(len=64)    :: fn_nml          !< namelist filename for surface data cycling
      real(kind=kind_phys) :: dtp             !< physics timestep in seconds
      real(kind=kind_phys) :: fhour           !< current forecast hour
      character(len=:), pointer, dimension(:) :: input_nml_file => null() !<character string containing full namelist
                                                                          !< for use with internal file reads
      integer              :: input_nml_file_length    !<length (number of lines) in namelist for internal reads
 

      ! Additional GFS_Control vars (not used currently)      
      ! integer              :: iau_offset
      ! integer              :: tile_num
      ! integer              :: nblks           !< for explicit data blocking: number of blocks
      ! integer,     pointer :: blksz(:)        !< for explicit data blocking: block sizes of all blocks
      ! integer              :: ncols           !< total number of columns for all blocks   
      ! integer              :: communicator    !< MPI communicator
      ! integer              :: ntasks          !< MPI size in communicator
      ! integer              :: nthreads        !< OpenMP threads available for physics
      ! integer              :: nlunit          !< unit for namelist
      !  character(len=:), pointer, dimension(:) :: input_nml_file => null() !< character string containing full namelist
      !  integer              :: logunit
      !--- calendars and time parameters and activation triggers      
      ! real(kind=kind_phys) :: dtf             !< dynamics timestep in seconds    
      ! integer              :: idat(1:8)       !< initialization date and time
      !                                        !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)
      ! integer              :: jdat(1:8)       !< current forecast date and time
      !                                        !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)
      ! real(kind=kind_phys)          :: sec    !< seconds since model initialization
      ! real(kind=kind_phys) :: phour           !< previous forecast hour      
      ! real(kind=kind_phys) :: zhour           !< previous hour diagnostic buckets emptied
      ! integer              :: kdt             !< current forecast iteration
      ! logical              :: first_time_step !< flag signaling first time step for time integration routine
  end type lnd_iau_control_type

  type(iau_state_type) :: IAU_state
  public lnd_iau_control_type, lnd_iau_external_data_type, lnd_iau_mod_set_control, lnd_iau_mod_init, lnd_iau_mod_getiauforcing

contains

subroutine lnd_iau_mod_set_control(LND_IAU_Control,fn_nml,input_nml_file_i,me, mpi_root, isc, jsc, nx, ny, nblks, blksz, &
                                   lsoil, lsnow_lsm, dtp, fhour, errmsg, errflg)          !nlunit

   type (lnd_iau_control_type), intent(inout) :: LND_IAU_Control
   character(*), intent(in)                   :: fn_nml          !< namelist filename for surface data cycling
   character(len=:), intent(in), dimension(:), pointer :: input_nml_file_i
   integer, intent(in)                        :: me, mpi_root          !< MPI rank of master atmosphere processor   
   integer, intent(in)                        :: isc, jsc, nx, ny, nblks, lsoil, lsnow_lsm
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


   !> 3.9.24 these are not available through the CCPP interface so need to read them from namelist file
   !> vars to read from namelist
   logical               :: do_lnd_iau_inc               = .false.
   real(kind=kind_phys)  :: lnd_iau_delthrs              = 0           !< iau time interval (to scale increments)
   character(len=240)    :: lnd_iau_inc_files(7)         = ''          !< list of increment files
   real(kind=kind_phys)  :: lnd_iaufhrs(7)               = -1          !< forecast hours associated with increment files
   logical               :: lnd_iau_filter_increments    = .false.     !< filter IAU increments
  
   NAMELIST /lnd_iau_nml/ do_lnd_iau_inc, lnd_iau_delthrs, lnd_iau_inc_files, lnd_iaufhrs, lnd_iau_filter_increments  !, lnd_iau_drymassfixer                                          &
   
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
      LND_IAU_Control%fn_nml = trim(fn_nml)   ! maynot need this
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
   
   LND_IAU_Control%do_lnd_iau_inc = do_lnd_iau_inc
   LND_IAU_Control%iau_delthrs = lnd_iau_delthrs
   LND_IAU_Control%iau_inc_files = lnd_iau_inc_files
   LND_IAU_Control%iaufhrs = lnd_iaufhrs   
   LND_IAU_Control%iau_filter_increments = lnd_iau_filter_increments
   ! LND_IAU_Control%iau_drymassfixer = lnd_iau_drymassfixer
   LND_IAU_Control%me = me
   LND_IAU_Control%mpi_root = mpi_root
   LND_IAU_Control%isc = isc
   LND_IAU_Control%jsc = jsc
   LND_IAU_Control%nx = nx
   LND_IAU_Control%ny = ny
   LND_IAU_Control%nblks = nblks
   LND_IAU_Control%lsoil = lsoil
   LND_IAU_Control%lsnow_lsm = lsnow_lsm
   LND_IAU_Control%dtp = dtp
   LND_IAU_Control%fhour = fhour

   LND_IAU_Control%input_nml_file = input_nml_file
   LND_IAU_Control%input_nml_file_length = input_nml_file_length

   allocate(LND_IAU_Control%blksz(nblks))
   allocate(LND_IAU_Control%blk_strt_indx(nblks))
   !start index of each block, for flattened (ncol=nx*ny) arrays 
   ! required in noahmpdriv_run to get subsection of the stc array for each
   ! proc/thread
   ix = 1
   do nb=1, nblks
      LND_IAU_Control%blksz(nb) = blksz(nb)
      LND_IAU_Control%blk_strt_indx(nb) = ix
      ix = ix + blksz(nb)
   enddo

end subroutine lnd_iau_mod_set_control

subroutine lnd_iau_mod_init (LND_IAU_Control, LND_IAU_Data, xlon, xlat, errmsg, errflg)     !nlunit, ncols, IPD_Data,,Init_parm)
   ! integer,                              intent(in) :: me, mpi_root
   type (lnd_iau_control_type),          intent(in) :: LND_IAU_Control
   type (lnd_iau_external_data_type), intent(inout) :: LND_IAU_Data
   ! type (IPD_init_type),    intent(in) :: Init_parm
   ! type (IPD_Data_type), dimension(:),    intent(in) :: IPD_Data
   ! integer, intent(in)                          :: ncols   
   real(kind_phys), dimension(:), intent(in)  :: xlon    ! longitude  !GFS_Data(cdata%blk_no)%Grid%xlon
   real(kind_phys), dimension(:), intent(in)  :: xlat    ! latitude
   character(len=*),              intent(out) :: errmsg
   integer,                       intent(out) :: errflg

   ! local
   character(len=128) :: fname
   ! real, dimension(:,:,:), allocatable:: u_inc, v_inc
   real(kind=kind_dyn), allocatable:: lat(:), lon(:),agrid(:,:,:)
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

   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0

   do_lnd_iau_inc = LND_IAU_Control%do_lnd_iau_inc
   n_soill = LND_IAU_Control%lsoil     !4  for sfc updates
!  n_snowl = LND_IAU_Control%lsnowl 
   npz = LND_IAU_Control%lsoil
   
   is  = LND_IAU_Control%isc
   ie  = is + LND_IAU_Control%nx-1
   js  = LND_IAU_Control%jsc
   je  = js + LND_IAU_Control%ny-1
   nlon = LND_IAU_Control%nx
   nlat = LND_IAU_Control%ny
   !nblks = LND_IAU_Control%nblks
   !blksz = LND_IAU_Control%blksz(1)

   allocate(Init_parm_xlon(nlon,nlat), Init_parm_xlat(nlon,nlat))
   ib = 1
   do j = 1, nlat  !ny
      ! do i = 1, nx      
         Init_parm_xlon (:,j) = xlon(ib:ib+nlon-1)  
         Init_parm_xlat (:,j) = xlat(ib:ib+nlon-1) 
         ib = ib+nlon
      ! enddo
   enddo
   !  call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   !  allocate (tracer_names(ntracers))
   !  allocate (tracer_indicies(ntracers))
   !  do i = 1, ntracers
   !     call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
   !     tracer_indicies(i)  = get_tracer_index(MODEL_ATMOS,tracer_names(i))
   !  enddo
   allocate(s2c(is:ie,js:je,4))
   allocate(id1(is:ie,js:je))
   allocate(id2(is:ie,js:je))
   allocate(jdc(is:ie,js:je))
   allocate(agrid(is:ie,js:je,2))
! determine number of increment files to read, and the valid forecast hours

   nfilesall = size(LND_IAU_Control%iau_inc_files)
   nfiles = 0
   if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print*,'in lnd_iau_init incfile1 iaufhr1 ', &
                                 trim(LND_IAU_Control%iau_inc_files(1)),LND_IAU_Control%iaufhrs(1)
   do k=1,nfilesall
      if (trim(LND_IAU_Control%iau_inc_files(k)) .eq. '' .or. LND_IAU_Control%iaufhrs(k) .lt. 0) exit   
      if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) then
         print *,k, " ", trim(adjustl(LND_IAU_Control%iau_inc_files(k)))
      endif
      nfiles = nfiles + 1
   enddo
   if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'nfiles = ',nfiles
   if (nfiles < 1) then
      return
   endif
   if (nfiles > 1) then
      allocate(idt(nfiles-1))
      idt = LND_IAU_Control%iaufhrs(2:nfiles)-LND_IAU_Control%iaufhrs(1:nfiles-1)
      do k=1,nfiles-1
         if (idt(k) .ne. LND_IAU_Control%iaufhrs(2)-LND_IAU_Control%iaufhrs(1)) then
           print *,'in lnd_iau_init: forecast intervals in iaufhrs must be constant'
         !   call mpp_error (FATAL,' forecast intervals in iaufhrs must be constant')
           errmsg = 'Fatal error in lnd_iau_init. forecast intervals in iaufhrs must be constant'
           errflg = 1
           return
         endif
      enddo
      deallocate(idt)
   endif
   if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'lnd_iau interval = ',LND_IAU_Control%iau_delthrs,' hours'
   dt = (LND_IAU_Control%iau_delthrs*3600.)
   rdt = 1.0/dt

!  set up interpolation weights to go from GSI's gaussian grid to cubed sphere
   deg2rad = pi/180.

   !  npz = LND_IAU_Control%levs
   fname = 'INPUT/'//trim(LND_IAU_Control%iau_inc_files(1))    
   inquire (file=trim(fname), exist=exists)    
   if (exists) then
   !  if( file_exist(fname) ) then
      call open_ncfile( fname, ncid )        ! open the file
!TODO !change to Latitude
      call get_ncdim1( ncid, 'longitude',   im)    
      call get_ncdim1( ncid, 'latitude',   jm)
      ! call get_ncdim1( ncid, 'nsoill',   km)
      km = n_soill
      ! if (km.ne.npz) then
      !   if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *, 'km = ', km
      ! !   call mpp_error(FATAL, '==> Error in IAU_initialize: km is not equal to npz')
      !   errmsg = 'Fatal Error in IAU_initialize: km is not equal to npz'
      !   errflg = 1
      !   return
      ! endif
      if(LND_IAU_Control%me == LND_IAU_Control%mpi_root)  write(*,*) fname, ' DA increment dimensions:', im,jm,km

      allocate (  lon(im) )
      allocate (  lat(jm) )

      call _GET_VAR1 (ncid, 'longitude', im, lon )
      call _GET_VAR1 (ncid, 'latitude', jm, lat )
      call close_ncfile(ncid)

      ! Convert to radians
      do i=1,im
        lon(i) = lon(i) * deg2rad
      enddo
      do j=1,jm
        lat(j) = lat(j) * deg2rad
      enddo
   else
      ! call mpp_error(FATAL,'==> Error in IAU_initialize: Expected file '&
      !     //trim(fname)//' for DA increment does not exist')
      errmsg = 'FATAL Error in IAU_initialize: Expected file '// trim(fname)//' for DA increment does not exist'
      errflg = 1
      return
   endif

! Initialize lat-lon to Cubed bi-linear interpolation coeff:
! populate agrid
!    print*,'is,ie,js,je=',is,ie,js,ie
!    print*,'size xlon=',size(Init_parm%xlon(:,1)),size(Init_parm%xlon(1,:))
!    print*,'size agrid=',size(agrid(:,1,1)),size(agrid(1,:,1)),size(agrid(1,1,:))
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
   deallocate ( lon, lat,agrid )
   if (allocated(Init_parm_xlon)) deallocate(Init_parm_xlon)
   if (allocated(Init_parm_xlat)) deallocate(Init_parm_xlat)

   !  allocate(LND_IAU_Data%ua_inc(is:ie, js:je, km))
   !  allocate(LND_IAU_Data%va_inc(is:ie, js:je, km))
   !  allocate(LND_IAU_Data%temp_inc(is:ie, js:je, km))
   !  allocate(LND_IAU_Data%delp_inc(is:ie, js:je, km))
   !  allocate(LND_IAU_Data%delz_inc(is:ie, js:je, km))
   !  allocate(LND_IAU_Data%tracer_inc(is:ie, js:je, km,ntracers))
   allocate(LND_IAU_Data%stc_inc(is:ie, js:je, km))
   allocate(LND_IAU_Data%slc_inc(is:ie, js:je, km))
   allocate(LND_IAU_Data%tmp2m_inc(is:ie, js:je, 1))
   allocate(LND_IAU_Data%spfh2m_inc(is:ie, js:je, 1))
! allocate arrays that will hold iau state
   allocate (iau_state%inc1%stc_inc(is:ie, js:je, km))
   allocate (iau_state%inc1%slc_inc(is:ie, js:je, km))
   allocate (iau_state%inc1%tmp2m_inc(is:ie, js:je, 1))
   allocate (iau_state%inc1%spfh2m_inc (is:ie, js:je, 1))
   iau_state%hr1=LND_IAU_Control%iaufhrs(1)
   iau_state%wt = 1.0 ! IAU increment filter weights (default 1.0)
   iau_state%wt_normfact = 1.0
   if (LND_IAU_Control%iau_filter_increments) then
      ! compute increment filter weights, sum to obtain normalization factor
      dtp=LND_IAU_Control%dtp
      nstep = 0.5*LND_IAU_Control%iau_delthrs*3600/dtp
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
         if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'filter wts',k,kstep,wt
      enddo
      iau_state%wt_normfact = (2*nstep+1)/normfact
   endif

!3.22.24 MB wants to read all increments files at iau init 
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
 
   ! call read_iau_forcing(LND_IAU_Control,iau_state%inc1,'INPUT/'//trim(LND_IAU_Control%iau_inc_files(1)), errmsg, errflg)  
   allocate (wk3_stc(nfiles, 1:im,jbeg:jend, 1:km))
   allocate (wk3_slc(nfiles, 1:im,jbeg:jend, 1:km))
   allocate (wk3_t2m(nfiles, 1:im,jbeg:jend, 1:1))
   allocate (wk3_q2m(nfiles, 1:im,jbeg:jend, 1:1))
   do k=1, nfiles
      call read_iau_forcing_all_timesteps(LND_IAU_Control, 'INPUT/'//trim(LND_IAU_Control%iau_inc_files(k)), errmsg, errflg, &
                                          wk3_stc(k, :, :, :), wk3_slc(k, :, :, :), wk3_t2m(k, :, :, :), wk3_q2m(k, :, :, :)) 
   enddo
   ! call interp_inc(LND_IAU_Control, 'soilt1_inc',increments%stc_inc(:,:,1),jbeg,jend)
   ! call interp_inc(LND_IAU_Control, 'tmp2m_inc',increments%tmp2m_inc(:,:,1),jbeg,jend)
   call interp_inc_at_timestep(LND_IAU_Control, km, wk3_stc(1, :, :, :), iau_state%inc1%stc_inc, errmsg, errflg)
   call interp_inc_at_timestep(LND_IAU_Control, km, wk3_slc(1, :, :, :), iau_state%inc1%slc_inc, errmsg, errflg)
   call interp_inc_at_timestep(LND_IAU_Control, 1, wk3_t2m(1, :, :, :), iau_state%inc1%tmp2m_inc, errmsg, errflg)
   call interp_inc_at_timestep(LND_IAU_Control, 1, wk3_q2m(1, :, :, :), iau_state%inc1%spfh2m_inc, errmsg, errflg)

   if (nfiles.EQ.1) then  ! only need to get incrments once since constant forcing over window
      call setiauforcing(LND_IAU_Control, LND_IAU_Data, iau_state%wt)
   endif
   if (nfiles.GT.1) then  !have multiple files, but only read in 2 at a time and interpoalte between them
      allocate (iau_state%inc2%stc_inc(is:ie, js:je, km))
      allocate (iau_state%inc2%slc_inc(is:ie, js:je, km))
      allocate (iau_state%inc2%tmp2m_inc(is:ie, js:je, 1))
      allocate (iau_state%inc2%spfh2m_inc(is:ie, js:je, 1))
      iau_state%hr2=LND_IAU_Control%iaufhrs(2)

      ! call read_iau_forcing(LND_IAU_Control,iau_state%inc2,'INPUT/'//trim(LND_IAU_Control%iau_inc_files(2)), errmsg, errflg)
      call interp_inc_at_timestep(LND_IAU_Control, km, wk3_stc(2, :, :, :), iau_state%inc2%stc_inc, errmsg, errflg)
      call interp_inc_at_timestep(LND_IAU_Control, km, wk3_slc(2, :, :, :), iau_state%inc2%slc_inc, errmsg, errflg)
      call interp_inc_at_timestep(LND_IAU_Control, 1, wk3_t2m(2, :, :, :), iau_state%inc2%tmp2m_inc, errmsg, errflg)
      call interp_inc_at_timestep(LND_IAU_Control, 1, wk3_q2m(2, :, :, :), iau_state%inc2%spfh2m_inc, errmsg, errflg)         
   endif
!   print*,'end of IAU init',dt,rdt

end subroutine lnd_iau_mod_init

subroutine lnd_iau_mod_finalize()

   implicit none

   if (allocated (wk3_stc)) deallocate (wk3_stc)
   if (allocated (wk3_slc)) deallocate (wk3_slc)
   if (allocated (wk3_t2m)) deallocate (wk3_t2m)
   if (allocated (wk3_q2m)) deallocate (wk3_q2m)

   if (allocated(LND_IAU_Data%stc_inc)) deallocate (LND_IAU_Data%stc_inc)
   if (allocated(LND_IAU_Data%slc_inc)) deallocate (LND_IAU_Data%slc_inc)
   if (allocated(LND_IAU_Data%tmp2m_inc)) deallocate (LND_IAU_Data%tmp2m_inc)
   if (allocated(LND_IAU_Data%spfh2m_inc)) deallocate (LND_IAU_Data%spfh2m_inc)

   if (allocated(iau_state%inc1%stc_inc)) deallocate(iau_state%inc1%stc_inc)
   if (allocated(iau_state%inc1%slc_inc)) deallocate(iau_state%inc1%slc_inc)
   if (allocated(iau_state%inc1%tmp2m_inc)) deallocate(iau_state%inc1%tmp2m_inc)
   if (allocated(iau_state%inc1%spfh2m_inc)) deallocate(iau_state%inc1%spfh2m_inc)

   if (allocated(iau_state%inc2%stc_inc)) deallocate(iau_state%inc2%stc_inc)
   if (allocated(iau_state%inc2%slc_inc)) deallocate(iau_state%inc2%slc_inc)
   if (allocated(iau_state%inc2%tmp2m_inc)) deallocate(iau_state%inc2%tmp2m_inc)
   if (allocated(iau_state%inc2%spfh2m_inc)) deallocate(iau_state%inc2%spfh2m_inc)

end subroutine lnd_iau_mod_finalize

 subroutine lnd_iau_mod_getiauforcing(LND_IAU_Control, LND_IAU_Data, errmsg, errflg)

   implicit none
   type (LND_IAU_Control_type),          intent(in) :: LND_IAU_Control
   type(lnd_iau_external_data_type),  intent(inout) :: LND_IAU_Data
   character(len=*),              intent(out) :: errmsg
   integer,                       intent(out) :: errflg
   real(kind=kind_phys) t1,t2,sx,wx,wt,dtp
   integer n,i,j,k,sphum,kstep,nstep,itnext

   LND_IAU_Data%in_interval=.false.
   if (nfiles.LE.0) then
       return
   endif

   if (nfiles .eq. 1) then 
       t1 = LND_IAU_Control%iaufhrs(1)-0.5*LND_IAU_Control%iau_delthrs
       t2 = LND_IAU_Control%iaufhrs(1)+0.5*LND_IAU_Control%iau_delthrs
   else
       t1 = LND_IAU_Control%iaufhrs(1)
       t2 = LND_IAU_Control%iaufhrs(nfiles)
   endif
   if (LND_IAU_Control%iau_filter_increments) then
      ! compute increment filter weight
      ! t1 is beginning of window, t2 end of window
      ! LND_IAU_Control%fhour current time
      ! in window kstep=-nstep,nstep (2*nstep+1 total)
      ! time step LND_IAU_Control%dtp
      dtp=LND_IAU_Control%dtp
      nstep = 0.5*LND_IAU_Control%iau_delthrs*3600/dtp
      ! compute normalized filter weight
      kstep = ((LND_IAU_Control%fhour-t1) - 0.5*LND_IAU_Control%iau_delthrs)*3600./dtp
      if (LND_IAU_Control%fhour >= t1 .and. LND_IAU_Control%fhour < t2) then
         sx     = acos(-1.)*kstep/nstep
         wx     = acos(-1.)*kstep/(nstep+1)
         if (kstep .ne. 0) then
            wt = (sin(wx)/wx*sin(sx)/sx)
         else
            wt = 1.
         endif
         iau_state%wt = iau_state%wt_normfact*wt
         !if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'kstep,t1,t,t2,filter wt=',kstep,t1,LND_IAU_Control%fhour,t2,iau_state%wt/iau_state%wt_normfact
      else
         iau_state%wt = 0.
      endif
   endif

   if (nfiles.EQ.1) then
!  on check to see if we are in the IAU window,  no need to update the
!  tendencies since they are fixed over the window
      if ( LND_IAU_Control%fhour < t1 .or. LND_IAU_Control%fhour >= t2 ) then
!         if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'no iau forcing',t1,LND_IAU_Control%fhour,t2
         LND_IAU_Data%in_interval=.false.
      else
         if (LND_IAU_Control%iau_filter_increments) call setiauforcing(LND_IAU_Control,LND_IAU_Data,iau_state%wt)
         if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'apply lnd iau forcing t1,t,t2,filter wt= ',t1,LND_IAU_Control%fhour,t2,iau_state%wt/iau_state%wt_normfact
         LND_IAU_Data%in_interval=.true.
      endif
      return
   endif

   if (nfiles > 1) then
      itnext=2
      if (LND_IAU_Control%fhour < t1 .or. LND_IAU_Control%fhour >= t2) then
!         if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'no iau forcing',LND_IAU_Control%iaufhrs(1),LND_IAU_Control%fhour,LND_IAU_Control%iaufhrs(nfiles)
         LND_IAU_Data%in_interval=.false.
      else
         if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'apply lnd iau forcing t1,t,t2,filter wt= ',t1,LND_IAU_Control%fhour,t2,iau_state%wt/iau_state%wt_normfact
         LND_IAU_Data%in_interval=.true.
         do k=nfiles, 1, -1
            if (LND_IAU_Control%iaufhrs(k) > LND_IAU_Control%fhour) then
               itnext=k
            endif
         enddo
!         if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'itnext=',itnext
         if (LND_IAU_Control%fhour >= iau_state%hr2) then ! need to read in next increment file
            iau_state%hr1=iau_state%hr2
            iau_state%hr2=LND_IAU_Control%iaufhrs(itnext)
            iau_state%inc1=iau_state%inc2
     
            ! if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'reading next lnd iau increment file',trim(LND_IAU_Control%iau_inc_files(itnext))
            ! call read_iau_forcing(LND_IAU_Control,iau_state%inc2,'INPUT/'//trim(LND_IAU_Control%iau_inc_files(itnext)), errmsg, errflg)            
            if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'interpolating next lnd iau increment ', itnext  !trim(LND_IAU_Control%iau_inc_files(itnext))
            call interp_inc_at_timestep(LND_IAU_Control, km, wk3_stc(itnext, :, :, :), iau_state%inc2%stc_inc, errmsg, errflg)
            call interp_inc_at_timestep(LND_IAU_Control, km, wk3_slc(itnext, :, :, :), iau_state%inc2%slc_inc, errmsg, errflg)
            call interp_inc_at_timestep(LND_IAU_Control, 1, wk3_t2m(itnext, :, :, :), iau_state%inc2%tmp2m_inc, errmsg, errflg)
            call interp_inc_at_timestep(LND_IAU_Control, 1, wk3_q2m(itnext, :, :, :), iau_state%inc2%spfh2m_inc, errmsg, errflg)
         endif
         call updateiauforcing(LND_IAU_Control,LND_IAU_Data,iau_state%wt)
      endif
   endif
   ! sphum=get_tracer_index(MODEL_ATMOS,'sphum')

 end subroutine lnd_iau_mod_getiauforcing

subroutine updateiauforcing(LND_IAU_Control, LND_IAU_Data, wt)

   implicit none
   type (LND_IAU_Control_type),        intent(in) :: LND_IAU_Control
   type(lnd_iau_external_data_type),  intent(inout) :: LND_IAU_Data
   real(kind_phys) delt, wt
   integer i,j,k,l

!   if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'in updateiauforcing',nfiles,LND_IAU_Control%iaufhrs(1:nfiles)
   delt = (iau_state%hr2-(LND_IAU_Control%fhour))/(IAU_state%hr2-IAU_state%hr1)
   do j = js,je
      do i = is,ie
         do k = 1,npz
         ! do k = 1,n_soill    !         
            LND_IAU_Data%stc_inc(i,j,k)  =(delt*IAU_state%inc1%stc_inc(i,j,k)  + (1.-delt)* IAU_state%inc2%stc_inc(i,j,k))*rdt*wt
            LND_IAU_Data%slc_inc(i,j,k)  =(delt*IAU_state%inc1%slc_inc(i,j,k)  + (1.-delt)* IAU_state%inc2%slc_inc(i,j,k))*rdt*wt
         end do
         LND_IAU_Data%tmp2m_inc(i,j,1)  =(delt*IAU_state%inc1%tmp2m_inc(i,j,1)  + (1.-delt)* IAU_state%inc2%tmp2m_inc(i,j,1))*rdt*wt
         LND_IAU_Data%spfh2m_inc(i,j,1)  =(delt*IAU_state%inc1%spfh2m_inc(i,j,1)  + (1.-delt)* IAU_state%inc2%spfh2m_inc(i,j,1))*rdt*wt
       enddo
   enddo
 end subroutine updateiauforcing


 subroutine setiauforcing(LND_IAU_Control, LND_IAU_Data, wt)

   implicit none
   type (LND_IAU_Control_type),        intent(in)   :: LND_IAU_Control
   type(lnd_iau_external_data_type),  intent(inout) :: LND_IAU_Data
   real(kind_phys) delt, dt,wt
   integer i,j,k,l,sphum
   !  this is only called if using 1 increment file
   if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'in lnd_iau setiauforcing rdt = ',rdt
   do j = js,je
      do i = is,ie
         do k = 1,npz
         !  do k = 1,n_soill    !         
            LND_IAU_Data%stc_inc(i,j,k) = wt*IAU_state%inc1%stc_inc(i,j,k)*rdt
            LND_IAU_Data%slc_inc(i,j,k) = wt*IAU_state%inc1%slc_inc(i,j,k)*rdt
         end do
         LND_IAU_Data%tmp2m_inc(i,j,1) = wt*IAU_state%inc1%tmp2m_inc(i,j,1)*rdt
         LND_IAU_Data%spfh2m_inc(i,j,1) = wt*IAU_state%inc1%spfh2m_inc(i,j,1)*rdt
      enddo
   enddo
   !  sphum=get_tracer_index(MODEL_ATMOS,'sphum')

 end subroutine setiauforcing

subroutine read_iau_forcing_all_timesteps(LND_IAU_Control, fname, errmsg, errflg, &
                                          wk3_out_stc, wk3_out_slc, wk3_out_t2m, wk3_out_q2m)   !, fname_sfc) is, ie, js, je, ks,ke, 
   type (LND_IAU_Control_type),   intent(in) :: LND_IAU_Control
   character(len=*),              intent(in) :: fname
   character(len=*),             intent(out) :: errmsg
   integer,                      intent(out) :: errflg
   ! integer,                       intent(in) :: is, ie, js, je, ks,ke
   ! real(kind=4),                 intent(out) :: wk3_out(is:ie,js:je,ks:ke)
   real(kind=4),                 intent(out) :: wk3_out_stc(1:im, jbeg:jend, 1:km)
   real(kind=4),                 intent(out) :: wk3_out_slc(1:im, jbeg:jend, 1:km)
   real(kind=4),                 intent(out) :: wk3_out_t2m(1:im, jbeg:jend, 1:1)
   real(kind=4),                 intent(out) :: wk3_out_q2m(1:im, jbeg:jend, 1:1)
   
   integer:: i, j, k, l, npz
   integer:: i1, i2, j1
   logical  :: exists
   integer  :: ncid

   character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
   character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']
   character(len=32), :: t2m_vars =  'tmp2m_inc'
   character(len=32), :: q2m_vars =  'spfh2m_inc'

   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0
   
   inquire (file=trim(fname), exist=exists)    
   if (exists) then
!  if( file_exist(fname) ) then
      call open_ncfile( fname, ncid )        ! open the file
   else
      ! call mpp_error(FATAL,'==> Error in read_iau_forcing: Expected file '&
      !     //trim(fname)//' for DA increment does not exist')
      errmsg = 'FATAL Error in read_iau_forcing: Expected file '//trim(fname)//' for DA increment does not exist'
      errflg = 1
      return
   endif

   do i = 1, size(stc_vars)
      print *, trim(stc_vars(i))
      call check_var_exists(ncid, trim(stc_vars(i)), ierr)
      if (ierr == 0) then
         ! call get_var3_r4( ncid, field_name, 1,im, jbeg,jend, 1,km, wk3 )
         call get_var3_r4( ncid, trim(stc_vars(i)), 1,im, jbeg,jend, 1,1, wk3_out_stc(:, :, i) )
      else
         if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'warning: no increment for ',trim(stc_vars(i)),' found, assuming zero'
         wk3_out = 0.
      endif
   enddo
   do i = 1, size(slc_vars)
      print *, trim(slc_vars(i))
      call check_var_exists(ncid, trim(slc_vars(i)), ierr)
      if (ierr == 0) then
         ! call get_var3_r4( ncid, field_name, 1,im, jbeg,jend, 1,km, wk3 )
         call get_var3_r4( ncid, trim(slc_vars(i)), 1,im, jbeg,jend, 1,1, wk3_out_slc(:, :, i) )
      else
         if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'warning: no increment for ',trim(slc_vars(i)),' found, assuming zero'
         wk3_out = 0.
      endif
   enddo
   print *, trim(t2m_vars)
   call check_var_exists(ncid, trim(t2m_vars), ierr)
   if (ierr == 0) then
      ! call get_var3_r4( ncid, field_name, 1,im, jbeg,jend, 1,km, wk3 )
      call get_var3_r4( ncid, trim(t2m_vars), 1,im, jbeg,jend, 1,1, wk3_out_t2m(:, :, :) )
   else
      if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'warning: no increment for ',trim(t2m_vars),' found, assuming zero'
      wk3_out = 0.
   endif
   print *, trim(q2m_vars)
   call check_var_exists(ncid, trim(q2m_vars), ierr)
   if (ierr == 0) then
      ! call get_var3_r4( ncid, field_name, 1,im, jbeg,jend, 1,km, wk3 )
      call get_var3_r4( ncid, trim(q2m_vars), 1,im, jbeg,jend, 1,1, wk3_out_q2m(:, :, :) )
   else
      if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'warning: no increment for ',trim(q2m_vars),' found, assuming zero'
      wk3_out = 0.
   endif

   call close_ncfile(ncid)
   
end subroutine read_iau_forcing_all_timesteps

subroutine interp_inc_at_timestep(LND_IAU_Control, km_in, wk3_in, var, errmsg, errflg)   !field_name, , jbeg, jend)
   ! interpolate increment from GSI gaussian grid to cubed sphere
   ! everying is on the A-grid, earth relative
   type (LND_IAU_Control_type), intent(in) :: LND_IAU_Control
   ! character(len=*), intent(in) :: field_name
   integer,                                intent(in) :: km_in        !jbeg,jend
   real(kind=4),                           intent(in) :: wk3_in(1:im,jbeg:jend, 1:km_in) 
   real, dimension(is:ie, js:je, 1:km), intent(inout) :: var
   
   character(len=*),              intent(out) :: errmsg
   integer,                       intent(out) :: errflg
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

subroutine read_iau_forcing(LND_IAU_Control, increments, fname, errmsg, errflg)   !, fname_sfc)
   type (LND_IAU_Control_type), intent(in) :: LND_IAU_Control
   type(iau_internal_data_type), intent(inout):: increments
   character(len=*),  intent(in) :: fname
!  character(len=*),  intent(in), optional :: fname_sfc
   character(len=*),              intent(out) :: errmsg
   integer,                       intent(out) :: errflg
!locals
!  real, dimension(:,:,:), allocatable:: u_inc, v_inc

   integer:: i, j, k, l, npz
   integer:: i1, i2, j1
   integer:: jbeg, jend
!  real(kind=R_GRID), dimension(2):: p1, p2, p3
!  real(kind=R_GRID), dimension(3):: e1, e2, ex, ey

!  logical  :: found
   integer  :: is,  ie,  js,  je, km_store
   logical  :: exists

   !Errors messages handled through CCPP error handling variables
   errmsg = ''
   errflg = 0

   is  = LND_IAU_Control%isc
   ie  = is + LND_IAU_Control%nx-1
   js  = LND_IAU_Control%jsc
   je  = js + LND_IAU_Control%ny-1

   deg2rad = pi/180.

   npz = LND_IAU_Control%lsoil
   
   inquire (file=trim(fname), exist=exists)    
   if (exists) then
!  if( file_exist(fname) ) then
   call open_ncfile( fname, ncid )        ! open the file
   else
   ! call mpp_error(FATAL,'==> Error in read_iau_forcing: Expected file '&
   !     //trim(fname)//' for DA increment does not exist')
   errmsg = 'FATAL Error in read_iau_forcing: Expected file '//trim(fname)//' for DA increment does not exist'
   errflg = 1
   return
   endif

   ! Find bounding latitudes:
   jbeg = jm-1;         jend = 2
   do j=js,je
   do i=is,ie
         j1 = jdc(i,j)
      jbeg = min(jbeg, j1)
      jend = max(jend, j1+1)
   enddo
   enddo

   km_store = km
   km = 1     ! n_soill Currently each soil layer increment is saved separately
   allocate ( wk3(1:im,jbeg:jend, 1:km) )
   ! call interp_inc('stc_inc',increments%stc_inc(:,:,:),jbeg,jend)    !TODO check var name
   call interp_inc(LND_IAU_Control, 'soilt1_inc',increments%stc_inc(:,:,1),jbeg,jend)
   call interp_inc(LND_IAU_Control, 'soilt2_inc',increments%stc_inc(:,:,2),jbeg,jend)
   call interp_inc(LND_IAU_Control, 'soilt3_inc',increments%stc_inc(:,:,3),jbeg,jend)
   call interp_inc(LND_IAU_Control, 'soilt4_inc',increments%stc_inc(:,:,4),jbeg,jend)

   call interp_inc(LND_IAU_Control, 'slc1_inc',increments%slc_inc(:,:,1),jbeg,jend)
   call interp_inc(LND_IAU_Control, 'slc2_inc',increments%slc_inc(:,:,2),jbeg,jend)
   call interp_inc(LND_IAU_Control, 'slc3_inc',increments%slc_inc(:,:,3),jbeg,jend)
   call interp_inc(LND_IAU_Control, 'slc4_inc',increments%slc_inc(:,:,4),jbeg,jend)

   call interp_inc(LND_IAU_Control, 'tmp2m_inc',increments%tmp2m_inc(:,:,1),jbeg,jend)
   call interp_inc(LND_IAU_Control, 'spfh2m_inc',increments%spfh2m_inc(:,:,1),jbeg,jend)
!  call interp_inc_sfc('stc_inc',increments%stc_inc(:,:,:),jbeg,jend, n_soill)    
   call close_ncfile(ncid)
   deallocate (wk3)
   km = km_store

end subroutine read_iau_forcing

subroutine interp_inc(LND_IAU_Control, field_name, var, jbeg, jend)
! interpolate increment from GSI gaussian grid to cubed sphere
! everying is on the A-grid, earth relative
 type (LND_IAU_Control_type), intent(in) :: LND_IAU_Control
 character(len=*), intent(in) :: field_name
 real, dimension(is:ie,js:je,1:km), intent(inout) :: var
 integer, intent(in) :: jbeg,jend
 integer:: i1, i2, j1, k,j,i,ierr
 call check_var_exists(ncid, field_name, ierr)
 if (ierr == 0) then
    call get_var3_r4( ncid, field_name, 1,im, jbeg,jend, 1,km, wk3 )
 else
    if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'warning: no increment for ',trim(field_name),' found, assuming zero'
    wk3 = 0.
 endif
 do k=1,km
    do j=js,je
       do i=is,ie
          i1 = id1(i,j)
          i2 = id2(i,j)
          j1 = jdc(i,j)
          var(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k)+&
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
       enddo
    enddo
 enddo
end subroutine interp_inc

!> This routine is copied from 'fv_treat_da_inc.F90 by Xi.Chen <xi.chen@noaa.gov>
! copying it here, due to inability to 'include' from the original module when the land iau mod is called through CCPP frameowrk
!
!> @author Xi.Chen <xi.chen@noaa.gov>
!> @date 02/12/2016
!
!  REVISION HISTORY:
!  02/12/2016 - Initial Version
  !=============================================================================
  !>@brief The subroutine 'remap_coef' calculates the coefficients for horizonal regridding.

  subroutine remap_coef( is, ie, js, je, isd, ied, jsd, jed, &
      im, jm, lon, lat, id1, id2, jdc, s2c, agrid )

    integer, intent(in):: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(in):: im, jm
    real(kind=kind_dyn),    intent(in):: lon(im), lat(jm)
    real,    intent(out):: s2c(is:ie,js:je,4)
    integer, intent(out), dimension(is:ie,js:je):: id1, id2, jdc
    real(kind=kind_dyn),    intent(in):: agrid(isd:ied,jsd:jed,2)
    ! local:
    real :: rdlon(im)
    real :: rdlat(jm)
    real:: a1, b1
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

! subroutine interp_inc_sfc(LND_IAU_Control, field_name,var,jbeg,jend, k_lv)  !is_land_in)
! ! interpolate increment from GSI gaussian grid to cubed sphere
! ! everying is on the A-grid, earth relative
!  type (LND_IAU_Control_type), intent(in) :: LND_IAU_Control
!  character(len=*), intent(in) :: field_name
!  integer, intent(in) :: jbeg, jend, k_lv
!  real, dimension(is:ie,js:je,1:k_lv), intent(inout) :: var
! !  logical, intent(in), optional :: is_land_in
! !  logical :: is_land
!  integer:: i1, i2, j1, k,j,i,ierr 
! !  k_lv = km
! !  is_land = .false.
! !  if ( present(is_land_in) ) is_land = is_land_in
! !  if (is_land) k_lv = n_soill 
!  call check_var_exists(ncid, field_name, ierr)
!  if (ierr == 0) then
!     call get_var3_r4( ncid, field_name, 1,im, jbeg,jend, 1,k_lv, wk3 )   !k, wk3 )
!  else
!     if (LND_IAU_Control%me == LND_IAU_Control%mpi_root) print *,'warning: no increment for ',trim(field_name),' found, assuming zero'
!     wk3 = 0.
!  endif
 
!  do k=1,k_lv       !km
!     do j=js,je
!        do i=is,ie
!           i1 = id1(i,j)
!           i2 = id2(i,j)
!           j1 = jdc(i,j)
!           var(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k)+&
!                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
!        enddo
!     enddo
!  enddo
 
! end subroutine interp_inc_sfc
  
end module lnd_iau_mod



!> \file GFS_debug.F90


!!
!! This is the place to switch between different debug outputs.
!! - The default behavior for Intel (or any compiler other than GNU)
!!   is to print mininmum, maximum and 32-bit Adler checksum for arrays.
!! - The default behavior for GNU is to mininmum, maximum and
!!   mean value of arrays, because calculating the checksum leads
!!   to segmentation faults with gfortran (bug in malloc?).
!! - If none of the #define preprocessor statements is used,
!!   arrays are printed in full (this is often unpractical).
!! - All output to stdout/stderr from these routines are prefixed
!!   with 'XXX: ' so that they can be easily removed from the log files
!!   using "grep -ve 'XXX: ' ..." if needed.
!! - Only one #define statement can be active at any time
!!
!! Available options for debug output:
!!
!!   #define PRINT_SUM: print mininmum, maximum and mean value of arrays
!!
!!   #define PRINT_CHKSUM: mininmum, maximum and 32-bit Adler checksum for arrays
!!

#ifdef __GFORTRAN__
#define PRINT_SUM
#else
#define PRINT_CHKSUM
#endif

!!
!!
!!

    module print_var_chksum

      use machine, only: kind_phys

      implicit none

      private

      public chksum_int, chksum_real, print_var

      interface print_var
        module procedure print_logic_0d
        module procedure print_logic_1d
        module procedure print_logic_2d
        module procedure print_int_0d
        module procedure print_int_1d
        module procedure print_int_2d
        module procedure print_real_0d
        module procedure print_real_1d
        module procedure print_real_2d
        module procedure print_real_3d
        module procedure print_real_4d
      end interface

    contains

      subroutine print_logic_0d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          logical, intent(in) :: var

          write(0,'(2a,3i6,1x,l)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_logic_0d

      subroutine print_logic_1d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          logical, intent(in) :: var(:)

          integer :: i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,2i8)') 'XXX: ', trim(name), mpirank, omprank, blkno, size(var), count(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,2i8)') 'XXX: ', trim(name), mpirank, omprank, blkno, size(var), count(var)
#else
          do i=lbound(var,1),ubound(var,1)
              write(0,'(2a,3i6,i6,2e16.7,1x,l)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, lat_d(i), lon_d(i), var(i)
          end do
#endif

      end subroutine print_logic_1d

      subroutine print_logic_2d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          logical, intent(in) :: var(:,:)

          integer :: i, k

#ifdef PRINT_SUM
          write(0,'(2a,3i6,2i8)') 'XXX: ', trim(name), mpirank, omprank, blkno, size(var), count(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,2i8)') 'XXX: ', trim(name), mpirank, omprank, blkno, size(var), count(var)
#else
          do i=lbound(var,1),ubound(var,1)
              do k=lbound(var,2),ubound(var,2)
                  write(0,'(2a,3i6,2i6,2e16.7,1x,l)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, lat_d(i), lon_d(i), var(i,k)
              end do
          end do
#endif

      end subroutine print_logic_2d

      subroutine print_int_0d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          integer, intent(in) :: var

          write(0,'(2a,3i6,i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_int_0d

      subroutine print_int_1d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          integer, intent(in) :: var(:)

          integer :: i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_int(size(var),var), minval(var), maxval(var)
#else
          do i=lbound(var,1),ubound(var,1)
              write(0,'(2a,3i6,i6,2e16.7,i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, lat_d(i), lon_d(i), var(i)
          end do
#endif

      end subroutine print_int_1d

      subroutine print_int_2d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          integer, intent(in) :: var(:,:)

          integer :: i, k

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_int(size(var),var), minval(var), maxval(var)
#else
          do i=lbound(var,1),ubound(var,1)
              do k=lbound(var,2),ubound(var,2)
                  write(0,'(2a,3i6,2i6,2e16.7,i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, lat_d(i), lon_d(i), var(i,k)
              end do
          end do
#endif

      end subroutine print_int_2d

      subroutine print_real_0d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var

          write(0,'(2a,3i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_real_0d

      subroutine print_real_1d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:)

          integer :: i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),var), minval(var), maxval(var)
#else
          do i=lbound(var,1),ubound(var,1)
              write(0,'(2a,3i6,i6,2e16.7,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, lat_d(i), lon_d(i), var(i)
          end do
#endif

      end subroutine print_real_1d

      subroutine print_real_2d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:)

          integer :: k, i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),reshape(var,(/size(var)/))), minval(var), maxval(var)
#else
          do i=lbound(var,1),ubound(var,1)
              do k=lbound(var,2),ubound(var,2)
                  write(0,'(2a,3i6,2i6,2e16.7,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, lat_d(i), lon_d(i), var(i,k)
              end do
          end do
#endif

      end subroutine print_real_2d

      subroutine print_real_3d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:,:)

          integer :: k, i, l

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),reshape(var,(/size(var)/))), minval(var), maxval(var)
#else
          do i=lbound(var,1),ubound(var,1)
              do k=lbound(var,2),ubound(var,2)
                  do l=lbound(var,3),ubound(var,3)
                      write(0,'(2a,3i6,3i6,2e16.7,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, l, lat_d(i), lon_d(i), var(i,k,l)
                  end do
              end do
          end do
#endif

      end subroutine print_real_3d

      subroutine print_real_4d(mpirank, omprank, blkno, lat_d, lon_d, name, var)

          integer, intent(in) :: mpirank, omprank, blkno
          real(kind_phys), intent(in) :: lat_d(:), lon_d(:)
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:,:,:)

          integer :: k, i, l, m

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),reshape(var,(/size(var)/))), minval(var), maxval(var)
#else
          do i=lbound(var,1),ubound(var,1)
              do k=lbound(var,2),ubound(var,2)
                  do l=lbound(var,3),ubound(var,3)
                      do m=lbound(var,4),ubound(var,4)
                          write(0,'(2a,3i6,4i6,2e16.7,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, l, m, lat_d(i), lon_d(i), var(i,k,l,m)
                      end do
                  end do
              end do
          end do
#endif

      end subroutine print_real_4d

      function chksum_int(N, var) result(hash)

          integer, intent(in) :: N
          integer, dimension(1:N), intent(in) :: var
          integer*8, dimension(1:N) :: int_var
          integer*8 :: a, b, i, hash
          integer*8, parameter :: mod_adler=65521

          a=1
          b=0
          i=1
          hash = 0
          int_var = TRANSFER(var, a, N)

          do i= 1, N
              a = MOD(a + int_var(i), mod_adler)
              b = MOD(b+a, mod_adler)
          end do

          hash = ior(b * 65536, a)

      end function chksum_int

      function chksum_real(N, var) result(hash)

          integer, intent(in) :: N
          real(kind_phys), dimension(1:N), intent(in) :: var
          integer*8, dimension(1:N) :: int_var
          integer*8 :: a, b, i, hash
          integer*8, parameter :: mod_adler=65521

          a=1
          b=0
          i=1
          hash = 0
          int_var = TRANSFER(var, a, N)

          do i= 1, N
              a = MOD(a + int_var(i), mod_adler)
              b = MOD(b+a, mod_adler)
          end do

          hash = ior(b * 65536, a)

      end function chksum_real

    end module print_var_chksum

    module GFS_diagtoscreen

      use print_var_chksum, only: print_var

      implicit none

      private

      public GFS_diagtoscreen_init, GFS_diagtoscreen_run, GFS_diagtoscreen_finalize

      contains

!> \section arg_table_GFS_diagtoscreen_init Argument Table
!! \htmlinclude GFS_diagtoscreen_init.html
!!
      subroutine GFS_diagtoscreen_init (Model, Data, Suite_Interstitial, Interstitial, errmsg, errflg)

         use GFS_typedefs,          only: GFS_control_type, GFS_data_type, &
                                          GFS_suite_interstitial_type, GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),            intent(in)  :: Model
         type(GFS_data_type),               intent(in)  :: Data(:)
         type(GFS_suite_interstitial_type), intent(in)  :: Suite_Interstitial(:)
         type(GFS_interstitial_type),       intent(in)  :: Interstitial(:)
         character(len=*),                  intent(out) :: errmsg
         integer,                           intent(out) :: errflg

         !--- local variables
         integer :: i

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         do i=1,size(Data)
           call GFS_diagtoscreen_run (Model, Data(i)%Statein, Data(i)%Stateout, Data(i)%Sfcprop,    &
                                      Data(i)%Coupling, Data(i)%Grid, Data(i)%Tbd, Data(i)%Cldprop, &
                                      Data(i)%Radtend, Data(i)%Intdiag, Suite_Interstitial(i),      &
                                      Interstitial(1), size(Interstitial), i, errmsg, errflg)
         end do

      end subroutine GFS_diagtoscreen_init

      subroutine GFS_diagtoscreen_finalize ()
      end subroutine GFS_diagtoscreen_finalize

!> \section arg_table_GFS_diagtoscreen_run Argument Table
!! \htmlinclude GFS_diagtoscreen_run.html
!!
      subroutine GFS_diagtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling,           &
                                       Grid, Tbd, Cldprop, Radtend, Diag, Suite_Interstitial, &
                                       Interstitial, nthreads, blkno, errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use GFS_typedefs,          only: GFS_control_type, GFS_statein_type,  &
                                          GFS_stateout_type, GFS_sfcprop_type, &
                                          GFS_coupling_type, GFS_grid_type,    &
                                          GFS_tbd_type, GFS_cldprop_type,      &
                                          GFS_radtend_type, GFS_diag_type,     &
                                          GFS_suite_interstitial_type,         &
                                          GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),            intent(in)  :: Model
         type(GFS_statein_type),            intent(in)  :: Statein
         type(GFS_stateout_type),           intent(in)  :: Stateout
         type(GFS_sfcprop_type),            intent(in)  :: Sfcprop
         type(GFS_coupling_type),           intent(in)  :: Coupling
         type(GFS_grid_type),               intent(in)  :: Grid
         type(GFS_tbd_type),                intent(in)  :: Tbd
         type(GFS_cldprop_type),            intent(in)  :: Cldprop
         type(GFS_radtend_type),            intent(in)  :: Radtend
         type(GFS_diag_type),               intent(in)  :: Diag
         type(GFS_suite_interstitial_type), intent(in)  :: Suite_Interstitial
         type(GFS_interstitial_type),       intent(in)  :: Interstitial
         integer,                           intent(in)  :: nthreads
         integer,                           intent(in)  :: blkno
         character(len=*),                  intent(out) :: errmsg
         integer,                           intent(out) :: errflg

         !--- local variables
         integer :: impi, iomp, ierr, n
         integer :: mpirank, mpisize, mpicomm
         integer :: omprank, ompsize

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef MPI
         mpicomm = Model%communicator
         mpirank = Model%me
         mpisize = Model%ntasks
#else
         mpirank = 0
         mpisize = 1
         mpicomm = 0
#endif
#ifdef OPENMP
         omprank = OMP_GET_THREAD_NUM()
         ompsize = nthreads
#else
         omprank = 0
         ompsize = 1
#endif

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
!         call MPI_BARRIER(mpicomm,ierr)
#endif

         do impi=0,mpisize-1
             do iomp=0,ompsize-1
                 if (mpirank==impi .and. omprank==iomp) then
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Model%kdt'        , Model%kdt)
                     ! Sfcprop
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%slmsk'    , Sfcprop%slmsk)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%oceanfrac', Sfcprop%oceanfrac)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%landfrac' , Sfcprop%landfrac)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%lakefrac' , Sfcprop%lakefrac)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tsfc'     , Sfcprop%tsfc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tsfco'    , Sfcprop%tsfco)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tsfcl'    , Sfcprop%tsfcl)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tisfc'    , Sfcprop%tisfc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%snowd'    , Sfcprop%snowd)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zorl'     , Sfcprop%zorl)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zorlw'    , Sfcprop%zorlw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zorll'    , Sfcprop%zorll)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zorli'    , Sfcprop%zorli)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zorlwav'  , Sfcprop%zorlwav)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%fice'     , Sfcprop%fice)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%hprime'   , Sfcprop%hprime)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%sncovr'   , Sfcprop%sncovr)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%snoalb'   , Sfcprop%snoalb)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%alvsf'    , Sfcprop%alvsf)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%alnsf'    , Sfcprop%alnsf)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%alvwf'    , Sfcprop%alvwf)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%alnwf'    , Sfcprop%alnwf)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%facsf'    , Sfcprop%facsf)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%facwf'    , Sfcprop%facwf)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%slope'    , Sfcprop%slope)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%shdmin'   , Sfcprop%shdmin)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%shdmax'   , Sfcprop%shdmax)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tg3'      , Sfcprop%tg3)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%vfrac'    , Sfcprop%vfrac)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%vtype'    , Sfcprop%vtype)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%stype'    , Sfcprop%stype)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%uustar'   , Sfcprop%uustar)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%oro'      , Sfcprop%oro)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%oro_uf'   , Sfcprop%oro_uf)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%hice'     , Sfcprop%hice)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%weasd'    , Sfcprop%weasd)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%canopy'   , Sfcprop%canopy)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%ffmm'     , Sfcprop%ffmm)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%ffhh'     , Sfcprop%ffhh)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%f10m'     , Sfcprop%f10m)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tprcp'    , Sfcprop%tprcp)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%srflag'   , Sfcprop%srflag)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%slc'      , Sfcprop%slc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%smc'      , Sfcprop%smc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%stc'      , Sfcprop%stc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%t2m'      , Sfcprop%t2m)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%q2m'      , Sfcprop%q2m)
                     if (Model%nstf_name(1)>0) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tref    ', Sfcprop%tref)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%z_c     ', Sfcprop%z_c)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%c_0     ', Sfcprop%c_0)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%c_d     ', Sfcprop%c_d)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%w_0     ', Sfcprop%w_0)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%w_d     ', Sfcprop%w_d)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%xt      ', Sfcprop%xt)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%xs      ', Sfcprop%xs)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%xu      ', Sfcprop%xu)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%xv      ', Sfcprop%xv)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%xz      ', Sfcprop%xz)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zm      ', Sfcprop%zm)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%xtts    ', Sfcprop%xtts)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%xzts    ', Sfcprop%xzts)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%d_conv  ', Sfcprop%d_conv)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%ifd     ', Sfcprop%ifd)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%dt_cool ', Sfcprop%dt_cool)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%qrain   ', Sfcprop%qrain)
                     end if
                     ! CCPP/RUC only
                     if (Model%lsm == Model%lsm_ruc) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%sh2o',            Sfcprop%sh2o)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%smois',           Sfcprop%smois)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tslb',            Sfcprop%tslb)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%clw_surf_land',   Sfcprop%clw_surf_land)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%clw_surf_ice',    Sfcprop%clw_surf_ice)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%qwv_surf_land',   Sfcprop%qwv_surf_land)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%qwv_surf_ice',    Sfcprop%qwv_surf_ice)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%flag_frsoil',     Sfcprop%flag_frsoil)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%rhofr',           Sfcprop%rhofr)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tsnow_land',      Sfcprop%tsnow_land)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tsnow_ice',       Sfcprop%tsnow_ice)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%snowfallac_land', Sfcprop%snowfallac_land)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%snowfallac_ice',  Sfcprop%snowfallac_ice)
                     end if
                     ! Radtend
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfcfsw%upfxc', Radtend%sfcfsw(:)%upfxc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfcfsw%dnfxc', Radtend%sfcfsw(:)%dnfxc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfcfsw%upfx0', Radtend%sfcfsw(:)%upfx0)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfcfsw%dnfx0', Radtend%sfcfsw(:)%dnfx0)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfcflw%upfxc', Radtend%sfcflw(:)%upfxc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfcflw%upfx0', Radtend%sfcflw(:)%upfx0)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfcflw%dnfxc', Radtend%sfcflw(:)%dnfxc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfcflw%dnfx0', Radtend%sfcflw(:)%dnfx0)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%htrsw',        Radtend%htrsw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%htrlw',        Radtend%htrlw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%sfalb',        Radtend%sfalb)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%coszen',       Radtend%coszen)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%tsflw',        Radtend%tsflw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%semis',        Radtend%semis)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%coszdg',       Radtend%coszdg)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%swhc',         Radtend%swhc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%lwhc',         Radtend%lwhc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Radtend%lwhd',         Radtend%lwhd)
                     ! Tbd
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%icsdsw'          , Tbd%icsdsw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%icsdlw'          , Tbd%icsdlw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%ozpl'            , Tbd%ozpl)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%h2opl'           , Tbd%h2opl)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%rann'            , Tbd%rann)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%acv'             , Tbd%acv)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%acvb'            , Tbd%acvb)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%acvt'            , Tbd%acvt)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%hpbl'            , Tbd%hpbl)
                     if (Model%do_sppt) then
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%dtdtnp'        , Tbd%dtdtnp)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%dtotprcp'      , Tbd%dtotprcp)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%dcnvprcp'      , Tbd%dcnvprcp)
                     end if
                     if (Model%cplflx .or. Model%cplchm) then
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%drain_cpl'     , Tbd%drain_cpl)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%dsnow_cpl'     , Tbd%dsnow_cpl)
                     end if
                     if (Model%nctp > 0 .and. Model%cscnv) then
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%phy_fctd'      , Tbd%phy_fctd)
                     end if
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%phy_f2d'         , Tbd%phy_f2d)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%phy_f3d'         , Tbd%phy_f3d)
                     do n=1,size(Tbd%phy_f3d(1,1,:))
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%phy_f3d_n'   , Tbd%phy_f3d(:,:,n))
                     end do
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%in_nm'           , Tbd%in_nm)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%ccn_nm'          , Tbd%ccn_nm)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%aer_nm'          , Tbd%aer_nm)
                     ! Diag
                     !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%fluxr       ',    Diag%fluxr)
                     !do n=1,size(Diag%fluxr(1,:))
                     !    call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%fluxr_n ',    Diag%fluxr(:,n))
                     !end do
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%srunoff     ',    Diag%srunoff)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%evbsa       ',    Diag%evbsa)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%evcwa       ',    Diag%evcwa)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%snohfa      ',    Diag%snohfa)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%transa      ',    Diag%transa)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%sbsnoa      ',    Diag%sbsnoa)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%snowca      ',    Diag%snowca)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%soilm       ',    Diag%soilm)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%tmpmin      ',    Diag%tmpmin)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%tmpmax      ',    Diag%tmpmax)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dusfc       ',    Diag%dusfc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dvsfc       ',    Diag%dvsfc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dtsfc       ',    Diag%dtsfc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dqsfc       ',    Diag%dqsfc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%totprcp     ',    Diag%totprcp)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%totice      ',    Diag%totice)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%totsnw      ',    Diag%totsnw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%totgrp      ',    Diag%totgrp)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%totprcpb    ',    Diag%totprcpb)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%toticeb     ',    Diag%toticeb)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%totsnwb     ',    Diag%totsnwb)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%totgrpb     ',    Diag%totgrpb)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%suntim      ',    Diag%suntim)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%runoff      ',    Diag%runoff)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%ep          ',    Diag%ep)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%cldwrk      ',    Diag%cldwrk)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dugwd       ',    Diag%dugwd)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dvgwd       ',    Diag%dvgwd)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%psmean      ',    Diag%psmean)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%cnvprcp     ',    Diag%cnvprcp)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%cnvprcpb    ',    Diag%cnvprcpb)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%spfhmin     ',    Diag%spfhmin)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%spfhmax     ',    Diag%spfhmax)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%u10mmax     ',    Diag%u10mmax)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%v10mmax     ',    Diag%v10mmax)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%wind10mmax  ',    Diag%wind10mmax)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%rain        ',    Diag%rain)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%rainc       ',    Diag%rainc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%ice         ',    Diag%ice)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%snow        ',    Diag%snow)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%graupel     ',    Diag%graupel)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%u10m        ',    Diag%u10m)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%v10m        ',    Diag%v10m)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dpt2m       ',    Diag%dpt2m)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%zlvl        ',    Diag%zlvl)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%psurf       ',    Diag%psurf)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%pwat        ',    Diag%pwat)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%t1          ',    Diag%t1)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%q1          ',    Diag%q1)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%u1          ',    Diag%u1)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%v1          ',    Diag%v1)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%chh         ',    Diag%chh)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%cmm         ',    Diag%cmm)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%epi         ',    Diag%epi)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%smcwlt2     ',    Diag%smcwlt2)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%smcref2     ',    Diag%smcref2)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%sr          ',    Diag%sr)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%tdomr       ',    Diag%tdomr)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%tdomzr      ',    Diag%tdomzr)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%tdomip      ',    Diag%tdomip)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%tdoms       ',    Diag%tdoms)
                     ! CCPP/RUC only
                     if (Model%lsm == Model%lsm_ruc) then
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%wet1        ',  Sfcprop%wetness)
                     else
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%wet1        ',  Diag%wet1)
                     end if
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%skebu_wts   ',    Diag%skebu_wts)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%skebv_wts   ',    Diag%skebv_wts)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%sppt_wts    ',    Diag%sppt_wts)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%shum_wts    ',    Diag%shum_wts)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%zmtnblck    ',    Diag%zmtnblck)
                     if (Model%ldiag3d) then
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%du3dt       ',    Diag%du3dt)
                       do n=1,size(Diag%du3dt(1,1,:))
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%du3dt_n     ',  Diag%du3dt(:,:,n))
                       end do
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dv3dt       ',    Diag%dv3dt)
                       do n=1,size(Diag%dv3dt(1,1,:))
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dv3dt_n     ',  Diag%dv3dt(:,:,n))
                       end do
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dt3dt       ',    Diag%dt3dt)
                       do n=1,size(Diag%dt3dt(1,1,:))
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dt3dt_n     ',  Diag%dt3dt(:,:,n))
                       end do
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dq3dt       ',    Diag%dq3dt)
                       do n=1,size(Diag%dq3dt(1,1,:))
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dq3dt_n     ',  Diag%dq3dt(:,:,n))
                       end do
                       !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%upd_mf      ',    Diag%upd_mf)
                       !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dwn_mf      ',    Diag%dwn_mf)
                       !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%det_mf      ',    Diag%det_mf)
                     end if
                     if(Model%lradar) then
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%refl_10cm   ',  Diag%refl_10cm)
                     end if
                     ! CCPP/MYNNPBL only
                     if (Model%do_mynnedmf) then
                       if (Model%bl_mynn_output .ne. 0) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%edmf_a      ',  Diag%edmf_a)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%edmf_w      ',  Diag%edmf_w)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%edmf_qt     ',  Diag%edmf_qt)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%edmf_thl    ',  Diag%edmf_thl)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%edmf_ent    ',  Diag%edmf_ent)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%edmf_qc     ',  Diag%edmf_qc)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%sub_thl     ',  Diag%sub_thl)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%sub_sqv     ',  Diag%sub_sqv)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%det_thl     ',  Diag%det_thl)
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%det_sqv     ',  Diag%det_sqv)
                       end if
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%nupdraft    ',  Diag%nupdraft)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%maxMF       ',  Diag%maxMF)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%ktop_plume  ',  Diag%ktop_plume)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%exch_h      ',  Diag%exch_h)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%exch_m      ',  Diag%exch_m)
                     end if
                     ! UGWP - incomplete list
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dudt_gw       ',  Diag%dudt_gw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dvdt_gw       ',  Diag%dvdt_gw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dtdt_gw       ',  Diag%dtdt_gw)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%kdis_gw       ',  Diag%kdis_gw)
                     if (Model%do_ugwp_v1 .or. Model%gwd_opt==33 .or. Model%gwd_opt==22) then
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dudt_ogw      ',  Diag%dudt_ogw )
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dvdt_ogw      ',  Diag%dvdt_ogw )
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dudt_obl      ',  Diag%dudt_obl )
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dvdt_obl      ',  Diag%dvdt_obl )
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dudt_oss      ',  Diag%dudt_oss )
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dvdt_oss      ',  Diag%dvdt_oss )
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dudt_ofd      ',  Diag%dudt_ofd )
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dvdt_ofd      ',  Diag%dvdt_ofd )
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%du_ogwcol     ',  Diag%du_ogwcol)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dv_ogwcol     ',  Diag%dv_ogwcol)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%du_oblcol     ',  Diag%du_oblcol)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dv_oblcol     ',  Diag%dv_oblcol)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%du_osscol     ',  Diag%du_osscol)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dv_osscol     ',  Diag%dv_osscol)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%du_ofdcol     ',  Diag%du_ofdcol)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dv_ofdcol     ',  Diag%dv_ofdcol)
                     else
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dudt_ogw      ',  Diag%dudt_ogw)
                     end if
                     ! Statein
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%phii'    ,     Statein%phii)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%prsi'    ,     Statein%prsi)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%prsik'   ,     Statein%prsik)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%phil'    ,     Statein%phil)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%prsl'    ,     Statein%prsl)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%prslk'   ,     Statein%prslk)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%pgr'     ,     Statein%pgr)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%ugrs'    ,     Statein%ugrs)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%vgrs'    ,     Statein%vgrs)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%vvl'     ,     Statein%vvl)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%tgrs'    ,     Statein%tgrs)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%qgrs'    ,     Statein%qgrs)
                     do n=1,size(Statein%qgrs(1,1,:))
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%qgrs_n',    Statein%qgrs(:,:,n))
                     end do
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%diss_est',     Statein%diss_est)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%smc'     ,     Statein%smc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%stc'     ,     Statein%stc)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Statein%slc'     ,     Statein%slc)
                     ! Stateout
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Stateout%gu0',         Stateout%gu0)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Stateout%gv0',         Stateout%gv0)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Stateout%gt0',         Stateout%gt0)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Stateout%gq0',         Stateout%gq0)
                     do n=1,size(Stateout%gq0(1,1,:))
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Stateout%gq0_n',    Stateout%gq0(:,:,n))
                     end do
                     ! Coupling
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nirbmdi', Coupling%nirbmdi)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nirdfdi', Coupling%nirdfdi)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%visbmdi', Coupling%visbmdi)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%visdfdi', Coupling%visdfdi)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nirbmui', Coupling%nirbmui)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nirdfui', Coupling%nirdfui)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%visbmui', Coupling%visbmui)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%visdfui', Coupling%visdfui)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%sfcdsw ', Coupling%sfcdsw )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%sfcnsw ', Coupling%sfcnsw )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%sfcdlw ', Coupling%sfcdlw )
                     if (Model%cplflx .or. Model%do_sppt .or. Model%cplchm) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%rain_cpl', Coupling%rain_cpl)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%snow_cpl', Coupling%snow_cpl)
                     end if
!                    if (Model%cplwav2atm) then
!                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%zorlwav_cpl' , Coupling%zorlwav_cpl  )
!                    end if
                     if (Model%cplflx) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%oro_cpl'     , Coupling%oro_cpl      )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%slmsk_cpl'   , Coupling%slmsk_cpl    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%slimskin_cpl', Coupling%slimskin_cpl )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dusfcin_cpl ', Coupling%dusfcin_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dvsfcin_cpl ', Coupling%dvsfcin_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dtsfcin_cpl ', Coupling%dtsfcin_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dqsfcin_cpl ', Coupling%dqsfcin_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ulwsfcin_cpl', Coupling%ulwsfcin_cpl )
!                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%tseain_cpl  ', Coupling%tseain_cpl   )
!                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%tisfcin_cpl ', Coupling%tisfcin_cpl  )
!                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ficein_cpl  ', Coupling%ficein_cpl   )
!                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%hicein_cpl  ', Coupling%hicein_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%hsnoin_cpl  ', Coupling%hsnoin_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dusfc_cpl   ', Coupling%dusfc_cpl    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dvsfc_cpl   ', Coupling%dvsfc_cpl    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dtsfc_cpl   ', Coupling%dtsfc_cpl    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dqsfc_cpl   ', Coupling%dqsfc_cpl    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dlwsfc_cpl  ', Coupling%dlwsfc_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dswsfc_cpl  ', Coupling%dswsfc_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dnirbm_cpl  ', Coupling%dnirbm_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dnirdf_cpl  ', Coupling%dnirdf_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dvisbm_cpl  ', Coupling%dvisbm_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dvisdf_cpl  ', Coupling%dvisdf_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nlwsfc_cpl  ', Coupling%nlwsfc_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nswsfc_cpl  ', Coupling%nswsfc_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nnirbm_cpl  ', Coupling%nnirbm_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nnirdf_cpl  ', Coupling%nnirdf_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nvisbm_cpl  ', Coupling%nvisbm_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nvisdf_cpl  ', Coupling%nvisdf_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dusfci_cpl  ', Coupling%dusfci_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dvsfci_cpl  ', Coupling%dvsfci_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dtsfci_cpl  ', Coupling%dtsfci_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dqsfci_cpl  ', Coupling%dqsfci_cpl   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dlwsfci_cpl ', Coupling%dlwsfci_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dswsfci_cpl ', Coupling%dswsfci_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dnirbmi_cpl ', Coupling%dnirbmi_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dnirdfi_cpl ', Coupling%dnirdfi_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dvisbmi_cpl ', Coupling%dvisbmi_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dvisdfi_cpl ', Coupling%dvisdfi_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nlwsfci_cpl ', Coupling%nlwsfci_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nswsfci_cpl ', Coupling%nswsfci_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nnirbmi_cpl ', Coupling%nnirbmi_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nnirdfi_cpl ', Coupling%nnirdfi_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nvisbmi_cpl ', Coupling%nvisbmi_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nvisdfi_cpl ', Coupling%nvisdfi_cpl  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%t2mi_cpl    ', Coupling%t2mi_cpl     )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%q2mi_cpl    ', Coupling%q2mi_cpl     )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%u10mi_cpl   ', Coupling%u10mi_cpl    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%v10mi_cpl   ', Coupling%v10mi_cpl    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%tsfci_cpl   ', Coupling%tsfci_cpl    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%psurfi_cpl  ', Coupling%psurfi_cpl   )
                     end if
                     if (Model%cplchm) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%rainc_cpl', Coupling%rainc_cpl)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ushfsfci ', Coupling%ushfsfci )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dkt      ', Coupling%dkt      )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%dqdti    ', Coupling%dqdti    )
                     end if
                     if (Model%do_sppt) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%sppt_wts', Coupling%sppt_wts)
                     end if
                     if (Model%do_shum) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%shum_wts', Coupling%shum_wts)
                     end if
                     if (Model%do_skeb) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%skebu_wts', Coupling%skebu_wts )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%skebv_wts', Coupling%skebv_wts )
                     end if
                     if (Model%lndp_type /= 0) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%sfc_wts'  , Coupling%sfc_wts   )
                     end if
                     if (Model%do_ca) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ca1      ', Coupling%ca1       )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ca_deep  ', Coupling%ca_deep   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ca_turb  ', Coupling%ca_turb   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ca_shal  ', Coupling%ca_shal   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ca_rad   ', Coupling%ca_rad    )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%ca_micro ', Coupling%ca_micro  )
                     end if
                     if(Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nwfa2d', Coupling%nwfa2d)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Coupling%nifa2d', Coupling%nifa2d)
                     end if
                     ! Grid
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%xlon  ', Grid%xlon  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%xlat  ', Grid%xlat  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%xlat_d', Grid%xlat_d)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%sinlat', Grid%sinlat)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%coslat', Grid%coslat)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%area  ', Grid%area  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%dx    ', Grid%dx    )
                     if (Model%ntoz > 0) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%ddy_o3   ', Grid%ddy_o3   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%jindx1_o3', Grid%jindx1_o3)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%jindx2_o3', Grid%jindx2_o3)
                     endif
                     if (Model%h2o_phys) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%ddy_h   ', Grid%ddy_h   )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%jindx1_h', Grid%jindx1_h)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%jindx2_h', Grid%jindx2_h)
                     endif
                     if (Model%do_ugwp_v1) then
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%ddy_j1tau ', Grid%ddy_j1tau  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%ddy_j2tau ', Grid%ddy_j2tau  )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%jindx1_tau', Grid%jindx1_tau )
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Grid%jindx2_tau', Grid%jindx2_tau )
                     endif
                 end if
#ifdef OPENMP
!$OMP BARRIER
#endif
             end do
#ifdef MPI
!             call MPI_BARRIER(mpicomm,ierr)
#endif
         end do

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
!         call MPI_BARRIER(mpicomm,ierr)
#endif

      end subroutine GFS_diagtoscreen_run

    end module GFS_diagtoscreen


    module GFS_suiteinterstitialtoscreen

      use print_var_chksum, only: print_var

      implicit none

      private

      public GFS_suiteinterstitialtoscreen_init, GFS_suiteinterstitialtoscreen_run, GFS_suiteinterstitialtoscreen_finalize

      contains

      subroutine GFS_suiteinterstitialtoscreen_init (Model, Data, Suite_Interstitial, Interstitial, errmsg, errflg)

         use GFS_typedefs,          only: GFS_control_type, GFS_data_type, &
                                          GFS_suite_interstitial_type, GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),            intent(in)  :: Model
         type(GFS_data_type),               intent(in)  :: Data(:)
         type(GFS_suite_interstitial_type), intent(in)  :: Suite_Interstitial(:)
         type(GFS_interstitial_type),       intent(in)  :: Interstitial(:)
         character(len=*),                  intent(out) :: errmsg
         integer,                           intent(out) :: errflg

         !--- local variables
         integer :: i

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0


         do i=1,size(Interstitial)
           call GFS_suiteinterstitialtoscreen_run (Model, Data(1)%Statein, Data(1)%Stateout, Data(1)%Sfcprop,    &
                                              Data(1)%Coupling, Data(1)%Grid, Data(1)%Tbd, Data(1)%Cldprop, &
                                              Data(1)%Radtend, Data(1)%Intdiag, Suite_Interstitial(1),      &
                                              Interstitial(i), size(Interstitial), -999, errmsg, errflg)
         end do

      end subroutine GFS_suiteinterstitialtoscreen_init

      subroutine GFS_suiteinterstitialtoscreen_finalize ()
      end subroutine GFS_suiteinterstitialtoscreen_finalize

!> \section arg_table_GFS_suiteinterstitialtoscreen_run Argument Table
!! \htmlinclude GFS_suiteinterstitialtoscreen_run.html
!!
      subroutine GFS_suiteinterstitialtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling,       &
                                           Grid, Tbd, Cldprop, Radtend, Diag, Suite_Interstitial, &
                                           Interstitial, nthreads, blkno, errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type, GFS_statein_type,  &
                                          GFS_stateout_type, GFS_sfcprop_type, &
                                          GFS_coupling_type, GFS_grid_type,    &
                                          GFS_tbd_type, GFS_cldprop_type,      &
                                          GFS_radtend_type, GFS_diag_type,     &
                                          GFS_suite_interstitial_type,         &
                                          GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),            intent(in)  :: Model
         type(GFS_statein_type),            intent(in)  :: Statein
         type(GFS_stateout_type),           intent(in)  :: Stateout
         type(GFS_sfcprop_type),            intent(in)  :: Sfcprop
         type(GFS_coupling_type),           intent(in)  :: Coupling
         type(GFS_grid_type),               intent(in)  :: Grid
         type(GFS_tbd_type),                intent(in)  :: Tbd
         type(GFS_cldprop_type),            intent(in)  :: Cldprop
         type(GFS_radtend_type),            intent(in)  :: Radtend
         type(GFS_diag_type),               intent(in)  :: Diag
         type(GFS_suite_interstitial_type), intent(in)  :: Suite_Interstitial
         type(GFS_interstitial_type),       intent(in)  :: Interstitial
         integer,                           intent(in)  :: nthreads
         integer,                           intent(in)  :: blkno
         character(len=*),                  intent(out) :: errmsg
         integer,                           intent(out) :: errflg

         !--- local variables
         integer :: impi, iomp, ierr
         integer :: mpirank, mpisize, mpicomm
         integer :: omprank, ompsize
         integer :: istart, iend, kstart, kend

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef MPI
         mpicomm = Model%communicator
         mpirank = Model%me
         call MPI_COMM_SIZE(mpicomm, mpisize, ierr)
#else
         mpirank = 0
         mpisize = 1
         mpicomm = 0
#endif
#ifdef OPENMP
         omprank = OMP_GET_THREAD_NUM()
         ompsize = nthreads
#else
         omprank = 0
         ompsize = 1
#endif

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
!         call MPI_BARRIER(mpicomm,ierr)
#endif

         do impi=0,mpisize-1
             do iomp=0,ompsize-1
                 if (mpirank==impi .and. omprank==iomp) then
                     ! Print static variables
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ipr                 ', Suite_Interstitial%ipr                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%itc                 ', Suite_Interstitial%itc                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%latidxprnt          ', Suite_Interstitial%latidxprnt              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%levi                ', Suite_Interstitial%levi                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%lmk                 ', Suite_Interstitial%lmk                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%lmp                 ', Suite_Interstitial%lmp                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nbdlw               ', Suite_Interstitial%nbdlw                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nbdsw               ', Suite_Interstitial%nbdsw                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nf_aelw             ', Suite_Interstitial%nf_aelw                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nf_aesw             ', Suite_Interstitial%nf_aesw                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nsamftrac           ', Suite_Interstitial%nsamftrac               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nscav               ', Suite_Interstitial%nscav                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nspc1               ', Suite_Interstitial%nspc1                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ntiwx               ', Suite_Interstitial%ntiwx                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nvdiff              ', Suite_Interstitial%nvdiff                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%phys_hydrostatic    ', Suite_Interstitial%phys_hydrostatic        )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%skip_macro          ', Suite_Interstitial%skip_macro              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%trans_aero          ', Suite_Interstitial%trans_aero              )
                     ! Print all other variables
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjsfculw_land      ', Suite_Interstitial%adjsfculw_land          )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjsfculw_ice       ', Suite_Interstitial%adjsfculw_ice           )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjsfculw_water     ', Suite_Interstitial%adjsfculw_water         )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjnirbmd           ', Suite_Interstitial%adjnirbmd               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjnirbmu           ', Suite_Interstitial%adjnirbmu               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjnirdfd           ', Suite_Interstitial%adjnirdfd               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjnirdfu           ', Suite_Interstitial%adjnirdfu               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjvisbmd           ', Suite_Interstitial%adjvisbmd               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjvisbmu           ', Suite_Interstitial%adjvisbmu               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjvisdfu           ', Suite_Interstitial%adjvisdfu               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%adjvisdfd           ', Suite_Interstitial%adjvisdfd               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%aerodp              ', Suite_Interstitial%aerodp                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%alb1d               ', Suite_Interstitial%alb1d                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%bexp1d              ', Suite_Interstitial%bexp1d                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cd                  ', Suite_Interstitial%cd                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cd_ice              ', Suite_Interstitial%cd_ice                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cd_land             ', Suite_Interstitial%cd_land                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cd_water            ', Suite_Interstitial%cd_water                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cdq                 ', Suite_Interstitial%cdq                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cdq_ice             ', Suite_Interstitial%cdq_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cdq_land            ', Suite_Interstitial%cdq_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cdq_water           ', Suite_Interstitial%cdq_water               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%chh_ice             ', Suite_Interstitial%chh_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%chh_land            ', Suite_Interstitial%chh_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%chh_water           ', Suite_Interstitial%chh_water               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cldf                ', Suite_Interstitial%cldf                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cldsa               ', Suite_Interstitial%cldsa                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cld1d               ', Suite_Interstitial%cld1d                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%clw                 ', Suite_Interstitial%clw                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%clx                 ', Suite_Interstitial%clx                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%clouds              ', Suite_Interstitial%clouds                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cmm_ice             ', Suite_Interstitial%cmm_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cmm_land            ', Suite_Interstitial%cmm_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cmm_water           ', Suite_Interstitial%cmm_water               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cnvc                ', Suite_Interstitial%cnvc                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cnvw                ', Suite_Interstitial%cnvw                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ctei_r              ', Suite_Interstitial%ctei_r                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ctei_rml            ', Suite_Interstitial%ctei_rml                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cumabs              ', Suite_Interstitial%cumabs                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dd_mf               ', Suite_Interstitial%dd_mf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%de_lgth             ', Suite_Interstitial%de_lgth                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%del                 ', Suite_Interstitial%del                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%del_gz              ', Suite_Interstitial%del_gz                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%delr                ', Suite_Interstitial%delr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dkt                 ', Suite_Interstitial%dkt                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dlength             ', Suite_Interstitial%dlength                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dqdt                ', Suite_Interstitial%dqdt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dqsfc1              ', Suite_Interstitial%dqsfc1                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%drain               ', Suite_Interstitial%drain                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dtdt                ', Suite_Interstitial%dtdt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dtsfc1              ', Suite_Interstitial%dtsfc1                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dtzm                ', Suite_Interstitial%dtzm                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dt_mf               ', Suite_Interstitial%dt_mf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dudt                ', Suite_Interstitial%dudt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dusfcg              ', Suite_Interstitial%dusfcg                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dusfc1              ', Suite_Interstitial%dusfc1                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dvdftra             ', Suite_Interstitial%dvdftra                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dvdt                ', Suite_Interstitial%dvdt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dvsfcg              ', Suite_Interstitial%dvsfcg                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dvsfc1              ', Suite_Interstitial%dvsfc1                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dzlyr               ', Suite_Interstitial%dzlyr                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%elvmax              ', Suite_Interstitial%elvmax                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ep1d                ', Suite_Interstitial%ep1d                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ep1d_ice            ', Suite_Interstitial%ep1d_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ep1d_land           ', Suite_Interstitial%ep1d_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ep1d_water          ', Suite_Interstitial%ep1d_water              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%evapq               ', Suite_Interstitial%evapq                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%evap_ice            ', Suite_Interstitial%evap_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%evap_land           ', Suite_Interstitial%evap_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%evap_water          ', Suite_Interstitial%evap_water              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%evbs                ', Suite_Interstitial%evbs                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%evcw                ', Suite_Interstitial%evcw                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%faerlw              ', Suite_Interstitial%faerlw                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%faersw              ', Suite_Interstitial%faersw                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ffhh_ice            ', Suite_Interstitial%ffhh_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ffhh_land           ', Suite_Interstitial%ffhh_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ffhh_water          ', Suite_Interstitial%ffhh_water              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fh2                 ', Suite_Interstitial%fh2                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fh2_ice             ', Suite_Interstitial%fh2_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fh2_land            ', Suite_Interstitial%fh2_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fh2_water           ', Suite_Interstitial%fh2_water               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%flag_cice           ', Suite_Interstitial%flag_cice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%flag_guess          ', Suite_Interstitial%flag_guess              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%flag_iter           ', Suite_Interstitial%flag_iter               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ffmm_ice            ', Suite_Interstitial%ffmm_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ffmm_land           ', Suite_Interstitial%ffmm_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ffmm_water          ', Suite_Interstitial%ffmm_water              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fm10                ', Suite_Interstitial%fm10                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fm10_ice            ', Suite_Interstitial%fm10_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fm10_land           ', Suite_Interstitial%fm10_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fm10_water          ', Suite_Interstitial%fm10_water              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%frain               ', Suite_Interstitial%frain                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%frland              ', Suite_Interstitial%frland                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fscav               ', Suite_Interstitial%fscav                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%fswtr               ', Suite_Interstitial%fswtr                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gabsbdlw            ', Suite_Interstitial%gabsbdlw                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gabsbdlw_ice        ', Suite_Interstitial%gabsbdlw_ice            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gabsbdlw_land       ', Suite_Interstitial%gabsbdlw_land           )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gabsbdlw_water      ', Suite_Interstitial%gabsbdlw_water          )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gamma               ', Suite_Interstitial%gamma                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gamq                ', Suite_Interstitial%gamq                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gamt                ', Suite_Interstitial%gamt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gasvmr              ', Suite_Interstitial%gasvmr                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gflx                ', Suite_Interstitial%gflx                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gflx_ice            ', Suite_Interstitial%gflx_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gflx_land           ', Suite_Interstitial%gflx_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gflx_water          ', Suite_Interstitial%gflx_water              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gwdcu               ', Suite_Interstitial%gwdcu                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%gwdcv               ', Suite_Interstitial%gwdcv                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%hefac               ', Suite_Interstitial%hefac                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%hffac               ', Suite_Interstitial%hffac                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%hflxq               ', Suite_Interstitial%hflxq                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%hflx_ice            ', Suite_Interstitial%hflx_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%hflx_land           ', Suite_Interstitial%hflx_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%hflx_water          ', Suite_Interstitial%hflx_water              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%htlwc               ', Suite_Interstitial%htlwc                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%htlw0               ', Suite_Interstitial%htlw0                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%htswc               ', Suite_Interstitial%htswc                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%htsw0               ', Suite_Interstitial%htsw0                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dry                 ', Suite_Interstitial%dry                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%idxday              ', Suite_Interstitial%idxday                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%icy                 ', Suite_Interstitial%icy                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%lake                ', Suite_Interstitial%lake                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ocean               ', Suite_Interstitial%ocean                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%islmsk              ', Suite_Interstitial%islmsk                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%islmsk_cice         ', Suite_Interstitial%islmsk_cice             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%wet                 ', Suite_Interstitial%wet                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%kb                  ', Suite_Interstitial%kb                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%kbot                ', Suite_Interstitial%kbot                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%kcnv                ', Suite_Interstitial%kcnv                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%kd                  ', Suite_Interstitial%kd                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%kinver              ', Suite_Interstitial%kinver                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%kpbl                ', Suite_Interstitial%kpbl                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%kt                  ', Suite_Interstitial%kt                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ktop                ', Suite_Interstitial%ktop                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%mbota               ', Suite_Interstitial%mbota                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%mtopa               ', Suite_Interstitial%mtopa                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nday                ', Suite_Interstitial%nday                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%oa4                 ', Suite_Interstitial%oa4                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%oc                  ', Suite_Interstitial%oc                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%olyr                ', Suite_Interstitial%olyr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%plvl                ', Suite_Interstitial%plvl                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%plyr                ', Suite_Interstitial%plyr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%prcpmp              ', Suite_Interstitial%prcpmp                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%prnum               ', Suite_Interstitial%prnum                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qlyr                ', Suite_Interstitial%qlyr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qss_ice             ', Suite_Interstitial%qss_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qss_land            ', Suite_Interstitial%qss_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qss_water           ', Suite_Interstitial%qss_water               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%radar_reset         ', Suite_Interstitial%radar_reset             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%raddt               ', Suite_Interstitial%raddt                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%raincd              ', Suite_Interstitial%raincd                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%raincs              ', Suite_Interstitial%raincs                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rainmcadj           ', Suite_Interstitial%rainmcadj               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rainp               ', Suite_Interstitial%rainp                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rb                  ', Suite_Interstitial%rb                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rb_ice              ', Suite_Interstitial%rb_ice                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rb_land             ', Suite_Interstitial%rb_land                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rb_water            ', Suite_Interstitial%rb_water                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%reset               ', Suite_Interstitial%reset                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rhc                 ', Suite_Interstitial%rhc                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%runoff              ', Suite_Interstitial%runoff                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%save_q              ', Suite_Interstitial%save_q                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%save_t              ', Suite_Interstitial%save_t                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%save_tcp            ', Suite_Interstitial%save_tcp                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%save_u              ', Suite_Interstitial%save_u                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%save_v              ', Suite_Interstitial%save_v                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%sbsno               ', Suite_Interstitial%sbsno                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%uvbfc        ', Suite_Interstitial%scmpsw%uvbfc            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%uvbf0        ', Suite_Interstitial%scmpsw%uvbf0            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%nirbm        ', Suite_Interstitial%scmpsw%nirbm            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%nirdf        ', Suite_Interstitial%scmpsw%nirdf            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%visbm        ', Suite_Interstitial%scmpsw%visbm            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%visdf        ', Suite_Interstitial%scmpsw%visdf            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%semis_ice           ', Suite_Interstitial%semis_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%semis_land          ', Suite_Interstitial%semis_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%semis_water         ', Suite_Interstitial%semis_water             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%sfcalb              ', Suite_Interstitial%sfcalb                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%sigma               ', Suite_Interstitial%sigma                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%sigmaf              ', Suite_Interstitial%sigmaf                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%sigmafrac           ', Suite_Interstitial%sigmafrac               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%sigmatot            ', Suite_Interstitial%sigmatot                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%slopetype           ', Suite_Interstitial%slopetype               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snowc               ', Suite_Interstitial%snowc                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snowd_ice           ', Suite_Interstitial%snowd_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snowd_land          ', Suite_Interstitial%snowd_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snowd_water         ', Suite_Interstitial%snowd_water             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snohf               ', Suite_Interstitial%snohf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snowmt              ', Suite_Interstitial%snowmt                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%soiltype            ', Suite_Interstitial%soiltype                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%stress              ', Suite_Interstitial%stress                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%stress_ice          ', Suite_Interstitial%stress_ice              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%stress_land         ', Suite_Interstitial%stress_land             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%stress_water        ', Suite_Interstitial%stress_water            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%theta               ', Suite_Interstitial%theta                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tice                ', Suite_Interstitial%tice                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tlvl                ', Suite_Interstitial%tlvl                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tlyr                ', Suite_Interstitial%tlyr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tprcp_ice           ', Suite_Interstitial%tprcp_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tprcp_land          ', Suite_Interstitial%tprcp_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tprcp_water         ', Suite_Interstitial%tprcp_water             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%trans               ', Suite_Interstitial%trans                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tseal               ', Suite_Interstitial%tseal                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsfa                ', Suite_Interstitial%tsfa                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsfc_ice            ', Suite_Interstitial%tsfc_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsfc_land           ', Suite_Interstitial%tsfc_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsfc_water          ', Suite_Interstitial%tsfc_water              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsfg                ', Suite_Interstitial%tsfg                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsurf               ', Suite_Interstitial%tsurf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsurf_ice           ', Suite_Interstitial%tsurf_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsurf_land          ', Suite_Interstitial%tsurf_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsurf_water         ', Suite_Interstitial%tsurf_water             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ud_mf               ', Suite_Interstitial%ud_mf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%uustar_ice          ', Suite_Interstitial%uustar_ice              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%uustar_land         ', Suite_Interstitial%uustar_land             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%uustar_water        ', Suite_Interstitial%uustar_water            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%vdftra              ', Suite_Interstitial%vdftra                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%vegf1d              ', Suite_Interstitial%vegf1d                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%vegtype             ', Suite_Interstitial%vegtype                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%wcbmax              ', Suite_Interstitial%wcbmax                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%weasd_ice           ', Suite_Interstitial%weasd_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%weasd_land          ', Suite_Interstitial%weasd_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%weasd_water         ', Suite_Interstitial%weasd_water             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%wind                ', Suite_Interstitial%wind                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%work1               ', Suite_Interstitial%work1                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%work2               ', Suite_Interstitial%work2                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%work3               ', Suite_Interstitial%work3                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%xcosz               ', Suite_Interstitial%xcosz                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%xlai1d              ', Suite_Interstitial%xlai1d                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%xmu                 ', Suite_Interstitial%xmu                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%z01d                ', Suite_Interstitial%z01d                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%zt1d                ', Suite_Interstitial%zt1d                    )
                     ! UGWP
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tau_mtb             ', Suite_Interstitial%tau_mtb                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tau_ogw             ', Suite_Interstitial%tau_ogw                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tau_tofd            ', Suite_Interstitial%tau_tofd                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tau_ngw             ', Suite_Interstitial%tau_ngw                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tau_oss             ', Suite_Interstitial%tau_oss                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dudt_mtb            ', Suite_Interstitial%dudt_mtb                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dudt_tms            ', Suite_Interstitial%dudt_tms                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%zmtb                ', Suite_Interstitial%zmtb                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%zlwb                ', Suite_Interstitial%zlwb                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%zogw                ', Suite_Interstitial%zogw                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%zngw                ', Suite_Interstitial%zngw                    )
                     ! UGWP v1
                     if (Model%do_ugwp_v1) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dudt_ngw            ', Suite_Interstitial%dudt_ngw                )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dvdt_ngw            ', Suite_Interstitial%dvdt_ngw                )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dtdt_ngw            ', Suite_Interstitial%dtdt_ngw                )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%kdis_ngw            ', Suite_Interstitial%kdis_ngw                )
                     end if
                     !-- GSD drag suite
                     if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
                         Model%gwd_opt==2 .or. Model%gwd_opt==22) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%varss               ', Suite_Interstitial%varss                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ocss                ', Suite_Interstitial%ocss                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%oa4ss               ', Suite_Interstitial%oa4ss                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%clxss               ', Suite_Interstitial%clxss                   )
                     end if
                     ! GFDL and Thompson MP
                     if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%graupelmp           ', Suite_Interstitial%graupelmp               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%icemp               ', Suite_Interstitial%icemp                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rainmp              ', Suite_Interstitial%rainmp                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snowmp              ', Suite_Interstitial%snowmp                  )
                     ! Ferrier-Aligo
                     else if (Model%imp_physics == Model%imp_physics_fer_hires) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%f_ice               ', Suite_Interstitial%f_ice                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%f_rain              ', Suite_Interstitial%f_rain                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%f_rimef             ', Suite_Interstitial%f_rimef                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cwm                 ', Suite_Interstitial%cwm                     )
                     ! Morrison-Gettelman
                     else if (Model%imp_physics == Model%imp_physics_mg) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ncgl                     ', Suite_Interstitial%ncgl                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ncpr                     ', Suite_Interstitial%ncpr                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ncps                     ', Suite_Interstitial%ncps                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qgl                      ', Suite_Interstitial%qgl                              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qrn                      ', Suite_Interstitial%qrn                              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qsnw                     ', Suite_Interstitial%qsnw                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qlcn                     ', Suite_Interstitial%qlcn                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qicn                     ', Suite_Interstitial%qicn                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%w_upi                    ', Suite_Interstitial%w_upi                            )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cf_upi                   ', Suite_Interstitial%cf_upi                           )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cnv_mfd                  ', Suite_Interstitial%cnv_mfd                          )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cnv_dqldt                ', Suite_Interstitial%cnv_dqldt                        )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%clcn                     ', Suite_Interstitial%clcn                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cnv_fice                 ', Suite_Interstitial%cnv_fice                         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cnv_ndrop                ', Suite_Interstitial%cnv_ndrop                        )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cnv_nice                 ', Suite_Interstitial%cnv_nice                         )
                     end if
                     ! SHOC
                     if (Model%do_shoc) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ncgl                     ', Suite_Interstitial%ncgl                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qrn                      ', Suite_Interstitial%qrn                              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qsnw                     ', Suite_Interstitial%qsnw                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qgl                      ', Suite_Interstitial%qgl                              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ncpi                     ', Suite_Interstitial%ncpi                             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ncpl                     ', Suite_Interstitial%ncpl                             )
                     end if
                     ! Noah MP
                     if (Model%lsm == Model%lsm_noahmp) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%t2mmp                        ', Suite_Interstitial%t2mmp                        )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%q2mp                         ', Suite_Interstitial%q2mp                         )
                     end if
                     if (.not. Model%do_RRTMGP) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%alpha                        ', Suite_Interstitial%alpha                        )
                     end if
                     if (Model%lsm == Model%lsm_noah_wrfv4) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%canopy_save                  ', Suite_Interstitial%canopy_save                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%chk_land                     ', Suite_Interstitial%chk_land                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%cmc                          ', Suite_Interstitial%cmc                          )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%dqsdt2                       ', Suite_Interstitial%dqsdt2                       )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%drain_in_m_sm1               ', Suite_Interstitial%drain_in_m_sm1               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%flag_lsm                     ', Suite_Interstitial%flag_lsm                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%flag_lsm_glacier             ', Suite_Interstitial%flag_lsm_glacier             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qs1                          ', Suite_Interstitial%qs1                          )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qv1                          ', Suite_Interstitial%qv1                          )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%rho1                         ', Suite_Interstitial%rho1                         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%runoff_in_m_sm1              ', Suite_Interstitial%runoff_in_m_sm1              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%slc_save                     ', Suite_Interstitial%slc_save                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%smcmax                       ', Suite_Interstitial%smcmax                       )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%smc_save                     ', Suite_Interstitial%smc_save                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snowd_land_save              ', Suite_Interstitial%snowd_land_save              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snow_depth                   ', Suite_Interstitial%snow_depth                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snohf_snow                   ', Suite_Interstitial%snohf_snow                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snohf_frzgra                 ', Suite_Interstitial%snohf_frzgra                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%snohf_snowmelt               ', Suite_Interstitial%snohf_snowmelt               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%soilm_in_m                   ', Suite_Interstitial%soilm_in_m                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%stc_save                     ', Suite_Interstitial%stc_save                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%th1                          ', Suite_Interstitial%th1                          )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tprcp_rate_land              ', Suite_Interstitial%tprcp_rate_land              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tsfc_land_save               ', Suite_Interstitial%tsfc_land_save               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%weasd_land_save              ', Suite_Interstitial%weasd_land_save              )
                     end if
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%mg3_as_mg2                   ', Suite_Interstitial%mg3_as_mg2                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ncstrac                      ', Suite_Interstitial%ncstrac                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nn                           ', Suite_Interstitial%nn                           )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%nncl                         ', Suite_Interstitial%nncl                         )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ntk                          ', Suite_Interstitial%ntk                          )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%ntkev                        ', Suite_Interstitial%ntkev                        )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%otspt                        ', Suite_Interstitial%otspt                        )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%oz_coeffp5                   ', Suite_Interstitial%oz_coeffp5                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tracers_start_index          ', Suite_Interstitial%tracers_start_index          )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tracers_total                ', Suite_Interstitial%tracers_total                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%tracers_water                ', Suite_Interstitial%tracers_water                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%lndp_vgf                     ', Suite_Interstitial%lndp_vgf                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%uvbfc                 ', Suite_Interstitial%scmpsw%uvbfc                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%uvbf0                 ', Suite_Interstitial%scmpsw%uvbf0                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%nirbm                 ', Suite_Interstitial%scmpsw%nirbm                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%nirdf                 ', Suite_Interstitial%scmpsw%nirdf                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%visbm                 ', Suite_Interstitial%scmpsw%visbm                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%scmpsw%visdf                 ', Suite_Interstitial%scmpsw%visdf                 )
                     if (Model%imp_physics == Model%imp_physics_fer_hires) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qv_r                         ', Suite_Interstitial%qv_r                         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qc_r                         ', Suite_Interstitial%qc_r                         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qi_r                         ', Suite_Interstitial%qi_r                         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qr_r                         ', Suite_Interstitial%qr_r                         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qs_r                         ', Suite_Interstitial%qs_r                         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Suite_Interstitial%qg_r                         ', Suite_Interstitial%qg_r                         )
                     end if
                 end if
#ifdef OPENMP
!$OMP BARRIER
#endif
             end do
#ifdef MPI
!             call MPI_BARRIER(mpicomm,ierr)
#endif
         end do

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
!         call MPI_BARRIER(mpicomm,ierr)
#endif

      end subroutine GFS_suiteinterstitialtoscreen_run

    end module GFS_suiteinterstitialtoscreen

    module GFS_interstitialtoscreen

      use print_var_chksum, only: print_var

      implicit none

      private

      public GFS_interstitialtoscreen_init, GFS_interstitialtoscreen_run, GFS_interstitialtoscreen_finalize

      contains

      subroutine GFS_interstitialtoscreen_init (Model, Data, Suite_Interstitial, Interstitial, errmsg, errflg)

         use GFS_typedefs,          only: GFS_control_type, GFS_data_type, &
                                          GFS_suite_interstitial_type, GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),            intent(in)  :: Model
         type(GFS_data_type),               intent(in)  :: Data(:)
         type(GFS_suite_interstitial_type), intent(in)  :: Suite_Interstitial(:)
         type(GFS_interstitial_type),       intent(in)  :: Interstitial(:)
         character(len=*),                  intent(out) :: errmsg
         integer,                           intent(out) :: errflg

         !--- local variables
         integer :: i

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0


         do i=1,size(Interstitial)
           call GFS_interstitialtoscreen_run (Model, Data(1)%Statein, Data(1)%Stateout, Data(1)%Sfcprop,    &
                                              Data(1)%Coupling, Data(1)%Grid, Data(1)%Tbd, Data(1)%Cldprop, &
                                              Data(1)%Radtend, Data(1)%Intdiag, Suite_Interstitial(1),      &
                                              Interstitial(i), size(Interstitial), -999, errmsg, errflg)
         end do

      end subroutine GFS_interstitialtoscreen_init

      subroutine GFS_interstitialtoscreen_finalize ()
      end subroutine GFS_interstitialtoscreen_finalize

!> \section arg_table_GFS_interstitialtoscreen_run Argument Table
!! \htmlinclude GFS_interstitialtoscreen_run.html
!!
      subroutine GFS_interstitialtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling,       &
                                           Grid, Tbd, Cldprop, Radtend, Diag, Suite_Interstitial, &
                                           Interstitial, nthreads, blkno, errmsg, errflg)

#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type, GFS_statein_type,  &
                                          GFS_stateout_type, GFS_sfcprop_type, &
                                          GFS_coupling_type, GFS_grid_type,    &
                                          GFS_tbd_type, GFS_cldprop_type,      &
                                          GFS_radtend_type, GFS_diag_type,     &
                                          GFS_suite_interstitial_type,         &
                                          GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),            intent(in)  :: Model
         type(GFS_statein_type),            intent(in)  :: Statein
         type(GFS_stateout_type),           intent(in)  :: Stateout
         type(GFS_sfcprop_type),            intent(in)  :: Sfcprop
         type(GFS_coupling_type),           intent(in)  :: Coupling
         type(GFS_grid_type),               intent(in)  :: Grid
         type(GFS_tbd_type),                intent(in)  :: Tbd
         type(GFS_cldprop_type),            intent(in)  :: Cldprop
         type(GFS_radtend_type),            intent(in)  :: Radtend
         type(GFS_diag_type),               intent(in)  :: Diag
         type(GFS_suite_interstitial_type), intent(in)  :: Suite_Interstitial
         type(GFS_interstitial_type),       intent(in)  :: Interstitial
         integer,                           intent(in)  :: nthreads
         integer,                           intent(in)  :: blkno
         character(len=*),                  intent(out) :: errmsg
         integer,                           intent(out) :: errflg

         !--- local variables
         integer :: impi, iomp, ierr
         integer :: mpirank, mpisize, mpicomm
         integer :: omprank, ompsize
         integer :: istart, iend, kstart, kend

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

#ifdef MPI
         mpicomm = Model%communicator
         mpirank = Model%me
         call MPI_COMM_SIZE(mpicomm, mpisize, ierr)
#else
         mpirank = 0
         mpisize = 1
         mpicomm = 0
#endif
#ifdef OPENMP
         omprank = OMP_GET_THREAD_NUM()
         ompsize = nthreads
#else
         omprank = 0
         ompsize = 1
#endif

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
!         call MPI_BARRIER(mpicomm,ierr)
#endif

         do impi=0,mpisize-1
             do iomp=0,ompsize-1
                 if (mpirank==impi .and. omprank==iomp) then
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cldtaulw            ', Interstitial%cldtaulw                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cldtausw            ', Interstitial%cldtausw                )
                     ! RRTMGP
                     if (Model%do_RRTMGP) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%aerosolslw          ', Interstitial%aerosolslw              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%aerosolssw          ', Interstitial%aerosolssw              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_frac            ', Interstitial%cld_frac                )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_lwp             ', Interstitial%cld_lwp                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_reliq           ', Interstitial%cld_reliq               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_iwp             ', Interstitial%cld_iwp                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_reice           ', Interstitial%cld_reice               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_swp             ', Interstitial%cld_swp                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_resnow          ', Interstitial%cld_resnow              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_rwp             ', Interstitial%cld_rwp                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld_rerain          ', Interstitial%cld_rerain              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%precip_frac         ', Interstitial%precip_frac             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%icseed_lw           ', Interstitial%icseed_lw               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%icseed_sw           ', Interstitial%icseed_sw               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxlwUP_allsky     ', Interstitial%fluxlwUP_allsky         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxlwDOWN_allsky   ', Interstitial%fluxlwDOWN_allsky       )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxlwUP_clrsky     ', Interstitial%fluxlwUP_clrsky         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxlwDOWN_clrsky   ', Interstitial%fluxlwDOWN_clrsky       )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxswUP_allsky     ', Interstitial%fluxswUP_allsky         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxswDOWN_allsky   ', Interstitial%fluxswDOWN_allsky       )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxswUP_clrsky     ', Interstitial%fluxswUP_clrsky         )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxswDOWN_clrsky   ', Interstitial%fluxswDOWN_clrsky       )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%relhum              ', Interstitial%relhum                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%q_lay               ', Interstitial%q_lay                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qs_lay              ', Interstitial%qs_lay                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%deltaZ              ', Interstitial%deltaZ                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%p_lay               ', Interstitial%p_lay                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%p_lev               ', Interstitial%p_lev                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%t_lay               ', Interstitial%t_lay                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%t_lev               ', Interstitial%t_lev                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tv_lay              ', Interstitial%tv_lay                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cloud_overlap_param ', Interstitial%cloud_overlap_param     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%precip_overlap_param', Interstitial%precip_overlap_param    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%minGPpres                    ', Interstitial%minGPpres                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%minGPtemp                    ', Interstitial%minGPtemp                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ipsdlw0                      ', Interstitial%ipsdlw0                      )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ipsdsw0                      ', Interstitial%ipsdsw0                      )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tracer                       ', Interstitial%tracer                       )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxlwUP_jac                 ', Interstitial%fluxlwUP_jac                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fluxlwDOWN_jac               ', Interstitial%fluxlwDOWN_jac               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sfc_emiss_byband             ', Interstitial%sfc_emiss_byband             )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sec_diff_byband              ', Interstitial%sec_diff_byband              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sfc_alb_nir_dir              ', Interstitial%sfc_alb_nir_dir              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sfc_alb_nir_dif              ', Interstitial%sfc_alb_nir_dif              )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sfc_alb_uvvis_dir            ', Interstitial%sfc_alb_uvvis_dir            )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sfc_alb_uvvis_dif            ', Interstitial%sfc_alb_uvvis_dif            )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%toa_src_lw                   ', Interstitial%toa_src_lw                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%toa_src_sw                   ', Interstitial%toa_src_sw                   )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%active_gases_array           ', Interstitial%active_gases_array           )
                         ! These DDTs do not have print routines
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%flxprf_lw                    ', Interstitial%flxprf_lw                    )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%flxprf_sw                    ', Interstitial%flxprf_sw                    )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lw_optical_props_cloudsByBand', Interstitial%lw_optical_props_cloudsByBand)
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lw_optical_props_clouds      ', Interstitial%lw_optical_props_clouds      )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lw_optical_props_precipByBand', Interstitial%lw_optical_props_precipByBand)
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lw_optical_props_precip      ', Interstitial%lw_optical_props_precip      )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lw_optical_props_clrsky      ', Interstitial%lw_optical_props_clrsky      )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lw_optical_props_aerosol     ', Interstitial%lw_optical_props_aerosol     )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sw_optical_props_cloudsByBand', Interstitial%sw_optical_props_cloudsByBand)
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sw_optical_props_clouds      ', Interstitial%sw_optical_props_clouds      )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sw_optical_props_precipByBand', Interstitial%sw_optical_props_precipByBand)
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sw_optical_props_precip      ', Interstitial%sw_optical_props_precip      )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sw_optical_props_clrsky      ', Interstitial%sw_optical_props_clrsky      )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sw_optical_props_aerosol     ', Interstitial%sw_optical_props_aerosol     )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gas_concentrations           ', Interstitial%gas_concentrations           )
                         !call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sources                      ', Interstitial%sources                      )
                     end if
                 end if
#ifdef OPENMP
!$OMP BARRIER
#endif
             end do
#ifdef MPI
!             call MPI_BARRIER(mpicomm,ierr)
#endif
         end do

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
!         call MPI_BARRIER(mpicomm,ierr)
#endif

      end subroutine GFS_interstitialtoscreen_run

    end module GFS_interstitialtoscreen

    module GFS_abort

      private

      public GFS_abort_init, GFS_abort_run, GFS_abort_finalize

      contains

      subroutine GFS_abort_init ()
      end subroutine GFS_abort_init

      subroutine GFS_abort_finalize ()
      end subroutine GFS_abort_finalize

!> \section arg_table_GFS_abort_run Argument Table
!! \htmlinclude GFS_abort_run.html
!!
      subroutine GFS_abort_run (Model, blkno, errmsg, errflg)

         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in   ) :: Model
         integer,                  intent(in   ) :: blkno
         character(len=*),         intent(  out) :: errmsg
         integer,                  intent(  out) :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (Model%kdt==1 .and. blkno==size(Model%blksz)) then
             if (Model%me==Model%master) write(0,*) "GFS_abort_run: ABORTING MODEL"
             call sleep(10)
             stop
         end if

      end subroutine GFS_abort_run

    end module GFS_abort

    module GFS_checkland

      private

      public GFS_checkland_init, GFS_checkland_run, GFS_checkland_finalize

      contains

      subroutine GFS_checkland_init ()
      end subroutine GFS_checkland_init

      subroutine GFS_checkland_finalize ()
      end subroutine GFS_checkland_finalize

!> \section arg_table_GFS_checkland_run Argument Table
!! \htmlinclude GFS_checkland_run.html
!!
      subroutine GFS_checkland_run (me, master, blkno, im, kdt, iter, flag_iter, flag_guess, &
              flag_init, flag_restart, frac_grid, isot, ivegsrc, stype, vtype, slope,        &
              soiltyp, vegtype, slopetyp, dry, icy, wet, lake, ocean,                        &
              oceanfrac, landfrac, lakefrac, slmsk, islmsk,                                  &
              zorl, zorlw, zorll, zorli, fice, errmsg, errflg )

         use machine, only: kind_phys

         implicit none

         ! Interface variables
         integer,          intent(in   ) :: me
         integer,          intent(in   ) :: master
         integer,          intent(in   ) :: blkno
         integer,          intent(in   ) :: im
         integer,          intent(in   ) :: kdt
         integer,          intent(in   ) :: iter
         logical,          intent(in   ) :: flag_iter(im)
         logical,          intent(in   ) :: flag_guess(im)
         logical,          intent(in   ) :: flag_init
         logical,          intent(in   ) :: flag_restart
         logical,          intent(in   ) :: frac_grid
         integer,          intent(in   ) :: isot
         integer,          intent(in   ) :: ivegsrc
         real(kind_phys),  intent(in   ) :: stype(im)
         real(kind_phys),  intent(in   ) :: vtype(im)
         real(kind_phys),  intent(in   ) :: slope(im)
         integer,          intent(in   ) :: soiltyp(im)
         integer,          intent(in   ) :: vegtype(im)
         integer,          intent(in   ) :: slopetyp(im)
         logical,          intent(in   ) :: dry(im)
         logical,          intent(in   ) :: icy(im)
         logical,          intent(in   ) :: wet(im)
         logical,          intent(in   ) :: lake(im)
         logical,          intent(in   ) :: ocean(im)
         real(kind_phys),  intent(in   ) :: oceanfrac(im)
         real(kind_phys),  intent(in   ) :: landfrac(im)
         real(kind_phys),  intent(in   ) :: lakefrac(im)
         real(kind_phys),  intent(in   ) :: slmsk(im)
         integer,          intent(in   ) :: islmsk(im)
         real(kind_phys),  intent(in   ) :: zorl(im)
         real(kind_phys),  intent(in   ) :: zorlw(im)
         real(kind_phys),  intent(in   ) :: zorll(im)
         real(kind_phys),  intent(in   ) :: zorli(im)
         real(kind_phys),  intent(in   ) :: fice(im)
         character(len=*), intent(  out) :: errmsg
         integer,          intent(  out) :: errflg

         ! Local variables
         integer :: i

         errflg = 0
         errmsg = ''

         write(0,'(a,i5)')   'YYY: me           :', me
         write(0,'(a,i5)')   'YYY: master       :', master
         write(0,'(a,i5)')   'YYY: blkno        :', blkno
         write(0,'(a,i5)')   'YYY: im           :', im
         write(0,'(a,i5)')   'YYY: kdt          :', kdt
         write(0,'(a,i5)')   'YYY: iter         :', iter
         write(0,'(a,1x,l)') 'YYY: flag_init    :', flag_init
         write(0,'(a,1x,l)') 'YYY: flag_restart :', flag_restart
         write(0,'(a,1x,l)') 'YYY: frac_grid    :', frac_grid
         write(0,'(a,i5)')   'YYY: isot         :', isot
         write(0,'(a,i5)')   'YYY: ivegsrc      :', ivegsrc

         do i=1,im
           !if (fice(i)>0.999) then
           !if (vegtype(i)==15) then
             write(0,'(a,2i5,1x,1x,l)') 'YYY: i, blk, flag_iter(i)  :', i, blkno, flag_iter(i)
             write(0,'(a,2i5,1x,1x,l)') 'YYY: i, blk, flag_guess(i) :', i, blkno, flag_guess(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, stype(i)      :', i, blkno, stype(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, vtype(i)      :', i, blkno, vtype(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, slope(i)      :', i, blkno, slope(i)
             write(0,'(a,2i5,1x,i5)')   'YYY: i, blk, soiltyp(i)    :', i, blkno, soiltyp(i)
             write(0,'(a,2i5,1x,i5)')   'YYY: i, blk, vegtype(i)    :', i, blkno, vegtype(i)
             write(0,'(a,2i5,1x,i5)')   'YYY: i, blk, slopetyp(i)   :', i, blkno, slopetyp(i)
             write(0,'(a,2i5,1x,1x,l)') 'YYY: i, blk, dry(i)        :', i, blkno, dry(i)
             write(0,'(a,2i5,1x,1x,l)') 'YYY: i, blk, icy(i)        :', i, blkno, icy(i)
             write(0,'(a,2i5,1x,1x,l)') 'YYY: i, blk, wet(i)        :', i, blkno, wet(i)
             write(0,'(a,2i5,1x,1x,l)') 'YYY: i, blk, lake(i)       :', i, blkno, lake(i)
             write(0,'(a,2i5,1x,1x,l)') 'YYY: i, blk, ocean(i)      :', i, blkno, ocean(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, oceanfrac(i)  :', i, blkno, oceanfrac(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, landfrac(i)   :', i, blkno, landfrac(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, lakefrac(i)   :', i, blkno, lakefrac(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, fice(i)       :', i, blkno, fice(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, slmsk(i)      :', i, blkno, slmsk(i)
             write(0,'(a,2i5,1x,i5)')   'YYY: i, blk, islmsk(i)     :', i, blkno, islmsk(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, zorl(i)       :', i, blkno, zorl(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, zorlw(i)      :', i, blkno, zorlw(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, zorli(i)      :', i, blkno, zorli(i)
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, zorll(i)      :', i, blkno, zorll(i)
           !end if
         end do

      end subroutine GFS_checkland_run

    end module GFS_checkland

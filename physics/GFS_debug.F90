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

      subroutine GFS_diagtoscreen_init ()
      end subroutine GFS_diagtoscreen_init

      subroutine GFS_diagtoscreen_finalize ()
      end subroutine GFS_diagtoscreen_finalize

!> \section arg_table_GFS_diagtoscreen_run Argument Table
!! \htmlinclude GFS_diagtoscreen_run.html
!!
      subroutine GFS_diagtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling,     &
                                       Grid, Tbd, Cldprop, Radtend, Diag, Interstitial, &
                                       nthreads, blkno, errmsg, errflg)

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
                                          GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in   ) :: Model
         type(GFS_statein_type),   intent(in   ) :: Statein
         type(GFS_stateout_type),  intent(in   ) :: Stateout
         type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
         type(GFS_coupling_type),  intent(in   ) :: Coupling
         type(GFS_grid_type),      intent(in   ) :: Grid
         type(GFS_tbd_type),       intent(in   ) :: Tbd
         type(GFS_cldprop_type),   intent(in   ) :: Cldprop
         type(GFS_radtend_type),   intent(in   ) :: Radtend
         type(GFS_diag_type),      intent(in   ) :: Diag
         type(GFS_interstitial_type), intent(in) :: Interstitial
         integer,                  intent(in   ) :: nthreads
         integer,                  intent(in   ) :: blkno
         character(len=*),           intent(out) :: errmsg
         integer,                    intent(out) :: errflg

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
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zorlo'    , Sfcprop%zorlo)
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zorll'    , Sfcprop%zorll)
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
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%sh2o',         Sfcprop%sh2o)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%smois',        Sfcprop%smois)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tslb',         Sfcprop%tslb)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%zs',           Sfcprop%zs)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%clw_surf',     Sfcprop%clw_surf)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%qwv_surf',     Sfcprop%qwv_surf)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%cndm_surf',    Sfcprop%cndm_surf)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%flag_frsoil',  Sfcprop%flag_frsoil)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%rhofr',        Sfcprop%rhofr)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%tsnow',        Sfcprop%tsnow)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%snowfallac  ', Sfcprop%snowfallac)
                        call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Sfcprop%acsnow      ', Sfcprop%acsnow)
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
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Tbd%dtdtr'         , Tbd%dtdtr)
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
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%upd_mf      ',    Diag%upd_mf)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%dwn_mf      ',    Diag%dwn_mf)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%det_mf      ',    Diag%det_mf)
                       call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Diag%cldcov      ',    Diag%cldcov)
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
                     ! Model/Control
                     ! not yet
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


    module GFS_interstitialtoscreen

      use print_var_chksum, only: print_var

      implicit none

      private

      public GFS_interstitialtoscreen_init, GFS_interstitialtoscreen_run, GFS_interstitialtoscreen_finalize

      contains

      subroutine GFS_interstitialtoscreen_init ()
      end subroutine GFS_interstitialtoscreen_init

      subroutine GFS_interstitialtoscreen_finalize ()
      end subroutine GFS_interstitialtoscreen_finalize

!> \section arg_table_GFS_interstitialtoscreen_run Argument Table
!! \htmlinclude GFS_interstitialtoscreen_run.html
!!
      subroutine GFS_interstitialtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling, &
                                           Grid, Tbd, Cldprop, Radtend, Diag, Interstitial, &
                                           nthreads, blkno, errmsg, errflg)

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
                                          GFS_interstitial_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in   ) :: Model
         type(GFS_statein_type),   intent(in   ) :: Statein
         type(GFS_stateout_type),  intent(in   ) :: Stateout
         type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
         type(GFS_coupling_type),  intent(in   ) :: Coupling
         type(GFS_grid_type),      intent(in   ) :: Grid
         type(GFS_tbd_type),       intent(in   ) :: Tbd
         type(GFS_cldprop_type),   intent(in   ) :: Cldprop
         type(GFS_radtend_type),   intent(in   ) :: Radtend
         type(GFS_diag_type),      intent(in   ) :: Diag
         type(GFS_interstitial_type), intent(in) :: Interstitial
         integer,                  intent(in   ) :: nthreads
         integer,                  intent(in   ) :: blkno
         character(len=*),         intent(  out) :: errmsg
         integer,                  intent(  out) :: errflg

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
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%h2o_coeff           ', Interstitial%h2o_coeff               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%h2o_pres            ', Interstitial%h2o_pres                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ipr                 ', Interstitial%ipr                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%itc                 ', Interstitial%itc                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%latidxprnt          ', Interstitial%latidxprnt              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%levi                ', Interstitial%levi                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%levh2o              ', Interstitial%levh2o                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%levozp              ', Interstitial%levozp                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lmk                 ', Interstitial%lmk                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lmp                 ', Interstitial%lmp                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nbdlw               ', Interstitial%nbdlw                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nbdsw               ', Interstitial%nbdsw                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nf_aelw             ', Interstitial%nf_aelw                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nf_aesw             ', Interstitial%nf_aesw                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nsamftrac           ', Interstitial%nsamftrac               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nscav               ', Interstitial%nscav                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nspc1               ', Interstitial%nspc1                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ntiwx               ', Interstitial%ntiwx                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nvdiff              ', Interstitial%nvdiff                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%oz_coeff            ', Interstitial%oz_coeff                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'sum(Interstitial%oz_pres)        ', Interstitial%oz_pres                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%phys_hydrostatic    ', Interstitial%phys_hydrostatic        )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%skip_macro          ', Interstitial%skip_macro              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%trans_aero          ', Interstitial%trans_aero              )
                     ! Print all other variables
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjsfculw_land      ', Interstitial%adjsfculw_land          )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjsfculw_ice       ', Interstitial%adjsfculw_ice           )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjsfculw_ocean     ', Interstitial%adjsfculw_ocean         )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjnirbmd           ', Interstitial%adjnirbmd               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjnirbmu           ', Interstitial%adjnirbmu               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjnirdfd           ', Interstitial%adjnirdfd               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjnirdfu           ', Interstitial%adjnirdfu               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjvisbmd           ', Interstitial%adjvisbmd               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjvisbmu           ', Interstitial%adjvisbmu               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjvisdfu           ', Interstitial%adjvisdfu               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%adjvisdfd           ', Interstitial%adjvisdfd               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%aerodp              ', Interstitial%aerodp                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%alb1d               ', Interstitial%alb1d                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%bexp1d              ', Interstitial%bexp1d                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cd                  ', Interstitial%cd                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cd_ice              ', Interstitial%cd_ice                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cd_land             ', Interstitial%cd_land                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cd_ocean            ', Interstitial%cd_ocean                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cdq                 ', Interstitial%cdq                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cdq_ice             ', Interstitial%cdq_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cdq_land            ', Interstitial%cdq_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cdq_ocean           ', Interstitial%cdq_ocean               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%chh_ice             ', Interstitial%chh_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%chh_land            ', Interstitial%chh_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%chh_ocean           ', Interstitial%chh_ocean               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cldf                ', Interstitial%cldf                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cldsa               ', Interstitial%cldsa                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cldtaulw            ', Interstitial%cldtaulw                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cldtausw            ', Interstitial%cldtausw                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cld1d               ', Interstitial%cld1d                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%clw                 ', Interstitial%clw                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%clx                 ', Interstitial%clx                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%clouds              ', Interstitial%clouds                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cmm_ice             ', Interstitial%cmm_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cmm_land            ', Interstitial%cmm_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cmm_ocean           ', Interstitial%cmm_ocean               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cnvc                ', Interstitial%cnvc                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cnvw                ', Interstitial%cnvw                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ctei_r              ', Interstitial%ctei_r                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ctei_rml            ', Interstitial%ctei_rml                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cumabs              ', Interstitial%cumabs                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dd_mf               ', Interstitial%dd_mf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%de_lgth             ', Interstitial%de_lgth                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%del                 ', Interstitial%del                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%del_gz              ', Interstitial%del_gz                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%delr                ', Interstitial%delr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dkt                 ', Interstitial%dkt                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dlength             ', Interstitial%dlength                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dqdt                ', Interstitial%dqdt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dqsfc1              ', Interstitial%dqsfc1                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%drain               ', Interstitial%drain                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dtdt                ', Interstitial%dtdt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dtdtc               ', Interstitial%dtdtc                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dtsfc1              ', Interstitial%dtsfc1                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dtzm                ', Interstitial%dtzm                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dt_mf               ', Interstitial%dt_mf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dudt                ', Interstitial%dudt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dusfcg              ', Interstitial%dusfcg                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dusfc1              ', Interstitial%dusfc1                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dvdftra             ', Interstitial%dvdftra                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dvdt                ', Interstitial%dvdt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dvsfcg              ', Interstitial%dvsfcg                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dvsfc1              ', Interstitial%dvsfc1                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dzlyr               ', Interstitial%dzlyr                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%elvmax              ', Interstitial%elvmax                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ep1d                ', Interstitial%ep1d                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ep1d_ice            ', Interstitial%ep1d_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ep1d_land           ', Interstitial%ep1d_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ep1d_ocean          ', Interstitial%ep1d_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%evapq               ', Interstitial%evapq                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%evap_ice            ', Interstitial%evap_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%evap_land           ', Interstitial%evap_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%evap_ocean          ', Interstitial%evap_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%evbs                ', Interstitial%evbs                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%evcw                ', Interstitial%evcw                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%faerlw              ', Interstitial%faerlw                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%faersw              ', Interstitial%faersw                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ffhh_ice            ', Interstitial%ffhh_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ffhh_land           ', Interstitial%ffhh_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ffhh_ocean          ', Interstitial%ffhh_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fh2                 ', Interstitial%fh2                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fh2_ice             ', Interstitial%fh2_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fh2_land            ', Interstitial%fh2_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fh2_ocean           ', Interstitial%fh2_ocean               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%flag_cice           ', Interstitial%flag_cice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%flag_guess          ', Interstitial%flag_guess              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%flag_iter           ', Interstitial%flag_iter               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ffmm_ice            ', Interstitial%ffmm_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ffmm_land           ', Interstitial%ffmm_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ffmm_ocean          ', Interstitial%ffmm_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fm10                ', Interstitial%fm10                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fm10_ice            ', Interstitial%fm10_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fm10_land           ', Interstitial%fm10_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fm10_ocean          ', Interstitial%fm10_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%frain               ', Interstitial%frain                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%frland              ', Interstitial%frland                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fscav               ', Interstitial%fscav                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%fswtr               ', Interstitial%fswtr                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gabsbdlw            ', Interstitial%gabsbdlw                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gabsbdlw_ice        ', Interstitial%gabsbdlw_ice            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gabsbdlw_land       ', Interstitial%gabsbdlw_land           )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gabsbdlw_ocean      ', Interstitial%gabsbdlw_ocean          )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gamma               ', Interstitial%gamma                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gamq                ', Interstitial%gamq                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gamt                ', Interstitial%gamt                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gasvmr              ', Interstitial%gasvmr                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gflx                ', Interstitial%gflx                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gflx_ice            ', Interstitial%gflx_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gflx_land           ', Interstitial%gflx_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gflx_ocean          ', Interstitial%gflx_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gwdcu               ', Interstitial%gwdcu                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gwdcv               ', Interstitial%gwdcv                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%hefac               ', Interstitial%hefac                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%hffac               ', Interstitial%hffac                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%hflxq               ', Interstitial%hflxq                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%hflx_ice            ', Interstitial%hflx_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%hflx_land           ', Interstitial%hflx_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%hflx_ocean          ', Interstitial%hflx_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%htlwc               ', Interstitial%htlwc                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%htlw0               ', Interstitial%htlw0                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%htswc               ', Interstitial%htswc                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%htsw0               ', Interstitial%htsw0                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dry                 ', Interstitial%dry                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%idxday              ', Interstitial%idxday                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%icy                 ', Interstitial%icy                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%lake                ', Interstitial%lake                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ocean               ', Interstitial%ocean                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%islmsk              ', Interstitial%islmsk                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%islmsk_cice         ', Interstitial%islmsk_cice             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%wet                 ', Interstitial%wet                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%kb                  ', Interstitial%kb                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%kbot                ', Interstitial%kbot                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%kcnv                ', Interstitial%kcnv                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%kd                  ', Interstitial%kd                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%kinver              ', Interstitial%kinver                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%kpbl                ', Interstitial%kpbl                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%kt                  ', Interstitial%kt                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ktop                ', Interstitial%ktop                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%mbota               ', Interstitial%mbota                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%mtopa               ', Interstitial%mtopa                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%nday                ', Interstitial%nday                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%oa4                 ', Interstitial%oa4                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%oc                  ', Interstitial%oc                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%olyr                ', Interstitial%olyr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%plvl                ', Interstitial%plvl                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%plyr                ', Interstitial%plyr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%prcpmp              ', Interstitial%prcpmp                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%prnum               ', Interstitial%prnum                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qlyr                ', Interstitial%qlyr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qss_ice             ', Interstitial%qss_ice                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qss_land            ', Interstitial%qss_land                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qss_ocean           ', Interstitial%qss_ocean               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%radar_reset         ', Interstitial%radar_reset             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%raddt               ', Interstitial%raddt                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%raincd              ', Interstitial%raincd                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%raincs              ', Interstitial%raincs                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%rainmcadj           ', Interstitial%rainmcadj               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%rainp               ', Interstitial%rainp                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%rb                  ', Interstitial%rb                      )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%rb_ice              ', Interstitial%rb_ice                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%rb_land             ', Interstitial%rb_land                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%rb_ocean            ', Interstitial%rb_ocean                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%reset               ', Interstitial%reset                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%rhc                 ', Interstitial%rhc                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%runoff              ', Interstitial%runoff                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%save_q              ', Interstitial%save_q                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%save_t              ', Interstitial%save_t                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%save_tcp            ', Interstitial%save_tcp                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%save_u              ', Interstitial%save_u                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%save_v              ', Interstitial%save_v                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sbsno               ', Interstitial%sbsno                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%scmpsw%uvbfc        ', Interstitial%scmpsw%uvbfc            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%scmpsw%uvbf0        ', Interstitial%scmpsw%uvbf0            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%scmpsw%nirbm        ', Interstitial%scmpsw%nirbm            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%scmpsw%nirdf        ', Interstitial%scmpsw%nirdf            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%scmpsw%visbm        ', Interstitial%scmpsw%visbm            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%scmpsw%visdf        ', Interstitial%scmpsw%visdf            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%semis_ice           ', Interstitial%semis_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%semis_land          ', Interstitial%semis_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%semis_ocean         ', Interstitial%semis_ocean             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sfcalb              ', Interstitial%sfcalb                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sigma               ', Interstitial%sigma                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sigmaf              ', Interstitial%sigmaf                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sigmafrac           ', Interstitial%sigmafrac               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%sigmatot            ', Interstitial%sigmatot                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%slopetype           ', Interstitial%slopetype               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%snowc               ', Interstitial%snowc                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%snowd_ice           ', Interstitial%snowd_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%snowd_land          ', Interstitial%snowd_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%snowd_ocean         ', Interstitial%snowd_ocean             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%snohf               ', Interstitial%snohf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%snowmt              ', Interstitial%snowmt                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%soiltype            ', Interstitial%soiltype                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%stress              ', Interstitial%stress                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%stress_ice          ', Interstitial%stress_ice              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%stress_land         ', Interstitial%stress_land             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%stress_ocean        ', Interstitial%stress_ocean            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%theta               ', Interstitial%theta                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tice                ', Interstitial%tice                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tlvl                ', Interstitial%tlvl                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tlyr                ', Interstitial%tlyr                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tprcp_ice           ', Interstitial%tprcp_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tprcp_land          ', Interstitial%tprcp_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tprcp_ocean         ', Interstitial%tprcp_ocean             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%trans               ', Interstitial%trans                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tseal               ', Interstitial%tseal                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsfa                ', Interstitial%tsfa                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsfc_ice            ', Interstitial%tsfc_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsfc_land           ', Interstitial%tsfc_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsfc_ocean          ', Interstitial%tsfc_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsfg                ', Interstitial%tsfg                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsurf               ', Interstitial%tsurf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsurf_ice           ', Interstitial%tsurf_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsurf_land          ', Interstitial%tsurf_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tsurf_ocean         ', Interstitial%tsurf_ocean             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ud_mf               ', Interstitial%ud_mf                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%uustar_ice          ', Interstitial%uustar_ice              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%uustar_land         ', Interstitial%uustar_land             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%uustar_ocean        ', Interstitial%uustar_ocean            )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%vdftra              ', Interstitial%vdftra                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%vegf1d              ', Interstitial%vegf1d                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%vegtype             ', Interstitial%vegtype                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%wcbmax              ', Interstitial%wcbmax                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%weasd_ice           ', Interstitial%weasd_ice               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%weasd_land          ', Interstitial%weasd_land              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%weasd_ocean         ', Interstitial%weasd_ocean             )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%wind                ', Interstitial%wind                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%work1               ', Interstitial%work1                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%work2               ', Interstitial%work2                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%work3               ', Interstitial%work3                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%xcosz               ', Interstitial%xcosz                   )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%xlai1d              ', Interstitial%xlai1d                  )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%xmu                 ', Interstitial%xmu                     )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%z01d                ', Interstitial%z01d                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%zorl_ice            ', Interstitial%zorl_ice                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%zorl_land           ', Interstitial%zorl_land               )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%zorl_ocean          ', Interstitial%zorl_ocean              )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%zt1d                ', Interstitial%zt1d                    )
                     ! CIRES UGWP v0
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gw_dudt             ', Interstitial%gw_dudt                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gw_dvdt             ', Interstitial%gw_dvdt                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gw_dtdt             ', Interstitial%gw_dtdt                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%gw_kdis             ', Interstitial%gw_kdis                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tau_mtb             ', Interstitial%tau_mtb                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tau_ogw             ', Interstitial%tau_ogw                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tau_tofd            ', Interstitial%tau_tofd                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%tau_ngw             ', Interstitial%tau_ngw                 )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%zmtb                ', Interstitial%zmtb                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%zlwb                ', Interstitial%zlwb                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%zogw                ', Interstitial%zogw                    )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dudt_mtb            ', Interstitial%dudt_mtb                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dudt_ogw            ', Interstitial%dudt_ogw                )
                     call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%dudt_tms            ', Interstitial%dudt_tms                )
                     !-- GSD drag suite
                     if (Model%gwd_opt==3 .or. Model%gwd_opt==33) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%varss               ', Interstitial%varss                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ocss                ', Interstitial%ocss                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%oa4ss               ', Interstitial%oa4ss                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%clxss               ', Interstitial%clxss                   )
                     end if
                     ! GFDL and Thompson MP
                     if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%graupelmp           ', Interstitial%graupelmp               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%icemp               ', Interstitial%icemp                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%rainmp              ', Interstitial%rainmp                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%snowmp              ', Interstitial%snowmp                  )
                     ! Ferrier-Aligo
                     else if (Model%imp_physics == Model%imp_physics_fer_hires) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%f_ice               ', Interstitial%f_ice                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%f_rain              ', Interstitial%f_rain                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%f_rimef             ', Interstitial%f_rimef                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cwm                 ', Interstitial%cwm                     )
                     ! Morrison-Gettelman
                     else if (Model%imp_physics == Model%imp_physics_mg) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ncgl                ', Interstitial%ncgl                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ncpr                ', Interstitial%ncpr                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ncps                ', Interstitial%ncps                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qgl                 ', Interstitial%qgl                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qrn                 ', Interstitial%qrn                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qsnw                ', Interstitial%qsnw                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qlcn                ', Interstitial%qlcn                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qicn                ', Interstitial%qicn                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%w_upi               ', Interstitial%w_upi                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cf_upi              ', Interstitial%cf_upi                  )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cnv_mfd             ', Interstitial%cnv_mfd                 )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cnv_dqldt           ', Interstitial%cnv_dqldt               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%clcn                ', Interstitial%clcn                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cnv_fice            ', Interstitial%cnv_fice                )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cnv_ndrop           ', Interstitial%cnv_ndrop               )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%cnv_nice            ', Interstitial%cnv_nice                )
                     end if
                     ! SHOC
                     if (Model%do_shoc) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ncgl                ', Interstitial%ncgl                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qrn                 ', Interstitial%qrn                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qsnw                ', Interstitial%qsnw                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%qgl                 ', Interstitial%qgl                     )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ncpi                ', Interstitial%ncpi                    )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%ncpl                ', Interstitial%ncpl                    )
                     end if
                     ! Noah MP
                     if (Model%lsm == Model%lsm_noahmp) then
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%t2mmp               ', Interstitial%t2mmp                   )
                         call print_var(mpirank, omprank, blkno, Grid%xlat_d, Grid%xlon_d, 'Interstitial%q2mp                ', Interstitial%q2mp                    )
                     end if
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
              oceanfrac, landfrac, lakefrac, slmsk, islmsk, errmsg, errflg )

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
             write(0,'(a,2i5,1x,e16.7)')'YYY: i, blk, slmsk(i)      :', i, blkno, slmsk(i)
             write(0,'(a,2i5,1x,i5)')   'YYY: i, blk, islmsk(i)     :', i, blkno, islmsk(i)
           !end if
         end do

      end subroutine GFS_checkland_run

    end module GFS_checkland

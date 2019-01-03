!> \file GFS_debug.F90

    module GFS_diagtoscreen

      private
 
      public GFS_diagtoscreen_init, GFS_diagtoscreen_run, GFS_diagtoscreen_finalize

      public print_my_stuff, chksum_int, chksum_real

#define PRINT_CHKSUM
!#define PRINT_SUM

      interface print_var
        module procedure print_logic_0d
        module procedure print_int_0d
        module procedure print_int_1d
        module procedure print_real_0d
        module procedure print_real_1d
        module procedure print_real_2d
        module procedure print_real_3d
      end interface

      integer, parameter :: ISTART = 1
      integer, parameter :: IEND = 9999999

      integer, parameter :: KSTART = 1
      integer, parameter :: KEND = 9999999

      contains

      subroutine GFS_diagtoscreen_init ()
      end subroutine GFS_diagtoscreen_init

      subroutine GFS_diagtoscreen_finalize ()
      end subroutine GFS_diagtoscreen_finalize

!> \section arg_table_GFS_diagtoscreen_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | derived type GFS_control_type in FV3                    | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                   | derived type GFS_statein_type in FV3                    | DDT           |    0 | GFS_statein_type      |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | derived type GFS_stateout_type in FV3                   | DDT           |    0 | GFS_stateout_type     |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                   | derived type GFS_sfcprop_type in FV3                    | DDT           |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Coupling       | FV3-GFS_Coupling_type                                  | derived type GFS_coupling_type in FV3                   | DDT           |    0 | GFS_coupling_type     |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | derived type GFS_grid_type in FV3                       | DDT           |    0 | GFS_grid_type         |           | in     | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                       | derived type GFS_tbd_type in FV3                        | DDT           |    0 | GFS_tbd_type          |           | in     | F        |
!! | Cldprop        | FV3-GFS_Cldprop_type                                   | derived type GFS_cldprop_type in FV3                    | DDT           |    0 | GFS_cldprop_type      |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                   | derived type GFS_radtend_type in FV3                    | DDT           |    0 | GFS_radtend_type      |           | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                      | derived type GFS_diag_type in FV3                       | DDT           |    0 | GFS_diag_type         |           | in     | F        |
!! | Interstitial   | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | in     | F        |
!! | nthreads       | omp_threads                                            | number of OpenMP threads or fast physics schemes        | count         |    0 | integer               |           | in     | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_diagtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling,     &
                                       Grid, Tbd, Cldprop, Radtend, Diag, Interstitial, &
                                       nthreads, errmsg, errflg)

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
                     ! Sfcprop
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%slmsk'    , Sfcprop%slmsk)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%lakemsk'  , Sfcprop%lakemsk)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%tsfc'     , Sfcprop%tsfc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%tisfc'    , Sfcprop%tisfc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%snowd'    , Sfcprop%snowd)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%zorl'     , Sfcprop%zorl)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%fice'     , Sfcprop%fice)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%hprim'    , Sfcprop%hprim)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%hprime'   , Sfcprop%hprime)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%sncovr'   , Sfcprop%sncovr)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%snoalb'   , Sfcprop%snoalb)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%alvsf'    , Sfcprop%alvsf)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%alnsf'    , Sfcprop%alnsf)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%alvwf'    , Sfcprop%alvwf)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%alnwf'    , Sfcprop%alnwf)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%facsf'    , Sfcprop%facsf)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%facwf'    , Sfcprop%facwf)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%slope'    , Sfcprop%slope)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%shdmin'   , Sfcprop%shdmin)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%shdmax'   , Sfcprop%shdmax)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%tg3'      , Sfcprop%tg3)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%vfrac'    , Sfcprop%vfrac)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%vtype'    , Sfcprop%vtype)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%stype'    , Sfcprop%stype)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%uustar'   , Sfcprop%uustar)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%oro'      , Sfcprop%oro)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%oro_uf'   , Sfcprop%oro_uf)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%hice'     , Sfcprop%hice)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%weasd'    , Sfcprop%weasd)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%canopy'   , Sfcprop%canopy)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%ffmm'     , Sfcprop%ffmm)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%ffhh'     , Sfcprop%ffhh)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%f10m'     , Sfcprop%f10m)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%tprcp'    , Sfcprop%tprcp)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%srflag'   , Sfcprop%srflag)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%slc'      , Sfcprop%slc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%smc'      , Sfcprop%smc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%stc'      , Sfcprop%stc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%t2m'      , Sfcprop%t2m)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%q2m'      , Sfcprop%q2m)
                     if (Model%nstf_name(1)>0) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%tref    ', Sfcprop%tref)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%z_c     ', Sfcprop%z_c)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%c_0     ', Sfcprop%c_0)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%c_d     ', Sfcprop%c_d)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%w_0     ', Sfcprop%w_0)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%w_d     ', Sfcprop%w_d)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%xt      ', Sfcprop%xt)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%xs      ', Sfcprop%xs)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%xu      ', Sfcprop%xu)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%xv      ', Sfcprop%xv)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%xz      ', Sfcprop%xz)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%zm      ', Sfcprop%zm)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%xtts    ', Sfcprop%xtts)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%xzts    ', Sfcprop%xzts)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%d_conv  ', Sfcprop%d_conv)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%ifd     ', Sfcprop%ifd)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%dt_cool ', Sfcprop%dt_cool)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%qrain   ', Sfcprop%qrain)
                     end if
                     ! CCPP/RUC only
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%sh2o',        Sfcprop%sh2o)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%smois',       Sfcprop%smois)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%tslb',        Sfcprop%tslb)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%zs',          Sfcprop%zs)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%clw_surf',    Sfcprop%clw_surf)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%cndm_surf',   Sfcprop%cndm_surf)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%flag_frsoil', Sfcprop%flag_frsoil)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%rhofr',       Sfcprop%rhofr)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Sfcprop%tsnow',       Sfcprop%tsnow)
                     ! Radtend
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfcfsw%upfxc', Radtend%sfcfsw(:)%upfxc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfcfsw%dnfxc', Radtend%sfcfsw(:)%dnfxc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfcfsw%upfx0', Radtend%sfcfsw(:)%upfx0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfcfsw%dnfx0', Radtend%sfcfsw(:)%dnfx0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfcflw%upfxc', Radtend%sfcflw(:)%upfxc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfcflw%upfx0', Radtend%sfcflw(:)%upfx0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfcflw%dnfxc', Radtend%sfcflw(:)%dnfxc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfcflw%dnfx0', Radtend%sfcflw(:)%dnfx0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%htrsw',        Radtend%htrsw)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%htrlw',        Radtend%htrlw)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%sfalb',        Radtend%sfalb)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%coszen',       Radtend%coszen)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%tsflw',        Radtend%tsflw)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%semis',        Radtend%semis)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%coszdg',       Radtend%coszdg)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%swhc',         Radtend%swhc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%lwhc',         Radtend%lwhc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Radtend%lwhd',         Radtend%lwhd)
                     ! Tbd
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%icsdsw'          , Tbd%icsdsw)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%icsdlw'          , Tbd%icsdlw)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%ozpl'            , Tbd%ozpl)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%h2opl'           , Tbd%h2opl)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%rann'            , Tbd%rann)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%acv'             , Tbd%acv)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%acvb'            , Tbd%acvb)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%acvt'            , Tbd%acvt)
                     if (Model%do_sppt) then
                       call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%dtdtr'         , Tbd%dtdtr)
                       call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%dtotprcp'      , Tbd%dtotprcp)
                       call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%dcnvprcp'      , Tbd%dcnvprcp)
                       call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%drain_cpl'     , Tbd%drain_cpl)
                       call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%dsnow_cpl'     , Tbd%dsnow_cpl)
                     end if
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%phy_fctd'        , Tbd%phy_fctd)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%phy_f2d'         , Tbd%phy_f2d)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%phy_f3d'         , Tbd%phy_f3d)
                     do n=1,size(Tbd%phy_f3d(1,1,:))
                         call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%phy_f3d_n'   , Tbd%phy_f3d(:,:,n))
                     end do
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%in_nm'           , Tbd%in_nm)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%ccn_nm'          , Tbd%ccn_nm)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%aer_nm'          , Tbd%aer_nm)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Tbd%blkno'           , Tbd%blkno)
                     ! Diag (incomplete)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%topfsw%upfxc',    Diag%topfsw%upfxc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%topfsw%dnfxc',    Diag%topfsw%dnfxc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%topfsw%upfx0',    Diag%topfsw%upfx0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%topflw%upfxc',    Diag%topflw%upfxc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%topflw%upfx0',    Diag%topflw%upfx0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%dswsfci',         Diag%dswsfci)
                     !call print_var(mpirank,omprank, Tbd%blkno, 'Diag%nswsfci',         Diag%nswsfci)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%uswsfci',         Diag%uswsfci)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%dlwsfci',         Diag%dlwsfci)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%ulwsfci',         Diag%ulwsfci)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%dusfci',          Diag%dusfci)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%dvsfci',          Diag%dvsfci)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%dtsfci',          Diag%dtsfci)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%dqsfci',          Diag%dqsfci)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%gfluxi',          Diag%gfluxi)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%gflux',           Diag%gflux)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%epi'    ,         Diag%epi)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%gfluxi' ,         Diag%gfluxi)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%t1'     ,         Diag%t1)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%q1'     ,         Diag%q1)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%u1'     ,         Diag%u1)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%v1'     ,         Diag%v1)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%dlwsfc',          Diag%dlwsfc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%ulwsfc',          Diag%ulwsfc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Diag%psmean',          Diag%psmean)
                     ! Statein
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%phii'    ,     Statein%phii)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%prsi'    ,     Statein%prsi)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%prsik'   ,     Statein%prsik)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%phil'    ,     Statein%phil)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%prsl'    ,     Statein%prsl)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%prslk'   ,     Statein%prslk)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%pgr'     ,     Statein%pgr)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%ugrs'    ,     Statein%ugrs)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%vgrs'    ,     Statein%vgrs)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%vvl'     ,     Statein%vvl)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%tgrs'    ,     Statein%tgrs)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%qgrs'    ,     Statein%qgrs)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%qgrs-qv' ,     Statein%qgrs(:,:,1))
                     if (Model%ntoz>0) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Statein%qgrs-o3' ,     Statein%qgrs(:,:,Model%ntoz))
                     end if
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%diss_est',     Statein%diss_est)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%smc'     ,     Statein%smc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%stc'     ,     Statein%stc)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Statein%slc'     ,     Statein%slc)
                     ! Stateout
                     call print_var(mpirank,omprank, Tbd%blkno, 'Stateout%gu0',         Stateout%gu0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Stateout%gv0',         Stateout%gv0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Stateout%gt0',         Stateout%gt0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Stateout%gq0',         Stateout%gq0)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Stateout%gq0-qv' ,     Stateout%gq0(:,:,1))
                     if (Model%ntoz>0) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Stateout%gq0-o3' ,     Stateout%gq0(:,:,Model%ntoz))
                     end if
                     ! Coupling
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nirbmdi', Coupling%nirbmdi)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nirdfdi', Coupling%nirdfdi)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%visbmdi', Coupling%visbmdi)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%visdfdi', Coupling%visdfdi)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nirbmui', Coupling%nirbmui)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nirdfui', Coupling%nirdfui)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%visbmui', Coupling%visbmui)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%visdfui', Coupling%visdfui)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%sfcdsw ', Coupling%sfcdsw )
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%sfcnsw ', Coupling%sfcnsw )
                     call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%sfcdlw ', Coupling%sfcdlw )
                     if (Model%cplflx .or. Model%do_sppt) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%rain_cpl', Coupling%rain_cpl)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%snow_cpl', Coupling%snow_cpl)
                     end if
                     if (Model%cplflx) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%slimskin_cpl', Coupling%slimskin_cpl )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dusfcin_cpl ', Coupling%dusfcin_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dvsfcin_cpl ', Coupling%dvsfcin_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dtsfcin_cpl ', Coupling%dtsfcin_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dqsfcin_cpl ', Coupling%dqsfcin_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%ulwsfcin_cpl', Coupling%ulwsfcin_cpl )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%tseain_cpl  ', Coupling%tseain_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%tisfcin_cpl ', Coupling%tisfcin_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%ficein_cpl  ', Coupling%ficein_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%hicein_cpl  ', Coupling%hicein_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%hsnoin_cpl  ', Coupling%hsnoin_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dusfc_cpl   ', Coupling%dusfc_cpl    )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dvsfc_cpl   ', Coupling%dvsfc_cpl    )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dtsfc_cpl   ', Coupling%dtsfc_cpl    )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dqsfc_cpl   ', Coupling%dqsfc_cpl    )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dlwsfc_cpl  ', Coupling%dlwsfc_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dswsfc_cpl  ', Coupling%dswsfc_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dnirbm_cpl  ', Coupling%dnirbm_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dnirdf_cpl  ', Coupling%dnirdf_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dvisbm_cpl  ', Coupling%dvisbm_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dvisdf_cpl  ', Coupling%dvisdf_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nlwsfc_cpl  ', Coupling%nlwsfc_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nswsfc_cpl  ', Coupling%nswsfc_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nnirbm_cpl  ', Coupling%nnirbm_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nnirdf_cpl  ', Coupling%nnirdf_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nvisbm_cpl  ', Coupling%nvisbm_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nvisdf_cpl  ', Coupling%nvisdf_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dusfci_cpl  ', Coupling%dusfci_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dvsfci_cpl  ', Coupling%dvsfci_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dtsfci_cpl  ', Coupling%dtsfci_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dqsfci_cpl  ', Coupling%dqsfci_cpl   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dlwsfci_cpl ', Coupling%dlwsfci_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dswsfci_cpl ', Coupling%dswsfci_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dnirbmi_cpl ', Coupling%dnirbmi_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dnirdfi_cpl ', Coupling%dnirdfi_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dvisbmi_cpl ', Coupling%dvisbmi_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dvisdfi_cpl ', Coupling%dvisdfi_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nlwsfci_cpl ', Coupling%nlwsfci_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nswsfci_cpl ', Coupling%nswsfci_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nnirbmi_cpl ', Coupling%nnirbmi_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nnirdfi_cpl ', Coupling%nnirdfi_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nvisbmi_cpl ', Coupling%nvisbmi_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nvisdfi_cpl ', Coupling%nvisdfi_cpl  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%t2mi_cpl    ', Coupling%t2mi_cpl     )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%q2mi_cpl    ', Coupling%q2mi_cpl     )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%u10mi_cpl   ', Coupling%u10mi_cpl    )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%v10mi_cpl   ', Coupling%v10mi_cpl    )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%tsfci_cpl   ', Coupling%tsfci_cpl    )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%psurfi_cpl  ', Coupling%psurfi_cpl   )
                     end if
                     if (Model%cplchm) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%rain_cpl ', Coupling%rain_cpl )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%rainc_cpl', Coupling%rainc_cpl)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%ushfsfci ', Coupling%ushfsfci )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dkt      ', Coupling%dkt      )
                     end if
                     if (Model%do_sppt) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%sppt_wts', Coupling%sppt_wts)
                     end if
                     if (Model%do_shum) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%shum_wts', Coupling%shum_wts)
                     end if
                     if (Model%do_skeb) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%skebu_wts', Coupling%skebu_wts)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%skebv_wts', Coupling%skebv_wts)
                     end if
                     if (Model%do_sfcperts) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%sfc_wts', Coupling%sfc_wts)
                     end if
                     if (Model%lgocart .or. Model%ldiag3d) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dqdti  ', Coupling%dqdti  )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%cnvqci ', Coupling%cnvqci )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%upd_mfi', Coupling%upd_mfi)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%dwn_mfi', Coupling%dwn_mfi)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%det_mfi', Coupling%det_mfi)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%cldcovi', Coupling%cldcovi)
                     end if
                     if(Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nwfa2d', Coupling%nwfa2d)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Coupling%nifa2d', Coupling%nifa2d)
                     end if
                     ! Grid
                     call print_var(mpirank,omprank, Tbd%blkno, 'Grid%xlon  ', Grid%xlon  )
                     call print_var(mpirank,omprank, Tbd%blkno, 'Grid%xlat  ', Grid%xlat  )
                     call print_var(mpirank,omprank, Tbd%blkno, 'Grid%xlat_d', Grid%xlat_d)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Grid%sinlat', Grid%sinlat)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Grid%coslat', Grid%coslat)
                     call print_var(mpirank,omprank, Tbd%blkno, 'Grid%area  ', Grid%area  )
                     call print_var(mpirank,omprank, Tbd%blkno, 'Grid%dx    ', Grid%dx    )
                     if (Model%ntoz > 0) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Grid%ddy_o3   ', Grid%ddy_o3   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Grid%jindx1_o3', Grid%jindx1_o3)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Grid%jindx2_o3', Grid%jindx2_o3)
                     endif
                     if (Model%h2o_phys) then
                        call print_var(mpirank,omprank, Tbd%blkno, 'Grid%ddy_h   ', Grid%ddy_h   )
                        call print_var(mpirank,omprank, Tbd%blkno, 'Grid%jindx1_h', Grid%jindx1_h)
                        call print_var(mpirank,omprank, Tbd%blkno, 'Grid%jindx2_h', Grid%jindx2_h)
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

      subroutine print_logic_0d(mpirank,omprank,blkno,name,var)

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          logical, intent(in) :: var

          write(0,'(2a,3i6,1x,l)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_logic_0d

      subroutine print_int_0d(mpirank,omprank,blkno,name,var)

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          integer, intent(in) :: var

          write(0,'(2a,3i6,i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_int_0d

      subroutine print_int_1d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          integer, intent(in) :: var(:)

          integer :: i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_int(size(var),var), minval(var), maxval(var)
#else
          do i=ISTART,min(IEND,size(var(:)))
              write(0,'(2a,3i6,i6,i15)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, var(i)
          end do
#endif

      end subroutine print_int_1d

      subroutine print_real_0d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var

          write(0,'(2a,3i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, var

      end subroutine print_real_0d

      subroutine print_real_1d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:)

          integer :: i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),var), minval(var), maxval(var)
#else
          do i=ISTART,min(IEND,size(var(:)))
              write(0,'(2a,3i6,i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, var(i)
          end do
#endif

      end subroutine print_real_1d

      subroutine print_real_2d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:)
          
          integer :: k, i

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),reshape(var,(/size(var)/))), minval(var), maxval(var)
#else
          do i=ISTART,min(IEND,size(var(:,1)))
              do k=KSTART,min(KEND,size(var(1,:)))
                  write(0,'(2a,3i6,2i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, var(i,k)
              end do
          end do
#endif

      end subroutine print_real_2d

      subroutine print_real_3d(mpirank,omprank,blkno,name,var)

          use machine,               only: kind_phys

          implicit none

          integer, intent(in) :: mpirank, omprank, blkno
          character(len=*), intent(in) :: name
          real(kind_phys), intent(in) :: var(:,:,:)

          integer :: k, i, l

#ifdef PRINT_SUM
          write(0,'(2a,3i6,3e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, sum(var), minval(var), maxval(var)
#elif defined(PRINT_CHKSUM)
          write(0,'(2a,3i6,i20,2e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, chksum_real(size(var),reshape(var,(/size(var)/))), minval(var), maxval(var)
#else
          do i=ISTART,min(IEND,size(var(:,1,1)))
              do k=KSTART,min(KEND,size(var(1,:,1)))
                  do l=1,size(var(1,1,:))
                      write(0,'(2a,3i6,3i6,e35.25)') 'XXX: ', trim(name), mpirank, omprank, blkno, i, k, l, var(i,k,l)
                  end do
              end do
          end do
#endif

      end subroutine print_real_3d

      function chksum_int(N, var) result(hash)
          implicit none
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
          use machine,               only: kind_phys
          implicit none
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

     function print_my_stuff(mpitoprint,omptoprint) result(flag)
#ifdef MPI
         use mpi
#endif
#ifdef OPENMP
         use omp_lib
#endif
         implicit none
         integer, intent(in) :: mpitoprint, omptoprint
         logical :: flag
         integer :: ompthread, mpirank, ierr
#ifdef MPI
         call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, ierr)
#else
         mpirank = 0
#endif
#ifdef OPENMP
         ompthread = OMP_GET_THREAD_NUM()
#else
         ompthread = 0
#endif

         if (mpitoprint==mpirank .and. omptoprint==ompthread) then
             flag = .true.
         else
             flag = .false.
         end if
     end function print_my_stuff

    end module GFS_diagtoscreen


    module GFS_interstitialtoscreen

      private
 
      public GFS_interstitialtoscreen_init, GFS_interstitialtoscreen_run, GFS_interstitialtoscreen_finalize

      contains

      subroutine GFS_interstitialtoscreen_init ()
      end subroutine GFS_interstitialtoscreen_init

      subroutine GFS_interstitialtoscreen_finalize ()
      end subroutine GFS_interstitialtoscreen_finalize

!> \section arg_table_GFS_interstitialtoscreen_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | derived type GFS_control_type in FV3                    | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                   | derived type GFS_statein_type in FV3                    | DDT           |    0 | GFS_statein_type      |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | derived type GFS_stateout_type in FV3                   | DDT           |    0 | GFS_stateout_type     |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                   | derived type GFS_sfcprop_type in FV3                    | DDT           |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Coupling       | FV3-GFS_Coupling_type                                  | derived type GFS_coupling_type in FV3                   | DDT           |    0 | GFS_coupling_type     |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | derived type GFS_grid_type in FV3                       | DDT           |    0 | GFS_grid_type         |           | in     | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                       | derived type GFS_tbd_type in FV3                        | DDT           |    0 | GFS_tbd_type          |           | in     | F        |
!! | Cldprop        | FV3-GFS_Cldprop_type                                   | derived type GFS_cldprop_type in FV3                    | DDT           |    0 | GFS_cldprop_type      |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                   | derived type GFS_radtend_type in FV3                    | DDT           |    0 | GFS_radtend_type      |           | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                      | derived type GFS_diag_type in FV3                       | DDT           |    0 | GFS_diag_type         |           | in     | F        |
!! | Interstitial   | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | in     | F        |
!! | nthreads       | omp_threads                                            | number of OpenMP threads or fast physics schemes        | count         |    0 | integer               |           | in     | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_interstitialtoscreen_run (Model, Statein, Stateout, Sfcprop, Coupling, &
                                           Grid, Tbd, Cldprop, Radtend, Diag, Interstitial, &
                                           nthreads, errmsg, errflg)

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
         character(len=*),           intent(out) :: errmsg
         integer,                    intent(out) :: errflg

         !--- local variables
         integer :: impi, iomp, ierr
         integer :: mpirank, mpisize, mpicomm
         integer :: omprank, ompsize

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
         call MPI_BARRIER(mpicomm,ierr)
#endif

         do impi=0,mpisize-1
             do iomp=0,ompsize-1
                 if (mpirank==impi .and. omprank==iomp) then
                     call Interstitial%mprint(Model,mpirank,omprank,Tbd%blkno)
                 end if
#ifdef OPENMP
!$OMP BARRIER
#endif
             end do
#ifdef MPI
             call MPI_BARRIER(mpicomm,ierr)
#endif
         end do

#ifdef OPENMP
!$OMP BARRIER
#endif
#ifdef MPI
         call MPI_BARRIER(mpicomm,ierr)
#endif

      end subroutine GFS_interstitialtoscreen_run

    end module GFS_interstitialtoscreen

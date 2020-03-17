!> \file GFS_debug.F90

    module GFS_diagtoscreen

      private

      public GFS_diagtoscreen_init, GFS_diagtoscreen_run, GFS_diagtoscreen_finalize

      public print_my_stuff, chksum_int, chksum_real, print_var

! Calculating the checksum leads to segmentation faults with gfortran (bug in malloc?),
! thus print the sum of the array instead of the checksum.
#ifdef __GFORTRAN__
#define PRINT_SUM
#else
#define PRINT_CHKSUM
#endif

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
                     ! Sfcprop
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%slmsk'    , Sfcprop%slmsk)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%oceanfrac', Sfcprop%oceanfrac)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%landfrac' , Sfcprop%landfrac)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%lakefrac' , Sfcprop%lakefrac)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%tsfc'     , Sfcprop%tsfc)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%tsfco'    , Sfcprop%tsfco)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%tsfcl'    , Sfcprop%tsfcl)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%tisfc'    , Sfcprop%tisfc)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%snowd'    , Sfcprop%snowd)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%zorl'     , Sfcprop%zorl)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%zorlo'    , Sfcprop%zorlo)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%zorll'    , Sfcprop%zorll)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%fice'     , Sfcprop%fice)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%hprime'   , Sfcprop%hprime)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%sncovr'   , Sfcprop%sncovr)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%snoalb'   , Sfcprop%snoalb)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%alvsf'    , Sfcprop%alvsf)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%alnsf'    , Sfcprop%alnsf)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%alvwf'    , Sfcprop%alvwf)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%alnwf'    , Sfcprop%alnwf)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%facsf'    , Sfcprop%facsf)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%facwf'    , Sfcprop%facwf)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%slope'    , Sfcprop%slope)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%shdmin'   , Sfcprop%shdmin)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%shdmax'   , Sfcprop%shdmax)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%tg3'      , Sfcprop%tg3)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%vfrac'    , Sfcprop%vfrac)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%vtype'    , Sfcprop%vtype)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%stype'    , Sfcprop%stype)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%uustar'   , Sfcprop%uustar)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%oro'      , Sfcprop%oro)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%oro_uf'   , Sfcprop%oro_uf)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%hice'     , Sfcprop%hice)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%weasd'    , Sfcprop%weasd)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%canopy'   , Sfcprop%canopy)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%ffmm'     , Sfcprop%ffmm)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%ffhh'     , Sfcprop%ffhh)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%f10m'     , Sfcprop%f10m)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%tprcp'    , Sfcprop%tprcp)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%srflag'   , Sfcprop%srflag)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%slc'      , Sfcprop%slc)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%smc'      , Sfcprop%smc)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%stc'      , Sfcprop%stc)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%t2m'      , Sfcprop%t2m)
                     call print_var(mpirank,omprank, blkno, 'Sfcprop%q2m'      , Sfcprop%q2m)
                     if (Model%nstf_name(1)>0) then
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%tref    ', Sfcprop%tref)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%z_c     ', Sfcprop%z_c)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%c_0     ', Sfcprop%c_0)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%c_d     ', Sfcprop%c_d)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%w_0     ', Sfcprop%w_0)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%w_d     ', Sfcprop%w_d)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%xt      ', Sfcprop%xt)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%xs      ', Sfcprop%xs)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%xu      ', Sfcprop%xu)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%xv      ', Sfcprop%xv)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%xz      ', Sfcprop%xz)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%zm      ', Sfcprop%zm)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%xtts    ', Sfcprop%xtts)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%xzts    ', Sfcprop%xzts)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%d_conv  ', Sfcprop%d_conv)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%ifd     ', Sfcprop%ifd)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%dt_cool ', Sfcprop%dt_cool)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%qrain   ', Sfcprop%qrain)
                     end if
                     ! CCPP/RUC only
                     if (Model%lsm == Model%lsm_ruc) then
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%sh2o',         Sfcprop%sh2o)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%smois',        Sfcprop%smois)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%tslb',         Sfcprop%tslb)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%zs',           Sfcprop%zs)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%clw_surf',     Sfcprop%clw_surf)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%qwv_surf',     Sfcprop%qwv_surf)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%cndm_surf',    Sfcprop%cndm_surf)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%flag_frsoil',  Sfcprop%flag_frsoil)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%rhofr',        Sfcprop%rhofr)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%tsnow',        Sfcprop%tsnow)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%snowfallac  ', Sfcprop%snowfallac)
                        call print_var(mpirank,omprank, blkno, 'Sfcprop%acsnow      ', Sfcprop%acsnow)
                     end if
                     ! Radtend
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfcfsw%upfxc', Radtend%sfcfsw(:)%upfxc)
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfcfsw%dnfxc', Radtend%sfcfsw(:)%dnfxc)
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfcfsw%upfx0', Radtend%sfcfsw(:)%upfx0)
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfcfsw%dnfx0', Radtend%sfcfsw(:)%dnfx0)
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfcflw%upfxc', Radtend%sfcflw(:)%upfxc)
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfcflw%upfx0', Radtend%sfcflw(:)%upfx0)
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfcflw%dnfxc', Radtend%sfcflw(:)%dnfxc)
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfcflw%dnfx0', Radtend%sfcflw(:)%dnfx0)
                     call print_var(mpirank,omprank, blkno, 'Radtend%htrsw',        Radtend%htrsw)
                     call print_var(mpirank,omprank, blkno, 'Radtend%htrlw',        Radtend%htrlw)
                     call print_var(mpirank,omprank, blkno, 'Radtend%sfalb',        Radtend%sfalb)
                     call print_var(mpirank,omprank, blkno, 'Radtend%coszen',       Radtend%coszen)
                     call print_var(mpirank,omprank, blkno, 'Radtend%tsflw',        Radtend%tsflw)
                     call print_var(mpirank,omprank, blkno, 'Radtend%semis',        Radtend%semis)
                     call print_var(mpirank,omprank, blkno, 'Radtend%coszdg',       Radtend%coszdg)
                     call print_var(mpirank,omprank, blkno, 'Radtend%swhc',         Radtend%swhc)
                     call print_var(mpirank,omprank, blkno, 'Radtend%lwhc',         Radtend%lwhc)
                     call print_var(mpirank,omprank, blkno, 'Radtend%lwhd',         Radtend%lwhd)
                     ! Tbd
                     call print_var(mpirank,omprank, blkno, 'Tbd%icsdsw'          , Tbd%icsdsw)
                     call print_var(mpirank,omprank, blkno, 'Tbd%icsdlw'          , Tbd%icsdlw)
                     call print_var(mpirank,omprank, blkno, 'Tbd%ozpl'            , Tbd%ozpl)
                     call print_var(mpirank,omprank, blkno, 'Tbd%h2opl'           , Tbd%h2opl)
                     call print_var(mpirank,omprank, blkno, 'Tbd%rann'            , Tbd%rann)
                     call print_var(mpirank,omprank, blkno, 'Tbd%acv'             , Tbd%acv)
                     call print_var(mpirank,omprank, blkno, 'Tbd%acvb'            , Tbd%acvb)
                     call print_var(mpirank,omprank, blkno, 'Tbd%acvt'            , Tbd%acvt)
                     call print_var(mpirank,omprank, blkno, 'Tbd%hpbl'            , Tbd%hpbl)
                     if (Model%do_sppt) then
                       call print_var(mpirank,omprank, blkno, 'Tbd%dtdtr'         , Tbd%dtdtr)
                       call print_var(mpirank,omprank, blkno, 'Tbd%dtotprcp'      , Tbd%dtotprcp)
                       call print_var(mpirank,omprank, blkno, 'Tbd%dcnvprcp'      , Tbd%dcnvprcp)
                       call print_var(mpirank,omprank, blkno, 'Tbd%drain_cpl'     , Tbd%drain_cpl)
                       call print_var(mpirank,omprank, blkno, 'Tbd%dsnow_cpl'     , Tbd%dsnow_cpl)
                     end if
                     if (Model%nctp > 0 .and. Model%cscnv) then
                       call print_var(mpirank,omprank, blkno, 'Tbd%phy_fctd'      , Tbd%phy_fctd)
                     end if
                     call print_var(mpirank,omprank, blkno, 'Tbd%phy_f2d'         , Tbd%phy_f2d)
                     call print_var(mpirank,omprank, blkno, 'Tbd%phy_f3d'         , Tbd%phy_f3d)
                     do n=1,size(Tbd%phy_f3d(1,1,:))
                         call print_var(mpirank,omprank, blkno, 'Tbd%phy_f3d_n'   , Tbd%phy_f3d(:,:,n))
                     end do
                     call print_var(mpirank,omprank, blkno, 'Tbd%in_nm'           , Tbd%in_nm)
                     call print_var(mpirank,omprank, blkno, 'Tbd%ccn_nm'          , Tbd%ccn_nm)
                     call print_var(mpirank,omprank, blkno, 'Tbd%aer_nm'          , Tbd%aer_nm)
                     ! Diag
                     !call print_var(mpirank,omprank, blkno, 'Diag%fluxr       ',    Diag%fluxr)
                     !do n=1,size(Diag%fluxr(1,:))
                     !    call print_var(mpirank,omprank, blkno, 'Diag%fluxr_n ',    Diag%fluxr(:,n))
                     !end do
                     call print_var(mpirank,omprank, blkno, 'Diag%srunoff     ',    Diag%srunoff)
                     call print_var(mpirank,omprank, blkno, 'Diag%evbsa       ',    Diag%evbsa)
                     call print_var(mpirank,omprank, blkno, 'Diag%evcwa       ',    Diag%evcwa)
                     call print_var(mpirank,omprank, blkno, 'Diag%snohfa      ',    Diag%snohfa)
                     call print_var(mpirank,omprank, blkno, 'Diag%transa      ',    Diag%transa)
                     call print_var(mpirank,omprank, blkno, 'Diag%sbsnoa      ',    Diag%sbsnoa)
                     call print_var(mpirank,omprank, blkno, 'Diag%snowca      ',    Diag%snowca)
                     call print_var(mpirank,omprank, blkno, 'Diag%soilm       ',    Diag%soilm)
                     call print_var(mpirank,omprank, blkno, 'Diag%tmpmin      ',    Diag%tmpmin)
                     call print_var(mpirank,omprank, blkno, 'Diag%tmpmax      ',    Diag%tmpmax)
                     call print_var(mpirank,omprank, blkno, 'Diag%dusfc       ',    Diag%dusfc)
                     call print_var(mpirank,omprank, blkno, 'Diag%dvsfc       ',    Diag%dvsfc)
                     call print_var(mpirank,omprank, blkno, 'Diag%dtsfc       ',    Diag%dtsfc)
                     call print_var(mpirank,omprank, blkno, 'Diag%dqsfc       ',    Diag%dqsfc)
                     call print_var(mpirank,omprank, blkno, 'Diag%totprcp     ',    Diag%totprcp)
                     call print_var(mpirank,omprank, blkno, 'Diag%totice      ',    Diag%totice)
                     call print_var(mpirank,omprank, blkno, 'Diag%totsnw      ',    Diag%totsnw)
                     call print_var(mpirank,omprank, blkno, 'Diag%totgrp      ',    Diag%totgrp)
                     call print_var(mpirank,omprank, blkno, 'Diag%totprcpb    ',    Diag%totprcpb)
                     call print_var(mpirank,omprank, blkno, 'Diag%toticeb     ',    Diag%toticeb)
                     call print_var(mpirank,omprank, blkno, 'Diag%totsnwb     ',    Diag%totsnwb)
                     call print_var(mpirank,omprank, blkno, 'Diag%totgrpb     ',    Diag%totgrpb)
                     call print_var(mpirank,omprank, blkno, 'Diag%suntim      ',    Diag%suntim)
                     call print_var(mpirank,omprank, blkno, 'Diag%runoff      ',    Diag%runoff)
                     call print_var(mpirank,omprank, blkno, 'Diag%ep          ',    Diag%ep)
                     call print_var(mpirank,omprank, blkno, 'Diag%cldwrk      ',    Diag%cldwrk)
                     call print_var(mpirank,omprank, blkno, 'Diag%dugwd       ',    Diag%dugwd)
                     call print_var(mpirank,omprank, blkno, 'Diag%dvgwd       ',    Diag%dvgwd)
                     call print_var(mpirank,omprank, blkno, 'Diag%psmean      ',    Diag%psmean)
                     call print_var(mpirank,omprank, blkno, 'Diag%cnvprcp     ',    Diag%cnvprcp)
                     call print_var(mpirank,omprank, blkno, 'Diag%cnvprcpb    ',    Diag%cnvprcpb)
                     call print_var(mpirank,omprank, blkno, 'Diag%spfhmin     ',    Diag%spfhmin)
                     call print_var(mpirank,omprank, blkno, 'Diag%spfhmax     ',    Diag%spfhmax)
                     call print_var(mpirank,omprank, blkno, 'Diag%u10mmax     ',    Diag%u10mmax)
                     call print_var(mpirank,omprank, blkno, 'Diag%v10mmax     ',    Diag%v10mmax)
                     call print_var(mpirank,omprank, blkno, 'Diag%wind10mmax  ',    Diag%wind10mmax)
                     call print_var(mpirank,omprank, blkno, 'Diag%rain        ',    Diag%rain)
                     call print_var(mpirank,omprank, blkno, 'Diag%rainc       ',    Diag%rainc)
                     call print_var(mpirank,omprank, blkno, 'Diag%ice         ',    Diag%ice)
                     call print_var(mpirank,omprank, blkno, 'Diag%snow        ',    Diag%snow)
                     call print_var(mpirank,omprank, blkno, 'Diag%graupel     ',    Diag%graupel)
                     call print_var(mpirank,omprank, blkno, 'Diag%u10m        ',    Diag%u10m)
                     call print_var(mpirank,omprank, blkno, 'Diag%v10m        ',    Diag%v10m)
                     call print_var(mpirank,omprank, blkno, 'Diag%dpt2m       ',    Diag%dpt2m)
                     call print_var(mpirank,omprank, blkno, 'Diag%zlvl        ',    Diag%zlvl)
                     call print_var(mpirank,omprank, blkno, 'Diag%psurf       ',    Diag%psurf)
                     call print_var(mpirank,omprank, blkno, 'Diag%pwat        ',    Diag%pwat)
                     call print_var(mpirank,omprank, blkno, 'Diag%t1          ',    Diag%t1)
                     call print_var(mpirank,omprank, blkno, 'Diag%q1          ',    Diag%q1)
                     call print_var(mpirank,omprank, blkno, 'Diag%u1          ',    Diag%u1)
                     call print_var(mpirank,omprank, blkno, 'Diag%v1          ',    Diag%v1)
                     call print_var(mpirank,omprank, blkno, 'Diag%chh         ',    Diag%chh)
                     call print_var(mpirank,omprank, blkno, 'Diag%cmm         ',    Diag%cmm)
                     call print_var(mpirank,omprank, blkno, 'Diag%epi         ',    Diag%epi)
                     call print_var(mpirank,omprank, blkno, 'Diag%smcwlt2     ',    Diag%smcwlt2)
                     call print_var(mpirank,omprank, blkno, 'Diag%smcref2     ',    Diag%smcref2)
                     call print_var(mpirank,omprank, blkno, 'Diag%sr          ',    Diag%sr)
                     call print_var(mpirank,omprank, blkno, 'Diag%tdomr       ',    Diag%tdomr)
                     call print_var(mpirank,omprank, blkno, 'Diag%tdomzr      ',    Diag%tdomzr)
                     call print_var(mpirank,omprank, blkno, 'Diag%tdomip      ',    Diag%tdomip)
                     call print_var(mpirank,omprank, blkno, 'Diag%tdoms       ',    Diag%tdoms)
                     ! CCPP/RUC only
                     if (Model%lsm == Model%lsm_ruc) then
                       call print_var(mpirank,omprank, blkno, 'Diag%wet1        ',  Sfcprop%wetness)
                     else
                       call print_var(mpirank,omprank, blkno, 'Diag%wet1        ',  Diag%wet1)
                     end if
                     call print_var(mpirank,omprank, blkno, 'Diag%skebu_wts   ',    Diag%skebu_wts)
                     call print_var(mpirank,omprank, blkno, 'Diag%skebv_wts   ',    Diag%skebv_wts)
                     call print_var(mpirank,omprank, blkno, 'Diag%sppt_wts    ',    Diag%sppt_wts)
                     call print_var(mpirank,omprank, blkno, 'Diag%shum_wts    ',    Diag%shum_wts)
                     call print_var(mpirank,omprank, blkno, 'Diag%zmtnblck    ',    Diag%zmtnblck)
                     if (Model%ldiag3d) then
                       call print_var(mpirank,omprank, blkno, 'Diag%du3dt       ',    Diag%du3dt)
                       do n=1,size(Diag%du3dt(1,1,:))
                         call print_var(mpirank,omprank, blkno, 'Diag%du3dt_n     ',  Diag%du3dt(:,:,n))
                       end do
                       call print_var(mpirank,omprank, blkno, 'Diag%dv3dt       ',    Diag%dv3dt)
                       do n=1,size(Diag%dv3dt(1,1,:))
                         call print_var(mpirank,omprank, blkno, 'Diag%dv3dt_n     ',  Diag%dv3dt(:,:,n))
                       end do
                       call print_var(mpirank,omprank, blkno, 'Diag%dt3dt       ',    Diag%dt3dt)
                       do n=1,size(Diag%dt3dt(1,1,:))
                         call print_var(mpirank,omprank, blkno, 'Diag%dt3dt_n     ',  Diag%dt3dt(:,:,n))
                       end do
                       call print_var(mpirank,omprank, blkno, 'Diag%dq3dt       ',    Diag%dq3dt)
                       do n=1,size(Diag%dq3dt(1,1,:))
                         call print_var(mpirank,omprank, blkno, 'Diag%dq3dt_n     ',  Diag%dq3dt(:,:,n))
                       end do
                       call print_var(mpirank,omprank, blkno, 'Diag%upd_mf      ',    Diag%upd_mf)
                       call print_var(mpirank,omprank, blkno, 'Diag%dwn_mf      ',    Diag%dwn_mf)
                       call print_var(mpirank,omprank, blkno, 'Diag%det_mf      ',    Diag%det_mf)
                       call print_var(mpirank,omprank, blkno, 'Diag%cldcov      ',    Diag%cldcov)
                     end if
                     if(Model%lradar) then
                       call print_var(mpirank,omprank, blkno, 'Diag%refl_10cm   ',  Diag%refl_10cm)
                     end if
                     ! CCPP/MYNNPBL only
                     if (Model%do_mynnedmf) then
                       call print_var(mpirank,omprank, blkno, 'Diag%edmf_a      ',  Diag%edmf_a)
                       call print_var(mpirank,omprank, blkno, 'Diag%edmf_w      ',  Diag%edmf_w)
                       call print_var(mpirank,omprank, blkno, 'Diag%edmf_qt     ',  Diag%edmf_qt)
                       call print_var(mpirank,omprank, blkno, 'Diag%edmf_thl    ',  Diag%edmf_thl)
                       call print_var(mpirank,omprank, blkno, 'Diag%edmf_ent    ',  Diag%edmf_ent)
                       call print_var(mpirank,omprank, blkno, 'Diag%edmf_qc     ',  Diag%edmf_qc)
                       call print_var(mpirank,omprank, blkno, 'Diag%nupdraft    ',  Diag%nupdraft)
                       call print_var(mpirank,omprank, blkno, 'Diag%maxMF       ',  Diag%maxMF)
                       call print_var(mpirank,omprank, blkno, 'Diag%ktop_shallow',  Diag%ktop_shallow)
                       call print_var(mpirank,omprank, blkno, 'Diag%exch_h      ',  Diag%exch_h)
                       call print_var(mpirank,omprank, blkno, 'Diag%exch_m      ',  Diag%exch_m)
                     end if
                     ! Statein
                     call print_var(mpirank,omprank, blkno, 'Statein%phii'    ,     Statein%phii)
                     call print_var(mpirank,omprank, blkno, 'Statein%prsi'    ,     Statein%prsi)
                     call print_var(mpirank,omprank, blkno, 'Statein%prsik'   ,     Statein%prsik)
                     call print_var(mpirank,omprank, blkno, 'Statein%phil'    ,     Statein%phil)
                     call print_var(mpirank,omprank, blkno, 'Statein%prsl'    ,     Statein%prsl)
                     call print_var(mpirank,omprank, blkno, 'Statein%prslk'   ,     Statein%prslk)
                     call print_var(mpirank,omprank, blkno, 'Statein%pgr'     ,     Statein%pgr)
                     call print_var(mpirank,omprank, blkno, 'Statein%ugrs'    ,     Statein%ugrs)
                     call print_var(mpirank,omprank, blkno, 'Statein%vgrs'    ,     Statein%vgrs)
                     call print_var(mpirank,omprank, blkno, 'Statein%vvl'     ,     Statein%vvl)
                     call print_var(mpirank,omprank, blkno, 'Statein%tgrs'    ,     Statein%tgrs)
                     call print_var(mpirank,omprank, blkno, 'Statein%qgrs'    ,     Statein%qgrs)
                     do n=1,size(Statein%qgrs(1,1,:))
                        call print_var(mpirank,omprank, blkno, 'Statein%qgrs_n',    Statein%qgrs(:,:,n))
                     end do
                     call print_var(mpirank,omprank, blkno, 'Statein%diss_est',     Statein%diss_est)
                     call print_var(mpirank,omprank, blkno, 'Statein%smc'     ,     Statein%smc)
                     call print_var(mpirank,omprank, blkno, 'Statein%stc'     ,     Statein%stc)
                     call print_var(mpirank,omprank, blkno, 'Statein%slc'     ,     Statein%slc)
                     ! Stateout
                     call print_var(mpirank,omprank, blkno, 'Stateout%gu0',         Stateout%gu0)
                     call print_var(mpirank,omprank, blkno, 'Stateout%gv0',         Stateout%gv0)
                     call print_var(mpirank,omprank, blkno, 'Stateout%gt0',         Stateout%gt0)
                     call print_var(mpirank,omprank, blkno, 'Stateout%gq0',         Stateout%gq0)
                     do n=1,size(Stateout%gq0(1,1,:))
                        call print_var(mpirank,omprank, blkno, 'Stateout%gq0_n',    Stateout%gq0(:,:,n))
                     end do
                     ! Coupling
                     call print_var(mpirank,omprank, blkno, 'Coupling%nirbmdi', Coupling%nirbmdi)
                     call print_var(mpirank,omprank, blkno, 'Coupling%nirdfdi', Coupling%nirdfdi)
                     call print_var(mpirank,omprank, blkno, 'Coupling%visbmdi', Coupling%visbmdi)
                     call print_var(mpirank,omprank, blkno, 'Coupling%visdfdi', Coupling%visdfdi)
                     call print_var(mpirank,omprank, blkno, 'Coupling%nirbmui', Coupling%nirbmui)
                     call print_var(mpirank,omprank, blkno, 'Coupling%nirdfui', Coupling%nirdfui)
                     call print_var(mpirank,omprank, blkno, 'Coupling%visbmui', Coupling%visbmui)
                     call print_var(mpirank,omprank, blkno, 'Coupling%visdfui', Coupling%visdfui)
                     call print_var(mpirank,omprank, blkno, 'Coupling%sfcdsw ', Coupling%sfcdsw )
                     call print_var(mpirank,omprank, blkno, 'Coupling%sfcnsw ', Coupling%sfcnsw )
                     call print_var(mpirank,omprank, blkno, 'Coupling%sfcdlw ', Coupling%sfcdlw )
                     if (Model%cplflx .or. Model%do_sppt .or. Model%cplchm) then
                        call print_var(mpirank,omprank, blkno, 'Coupling%rain_cpl', Coupling%rain_cpl)
                        call print_var(mpirank,omprank, blkno, 'Coupling%snow_cpl', Coupling%snow_cpl)
                     end if
                     if (Model%cplflx) then
                        call print_var(mpirank,omprank, blkno, 'Coupling%slimskin_cpl', Coupling%slimskin_cpl )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dusfcin_cpl ', Coupling%dusfcin_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dvsfcin_cpl ', Coupling%dvsfcin_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dtsfcin_cpl ', Coupling%dtsfcin_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dqsfcin_cpl ', Coupling%dqsfcin_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%ulwsfcin_cpl', Coupling%ulwsfcin_cpl )
                        call print_var(mpirank,omprank, blkno, 'Coupling%tseain_cpl  ', Coupling%tseain_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%tisfcin_cpl ', Coupling%tisfcin_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%ficein_cpl  ', Coupling%ficein_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%hicein_cpl  ', Coupling%hicein_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%hsnoin_cpl  ', Coupling%hsnoin_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dusfc_cpl   ', Coupling%dusfc_cpl    )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dvsfc_cpl   ', Coupling%dvsfc_cpl    )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dtsfc_cpl   ', Coupling%dtsfc_cpl    )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dqsfc_cpl   ', Coupling%dqsfc_cpl    )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dlwsfc_cpl  ', Coupling%dlwsfc_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dswsfc_cpl  ', Coupling%dswsfc_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dnirbm_cpl  ', Coupling%dnirbm_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dnirdf_cpl  ', Coupling%dnirdf_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dvisbm_cpl  ', Coupling%dvisbm_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dvisdf_cpl  ', Coupling%dvisdf_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nlwsfc_cpl  ', Coupling%nlwsfc_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nswsfc_cpl  ', Coupling%nswsfc_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nnirbm_cpl  ', Coupling%nnirbm_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nnirdf_cpl  ', Coupling%nnirdf_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nvisbm_cpl  ', Coupling%nvisbm_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nvisdf_cpl  ', Coupling%nvisdf_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dusfci_cpl  ', Coupling%dusfci_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dvsfci_cpl  ', Coupling%dvsfci_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dtsfci_cpl  ', Coupling%dtsfci_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dqsfci_cpl  ', Coupling%dqsfci_cpl   )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dlwsfci_cpl ', Coupling%dlwsfci_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dswsfci_cpl ', Coupling%dswsfci_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dnirbmi_cpl ', Coupling%dnirbmi_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dnirdfi_cpl ', Coupling%dnirdfi_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dvisbmi_cpl ', Coupling%dvisbmi_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dvisdfi_cpl ', Coupling%dvisdfi_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nlwsfci_cpl ', Coupling%nlwsfci_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nswsfci_cpl ', Coupling%nswsfci_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nnirbmi_cpl ', Coupling%nnirbmi_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nnirdfi_cpl ', Coupling%nnirdfi_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nvisbmi_cpl ', Coupling%nvisbmi_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%nvisdfi_cpl ', Coupling%nvisdfi_cpl  )
                        call print_var(mpirank,omprank, blkno, 'Coupling%t2mi_cpl    ', Coupling%t2mi_cpl     )
                        call print_var(mpirank,omprank, blkno, 'Coupling%q2mi_cpl    ', Coupling%q2mi_cpl     )
                        call print_var(mpirank,omprank, blkno, 'Coupling%u10mi_cpl   ', Coupling%u10mi_cpl    )
                        call print_var(mpirank,omprank, blkno, 'Coupling%v10mi_cpl   ', Coupling%v10mi_cpl    )
                        call print_var(mpirank,omprank, blkno, 'Coupling%tsfci_cpl   ', Coupling%tsfci_cpl    )
                        call print_var(mpirank,omprank, blkno, 'Coupling%psurfi_cpl  ', Coupling%psurfi_cpl   )
                     end if
                     if (Model%cplchm) then
                        call print_var(mpirank,omprank, blkno, 'Coupling%rainc_cpl', Coupling%rainc_cpl)
                        call print_var(mpirank,omprank, blkno, 'Coupling%ushfsfci ', Coupling%ushfsfci )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dkt      ', Coupling%dkt      )
                        call print_var(mpirank,omprank, blkno, 'Coupling%dqdti    ', Coupling%dqdti    )
                     end if
                     if (Model%do_sppt) then
                        call print_var(mpirank,omprank, blkno, 'Coupling%sppt_wts', Coupling%sppt_wts)
                     end if
                     if (Model%do_shum) then
                        call print_var(mpirank,omprank, blkno, 'Coupling%shum_wts', Coupling%shum_wts)
                     end if
                     if (Model%do_skeb) then
                        call print_var(mpirank,omprank, blkno, 'Coupling%skebu_wts', Coupling%skebu_wts)
                        call print_var(mpirank,omprank, blkno, 'Coupling%skebv_wts', Coupling%skebv_wts)
                     end if
                     if (Model%do_sfcperts) then
                        call print_var(mpirank,omprank, blkno, 'Coupling%sfc_wts', Coupling%sfc_wts)
                     end if
                     if(Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
                        call print_var(mpirank,omprank, blkno, 'Coupling%nwfa2d', Coupling%nwfa2d)
                        call print_var(mpirank,omprank, blkno, 'Coupling%nifa2d', Coupling%nifa2d)
                     end if
                     ! Grid
                     call print_var(mpirank,omprank, blkno, 'Grid%xlon  ', Grid%xlon  )
                     call print_var(mpirank,omprank, blkno, 'Grid%xlat  ', Grid%xlat  )
                     call print_var(mpirank,omprank, blkno, 'Grid%xlat_d', Grid%xlat_d)
                     call print_var(mpirank,omprank, blkno, 'Grid%sinlat', Grid%sinlat)
                     call print_var(mpirank,omprank, blkno, 'Grid%coslat', Grid%coslat)
                     call print_var(mpirank,omprank, blkno, 'Grid%area  ', Grid%area  )
                     call print_var(mpirank,omprank, blkno, 'Grid%dx    ', Grid%dx    )
                     if (Model%ntoz > 0) then
                        call print_var(mpirank,omprank, blkno, 'Grid%ddy_o3   ', Grid%ddy_o3   )
                        call print_var(mpirank,omprank, blkno, 'Grid%jindx1_o3', Grid%jindx1_o3)
                        call print_var(mpirank,omprank, blkno, 'Grid%jindx2_o3', Grid%jindx2_o3)
                     endif
                     if (Model%h2o_phys) then
                        call print_var(mpirank,omprank, blkno, 'Grid%ddy_h   ', Grid%ddy_h   )
                        call print_var(mpirank,omprank, blkno, 'Grid%jindx1_h', Grid%jindx1_h)
                        call print_var(mpirank,omprank, blkno, 'Grid%jindx2_h', Grid%jindx2_h)
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
                     call Interstitial%mprint(Model,mpirank,omprank,blkno)
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

         if (Model%kdt==1 .and. blkno==4) then
             if (Model%me==0) write(0,*) "GFS_abort_run: ABORTING MODEL"
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

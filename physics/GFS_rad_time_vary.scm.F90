!>\file GFS_rad_time_vary.F90
!!  Contains code related to GFS physics suite setup (radiation part of time_vary_step)
      module GFS_rad_time_vary

      implicit none

      private

      public GFS_rad_time_vary_init, GFS_rad_time_vary_run, GFS_rad_time_vary_finalize

      contains

!>\defgroup GFS_rad_time_vary GFS RRTMG Update
!!\ingroup RRTMG
!! @{
      subroutine GFS_rad_time_vary_init
      end subroutine GFS_rad_time_vary_init

!> \section arg_table_GFS_rad_time_vary_run Argument Table
!! \htmlinclude GFS_rad_time_vary_run.html
!!
      subroutine GFS_rad_time_vary_run (Model, Statein, Tbd, errmsg, errflg)

      use physparam,                 only: ipsd0, ipsdlim, iaerflg
      use mersenne_twister,          only: random_setseed, random_index, random_stat
      use machine,                   only: kind_phys
      use GFS_typedefs,              only: GFS_statein_type,   &
                                           GFS_control_type,   &
                                           GFS_grid_type,      &
                                           GFS_tbd_type
      use radcons,                   only: qmin, con_100

      implicit none

      type(GFS_control_type), intent(inout) :: Model
      type(GFS_statein_type), intent(in)    :: Statein
      type(GFS_tbd_type),     intent(inout) :: Tbd
      character(len=*),       intent(out) :: errmsg
      integer,                intent(out) :: errflg

      !--- local variables
      type (random_stat) :: stat
      integer :: ix, nb, j, i, nblks, ipseed
      integer :: numrdm(Model%cnx*Model%cny*2)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      nb = 1

      if (Model%lsswr .or. Model%lslwr) then

        !--- call to GFS_radupdate_run is now in GFS_rrtmg_setup_run

        !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
        if ((Model%isubc_lw==2) .or. (Model%isubc_sw==2)) then
          ipseed = mod(nint(con_100*sqrt(Model%sec)), ipsdlim) + 1 + ipsd0
          call random_setseed (ipseed, stat)
          call random_index (ipsdlim, numrdm, stat)

          !--- set the random seeds for each column in a reproducible way
          do ix=1,Model%blksz(nb)
             j = Tbd%jmap(ix)
             i = Tbd%imap(ix)
             !--- for testing purposes, replace numrdm with '100'
             Tbd%icsdsw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx)
             Tbd%icsdlw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + Model%cnx*Model%cny)
          enddo
        endif  ! isubc_lw and isubc_sw

        if (Model%imp_physics == 99) then
          if (Model%kdt == 1) then
            Tbd%phy_f3d(:,:,1) = Statein%tgrs
            Tbd%phy_f3d(:,:,2) = max(qmin,Statein%qgrs(:,:,1))
            Tbd%phy_f3d(:,:,3) = Statein%tgrs
            Tbd%phy_f3d(:,:,4) = max(qmin,Statein%qgrs(:,:,1))
            Tbd%phy_f2d(:,1)   = Statein%prsi(:,1)
            Tbd%phy_f2d(:,2)   = Statein%prsi(:,1)
          endif
        endif

      endif

  end subroutine GFS_rad_time_vary_run

  subroutine GFS_rad_time_vary_finalize()
  end subroutine GFS_rad_time_vary_finalize
!! @}
  end module GFS_rad_time_vary

!>\file GFS_rad_time_vary.f90
!! This file contains
      module GFS_rad_time_vary
      contains

!>\defgroup GFS_rad_time_vary GFS RRTMG Update 
!! @{
!! \section arg_table_GFS_rad_time_vary_init Argument Table
!!
      subroutine GFS_rad_time_vary_init
      end subroutine GFS_rad_time_vary_init

!> \section arg_table_GFS_rad_time_vary_run Argument Table
!! | local var name    | longname                                                      | description                                                                   | units    | rank |  type                         |   kind    | intent | optional |
!! |-------------------|---------------------------------------------------------------|-------------------------------------------------------------------------------|----------|------|-------------------------------|-----------|--------|----------|
!! |   Model           | FV3-GFS_Control_type                                          | Fortran DDT containing FV3-GFS model control parameters                       | DDT      |  0   | GFS_typedefs%GFS_control_type |           | in     | F        |
!! |   Statein         | FV3-GFS_Statein_type                                          | Fortran DDT containing FV3-GFS prognostic state data in from dycore           | DDT      |  0   | GFS_typedefs%GFS_statein_type |           | in     | F        |
!! |   Tbd             | FV3-GFS_Tbd_type                                              | Fortran DDT containing FV3-GFS data not yet assigned to a defined container   | DDT      |  0   | GFS_typedefs%GFS_tbd_type     |           | in     | F        |
!! |   blksz           | horizontal_block_size                                         | horizontal block size for explicit data blocking                              | none     |  1   | integer                       |           | in     | F        |
!! |   sec             | seconds_elapsed_since_model_initialization                    | seconds elapsed since model initialization                                    | s        |  0   | real                          | kind_phys | in     | F        |
!! |   ictmflg         | flag_for_initial_time-date_control                            | flag for initial time/date control                                            | none     |  0   | integer                       |           | in     | F        |
!! |   isolar          | flag_for_solar_constant                                       | solar constant control flag                                                   | none     |  0   | integer                       |           | in     | F        |
!!
      subroutine GFS_rad_time_vary_run (Model, Statein, Tbd, blksz, sec, ictmflg, isolar)

      use physparam,                 only: ipsd0, ipsdlim, iaerflg
      use mersenne_twister,          only: random_setseed, random_index, random_stat
      use machine,                   only: kind_phys
      use GFS_typedefs,              only: GFS_statein_type,   &
                                           GFS_control_type,   &
                                           GFS_grid_type,      &
                                           GFS_tbd_type     
      !use module_radiation_driver,   only: radupdate
      use GFS_radupdate,              only: GFS_radupdate_run
      use radcons,                    only: qmin, con_100
      


      implicit none

      type(GFS_control_type), intent(inout) :: Model
      type(GFS_statein_type), intent(in)    :: Statein(:)
      type(GFS_tbd_type),     intent(inout) :: Tbd(:)
      real(kind=kind_phys),   intent(in)    :: sec
      integer, intent(in)  :: ictmflg, isolar
      integer, allocatable :: blksz(:)

      !--- local variables
      type (random_stat) :: stat
      integer :: ix, nb, j, i, nblks, ipseed
      integer :: numrdm(Model%cnx*Model%cny*2)

      
     if (Model%lsswr .or. Model%lslwr) then

      nblks = size(blksz,1)
!
      call GFS_radupdate_run (Model%idat, Model%jdat, Model%fhswr, Model%dtf, Model%lsswr, &
                    Model%me, Model%slag, Model%sdec, Model%cdec, Model%solcon,            &
                    ictmflg, isolar )

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

    endif

  end subroutine GFS_rad_time_vary_run
 
!> \section arg_table_GFS_rad_time_vary_finalize Argument Table
!!
  subroutine GFS_rad_time_vary_finalize()
  end subroutine GFS_rad_time_vary_finalize
!! @}
  end module GFS_rad_time_vary

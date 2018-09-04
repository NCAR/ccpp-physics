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
!! \section arg_table_GFS_rad_time_vary_init Argument Table
!!
      subroutine GFS_rad_time_vary_init
      end subroutine GFS_rad_time_vary_init

!> \section arg_table_GFS_rad_time_vary_run Argument Table
!! | local_name        | standard_name                                          | long_name                                                                     | units    | rank |  type                 |   kind    | intent | optional |
!! |-------------------|--------------------------------------------------------|-------------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model             | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                       | DDT      |    0 | GFS_control_type      |           | inout  | F        |
!! | Data              | FV3-GFS_Data_type_all_blocks                           | Fortran DDT containing FV3-GFS data                                           | DDT      |    1 | GFS_data_type         |           | inout  | F        |
!! | errmsg            | ccpp_error_message                                     | error message for error handling in CCPP                                      | none     |    0 | character             | len=*     | out    | F        |
!! | errflg            | ccpp_error_flag                                        | error flag for error handling in CCPP                                         | flag     |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_rad_time_vary_run (Model, Data, errmsg, errflg)

         use physparam,                 only: ipsd0, ipsdlim, iaerflg
         use mersenne_twister,          only: random_setseed, random_index, random_stat
         use machine,                   only: kind_phys
         use GFS_typedefs,              only: GFS_control_type, &
                                              GFS_data_type
         use radcons,                   only: qmin, con_100

         implicit none

         type(GFS_control_type), intent(inout) :: Model
         type(GFS_data_type),    intent(inout) :: Data(:)
         character(len=*),       intent(out) :: errmsg
         integer,                intent(out) :: errflg

         !--- local variables
         type (random_stat) :: stat
         integer :: ix, nb, j, i, nblks, ipseed
         integer :: numrdm(Model%cnx*Model%cny*2)

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (Model%lsswr .or. Model%lslwr) then

           nblks = size(Model%blksz)

           !--- call to GFS_radupdate_run is now in GFS_rrtmg_setup_run

           !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
           if ((Model%isubc_lw==2) .or. (Model%isubc_sw==2)) then
             ipseed = mod(nint(con_100*sqrt(Model%sec)), ipsdlim) + 1 + ipsd0
             call random_setseed (ipseed, stat)
             call random_index (ipsdlim, numrdm, stat)
    
             !--- set the random seeds for each column in a reproducible way
             ix = 0
             nb = 1
             ! DH* TODO - this could be sped up by saving jsc, jec, isc, iec in Tbd (for example)
             ! and looping just over them; ix would then run from 1 to blksz(nb); one could also
             ! use OpenMP to speed up this loop *DH
             do j = 1,Model%ny
               do i = 1,Model%nx
                 ix = ix + 1
                 if (ix .gt. Model%blksz(nb)) then
                   ix = 1
                   nb = nb + 1
                 endif
              
                 !--- for testing purposes, replace numrdm with '100'
                 Data(nb)%Tbd%icsdsw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx)
                 Data(nb)%Tbd%icsdlw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + Model%cnx*Model%cny)
               enddo
             enddo
           endif  ! isubc_lw and isubc_sw

           if (Model%imp_physics == 99) then
             if (Model%kdt == 1) then
   ! DH* OpenMP?
               do nb = 1,nblks
                 Data(nb)%Tbd%phy_f3d(:,:,1) = Data(nb)%Statein%tgrs
                 Data(nb)%Tbd%phy_f3d(:,:,2) = max(qmin,Data(nb)%Statein%qgrs(:,:,1))
                 Data(nb)%Tbd%phy_f3d(:,:,3) = Data(nb)%Statein%tgrs
                 Data(nb)%Tbd%phy_f3d(:,:,4) = max(qmin,Data(nb)%Statein%qgrs(:,:,1))
                 Data(nb)%Tbd%phy_f2d(:,1)   = Data(nb)%Statein%prsi(:,1)
                 Data(nb)%Tbd%phy_f2d(:,2)   = Data(nb)%Statein%prsi(:,1)
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

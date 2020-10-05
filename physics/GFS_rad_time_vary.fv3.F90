!>\file GFS_rad_time_vary.F90
!!  Contains code related to GFS physics suite setup (radiation part of time_vary_step)
   module GFS_rad_time_vary

      implicit none

      private

      public GFS_rad_time_vary_init, GFS_rad_time_vary_run, GFS_rad_time_vary_finalize

      contains

      subroutine GFS_rad_time_vary_init
      end subroutine GFS_rad_time_vary_init

!>\defgroup mod_GFS_rad_time_vary GFS Radiation Time Update
!> @{
!> \section arg_table_GFS_rad_time_vary_run Argument Table
!! \htmlinclude GFS_rad_time_vary_run.html
!!
      subroutine GFS_rad_time_vary_run (Model, Data, nthrds, errmsg, errflg)

         use physparam,                 only: ipsd0, ipsdlim, iaerflg
         use mersenne_twister,          only: random_setseed, random_index, random_stat
         use machine,                   only: kind_phys
         use GFS_typedefs,              only: GFS_control_type, &
                                              GFS_data_type
         use radcons,                   only: qmin, con_100

         implicit none

         type(GFS_control_type), intent(inout) :: Model
         type(GFS_data_type),    intent(inout) :: Data(:)
         integer,                intent(in)    :: nthrds
         character(len=*),       intent(out)   :: errmsg
         integer,                intent(out)   :: errflg

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

!$OMP parallel num_threads(nthrds) default(none)        &
!$OMP          private (nb,ix,i,j)                      &
!$OMP          shared (Model,Data,ipsdlim,ipsd0,ipseed) &
!$OMP          shared (numrdm,stat,nblks)

           !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
           if ((Model%isubc_lw==2) .or. (Model%isubc_sw==2)) then
!$OMP single
             ipseed = mod(nint(con_100*sqrt(Model%sec)), ipsdlim) + 1 + ipsd0
             call random_setseed (ipseed, stat)
             call random_index (ipsdlim, numrdm, stat)
!$OMP end single

!$OMP do schedule (dynamic,1)
             do nb=1,nblks
               do ix=1,Model%blksz(nb)
                 j = Data(nb)%Tbd%jmap(ix)
                 i = Data(nb)%Tbd%imap(ix)
                 !--- for testing purposes, replace numrdm with '100'
                 Data(nb)%Tbd%icsdsw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx)
                 Data(nb)%Tbd%icsdlw(ix) = numrdm(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + Model%cnx*Model%cny)
               enddo
             enddo
!$OMP end do
           endif  ! isubc_lw and isubc_sw

           if (Model%imp_physics == 99) then
             if (Model%kdt == 1) then
!$OMP do schedule (dynamic,1)
               do nb = 1,nblks
                 Data(nb)%Tbd%phy_f3d(:,:,1) = Data(nb)%Statein%tgrs
                 Data(nb)%Tbd%phy_f3d(:,:,2) = max(qmin,Data(nb)%Statein%qgrs(:,:,1))
                 Data(nb)%Tbd%phy_f3d(:,:,3) = Data(nb)%Statein%tgrs
                 Data(nb)%Tbd%phy_f3d(:,:,4) = max(qmin,Data(nb)%Statein%qgrs(:,:,1))
                 Data(nb)%Tbd%phy_f2d(:,1)   = Data(nb)%Statein%prsi(:,1)
                 Data(nb)%Tbd%phy_f2d(:,2)   = Data(nb)%Statein%prsi(:,1)
               enddo
!$OMP end do
             endif
           endif

!$OMP end parallel

         endif

      end subroutine GFS_rad_time_vary_run
!> @}
 
      subroutine GFS_rad_time_vary_finalize()
      end subroutine GFS_rad_time_vary_finalize

   end module GFS_rad_time_vary

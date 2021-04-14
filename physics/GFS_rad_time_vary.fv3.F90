!>\file GFS_rad_time_vary.F90
!!  Contains code related to GFS physics suite setup (radiation part of time_vary_step)
   module GFS_rad_time_vary

      implicit none

      private

      public GFS_rad_time_vary_timestep_init

      contains

!>\defgroup mod_GFS_rad_time_vary GFS Radiation Time Update
!> @{
!> \section arg_table_GFS_rad_time_vary_timestep_init Argument Table
!! \htmlinclude GFS_rad_time_vary_timestep_init.html
!!
      subroutine GFS_rad_time_vary_timestep_init (                                     &
              lslwr, lsswr, isubc_lw, isubc_sw, icsdsw, icsdlw, cnx, cny, isc, jsc,    &
              imap, jmap, sec, kdt, imp_physics, imp_physics_zhao_carr, ps_2delt,      &
              ps_1delt, t_2delt, t_1delt, qv_2delt, qv_1delt, t, qv, ps, errmsg, errflg)

         use physparam,                 only: ipsd0, ipsdlim, iaerflg
         use mersenne_twister,          only: random_setseed, random_index, random_stat
         use machine,                   only: kind_phys
         use radcons,                   only: qmin, con_100

         implicit none

         ! Interface variables
         integer,                intent(in)    :: isubc_lw, isubc_sw, cnx, cny, isc, jsc, kdt
         integer,                intent(in)    :: imp_physics, imp_physics_zhao_carr
         logical,                intent(in)    :: lslwr, lsswr
         integer,                intent(inout) :: icsdsw(:), icsdlw(:)
         integer,                intent(in)    :: imap(:), jmap(:)
         real(kind_phys),        intent(in)    :: sec
         real(kind_phys),        intent(inout) :: ps_2delt(:)
         real(kind_phys),        intent(inout) :: ps_1delt(:)
         real(kind_phys),        intent(inout) :: t_2delt(:,:)
         real(kind_phys),        intent(inout) :: t_1delt(:,:)
         real(kind_phys),        intent(inout) :: qv_2delt(:,:)
         real(kind_phys),        intent(inout) :: qv_1delt(:,:)
         real(kind_phys),        intent(in)    :: t(:,:), qv(:,:), ps(:)
         character(len=*),       intent(out)   :: errmsg
         integer,                intent(out)   :: errflg

         ! Local variables
         type (random_stat) :: stat
         integer :: ix, j, i, nblks, ipseed
         integer :: numrdm(cnx*cny*2)

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (lsswr .or. lslwr) then

           !--- call to GFS_radupdate_timestep_init is now in GFS_rrtmg_setup_timestep_init

           !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
           if ((isubc_lw==2) .or. (isubc_sw==2)) then
             ipseed = mod(nint(con_100*sqrt(sec)), ipsdlim) + 1 + ipsd0
             call random_setseed (ipseed, stat)
             call random_index (ipsdlim, numrdm, stat)

             do ix=1,size(jmap)
               j = jmap(ix)
               i = imap(ix)
               !--- for testing purposes, replace numrdm with '100'
               icsdsw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx)
               icsdlw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx + cnx*cny)
             enddo

           endif  ! isubc_lw and isubc_sw

           if (imp_physics == imp_physics_zhao_carr) then
             if (kdt == 1) then
               t_2delt  = t
               t_1delt  = t
               qv_2delt = qv
               qv_1delt = qv
               ps_2delt = ps
               ps_1delt = ps
             endif
           endif

         endif

      end subroutine GFS_rad_time_vary_timestep_init
!> @}

   end module GFS_rad_time_vary

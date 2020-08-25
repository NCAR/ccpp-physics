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
!! \htmlinclude GFS_rad_time_vary_run.html
!!
      subroutine GFS_rad_time_vary_run (cnx, cny, lsswr, lslwr, isubc_sw, &
          isubc_lw, sec, nblks, blksz, isc, jsc, imp_physics,             &
          imp_physics_zhao_carr, kdt, tgrs, qgrs_wv, prsi, imap, jmap,    &
          icsdsw, icsdlw, t_minus_two_delt, qv_minus_two_delt,            &
          t_minus_delt, qv_minus_delt, ps_minus_two_delt, ps_minus_delt,  &
          errmsg, errflg)

      use physparam,                 only: ipsd0, ipsdlim, iaerflg
      use mersenne_twister,          only: random_setseed, random_index, random_stat
      use machine,                   only: kind_phys
      use radcons,                   only: qmin, con_100

      implicit none

      integer,                intent(in)    :: cnx, cny, isubc_sw, isubc_lw, &
                                                nblks, isc, jsc, imp_physics,&
                                                imp_physics_zhao_carr, kdt
      logical,                intent(in)    :: lsswr, lslwr
      real(kind=kind_phys),   intent(in)    :: sec
      
      integer, dimension(nblks), intent(in) :: blksz
      integer, dimension(:),     intent(in) :: imap, jmap
      
      integer, dimension(:),  intent(inout) :: icsdsw, icsdlw
      
      real(kind=kind_phys), dimension(:,:), intent(in) :: tgrs, qgrs_wv
      real(kind=kind_phys), dimension(:,:), intent(in) :: prsi
      
      real(kind=kind_phys), dimension(:,:), intent(inout) :: t_minus_two_delt, &
                                qv_minus_two_delt, t_minus_delt, qv_minus_delt
      real(kind=kind_phys), dimension(:),   intent(inout) :: ps_minus_two_delt,&
                                ps_minus_delt
      
      character(len=*),       intent(out) :: errmsg
      integer,                intent(out) :: errflg

      !--- local variables
      type (random_stat) :: stat
      integer :: ix, nb, j, i, ipseed
      integer :: numrdm(cnx*cny*2)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      nb = 1

      if (lsswr .or. lslwr) then

        !--- call to GFS_radupdate_run is now in GFS_rrtmg_setup_run

        !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
        if ((isubc_lw==2) .or. (isubc_sw==2)) then
          ipseed = mod(nint(con_100*sqrt(sec)), ipsdlim) + 1 + ipsd0
          call random_setseed (ipseed, stat)
          call random_index (ipsdlim, numrdm, stat)

          !--- set the random seeds for each column in a reproducible way
          do ix=1,blksz(nb)
             j = jmap(ix)
             i = imap(ix)
             !--- for testing purposes, replace numrdm with '100'
             icsdsw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx)
             icsdlw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx + cnx*cny)
          enddo
        endif  ! isubc_lw and isubc_sw

        if (imp_physics == imp_physics_zhao_carr) then
          if (kdt == 1) then
            t_minus_two_delt(:,:)  = tgrs
            qv_minus_two_delt(:,:) = max(qmin,qgrs_wv(:,:))
            t_minus_delt(:,:)      = tgrs
            qv_minus_delt(:,:)     = max(qmin,qgrs_wv(:,:))
            ps_minus_two_delt(:)   = prsi(:,1)
            ps_minus_delt(:)       = prsi(:,1)
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

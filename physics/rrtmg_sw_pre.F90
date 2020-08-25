!>\file rrtmg_sw_pre.f90
!! This file contains a subroutine to module_radiation_surface::setalb() to
!! setup surface albedo for SW radiation.
      module rrtmg_sw_pre
      contains

!>\defgroup rrtmg_sw_pre GFS RRTMG scheme Pre
!! @{
!> \section arg_table_rrtmg_sw_pre_init Argument Table
!!
      subroutine rrtmg_sw_pre_init ()
      end subroutine rrtmg_sw_pre_init

!> \section arg_table_rrtmg_sw_pre_run Argument Table
!! \htmlinclude rrtmg_sw_pre_run.html
!!
      subroutine rrtmg_sw_pre_run (Model, Grid, Sfcprop, Radtend, im, &
        nday, idxday, tsfg, tsfa, sfcalb1, sfcalb2, sfcalb3, sfcalb4, &
        alb1d, errmsg, errflg)

      use machine,                   only: kind_phys

      use GFS_typedefs,              only: GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_sfcprop_type
      use module_radiation_surface,  only: NF_ALBD, setalb

      implicit none

      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
      type(GFS_grid_type),            intent(in)    :: Grid
      integer,                        intent(in)    :: im
      integer,                        intent(out)   :: nday
      integer, dimension(size(Grid%xlon,1)), intent(out) :: idxday
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)  :: tsfa, tsfg
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(out) :: sfcalb1, sfcalb2, sfcalb3, sfcalb4
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)  :: alb1d
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer :: i
      real(kind=kind_phys), dimension(size(Grid%xlon,1),NF_ALBD) :: sfcalb

      real(kind=kind_phys) :: lndp_alb

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- ...  start radiation calculations
!           remember to set heating rate unit to k/sec!
!> -# Start SW radiation calculations
      if (Model%lsswr) then

!>  - Check for daytime points for SW radiation.
        nday = 0
        idxday = 0
        do i = 1, IM
          if (Radtend%coszen(i) >= 0.0001) then
            nday = nday + 1
            idxday(nday) = i
          endif
        enddo

! set albedo pert, if requested.
        lndp_alb = -999.
        if (Model%lndp_type==1) then
          do i =1,Model%n_var_lndp
            if (Model%lndp_var_list(i) == 'alb') then
                lndp_alb = Model%lndp_prt_list(i)
            endif
          enddo
        endif

!>  - Call module_radiation_surface::setalb() to setup surface albedo.
!!  for SW radiation.

        call setalb (Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,   &  !  ---  inputs:
                     Sfcprop%snoalb, Sfcprop%zorl, Radtend%coszen,   &
                     tsfg, tsfa, Sfcprop%hprime(:,1), Sfcprop%alvsf, &
                     Sfcprop%alnsf, Sfcprop%alvwf, Sfcprop%alnwf,    &
                     Sfcprop%facsf, Sfcprop%facwf, Sfcprop%fice,     &
                     Sfcprop%tisfc, IM,                              &
                     alb1d, lndp_alb,                           &  !  mg, sfc-perts
                     sfcalb)                                           !  ---  outputs

!> -# Approximate mean surface albedo from vis- and nir-  diffuse values.
        Radtend%sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
      else
        nday   = 0
        idxday = 0
        sfcalb = 0.0
      endif

      do i = 1, im
        sfcalb1(i) = sfcalb(i,1)
        sfcalb2(i) = sfcalb(i,2)
        sfcalb3(i) = sfcalb(i,3)
        sfcalb4(i) = sfcalb(i,4)
      enddo

      end subroutine rrtmg_sw_pre_run

!> \section arg_table_rrtmg_sw_pre_finalize Argument Table
!!
      subroutine rrtmg_sw_pre_finalize ()
      end subroutine rrtmg_sw_pre_finalize

!! @}
      end module rrtmg_sw_pre

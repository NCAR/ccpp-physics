!> \file fix_sys_bias_sfc_time.F90
!!  Modifies surface fluxes used in GFS-based PBL schemes.

!> This module contains the CCPP-compliant "fix_sys_bias_sfc_time" scheme.
    module fix_sys_bias_sfc_time
      
      use machine , only : kind_phys
      
      contains

      subroutine fix_sys_bias_sfc_time_init()
      end subroutine fix_sys_bias_sfc_time_init
      
      subroutine fix_sys_bias_sfc_time_finalize()
      end subroutine fix_sys_bias_sfc_time_finalize
      
!>  \brief  This subroutine contains all of the logic for the
!! fix_sys_bias_sfc_time scheme used in the CCPP-SCM online tutorial.
!!
!> \section arg_table_fix_sys_bias_sfc_time_run Argument Table
!! \htmlinclude fix_sys_bias_sfc_time_run.html
!!
      subroutine fix_sys_bias_sfc_time_run (im, con_cp, con_rd, con_hvap, p1, t1, solhr, hflx_r, qflx_r, errmsg, errflg)

      implicit none

!     arguments

      integer, intent(in) :: im
      
      real(kind=kind_phys), intent(in) :: con_cp, con_rd, con_hvap
      real(kind=kind_phys), intent(in) :: solhr
      real(kind=kind_phys), intent(in) :: p1(:), t1(:)
      real(kind=kind_phys), intent(inout) :: hflx_r(:), qflx_r(:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!     locals

      integer :: i
      real(kind=kind_phys) :: rho
      real(kind=kind_phys), parameter :: sens_mod_factor = 0 !W m-2
      real(kind=kind_phys), parameter :: lat_mod_factor = 40 !W m-2

!     Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      !Darwin, Australia has a local time offset of UTC + 9.5 hours;
      !If solhr (forecast UTC time) is outside of the Darwin local 6AM - 6PM, don't apply modified surface fluxes
      if (solhr > 8.5 .and. solhr < 20.5) return
      
      do i=1, im
        rho = p1(i)/(con_rd*t1(i))
        !convert mod_factors to kinematic units and add to heat fluxes
        hflx_r(i) = MAX(sens_mod_factor/(rho*con_cp) + hflx_r(i), 0.0)
        qflx_r(i) = MAX(lat_mod_factor/(rho*con_hvap) + qflx_r(i), 0.0)
      end do

    end subroutine fix_sys_bias_sfc_time_run

  end module fix_sys_bias_sfc_time

!>  \file sfc_noahmp_pre.F90
!!  This file contains data preparation for the NoahMP LSM for use in the GFS physics suite.

!> This module contains the CCPP-compliant data preparation for NoahMP LSM.
    module sfc_noahmp_pre

      implicit none

      private

      public :: sfc_noahmp_pre_init, sfc_noahmp_pre_run, sfc_noahmp_pre_finalize

      contains

      subroutine sfc_noahmp_pre_init()      
      end subroutine sfc_noahmp_pre_init

      subroutine sfc_noahmp_pre_finalize
      end subroutine sfc_noahmp_pre_finalize

!> \section arg_table_sfc_noahmp_pre_run Argument Table
!! \htmlinclude sfc_noahmp_pre_run.html
!!
!-----------------------------------
      subroutine sfc_noahmp_pre_run (im, lsm, lsm_noahmp, imp_physics,  &
        imp_physics_gfdl, imp_physics_mg, dtp, rain, rainc, ice, snow,  &
        graupel, rainn_mp, rainc_mp, ice_mp, snow_mp, graupel_mp,       &
        errmsg, errflg)     

        use machine ,   only : kind_phys

        implicit none
      
        integer,               intent(in) :: im, lsm, lsm_noahmp,       &
          imp_physics, imp_physics_gfdl, imp_physics_mg
        real (kind=kind_phys), intent(in) :: dtp
        real (kind=kind_phys), dimension(im), intent(in) :: rain, rainc,&
          ice, snow, graupel
        real (kind=kind_phys), dimension(:), intent(inout) :: rainn_mp, &
          rainc_mp, ice_mp, snow_mp, graupel_mp

        ! error messages
        character(len=*), intent(out)    :: errmsg
        integer,          intent(out)    :: errflg

        !  ---  locals:
        integer :: i
        real(kind=kind_phys) :: tem
        real(kind=kind_phys), parameter :: con_p001= 0.001d0

        !---  get the amount of different precip type for Noah MP
        !  ---  convert from m/dtp to mm/s
        if (lsm ==  lsm_noahmp .and. (imp_physics == imp_physics_mg .or. imp_physics == imp_physics_gfdl)) then
          tem = 1.0 / (dtp*con_p001)
          do  i=1,im
            rainn_mp(i)   = tem * (rain(i)-rainc(i))
            rainc_mp(i)   = tem * rainc(i)
            snow_mp(i)    = tem * snow(i)
            graupel_mp(i) = tem * graupel(i)
            ice_mp(i)     = tem * ice(i)
          enddo
        endif
        
      end subroutine sfc_noahmp_pre_run
    end module sfc_noahmp_pre

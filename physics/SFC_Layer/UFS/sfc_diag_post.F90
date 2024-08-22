!> \file sfc_diag_post.F90
!!  Contains code related to the surface diagnostic scheme.

      module sfc_diag_post

      contains

!>\defgroup sfc_diag_post_mod GFS sfc_diag_post Module
!! This module contains code related to the surface diagnostic scheme.
!> @{
#if 0
!> \section arg_table_sfc_diag_post_run Argument Table
!! \htmlinclude sfc_diag_post_run.html
!!
#endif
      subroutine sfc_diag_post_run (im, lsm, lsm_noahmp, opt_diag, dry, lssav, dtf, con_eps, con_epsm1, pgr,&
                 t2mmp,q2mp, t2m, q2m, u10m, v10m, tmpmin, tmpmax, spfhmin, spfhmax,                  &
                         wind10mmax, u10mmax, v10mmax, dpt2m, errmsg, errflg)

        use machine,               only: kind_phys, kind_dbl_prec

        implicit none

        integer,                             intent(in) :: im, lsm, lsm_noahmp,opt_diag
        logical,                             intent(in) :: lssav
        real(kind=kind_phys),                intent(in) :: dtf, con_eps, con_epsm1
        logical             , dimension(:),  intent(in) :: dry
        real(kind=kind_phys), dimension(:),  intent(in) :: pgr, u10m, v10m
        real(kind=kind_phys), dimension(:),  intent(inout) :: t2m, q2m, tmpmin, tmpmax, spfhmin, spfhmax
        real(kind=kind_phys), dimension(:),  intent(in), optional :: t2mmp, q2mp
        real(kind=kind_phys), dimension(:),  intent(inout) :: wind10mmax, u10mmax, v10mmax, dpt2m

        character(len=*),                     intent(out) :: errmsg
        integer,                              intent(out) :: errflg

        integer :: i
        real(kind=kind_dbl_prec) :: tem    ! made dbl prec always, JM 20211104

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (lsm == lsm_noahmp) then
         if (opt_diag == 2 .or. opt_diag == 3) then
          do i=1,im
            if(dry(i)) then
              t2m(i) = t2mmp(i)
              q2m(i) = q2mp(i)
             endif
          enddo
         endif
        endif

        if (lssav) then
          do i=1,im
            tmpmax(i)  = max(tmpmax(i),t2m(i))
            tmpmin(i)  = min(tmpmin(i),t2m(i))
            spfhmax(i) = max(spfhmax(i),q2m(i))
            spfhmin(i) = min(spfhmin(i),q2m(i))
          enddo
          ! Find max wind speed then decompose
          do i=1, im
             tem = sqrt(u10m(i)*u10m(i) + v10m(i)*v10m(i))
             if (tem > wind10mmax(i)) then
                wind10mmax(i) = tem
                u10mmax(i)    = u10m(i)
                v10mmax(i)    = v10m(i)
             endif
             ! Compute dew point, first using vapor pressure
             tem = max(pgr(i) * q2m(i) / ( con_eps - con_epsm1 *q2m(i)), 1.d-8)
             dpt2m(i) = 243.5_kind_dbl_prec / &
                ( ( 17.67_kind_dbl_prec / log(tem/611.2_kind_dbl_prec) ) - 1.) + 273.14
          enddo
        endif

      end subroutine sfc_diag_post_run
!> @}
      end module sfc_diag_post

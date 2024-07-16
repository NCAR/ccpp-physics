!>  \file sfc_nst_pre.f90
!!  This file contains preparation for the GFS NSST model.

module sfc_nst_pre

  use machine               , only : kind_phys
  use module_nst_water_prop , only : get_dtzm_2d
  use module_nst_parameters , only : zero, one

  implicit none

contains

  !> \defgroup GFS_NSST_PRE GFS Near-Surface Sea Temperature Pre
  !!
  !! The NSST scheme is one of the three schemes used to represent the
  !! surface in the GFS physics suite. The other two are the Noah land
  !! surface model and the sice simplified ice model.
  !!
  !! \section arg_table_sfc_nst_pre_run Argument Table
  !! \htmlinclude sfc_nst_pre_run.html
  !!
  !> \section NSST_general_pre_algorithm General Algorithm
  subroutine sfc_nst_pre_run                      &
       (im, wet, tgice, tsfco, tsurf_wat,         &
       tseal, xt, xz, dt_cool, z_c, tref, cplflx, &
       oceanfrac, nthreads, errmsg, errflg)

    !  ---  inputs:
    integer, intent(in) :: im, nthreads
    logical, dimension(:), intent(in) :: wet
    real (kind=kind_phys), intent(in) :: tgice
    real (kind=kind_phys), dimension(:), intent(in) :: tsfco, oceanfrac
    real (kind=kind_phys), dimension(:), intent(in), optional :: xt, xz, dt_cool, z_c
    logical, intent(in) :: cplflx

    !  ---  input/outputs:
    real (kind=kind_phys), dimension(:), intent(inout) :: tsurf_wat, tseal
    real (kind=kind_phys), dimension(:), intent(inout), optional :: tref

    !  ---  outputs:
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !  ---  locals
    integer :: i
    real(kind=kind_phys), parameter :: omz1 = 2.0_kind_phys
    real(kind=kind_phys) :: tem2, dnsst
    real(kind=kind_phys), dimension(im) :: dtzm, z_c_0

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    do i=1,im
       if (wet(i) .and. oceanfrac(i) > 0.0) then
          !          tem         = (oro(i)-oro_uf(i)) * rlapse
          ! DH* 20190927 simplyfing this code because tem is zero
          !tem          = zero
          !tseal(i)     = tsfco(i)  + tem
          tseal(i)      = tsfco(i)
          !tsurf_wat(i) = tsurf_wat(i) + tem
          ! *DH
       endif
    enddo
    !
    !   update tsfc & tref with T1 from OGCM & NSST Profile if coupled
    !
    if (cplflx) then
       z_c_0 = zero
       call get_dtzm_2d (xt, xz, dt_cool, z_c_0, wet, zero, omz1, im, 1, nthreads, dtzm)
       do i=1,im
          if (wet(i) .and. oceanfrac(i) > zero ) then
             !           dnsst   = tsfc_wat(i) - tref(i)           !  retrive/get difference of Ts and Tf
             tref(i) = max(tgice, tsfco(i) - dtzm(i))              !  update Tf with T1 and NSST T-Profile
             !           tsfc_wat(i) = max(271.2,tref(i) + dnsst)  !  get Ts updated due to Tf update
             !           tseal(i)    = tsfc_wat(i)
             if (abs(xz(i)) > zero) then
                tem2 = one / xz(i)
             else
                tem2 = zero
             endif
             tseal(i)     = tref(i) + (xt(i)+xt(i)) * tem2 - dt_cool(i)
             tsurf_wat(i) = tseal(i)
          endif
       enddo
    endif

    return
  end subroutine sfc_nst_pre_run
end module sfc_nst_pre

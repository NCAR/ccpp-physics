!>  \file sfc_nst_post.f
!!  This file contains code to be executed after the GFS NSST model.

      module sfc_nst_post

      contains

! \defgroup GFS_NSST_POST GFS Near-Surface Sea Temperature Post

!> \section arg_table_sfc_nst_post_run Argument Table
!! \htmlinclude sfc_nst_post_run.html
!!
! \section NSST_general_post_algorithm General Algorithm
!
! \section NSST_detailed_post_algorithm Detailed Algorithm
! @{
      subroutine sfc_nst_post_run                                       &
     &     ( im, kdt, rlapse, tgice, wet, use_lake_model, icy, oro,     &
     &       oro_uf, nstf_name1,                                        &
     &       nstf_name4, nstf_name5, xt, xz, dt_cool, z_c, tref, xlon,  &
     &       tsurf_wat, tsfc_wat, nthreads, dtzm, errmsg, errflg        &
     &     )

      use machine , only : kind_phys
      use module_nst_water_prop, only: get_dtzm_2d

      implicit none

      integer, parameter :: kp = kind_phys

!  ---  inputs:
      integer, intent(in) :: im, kdt, nthreads
      logical, dimension(:), intent(in) :: wet, icy
      integer, dimension(:), intent(in) :: use_lake_model
      real (kind=kind_phys), intent(in) :: rlapse, tgice
      real (kind=kind_phys), dimension(:), intent(in) :: oro, oro_uf
      integer, intent(in) :: nstf_name1, nstf_name4, nstf_name5
      real (kind=kind_phys), dimension(:), intent(in) :: xt, xz,        &
     &      dt_cool, z_c, tref, xlon

!  ---  input/outputs:
      real (kind=kind_phys), dimension(:), intent(inout) :: tsurf_wat,  &
     &      tsfc_wat

!  ---  outputs:
      real (kind=kind_phys), dimension(:), intent(out) :: dtzm

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals
      integer :: i
      real(kind=kind_phys) :: zsea1, zsea2

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!     if (lprnt) print *,' tseaz2=',tseal(ipr),' tref=',tref(ipr),
!    &     ' dt_cool=',dt_cool(ipr),' dt_warm=',2.0*xt(ipr)/xz(ipr),
!    &     ' kdt=',kdt

!      do i = 1, im
!        if (wet(i) .and. .not. icy(i)) then
!          tsurf_wat(i) = tsurf_wat(i) - (oro(i)-oro_uf(i)) * rlapse
!        endif
!      enddo

!  --- ...  run nsst model  ... ---

      if (nstf_name1 > 1) then
        zsea1 = 0.001_kp*real(nstf_name4)
        zsea2 = 0.001_kp*real(nstf_name5)
        call get_dtzm_2d (xt, xz, dt_cool, z_c, wet, zsea1, zsea2,      &
     &                    im, 1, nthreads, dtzm)
        do i = 1, im
!         if (wet(i) .and. .not.icy(i)) then
!         if (wet(i) .and. (frac_grid .or. .not. icy(i))) then
          if (wet(i) .and. use_lake_model(i) /=1) then
            tsfc_wat(i) = max(tgice, tref(i) + dtzm(i))
!           tsfc_wat(i) = max(271.2, tref(i) + dtzm(i)) -  &
!                           (oro(i)-oro_uf(i))*rlapse
          endif
        enddo
      endif

!     if (lprnt) print *,' tseaz2=',tsea(ipr),' tref=',tref(ipr),   &
!    &    ' dt_cool=',dt_cool(ipr),' dt_warm=',dt_warm(ipr),' kdt=',kdt

      return
      end subroutine sfc_nst_post_run

      end module sfc_nst_post

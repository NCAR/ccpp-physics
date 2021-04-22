!> \file GFS_DCNV_generic.F90
!!  Contains code related to deep convective schemes to be used within the GFS physics suite.

      module GFS_DCNV_generic_pre

      contains

      subroutine GFS_DCNV_generic_pre_init ()
      end subroutine GFS_DCNV_generic_pre_init

      subroutine GFS_DCNV_generic_pre_finalize()
      end subroutine GFS_DCNV_generic_pre_finalize

!> \brief Interstitial scheme called prior to any deep convective scheme to save state variables for calculating tendencies after the deep convective scheme is executed
!! \section arg_table_GFS_DCNV_generic_pre_run Argument Table
!! \htmlinclude GFS_DCNV_generic_pre_run.html
!!
    subroutine GFS_DCNV_generic_pre_run (im, levs, ldiag3d, qdiag3d, do_cnvgwd, cplchm,&
                                         gu0, gv0, gt0, gq0_water_vapor,               &
                                         save_u, save_v, save_t, save_qv, dqdti,       &
                                         errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      integer, intent(in) :: im, levs
      logical, intent(in) :: ldiag3d, qdiag3d, do_cnvgwd, cplchm
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: gu0
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: gv0
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: gt0
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: gq0_water_vapor
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_u
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_v
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_t
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_qv
      ! dqdti only allocated if cplchm is .true.
      real(kind=kind_phys), dimension(:,:),     intent(inout) :: dqdti
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      real(kind=kind_phys), parameter :: zero    = 0.0d0
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (ldiag3d) then
        do k=1,levs
          do i=1,im
            save_t(i,k) = gt0(i,k)
            save_u(i,k) = gu0(i,k)
            save_v(i,k) = gv0(i,k)
          enddo
        enddo
      elseif (do_cnvgwd) then
        do k=1,levs
          do i=1,im
            save_t(i,k) = gt0(i,k)
          enddo
        enddo
      endif

      if ((ldiag3d.and.qdiag3d) .or. cplchm) then
        do k=1,levs
          do i=1,im
            save_qv(i,k) = gq0_water_vapor(i,k)
          enddo
        enddo
      endif

      if (cplchm) then
        dqdti = zero
      endif

    end subroutine GFS_DCNV_generic_pre_run

    end module GFS_DCNV_generic_pre

    module GFS_DCNV_generic_post

    contains

    subroutine GFS_DCNV_generic_post_init ()
    end subroutine GFS_DCNV_generic_post_init

    subroutine GFS_DCNV_generic_post_finalize ()
    end subroutine GFS_DCNV_generic_post_finalize

!> \section arg_table_GFS_DCNV_generic_post_run Argument Table
!! \htmlinclude GFS_DCNV_generic_post_run.html
!!
    subroutine GFS_DCNV_generic_post_run (im, levs, lssav, ldiag3d, qdiag3d, ras, cscnv, &
      frain, rain1, dtf, cld1d, save_u, save_v, save_t, save_qv, gu0, gv0, gt0,          &
      gq0_water_vapor, ud_mf, dd_mf, dt_mf, con_g, npdf3d, num_p3d, ncnvcld3d,           &
      rainc, cldwrk, dt3dt, dq3dt, du3dt, dv3dt, upd_mf, dwn_mf, det_mf,                 &
      cnvw, cnvc, cnvw_phy_f3d, cnvc_phy_f3d, flag_for_dcnv_generic_tend, errmsg, errflg)


      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs
      logical, intent(in) :: lssav, ldiag3d, qdiag3d, ras, cscnv
      logical, intent(in) :: flag_for_dcnv_generic_tend

      real(kind=kind_phys), intent(in) :: frain, dtf
      real(kind=kind_phys), dimension(im), intent(in) :: rain1, cld1d
      real(kind=kind_phys), dimension(im,levs), intent(in) :: save_u, save_v, save_t, save_qv
      real(kind=kind_phys), dimension(im,levs), intent(in) :: gu0, gv0, gt0, gq0_water_vapor
      real(kind=kind_phys), dimension(im,levs), intent(in) :: ud_mf, dd_mf, dt_mf
      real(kind=kind_phys), intent(in) :: con_g
      integer, intent(in) :: npdf3d, num_p3d, ncnvcld3d

      real(kind=kind_phys), dimension(im), intent(inout) :: rainc, cldwrk
      ! dt3dt, dq3dt, du3dt, dv3dt upd_mf, dwn_mf, det_mf only allocated if ldiag3d == .true.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dt3dt, dq3dt, du3dt, dv3dt
      real(kind=kind_phys), dimension(:,:), intent(inout) :: upd_mf, dwn_mf, det_mf
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: cnvw, cnvc
      ! The following arrays may not be allocated, depending on certain flags and microphysics schemes.
      ! Since Intel 15 crashes when passing unallocated arrays to arrays defined with explicit shape,
      ! use assumed-shape arrays. Note that Intel 18 and GNU 6.2.0-8.1.0 tolerate explicit-shape arrays
      ! as long as these do not get used when not allocated (it is still invalid Fortran code, though).
      real(kind=kind_phys), dimension(:,:), intent(inout) :: cnvw_phy_f3d, cnvc_phy_f3d

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (.not. ras .and. .not. cscnv) then
        if (npdf3d == 3 .and. num_p3d == 4) then
          do k=1,levs
            do i=1,im
              cnvw_phy_f3d(i,k) = cnvw(i,k)
              cnvc_phy_f3d(i,k) = cnvc(i,k)
              cnvw(i,k)         = 0.0
              cnvc(i,k)         = 0.0
            enddo
          enddo
        elseif (npdf3d == 0 .and. ncnvcld3d == 1) then
          do k=1,levs
            do i=1,im
              cnvw_phy_f3d(i,k) = cnvw(i,k)
              cnvw(i,k)         = 0.0
            enddo
          enddo
        endif
      endif ! if (.not. ras .and. .not. cscnv)

      do i=1,im
        rainc(i) = frain * rain1(i)
      enddo
!
      if (lssav) then
        do i=1,im
          cldwrk (i)  = cldwrk (i)  + cld1d(i) * dtf
        enddo

        if (ldiag3d .and. flag_for_dcnv_generic_tend) then
          do k=1,levs
            do i=1,im
              dt3dt(i,k) = dt3dt(i,k) + (gt0(i,k)-save_t(i,k)) * frain
              du3dt(i,k) = du3dt(i,k) + (gu0(i,k)-save_u(i,k)) * frain
              dv3dt(i,k) = dv3dt(i,k) + (gv0(i,k)-save_v(i,k)) * frain
              ! convective mass fluxes
              upd_mf(i,k)  = upd_mf(i,k)  + ud_mf(i,k) * (con_g*frain)
              dwn_mf(i,k)  = dwn_mf(i,k)  + dd_mf(i,k) * (con_g*frain)
              det_mf(i,k)  = det_mf(i,k)  + dt_mf(i,k) * (con_g*frain)
            enddo
          enddo
          if(qdiag3d) then
             do k=1,levs
                do i=1,im
                   dq3dt(i,k) = dq3dt(i,k) + (gq0_water_vapor(i,k)-save_qv(i,k)) * frain
                enddo
             enddo
          endif
        endif ! if (ldiag3d)

      endif ! if (lssav)

    end subroutine GFS_DCNV_generic_post_run

    end module GFS_DCNV_generic_post

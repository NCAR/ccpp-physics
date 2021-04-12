!> \file GFS_SCNV_generic.F90
!!  Contains code related to shallow convective schemes to be used within the GFS physics suite.

      module GFS_SCNV_generic_pre

      contains

      subroutine GFS_SCNV_generic_pre_init ()
      end subroutine GFS_SCNV_generic_pre_init

      subroutine GFS_SCNV_generic_pre_finalize()
      end subroutine GFS_SCNV_generic_pre_finalize

!> \section arg_table_GFS_SCNV_generic_pre_run Argument Table
!! \htmlinclude GFS_SCNV_generic_pre_run.html
!!
      subroutine GFS_SCNV_generic_pre_run (im, levs, ldiag3d, qdiag3d, gu0, gv0, gt0, gq0, &
        save_u, save_v, save_t, save_q, ntqv, nsamftrac, flag_for_scnv_generic_tend,       &
        dtidx, index_of_process_scnv, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none

        integer, intent(in) :: im, levs, ntqv, nsamftrac, index_of_process_scnv, dtidx(:,:)
        logical, intent(in) :: ldiag3d, qdiag3d, flag_for_scnv_generic_tend
        real(kind=kind_phys), dimension(im,levs), intent(in) :: gu0, gv0, gt0
        real(kind=kind_phys), intent(in) :: gq0(:,:,:)
        real(kind=kind_phys), intent(inout) :: save_q(:,:,:)
        real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_u, save_v, save_t
        character(len=*),                 intent(out) :: errmsg
        integer,                          intent(out) :: errflg

        integer :: i, k, n

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (ldiag3d .and. flag_for_scnv_generic_tend) then
          do k=1,levs
            do i=1,im
              save_u(i,k)   = gu0(i,k)
              save_v(i,k)   = gv0(i,k)
              save_t(i,k)   = gt0(i,k)
            enddo
          enddo
          if (qdiag3d) then
            if(nsamftrac>0) then
              do n=1,nsamftrac
                if(n==ntqv .or. dtidx(ntqv,index_of_process_scnv)>=1) then
                  save_q(:,:,n) = gq0(:,:,n)
                endif
              enddo
            else
              save_q(:,:,ntqv) = gq0(:,:,ntqv)
            endif
          endif
       endif

    end subroutine GFS_SCNV_generic_pre_run


    end module GFS_SCNV_generic_pre

    module GFS_SCNV_generic_post

    contains

    subroutine GFS_SCNV_generic_post_init ()
    end subroutine GFS_SCNV_generic_post_init

    subroutine GFS_SCNV_generic_post_finalize ()
    end subroutine GFS_SCNV_generic_post_finalize

!> \section arg_table_GFS_SCNV_generic_post_run Argument Table
!! \htmlinclude GFS_SCNV_generic_post_run.html
!!
      subroutine GFS_SCNV_generic_post_run (im, levs, nn, lssav, ldiag3d, qdiag3d, &
        cplchm, frain, gu0, gv0, gt0, gq0, save_u, save_v, save_t, save_q, dqdti,  &
        clw, shcnvcw, rain1, npdf3d, num_p3d, ncnvcld3d, cnvc, cnvw, nsamftrac,    &
        rainc, cnvprcp, cnvprcpb, cnvw_phy_f3d, cnvc_phy_f3d,                      &
        dtend, dtidx, index_of_temperature, index_of_x_wind, index_of_y_wind,      &
        index_of_process_scnv, ntqv, flag_for_scnv_generic_tend,                   &
        imfshalcnv, imfshalcnv_sas, imfshalcnv_samf, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, nn, ntqv, nsamftrac
      logical, intent(in) :: lssav, ldiag3d, qdiag3d, cplchm, flag_for_scnv_generic_tend
      real(kind=kind_phys),                     intent(in) :: frain
      real(kind=kind_phys), dimension(im,levs), intent(in) :: gu0, gv0, gt0
      real(kind=kind_phys), dimension(im,levs), intent(in) :: save_u, save_v, save_t
      real(kind=kind_phys), dimension(:,:,:),   intent(in) :: save_q, gq0

      ! dtend only allocated if ldiag3d == .true.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dqdti
      real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:)
      integer, intent(in) :: index_of_temperature, index_of_x_wind, index_of_y_wind, index_of_process_scnv
      real(kind=kind_phys), dimension(im,levs,nn), intent(inout) :: clw

      ! Post code for SAS/SAMF
      integer, intent(in) :: npdf3d, num_p3d, ncnvcld3d
      logical, intent(in) :: shcnvcw
      real(kind=kind_phys), dimension(im), intent(in) :: rain1
      real(kind=kind_phys), dimension(im,levs), intent(in) :: cnvw, cnvc
      real(kind=kind_phys), dimension(im), intent(inout) :: rainc, cnvprcp, cnvprcpb
      ! The following arrays may not be allocated, depending on certain flags and microphysics schemes.
      ! Since Intel 15 crashes when passing unallocated arrays to arrays defined with explicit shape,
      ! use assumed-shape arrays. Note that Intel 18 and GNU 6.2.0-8.1.0 tolerate explicit-shape arrays
      ! as long as these do not get used when not allocated.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: cnvw_phy_f3d, cnvc_phy_f3d
      integer, intent(in) :: imfshalcnv, imfshalcnv_sas, imfshalcnv_samf

      character(len=*),              intent(out) :: errmsg
      integer,                       intent(out) :: errflg

      integer :: i, k, n, idtend
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (imfshalcnv==imfshalcnv_sas .or. imfshalcnv==imfshalcnv_samf) then
        do i=1,im
          rainc(i) = rainc(i) + frain * rain1(i)
        enddo
! 'cnvw' and 'cnvc' are set to zero before computation starts:
        if (shcnvcw .and. num_p3d == 4 .and. npdf3d == 3) then
          do k=1,levs
            do i=1,im
              cnvw_phy_f3d(i,k) = cnvw_phy_f3d(i,k) + cnvw(i,k)
              cnvc_phy_f3d(i,k) = cnvc_phy_f3d(i,k) + cnvc(i,k)
            enddo
          enddo
        elseif (npdf3d == 0 .and. ncnvcld3d == 1) then
          do k=1,levs
            do i=1,im
              cnvw_phy_f3d(i,k) = cnvw_phy_f3d(i,k) +  cnvw(i,k)
            enddo
          enddo
        endif
      endif

      if (lssav .and. flag_for_scnv_generic_tend) then
        if (ldiag3d) then
          idtend = dtidx(index_of_temperature, index_of_process_scnv)
          if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + (gt0 - save_t) * frain
          endif

          idtend = dtidx(index_of_x_wind, index_of_process_scnv)
          if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + (gu0 - save_u) * frain
          endif

          idtend = dtidx(index_of_y_wind, index_of_process_scnv)
          if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + (gv0 - save_v) * frain
          endif

          if(nsamftrac>0) then
            do n=1,nsamftrac
              idtend = dtidx(100+n, index_of_process_scnv)
              if(idtend>=1) then
                dtend(:,:,idtend) = dtend(:,:,idtend) + (gq0(:,:,n) - save_q(:,:,n)) * frain
              endif
            enddo
          else
            idtend = dtidx(100+ntqv, index_of_process_scnv)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + (gq0(:,:,ntqv) - save_q(:,:,ntqv)) * frain
            endif
          endif
        endif
      endif
!
      if (cplchm) then
        do k=1,levs
          do i=1,im
            tem  = (gq0(i,k,ntqv)-save_q(i,k,ntqv)) * frain
            dqdti(i,k) = dqdti(i,k) + tem
          enddo
        enddo
      endif
!
      do k=1,levs
        do i=1,im
          if (clw(i,k,2) <= -999.0) clw(i,k,2) = 0.0
        enddo
      enddo

      end subroutine GFS_SCNV_generic_post_run

      end module GFS_SCNV_generic_post

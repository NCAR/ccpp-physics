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
    subroutine GFS_DCNV_generic_pre_run (im, levs, ldiag3d, qdiag3d, do_cnvgwd, cplchm,  &
                                         gu0, gv0, gt0, gq0, nsamftrac, ntqv,            &
                                         save_u, save_v, save_t, save_q, dqdti, clw,     &
                                         ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,   &
                                         ntgnc, cscnv, satmedmf, trans_trac, ras, ntrac, &
                                         dtidx, index_of_process_dcnv, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, nsamftrac, ntqv, index_of_process_dcnv, dtidx(:,:), &
           ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,ntrac,ntgnc
      logical, intent(in) :: ldiag3d, qdiag3d, do_cnvgwd, cplchm
      real(kind=kind_phys), dimension(:,:),   intent(in)    :: gu0
      real(kind=kind_phys), dimension(:,:),   intent(in)    :: gv0
      real(kind=kind_phys), dimension(:,:),   intent(in)    :: gt0
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: gq0
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: save_u
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: save_v
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: save_t
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: save_q
      ! dqdti only allocated if cplchm is .true.
      real(kind=kind_phys), dimension(:,:),     intent(inout) :: dqdti
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg
      logical, intent(in) :: cscnv, satmedmf, trans_trac, ras
      real(kind=kind_phys), parameter :: zero    = 0.0d0
      real(kind=kind_phys), dimension(:,:,:), intent(in) :: clw

      integer :: i, k, n, tracers

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
         if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
            tracers = 2
            do n=2,ntrac
               if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
                    n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
                    n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc) then
                  tracers = tracers + 1
                  if(dtidx(100+n,index_of_process_dcnv)>0) then
                     save_q(:,:,n) = clw(:,:,tracers)
                  endif
               endif
            enddo
         else
            do n=2,ntrac
               if(dtidx(100+n,index_of_process_dcnv)>0) then
                  save_q(:,:,n) = gq0(:,:,n)
               endif
            enddo
         endif ! end if_ras or cfscnv or samf
         save_q(:,:,ntqv) = gq0(:,:,ntqv)
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
    subroutine GFS_DCNV_generic_post_run (im, levs, lssav, ldiag3d, qdiag3d, ras, &
      cscnv, frain, rain1, dtf, cld1d, save_u, save_v, save_t, gu0, gv0, gt0,     &
      ud_mf, dd_mf, dt_mf, con_g, npdf3d, num_p3d, ncnvcld3d, nsamftrac,          &
      rainc, cldwrk, upd_mf, dwn_mf, det_mf, dtend, dtidx, index_of_process_dcnv, &
      index_of_temperature, index_of_x_wind, index_of_y_wind, ntqv, gq0, save_q,  &
      cnvw, cnvc, cnvw_phy_f3d, cnvc_phy_f3d, flag_for_dcnv_generic_tend,         &
      ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,ntgnc, ntrac,clw,              &
      satmedmf, trans_trac, errmsg, errflg)


      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, nsamftrac
      logical, intent(in) :: lssav, ldiag3d, qdiag3d, ras, cscnv
      logical, intent(in) :: flag_for_dcnv_generic_tend

      real(kind=kind_phys), intent(in) :: frain, dtf
      real(kind=kind_phys), dimension(:),     intent(in) :: rain1, cld1d
      real(kind=kind_phys), dimension(:,:),   intent(in) :: save_u, save_v, save_t
      real(kind=kind_phys), dimension(:,:),   intent(in) :: gu0, gv0, gt0
      real(kind=kind_phys), dimension(:,:,:), intent(in) :: gq0, save_q
      real(kind=kind_phys), dimension(:,:),   intent(in) :: ud_mf, dd_mf, dt_mf
      real(kind=kind_phys), intent(in) :: con_g
      integer, intent(in) :: npdf3d, num_p3d, ncnvcld3d
      logical, intent(in) :: satmedmf, trans_trac

      real(kind=kind_phys), dimension(:),   intent(inout) :: rainc, cldwrk
      real(kind=kind_phys), dimension(:,:), intent(inout) :: upd_mf, dwn_mf, det_mf
      real(kind=kind_phys), dimension(:,:), intent(inout) :: cnvw, cnvc

      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: dtend
      integer, intent(in) :: dtidx(:,:), index_of_process_dcnv, index_of_temperature, &
           index_of_x_wind, index_of_y_wind, ntqv
      integer, intent(in) :: ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,ntrac,ntgnc
      real(kind=kind_phys), dimension(:,:,:), intent(in) :: clw


      real(kind=kind_phys), dimension(:,:), intent(inout) :: cnvw_phy_f3d, cnvc_phy_f3d

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k, n, idtend, tracers

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
          idtend=dtidx(index_of_temperature,index_of_process_dcnv)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + (gt0-save_t)*frain
          endif

          idtend=dtidx(index_of_x_wind,index_of_process_dcnv)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + (gu0-save_u)*frain
          endif

          idtend=dtidx(index_of_y_wind,index_of_process_dcnv)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + (gv0-save_v)*frain
          endif

          if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
             tracers = 2
             do n=2,ntrac
                if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
                     n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
                     n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc) then
                   tracers = tracers + 1
                   idtend = dtidx(100+n,index_of_process_dcnv)
                   if(idtend>0) then
                      dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,tracers)-save_q(:,:,n) * frain
                   endif
                endif
             enddo
          else
            do n=2,ntrac
               idtend = dtidx(100+n,index_of_process_dcnv)
               if(idtend>0) then
                  dtend(:,:,idtend) = dtend(:,:,idtend) + (gq0(:,:,n)-save_q(:,:,n))*frain
               endif
            enddo
          endif
          idtend = dtidx(100+ntqv, index_of_process_dcnv)
          if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + (gq0(:,:,ntqv) - save_q(:,:,ntqv)) * frain
          endif

          ! convective mass fluxes
          if(qdiag3d) then
            do k=1,levs
              do i=1,im
                upd_mf(i,k)  = upd_mf(i,k)  + ud_mf(i,k) * (con_g*frain)
                dwn_mf(i,k)  = dwn_mf(i,k)  + dd_mf(i,k) * (con_g*frain)
                det_mf(i,k)  = det_mf(i,k)  + dt_mf(i,k) * (con_g*frain)
              enddo
            enddo
          endif
        endif ! if (ldiag3d)

      endif ! if (lssav)

    end subroutine GFS_DCNV_generic_post_run

    end module GFS_DCNV_generic_post

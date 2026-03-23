!> \file GFS_SCNV_generic_post.F90
!!  Contains code related to shallow convective schemes to be used after shallow convection for GFS-based physics suites.

    module GFS_SCNV_generic_post

    contains

!> \section arg_table_GFS_SCNV_generic_post_run Argument Table
!! \htmlinclude GFS_SCNV_generic_post_run.html
!!
      subroutine GFS_SCNV_generic_post_run (im, levs, tracers_total, otsptflag, imp_physics, imp_physics_gfdl, tend_opt_scnv, lssav, ldiag3d, qdiag3d, &
        frain, gu0, gv0, gt0, gq0, dudt, dvdt, dtdt, dqdt, ten_t, ten_u, ten_v, ten_q, delt,              &
        clw, dclw, shcnvcw, rain1, npdf3d, num_p3d, ncnvcld3d, cnvc, cnvw, nsamftrac,    &
        rainc, cnvprcp, cnvw_phy_f3d, cnvc_phy_f3d,                      &
        dtend, dtidx, index_of_temperature, index_of_x_wind, index_of_y_wind,      &
        index_of_process_scnv, ntqv, flag_for_scnv_generic_tend,                   &
        ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,ntgnc,ntsigma,                &
        imfshalcnv, imfshalcnv_sas, imfshalcnv_samf, ntrac,                        &
        cscnv, satmedmf, trans_trac, ras, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, ntqv, nsamftrac, tracers_total, tend_opt_scnv
      integer, intent(in) :: imp_physics, imp_physics_gfdl
      integer, intent(in) :: ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,ntgnc,ntsigma,ntrac
      logical, intent(in) :: lssav, ldiag3d, qdiag3d, flag_for_scnv_generic_tend
      logical, dimension(:), intent(in) :: otsptflag
      real(kind=kind_phys),                     intent(in) :: frain
      real(kind=kind_phys), dimension(:,:), intent(inout) :: gu0, gv0, gt0
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: gq0

      ! dtend only allocated if ldiag3d == .true.
      real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:)
      integer, intent(in) :: index_of_temperature, index_of_x_wind, index_of_y_wind, index_of_process_scnv
      real(kind=kind_phys), dimension(:,:,:), intent(in) :: clw, dclw

      ! Post code for SAS/SAMF
      integer, intent(in) :: npdf3d, num_p3d, ncnvcld3d
      logical, intent(in) :: shcnvcw
      real(kind=kind_phys), dimension(:), intent(in) :: rain1
      real(kind=kind_phys), dimension(:, :), intent(in) :: cnvw, cnvc
      real(kind=kind_phys), dimension(:), intent(inout) :: rainc, cnvprcp
      ! The following arrays may not be allocated, depending on certain flags and microphysics schemes.
      ! Since Intel 15 crashes when passing unallocated arrays to arrays defined with explicit shape,
      ! use assumed-shape arrays. Note that Intel 18 and GNU 6.2.0-8.1.0 tolerate explicit-shape arrays
      ! as long as these do not get used when not allocated.
      real(kind=kind_phys), dimension(:,:), intent(inout), optional :: cnvw_phy_f3d, cnvc_phy_f3d
      integer, intent(in) :: imfshalcnv, imfshalcnv_sas, imfshalcnv_samf
      logical, intent(in) :: cscnv, satmedmf, trans_trac, ras

      character(len=*),              intent(out) :: errmsg
      integer,                       intent(out) :: errflg

      integer :: i, k, n, idtend, tracers
      real(kind=kind_phys) :: tem

      real(kind=kind_phys), dimension(:,:), intent(in) :: ten_t, ten_u, ten_v
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: ten_q
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dudt, dvdt, dtdt 
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: dqdt
      real(kind=kind_phys), intent(in) ::  delt

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      !ten_q(:,:,1) already has a value from the shallow convection scheme
      if (tracers_total > 0) then
        tracers = 2
        do n=2,ntrac
          if ( otsptflag(n) ) then                                                   
            tracers = tracers + 1
            ten_q(1:im,:,n) = dclw(1:im,:,tracers)
          endif
        enddo
      endif
      if (ntcw > 0) then
        if (imp_physics == imp_physics_gfdl) then
           ten_q(1:im,:,ntcw) = dclw(1:im,:,1) + dclw(1:im,:,2)
        elseif (ntiw > 0) then
          ten_q(1:im,:,ntiw) = dclw(1:im,:,1)
          ten_q(1:im,:,ntcw) = dclw(1:im,:,2)
        else
          ten_q(1:im,:,ntcw) = dclw(1:im,:,1) + dclw(1:im,:,2)
        endif   ! end if_ntiw
      endif   ! end if_ntcw
      
      case_SCNV_ten: select case (tend_opt_scnv)
        case (1) !immediately apply tendencies
                  !Current state = current state + dt*current tendency
                  !Accumulated tendency unchanged
          do k=1,levs
            do i=1,im
              gt0(i,k) = gt0(i,k) + delt*ten_t(i,k)
              gu0(i,k) = gu0(i,k) + delt*ten_u(i,k)
              gv0(i,k) = gv0(i,k) + delt*ten_v(i,k)
              do n = 1, ntrac
                gq0(i,k,n) = gq0(i,k,n) + delt*ten_q(i,k,n)
              end do
            end do
          end do
        case (2) !add tendencies to sum
                  !Accumulated tendency = accumulated tendency + current tendency
                  !Current state unchanged
          do k=1,levs
            do i=1,im
              dtdt(i,k) = dtdt(i,k) + ten_t(i,k)
              dudt(i,k) = dudt(i,k) + ten_u(i,k)
              dvdt(i,k) = dvdt(i,k) + ten_v(i,k)
              do n = 1, ntrac
                dqdt(i,k,n) = dqdt(i,k,n) + ten_q(i,k,n)
              end do
            end do
          end do
        case (3) !add tendencies to sum and apply
                  !Current state = current state + dt*(accumulated tendency + current tendency)
                  !Accumulated tendency = 0
          do k=1,levs
            do i=1,im
              gt0(i,k) = gt0(i,k) + delt*(dtdt(i,k) + ten_t(i,k))
              dtdt(i,k) = 0.0
              gu0(i,k) = gu0(i,k) + delt*(dudt(i,k) + ten_u(i,k))
              dudt(i,k) = 0.0
              gv0(i,k) = gv0(i,k) + delt*(dvdt(i,k) + ten_v(i,k))
              dvdt(i,k) = 0.0
              do n = 1, ntrac
                gq0(i,k,n) = gq0(i,k,n) + delt*(dqdt(i,k,n) + ten_q(i,k,n))
                dqdt(i,k,n) = 0.0
              end do
            end do
          end do
        case (4) !Current state unchanged
                  !Accumulated tendency unchanged
                  !Current tendency unchanged (but will be overwritten during next primary scheme)
          exit case_SCNV_ten
        case default
          errflg = 1
          errmsg = 'A tendency application control was outside of the acceptable range (1-4)'
          return
      end select case_SCNV_ten

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
             dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_t * delt) * frain
          endif

          idtend = dtidx(index_of_x_wind, index_of_process_scnv)
          if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_u * delt) * frain
          endif

          idtend = dtidx(index_of_y_wind, index_of_process_scnv)
          if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_v * delt) * frain
          endif

          if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
             tracers = 2
             do n=2,ntrac
                if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
                     n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
                     n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc   .and. n /= ntsigma) then
                   tracers = tracers + 1
                   idtend = dtidx(100+n,index_of_process_scnv)
                   if(idtend>0) then
                      dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_q(:,:,n)*delt) * frain
                   endif
                endif
             enddo
          else
            do n=2,ntrac
               idtend = dtidx(100+n,index_of_process_scnv)
               if(idtend>0) then
                  dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_q(:,:,n)*delt)*frain
               endif
            enddo
          endif
          idtend = dtidx(100+ntqv, index_of_process_scnv)
          if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_q(:,:,ntqv)*delt) * frain
          endif
        endif
      endif

      end subroutine GFS_SCNV_generic_post_run

      end module GFS_SCNV_generic_post

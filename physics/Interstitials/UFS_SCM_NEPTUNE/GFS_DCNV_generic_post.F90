!> \file GFS_DCNV_generic_post.F90
!!  Contains code related to deep convective schemes to be used within the GFS physics suite.

    module GFS_DCNV_generic_post

    contains

!> \section arg_table_GFS_DCNV_generic_post_run Argument Table
!! \htmlinclude GFS_DCNV_generic_post_run.html
!!
    subroutine GFS_DCNV_generic_post_run (im, levs, tracers_total, otsptflag,     &
      imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_nssl,      &
      imp_physics_wsm6, imp_physics_mg, imp_physics_fer_hires, tend_opt_dcnv,     &
      lssav, ldiag3d, qdiag3d, ras, cscnv, frain, rain1, dtf, cld1d, gu0, gv0,    &
      gt0, ten_t, ten_u, ten_v, ten_q, dudt, dvdt, dtdt, dqdt,                    &
      delt, ud_mf, dd_mf, dt_mf, con_g, npdf3d, num_p3d, ncnvcld3d, nsamftrac,    &
      rainc, cldwrk, upd_mf, dwn_mf, det_mf, dtend, dtidx, index_of_process_dcnv, &
      index_of_temperature, index_of_x_wind, index_of_y_wind, ntqv, gq0,          &
      cnvw, cnvc, cnvw_phy_f3d, cnvc_phy_f3d, flag_for_dcnv_generic_tend,         &
      ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,                               &
      ntgnc, nthl, nthnc, nthv, ntgv, ntrz, ntgz, nthz, ntsigma, ntomega, ntrac,  &
      clw, dclw, satmedmf, trans_trac, errmsg, errflg)


      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, nsamftrac, tracers_total, tend_opt_dcnv
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_nssl, imp_physics_wsm6, imp_physics_mg, imp_physics_fer_hires
      logical, intent(in) :: lssav, ldiag3d, qdiag3d, ras, cscnv
      logical, intent(in) :: flag_for_dcnv_generic_tend
      logical, dimension(:), intent(in) :: otsptflag

      real(kind=kind_phys), intent(in) :: frain, dtf
      real(kind=kind_phys), dimension(:),     intent(in) :: rain1, cld1d
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: gu0, gv0, gt0
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: gq0
      real(kind=kind_phys), dimension(:,:),   intent(in) :: dd_mf, dt_mf
      real(kind=kind_phys), dimension(:,:),   intent(in), optional :: ud_mf
      real(kind=kind_phys), intent(in) :: con_g
      integer, intent(in) :: npdf3d, num_p3d, ncnvcld3d
      logical, intent(in) :: satmedmf, trans_trac

      real(kind=kind_phys), dimension(:),   intent(inout) :: rainc, cldwrk
      real(kind=kind_phys), dimension(:,:), intent(inout), optional :: upd_mf, dwn_mf, det_mf
      real(kind=kind_phys), dimension(:,:), intent(inout) :: cnvw, cnvc

      real(kind=kind_phys), dimension(:,:,:), intent(inout), optional :: dtend
      integer, intent(in) :: dtidx(:,:), index_of_process_dcnv, index_of_temperature, &
           index_of_x_wind, index_of_y_wind, ntqv
      integer, intent(in) :: ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,     &
                             ntgnc, nthl, nthnc, nthv, ntgv, ntrz, ntgz, nthz, &
                             ntsigma, ntomega, ntrac
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: clw
      real(kind=kind_phys), dimension(:,:,:), intent(in) :: dclw


      real(kind=kind_phys), dimension(:,:), intent(inout), optional :: cnvw_phy_f3d, cnvc_phy_f3d

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k, n, idtend, tracers

      real(kind=kind_phys), dimension(:,:), intent(in) :: ten_t, ten_u, ten_v
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: ten_q
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dudt, dvdt, dtdt 
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: dqdt
      real(kind=kind_phys), intent(in) ::  delt
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      !ten_q(:,:,1) already has a value from the deep convection scheme
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
      
      
      case_DCNV_ten: select case (tend_opt_dcnv)
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
          exit case_DCNV_ten
        case default
          errflg = 1
          errmsg = 'A tendency application control was outside of the acceptable range (1-4)'
          return
      end select case_DCNV_ten      
      
      if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
        tracers = 2
        do n=2,ntrac
!          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
!               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
!               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc) then
            IF ( otsptflag(n) ) THEN
            tracers = tracers + 1
            do k=1,levs
              do i=1,im
                clw(i,k,tracers) = gq0(i,k,n)
              enddo
            enddo
          endif
        enddo
      endif ! end if_ras or cfscnv or samf
      if (imp_physics == imp_physics_gfdl) then
        clw(1:im,:,1) = gq0(1:im,:,ntcw)
      elseif (imp_physics == imp_physics_thompson) then
        do k=1,levs
          do i=1,im
            clw(i,k,1)    = gq0(i,k,ntiw)                    ! ice
            clw(i,k,2)    = gq0(i,k,ntcw)                    ! water
          enddo
        enddo
      else if (imp_physics == imp_physics_nssl ) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntiw)                    ! cloud ice
            clw(i,k,2) = gq0(i,k,ntcw)                    ! cloud droplets
          enddo
        enddo
      elseif (imp_physics == imp_physics_wsm6 .or. imp_physics == imp_physics_mg .or. imp_physics == imp_physics_fer_hires) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = gq0(i,k,ntcw)                    ! water
          enddo
        enddo
      endif
      
      !shallow convection expects clw has already been updated
      
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
            dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_t * delt)*frain
          endif

          idtend=dtidx(index_of_x_wind,index_of_process_dcnv)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_u * delt)*frain
          endif

          idtend=dtidx(index_of_y_wind,index_of_process_dcnv)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_v * delt)*frain
          endif

          if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
             tracers = 2
             do n=2,ntrac
                if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
                     n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
                     n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc   .and. &
                     n /= nthl  .and. n /= nthnc .and. n /= nthv    .and. &
                     n /= ntrz  .and. n /= ntgz  .and. n /= nthz    .and. &
                     n /= ntgv  .and. n /= ntsigma .and. n /= ntomega) then
                   tracers = tracers + 1
                   idtend = dtidx(100+n,index_of_process_dcnv)
                   if(idtend>0) then
                      dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_q(:,:,n)*delt) * frain
                   endif
                endif
             enddo
          else
            do n=2,ntrac
               idtend = dtidx(100+n,index_of_process_dcnv)
               if(idtend>0) then
                  dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_q(:,:,n)*delt)*frain
               endif
            enddo
          endif
          idtend = dtidx(100+ntqv, index_of_process_dcnv)
          if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + (ten_q(:,:,ntqv)*delt) * frain
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

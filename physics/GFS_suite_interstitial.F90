!> \file GFS_suite_interstitial.f90
!!  Contains code related to more than one scheme in the GFS physics suite.

    module GFS_suite_interstitial_rad_reset

    contains

    subroutine GFS_suite_interstitial_rad_reset_init ()
    end subroutine GFS_suite_interstitial_rad_reset_init

    subroutine GFS_suite_interstitial_rad_reset_finalize()
    end subroutine GFS_suite_interstitial_rad_reset_finalize

!> \section arg_table_GFS_suite_interstitial_rad_reset_run Argument Table
!! \htmlinclude GFS_suite_interstitial_rad_reset_run.html
!!
    subroutine GFS_suite_interstitial_rad_reset_run (Interstitial, Model, errmsg, errflg)

      use machine,      only: kind_phys
      use GFS_typedefs, only: GFS_control_type, GFS_interstitial_type

      implicit none

      ! interface variables
      type(GFS_interstitial_type), intent(inout) :: Interstitial
      type(GFS_control_type),      intent(in)    :: Model
      character(len=*),            intent(out)   :: errmsg
      integer,                     intent(out)   :: errflg

      errmsg = ''
      errflg = 0

      call Interstitial%rad_reset(Model)

    end subroutine GFS_suite_interstitial_rad_reset_run

    end module GFS_suite_interstitial_rad_reset


    module GFS_suite_interstitial_phys_reset

    contains

    subroutine GFS_suite_interstitial_phys_reset_init ()
    end subroutine GFS_suite_interstitial_phys_reset_init

    subroutine GFS_suite_interstitial_phys_reset_finalize()
    end subroutine GFS_suite_interstitial_phys_reset_finalize

!> \section arg_table_GFS_suite_interstitial_phys_reset_run Argument Table
!! \htmlinclude GFS_suite_interstitial_phys_reset_run.html
!!
    subroutine GFS_suite_interstitial_phys_reset_run (Interstitial, Model, errmsg, errflg)

      use machine,      only: kind_phys
      use GFS_typedefs, only: GFS_control_type, GFS_interstitial_type

      implicit none

      ! interface variables
      type(GFS_interstitial_type), intent(inout) :: Interstitial
      type(GFS_control_type),      intent(in)    :: Model
      character(len=*),            intent(out)   :: errmsg
      integer,                     intent(out)   :: errflg

      errmsg = ''
      errflg = 0

      call Interstitial%phys_reset(Model)

    end subroutine GFS_suite_interstitial_phys_reset_run

    end module GFS_suite_interstitial_phys_reset


    module GFS_suite_interstitial_1

    contains

    subroutine GFS_suite_interstitial_1_init ()
    end subroutine GFS_suite_interstitial_1_init

    subroutine GFS_suite_interstitial_1_finalize()
    end subroutine GFS_suite_interstitial_1_finalize

!> \section arg_table_GFS_suite_interstitial_1_run Argument Table
!! \htmlinclude GFS_suite_interstitial_1_run.html
!!
    subroutine GFS_suite_interstitial_1_run (im, levs, ntrac, dtf, dtp, slmsk, area, dxmin, dxinv, pgr, &
      islmsk, work1, work2, psurf, dudt, dvdt, dtdt, dqdt, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! interface variables
      integer,              intent(in) :: im, levs, ntrac
      real(kind=kind_phys), intent(in) :: dtf, dtp, dxmin, dxinv
      real(kind=kind_phys), intent(in), dimension(im) :: slmsk, area, pgr

      integer,              intent(out), dimension(im) :: islmsk
      real(kind=kind_phys), intent(out), dimension(im) :: work1, work2, psurf
      real(kind=kind_phys), intent(out), dimension(im,levs) :: dudt, dvdt, dtdt
      real(kind=kind_phys), intent(out), dimension(im,levs,ntrac) ::  dqdt
      real(kind=kind_phys), parameter   :: zero = 0.0_kind_phys, one = 1.0_kind_phys
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! local variables
      integer :: i, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i = 1, im
        islmsk(i)   = nint(slmsk(i))

        work1(i) = (log(area(i)) - dxmin) * dxinv
        work1(i) = max(zero, min(one, work1(i)))
        work2(i) = one - work1(i)
        psurf(i) = pgr(i)
      end do

      do k=1,levs
        do i=1,im
          dudt(i,k)  = zero
          dvdt(i,k)  = zero
          dtdt(i,k)  = zero
        enddo
      enddo
      do n=1,ntrac
        do k=1,levs
          do i=1,im
            dqdt(i,k,n) = zero
          enddo
        enddo
      enddo

    end subroutine GFS_suite_interstitial_1_run

  end module GFS_suite_interstitial_1


  module GFS_suite_interstitial_2

  use machine, only: kind_phys
  real(kind=kind_phys), parameter :: one = 1.0_kind_phys
  logical :: linit_mod  = .false. 

  contains

    subroutine GFS_suite_interstitial_2_init ()
    end subroutine GFS_suite_interstitial_2_init

    subroutine GFS_suite_interstitial_2_finalize()
    end subroutine GFS_suite_interstitial_2_finalize

!> \section arg_table_GFS_suite_interstitial_2_run Argument Table
!! \htmlinclude GFS_suite_interstitial_2_run.html
!!
    subroutine GFS_suite_interstitial_2_run (im, levs, lssav, ldiag3d, lsidea, cplflx, flag_cice, shal_cnv, old_monin, mstrat,    &
      do_shoc, frac_grid, imfshalcnv, dtf, xcosz, adjsfcdsw, adjsfcdlw, cice, pgr, ulwsfc_cice, lwhd, htrsw, htrlw, xmu, ctei_rm, &
      work1, work2, prsi, tgrs, prsl, qgrs_water_vapor, qgrs_cloud_water, cp, hvap, prslk, suntim, adjsfculw, adjsfculw_lnd,      &
      adjsfculw_ice, adjsfculw_wat, dlwsfc, ulwsfc, psmean, dtend, dtidx, index_of_process_longwave, index_of_process_shortwave,  &
      index_of_process_pbl, index_of_process_dcnv, index_of_process_scnv, index_of_process_mp, index_of_temperature,              &
      ctei_rml, ctei_r, kinver, dry, icy, wet, frland, huge, use_LW_jacobian, errmsg, errflg)

      implicit none

      ! interface variables
      integer,              intent(in   ) :: im, levs, imfshalcnv
      logical,              intent(in   ) :: lssav, ldiag3d, lsidea, cplflx, shal_cnv
      logical,              intent(in   ) :: old_monin, mstrat, do_shoc, frac_grid, use_LW_jacobian
      real(kind=kind_phys), intent(in   ) :: dtf, cp, hvap

      logical,              intent(in   ), dimension(im) :: flag_cice
      real(kind=kind_phys), intent(in   ), dimension(2) :: ctei_rm
      real(kind=kind_phys), intent(in   ), dimension(im) :: xcosz, adjsfcdsw, adjsfcdlw, pgr, xmu, ulwsfc_cice, work1, work2
      real(kind=kind_phys), intent(in   ), dimension(im) :: cice
      real(kind=kind_phys), intent(in   ), dimension(im, levs) :: htrsw, htrlw, tgrs, prsl, qgrs_water_vapor, qgrs_cloud_water, prslk
      real(kind=kind_phys), intent(in   ), dimension(im, levs+1) :: prsi
      real(kind=kind_phys), intent(in   ), dimension(im, levs, 6) :: lwhd
      integer,              intent(inout), dimension(im) :: kinver
      real(kind=kind_phys), intent(inout), dimension(im) :: suntim, dlwsfc, ulwsfc, psmean, ctei_rml, ctei_r
      real(kind=kind_phys), intent(in   ), dimension(im) :: adjsfculw_lnd, adjsfculw_ice, adjsfculw_wat
      real(kind=kind_phys), intent(inout), dimension(im) :: adjsfculw

      ! dtend is only allocated if ldiag3d is .true.
      real(kind=kind_phys), optional, intent(inout), dimension(:,:,:) :: dtend
      integer,              intent(in),    dimension(:,:) :: dtidx
      integer, intent(in) :: index_of_process_longwave, index_of_process_shortwave, &
           index_of_process_pbl, index_of_process_dcnv, index_of_process_scnv,       &
           index_of_process_mp, index_of_temperature

      logical,              intent(in   ), dimension(im) :: dry, icy, wet
      real(kind=kind_phys), intent(in   ), dimension(im) :: frland
      real(kind=kind_phys), intent(in   ) :: huge

      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! local variables
      real(kind=kind_phys), parameter :: czmin   = 0.0001_kind_phys      ! cos(89.994)
      integer :: i, k, idtend
      real(kind=kind_phys) :: tem1, tem2, tem, hocp
      logical, dimension(im) :: invrsn
      real(kind=kind_phys), dimension(im) :: tx1, tx2

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys
      real(kind=kind_phys), parameter :: qmin = 1.0e-10_kind_phys, epsln=1.0e-10_kind_phys

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      hocp = hvap/cp

      if (lssav) then      !  --- ...  accumulate/save output variables

!  --- ...  sunshine duration time is defined as the length of time (in mdl output
!           interval) that solar radiation falling on a plane perpendicular to the
!           direction of the sun >= 120 w/m2

        do i = 1, im
          if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
            tem1 = adjsfcdsw(i) / xcosz(i)
            if ( tem1 >= 120.0_kind_phys ) then
              suntim(i) = suntim(i) + dtf
            endif
          endif
        enddo

!  --- ...  sfc lw fluxes used by atmospheric model are saved for output
        if (.not. use_LW_jacobian) then
        if (frac_grid) then
           do i=1,im
              tem = (one - frland(i)) * cice(i) ! tem = ice fraction wrt whole cell
              if (flag_cice(i)) then
                 adjsfculw(i) = adjsfculw_lnd(i) * frland(i)               &
                              + ulwsfc_cice(i)   * tem                     &
                              + adjsfculw_wat(i) * (one - frland(i) - tem)
              else
                 adjsfculw(i) = adjsfculw_lnd(i) * frland(i)               &
                              + adjsfculw_ice(i) * tem                     &
                              + adjsfculw_wat(i) * (one - frland(i) - tem)
              endif
           enddo
        else
           do i=1,im
              if (dry(i)) then                     ! all land
                 adjsfculw(i) = adjsfculw_lnd(i)
              elseif (icy(i)) then                 ! ice (and water)
                 tem = one - cice(i)
                 if (flag_cice(i)) then
                    if (wet(i) .and. adjsfculw_wat(i) /= huge) then
                       adjsfculw(i) = ulwsfc_cice(i)*cice(i) + adjsfculw_wat(i)*tem
                    else
                       adjsfculw(i) = ulwsfc_cice(i)
                    endif
                 else
                    if (wet(i) .and. adjsfculw_wat(i) /= huge) then
                       adjsfculw(i) = adjsfculw_ice(i)*cice(i) + adjsfculw_wat(i)*tem
                    else
                       adjsfculw(i) = adjsfculw_ice(i)
                    endif
                 endif
              else                                 ! all water
                 adjsfculw(i) = adjsfculw_wat(i)
              endif
           enddo
        endif
        endif

        do i=1,im
          dlwsfc(i) = dlwsfc(i) + adjsfcdlw(i)*dtf
          ulwsfc(i) = ulwsfc(i) + adjsfculw(i)*dtf
          psmean(i) = psmean(i) + pgr(i)*dtf        ! mean surface pressure
        end do

        if (ldiag3d) then
          if (lsidea) then
            idtend = dtidx(index_of_temperature,index_of_process_longwave)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + lwhd(:,:,1)*dtf
            endif

            idtend = dtidx(index_of_temperature,index_of_process_shortwave)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + lwhd(:,:,2)*dtf
            endif

            idtend = dtidx(index_of_temperature,index_of_process_pbl)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + lwhd(:,:,3)*dtf
            endif

            idtend = dtidx(index_of_temperature,index_of_process_dcnv)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + lwhd(:,:,4)*dtf
            endif

            idtend = dtidx(index_of_temperature,index_of_process_scnv)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + lwhd(:,:,5)*dtf
            endif

            idtend = dtidx(index_of_temperature,index_of_process_mp)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + lwhd(:,:,6)*dtf
            endif
          else
            idtend = dtidx(index_of_temperature,index_of_process_longwave)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + htrlw(:,:)*dtf
            endif

            idtend = dtidx(index_of_temperature,index_of_process_shortwave)
            if(idtend>=1) then
               do k=1,levs
                  do i=1,im
                     dtend(i,k,idtend) = dtend(i,k,idtend) + htrsw(i,k)*dtf*xmu(i)
                  enddo
               enddo
            endif
          endif
        endif
      endif    ! end if_lssav_block

      do i=1, im
        invrsn(i) = .false.
        tx1(i)    = zero
        tx2(i)    = 10.0_kind_phys
        ctei_r(i) = 10.0_kind_phys
      enddo

      if ((((imfshalcnv == 0 .and. shal_cnv) .or. old_monin) .and. mstrat) &
         .or. do_shoc) then
        ctei_rml(:) = ctei_rm(1)*work1(:) + ctei_rm(2)*work2(:)
        do k=1,levs/2
          do i=1,im
            if (prsi(i,1)-prsi(i,k+1) < 0.35_kind_phys*prsi(i,1)       &
                .and. (.not. invrsn(i))) then
              tem = (tgrs(i,k+1) - tgrs(i,k))  &
                  / (prsl(i,k)   - prsl(i,k+1))

              if (((tem > 0.0001_kind_phys) .and. (tx1(i) < zero)) .or.  &
                  ((tem-abs(tx1(i)) > zero) .and. (tx2(i) < zero))) then
                invrsn(i) = .true.

                if (qgrs_water_vapor(i,k) > qgrs_water_vapor(i,k+1)) then
                  tem1 = tgrs(i,k+1) + hocp*max(qgrs_water_vapor(i,k+1),qmin)
                  tem2 = tgrs(i,k)   + hocp*max(qgrs_water_vapor(i,k),qmin)

                  tem1 = tem1 / prslk(i,k+1) - tem2 / prslk(i,k)

!  --- ...  (cp/l)(deltathetae)/(deltatwater) > ctei_rm -> conditon for CTEI
                  ctei_r(i) = (one/hocp)*tem1/(qgrs_water_vapor(i,k+1)-qgrs_water_vapor(i,k)  &
                            + qgrs_cloud_water(i,k+1)-qgrs_cloud_water(i,k))
                else
                  ctei_r(i) = 10.0_kind_phys
                endif

                if ( ctei_rml(i) > ctei_r(i) ) then
                  kinver(i) = k
                else
                  kinver(i) = levs
                endif
              endif

              tx2(i) = tx1(i)
              tx1(i) = tem
            endif
          enddo
        enddo
      endif

    end subroutine GFS_suite_interstitial_2_run

  end module GFS_suite_interstitial_2


  module GFS_suite_stateout_reset

  contains

    subroutine GFS_suite_stateout_reset_init ()
    end subroutine GFS_suite_stateout_reset_init

    subroutine GFS_suite_stateout_reset_finalize()
    end subroutine GFS_suite_stateout_reset_finalize

!> \section arg_table_GFS_suite_stateout_reset_run Argument Table
!! \htmlinclude GFS_suite_stateout_reset_run.html
!!
    subroutine GFS_suite_stateout_reset_run (im, levs, ntrac,        &
                                             tgrs, ugrs, vgrs, qgrs, &
                                             gt0 , gu0 , gv0 , gq0 , &
                                             errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! interface variables
      integer, intent(in) :: im
      integer, intent(in) :: levs
      integer, intent(in) :: ntrac
      real(kind=kind_phys), dimension(im,levs),       intent(in)  :: tgrs, ugrs, vgrs
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(in)  :: qgrs
      real(kind=kind_phys), dimension(im,levs),       intent(out) :: gt0, gu0, gv0
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(out) :: gq0

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      gt0(:,:)   = tgrs(:,:)
      gu0(:,:)   = ugrs(:,:)
      gv0(:,:)   = vgrs(:,:)
      gq0(:,:,:) = qgrs(:,:,:)

    end subroutine GFS_suite_stateout_reset_run

  end module GFS_suite_stateout_reset


  module GFS_suite_stateout_update

  contains

    subroutine GFS_suite_stateout_update_init ()
    end subroutine GFS_suite_stateout_update_init

    subroutine GFS_suite_stateout_update_finalize()
    end subroutine GFS_suite_stateout_update_finalize

!> \section arg_table_GFS_suite_stateout_update_run Argument Table
!! \htmlinclude GFS_suite_stateout_update_run.html
!!
    subroutine GFS_suite_stateout_update_run (im, levs, ntrac, dtp,  &
                     tgrs, ugrs, vgrs, qgrs, dudt, dvdt, dtdt, dqdt, &
                     gt0, gu0, gv0, gq0, ntiw, nqrimef, imp_physics, &
                     imp_physics_fer_hires, epsq, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! Interface variables
      integer,              intent(in) :: im
      integer,              intent(in) :: levs
      integer,              intent(in) :: ntrac
      integer,              intent(in) :: imp_physics,imp_physics_fer_hires
      integer,              intent(in) :: ntiw, nqrimef
      real(kind=kind_phys), intent(in) :: dtp, epsq

      real(kind=kind_phys), dimension(im,levs),       intent(in)  :: tgrs, ugrs, vgrs
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(in)  :: qgrs
      real(kind=kind_phys), dimension(im,levs),       intent(in)  :: dudt, dvdt, dtdt
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(in)  :: dqdt
      real(kind=kind_phys), dimension(im,levs),       intent(out) :: gt0, gu0, gv0
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(out) :: gq0

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer                       :: i, k
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      gt0(:,:)   = tgrs(:,:)   + dtdt(:,:)   * dtp
      gu0(:,:)   = ugrs(:,:)   + dudt(:,:)   * dtp
      gv0(:,:)   = vgrs(:,:)   + dvdt(:,:)   * dtp
      gq0(:,:,:) = qgrs(:,:,:) + dqdt(:,:,:) * dtp
      
      if (imp_physics == imp_physics_fer_hires) then
       do k=1,levs
         do i=1,im
           if(gq0(i,k,ntiw) > epsq) then
             gq0(i,k,nqrimef) = max(1., gq0(i,k,nqrimef)/gq0(i,k,ntiw))
           else
             gq0(i,k,nqrimef) = 1.
           end if
         end do
       end do
      end if

    end subroutine GFS_suite_stateout_update_run

  end module GFS_suite_stateout_update


  module GFS_suite_interstitial_3

  contains

    subroutine GFS_suite_interstitial_3_init ()
    end subroutine GFS_suite_interstitial_3_init

    subroutine GFS_suite_interstitial_3_finalize()
    end subroutine GFS_suite_interstitial_3_finalize

!> \section arg_table_GFS_suite_interstitial_3_run Argument Table
!! \htmlinclude GFS_suite_interstitial_3_run.html
!!
    subroutine GFS_suite_interstitial_3_run (im, levs, nn, cscnv,       &
               satmedmf, trans_trac, do_shoc, ltaerosol, ntrac, ntcw,   &
               ntiw, ntclamt, ntrw, ntsw, ntrnc, ntsnc, ntgl, ntgnc,    &
               xlon, xlat, gt0, gq0, imp_physics, imp_physics_mg,       &
               imp_physics_zhao_carr, imp_physics_zhao_carr_pdf,        &
               imp_physics_gfdl, imp_physics_thompson,                  &
               imp_physics_wsm6, imp_physics_fer_hires, prsi,           &
               prsl, prslk, rhcbot,rhcpbl, rhctop, rhcmax, islmsk,      &
               work1, work2, kpbl, kinver, ras, me,                     &
               clw, rhc, save_qc, save_qi, save_tcp, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! interface variables
      integer,                                          intent(in) :: im, levs, nn, ntrac, ntcw, ntiw, ntclamt, ntrw,   &
        ntsw, ntrnc, ntsnc, ntgl, ntgnc, imp_physics, imp_physics_mg, imp_physics_zhao_carr, imp_physics_zhao_carr_pdf, &
        imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,imp_physics_fer_hires, me
      integer, dimension(im),                           intent(in) :: islmsk, kpbl, kinver
      logical,                                          intent(in) :: cscnv, satmedmf, trans_trac, do_shoc, ltaerosol, ras

      real(kind=kind_phys),                             intent(in) :: rhcbot, rhcmax, rhcpbl, rhctop
      real(kind=kind_phys), dimension(im),              intent(in) :: work1, work2
      real(kind=kind_phys), dimension(im, levs),        intent(in) :: prsl, prslk
      real(kind=kind_phys), dimension(im, levs+1),      intent(in) :: prsi
      real(kind=kind_phys), dimension(im),              intent(in) :: xlon, xlat
      real(kind=kind_phys), dimension(im, levs),        intent(in) :: gt0
      real(kind=kind_phys), dimension(im, levs, ntrac), intent(in) :: gq0

      real(kind=kind_phys), dimension(im, levs),      intent(inout) :: rhc, save_qc
      ! save_qi is not allocated for Zhao-Carr MP
      real(kind=kind_phys), dimension(:, :),          intent(inout) :: save_qi
      real(kind=kind_phys), dimension(:, :),          intent(inout) :: save_tcp
      real(kind=kind_phys), dimension(im, levs, nn),  intent(inout) :: clw

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! local variables
      integer :: i,k,n,tracers,kk
      real(kind=kind_phys) :: tem, tem1, tem2
      real(kind=kind_phys), dimension(im) :: tx1, tx2, tx3, tx4

      !real(kind=kind_phys),parameter :: slope_mg = 0.02, slope_upmg = 0.04,  &
      !                   turnrhcrit = 0.900, turnrhcrit_upper = 0.150
      ! in the following inverse of slope_mg and slope_upmg are specified
      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys
      real(kind=kind_phys), parameter :: slope_mg   = 50.0_kind_phys,   &
                                         slope_upmg = 25.0_kind_phys

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
        tracers = 2
        do n=2,ntrac
          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc) then
            tracers = tracers + 1
            do k=1,levs
              do i=1,im
                clw(i,k,tracers) = gq0(i,k,n)
              enddo
            enddo
          endif
        enddo
      endif ! end if_ras or cfscnv or samf

      if (ntcw > 0) then
        if (imp_physics == imp_physics_mg .and. rhcpbl < 0.5_kind_phys) then ! compute rhc for GMAO macro physics cloud pdf
          do i=1,im
            tx1(i) = one / prsi(i,1)
            tx2(i) = one - rhcmax*work1(i)-rhcbot*work2(i)

            kk     = min(kinver(i), max(2,kpbl(i)))
            tx3(i) = prsi(i,kk)*tx1(i)
            tx4(i) = rhcpbl - rhctop*abs(cos(xlat(i)))
          enddo
          do k = 1, levs
            do i = 1, im
              tem  = prsl(i,k) * tx1(i)
              tem1 = min(max((tem-tx3(i))*slope_mg, -20.0_kind_phys), 20.0_kind_phys)
              ! Using rhcpbl and rhctop from the namelist instead of 0.3 and 0.2
              ! and rhcbot represents pbl top critical relative humidity
              tem2 = min(max((tx4(i)-tem)*slope_upmg, -20.0_kind_phys), 20.0_kind_phys) ! Anning
              if (islmsk(i) > 0) then
                tem1 = one / (one+exp(tem1+tem1))
              else
                tem1 = 2.0_kind_phys / (one+exp(tem1+tem1))
              endif
              tem2 = one / (one+exp(tem2))

              rhc(i,k) = min(rhcmax, max(0.7_kind_phys, one-tx2(i)*tem1*tem2))
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              kk = max(10,kpbl(i))
              if (k < kk) then
                tem    = rhcbot - (rhcbot-rhcpbl) * (one-prslk(i,k)) / (one-prslk(i,kk))
              else
                tem    = rhcpbl - (rhcpbl-rhctop) * (prslk(i,kk)-prslk(i,k)) / prslk(i,kk)
              endif
              tem      = rhcmax * work1(i) + tem * work2(i)
              rhc(i,k) = max(zero, min(one,tem))
            enddo
          enddo
        endif
      else
        rhc(:,:) = 1.0
      endif

      if (imp_physics == imp_physics_zhao_carr .or. imp_physics == imp_physics_zhao_carr_pdf) then   ! zhao-carr microphysics
        !GF* move to GFS_MP_generic_pre (from gscond/precpd)
        ! do i=1,im
        !   psautco_l(i) = Model%psautco(1)*work1(i) + Model%psautco(2)*work2(i)
        !   prautco_l(i) = Model%prautco(1)*work1(i) + Model%prautco(2)*work2(i)
        ! enddo
        !*GF
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntcw)
          enddo
        enddo
      elseif (imp_physics == imp_physics_gfdl) then
        clw(1:im,:,1) = gq0(1:im,:,ntcw)
      elseif (imp_physics == imp_physics_thompson) then
        do k=1,levs
          do i=1,im
            clw(i,k,1)    = gq0(i,k,ntiw)                    ! ice
            clw(i,k,2)    = gq0(i,k,ntcw)                    ! water
            save_tcp(i,k) = gt0(i,k)
          enddo
        enddo
        if(ltaerosol) then
          save_qi(:,:) = clw(:,:,1)
          save_qc(:,:) = clw(:,:,2)
        else
          save_qi(:,:) = clw(:,:,1)
        endif
      elseif (imp_physics == imp_physics_wsm6 .or. imp_physics == imp_physics_mg .or. imp_physics == imp_physics_fer_hires) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = gq0(i,k,ntcw)                    ! water
          enddo
        enddo
      endif

    end subroutine GFS_suite_interstitial_3_run

  end module GFS_suite_interstitial_3

  module GFS_suite_interstitial_4

  contains

    subroutine GFS_suite_interstitial_4_init ()
    end subroutine GFS_suite_interstitial_4_init

    subroutine GFS_suite_interstitial_4_finalize()
    end subroutine GFS_suite_interstitial_4_finalize

!> \section arg_table_GFS_suite_interstitial_4_run Argument Table
!! \htmlinclude GFS_suite_interstitial_4_run.html
!!
    subroutine GFS_suite_interstitial_4_run (im, levs, ltaerosol, cplchm, tracers_total, ntrac, ntcw, ntiw, ntclamt, &
      ntrw, ntsw, ntrnc, ntsnc, ntgl, ntgnc, ntlnc, ntinc, nn, imp_physics, imp_physics_gfdl, imp_physics_thompson,  &
      imp_physics_zhao_carr, imp_physics_zhao_carr_pdf, convert_dry_rho, dtf, save_qc, save_qi, con_pi, dtidx, dtend,&
      index_of_process_conv_trans, gq0, clw, prsl, save_tcp, con_rd, con_eps, nwfa, spechum, dqdti, ldiag3d,         &
      ntk, ntke, errmsg, errflg)

      use machine,               only: kind_phys
      use module_mp_thompson_make_number_concentrations, only: make_IceNumber, make_DropletNumber

      implicit none

      ! interface variables

      integer,                                  intent(in) :: im, levs, tracers_total, ntrac, ntcw, ntiw, ntclamt, ntrw,  &
        ntsw, ntrnc, ntsnc, ntgl, ntgnc, ntlnc, ntinc, nn, imp_physics, imp_physics_gfdl, imp_physics_thompson,           &
        imp_physics_zhao_carr, imp_physics_zhao_carr_pdf

      logical,                                  intent(in) :: ltaerosol, cplchm, convert_dry_rho

      real(kind=kind_phys),                     intent(in) :: con_pi, dtf
      real(kind=kind_phys), dimension(im,levs), intent(in) :: save_qc
      ! save_qi is not allocated for Zhao-Carr MP
      real(kind=kind_phys), dimension(:, :),    intent(in) :: save_qi

      ! dtend and dtidx are only allocated if ldiag3d
      logical, intent(in)                                     :: ldiag3d
      real(kind=kind_phys), dimension(:,:,:),   intent(inout) :: dtend
      integer,              dimension(:,:),     intent(in)    :: dtidx
      integer,                                  intent(in)    :: index_of_process_conv_trans,ntk,ntke

      real(kind=kind_phys), dimension(im,levs,ntrac), intent(inout) :: gq0
      real(kind=kind_phys), dimension(im,levs,nn),    intent(inout) :: clw
      real(kind=kind_phys), dimension(im,levs),       intent(in) :: prsl
      real(kind=kind_phys),                           intent(in) :: con_rd, con_eps
      real(kind=kind_phys), dimension(:,:),           intent(in) :: nwfa, save_tcp
      real(kind=kind_phys), dimension(im,levs),       intent(in) :: spechum

      ! dqdti may not be allocated
      real(kind=kind_phys), dimension(:,:),           intent(inout) :: dqdti

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! local variables
      integer :: i,k,n,tracers,idtend

      real(kind=kind_phys) :: rho, orho
      real(kind=kind_phys), dimension(im,levs) :: qv_mp !< kg kg-1 (dry mixing ratio)
      real(kind=kind_phys), dimension(im,levs) :: qc_mp !< kg kg-1 (dry mixing ratio)
      real(kind=kind_phys), dimension(im,levs) :: qi_mp !< kg kg-1 (dry mixing ratio)
      real(kind=kind_phys), dimension(im,levs) :: nc_mp !< kg-1 (dry mixing ratio)
      real(kind=kind_phys), dimension(im,levs) :: ni_mp !< kg-1 (dry mixing ratio)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if(ldiag3d) then
         if(ntk>0 .and. ntk<=size(clw,3)) then
            idtend=dtidx(100+ntke,index_of_process_conv_trans)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,ntk)
            endif
         endif
         if(ntclamt<=size(clw,3) .and. ntclamt>0) then
            idtend=dtidx(100+ntclamt,index_of_process_conv_trans)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,ntclamt)
            endif
         endif
      endif

      if(ldiag3d .and. ntk>0) then
         idtend=dtidx(100+ntke,index_of_process_conv_trans)
         if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,ntk)
         endif
      endif

!  --- update the tracers due to deep & shallow cumulus convective transport
!           (except for suspended water and ice)

      if (tracers_total > 0) then
        tracers = 2
        do n=2,ntrac
!         if ( n /= ntcw .and. n /= ntiw .and. n /= ntclamt) then
          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc ) then
              tracers = tracers + 1
            do k=1,levs
              do i=1,im
                gq0(i,k,n) = clw(i,k,tracers)
              enddo
            enddo
          endif
        enddo
      endif

      if (ntcw > 0) then

!  for microphysics
        if (imp_physics == imp_physics_zhao_carr     .or. &
            imp_physics == imp_physics_zhao_carr_pdf .or. &
            imp_physics == imp_physics_gfdl) then
           gq0(1:im,:,ntcw) = clw(1:im,:,1) + clw(1:im,:,2)

        elseif (ntiw > 0) then
          do k=1,levs
            do i=1,im
              gq0(i,k,ntiw) = clw(i,k,1)                     ! ice
              gq0(i,k,ntcw) = clw(i,k,2)                     ! water
            enddo
          enddo

          if (imp_physics == imp_physics_thompson .and. (ntlnc>0 .or. ntinc>0)) then
            if_convert_dry_rho: if (convert_dry_rho) then
              do k=1,levs
                do i=1,im
                  !> - Convert specific humidity to dry mixing ratio
                  qv_mp(i,k) = spechum(i,k) / (one-spechum(i,k))
                  !> - Density of air in kg m-3 and inverse density
                  rho = con_eps*prsl(i,k) / (con_rd*save_tcp(i,k)*(qv_mp(i,k)+con_eps))
                  orho = one/rho
                  if (ntlnc>0) then
                    !> - Convert moist mixing ratio to dry mixing ratio
                    qc_mp(i,k) = (clw(i,k,2)-save_qc(i,k)) / (one-spechum(i,k))
                    !> - Convert number concentration from moist to dry
                    nc_mp(i,k) = gq0(i,k,ntlnc) / (one-spechum(i,k))
                    nc_mp(i,k) = max(zero, nc_mp(i,k) + make_DropletNumber(qc_mp(i,k) * rho, nwfa(i,k)*rho) * orho)
                    !> - Convert number concentrations from dry to moist
                    gq0(i,k,ntlnc) = nc_mp(i,k) / (one+qv_mp(i,k))
                  endif
                  if (ntinc>0) then
                    !> - Convert moist mixing ratio to dry mixing ratio
                    qi_mp(i,k) = (clw(i,k,1)-save_qi(i,k)) / (one-spechum(i,k))
                    !> - Convert number concentration from moist to dry
                    ni_mp(i,k) = gq0(i,k,ntinc) / (one-spechum(i,k)) 
                    ni_mp(i,k) = max(zero, ni_mp(i,k) + make_IceNumber(qi_mp(i,k) * rho, save_tcp(i,k)) * orho)
                    !> - Convert number concentrations from dry to moist
                    gq0(i,k,ntinc) = ni_mp(i,k) / (one+qv_mp(i,k))
                  endif
                enddo
              enddo
            else
              do k=1,levs
                do i=1,im
                  !> - Density of air in kg m-3 and inverse density
                  rho = con_eps*prsl(i,k) / (con_rd*save_tcp(i,k)*(spechum(i,k)+con_eps))
                  orho = one/rho
                  if (ntlnc>0) then
                    !> - Update cloud water mixing ratio
                    qc_mp(i,k) = (clw(i,k,2)-save_qc(i,k))
                    !> - Update cloud water number concentration
                    gq0(i,k,ntlnc) = max(zero, gq0(i,k,ntlnc) + make_DropletNumber(qc_mp(i,k) * rho, nwfa(i,k)*rho) * orho)
                  endif
                  if (ntinc>0) then
                    !> - Update cloud ice mixing ratio
                    qi_mp(i,k) = (clw(i,k,1)-save_qi(i,k))
                    !> - Update cloud ice number concentration
                    gq0(i,k,ntinc) = max(zero, gq0(i,k,ntinc) + make_IceNumber(qi_mp(i,k) * rho, save_tcp(i,k)) * orho)
                  endif
                enddo
              enddo
            end if if_convert_dry_rho
          endif

        else
          do k=1,levs
            do i=1,im
              gq0(i,k,ntcw) = clw(i,k,1) + clw(i,k,2)
            enddo
          enddo
        endif   ! end if_ntiw

      else
        do k=1,levs
          do i=1,im
            clw(i,k,1) = clw(i,k,1) + clw(i,k,2)
          enddo
        enddo
      endif   ! end if_ntcw

! dqdt_v : instaneous moisture tendency (kg/kg/sec)
      if (cplchm) then
        do k=1,levs
          do i=1,im
            dqdti(i,k) = dqdti(i,k) * (one / dtf)
          enddo
        enddo
      endif

    end subroutine GFS_suite_interstitial_4_run

  end module GFS_suite_interstitial_4

  module GFS_suite_interstitial_5

  contains

    subroutine GFS_suite_interstitial_5_init ()
    end subroutine GFS_suite_interstitial_5_init

    subroutine GFS_suite_interstitial_5_finalize()
    end subroutine GFS_suite_interstitial_5_finalize

!> \section arg_table_GFS_suite_interstitial_5_run Argument Table
!! \htmlinclude GFS_suite_interstitial_5_run.html
!!
    subroutine GFS_suite_interstitial_5_run (im, levs, ntrac, ntcw, ntiw, nn, gq0, clw, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! interface variables
      integer,                                          intent(in)  :: im, levs, ntrac, ntcw, ntiw, nn

      real(kind=kind_phys), dimension(im, levs, ntrac), intent(in)  :: gq0

      real(kind=kind_phys), dimension(im, levs, nn),    intent(out) :: clw

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! local variables
      integer :: i,k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do k=1,levs
        do i=1,im
          clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
          clw(i,k,2) = gq0(i,k,ntcw)                    ! water
        enddo
      enddo

    end subroutine GFS_suite_interstitial_5_run

  end module GFS_suite_interstitial_5


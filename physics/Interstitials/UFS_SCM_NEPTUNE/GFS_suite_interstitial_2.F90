!> \file GFS_suite_interstitial_2.f90
!!  Contains code related used to calculate radiation-based and PBL-based diagnostics that are executed after radiation time interpolation and before the surface layer. 

  module GFS_suite_interstitial_2

  use machine, only: kind_phys
  real(kind=kind_phys), parameter :: one = 1.0_kind_phys
  logical :: linit_mod  = .false. 

  contains

!> \section arg_table_GFS_suite_interstitial_2_run Argument Table
!! \htmlinclude GFS_suite_interstitial_2_run.html
!!
    subroutine GFS_suite_interstitial_2_run (im, levs, lssav, ldiag3d, lsidea, flag_cice, shal_cnv, old_monin, mstrat,            &
      do_shoc, frac_grid, imfshalcnv, dtf, xcosz, adjsfcdsw, adjsfcdlw, cice, pgr, ulwsfc_cice, lwhd, htrsw, htrlw, xmu, ctei_rm, &
      work1, work2, prsi, tgrs, prsl, qgrs_water_vapor, qgrs_cloud_water, cp, hvap, prslk, suntim, adjsfculw, adjsfculw_lnd,      &
      adjsfculw_ice, adjsfculw_wat, dlwsfc, ulwsfc, psmean, dtend, dtidx, index_of_process_longwave, index_of_process_shortwave,  &
      index_of_process_pbl, index_of_process_dcnv, index_of_process_scnv, index_of_process_mp, index_of_temperature,              &
      ctei_rml, ctei_r, kinver, dry, icy, wet, frland, huge, use_LW_jacobian, htrlwu, errmsg, errflg)

      implicit none

      ! interface variables
      integer,              intent(in   )                   :: im, levs, imfshalcnv
      logical,              intent(in   )                   :: lssav, ldiag3d, lsidea, shal_cnv
      logical,              intent(in   )                   :: old_monin, mstrat, do_shoc, frac_grid, use_LW_jacobian
      real(kind=kind_phys), intent(in   )                   :: dtf, cp, hvap

      logical,              intent(in   ), dimension(:)     :: flag_cice
      real(kind=kind_phys), intent(in   ), dimension(:)     :: ctei_rm
      real(kind=kind_phys), intent(in   ), dimension(:)     :: xcosz, adjsfcdsw, adjsfcdlw, pgr, xmu, work1, work2
      real(kind=kind_phys), intent(in   ), dimension(:), optional :: ulwsfc_cice
      real(kind=kind_phys), intent(in   ), dimension(:)     :: cice
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: htrsw, htrlw, tgrs, prsl, qgrs_water_vapor, qgrs_cloud_water, prslk
      real(kind=kind_phys), intent(in   ), dimension(:,:), optional :: htrlwu
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: prsi
      real(kind=kind_phys), intent(in   ), dimension(:,:,:) :: lwhd
      integer,              intent(inout), dimension(:)     :: kinver
      real(kind=kind_phys), intent(inout), dimension(:)     :: suntim, dlwsfc, ulwsfc, psmean, ctei_rml, ctei_r
      real(kind=kind_phys), intent(in   ), dimension(:)     :: adjsfculw_lnd, adjsfculw_ice, adjsfculw_wat
      real(kind=kind_phys), intent(inout), dimension(:)     :: adjsfculw

      ! dtend is only allocated if ldiag3d is .true.
      real(kind=kind_phys), optional, intent(inout), dimension(:,:,:) :: dtend
      integer,              intent(in),    dimension(:,:) :: dtidx
      integer, intent(in) :: index_of_process_longwave, index_of_process_shortwave, &
           index_of_process_pbl, index_of_process_dcnv, index_of_process_scnv,       &
           index_of_process_mp, index_of_temperature

      logical,              intent(in   ), dimension(:)     :: dry, icy, wet
      real(kind=kind_phys), intent(in   ), dimension(:)     :: frland
      real(kind=kind_phys), intent(in   )                   :: huge

      character(len=*),     intent(  out)                   :: errmsg
      integer,              intent(  out)                   :: errflg

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

        do i=1,im
          dlwsfc(i) = dlwsfc(i) + adjsfcdlw(i)*dtf
          ulwsfc(i) = ulwsfc(i) + adjsfculw(i)*dtf
          psmean(i) = psmean(i) + pgr(i)*dtf        ! mean surface pressure
        enddo

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
              if (use_LW_jacobian) then
                dtend(:,:,idtend) = dtend(:,:,idtend) + htrlwu(:,:)*dtf
              else
                dtend(:,:,idtend) = dtend(:,:,idtend) + htrlw(:,:)*dtf
              endif
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
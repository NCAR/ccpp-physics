!> \file GFS_surface_generic.F90
!!  Contains code related to all GFS surface schemes.

      module GFS_surface_generic_pre

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_pre_init, GFS_surface_generic_pre_finalize, GFS_surface_generic_pre_run

      real(kind=kind_phys), parameter :: one = 1.0d0
      real(kind=kind_phys), parameter :: zero = 0.0d0

      contains

      subroutine GFS_surface_generic_pre_init ()
      end subroutine GFS_surface_generic_pre_init

      subroutine GFS_surface_generic_pre_finalize()
      end subroutine GFS_surface_generic_pre_finalize

!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! \htmlinclude GFS_surface_generic_pre_run.html
!!
      subroutine GFS_surface_generic_pre_run (im, levs, vfrac, islmsk, isot, ivegsrc, stype, vtype, slope, &
                          prsik_1, prslk_1, tsfc, phil, con_g,                                             &
                          sigmaf, soiltyp, vegtype, slopetyp, work3, tsurf, zlvl, do_sppt, dtdtr,          &
                          drain_cpl, dsnow_cpl, rain_cpl, snow_cpl, do_sfcperts, nsfcpert, sfc_wts,        &
                          pertz0, pertzt, pertshc, pertlai, pertvegf, z01d, zt1d, bexp1d, xlai1d, vegf1d,  &
                          cplflx, flag_cice, islmsk_cice,slimskin_cpl, dusfcin_cpl, dvsfcin_cpl,           &
                          dtsfcin_cpl, dqsfcin_cpl, ulwsfcin_cpl, ulwsfc_cice, dusfc_cice, dvsfc_cice,     &
                          dtsfc_cice, dqsfc_cice, tisfc, tsfco, fice, hice, dry, icy, wet,                 &
                          wind, u1, v1, cnvwind, errmsg, errflg)

        use surface_perturbation,  only: cdfnor

        implicit none

        ! Interface variables
        integer, intent(in) :: im, levs, isot, ivegsrc
        integer, dimension(im), intent(in) :: islmsk
        integer, dimension(im), intent(inout) :: soiltyp, vegtype, slopetyp
        logical, dimension(im), intent(in) :: dry, icy, wet

        real(kind=kind_phys), intent(in) :: con_g
        real(kind=kind_phys), dimension(im), intent(in) :: vfrac, stype, vtype, slope, prsik_1, prslk_1

        real(kind=kind_phys), dimension(im), intent(inout) :: tsfc
        real(kind=kind_phys), dimension(im,levs), intent(in) :: phil

        real(kind=kind_phys), dimension(im), intent(inout) :: sigmaf, work3, tsurf, zlvl

        ! Stochastic physics / surface perturbations
        logical, intent(in) :: do_sppt
        real(kind=kind_phys), dimension(im,levs),     intent(out) :: dtdtr
        real(kind=kind_phys), dimension(im),          intent(out) :: drain_cpl
        real(kind=kind_phys), dimension(im),          intent(out) :: dsnow_cpl
        real(kind=kind_phys), dimension(im),          intent(in)  :: rain_cpl
        real(kind=kind_phys), dimension(im),          intent(in)  :: snow_cpl
        logical, intent(in) :: do_sfcperts
        integer, intent(in) :: nsfcpert
        real(kind=kind_phys), dimension(im,nsfcpert), intent(in)  :: sfc_wts
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertz0
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertzt
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertshc
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertlai
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertvegf
        real(kind=kind_phys), dimension(im),          intent(out) :: z01d
        real(kind=kind_phys), dimension(im),          intent(out) :: zt1d
        real(kind=kind_phys), dimension(im),          intent(out) :: bexp1d
        real(kind=kind_phys), dimension(im),          intent(out) :: xlai1d
        real(kind=kind_phys), dimension(im),          intent(out) :: vegf1d

        logical, intent(in) :: cplflx
        real(kind=kind_phys), dimension(im), intent(in) :: slimskin_cpl
        logical, dimension(im), intent(inout) :: flag_cice
              integer, dimension(im), intent(out) :: islmsk_cice
        real(kind=kind_phys), dimension(im), intent(in) ::ulwsfcin_cpl, &
             dusfcin_cpl, dvsfcin_cpl, dtsfcin_cpl, dqsfcin_cpl, &
             tisfc, tsfco, fice, hice
        real(kind=kind_phys), dimension(im), intent(out) ::ulwsfc_cice, &
             dusfc_cice, dvsfc_cice, dtsfc_cice, dqsfc_cice

        real(kind=kind_phys), dimension(im), intent(out) :: wind
        real(kind=kind_phys), dimension(im), intent(in ) :: u1, v1
        ! surface wind enhancement due to convection
        real(kind=kind_phys), dimension(im), intent(in ) :: cnvwind

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer              :: i
        real(kind=kind_phys) :: onebg
        real(kind=kind_phys) :: cdfz

        ! Set constants
        onebg  = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Set initial quantities for stochastic physics deltas
        if (do_sppt) then
          dtdtr     = 0.0
          do i=1,im
            drain_cpl(i) = rain_cpl (i)
            dsnow_cpl(i) = snow_cpl (i)
          enddo
        endif

        ! Scale random patterns for surface perturbations with perturbation size
        ! Turn vegetation fraction pattern into percentile pattern
        if (do_sfcperts) then
          if (pertz0(1) > 0.) then
            z01d(:) = pertz0(1) * sfc_wts(:,1)
  !          if (me == 0) print*,'sfc_wts(:,1) min and max',minval(sfc_wts(:,1)),maxval(sfc_wts(:,1))
  !          if (me == 0) print*,'z01d min and max ',minval(z01d),maxval(z01d)
          endif
          if (pertzt(1) > 0.) then
            zt1d(:) = pertzt(1) * sfc_wts(:,2)
          endif
          if (pertshc(1) > 0.) then
            bexp1d(:) = pertshc(1) * sfc_wts(:,3)
          endif
          if (pertlai(1) > 0.) then
            xlai1d(:) = pertlai(1) * sfc_wts(:,4)
          endif
  ! --- do the albedo percentile calculation in GFS_radiation_driver instead --- !
  !        if (pertalb(1) > 0.) then
  !          do i=1,im
  !            call cdfnor(sfc_wts(i,5),cdfz)
  !            alb1d(i) = cdfz
  !          enddo
  !        endif
          if (pertvegf(1) > 0.) then
            do i=1,im
              call cdfnor(sfc_wts(i,6),cdfz)
              vegf1d(i) = cdfz
            enddo
          endif
        endif

        ! End of stochastic physics / surface perturbation

        do i=1,im
          sigmaf(i) = max(vfrac(i),0.01 )
          if (islmsk(i) == 2) then
            if (isot == 1) then
              soiltyp(i)  = 16
            else
              soiltyp(i)  = 9
            endif
            if (ivegsrc == 1) then
              vegtype(i)  = 15
            elseif(ivegsrc == 2) then
              vegtype(i)  = 13
            endif
            slopetyp(i) = 9
          else
            soiltyp(i)  = int( stype(i)+0.5 )
            vegtype(i)  = int( vtype(i)+0.5 )
            slopetyp(i) = int( slope(i)+0.5 )    !! clu: slope -> slopetyp
            if (soiltyp(i)  < 1) soiltyp(i)  = 14
            if (vegtype(i)  < 1) vegtype(i)  = 17
            if (slopetyp(i) < 1) slopetyp(i) = 1
          endif

          work3(i)   = prsik_1(i) / prslk_1(i)
        end do

        do i=1,im
          !tsurf(i) = tsfc(i)
          zlvl(i)  = phil(i,1) * onebg
          wind(i)  = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))   &
                         + max(zero, min(cnvwind(i), 30.0)), one)
          !wind(i)  = max(sqrt(Statein%ugrs(i,1)*Statein%ugrs(i,1) + &
          !                         Statein%vgrs(i,1)*Statein%vgrs(i,1))  &
          !              + max(zero, min(Tbd%phy_f2d(i,Model%num_p2d), 30.0)), one)
        end do


      if (cplflx) then
        do i=1,im
          islmsk_cice(i) = nint(slimskin_cpl(i))
          flag_cice(i)   = (islmsk_cice(i) == 4)

          if (flag_cice(i)) then
!           ulwsfc_cice(i) = ulwsfcin_cpl(i)
            dusfc_cice(i)  = dusfcin_cpl(i)
            dvsfc_cice(i)  = dvsfcin_cpl(i)
            dtsfc_cice(i)  = dtsfcin_cpl(i)
            dqsfc_cice(i)  = dqsfcin_cpl(i)
          endif
        enddo
      endif

      end subroutine GFS_surface_generic_pre_run

      end module GFS_surface_generic_pre


      module GFS_surface_generic_post

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_post_init, GFS_surface_generic_post_finalize, GFS_surface_generic_post_run

      real(kind=kind_phys), parameter :: one = 1.0d0
      real(kind=kind_phys), parameter :: zero = 0.0d0

      contains

      subroutine GFS_surface_generic_post_init ()
      end subroutine GFS_surface_generic_post_init

      subroutine GFS_surface_generic_post_finalize()
      end subroutine GFS_surface_generic_post_finalize

!> \section arg_table_GFS_surface_generic_post_run Argument Table
!! \htmlinclude GFS_surface_generic_post_run.html
!!
      subroutine GFS_surface_generic_post_run (im, cplflx, cplwav, lssav, icy, wet, dtf, ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1,&
        adjsfcdlw, adjsfcdsw, adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_ocn, adjnirbmu, adjnirdfu,           &
        adjvisbmu, adjvisdfu,t2m, q2m, u10m, v10m, tsfc, tsfc_ocn, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf,              &
        epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl,       &
        dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl,        &
        v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl,  &
        nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, evcwa, transa, sbsnoa, snowca, snohfa, ep,                                &
        runoff, srunoff, runof, drain, errmsg, errflg)

        implicit none

        integer,                                intent(in) :: im
        logical,                                intent(in) :: cplflx, cplwav, lssav
        logical, dimension(im),                 intent(in) :: icy, wet
        real(kind=kind_phys),                   intent(in) :: dtf

        real(kind=kind_phys), dimension(im),  intent(in)  :: ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1, adjsfcdlw, adjsfcdsw, &
          adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_ocn, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,    &
          t2m, q2m, u10m, v10m, tsfc, tsfc_ocn, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf

        real(kind=kind_phys), dimension(im),  intent(inout) :: epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, &
          dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl, dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, &
          nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl, v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, &
          nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl, nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, &
          evcwa, transa, sbsnoa, snowca, snohfa, ep

        real(kind=kind_phys), dimension(im), intent(inout) :: runoff, srunoff
        real(kind=kind_phys), dimension(im), intent(in)    :: drain, runof

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        real(kind=kind_phys), parameter :: albdf   = 0.06d0

        integer :: i
        real(kind=kind_phys) :: xcosz_loc, ocalnirdf_cpl, ocalnirbm_cpl, ocalvisdf_cpl, ocalvisbm_cpl

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        do i=1,im
          epi(i)     = ep1d(i)
          gfluxi(i)  = gflx(i)
          t1(i)      = tgrs_1(i)
          q1(i)      = qgrs_1(i)
          u1(i)      = ugrs_1(i)
          v1(i)      = vgrs_1(i)
        enddo

        if (cplflx .or. cplwav) then
          do i=1,im
            u10mi_cpl   (i) = u10m(i)
            v10mi_cpl   (i) = v10m(i)
          enddo
        endif

        if (cplflx) then
          do i=1,im
            dlwsfci_cpl (i) = adjsfcdlw(i)
            dswsfci_cpl (i) = adjsfcdsw(i)
            dlwsfc_cpl  (i) = dlwsfc_cpl(i) + adjsfcdlw(i)*dtf
            dswsfc_cpl  (i) = dswsfc_cpl(i) + adjsfcdsw(i)*dtf
            dnirbmi_cpl (i) = adjnirbmd(i)
            dnirdfi_cpl (i) = adjnirdfd(i)
            dvisbmi_cpl (i) = adjvisbmd(i)
            dvisdfi_cpl (i) = adjvisdfd(i)
            dnirbm_cpl  (i) = dnirbm_cpl(i) + adjnirbmd(i)*dtf
            dnirdf_cpl  (i) = dnirdf_cpl(i) + adjnirdfd(i)*dtf
            dvisbm_cpl  (i) = dvisbm_cpl(i) + adjvisbmd(i)*dtf
            dvisdf_cpl  (i) = dvisdf_cpl(i) + adjvisdfd(i)*dtf
            nlwsfci_cpl (i) = adjsfcdlw(i)  - adjsfculw(i)
            if (wet(i)) then
              nlwsfci_cpl(i) = adjsfcdlw(i) - adjsfculw_ocn(i)
            endif
            nlwsfc_cpl  (i) = nlwsfc_cpl(i) + nlwsfci_cpl(i)*dtf
            t2mi_cpl    (i) = t2m(i)
            q2mi_cpl    (i) = q2m(i)
!            tsfci_cpl   (i) = tsfc(i)
            tsfci_cpl   (i) = tsfc_ocn(i)
            psurfi_cpl  (i) = pgr(i)
          enddo

!  ---  estimate mean albedo for ocean point without ice cover and apply
!       them to net SW heat fluxes

          do i=1,im
!           if (Sfcprop%landfrac(i) < one) then ! Not 100% land
            if (wet(i)) then                    ! some open water 
!  ---  compute open water albedo
              xcosz_loc = max( 0.0, min( 1.0, xcosz(i) ))
              ocalnirdf_cpl = 0.06
              ocalnirbm_cpl = max(albdf, 0.026/(xcosz_loc**1.7+0.065)  &
       &                       + 0.15 * (xcosz_loc-0.1) * (xcosz_loc-0.5) &
       &                       * (xcosz_loc-1.0))
              ocalvisdf_cpl = 0.06
              ocalvisbm_cpl = ocalnirbm_cpl

              nnirbmi_cpl(i) = adjnirbmd(i) * (one-ocalnirbm_cpl)
              nnirdfi_cpl(i) = adjnirdfd(i) * (one-ocalnirdf_cpl)
              nvisbmi_cpl(i) = adjvisbmd(i) * (one-ocalvisbm_cpl)
              nvisdfi_cpl(i) = adjvisdfd(i) * (one-ocalvisdf_cpl)
            else
              nnirbmi_cpl(i) = adjnirbmd(i) - adjnirbmu(i)
              nnirdfi_cpl(i) = adjnirdfd(i) - adjnirdfu(i)
              nvisbmi_cpl(i) = adjvisbmd(i) - adjvisbmu(i)
              nvisdfi_cpl(i) = adjvisdfd(i) - adjvisdfu(i)
            endif
            nswsfci_cpl(i) = nnirbmi_cpl(i) + nnirdfi_cpl(i)   &
                            + nvisbmi_cpl(i) + nvisdfi_cpl(i)
            nswsfc_cpl(i)  = nswsfc_cpl(i)  + nswsfci_cpl(i)*dtf
            nnirbm_cpl(i)  = nnirbm_cpl(i)  + nnirbmi_cpl(i)*dtf
            nnirdf_cpl(i)  = nnirdf_cpl(i)  + nnirdfi_cpl(i)*dtf
            nvisbm_cpl(i)  = nvisbm_cpl(i)  + nvisbmi_cpl(i)*dtf
            nvisdf_cpl(i)  = nvisdf_cpl(i)  + nvisdfi_cpl(i)*dtf
          enddo
        endif

        if (lssav) then
          do i=1,im
            gflux(i)   = gflux(i)  + gflx(i)  * dtf
            evbsa(i)   = evbsa(i)  + evbs(i)  * dtf
            evcwa(i)   = evcwa(i)  + evcw(i)  * dtf
            transa(i)  = transa(i) + trans(i) * dtf
            sbsnoa(i)  = sbsnoa(i) + sbsno(i) * dtf
            snowca(i)  = snowca(i) + snowc(i) * dtf
            snohfa(i)  = snohfa(i) + snohf(i) * dtf
            ep(i)      = ep(i)     + ep1d(i)  * dtf
          enddo
        endif

!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters
        if (lssav) then
          do i=1,im
            runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * dtf
            srunoff(i) = srunoff(i) + runof(i) * dtf
          enddo
        endif

      end subroutine GFS_surface_generic_post_run

      end module GFS_surface_generic_post

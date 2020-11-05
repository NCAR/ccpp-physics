!> \file GFS_surface_generic.F90
!!  Contains code related to all GFS surface schemes.

      module GFS_surface_generic_pre

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_pre_init, GFS_surface_generic_pre_finalize, GFS_surface_generic_pre_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

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
                          sigmaf, soiltyp, vegtype, slopetyp, work3, tsurf, zlvl, do_sppt, ca_global,dtdtr,&
                          drain_cpl, dsnow_cpl, rain_cpl, snow_cpl, lndp_type, n_var_lndp, sfc_wts,        &
                          lndp_var_list, lndp_prt_list,                                                    &
                          z01d, zt1d, bexp1d, xlai1d, vegf1d, lndp_vgf,                                    &
                          cplflx, flag_cice, islmsk_cice, slimskin_cpl, tisfc, tsfco, fice, hice,          &
                          wind, u1, v1, cnvwind, smcwlt2, smcref2, errmsg, errflg)

        use surface_perturbation,  only: cdfnor

        implicit none

        ! Interface variables
        integer, intent(in) :: im, levs, isot, ivegsrc
        integer, dimension(im), intent(in) :: islmsk
        integer, dimension(im), intent(inout) :: soiltyp, vegtype, slopetyp

        real(kind=kind_phys), intent(in) :: con_g
        real(kind=kind_phys), dimension(im), intent(in) :: vfrac, stype, vtype, slope, prsik_1, prslk_1

        real(kind=kind_phys), dimension(im), intent(inout) :: tsfc
        real(kind=kind_phys), dimension(im,levs), intent(in) :: phil

        real(kind=kind_phys), dimension(im), intent(inout) :: sigmaf, work3, tsurf, zlvl

        ! Stochastic physics / surface perturbations
        logical, intent(in) :: do_sppt, ca_global
        real(kind=kind_phys), dimension(im,levs),     intent(out) :: dtdtr
        real(kind=kind_phys), dimension(im),          intent(out) :: drain_cpl
        real(kind=kind_phys), dimension(im),          intent(out) :: dsnow_cpl
        real(kind=kind_phys), dimension(im),          intent(in)  :: rain_cpl
        real(kind=kind_phys), dimension(im),          intent(in)  :: snow_cpl
        integer, intent(in) :: lndp_type
        integer, intent(in) :: n_var_lndp
        character(len=3), dimension(n_var_lndp),      intent(in)  :: lndp_var_list
        real(kind=kind_phys), dimension(n_var_lndp),  intent(in)  :: lndp_prt_list
        real(kind=kind_phys), dimension(im,n_var_lndp), intent(in)  :: sfc_wts
        real(kind=kind_phys), dimension(im),          intent(out) :: z01d
        real(kind=kind_phys), dimension(im),          intent(out) :: zt1d
        real(kind=kind_phys), dimension(im),          intent(out) :: bexp1d
        real(kind=kind_phys), dimension(im),          intent(out) :: xlai1d
        real(kind=kind_phys), dimension(im),          intent(out) :: vegf1d
        real(kind=kind_phys),                         intent(out) :: lndp_vgf

        logical, intent(in) :: cplflx
        real(kind=kind_phys), dimension(im), intent(in) :: slimskin_cpl
        logical, dimension(im), intent(inout) :: flag_cice
        integer, dimension(im), intent(out) :: islmsk_cice
        real(kind=kind_phys), dimension(im), intent(in) :: &
             tisfc, tsfco, fice, hice

        real(kind=kind_phys), dimension(im), intent(out) :: wind
        real(kind=kind_phys), dimension(im), intent(in ) :: u1, v1
        ! surface wind enhancement due to convection
        real(kind=kind_phys), dimension(im), intent(inout ) :: cnvwind
        !
        real(kind=kind_phys), dimension(im), intent(out) :: smcwlt2, smcref2

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer              :: i, k
        real(kind=kind_phys) :: onebg
        real(kind=kind_phys) :: cdfz

        ! Set constants
        onebg  = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Set initial quantities for stochastic physics deltas
        if (do_sppt .or. ca_global) then
          dtdtr     = 0.0
        endif

        ! Scale random patterns for surface perturbations with perturbation size
        ! Turn vegetation fraction pattern into percentile pattern
        lndp_vgf=-999.

        if (lndp_type==1) then
        do k =1,n_var_lndp
           select case(lndp_var_list(k))
           case ('rz0')
                z01d(:) = lndp_prt_list(k)* sfc_wts(:,k)
           case ('rzt')
                zt1d(:) = lndp_prt_list(k)* sfc_wts(:,k)
           case ('shc')
                bexp1d(:) = lndp_prt_list(k) * sfc_wts(:,k)
           case ('lai')
                xlai1d(:) = lndp_prt_list(k)* sfc_wts(:,k)
           case ('vgf')
        ! note that the pertrubed vegfrac is being used in sfc_drv, but not sfc_diff
              do i=1,im
                call cdfnor(sfc_wts(i,k),cdfz)
                vegf1d(i) = cdfz
              enddo
              lndp_vgf = lndp_prt_list(k)
           end select
        enddo
        endif

        ! End of stochastic physics / surface perturbation

        do i=1,im
          sigmaf(i) = max(vfrac(i), 0.01_kind_phys)
          islmsk_cice(i) = islmsk(i)
          if (islmsk(i) == 2) then
            if (isot == 1) then
              soiltyp(i)  = 16
            else
              soiltyp(i)  = 9
            endif
            if (ivegsrc == 0 .or. ivegsrc == 4) then
              vegtype(i)  = 24
            elseif (ivegsrc == 1) then
              vegtype(i)  = 15
            elseif (ivegsrc == 2) then
              vegtype(i)  = 13
            elseif (ivegsrc == 3 .or. ivegsrc == 5) then
              vegtype(i)  = 15
            endif
            slopetyp(i) = 9
          else
            soiltyp(i)  = int( stype(i)+0.5_kind_phys )
            vegtype(i)  = int( vtype(i)+0.5_kind_phys )
            slopetyp(i) = int( slope(i)+0.5_kind_phys )    !! clu: slope -> slopetyp
            if (soiltyp(i)  < 1) soiltyp(i)  = 14
            if (vegtype(i)  < 1) vegtype(i)  = 17
            if (slopetyp(i) < 1) slopetyp(i) = 1
          endif

          work3(i)   = prsik_1(i) / prslk_1(i)

          !tsurf(i) = tsfc(i)
          zlvl(i)    = phil(i,1) * onebg
          smcwlt2(i) = zero
          smcref2(i) = zero

          wind(i)  = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))   &
                         + max(zero, min(cnvwind(i), 30.0_kind_phys)), one)
          !wind(i)  = max(sqrt(Statein%ugrs(i,1)*Statein%ugrs(i,1) + &
          !                         Statein%vgrs(i,1)*Statein%vgrs(i,1))  &
          !              + max(zero, min(Tbd%phy_f2d(i,Model%num_p2d), 30.0)), one)
          cnvwind(i) = zero

        enddo

      if (cplflx) then
        do i=1,im
          islmsk_cice(i) = nint(slimskin_cpl(i))
          flag_cice(i)   = (islmsk_cice(i) == 4)
        enddo
      endif

      end subroutine GFS_surface_generic_pre_run

      end module GFS_surface_generic_pre


      module GFS_surface_generic_post

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_post_init, GFS_surface_generic_post_finalize, GFS_surface_generic_post_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      contains

      subroutine GFS_surface_generic_post_init ()
      end subroutine GFS_surface_generic_post_init

      subroutine GFS_surface_generic_post_finalize()
      end subroutine GFS_surface_generic_post_finalize

!> \section arg_table_GFS_surface_generic_post_run Argument Table
!! \htmlinclude GFS_surface_generic_post_run.html
!!
      subroutine GFS_surface_generic_post_run (im, cplflx, cplwav, lssav, icy, wet, dtf, ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1,&
        adjsfcdlw, adjsfcdsw, adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_wat, adjnirbmu, adjnirdfu,           &
        adjvisbmu, adjvisdfu,t2m, q2m, u10m, v10m, tsfc, tsfc_wat, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf,              &
        epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl,       &
        dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl,        &
        v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl,  &
        nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, evcwa, transa, sbsnoa, snowca, snohfa, ep,                                &
        runoff, srunoff, runof, drain, lheatstrg, z0fac, e0fac, zorl, hflx, evap, hflxq, evapq, hffac, hefac, errmsg, errflg)

        implicit none

        integer,                                intent(in) :: im
        logical,                                intent(in) :: cplflx, cplwav, lssav
        logical, dimension(im),                 intent(in) :: icy, wet
        real(kind=kind_phys),                   intent(in) :: dtf

        real(kind=kind_phys), dimension(im),  intent(in)  :: ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1, adjsfcdlw, adjsfcdsw, &
          adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_wat, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,    &
          t2m, q2m, u10m, v10m, tsfc, tsfc_wat, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf

        real(kind=kind_phys), dimension(im),  intent(inout) :: epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, &
          dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl, dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, &
          nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl, v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, &
          nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl, nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, &
          evcwa, transa, sbsnoa, snowca, snohfa, ep

        real(kind=kind_phys), dimension(im), intent(inout) :: runoff, srunoff
        real(kind=kind_phys), dimension(im), intent(in)    :: drain, runof

        ! For canopy heat storage
        logical, intent(in) :: lheatstrg
        real(kind=kind_phys), intent(in) :: z0fac, e0fac
        real(kind=kind_phys), dimension(im), intent(in)  :: zorl
        real(kind=kind_phys), dimension(im), intent(in)  :: hflx,  evap
        real(kind=kind_phys), dimension(im), intent(out) :: hflxq, evapq
        real(kind=kind_phys), dimension(im), intent(out) :: hffac, hefac

        ! CCPP error handling variables
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        real(kind=kind_phys), parameter :: albdf = 0.06_kind_phys

        ! Parameters for canopy heat storage parametrization
        real(kind=kind_phys), parameter :: z0min=0.2, z0max=1.0
        real(kind=kind_phys), parameter :: u10min=2.5, u10max=7.5

        integer :: i
        real(kind=kind_phys) :: xcosz_loc, ocalnirdf_cpl, ocalnirbm_cpl, ocalvisdf_cpl, ocalvisbm_cpl
        real(kind=kind_phys) :: tem, tem1, tem2

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        do i=1,im
          epi(i)    = ep1d(i)
          gfluxi(i) = gflx(i)
          t1(i)     = tgrs_1(i)
          q1(i)     = qgrs_1(i)
          u1(i)     = ugrs_1(i)
          v1(i)     = vgrs_1(i)
        enddo

        if (cplflx .or. cplwav) then
          do i=1,im
            u10mi_cpl(i) = u10m(i)
            v10mi_cpl(i) = v10m(i)
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
              nlwsfci_cpl(i) = adjsfcdlw(i) - adjsfculw_wat(i)
            endif
            nlwsfc_cpl  (i) = nlwsfc_cpl(i) + nlwsfci_cpl(i)*dtf
            t2mi_cpl    (i) = t2m(i)
            q2mi_cpl    (i) = q2m(i)
            tsfci_cpl   (i) = tsfc(i)
!           tsfci_cpl   (i) = tsfc_wat(i)
            psurfi_cpl  (i) = pgr(i)
          enddo

!  ---  estimate mean albedo for ocean point without ice cover and apply
!       them to net SW heat fluxes

          do i=1,im
!           if (Sfcprop%landfrac(i) < one) then ! Not 100% land
            if (wet(i)) then                    ! some open water 
!  ---  compute open water albedo
              xcosz_loc = max( zero, min( one, xcosz(i) ))
              ocalnirdf_cpl = 0.06_kind_phys
              ocalnirbm_cpl = max(albdf, 0.026_kind_phys/(xcosz_loc**1.7_kind_phys+0.065_kind_phys)     &
       &                       + 0.15_kind_phys * (xcosz_loc-0.1_kind_phys) * (xcosz_loc-0.5_kind_phys) &
       &                       * (xcosz_loc-one))
              ocalvisdf_cpl = 0.06_kind_phys
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

!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters
            runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * dtf
            srunoff(i) = srunoff(i) + runof(i) * dtf
          enddo
        endif

!  --- ...  Boundary Layer and Free atmospheic turbulence parameterization
!
!  in order to achieve heat storage within canopy layer, in the canopy heat
!    storage parameterization the kinematic sensible and latent heat fluxes
!    (hflx & evap) as surface boundary forcings to the pbl scheme are
!    reduced as a function of surface roughness
!
        do i=1,im
          hflxq(i) = hflx(i)
          evapq(i) = evap(i)
          hffac(i) = 1.0
          hefac(i) = 1.0
        enddo
        if (lheatstrg) then
          do i=1,im
            tem = 0.01 * zorl(i)     ! change unit from cm to m
            tem1 = (tem - z0min) / (z0max - z0min)
            hffac(i) = z0fac * min(max(tem1, 0.0), 1.0)
            tem = sqrt(u10m(i)**2+v10m(i)**2)
            tem1 = (tem - u10min) / (u10max - u10min)
            tem2 = 1.0 - min(max(tem1, 0.0), 1.0)
            hffac(i) = tem2 * hffac(i)
            hefac(i) = 1. + e0fac * hffac(i)
            hffac(i) = 1. + hffac(i)
            hflxq(i) = hflx(i) / hffac(i)
            evapq(i) = evap(i) / hefac(i)
          enddo
        endif

      end subroutine GFS_surface_generic_post_run

      end module GFS_surface_generic_post

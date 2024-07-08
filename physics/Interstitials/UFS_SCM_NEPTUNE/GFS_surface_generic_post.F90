!> \file GFS_surface_generic_post.F90
!!  Contains code related to all GFS surface schemes to be run afterward.

      module GFS_surface_generic_post

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_post_init, GFS_surface_generic_post_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      contains

!>\defgroup gfs_sfc_gen_post_mode GFS surface_generic_post Module
!! This module contains code related to all GFS surface schemes to be run afterward.
!> @{
!> \section arg_table_GFS_surface_generic_post_init Argument Table
!! \htmlinclude GFS_surface_generic_post_init.html
!!
      subroutine GFS_surface_generic_post_init (vtype, stype,scolor, slope, vtype_save, stype_save,scolor_save, slope_save, errmsg, errflg)

        integer, dimension(:), intent(in)  :: vtype_save, stype_save,scolor_save, slope_save
        integer, dimension(:), intent(out) :: vtype, stype, scolor,slope

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Restore vegetation, soil and slope type
        vtype(:) = vtype_save(:)
        stype(:) = stype_save(:)
        scolor(:) = scolor_save(:)
        slope(:) = slope_save(:)

      end subroutine GFS_surface_generic_post_init

!> \section arg_table_GFS_surface_generic_post_run Argument Table
!! \htmlinclude GFS_surface_generic_post_run.html
!!
      subroutine GFS_surface_generic_post_run (im, cplflx, cplaqm, cplchm, cplwav, cpllnd, lssav, dry, icy, wet,                    &
        lsm, lsm_noahmp, dtf, ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1,                                                           &
        adjsfcdlw, adjsfcdsw, adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_wat, adjnirbmu, adjnirdfu,           &
        adjvisbmu, adjvisdfu, t2m, q2m, u10m, v10m, tsfc, tsfc_wat, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf, pah, pahi,  &
        epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl,       &
        dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl,        &
        v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl,  &
        nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, evcwa, transa, sbsnoa, snowca, snohfa, paha, ep, ecan, etran, edir, waxy, &
        runoff, srunoff, runof, drain, tecan, tetran, tedir, twa, lheatstrg, h0facu, h0facs, zvfun, hflx, evap, hflxq, hffac,       &
        isot, ivegsrc, islmsk, vtype, stype,scolor, slope, vtype_save, stype_save,scolor_save, slope_save, errmsg, errflg)

        implicit none

        integer,                                intent(in) :: im
        logical,                                intent(in) :: cplflx, cplaqm, cplchm, cplwav, cpllnd, lssav
        logical, dimension(:),                  intent(in) :: dry, icy, wet
        integer,                                intent(in) :: lsm, lsm_noahmp
        real(kind=kind_phys),                   intent(in) :: dtf

        real(kind=kind_phys), dimension(:),  intent(in)  :: ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1, adjsfcdlw, adjsfcdsw,  &
          adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_wat, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,    &
          t2m, q2m, u10m, v10m, tsfc, tsfc_wat, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf, pah, ecan, etran, edir
        real(kind=kind_phys), dimension(:),  intent(in), optional  ::   &
          waxy

        real(kind=kind_phys), dimension(:),  intent(inout) :: epi, gfluxi, t1, q1, u1, v1,gflux, evbsa,       &
          evcwa, transa, sbsnoa, snowca, snohfa, ep, tecan, tetran, tedir
        real(kind=kind_phys), dimension(:),  intent(inout), optional :: pahi, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, &
          dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl, dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl,        &
          nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl, v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl,    &
          nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl, nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, paha, twa

        real(kind=kind_phys), dimension(:), intent(inout) :: runoff, srunoff
        real(kind=kind_phys), dimension(:), intent(in)    :: drain, runof

        ! For canopy heat storage
        logical, intent(in) :: lheatstrg
        real(kind=kind_phys), intent(in) :: h0facu, h0facs
        real(kind=kind_phys), dimension(:), intent(in)  :: zvfun
        real(kind=kind_phys), dimension(:), intent(in)  :: hflx,  evap
        real(kind=kind_phys), dimension(:), intent(out) :: hflxq
        real(kind=kind_phys), dimension(:), intent(out) :: hffac

        integer, intent(in) :: isot, ivegsrc, islmsk(:), vtype_save(:), stype_save(:),scolor_save(:), slope_save(:)
        integer, intent(out) :: vtype(:), stype(:),scolor(:), slope(:)

        ! CCPP error handling variables
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        real(kind=kind_phys), parameter :: albdf = 0.06_kind_phys

        integer :: i
        real(kind=kind_phys) :: xcosz_loc, ocalnirdf_cpl, ocalnirbm_cpl, ocalvisdf_cpl, ocalvisbm_cpl

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        do i=1,im
          epi(i)    = ep1d(i)
          gfluxi(i) = gflx(i)
          if (lsm == lsm_noahmp) then
            pahi(i)   = pah(i)
          endif
          t1(i)     = tgrs_1(i)
          q1(i)     = qgrs_1(i)
          u1(i)     = ugrs_1(i)
          v1(i)     = vgrs_1(i)
        enddo

        if (cplflx .or. cplchm .or. cplwav) then
          do i=1,im
            u10mi_cpl(i) = u10m(i)
            v10mi_cpl(i) = v10m(i)
          enddo
        endif

        if (cplflx .or. cplchm .or. cpllnd) then
          do i=1,im
            tsfci_cpl(i) = tsfc(i)
          enddo
        endif

        if (cplflx .or. cpllnd) then
          do i=1,im
            dlwsfci_cpl (i) = adjsfcdlw(i)
            dswsfci_cpl (i) = adjsfcdsw(i)
            dlwsfc_cpl  (i) = dlwsfc_cpl(i) + adjsfcdlw(i)*dtf
            dswsfc_cpl  (i) = dswsfc_cpl(i) + adjsfcdsw(i)*dtf
            psurfi_cpl  (i) = pgr(i)
          enddo
        endif

        if (cplflx) then
          do i=1,im
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
          enddo
        endif

!  ---  estimate mean albedo for ocean point without ice cover and apply
!       them to net SW heat fluxes

        if (cplflx .or. cpllnd) then
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

        if (cplaqm .and. .not.cplflx) then
          do i=1,im
            t2mi_cpl    (i) = t2m(i)
            q2mi_cpl    (i) = q2m(i)
            psurfi_cpl  (i) = pgr(i)
            if (wet(i)) then                    ! some open water
!  ---  compute open water albedo
              xcosz_loc = max( zero, min( one, xcosz(i) ))
              ocalnirdf_cpl = 0.06_kind_phys
              ocalnirbm_cpl = max(albdf, 0.026_kind_phys/(xcosz_loc**1.7_kind_phys+0.065_kind_phys)     &
       &                       + 0.15_kind_phys * (xcosz_loc-0.1_kind_phys) * (xcosz_loc-0.5_kind_phys) &
       &                       * (xcosz_loc-one))
              ocalvisdf_cpl = 0.06_kind_phys
              ocalvisbm_cpl = ocalnirbm_cpl

              nswsfci_cpl(i) = adjnirbmd(i) * (one-ocalnirbm_cpl) + &
                               adjnirdfd(i) * (one-ocalnirdf_cpl) + &
                               adjvisbmd(i) * (one-ocalvisbm_cpl) + &
                               adjvisdfd(i) * (one-ocalvisdf_cpl)
            else
              nswsfci_cpl(i) = adjnirbmd(i) - adjnirbmu(i) + &
                               adjnirdfd(i) - adjnirdfu(i) + &
                               adjvisbmd(i) - adjvisbmu(i) + &
                               adjvisdfd(i) - adjvisdfu(i)
            endif
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
            tecan(i)   = tecan(i)   + ecan(i) * dtf
            tetran(i)  = tetran(i)  + etran(i) * dtf
            tedir(i)   = tedir(i)   + edir(i) * dtf
            if (lsm == lsm_noahmp) then
             paha(i)    = paha(i)    + pah(i)   * dtf
             twa(i)     = waxy(i)
            endif
          enddo
        endif

!
!  in order to achieve heat storage within canopy layer, in the canopy
!    heat torage parameterization the kinematic sensible heat flux
!    (hflx) as surface boundary forcing to the pbl scheme is
!    reduced in a factor of hffac given as a function of surface roughness &
!    green vegetation fraction (zvfun)
!
        do i=1,im
          hflxq(i) = hflx(i)
          hffac(i) = 1.0
        enddo
        if (lheatstrg) then
          do i=1,im
            if (dry(i)) then
              if(hflx(i) > 0.) then
                hffac(i) = h0facu * zvfun(i)
              else
                hffac(i) = h0facs * zvfun(i)
              endif
              hffac(i) = 1. + hffac(i)
              hflxq(i) = hflx(i) / hffac(i)
            endif
          enddo
        endif

        ! Restore vegetation, soil and slope type
        vtype(:) = vtype_save(:)
        stype(:) = stype_save(:)
        scolor(:) = scolor_save(:)
        slope(:) = slope_save(:)

      end subroutine GFS_surface_generic_post_run
!> @}
      end module GFS_surface_generic_post

!> \file GFS_surface_generic.F90
!!  Contains code related to all GFS surface schemes.

!>\defgroup mod_GFS_surface_generic_pre GFS Surface Generic Pre module
      module GFS_surface_generic_pre

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_pre_init, GFS_surface_generic_pre_finalize, GFS_surface_generic_pre_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      contains

!> \section arg_table_GFS_surface_generic_pre_init Argument Table
!! \htmlinclude GFS_surface_generic_pre_init.html
!!
      subroutine GFS_surface_generic_pre_init (nthreads, im, slmsk, isot, ivegsrc, stype, vtype, slope, &
                                               vtype_save, stype_save, slope_save, errmsg, errflg)

        implicit none

        ! Interface variables
        integer,                       intent(in)    :: nthreads, im, isot, ivegsrc
        real(kind_phys), dimension(:), intent(in)    :: slmsk
        integer,         dimension(:), intent(inout) :: vtype, stype, slope
        integer,         dimension(:), intent(out)   :: vtype_save, stype_save, slope_save

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer, dimension(1:im) :: islmsk
        integer :: i

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        islmsk = nint(slmsk)

        ! Save current values of vegetation, soil and slope type
        vtype_save(:) = vtype(:)
        stype_save(:) = stype(:)
        slope_save(:) = slope(:)

        call update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype, slope)

      end subroutine GFS_surface_generic_pre_init

      subroutine GFS_surface_generic_pre_finalize()
      end subroutine GFS_surface_generic_pre_finalize

!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! \htmlinclude GFS_surface_generic_pre_run.html
!!
      subroutine GFS_surface_generic_pre_run (nthreads, im, levs, vfrac, islmsk, isot, ivegsrc, stype, vtype, slope, &
                          prsik_1, prslk_1, tsfc, phil, con_g, sigmaf, work3, zlvl,                        &
                          drain_cpl, dsnow_cpl, rain_cpl, snow_cpl, lndp_type, n_var_lndp, sfc_wts,        &
                          lndp_var_list, lndp_prt_list,                                                    &
                          z01d, zt1d, bexp1d, xlai1d, vegf1d, lndp_vgf,                                    &
                          cplflx, flag_cice, islmsk_cice, slimskin_cpl,                                    &
                          wind, u1, v1, cnvwind, smcwlt2, smcref2, vtype_save, stype_save, slope_save,     &
                          errmsg, errflg)

        use surface_perturbation,  only: cdfnor

        implicit none

        ! Interface variables
        integer, intent(in) :: nthreads, im, levs, isot, ivegsrc
        integer, dimension(:), intent(in) :: islmsk

        real(kind=kind_phys), intent(in) :: con_g
        real(kind=kind_phys), dimension(:), intent(in) :: vfrac, prsik_1, prslk_1
        integer, dimension(:), intent(inout) :: vtype, stype, slope
        integer, dimension(:), intent(out)   :: vtype_save(:), stype_save(:), slope_save(:)

        real(kind=kind_phys), dimension(:), intent(inout) :: tsfc
        real(kind=kind_phys), dimension(:,:), intent(in) :: phil

        real(kind=kind_phys), dimension(:), intent(inout) :: sigmaf, work3, zlvl

        ! Stochastic physics / surface perturbations
        real(kind=kind_phys), dimension(:),   intent(out) :: drain_cpl
        real(kind=kind_phys), dimension(:),   intent(out) :: dsnow_cpl
        real(kind=kind_phys), dimension(:),   intent(in)  :: rain_cpl
        real(kind=kind_phys), dimension(:),   intent(in)  :: snow_cpl
        integer,                              intent(in)  :: lndp_type, n_var_lndp
        character(len=3),     dimension(:),   intent(in)  :: lndp_var_list
        real(kind=kind_phys), dimension(:),   intent(in)  :: lndp_prt_list
        real(kind=kind_phys), dimension(:,:), intent(in)  :: sfc_wts
        real(kind=kind_phys), dimension(:),   intent(out) :: z01d
        real(kind=kind_phys), dimension(:),   intent(out) :: zt1d
        real(kind=kind_phys), dimension(:),   intent(out) :: bexp1d
        real(kind=kind_phys), dimension(:),   intent(out) :: xlai1d
        real(kind=kind_phys), dimension(:),   intent(out) :: vegf1d
        real(kind=kind_phys),                 intent(out) :: lndp_vgf

        logical,                              intent(in)    :: cplflx
        real(kind=kind_phys), dimension(:),   intent(in)    :: slimskin_cpl
        logical,              dimension(:),   intent(inout) :: flag_cice
        integer,              dimension(:),   intent(out)   :: islmsk_cice

        real(kind=kind_phys), dimension(:),   intent(out) :: wind
        real(kind=kind_phys), dimension(:),   intent(in ) :: u1, v1
        ! surface wind enhancement due to convection
        real(kind=kind_phys), dimension(:),   intent(inout ) :: cnvwind
        !
        real(kind=kind_phys), dimension(:),   intent(out)    :: smcwlt2, smcref2

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer              :: i, k
        real(kind=kind_phys) :: onebg, cdfz

        ! Set constants
        onebg  = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

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

        ! Save current values of vegetation, soil and slope type
        vtype_save(:) = vtype(:)
        stype_save(:) = stype(:)
        slope_save(:) = slope(:)

        call update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype, slope)

        do i=1,im
          sigmaf(i) = max(vfrac(i), 0.01_kind_phys)
          islmsk_cice(i) = islmsk(i)

          work3(i)   = prsik_1(i) / prslk_1(i)

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

      subroutine update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype, slope)

        implicit none

        integer, intent(in)    :: nthreads, im, isot, ivegsrc, islmsk(:)
        integer, intent(inout) :: vtype(:), stype(:), slope(:)
        integer :: i

!$OMP  parallel do num_threads(nthreads) default(none) private(i) &
!$OMP      shared(im, isot, ivegsrc, islmsk, vtype, stype, slope)
        do i=1,im
          if (islmsk(i) == 2) then
            if (isot == 1) then
              stype(i) = 16
            else
              stype(i) = 9
            endif
            if (ivegsrc == 0 .or. ivegsrc == 4) then
              vtype(i) = 24
            elseif (ivegsrc == 1) then
              vtype(i) = 15
            elseif (ivegsrc == 2) then
              vtype(i) = 13
            elseif (ivegsrc == 3 .or. ivegsrc == 5) then
              vtype(i) = 15
            endif
            slope(i)  = 9
          else
            if (vtype(i)  < 1) vtype(i)  = 17
            if (slope(i) < 1) slope(i) = 1
          endif
        enddo
!$OMP end parallel do

      end subroutine update_vegetation_soil_slope_type

      end module GFS_surface_generic_pre


      module GFS_surface_generic_post

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_post_init, GFS_surface_generic_post_finalize, GFS_surface_generic_post_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      contains

!> \section arg_table_GFS_surface_generic_post_init Argument Table
!! \htmlinclude GFS_surface_generic_post_init.html
!!
      subroutine GFS_surface_generic_post_init (vtype, stype, slope, vtype_save, stype_save, slope_save, errmsg, errflg)

        integer, dimension(:), intent(in)  :: vtype_save, stype_save, slope_save
        integer, dimension(:), intent(out) :: vtype, stype, slope

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Restore vegetation, soil and slope type
        vtype(:) = vtype_save(:)
        stype(:) = stype_save(:)
        slope(:) = slope_save(:)

      end subroutine GFS_surface_generic_post_init

      subroutine GFS_surface_generic_post_finalize()
      end subroutine GFS_surface_generic_post_finalize

!> \section arg_table_GFS_surface_generic_post_run Argument Table
!! \htmlinclude GFS_surface_generic_post_run.html
!!
      subroutine GFS_surface_generic_post_run (im, cplflx, cplaqm, cplchm, cplwav, lssav, dry, icy, wet,                            &
        lsm, lsm_noahmp, dtf, ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1,                                                           &
        adjsfcdlw, adjsfcdsw, adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_wat, adjnirbmu, adjnirdfu,           &
        adjvisbmu, adjvisdfu, t2m, q2m, u10m, v10m, tsfc, tsfc_wat, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf, pah, pahi,  &
        epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl,       &
        dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl,        &
        v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl,  &
        nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, evcwa, transa, sbsnoa, snowca, snohfa, paha, ep, ecan, etran, edir, waxy, &
        runoff, srunoff, runof, drain, tecan, tetran, tedir, twa, lheatstrg, h0facu, h0facs, zvfun, hflx, evap, hflxq, hffac,       &
        isot, ivegsrc, islmsk, vtype, stype, slope, vtype_save, stype_save, slope_save, errmsg, errflg)

        implicit none

        integer,                                intent(in) :: im
        logical,                                intent(in) :: cplflx, cplaqm, cplchm, cplwav, lssav
        logical, dimension(:),                  intent(in) :: dry, icy, wet
        integer,                                intent(in) :: lsm, lsm_noahmp
        real(kind=kind_phys),                   intent(in) :: dtf

        real(kind=kind_phys), dimension(:),  intent(in)  :: ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1, adjsfcdlw, adjsfcdsw,  &
          adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjsfculw_wat, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,    &
          t2m, q2m, u10m, v10m, tsfc, tsfc_wat, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf, pah, ecan, etran, edir,    &
          waxy

        real(kind=kind_phys), dimension(:),  intent(inout) :: epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, &
          dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl, dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl,        &
          nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl, v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl,    &
          nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl, nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa,       &
          evcwa, transa, sbsnoa, snowca, snohfa, ep, paha, tecan, tetran, tedir, twa, pahi

        real(kind=kind_phys), dimension(:), intent(inout) :: runoff, srunoff
        real(kind=kind_phys), dimension(:), intent(in)    :: drain, runof

        ! For canopy heat storage
        logical, intent(in) :: lheatstrg
        real(kind=kind_phys), intent(in) :: h0facu, h0facs
        real(kind=kind_phys), dimension(:), intent(in)  :: zvfun
        real(kind=kind_phys), dimension(:), intent(in)  :: hflx,  evap
        real(kind=kind_phys), dimension(:), intent(out) :: hflxq
        real(kind=kind_phys), dimension(:), intent(out) :: hffac

        integer, intent(in) :: isot, ivegsrc, islmsk(:), vtype_save(:), stype_save(:), slope_save(:)
        integer, intent(out) :: vtype(:), stype(:), slope(:)

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

        if (cplflx .or. cplchm) then
          do i=1,im
            tsfci_cpl(i) = tsfc(i)
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
        slope(:) = slope_save(:)

      end subroutine GFS_surface_generic_post_run

      end module GFS_surface_generic_post

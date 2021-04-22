!>  \file gfdl_sfc_layer.f
!! This file contains ...

!> This module contains the CCPP-compliant GFDL surface layer scheme.
      module gfdl_sfc_layer

      use machine , only : kind_phys

      implicit none

      public :: gfdl_sfc_layer_init, gfdl_sfc_layer_run, gfdl_sfc_layer_finalize

      private

      contains

!> \section arg_table_gfdl_sfc_layer_init Argument Table
!! \htmlinclude gfdl_sfc_layer_init.html
!!
      subroutine gfdl_sfc_layer_init (icoef_sf, cplwav, cplwav2atm, lcurr_sf,   &
        pert_cd, ntsflg, errmsg, errflg)

        implicit none

        integer, intent(in) :: icoef_sf, ntsflg
        logical, intent(in) :: cplwav, cplwav2atm, lcurr_sf, pert_cd

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

#if HWRF==1
        write(errmsg,'(*(a))') 'The GFDL surface layer scheme does not support '&
          //'use of the HWRF preprocessor flag in gfdl_sfc_layer.F90'
        errflg = 1
        return
#endif

        if (icoef_sf < 0 .or. icoef_sf > 8) then
          write(errmsg,'(*(a))') 'The value of icoef_sf is outside of the '     &
            //'supported range (0-8) in gfdl_sfc_layer.F90'
          errflg = 1
          return
        end if

        if (cplwav .or. cplwav2atm) then
          write(errmsg,'(*(a))') 'The GFDL surface layer scheme is not set up ' &
            //'to be coupled to waves in gfdl_sfc_layer.F90'
          errflg = 1
          return
        end if

        if (lcurr_sf) then
          write(errmsg,'(*(a))') 'The GFDL surface layer scheme is not set up ' &
            //'to be used with the lcurr_sf option in gfdl_sfc_layer.F90'
          errflg = 1
          return
        end if

        if (pert_cd) then
          write(errmsg,'(*(a))') 'The GFDL surface layer scheme is not set up ' &
            //'to be used with the pert_cd option in gfdl_sfc_layer.F90'
          errflg = 1
          return
        end if

        if (ntsflg > 0) then
          !GJF: In order to enable ntsflg > 0, the variable 'tstrc' passed into MFLUX2 should be set
          !  to the surface_skin_temperature_over_X_interstitial rather than the average of it and
          !  surface_skin_temperature_after_iteration_over_X
          write(errmsg,'(*(a))') 'Setting ntsflg > 0 is currently not supported'&
            //' in gfdl_sfc_layer.F90'
          errflg = 1
          return
        end if

        !GJF: Initialization notes: In WRF, the subroutine module_sf_myjsfc/myjsfcinit
        !     is called for initialization of the GFDL surface layer scheme from
        !     the module_physics_init subroutine. It contains the following
        !     initializations which should already have been done by other
        !     code in UFS-related host models:
        ! IF(.NOT.RESTART)THEN
        !   DO J=JTS,JTE
        !   DO I=ITS,ITF
        !     USTAR(I,J)=0.1
        !   ENDDO
        !   ENDDO
        ! ENDIF
        !also initialize surface roughness length

      end subroutine gfdl_sfc_layer_init

      subroutine gfdl_sfc_layer_finalize ()
      end subroutine gfdl_sfc_layer_finalize

!> \section arg_table_gfdl_sfc_layer_run Argument Table
!! \htmlinclude gfdl_sfc_layer_run.html
!!
      subroutine gfdl_sfc_layer_run (im, nsoil, km, xlat, xlon, flag_iter, lsm, &
        lsm_noah, lsm_noahmp, lsm_ruc, lsm_noah_wrfv4, icoef_sf, cplwav,        &
        cplwav2atm, lcurr_sf, pert_Cd, ntsflg, sfenth, z1, shdmax, ivegsrc,     &
        vegtype, sigmaf, dt, wet, dry, icy, isltyp, rd, grav, ep1, ep2, smois,  &
        psfc, prsl1, q1, t1, u1, v1, wspd, u10, v10, gsw, glw, tsurf_wat,       &
        tsurf_lnd, tsurf_ice, tskin_wat, tskin_lnd, tskin_ice, ustar_wat,       &
        ustar_lnd, ustar_ice, znt_wat, znt_lnd, znt_ice, cdm_wat, cdm_lnd,      &
        cdm_ice, stress_wat, stress_lnd, stress_ice, rib_wat, rib_lnd, rib_ice, &
        fm_wat, fm_lnd, fm_ice, fh_wat, fh_lnd, fh_ice, fh2_wat, fh2_lnd,       &
        fh2_ice, ch_wat, ch_lnd, ch_ice, fm10_wat, fm10_lnd, fm10_ice, qss_wat, &
        qss_lnd, qss_ice, errmsg, errflg)

        use funcphys, only: fpvs

        !####  GJF: temporarily grab parameters from LSM-specific modules -- should go through CCPP ####
        !           (fixing this involves replacing the functionality of set_soilveg and namelist_soilveg)
        use namelist_soilveg, only: maxsmc_noah => maxsmc, drysmc_noah => drysmc
        use namelist_soilveg_ruc, only: maxsmc_ruc => maxsmc, drysmc_ruc => drysmc
        use noahmp_tables, only: maxsmc_noahmp => smcmax_table, drysmc_noahmp => smcdry_table
        use module_sf_noahlsm, only: maxsmc_noah_wrfv4 => maxsmc, drysmc_noah_wrfv4 => drysmc
        !################################################################################################

        implicit none

        integer,                intent(in) :: im, nsoil, km, ivegsrc
        integer,                intent(in) :: lsm, lsm_noah, lsm_noahmp,        &
                                              lsm_ruc, lsm_noah_wrfv4, icoef_sf,&
                                              ntsflg
        logical,                intent(in) :: cplwav, cplwav2atm !GJF: this scheme has not been tested with these on
        logical,                intent(in) :: lcurr_sf           !GJF: this scheme has not been tested with this option turned on; the variables scurx and scury need to be input in order to use this
        logical,                intent(in) :: pert_Cd            !GJF: this scheme has not been tested with this option turned on; the variables ens_random_seed and ens_Cdamp need to be input in order to use this
        logical, dimension(im), intent(in) :: flag_iter, wet, dry, icy
        integer, dimension(im), intent(in) :: isltyp, vegtype
        real(kind=kind_phys),                      intent(in) :: dt, sfenth
        real(kind=kind_phys),                      intent(in) :: rd,grav,ep1,ep2
        real(kind=kind_phys), dimension(im,nsoil), intent(in) :: smois
        real(kind=kind_phys), dimension(im),       intent(in) :: psfc, prsl1,   &
            q1, t1, u1, v1, wspd, u10, v10, gsw, glw, z1, shdmax, sigmaf, xlat, &
            xlon, tsurf_wat, tsurf_lnd, tsurf_ice

        real(kind=kind_phys), intent(inout), dimension(im) :: tskin_wat,        &
            tskin_lnd, tskin_ice, ustar_wat, ustar_lnd, ustar_ice,              &
            znt_wat, znt_lnd, znt_ice, cdm_wat, cdm_lnd, cdm_ice,               &
            stress_wat, stress_lnd, stress_ice, rib_wat, rib_lnd, rib_ice,      &
            fm_wat, fm_lnd, fm_ice, fh_wat, fh_lnd, fh_ice, fh2_wat, fh2_lnd,   &
            fh2_ice, ch_wat, ch_lnd, ch_ice, fm10_wat, fm10_lnd, fm10_ice,      &
            qss_wat, qss_lnd, qss_ice

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        !local variables

        integer :: i, its, ite, ims, ime

        logical :: ch_bound_excursion

        !GJF: the vonKarman constant should come in through the CCPP and be defined by the host model
        real (kind=kind_phys), parameter :: karman = 0.4
        real (kind=kind_phys), parameter :: log01=log(0.01), log05=log(0.05),   &
            log07=log(0.07)

        !GJF: if the following variables will be used, they should be turned into intent(in) namelist options
        integer :: iwavecpl, ens_random_seed, issflx
        logical :: diag_wind10m, diag_qss
        real(kind=kind_phys) :: ens_Cdamp

        real(kind=kind_phys), dimension(im)   :: wetc, pspc, pkmax, tstrc, upc, &
            vpc, mznt, slwdc, wind10, qfx, qgh, zkmax, z1_cm, z0max, ztmax
        real(kind=kind_phys), dimension(im)   :: u10_lnd, u10_ocn, u10_ice,     &
            v10_lnd, v10_ocn, v10_ice

        !GJF: the following variables are identified as:
        !"SCURX"       "Surface Currents(X)"                    "m s-1"
        !"SCURY"       "Surface Currents(Y)"                    "m s-1
        !"CHARN"       "Charnock Coeff"                         " "
        !"MSANG"       "Wind/Stress Angle"                      "Radian"
        real(kind=kind_phys), dimension(im)   :: charn, msang, scurx, scury

        real(kind=kind_phys), dimension(im)   :: fxh, fxe, fxmx, fxmy, xxfh,    &
                                                 xxfh2, tzot
        real(kind=kind_phys), dimension(1:30) :: maxsmc, drysmc
        real(kind=kind_phys)                  :: smcmax, smcdry, zhalf, cd10,   &
            esat, fm_lnd_old, fh_lnd_old, tem1, tem2, czilc, cd_low_limit,      &
            cd_high_limit, ch_low_limit, ch_high_limit, fh2_fh_ratio

        !#### This block will become unnecessary when maxsmc and drysmc come through the CCPP ####
        if (lsm == lsm_noah) then
          maxsmc = maxsmc_noah
          drysmc = drysmc_noah
        else if (lsm == lsm_noahmp) then
          maxsmc = maxsmc_noahmp
          drysmc = drysmc_noahmp
        else if (lsm == lsm_ruc) then
          maxsmc = maxsmc_ruc
          drysmc = drysmc_ruc
        else if (lsm == lsm_noah_wrfv4) then
          maxsmc = maxsmc_noah_wrfv4
          drysmc = drysmc_noah_wrfv4
        else
          !GJF: These data were from the original GFDL surface layer scheme, but
          !     rather than being hard-coded here, they should be shared with the
          !     LSM. These data are kept for legacy purposes. Note that these only
          !     have nonzero values for 16 soil types vs 19 for other STAS datasets
          data maxsmc/0.339, 0.421, 0.434, 0.476, 0.476, 0.439,  &
                       0.404, 0.464, 0.465, 0.406, 0.468, 0.468,  &
                       0.439, 1.000, 0.200, 0.421, 0.000, 0.000,  &
                       0.000, 0.000, 0.000, 0.000, 0.000, 0.000,  &
                       0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
          data drysmc/0.010, 0.028, 0.047, 0.084, 0.084, 0.066,     &
                       0.067, 0.120, 0.103, 0.100, 0.126, 0.138,     &
                       0.066, 0.000, 0.006, 0.028, 0.000, 0.000,     &
                       0.000, 0.000, 0.000, 0.000, 0.000, 0.000,     &
                       0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
        end if
        !########################################################################

        !GJF: This code has not been tested with iwavecpl = 1; the variables 'charn' and 'msang' (and others?) need to be input in order to use this
        ! if (cplwav .or. cplwav2atm) then
        !   iwavecpl = 1
        ! else
        !   iwavecpl = 0
        ! end if
        iwavecpl = 0

        !GJF: temporary setting of variables that should be moved to namelist is they are used
        ens_random_seed = 0   !used for HWRF ensemble?
        ens_Cdamp = 0.0       !used for HWRF ensemble?

        issflx = 0              !GJF:  1 = calculate surface fluxes, 0 = don't
        diag_wind10m = .false.  !GJF: if one wants 10m wind speeds to come from this scheme, set this to True,
                                !  put [u,v]10_[lnd/ocn/ice] in the scheme argument list (and metadata), and modify
                                !  GFS_surface_compsites to receive the individual components and calculate an all-grid value
        diag_qss = .false.      !GJF: saturation specific humidities are calculated by LSM, sea surface, and sea ice schemes in
                                !  GFS-based suites

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        its = 1
        ims = 1
        ite = im
        ime = im

        do i=its, ite
          if (flag_iter(i)) then
            !GJF: Perform data preparation that is the same for all surface types

            pspc(i) = psfc(i)*10.     ! convert from Pa to cgs
            pkmax(i) = prsl1(i)*10.   ! convert from Pa to cgs

            upc(i) = u1(i)*100.       ! convert from m s-1 to cm s-1
            vpc(i) = v1(i)*100.       ! convert from m s-1 to cm s-1

            !Wang:  use previous u10 v10 to compute wind10, input to MFLUX2 to compute z0 (for first time step, u10 and v10 may be zero)
            wind10(i)=sqrt(u10(i)*u10(i)+v10(i)*v10(i)) !m s-1

            !Wang: calulate height of the first half level
            ! if (wind10(i) <= 1.0e-10 .or. wind10(i) > 150.0) then
            !   zhalf = -rd*t1(i)*alog(pkmax(i)/pspc(i))/grav !m
            ! endif

            !GJF: rather than calculate the height of the first half level, if it is precalculated
            !  in a different scheme, pass it in and use it; note that in FV3, calculating via the hypsometric equation
            !  occasionally produced values much shallower than those passed in
            !zkmax(i) = -rd*t1(i)*alog(pkmax(i)/pspc(i))/grav !m
            zkmax(i) = z1(i)
            z1_cm(i) = 100.0*z1(i)

            !GJF: these drag coefficient limits were suggested by Chunxi Zhang via his module_sf_sfclayrev.f90
            cd_low_limit  = 1.0e-5/zkmax(i)
            cd_high_limit = 0.1
            !GJF: use the lower of 0.1 from Chunxi Zhang or 0.05/wspd from WRF's module_sf_gfdl.F
            !    (this will always be the latter if wspd has a minimum of 1.0 m s-1 from above)
            ch_low_limit = cd_low_limit
            ch_high_limit = min(0.1,0.05/wspd(i))

            !slwdc... GFDL downward net flux in units of cal/(cm**2/min)
            !also divide by 10**4 to convert from /m**2 to /cm**2
            slwdc(i)=gsw(i)+glw(i)
            slwdc(i)=0.239*60.*slwdc(i)*1.e-4

            !GJF: these variables should be passed in if these options are used
            charn(i) = 0.0 !used with wave coupling (iwavecpl == 1)
            msang(i) = 0.0 !used with wave coupling (iwavecpl == 1)
            scurx(i) = 0.0 !used with ocean currents? (lcurr_sf == T)
            scury(i) = 0.0 !used with ocean currents? (lcurr_sf == T)

            if (diag_qss) then
              esat = fpvs(t1(i))
              qgh(i) = ep2*esat/(psfc(i)-esat)
            end if

            !GJF: these vars are not needed in a GFS-based suite
            !rho1(i)=prsl1(i)/(rd*t1(i)*(1.+ep1*q1(i)))
            !cpm(i)=cp*(1.+0.8*q1(i))

            !GJF: perform data preparation that depends on surface types and call the mflux2 subroutine for each surface type
            !  Note that this is different than the original WRF module_sf_gfdl.F where mflux2 is called once for all surface
            !  types, with negative roughness lengths denoting open ocean.
            if (dry(i)) then
              !GJF: from WRF's module_sf_gfdl.F
              smcdry=drysmc(isltyp(i))
              smcmax=maxsmc(isltyp(i))
              wetc(i)=(smois(i,1)-smcdry)/(smcmax-smcdry)
              wetc(i)=amin1(1.,amax1(wetc(i),0.))

              !GJF: the lower boundary temperature passed in to MFLUX2 either follows GFS:
              tstrc(i) = 0.5*(tskin_lnd(i) + tsurf_lnd(i)) !averaging tskin_lnd and tsurf_lnd as in GFS surface layer breaks ntsflg functionality
              !GJF: or WRF module_sf_gfdl.F:
              !tstrc(i) = tskin_lnd(i)

              !GJF: Roughness Length Limitation section
              !  The WRF version of module_sf_gfdl.F has no checks on the roughness lengths prior to entering MFLUX2.
              !  The following limits were placed on roughness lengths from the GFS surface layer scheme at the suggestion
              !  of Chunxi Zhang. Using the GFDL surface layer without such checks can lead to instability in the UFS.

              !znt_lnd is in cm, z0max/ztmax are in m at this point
              z0max(i) = max(1.0e-6, min(0.01 * znt_lnd(i), zkmax(i)))

              tem1  = 1.0 - shdmax(i)
              tem2  = tem1 * tem1
              tem1  = 1.0  - tem2

              if( ivegsrc == 1 ) then
                if (vegtype(i) == 10) then
                  z0max(i) = exp( tem2*log01 + tem1*log07 )
                elseif (vegtype(i) == 6) then
                  z0max(i) = exp( tem2*log01 + tem1*log05 )
                elseif (vegtype(i) == 7) then
  !               z0max(i) = exp( tem2*log01 + tem1*log01 )
                  z0max(i) = 0.01
                elseif (vegtype(i) == 16) then
  !               z0max(i) = exp( tem2*log01 + tem1*log01 )
                  z0max(i) = 0.01
                else
                  z0max(i) = exp( tem2*log01 + tem1*log(z0max(i)) )
                endif
              elseif (ivegsrc == 2 ) then
                if (vegtype(i) == 7) then
                  z0max(i) = exp( tem2*log01 + tem1*log07 )
                elseif (vegtype(i) == 8) then
                  z0max(i) = exp( tem2*log01 + tem1*log05 )
                elseif (vegtype(i) == 9) then
  !               z0max(i) = exp( tem2*log01 + tem1*log01 )
                  z0max(i) = 0.01
                elseif (vegtype(i) == 11) then
  !               z0max(i) = exp( tem2*log01 + tem1*log01 )
                  z0max(i) = 0.01
                else
                  z0max(i) = exp( tem2*log01 + tem1*log(z0max(i)) )
                endif
              endif

              z0max(i) = max(z0max(i), 1.0e-6)

  !           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
              czilc = 0.8

              tem1  = 1.0 - sigmaf(i)
              ztmax(i) = z0max(i)*exp( - tem1*tem1 &
       &                     * czilc*karman*sqrt(ustar_lnd(i)*(0.01/1.5e-05)))
              ztmax(i) = max(ztmax(i), 1.0e-6)

              !GJF: from WRF's module_sf_gfdl.F
              if (wind10(i) <= 1.0e-10 .or. wind10(i) > 150.0) then
                 wind10(i)=wspd(i)*alog(10.0/z0max(i))/alog(z1(i)/z0max(i)) !m s-1
              end if
              wind10(i)=wind10(i)*100.0   !convert from m/s to cm/s

              ztmax(i) = ztmax(i)*100.0   !convert from m to cm
              z0max(i) = z0max(i)*100.0   !convert from m to cm

              call mflux2 (fxh(i), fxe(i), fxmx(i), fxmy(i), cdm_lnd(i), rib_lnd(i), &
                xxfh(i), ztmax(i), z0max(i), tstrc(i),   &
                pspc(i), pkmax(i), wetc(i), slwdc(i), z1_cm(i), icoef_sf, iwavecpl, lcurr_sf, charn(i), msang(i), &
                scurx(i), scury(i), pert_Cd, ens_random_seed, ens_Cdamp, upc(i), vpc(i), t1(i), q1(i), &
                dt, wind10(i), xxfh2(i), ntsflg, sfenth, tzot(i), ep2, errmsg, &
                errflg)
                if (errflg /= 0) return

              !GJF: this is broken when tstrc is set to an average of two variables
              if (ntsflg==1) then
                tskin_lnd(i) = tstrc(i)      ! gopal's doing
              end if

              if (diag_wind10m) then
                u10_lnd(i) = u1(i)*(0.01*wind10(i)/wspd(i))
                v10_lnd(i) = v1(i)*(0.01*wind10(i)/wspd(i))
              end if

              !GJF: these variables are not needed in a GFS-based suite, but are found in WRF's module_sf_gfdl.F and kept in comments for legacy
              !gz1oz0(i) = alog(zkmax(i)/(0.01*znt_lnd(i)))
              !taux(i) = fxmx(i)/10.    ! gopal's doing for Ocean coupling
              !tauy(i) = fxmy(i)/10.    ! gopal's doing for Ocean coupling

              cdm_lnd(i) = max(cdm_lnd(i), cd_low_limit)
              cdm_lnd(i) = min(cdm_lnd(i), cd_high_limit)
              fm_lnd(i) = karman/sqrt(cdm_lnd(i))

              !1) try fh_lnd from MFLUX2
              fh_lnd(i) = karman*xxfh(i)

              !2) calc ch_lnd from fm_lnd and fh_lnd
              ch_lnd(i)  = karman*karman/(fm_lnd(i) * fh_lnd(i))

              !3) check if ch_lnd is out of bounds (if so, recalculate fh_lnd from bounded value)
              ch_bound_excursion = .false.
              if (ch_lnd(i) < ch_low_limit) then
                ch_bound_excursion = .true.
                ch_lnd(i) = ch_low_limit
              else if (ch_lnd(i) > ch_high_limit) then
                ch_bound_excursion = .true.
                ch_lnd(i) = ch_high_limit
              end if

              fh2_lnd(i) = karman*xxfh2(i)

              if (ch_bound_excursion) then
                fh2_fh_ratio = min(xxfh2(i)/xxfh(i), 1.0)
                fh_lnd(i) = karman*karman/(fm_lnd(i)*ch_lnd(i))
                fh2_lnd(i) = fh2_fh_ratio*fh_lnd(i)
              end if

              !GJF: Other CCPP schemes (PBL) ask for fm/fh instead of psim/psih
              !psim_lnd(i)=gz1oz0(i)-fm_lnd(i)
              !psih_lnd(i)=gz1oz0(i)-fh_lnd(i)

              !GJF: from WRF's module_sf_gfdl.F
              ustar_lnd(i) = 0.01*sqrt(cdm_lnd(i)*   &
                         (upc(i)*upc(i) + vpc(i)*vpc(i)))
              !GJF: from Chunxi Zhang's module_sf_sfclayrev.f90 (I'm not sure it's necessary.)
              ustar_lnd(i) = amax1(ustar_lnd(i),0.001)

              stress_lnd(i) = cdm_lnd(i)*wspd(i)*wspd(i)

              !GJF: from WRF's module_sf_gfdl.F
              ! convert cd, ch to values at 10m, for output
              cd10 = cdm_lnd(i)
              if ( wind10(i) .ge. 0.1 ) then
                cd10=cdm_lnd(i)* (wspd(i)/(0.01*wind10(i)) )**2
                !tmp9=0.01*abs(tzot(i))
                !ch_out(i)=ch_lnd(i)*(wspd(i)/(0.01*wind10(i)) ) * &
               !           (alog(zkmax(i)/tmp9)/alog(10.0/tmp9))
              end if
              fm10_lnd(i) = karman/sqrt(cd10)

              !GJF: conductances aren't used in other CCPP schemes, but this limit
              !  might be able to replace the limits on drag coefficients above

              !chs_lnd(i)=ch_lnd(i)*wspd (i) !conductance
              !chs2_lnd(i)=ustar_lnd(i)*karman/fh2_lnd(i) !2m conductance

              !!!2014-0922  cap CHS over land points
              ! chs_lnd(i)=amin1(chs_lnd(i), 0.05)
              ! chs2_lnd(i)=amin1(chs2_lnd(i), 0.05)
              ! if (chs2_lnd(i) < 0) chs2_lnd(i)=1.0e-6

              if (diag_qss) then
                esat = fpvs(tskin_lnd(i))
                qss_lnd(i) = ep2*esat/(psfc(i)-esat)
              end if

              !GJF: not used in CCPP
              !flhc_lnd(i)=cpm(i)*rho1(i)*chs_lnd(i)
              !flqc_lnd(i)=rho1(i)*chs_lnd(i)
              !cqs2_lnd(i)=chs2_lnd(i)
            end if !dry

            if (icy(i)) then
              !GJF: from WRF's module_sf_gfdl.F
              smcdry=drysmc(isltyp(i))
              smcmax=maxsmc(isltyp(i))
              wetc(i)=(smois(i,1)-smcdry)/(smcmax-smcdry)
              wetc(i)=amin1(1.,amax1(wetc(i),0.))


              !GJF: the lower boundary temperature passed in to MFLUX2 either follows GFS:
              tstrc(i) = 0.5*(tskin_ice(i) + tsurf_ice(i)) !averaging tskin_ice and tsurf_ice as in GFS surface layer breaks ntsflg functionality
              !GJF: or WRF module_sf_gfdl.F:
              !tstrc(i) = tskin_ice(i)
              !averaging tskin_ice and tsurf_ice as in GFS surface layer breaks ntsflg functionality

              !GJF: Roughness Length Limitation section
              !  The WRF version of module_sf_gfdl.F has no checks on the roughness lengths prior to entering MFLUX2.
              !  The following limits were placed on roughness lengths from the GFS surface layer scheme at the suggestion
              !  of Chunxi Zhang. Using the GFDL surface layer without such checks can lead to instability in the UFS.

              !znt_ice is in cm, z0max/ztmax are in m at this point
              z0max(i) = max(1.0e-6, min(0.01 * znt_ice(i), zkmax(i)))
  !** xubin's new z0  over land and sea ice
              tem1  = 1.0 - shdmax(i)
              tem2  = tem1 * tem1
              tem1  = 1.0  - tem2

              if( ivegsrc == 1 ) then
                z0max(i) = exp( tem2*log01 + tem1*log(z0max(i)) )
              elseif (ivegsrc == 2 ) then
                z0max(i) = exp( tem2*log01 + tem1*log(z0max(i)) )
              endif

              z0max(i) = max(z0max(i), 1.0e-6)

  !           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height
  !           dependance of czil
              czilc = 0.8

              tem1  = 1.0 - sigmaf(i)
              ztmax(i) = z0max(i)*exp( - tem1*tem1 &
       &                     * czilc*karman*sqrt(ustar_ice(i)*(0.01/1.5e-05)))
              ztmax(i) = max(ztmax(i), 1.0e-6)


              !GJF: from WRF's module_sf_gfdl.F
              if (wind10(i) <= 1.0e-10 .or. wind10(i) > 150.0) then
                 wind10(i)=wspd(i)*alog(10.0/z0max(i))/alog(z1(i)/z0max(i))
              end if
              wind10(i)=wind10(i)*100.0   !! m/s to cm/s

              ztmax(i) = ztmax(i)*100.0 !m to cm
              z0max(i) = z0max(i)*100.0 !m to cm

              call mflux2 (fxh(i), fxe(i), fxmx(i), fxmy(i), cdm_ice(i), rib_ice(i), &
                xxfh(i), ztmax(i), z0max(i), tstrc(i),   &
                pspc(i), pkmax(i), wetc(i), slwdc(i), z1_cm(i), icoef_sf, iwavecpl, lcurr_sf, charn(i), msang(i), &
                scurx(i), scury(i), pert_Cd, ens_random_seed, ens_Cdamp, upc(i), vpc(i), t1(i), q1(i), &
                dt, wind10(i), xxfh2(i), ntsflg, sfenth, tzot(i), ep2, errmsg, &
                errflg)
                if (errflg /= 0) return

              !GJF: this is broken when tstrc is set to an average of two variables
              if (ntsflg==1) then
                tskin_ice(i) = tstrc(i)      ! gopal's doing
              end if

              if (diag_wind10m) then
                u10_ice(i) = u1(i)*(0.01*wind10(i)/wspd(i))
                v10_ice(i) = v1(i)*(0.01*wind10(i)/wspd(i))
              end if

              !GJF: these variables are not needed in a GFS-based suite, but are found in WRF's module_sf_gfdl.F and kept in comments for legacy
              !gz1oz0(i) = alog(zkmax(i)/znt_ice(i))
              !taux(i) = fxmx(i)/10.    ! gopal's doing for Ocean coupling
              !tauy(i) = fxmy(i)/10.    ! gopal's doing for Ocean coupling

              cdm_ice(i) = max(cdm_ice(i), cd_low_limit)
              cdm_ice(i) = min(cdm_ice(i), cd_high_limit)
              fm_ice(i) = karman/sqrt(cdm_ice(i))

              !1) try fh_ice from MFLUX2
              fh_ice(i) = karman*xxfh(i)

              !2) calc ch_ice from fm_ice and fh_ice
              ch_ice(i)  = karman*karman/(fm_ice(i) * fh_ice(i))

              !3) check if ch_ice is out of bounds (if so, recalculate fh_ice from bounded value)
              ch_bound_excursion = .false.
              if (ch_ice(i) < ch_low_limit) then
                ch_bound_excursion = .true.
                ch_ice(i) = ch_low_limit
              else if (ch_ice(i) > ch_high_limit) then
                ch_bound_excursion = .true.
                ch_ice(i) = ch_high_limit
              end if

              fh2_ice(i) = karman*xxfh2(i)

              if (ch_bound_excursion) then
                fh2_fh_ratio = min(xxfh2(i)/xxfh(i), 1.0)
                fh_ice(i) = karman*karman/(fm_ice(i)*ch_ice(i))
                fh2_ice(i) = fh2_fh_ratio*fh_ice(i)
              end if

              !Other CCPP schemes (PBL) ask for fm/fh instead of psim/psih
              !psim_ice(i)=gz1oz0(i)-fm_ice(i)
              !psih_ice(i)=gz1oz0(i)-fh_ice(i)

              ustar_ice(i) = 0.01*sqrt(cdm_ice(i)*   &
                         (upc(i)*upc(i) + vpc(i)*vpc(i)))
              !GJF: from Chunxi Zhang's module_sf_sfclayrev.f90 (I'm not sure it's necessary.)
              ustar_ice(i) = amax1(ustar_ice(i),0.001)

              stress_ice(i) = cdm_ice(i)*wspd(i)*wspd(i)

              !GJF: from WRF's module_sf_gfdl.F
              !!! convert cd, ch to values at 10m, for output
              cd10 = cdm_ice(i)
              if ( wind10(i) .ge. 0.1 ) then
                cd10=cdm_ice(i)* (wspd(i)/(0.01*wind10(i)) )**2
                !tmp9=0.01*abs(tzot(i))
                !ch_out(i)=ch_ice(i)*(wspd(i)/(0.01*wind10(i)) ) * &
               !           (alog(zkmax(i)/tmp9)/alog(10.0/tmp9))
              end if
              fm10_ice(i) = karman/sqrt(cd10)

              !GJF: conductances aren't used in other CCPP schemes
              !chs_ice(i)=ch_ice(i)*wspd (i) !conductance
              !chs2_ice(i)=ustar_ice(i)*karman/fh2_ice(i) !2m conductance

              if (diag_qss) then
                esat = fpvs(tskin_ice(i))
                qss_ice(i) = ep2*esat/(psfc(i)-esat)
              end if

              !flhc_ice(i)=cpm(i)*rho1(i)*chs_ice(i)
              !flqc_ice(i)=rho1(i)*chs_ice(i)
              !cqs2_ice(i)=chs2_ice(i)
            end if !ice

            if (wet(i)) then
              wetc(i) = 1.0

              !GJF: the lower boundary temperature passed in to MFLUX2 either follows GFS:
              tstrc(i) = 0.5*(tskin_wat(i) + tsurf_wat(i)) !averaging tskin_wat and tsurf_wat as in GFS surface layer breaks ntsflg functionality
              !GJF: or WRF module_sf_gfdl.F:
              !tstrc(i) = tskin_wat(i)

              ! DH* 20201009: these bounds on ocean roughness lengths are from Chunxi Zhang's module_sf_sfclayrev.f90 (in cm)
              znt_wat(i)=min(2.85e-1,max(znt_wat(i),1.27e-5))

              !GJF: from WRF's module_sf_gfdl.F
              if (wind10(i) <= 1.0e-10 .or. wind10(i) > 150.0) then
                 wind10(i)=wspd(i)*alog(10.0/(0.01*znt_wat(i)))/alog(z1(i)/(0.01*znt_wat(i)))
              end if
              wind10(i)=wind10(i)*100.0   !! m/s to cm/s

              !GJF: mflux2 expects negative roughness length for ocean points
              znt_wat(i) = -znt_wat(i)

              call mflux2 (fxh(i), fxe(i), fxmx(i), fxmy(i), cdm_wat(i), rib_wat(i), &
                xxfh(i), znt_wat(i), mznt(i), tstrc(i),   &
                pspc(i), pkmax(i), wetc(i), slwdc(i), z1_cm(i), icoef_sf, iwavecpl, lcurr_sf, charn(i), msang(i), &
                scurx(i), scury(i), pert_Cd, ens_random_seed, ens_Cdamp, upc(i), vpc(i), t1(i), q1(i), &
                dt, wind10(i), xxfh2(i), ntsflg, sfenth, tzot(i), ep2, errmsg, &
                errflg)
                if (errflg /= 0) return

              !GJF: this is broken when tstrc is set to an average of two variables
              if (ntsflg==1) then
                tskin_wat(i) = tstrc(i)      ! gopal's doing
              end if

              znt_wat(i)= abs(znt_wat(i))
              mznt(i)= abs(mznt(i))

              !GJF: these bounds on ocean roughness lengths are from Chunxi Zhang's module_sf_sfclayrev.f90 (in cm)
              znt_wat(i)=min(2.85e-1,max(znt_wat(i),1.27e-5))

              if (diag_wind10m) then
                u10_ocn(i) = u1(i)*(0.01*wind10(i)/wspd(i))
                v10_ocn(i) = v1(i)*(0.01*wind10(i)/wspd(i))
              end if

              !GJF: these variables are not needed in a GFS-based suite, but are found in WRF's module_sf_gfdl.F and kept in comments for legacy
              !gz1oz0(i) = alog(zkmax(i)/znt_wat(i))
              !taux(i) = fxmx(i)/10.    ! gopal's doing for Ocean coupling
              !tauy(i) = fxmy(i)/10.    ! gopal's doing for Ocean coupling

              cdm_wat(i) = max(cdm_wat(i), cd_low_limit)
              cdm_wat(i) = min(cdm_wat(i), cd_high_limit)
              fm_wat(i) = karman/sqrt(cdm_wat(i))

              !1) try fh_wat from MFLUX2
              fh_wat(i) = karman*xxfh(i)

              !2) calc ch_wat from fm_wat and fh_wat
              ch_wat(i)  = karman*karman/(fm_wat(i) * fh_wat(i))

              !3) check if ch_lnd is out of bounds (if so, recalculate fh_lnd from bounded value)
              ch_bound_excursion = .false.
              if (ch_wat(i) < ch_low_limit) then
                ch_bound_excursion = .true.
                ch_wat(i) = ch_low_limit
              else if (ch_wat(i) > ch_high_limit) then
                ch_bound_excursion = .true.
                ch_wat(i) = ch_high_limit
              end if

              fh2_wat(i) = karman*xxfh2(i)

              if (ch_bound_excursion) then
                fh2_fh_ratio = min(xxfh2(i)/xxfh(i), 1.0)
                fh_wat(i) = karman*karman/(fm_wat(i)*ch_wat(i))
                fh2_wat(i) = fh2_fh_ratio*fh_wat(i)
              end if

              !Other CCPP schemes (PBL) ask for fm/fh instead of psim/psih
              !psim_ocn(i)=gz1oz0(i)-fm_wat(i)
              !psih_ocn(i)=gz1oz0(i)-fh_wat(i)

              ustar_wat(i) = 0.01*sqrt(cdm_wat(i)*   &
                         (upc(i)*upc(i) + vpc(i)*vpc(i)))
              !GJF: from Chunxi Zhang's module_sf_sfclayrev.f90 (I'm not sure it's necessary.)
              ustar_wat(i) = amax1(ustar_wat(i),0.001)

              stress_wat(i) = cdm_wat(i)*wspd(i)*wspd(i)

              !GJF: from WRF's module_sf_gfdl.F
              !!! convert cd, ch to values at 10m, for output
              cd10 = cdm_wat(i)
              if ( wind10(i) .ge. 0.1 ) then
                cd10=cdm_wat(i)* (wspd(i)/(0.01*wind10(i)) )**2
                !tmp9=0.01*abs(tzot(i))
                !ch_out(i)=ch_wat(i)*(wspd(i)/(0.01*wind10(i)) ) * &
               !           (alog(zkmax(i)/tmp9)/alog(10.0/tmp9))
              end if
              fm10_wat(i) = karman/sqrt(cd10)

              !GJF: conductances aren't used in other CCPP schemes
              !chs_ocn(i)=ch_wat(i)*wspd (i) !conductance
              !chs2_ocn(i)=ustar_wat(i)*karman/fh2_wat(i) !2m conductance

              if (diag_qss) then
                esat = fpvs(tskin_wat(i))
                qss_wat(i) = ep2*esat/(psfc(i)-esat)
              end if
            end if !wet

            !flhc_ocn(i)=cpm(i)*rho1(i)*chs_ocn(i)
            !flqc_ocn(i)=rho1(i)*chs_ocn(i)
            !cqs2_ocn(i)=chs2_ocn(i)
          end if !flag_iter
        end do

        !GJF: this code has not been updated since GFS suites don't require this; one would need to have different values of hfx, qfx, lh for each surface type
        ! if (isfflx.eq.0) then
        !   do i=its,ite
        !     hfx(i)=0.
        !     lh(i)=0.
        !     qfx(i)=0.
        !   enddo
        ! else
        !   do i=its,ite
        !     if(islmsk == 0) then
        !       !water
        !       hfx(i)= -10.*cp*fxh(i)
        !     else if (islmsk == 1) then
        !       hfx(i)= -10.*cp*fxh(i)
        !       hfx(i)=amax1(hfx(i),-250.)
        !     end if
        !     qfx(j)=-10.*fxe(i)
        !     qfx(i)=amax1(qfx(i),0.)
        !     lh(i)=xlv*qfx(i)
        !   enddo
        ! endif


        end subroutine gfdl_sfc_layer_run

!---------------------------------
!GJF (2020/04/21): The starting point for the MFLUX2 subroutine here was module_sf_gfdl.F in WRF
      SUBROUTINE MFLUX2( fxh,fxe,fxmx,fxmy,cdm,rib,xxfh,zoc,mzoc,tstrc,        &    !mzoc KWON
                         pspc,pkmax,wetc,slwdc,z1,                             &
                         icoef_sf,iwavecpl,lcurr_sf,alpha,gamma,xcur,ycur,     &
                         pert_Cd, ens_random_seed, ens_Cdamp,                  &
                         upc,vpc,tpc,rpc,dt,wind10,xxfh2,ntsflg,sfenth,        &
                         tzot, ep2, errmsg, errflg)

!------------------------------------------------------------------------
!
!     MFLUX2 computes surface fluxes of momentum, heat,and moisture
!     using monin-obukhov. the roughness length "z0" is prescribed
!     over land and over ocean "z0" is computed using charnocks formula.
!     the universal functions (from similarity theory approach) are
!     those of hicks. This is Bob's doing.
!
!------------------------------------------------------------------------

      USE module_sf_exchcoef
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     user interface variables
!-----------------------------------------------------------------------
      !GJF: This subroutine was converted to expect data from a single point instead of a horizontal array to accommodate a fractional landmask
      !integer,intent(in)  :: ims,ime
      !integer,intent(in)  :: its,ite
      integer, parameter :: ims = 1
      integer, parameter :: ime = 1
      integer, parameter :: its = 1
      integer, parameter :: ite = 1
      integer,intent(in)  :: ntsflg
      integer,intent(in)  :: icoef_sf
      integer,intent(in)  :: iwavecpl
      logical,intent(in)  :: lcurr_sf
      logical,intent(in)  :: pert_Cd
      integer,intent(in)  :: ens_random_seed
      real(kind=kind_phys),intent(in)     :: ens_Cdamp

      real(kind=kind_phys), intent (out), dimension (ims :ime ) :: fxh
      real(kind=kind_phys), intent (out), dimension (ims :ime ) :: fxe
      real(kind=kind_phys), intent (out), dimension (ims :ime ) :: fxmx
      real(kind=kind_phys), intent (out), dimension (ims :ime ) :: fxmy
      real(kind=kind_phys), intent (inout), dimension (ims :ime ) :: cdm
!       real, intent (out), dimension (ims :ime ) :: cdm2
      real(kind=kind_phys), intent (out), dimension (ims :ime ) :: rib
      real(kind=kind_phys), intent (out), dimension (ims :ime ) :: xxfh
      real(kind=kind_phys), intent (out), dimension (ims :ime ) :: xxfh2
      real(kind=kind_phys), intent (out), dimension (ims :ime ) :: wind10

      real(kind=kind_phys), intent ( inout), dimension (ims :ime ) :: zoc,mzoc    !KWON
      real(kind=kind_phys), intent ( inout), dimension (ims :ime ) :: tzot        !WANG
      real(kind=kind_phys), intent ( inout), dimension (ims :ime ) :: tstrc

      real(kind=kind_phys), intent ( in)                        :: dt
      real(kind=kind_phys), intent ( in)                        :: sfenth
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: pspc
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: pkmax
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: wetc
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: slwdc
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: alpha, gamma
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: xcur, ycur
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: z1

      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: upc
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: vpc
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: tpc
      real(kind=kind_phys), intent ( in), dimension (ims :ime ) :: rpc

      real(kind=kind_phys), intent ( in) :: ep2

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!-----------------------------------------------------------------------
!     internal variables
!-----------------------------------------------------------------------

      integer, parameter :: icntx = 30

      integer, dimension(1   :ime) :: ifz
      integer, dimension(1   :ime) :: indx
      integer, dimension(1   :ime) :: istb
      integer, dimension(1   :ime) :: it
      integer, dimension(1   :ime) :: iutb

      real(kind=kind_phys), dimension(1   :ime) :: aap
      real(kind=kind_phys), dimension(1   :ime) :: bq1
      real(kind=kind_phys), dimension(1   :ime) :: bq1p
      real(kind=kind_phys), dimension(1   :ime) :: delsrad
      real(kind=kind_phys), dimension(1   :ime) :: ecof
      real(kind=kind_phys), dimension(1   :ime) :: ecofp
      real(kind=kind_phys), dimension(1   :ime) :: estso
      real(kind=kind_phys), dimension(1   :ime) :: estsop
      real(kind=kind_phys), dimension(1   :ime) :: fmz1
      real(kind=kind_phys), dimension(1   :ime) :: fmz10
      real(kind=kind_phys), dimension(1   :ime) :: fmz2
      real(kind=kind_phys), dimension(1   :ime) :: fmzo1
      real(kind=kind_phys), dimension(1   :ime) :: foft
      real(kind=kind_phys), dimension(1   :ime) :: foftm
      real(kind=kind_phys), dimension(1   :ime) :: frac
      real(kind=kind_phys), dimension(1   :ime) :: land
      real(kind=kind_phys), dimension(1   :ime) :: pssp
      real(kind=kind_phys), dimension(1   :ime) :: qf
      real(kind=kind_phys), dimension(1   :ime) :: rdiff
      real(kind=kind_phys), dimension(1   :ime) :: rho
      real(kind=kind_phys), dimension(1   :ime) :: rkmaxp
      real(kind=kind_phys), dimension(1   :ime) :: rstso
      real(kind=kind_phys), dimension(1   :ime) :: rstsop
      real(kind=kind_phys), dimension(1   :ime) :: sf10
      real(kind=kind_phys), dimension(1   :ime) :: sf2
      real(kind=kind_phys), dimension(1   :ime) :: sfm
      real(kind=kind_phys), dimension(1   :ime) :: sfzo
      real(kind=kind_phys), dimension(1   :ime) :: sgzm
      real(kind=kind_phys), dimension(1   :ime) :: slwa
      real(kind=kind_phys), dimension(1   :ime) :: szeta
      real(kind=kind_phys), dimension(1   :ime) :: szetam
      real(kind=kind_phys), dimension(1   :ime) :: t1
      real(kind=kind_phys), dimension(1   :ime) :: t2
      real(kind=kind_phys), dimension(1   :ime) :: tab1
      real(kind=kind_phys), dimension(1   :ime) :: tab2
      real(kind=kind_phys), dimension(1   :ime) :: tempa1
      real(kind=kind_phys), dimension(1   :ime) :: tempa2
      real(kind=kind_phys), dimension(1   :ime) :: theta
      real(kind=kind_phys), dimension(1   :ime) :: thetap
      real(kind=kind_phys), dimension(1   :ime) :: tsg
      real(kind=kind_phys), dimension(1   :ime) :: tsm
      real(kind=kind_phys), dimension(1   :ime) :: tsp
      real(kind=kind_phys), dimension(1   :ime) :: tss
      real(kind=kind_phys), dimension(1   :ime) :: ucom
      real(kind=kind_phys), dimension(1   :ime) :: uf10
      real(kind=kind_phys), dimension(1   :ime) :: uf2
      real(kind=kind_phys), dimension(1   :ime) :: ufh
      real(kind=kind_phys), dimension(1   :ime) :: ufm
      real(kind=kind_phys), dimension(1   :ime) :: ufzo
      real(kind=kind_phys), dimension(1   :ime) :: ugzm
      real(kind=kind_phys), dimension(1   :ime) :: uzeta
      real(kind=kind_phys), dimension(1   :ime) :: uzetam
      real(kind=kind_phys), dimension(1   :ime) :: vcom
      real(kind=kind_phys), dimension(1   :ime) :: vrtkx
      real(kind=kind_phys), dimension(1   :ime) :: vrts
      real(kind=kind_phys), dimension(1   :ime) :: wind
      real(kind=kind_phys), dimension(1   :ime) :: windp
      real(kind=kind_phys), dimension(1   :ime) :: wind10p  !WANG, 10m wind previous step
      real(kind=kind_phys), dimension(1   :ime) :: uvs1
!      real(kind=kind_phys), dimension(1   :ime) :: xxfh
      real(kind=kind_phys), dimension(1   :ime) :: xxfm
      real(kind=kind_phys), dimension(1   :ime) :: xxsh
      real(kind=kind_phys), dimension(1   :ime) :: z10
      real(kind=kind_phys), dimension(1   :ime) :: z2
      real(kind=kind_phys), dimension(1   :ime) :: zeta
      real(kind=kind_phys), dimension(1   :ime) :: zkmax

      real(kind=kind_phys), dimension(1   :ime) :: pss
      real(kind=kind_phys), dimension(1   :ime) :: tstar
      real(kind=kind_phys), dimension(1   :ime) :: ukmax
      real(kind=kind_phys), dimension(1   :ime) :: vkmax
      real(kind=kind_phys), dimension(1   :ime) :: tkmax
      real(kind=kind_phys), dimension(1   :ime) :: rkmax
      real(kind=kind_phys), dimension(1   :ime) :: zot
      real(kind=kind_phys), dimension(1   :ime) :: fhzo1
      real(kind=kind_phys), dimension(1   :ime) :: sfh

      real(kind=kind_phys) :: ux13, yo, y,xo,x,ux21,ugzzo,ux11,ux12,uzetao,xnum,alll
      real(kind=kind_phys) :: ux1,ugz,x10,uzo,uq,ux2,ux3,xtan,xden,y10,uzet1o,ugz10
      real(kind=kind_phys) :: szet2, zal2,ugz2
      real(kind=kind_phys) :: rovcp,boycon,cmo2,psps1,zog,enrca,rca,cmo1,amask,en,ca,a,c
      real(kind=kind_phys) :: sgz,zal10,szet10,fmz,szo,sq,fmzo,rzeta1,zal1g,szetao,rzeta2,zal2g
      real(kind=kind_phys) :: hcap,xks,pith,teps,diffot,delten,alevp,psps2,alfus,nstep
      real(kind=kind_phys) :: shfx,sigt4,reflect
      real(kind=kind_phys) :: cor1,cor2,szetho,zal2gh,cons_p000001,cons_7,vis,ustar,restar,rat
      real(kind=kind_phys) :: wndm,ckg
      real(kind=kind_phys) :: windmks,znott,znotm
      real(kind=kind_phys) :: ubot, vbot
      integer:: i,j,ii,iq,nnest,icnt,ngd,ip

!-----------------------------------------------------------------------
!     internal variables
!-----------------------------------------------------------------------

      real(kind=kind_phys), dimension (223) :: tab
      real(kind=kind_phys), dimension (223) :: table
      real(kind=kind_phys), dimension (101) :: tab11
      real(kind=kind_phys), dimension (41) :: table4
      real(kind=kind_phys), dimension (42) :: tab3
      real(kind=kind_phys), dimension (54) :: table2
      real(kind=kind_phys), dimension (54) :: table3
      real(kind=kind_phys), dimension (74) :: table1
      real(kind=kind_phys), dimension (80) :: tab22

      character(len=255) :: message

      equivalence (tab(1),tab11(1))
      equivalence (tab(102),tab22(1))
      equivalence (tab(182),tab3(1))
      equivalence (table(1),table1(1))
      equivalence (table(75),table2(1))
      equivalence (table(129),table3(1))
      equivalence (table(183),table4(1))

      data amask/ -98.0/
!-----------------------------------------------------------------------
!     tables used to obtain the vapor pressures or saturated vapor
!     pressure
!-----------------------------------------------------------------------

      data tab11/21*0.01403,0.01719,0.02101,0.02561,0.03117,0.03784,      &
     &.04584,.05542,.06685,.08049,.09672,.1160,.1388,.1658,.1977,.2353,   &
     &.2796,.3316,.3925,.4638,.5472,.6444,.7577,.8894,1.042,1.220,1.425,  &
     &1.662,1.936,2.252,2.615,3.032,3.511,4.060,4.688,5.406,6.225,7.159,  &
     &8.223,9.432,10.80,12.36,14.13,16.12,18.38,20.92,23.80,27.03,30.67,  &
     &34.76,39.35,44.49,50.26,56.71,63.93,71.98,80.97,90.98,102.1,114.5,  &
     &128.3,143.6,160.6,179.4,200.2,223.3,248.8,276.9,307.9,342.1,379.8,  &
     &421.3,466.9,517.0,572.0,632.3,698.5,770.9,850.2,937.0,1032./

      data tab22/1146.6,1272.0,1408.1,1556.7,1716.9,1890.3,2077.6,2279.6  &
     &,2496.7,2729.8,2980.0,3247.8,3534.1,3839.8,4164.8,4510.5,4876.9,    &
     &5265.1,5675.2,6107.8,6566.2,7054.7,7575.3,8129.4,8719.2,9346.5,     &
     &10013.,10722.,11474.,12272.,13119.,14017.,14969.,15977.,17044.,     &
     &18173.,19367.,20630.,21964.,23373.,24861.,26430.,28086.,29831.,     &
     &31671.,33608.,35649.,37796.,40055.,42430.,44927.,47551.,50307.,     &
     &53200.,56236.,59422.,62762.,66264.,69934.,73777.,77802.,82015.,     &
     &86423.,91034.,95855.,100890.,106160.,111660.,117400.,123400.,       &
     &129650.,136170.,142980.,150070.,157460.,165160.,173180.,181530.,    &
     &190220.,199260./

      data tab3/208670.,218450.,228610.,239180.,250160.,261560.,273400.,  &
     &285700.,298450.,311690.,325420.,339650.,354410.,369710.,385560.,    &
     &401980.,418980.,436590.,454810.,473670.,493170.,513350.,534220.,    &
     &555800.,578090.,601130.,624940.,649530.,674920.,701130.,728190.,    &
     &756110.,784920.,814630.,845280.,876880.,909450.,943020.,977610.,    &
     &1013250.,1049940.,1087740./

      data table1/20*0.0,.3160e-02,.3820e-02,.4600e-02,.5560e-02,.6670e-02, &
     & .8000e-02,.9580e-02,.1143e-01,.1364e-01,.1623e-01,.1928e-01,       &
     &.2280e-01,.2700e-01,.3190e-01,.3760e-01,.4430e-01,.5200e-01,          &
     &.6090e-01,.7130e-01,.8340e-01,.9720e-01,.1133e+00,.1317e-00,          &
     &.1526e-00,.1780e-00,.2050e-00,.2370e-00,.2740e-00,.3160e-00,          &
     &.3630e-00,.4170e-00,.4790e-00,.5490e-00,.6280e-00,.7180e-00,          &
     &.8190e-00,.9340e-00,.1064e+01,.1209e+01,.1368e+01,.1560e+01,          &
     &.1770e+01,.1990e+01,.2260e+01,.2540e+01,.2880e+01,.3230e+01,          &
     &.3640e+01,.4090e+01,.4590e+01,.5140e+01,.5770e+01,.6450e+01,          &
     &.7220e+01/

      data table2/.8050e+01,.8990e+01,.1001e+02,.1112e+02,.1240e+02,      &
     &.1380e+02,.1530e+02,.1700e+02,.1880e+02,.2080e+02,.2310e+02,        &
     &.2550e+02,.2810e+02,.3100e+02,.3420e+02,.3770e+02,.4150e+02,        &
     &.4560e+02,.5010e+02,.5500e+02,.6030e+02,.6620e+02,.7240e+02,        &
     &.7930e+02,.8680e+02,.9500e+02,.1146e+03,.1254e+03,.1361e+03,        &
     &.1486e+03,.1602e+03,.1734e+03,.1873e+03,.2020e+03,.2171e+03,        &
     &.2331e+03,.2502e+03,.2678e+03,.2863e+03,.3057e+03,.3250e+03,        &
     &.3457e+03,.3664e+03,.3882e+03,.4101e+03,.4326e+03,.4584e+03,        &
     &.4885e+03,.5206e+03,.5541e+03,.5898e+03,.6273e+03,.6665e+03,        &
     &.7090e+03/

      data table3/.7520e+03,.7980e+03,.8470e+03,.8980e+03,.9520e+03,      &
     &.1008e+04,.1067e+04,.1129e+04,.1194e+04,.1263e+04,.1334e+04,        &
     &.1409e+04,.1488e+04,.1569e+04,.1656e+04,.1745e+04,.1840e+04,        &
     &.1937e+04,.2041e+04,.2147e+04,.2259e+04,.2375e+04,.2497e+04,        &
     &.2624e+04,.2756e+04,.2893e+04,.3036e+04,.3186e+04,.3340e+04,        &
     &.3502e+04,.3670e+04,.3843e+04,.4025e+04,.4213e+04,.4408e+04,        &
     &.4611e+04,.4821e+04,.5035e+04,.5270e+04,.5500e+04,.5740e+04,        &
     &.6000e+04,.6250e+04,.6520e+04,.6810e+04,.7090e+04,.7390e+04,        &
     &.7700e+04,.8020e+04,.8350e+04,.8690e+04,.9040e+04,.9410e+04,        &
     &.9780e+04/

      data table4/.1016e+05,.1057e+05,.1098e+05,.1140e+05,.1184e+05,      &
     &.1230e+05,.1275e+05,.1324e+05,.1373e+05,.1423e+05,.1476e+05,        &
     &.1530e+05,.1585e+05,.1642e+05,.1700e+05,.1761e+05,.1822e+05,        &
     &.1886e+05,.1950e+05,.2018e+05,.2087e+05,.2158e+05,.2229e+05,        &
     &.2304e+05,.2381e+05,.2459e+05,.2539e+05,.2621e+05,.2706e+05,        &
     &.2792e+05,.2881e+05,.2971e+05,.3065e+05,.3160e+05,.3257e+05,        &
     &.3357e+05,.3459e+05,.3564e+05,.3669e+05,.3780e+05,.0000e+00/
!
! spcify constants needed by MFLUX2
!
!GJF: should send through argument list, but these have nonstandard units
      real,parameter :: cp    = 1.00464e7
      real,parameter :: g     = 980.6
      real,parameter :: rgas  = 2.87e6
      real,parameter :: og    = 1./g
      integer :: ntstep = 0

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
#if HWRF==1
      real*8 :: gasdev,ran1          !zhang
      real :: rr                     !zhang
      logical,save :: pert_Cd_local            !zhang
      CHARACTER(len=3) :: env_memb,env_pp
      integer,save :: ens_random_seed_local,env_pp_local         !zhang
      integer :: ensda_physics_pert !zhang
      real,save :: ens_Cdamp_local         !zhang
      data ens_random_seed_local/0/
      data env_pp_local/0/
      if ( ens_random_seed_local .eq. 0 ) then
         CALL nl_get_ensda_physics_pert(1,ensda_physics_pert)
         ens_random_seed_local=ens_random_seed
         env_pp_local=ensda_physics_pert
         pert_Cd_local=.false.
         ens_Cdamp_local=0.0
! env_pp=1: do physics perturbations for ensda members, ens_random_seed must be 99
         if ( env_pp_local .eq. 1 ) then
            if ( ens_random_seed .ne. 99 ) then
               pert_Cd_local=.true.
               ens_Cdamp_local=ens_Cdamp
            else
! ens_random_seed=99 do physics perturbation for ensemble forecasts, env_pp must be zero
               ens_random_seed_local=ens_random_seed
               pert_Cd_local=pert_Cd
               ens_Cdamp_local=ens_Cdamp
            endif
         else
            ens_random_seed_local=ens_random_seed
            pert_Cd_local=pert_Cd
            ens_Cdamp_local=ens_Cdamp
         endif
      print*, "Cd ===", ens_random_seed_local,pert_Cd_local,ens_Cdamp_local,ensda_physics_pert
      endif
#endif

!     character*10 routine
!     routine = 'mflux2'
!
!------------------------------------------------------------------------
!     set water availability constant "ecof" and land mask "land".
!     limit minimum wind speed to 100 cm/s
!------------------------------------------------------------------------
!  constants for 10 m winds  (correction for knots
!
      cor1 = .120
      cor2 = 720.
! KWON : remove the artificial increase of 10m wind speed over 60kts
!        which comes from GFDL hurricane model
      cor1 = 0.
      cor2 = 0.
!

      do i = its,ite
        z10(i) = 1000.
        z2 (i) =  200.
        pss(i) = pspc(i)
        tstar(i) = tstrc(i)

        if ( lcurr_sf .and. zoc(i) .le. 0.0 ) then
          ubot = upc(i)  - xcur(i) * 100.0
          vbot = vpc(i)  - ycur(i) * 100.0
!          ubot = upc(i)
!          vbot = vpc(i)
        else
          ubot = upc(i)
          vbot = vpc(i)
        endif
        uvs1(i)= amax1( SQRT(ubot*ubot +    &
                             vbot*vbot), 100.0)
        if ( iwavecpl .eq. 1 .and. zoc(i) .le. 0.0 ) then
          ukmax(i) = ( ubot * cos(gamma(i))  -          &
                      vbot * sin(gamma(i)) )            &
                                  * cos(gamma(i))
          vkmax(i) = ( vbot * cos(gamma(i))  -          &
                      ubot * sin(gamma(i)) )            &
                                  * cos(gamma(i))

        else
          ukmax(i) = ubot
          vkmax(i) = vbot
        endif

!       ukmax(i) = upc(i)
!       vkmax(i) = vpc(i)
        tkmax(i) = tpc(i)
        rkmax(i) = rpc(i)
      enddo

      do i = its,ite
        windp(i) = SQRT(ukmax(i)*ukmax(i) + vkmax(i)*vkmax(i))
        wind (i) = amax1(windp(i),100.)

!! use wind10 previous step
        wind10p(i) = wind10(i)  !! cm/s
        wind10p(i) = amax1(wind10p(i),100.)
!!

        if (zoc(i) .LT. amask) zoc(i) = -0.0185*0.001*wind10p(i)*wind10p(i)*og
        if (zoc(i) .GT. 0.0) then
          ecof(i) = wetc(i)
          land(i) = 1.0
          zot (i) = zoc(i)
        else
          ecof(i) = wetc(i)
          land(i) = 0.0
          windmks=wind10p(i)*.01
          if ( iwavecpl .eq. 1 ) then
            call znot_wind10m(windmks,znott,znotm,icoef_sf)
            !Check if Charnock parameter ratio is received in a proper range.
            if ( alpha(i) .ge. 0.2 .and. alpha(i) .le. 5. ) then
              znotm = znotm*alpha(i)
            endif
            zoc(i)  = -100.*znotm
            zot(i) =  -100* znott
          else
            call znot_wind10m(windmks,znott,znotm,icoef_sf)
            zoc(i)  = -100.*znotm
            zot(i) =  -100* znott
          endif
        endif
!------------------------------------------------------------------------
!     where necessary modify zo values over ocean.
!------------------------------------------------------------------------
!
        mzoc(i) = zoc(i)                !FOR SAVE MOMENTUM Zo
        tzot(i) = zot(i)                 !output wang
      enddo

!------------------------------------------------------------------------
!     define constants:
!     a and c = constants used in evaluating universal function for
!               stable case
!     ca      = karmen constant
!     cm01    = constant part of vertical integral of universal
!               function; stable case ( 0.5 < zeta < or = 10.0)
!     cm02    = constant part of vertical integral of universal
!               function; stable case ( zeta > 10.0)
!------------------------------------------------------------------------

      en     = 2.
      c      = .76
      a      = 5.
      ca     = .4
      cmo1   = .5*a - 1.648
      cmo2   = 17.193 + .5*a - 10.*c
      boycon = .61
      rovcp=rgas/cp

      do i = its,ite
        theta(i) = tkmax(i)/((pkmax(i)/pspc(i))**rovcp)
        vrtkx(i) = 1.0 + boycon*rkmax(i)
        !zkmax(i) = -rgas*tkmax(i)*alog(pkmax(i)/pspc(i))*og
        zkmax(i) = z1(i) !use precalculated height of first model layer center
      enddo

!------------------------------------------------------------------------
!     get saturation mixing ratios at surface
!------------------------------------------------------------------------

      do i = its,ite
        tsg  (i) = tstar(i)
        tab1 (i) = tstar(i) - 153.16
        it   (i) = IFIX(tab1(i))
        tab2 (i) = tab1(i) - FLOAT(it(i))
        t1   (i) = tab(min(223,max(1,it(i) + 1)))
        t2   (i) = table(min(223,max(1,it(i) + 1)))
        estso(i) = t1(i) + tab2(i)*t2(i)
         psps1 = (pss(i) - estso(i))
          if(psps1 .EQ. 0.0)then
           psps1 = .1
          endif
        rstso(i) = ep2*estso(i)/psps1
        vrts (i) = 1. + boycon*ecof(i)*rstso(i)
      enddo

!------------------------------------------------------------------------
!     check if consideration of virtual temperature changes stability.
!     if so, set "dthetav" to near neutral value (1.0e-4). also check
!     for very small lapse rates; if ABS(tempa1) <1.0e-4 then
!     tempa1=1.0e-4
!------------------------------------------------------------------------

      do i = its,ite
        tempa1(i) = theta(i)*vrtkx(i) - tstar(i)*vrts(i)
        tempa2(i) = tempa1(i)*(theta(i) - tstar(i))
        if (tempa2(i) .LT. 0.) tempa1(i) = 1.0e-4
        tab1(i) = ABS(tempa1(i))
        if (tab1(i) .LT. 1.0e-4) tempa1(i) = 1.0e-4
!------------------------------------------------------------------------
!     compute bulk richardson number "rib" at each point. if "rib"
!     exceeds 95% of critical richardson number "tab1" then "rib = tab1"
!------------------------------------------------------------------------

        rib (i) = g*zkmax(i)*tempa1(i)/                             &
                                    (tkmax(i)*vrtkx(i)*wind(i)*wind(i))
        tab2(i) = ABS(zoc(i))
        tab1(i) = 0.95/(c*(1. - tab2(i)/zkmax(i)))
        if (rib(i) .GT. tab1(i)) rib(i) = tab1(i)
      enddo

      do i = its,ite
        zeta(i) = ca*rib(i)/0.03
      enddo

!------------------------------------------------------------------------
!     begin looping through points on line, solving wegsteins iteration
!     for zeta at each point, and using hicks functions
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!     set initial guess of zeta=non - dimensional height "szeta" for
!     stable points
!------------------------------------------------------------------------

      rca   = 1./ca
      enrca = en*rca
!     turn off interfacial layer by zeroing out enrca
      enrca = 0.0
      zog   = .0185*og

!------------------------------------------------------------------------
!     stable points
!------------------------------------------------------------------------

      ip    = 0
      do i = its,ite
        if (zeta(i) .GE. 0.0) then
          ip = ip + 1
          istb(ip) = i
        endif
      enddo

      if (ip .EQ. 0) go to 170
      do i = 1,ip
        szetam(i) = 1.0e+30
        sgzm(i)   = 0.0e+00
        szeta(i)  = zeta(istb(i))
        ifz(i)    = 1
      enddo

!------------------------------------------------------------------------
!     begin wegstein iteration for "zeta" at stable points using
!     hicks(1976)
!------------------------------------------------------------------------

      do icnt = 1,icntx
        do i = 1,ip
          if (ifz(i) .EQ. 0) go to 80
          zal1g = ALOG(szeta(i))
          if (szeta(i) .LE. 0.5) then
            fmz1(i) = (zal1g + a*szeta(i))*rca
          else if (szeta(i) .GT. 0.5 .AND. szeta(i) .LE. 10.) then
            rzeta1  = 1./szeta(i)
            fmz1(i) = (8.*zal1g + 4.25*rzeta1 - &
                                          0.5*rzeta1*rzeta1 + cmo1)*rca
          else if (szeta(i) .GT. 10.) then
            fmz1(i) = (c*szeta(i) + cmo2)*rca
          endif
          szetao = ABS(zoc(istb(i)))/zkmax(istb(i))*szeta(i)
          zal2g  = ALOG(szetao)
          if (szetao .LE. 0.5) then
            fmzo1(i) = (zal2g + a*szetao)*rca
            sfzo (i) = 1. + a*szetao
          else if (szetao .GT. 0.5 .AND. szetao .LE. 10.) then
            rzeta2   = 1./szetao
            fmzo1(i) = (8.*zal2g + 4.25*rzeta2 - &
                                          0.5*rzeta2*rzeta2 + cmo1)*rca
            sfzo (i) = 8.0 - 4.25*rzeta2 + rzeta2*rzeta2
          else if (szetao .GT. 10.) then
            fmzo1(i) = (c*szetao + cmo2)*rca
            sfzo (i) = c*szetao
          endif


!        compute heat & moisture parts of zot.. for calculation of sfh

          szetho = ABS(zot(istb(i)))/zkmax(istb(i))*szeta(i)
          zal2gh = ALOG(szetho)
          if (szetho .LE. 0.5) then
            fhzo1(i) = (zal2gh + a*szetho)*rca
            sfzo (i) = 1. + a*szetho
          else if (szetho .GT. 0.5 .AND. szetho .LE. 10.) then
            rzeta2   = 1./szetho
            fhzo1(i) = (8.*zal2gh + 4.25*rzeta2 -   &
                                          0.5*rzeta2*rzeta2 + cmo1)*rca
            sfzo (i) = 8.0 - 4.25*rzeta2 + rzeta2*rzeta2
          else if (szetho .GT. 10.) then
            fhzo1(i) = (c*szetho + cmo2)*rca
            sfzo (i) = c*szetho
          endif

!------------------------------------------------------------------------
!      compute universal function at 10 meters for diagnostic purposes
!------------------------------------------------------------------------

           szet10 = ABS(z10(istb(i)))/zkmax(istb(i))*szeta(i)
           zal10  = ALOG(szet10)
           if (szet10 .LE. 0.5) then
             fmz10(i) = (zal10 + a*szet10)*rca
           else if (szet10 .GT. 0.5 .AND. szet10 .LE. 10.) then
             rzeta2   = 1./szet10
             fmz10(i) = (8.*zal10 + 4.25*rzeta2 - &
                                         0.5*rzeta2*rzeta2 + cmo1)*rca
           else if (szet10 .GT. 10.) then
             fmz10(i) = (c*szet10 + cmo2)*rca
           endif
           sf10(i) = fmz10(i) - fmzo1(i)
!          compute 2m values for diagnostics in HWRF
           szet2  = ABS(z2 (istb(i)))/zkmax(istb(i))*szeta(i)
           zal2   = ALOG(szet2 )
           if (szet2  .LE. 0.5) then
             fmz2 (i) = (zal2  + a*szet2 )*rca
           else if (szet2  .GT. 0.5 .AND. szet2  .LE. 2.) then
             rzeta2   = 1./szet2
             fmz2 (i) = (8.*zal2  + 4.25*rzeta2 - &
                                         0.5*rzeta2*rzeta2 + cmo1)*rca
           else if (szet2  .GT. 2.) then
             fmz2 (i) = (c*szet2  + cmo2)*rca
           endif
           sf2 (i) = fmz2 (i) - fmzo1(i)

           sfm(i) = fmz1(i) - fmzo1(i)
           sfh(i) = fmz1(i) - fhzo1(i)
           sgz    = ca*rib(istb(i))*sfm(i)*sfm(i)/ &
                                                (sfh(i) + enrca*sfzo(i))
           fmz    = (sgz - szeta(i))/szeta(i)
           fmzo   = ABS(fmz)
           if (fmzo .GE. 5.0e-5) then
             sq        = (sgz - sgzm(i))/(szeta(i) - szetam(i))
             if(sq .EQ. 1) then
              write(errmsg,'(*(a))') 'NCO ERROR DIVIDE BY ZERO IN gfdl_sfc_layer.F90/MFLUX2 (STABLE CASE)'// &
                  'sq is 1 ',fmzo,sgz,sgzm(i),szeta(i),szetam(i)
              errflg = 1
              return
             endif
             szetam(i) = szeta(i)
             szeta (i) = (sgz - szeta(i)*sq)/(1.0 - sq)
             sgzm  (i) = sgz
           else
             ifz(i) = 0
           endif
80      continue
        enddo
      enddo

      do i = 1,ip
        if (ifz(i) .GE. 1) go to 110
      enddo

      go to 130

110   continue

      write(errmsg,'(*(a))') 'NON-CONVERGENCE FOR STABLE ZETA IN gfdl_sfc_layer.F90/MFLUX2'
      errflg = 1
      return
!     call MPI_CLOSE(1,routine)

!------------------------------------------------------------------------
!     update "zo" for ocean points.  "zo"cannot be updated within the
!     wegsteins iteration as the scheme (for the near neutral case)
!     can become unstable
!------------------------------------------------------------------------

130   continue
      do i = 1,ip
        szo = zoc(istb(i))
        if (szo .LT. 0.0)  then
        wndm=wind(istb(i))*0.01
          if(wndm.lt.15.0) then
          ckg=0.0185*og
          else
           ckg=(sfenth*(4*0.000308*wndm) + (1.-sfenth)*0.0185 )*og
          endif

        szo =  - ckg*wind(istb(i))*wind(istb(i))/  &
                            (sfm(i)*sfm(i))
        cons_p000001    =    .000001
        cons_7          =         7.
        vis             =     1.4E-1

        ustar    = sqrt( -szo / zog)
        restar = -ustar * szo    / vis
        restar = max(restar,cons_p000001)
!  Rat taken from Zeng,  Zhao and Dickinson 1997
        rat    = 2.67 * restar ** .25 - 2.57
        rat    = min(rat   ,cons_7)                      !constant
        rat=0.
        zot(istb(i)) = szo   * exp(-rat)
        else
         zot(istb(i)) = zoc(istb(i))
        endif

!      in hwrf thermal znot is loaded back into the zoc array for next step
            zoc(istb(i)) = szo
      enddo

      do i = 1,ip
        xxfm(istb(i)) = sfm(i)
        xxfh(istb(i)) = sfh(i)
        xxfh2(istb(i)) = sf2 (i)
        xxsh(istb(i)) = sfzo(i)
      enddo

!------------------------------------------------------------------------
!     obtain wind at 10 meters for diagnostic purposes
!------------------------------------------------------------------------

       do i = 1,ip
        wind10(istb(i)) = sf10(i)*uvs1(istb(i))/sfm(i)
        wind10(istb(i)) = wind10(istb(i)) * 1.944
          if(wind10(istb(i)) .GT. 6000.0) then
            wind10(istb(i))=wind10(istb(i))+wind10(istb(i))*cor1 &
                       - cor2
          endif
!     the above correction done by GFDL in centi-kts!!!-change back
        wind10(istb(i)) = wind10(istb(i)) / 1.944
       enddo

!------------------------------------------------------------------------
!     unstable points
!------------------------------------------------------------------------

170   continue

      iq = 0
      do i = its,ite
        if (zeta(i) .LT. 0.0) then
          iq       = iq + 1
          iutb(iq) = i
        endif
      enddo

      if (iq .EQ. 0) go to 290
      do i = 1,iq
        uzeta (i) = zeta(iutb(i))
        ifz   (i) = 1
        uzetam(i) = 1.0e+30
        ugzm  (i) = 0.0e+00
      enddo

!------------------------------------------------------------------------
!     begin wegstein iteration for "zeta" at unstable points using
!     hicks functions
!------------------------------------------------------------------------

      do icnt = 1,icntx
        do i = 1,iq
          if (ifz(i) .EQ. 0) go to 200
          ugzzo   = ALOG(zkmax(iutb(i))/ABS(zot(iutb(i))))
          uzetao  = ABS(zot(iutb(i)))/zkmax(iutb(i))*uzeta(i)
          ux11    = 1. - 16.*uzeta(i)
          ux12    = 1. - 16.*uzetao
          y       = SQRT(ux11)
          yo      = SQRT(ux12)
          ufzo(i) = 1./yo
          ux13    = (1. + y)/(1. + yo)
          ux21    = ALOG(ux13)
          ufh(i)  = (ugzzo - 2.*ux21)*rca
!          recompute scalers for ufm in terms of mom znot... zoc
          ugzzo   = ALOG(zkmax(iutb(i))/ABS(zoc(iutb(i))))
          uzetao  = ABS(zoc(iutb(i)))/zkmax(iutb(i))*uzeta(i)
          ux11    = 1. - 16.*uzeta(i)
          ux12    = 1. - 16.*uzetao
          y       = SQRT(ux11)
          yo      = SQRT(ux12)
          ux13    = (1. + y)/(1. + yo)
          ux21    = ALOG(ux13)
!           ufzo(i) = 1./yo
          x       = SQRT(y)
          xo      = SQRT(yo)
          xnum    = (x**2 + 1.)*((x + 1.)**2)
          xden    = (xo**2 + 1.)*((xo + 1.)**2)
          xtan    = ATAN(x) - ATAN(xo)
          ux3     = ALOG(xnum/xden)
          ufm(i)  = (ugzzo - ux3 + 2.*xtan)*rca

!------------------------------------------------------------------------
!     obtain ten meter winds for diagnostic purposes
!------------------------------------------------------------------------

          ugz10   = ALOG(z10(iutb(i))/ABS(zoc(iutb(i))))
          uzet1o  = ABS(z10(iutb(i)))/zkmax(iutb(i))*uzeta(i)
          uzetao  = ABS(zoc(iutb(i)))/zkmax(iutb(i))*uzeta(i)
          ux11    = 1. - 16.*uzet1o
          ux12    = 1. - 16.*uzetao
          y       = SQRT(ux11)
          y10     = SQRT(ux12)
          ux13    = (1. + y)/(1. + y10)
          ux21    = ALOG(ux13)
          x       = SQRT(y)
          x10     = SQRT(y10)
          xnum    = (x**2 + 1.)*((x + 1.)**2)
          xden    = (x10**2 + 1.)*((x10 + 1.)**2)
          xtan    = ATAN(x) - ATAN(x10)
          ux3     = ALOG(xnum/xden)
          uf10(i) = (ugz10 - ux3 + 2.*xtan)*rca

!   obtain 2m values for diagnostics...


          ugz2    = ALOG(z2   (iutb(i))/ABS(zoc(iutb(i))))
          uzet1o  = ABS(z2 (iutb(i)))/zkmax(iutb(i))*uzeta(i)
          uzetao  = ABS(zoc(iutb(i)))/zkmax(iutb(i))*uzeta(i)
          ux11    = 1. - 16.*uzet1o
          ux12    = 1. - 16.*uzetao
          y       = SQRT(ux11)
          yo      = SQRT(ux12)
          ux13    = (1. + y)/(1. + yo)
          ux21    = ALOG(ux13)
          uf2 (i)  = (ugzzo - 2.*ux21)*rca


          ugz = ca*rib(iutb(i))*ufm(i)*ufm(i)/(ufh(i) + enrca*ufzo(i))
          ux1 = (ugz - uzeta(i))/uzeta(i)
          ux2 = ABS(ux1)
          if (ux2 .GE. 5.0e-5) then
            uq        = (ugz - ugzm(i))/(uzeta(i) - uzetam(i))
            uzetam(i) = uzeta(i)
            if(uq .EQ. 1) then
             write(errmsg,'(*(a))') 'NCO ERROR DIVIDE BY ZERO IN gfdl_sfc_layer.F90/MFLUX2 (UNSTABLE CASE)'// &
                 'uq is 1 ',ux2,ugz,ugzm(i),uzeta(i),uzetam(i)
             errflg = 1
             return
            endif
            uzeta (i) = (ugz - uzeta(i)*uq)/(1.0 - uq)
            ugzm  (i) = ugz
          else
            ifz(i) = 0
          endif
200       continue
        enddo
      enddo


      do i = 1,iq
        if (ifz(i) .GE. 1) go to 230
      enddo

      go to 250

230   continue
      write(errmsg,'(*(a))') 'NON-CONVERGENCE FOR UNSTABLE ZETA IN ROW'// &
          'uq is 1 ',ux2,ugz,ugzm(i),uzeta(i),uzetam(i)
      errflg = 1
      return

!     call MPI_CLOSE(1,routine)

!------------------------------------------------------------------------
!     gather unstable values
!------------------------------------------------------------------------

250   continue

!------------------------------------------------------------------------
!     update "zo" for ocean points.  zo cannot be updated within the
!     wegsteins iteration as the scheme (for the near neutral case)
!     can become unstable.
!------------------------------------------------------------------------

      do i = 1,iq
        uzo = zoc(iutb(i))
        if (zoc(iutb(i)) .LT. 0.0)   then
        wndm=wind(iutb(i))*0.01
         if(wndm.lt.15.0) then
         ckg=0.0185*og
         else
            ckg=(4*0.000308*wndm)*og
            ckg=(sfenth*(4*0.000308*wndm) + (1.-sfenth)*0.0185 )*og
         endif
        uzo         =-ckg*wind(iutb(i))*wind(iutb(i))/(ufm(i)*ufm(i))
        cons_p000001    =    .000001
        cons_7          =         7.
        vis             =     1.4E-1

        ustar    = sqrt( -uzo / zog)
        restar = -ustar * uzo    / vis
        restar = max(restar,cons_p000001)
!   Rat taken from Zeng,  Zhao and Dickinson 1997
        rat    = 2.67 * restar ** .25 - 2.57
        rat    = min(rat   ,cons_7)                      !constant
        rat=0.0
        zot(iutb(i)) =  uzo   * exp(-rat)
        else
        zot(iutb(i)) = zoc(iutb(i))
       endif
!      in hwrf thermal znot is loaded back into the zoc array for next step
           zoc(iutb(i)) = uzo
      enddo

!------------------------------------------------------------------------
!     obtain wind at ten meters for diagnostic purposes
!------------------------------------------------------------------------
       do i = 1,iq
         wind10(iutb(i)) = uf10(i)*uvs1(iutb(i))/ufm(i)
         wind10(iutb(i)) = wind10(iutb(i)) * 1.944
          if(wind10(iutb(i)) .GT. 6000.0) then
        wind10(iutb(i))=wind10(iutb(i))+wind10(iutb(i))*cor1 &
                        - cor2
          endif
!     the above correction done by GFDL in centi-kts!!!-change back
         wind10(iutb(i)) = wind10(iutb(i)) / 1.944
       enddo

      do i = 1,iq
        xxfm(iutb(i)) = ufm(i)
        xxfh(iutb(i)) = ufh(i)
        xxfh2(iutb(i)) = uf2 (i)
        xxsh(iutb(i)) = ufzo(i)
      enddo

290   continue

      do i = its,ite
        ucom(i) = ukmax(i)
        vcom(i) = vkmax(i)
        if (windp(i) .EQ. 0.0) then
          windp(i) = 100.0
          ucom (i) = 100.0/SQRT(2.0)
          vcom (i) = 100.0/SQRT(2.0)
        endif
        rho(i) = pss(i)/(rgas*(tsg(i) + enrca*(theta(i) - &
                tsg(i))*xxsh(i)/(xxfh(i) + enrca*xxsh(i))))
        bq1(i) = wind(i)*rho(i)/(xxfm(i)*(xxfh(i) + enrca*xxsh(i)))
      enddo

!     do land sfc temperature prediction if ntsflg=1
!     ntsflg = 1                                    ! gopal's doing

      if (ntsflg .EQ. 0) go to 370
      alll = 600.
      xks   = 0.01
      hcap  = .5/2.39e-8
      pith  = SQRT(4.*ATAN(1.0))
      alfus = alll/2.39e-8
      teps  = 0.1
!     slwdc... in units of cal/min ????
!     slwa...  in units of ergs/sec/cm*2
!     1 erg=2.39e-8 cal
!------------------------------------------------------------------------
!     pack land and sea ice points
!------------------------------------------------------------------------

      ip    = 0
      do i = its,ite
        if (land(i) .EQ. 1) then
          ip = ip + 1
          indx   (ip) = i
!         slwa is defined as positive down....
          slwa   (ip) =    slwdc(i)/(2.39e-8*60.)
          tss    (ip) = tstar(i)
          thetap (ip) = theta(i)
          rkmaxp (ip) = rkmax(i)
          aap    (ip) = 5.673e-5
          pssp   (ip) = pss(i)
          ecofp  (ip) = ecof(i)
          estsop (ip) = estso(i)
          rstsop (ip) = rstso(i)
          bq1p   (ip) = bq1(i)
          bq1p   (ip) = amax1(bq1p(ip),0.1e-3)
          delsrad(ip) = dt   *pith/(hcap*SQRT(3600.*24.*xks))
        endif
      enddo

!------------------------------------------------------------------------
!     initialize variables for first pass of iteration
!------------------------------------------------------------------------

      do i = 1,ip
        ifz  (i) = 1
        tsm  (i) = tss(i)
        rdiff(i) = amin1(0.0,(rkmaxp(i) - rstsop(i)))

300   format(2X, ' SURFACE EQUILIBRIUM CALCULATION ')

          foftm(i) = tss(i) + delsrad(i)*(slwa(i) - aap(i)*tsm(i)**4 - &
           cp*bq1p(i)*(tsm(i) - thetap(i)) + ecofp(i)*alfus*bq1p(i)* &
           rdiff(i))
        tsp(i) = foftm(i)
      enddo

!------------------------------------------------------------------------
!     do iteration to determine "tstar" at new time level
!------------------------------------------------------------------------

      do icnt = 1,icntx
        do i = 1,ip
          if (ifz(i) .EQ. 0) go to 330
          tab1  (i) = tsp(i) - 153.16
          it    (i) = IFIX(tab1(i))
          tab2  (i) = tab1(i) - FLOAT(it(i))
          t1    (i) = tab(min(223,max(1,it(i) + 1)))
          t2    (i) = table(min(223,max(1,it(i) + 1)))
          estsop(i) = t1(i) + tab2(i)*t2(i)
             psps2 = (pssp(i) - estsop(i))
             if(psps2 .EQ. 0.0)then
               psps2 = .1
             endif
          rstsop(i) = ep2*estsop(i)/psps2
          rdiff (i) = amin1(0.0,(rkmaxp(i) - rstsop(i)))

            foft(i) = tss(i) + delsrad(i)*(slwa(i) - aap(i)*tsp(i)**4 - &
             cp*bq1p(i)*(tsp(i) - thetap(i)) + ecofp(i)*alfus*bq1p(i)* &
             rdiff(i))

          frac(i) = ABS((foft(i) - tsp(i))/tsp(i))

!------------------------------------------------------------------------
!      check for convergence of all points use wegstein iteration
!------------------------------------------------------------------------

         if (frac(i) .GE. teps) then
           qf   (i) = (foft(i) - foftm(i))/(tsp(i) - tsm(i))
           tsm  (i) = tsp(i)
           tsp  (i) = (foft(i) - tsp(i)*qf(i))/(1. - qf(i))
           foftm(i) = foft(i)
         else
           ifz(i) = 0
         endif
330       continue
        enddo
      enddo

!------------------------------------------------------------------------
!     check for convergence of "t star" prediction
!------------------------------------------------------------------------

      do i = 1,ip
        if (ifz(i) .EQ. 1) then
          write(errmsg,'(*(a))') 'NON-CONVERGENCE OF T* PREDICTED (T*,I) = ', &
              tsp(i), i
          errflg = 1
          return
!          call MPI_CLOSE(1,routine)
        endif
      enddo

      do i = 1,ip
        ii        = indx(i)
        tstrc(ii) = tsp (i)
      enddo

!------------------------------------------------------------------------
!     compute fluxes  and momentum drag coef
!------------------------------------------------------------------------

370   continue
      do i = its,ite
!!!
        if ( iwavecpl .eq. 1 .and. zoc(i) .le. 0.0 ) then
          windmks = wind10(i) * 0.01
          call znot_wind10m(windmks,znott,znotm,icoef_sf)
          !Check if Charnock parameter ratio is received in a proper range.
          if ( alpha(i) .ge. 0.2 .and. alpha(i) .le. 5. ) then
            znotm = znotm*alpha(i)
          endif
          zoc(i)  = -100.*znotm
          zot(i) =  -100* znott
        endif
!!!!
        fxh(i) = bq1(i)*(theta(i) - tsg(i))
        fxe(i) = ecof(i)*bq1(i)*(rkmax(i) - rstso(i))
        if (fxe(i) .GT. 0.0) fxe(i) = 0.0
        fxmx(i) = rho(i)/(xxfm(i)*xxfm(i))*wind(i)*wind(i)*ucom(i)/ &
                 windp(i)
        fxmy(i) = rho(i)/(xxfm(i)*xxfm(i))*wind(i)*wind(i)*vcom(i)/ &
        windp(i)
        cdm(i) = 1./(xxfm(i)*xxfm(i))
#if HWRF==1
! randomly perturb the Cd
!zzz        if( pert_Cd_local .and. ens_random_seed_local .gt. 0 ) then
       if( pert_Cd_local ) then
       ens_random_seed_local=ran1(-ens_random_seed_local)*1000
       rr=2.0*ens_Cdamp_local*ran1(-ens_random_seed_local)-ens_Cdamp_local
       cdm(i) = cdm(i) *(1.0+rr)
       endif
#endif

      enddo
      ntstep = ntstep + 1
      return
      end subroutine MFLUX2

      end module gfdl_sfc_layer

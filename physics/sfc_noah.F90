!>  \file sfc_noah.F90
!!  This file contains the Noah land surface scheme driver.

!> This module contains the CCPP-compliant Noah land surface scheme driver for 
!! the hurricane application.
      module sfc_noah

      implicit none

      private

      public :: sfc_noah_init, sfc_noah_run, sfc_noah_finalize

      contains

!> \ingroup Noah_LSM_hafs
!! \section arg_table_sfc_noah_init Argument Table
!! \htmlinclude sfc_noah_init.html
!!
      subroutine sfc_noah_init(lsm, lsm_noah_hafs, nsoil, ua_phys, fasdas, errmsg, errflg)

      use machine, only : kind_phys
      
      implicit none
      
      integer,              intent(in)  :: lsm, lsm_noah_hafs, nsoil, fasdas
      logical,              intent(in)  :: ua_phys

      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      if (lsm/=lsm_noah_hafs) then
        write(errmsg,'(*(a))') "Logic error: namelist choice of LSM is different from NOAH HAFS"
        errflg = 1
        return
      end if
      
      if (nsoil < 2) then
        write(errmsg,'(*(a))') "The NOAH HAFS scheme expects at least 2 soil layers."
        errflg = 1
        return
      end if
      
      if (ua_phys) then
        write(errmsg,'(*(a))') "The NOAH HAFS scheme has not been tested with ua_phys = T"
        errflg = 1
        return
      end if
      
      
      if (fasdas > 0) then
        write(errmsg,'(*(a))') "The NOAH HAFS scheme has not been tested with fasdas > 0"
        errflg = 1
        return
      end if
      
      !GJF: check for rdlai != F?
      !GJF: check for usemonalb != T?
      
      end subroutine sfc_noah_init


!! \section arg_table_sfc_noah_finalize Argument Table
!! \htmlinclude sfc_noah_finalize.html
!!
      subroutine sfc_noah_finalize(errmsg, errflg)

      implicit none

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      end subroutine sfc_noah_finalize


!> \defgroup Noah_LSM_hafs Noah LSM Model for the hurricane application
!! \section arg_table_sfc_noah_run Argument Table
!! \htmlinclude sfc_noah_run.html
!!
!> \section general_noah_hafs_drv GFS sfc_drv General Algorithm
!>  @{
      subroutine sfc_noah_run (im, flag_lsm, srflag, isurban, rdlai,    &
        ua_phys, usemonalb, aoasis, fasdas, dt, zlvl,                   &
        nsoil, sthick, lwdn, soldn, solnet, sfcprs, prcp, sfctmp, q1k,  &
        th1, qs1, dqsdt2, vegtyp, soiltyp, slopetyp, shdfac, shmin,     &
        shmax, albbrd, snoalb, tbot, z0brd, z0k, emissi, embrd, cmc, t1,&
        stc, smc, swc, snowhk, sneqv, chk, cp, rd, sigma, cph2o, cpice, &
        lsubf, sheat, eta_kinematic, ec, edir, ett, esnow, etp, ssoil,  &
        flx1, flx2, flx3, sncovr, runoff1, runoff2, soilm, qsurf, ribb, &
        smcwlt, smcref, smcmax, opt_thcnd, errmsg, errflg)
        
      use machine , only : kind_phys
      use module_sf_noahlsm, only: sflx, lutype, sltype

      implicit none
      
      integer,                             intent(in) :: im, isurban, nsoil, opt_thcnd, fasdas
      logical,                             intent(in) :: rdlai, ua_phys, usemonalb
      !GJF: usemonalb = True if the surface diffused shortwave albedo is EITHER read from input OR
      !  provided by a previous scheme (like radiation: as is done in GFS_rrtmgp_sw_pre)
      real(kind=kind_phys),                intent(in) :: aoasis
      
      real(kind=kind_phys),                intent(in) :: dt, cp, rd, sigma, cph2o, cpice, lsubf

      integer, dimension(im),              intent(in) :: vegtyp, soiltyp, slopetyp
      logical, dimension(im),              intent(in) :: flag_lsm
      real(kind=kind_phys), dimension(im), intent(in) :: srflag, zlvl, lwdn, soldn, solnet, &
                                                         sfcprs, prcp, sfctmp, q1k, th1, qs1, &
                                                         dqsdt2, shmin, shmax, snoalb, tbot
      real(kind=kind_phys), dimension(nsoil), intent(in) :: sthick

      real(kind=kind_phys), dimension(im), intent(inout) :: shdfac, albbrd, z0brd, z0k, emissi, &
                                                            cmc, t1, snowhk, sneqv, chk, flx1, &
                                                            flx2, flx3, ribb
      real(kind=kind_phys), dimension(im,nsoil), intent(inout) :: stc, smc, swc
      
      !variables that are intent(out) in module_sf_noahlsm, but are inout here due to being set within an IF statement
      real(kind=kind_phys), dimension(im), intent(inout) :: embrd, sheat, eta_kinematic, ec, &
                                                            edir, ett, esnow, etp, ssoil, sncovr, &
                                                            runoff1, runoff2, soilm, qsurf, smcwlt, &
                                                            smcref, smcmax

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      !GJF: There is some confusion regarding specific humidities vs mixing ratios in NOAH LSM.
      ! Looking at module_sf_noahlsm.F, sometimes the comments say mixing ratio and sometimes
      ! specific humidity. The WRF code (module_sf_noahdrv.F) specifically converts from mixing 
      ! ratio to specific humidity in preparation for calling SFLX, so I am assuming that
      ! all inputs/outputs into SFLX should be specific humidities, despite some comments in 
      ! module_sf_noahdrv.F describing arguments saying "mixing ratios". This applies to many
      ! arguments into SFLX (q1k, qs1, dqsdt2, eta_kinematic, qsurf, etc.).
      
!     local Variables
      integer :: i, k
      logical, parameter :: local = .false.  !(not actually used in SFLX) described in module_sf_noahlsm as:
      ! Flag for local-site simulation (where there is no maps for albedo, veg fraction, and roughness
      !   true:  all LSM parameters (inluding albedo, veg fraction and roughness length) will be defined by three tables

      real(kind=kind_phys) :: dummy 
      
      !GJF: The following variables are part of the interface to SFLX but not required as diagnostic
      ! output or otherwise outside of this subroutine (at least as part of a GFS-based suite). 
      ! If any of these variables are needed by other schemes or diagnostics, one needs to add it to 
      ! the host model and CCPP metadata. Alternatively, none of these variables NEED to be allocated
      ! and one could also just pass in dummy arguments. 
      !
      ! The variables descriptions are from module_sf_noahlsm.F:
      !
      ! albedok (output from SFLX): surface albedo including snow effect (unitless fraction)
      !                =snow-free albedo (alb) when sneqv=0, or
      !                =fct(msnoalb,alb,vegtyp,shdfac,shdmin) when sneqv>0
      ! eta (output from SFLX), eta_kinematic is what is passed out instead of eta
      ! fdown (output from SFLX) : Radiation forcing at the surface (W m-2) = SOLDN*(1-alb)+LWDN
      ! et (output from SFLX): plant transpiration from a particular root (soil) layer (W m-2)
      ! drip (output from SFLX): through-fall of precip and/or dew in excess of canopy water-holding capacity (m)
      ! dew (output from SFLX): dewfall (or frostfall for t<273.15) (m)
      ! beta (output from SFLX): ratio of actual/potential evap (dimensionless)
      ! snomlt (output from SFLX): snow melt (m) (water equivalent)
      ! runoff3 (output from SFLX): numerical trunctation in excess of porosity (smcmax) for a given soil layer at the end of a time step (m s-1).
      ! rc (output from SFLX): canopy resistance (s m-1)
      ! pc (output from SFLX): plant coefficient (unitless fraction, 0-1) where pc*etp = actual transp
      ! rsmin (output from SFLX): minimum canopy resistance (s m-1)
      ! xlai (output from SFLX): leaf area index (dimensionless)
      ! rcs (output from SFLX): incoming solar rc factor (dimensionless)
      ! rct (output from SFLX): air temperature rc factor (dimensionless)
      ! rcq (output from SFLX): atmos vapor pressure deficit rc factor (dimensionless)
      ! rcsoil (output from SFLX): soil moisture rc factor (dimensionless)
      ! soilw (output from SFLX): available soil moisture in root zone (unitless fraction between smcwlt and smcmax)
      ! smav (output from SFLX): soil moisture availability for each layer, as a fraction between smcwlt and smcmax.
      ! snotime1 (input/output from SFLX): no documentation in module_sf_noahlsm.F, but described as "initial number of timesteps since last snowfall" in module_sf_noahdrv.F; related to CCPP nondimensional_snow_age for NoahMP? Since inout, need to initialize here?
      ! smcdry (output from SFLX): dry soil moisture threshold where direct evap frm top layer ends (volumetric)
      ! smcmax (output from SFLX): porosity, i.e. saturated value of soil moisture (volumetric)
      ! nroot (output from SFLX): number of root layers, a function of veg type, determined in subroutine redprm.
      
      integer :: nroot
      real(kind=kind_phys) :: albedok, eta, fdown, drip, dew, beta, snomlt, &
                              runoff3, rc, pc, rsmin, xlai, rcs, rct, rcq,  &
                              rcsoil, soilw, snotime1, smcdry
      real (kind=kind_phys), dimension(nsoil) :: et, smav
      real(kind=kind_phys) :: sfcheadrt, infxsrt, etpnd1 !don't appear to be used unless WRF_HYDRO preprocessor directive is defined and no documentation
      real(kind=kind_phys) :: xsda_qfx, hfx_phy, qfx_phy, xqnorm, hcpct_fasdas !only used if fasdas = 1
      
      !variables associated with UA_PHYS (not used for now)
      real(kind=kind_phys) :: flx4, fvb, fbur, fgsn

      errmsg = ''
      errflg = 0
      
      do i=1, im
        if (flag_lsm(i)) then  
          !GJF: Why do LSMs want the dynamics time step instead of the physics time step?
          call sflx (i, 1, srflag(i), &
                   isurban, dt, zlvl(i), nsoil, sthick,              &    !c
                   local,                                            &    !L
                   lutype, sltype,                                   &    !CL
                   lwdn(i), soldn(i), solnet(i), sfcprs(i), prcp(i), &    !F
                   sfctmp(i), q1k(i), dummy, dummy, dummy, dummy,    &    !F
                   th1(i), qs1(i), dqsdt2(i),                        &    !I
                   vegtyp(i), soiltyp(i), slopetyp(i), shdfac(i),    &    !I
                   shmin(i), shmax(i),                               &    !I
                   albbrd(i), snoalb(i), tbot(i), z0brd(i), z0k(i),  &    !S
                   emissi(i), embrd(i),                              &    !S
                   cmc(i), t1(i), stc(i,:), smc(i,:), swc(i,:),      &    !H
                   snowhk(i), sneqv(i), albedok, chk(i), dummy,      &    !H
                   cp, rd, sigma, cph2o, cpice, lsubf,               &
                   eta, sheat(i), eta_kinematic(i), fdown,           &    !O
                   ec(i), edir(i), et, ett(i), esnow(i), drip, dew,  &    !O
                   beta, etp(i), ssoil(i), flx1(i), flx2(i), flx3(i),&    !O
                   flx4, fvb, fbur, fgsn, ua_phys,                   &    !UA 
                   snomlt, sncovr(i), runoff1(i), runoff2(i),runoff3,&    !O
                   rc, pc, rsmin, xlai, rcs, rct, rcq, rcsoil,       &    !O
                   soilw, soilm(i), qsurf(i), smav,                  &    !D
                   rdlai, usemonalb, snotime1, ribb(i),              &
                   smcwlt(i), smcdry, smcref(i), smcmax(i), nroot,   &
                   sfcheadrt, infxsrt, etpnd1, opt_thcnd, aoasis,    &
                   xsda_qfx, hfx_phy, qfx_phy, xqnorm, fasdas,       &    !fasdas
                   hcpct_fasdas,                                     &    !fasdas
                   errflg, errmsg)
          if (errflg > 0) return
        endif
      end do
      
      
      end subroutine sfc_noah_run
!> @}

end module sfc_noah

!>  \file ugwpv1_gsldrag.F90
!! This introduces two gravity wave drag schemes ugwpv1/CIRES and GSL/drag_suite.F90 under "ugwpv1_gsldrag" suite:
!!      1) The "V1 CIRES UGWP" scheme as tested in the FV3GFSv16-127L atmosphere model and workflow, which includes:
!!            a) the orograhic gravity wave drag, flow blocking scheme and TOFD (Beljaars et al, 2004).
!!            b) the v1 CIRE ugwp non-stationary GW scheme, new revision that generate realistic climate of FV3GFS-127L
!!               in the strato-mesosphere in the multi-year simulations (Annual cycles, SAO and QBO in th tropical dynamics).
!!      2) The GSL orographic drag suite (drag_suite.F90), as implmeneted in the RAP/HRRR, which includes:
!!            a) large-scale gravity wave drag and low-level flow blocking -- active at horizontal scales
!!               down to ~5km (Kim and Arakawa, 1995 \cite kim_and_arakawa_1995; Kim and Doyle, 2005 \cite kim_and_doyle_2005)
!!            b) small-scale gravity wave drag scheme -- active typically in stable PBL at horizontal grid resolutions down to ~1km
!!               (Steeneveld et al, 2008 \cite steeneveld_et_al_2008; Tsiringakis et al, 2017 \cite tsiringakis_et_al_2017)
!!            c) turbulent orographic form drag -- active at horizontal grid ersolutions down to ~1km
!!               (Beljaars et al, 2004 \cite beljaars_et_al_2004)
!! See Valery Yudin's presentation at 2020 UFS User's meeting (Jul 2020):
!! Gravity waves (GWs): Mesoscale GWs transport momentum, energy (heat) , and create eddy mixing in the whole atmosphere domain; Breaking and dissipating GWs deposit: (a) momentum; (b) heat (energy); and create (c) turbulent mixing of momentum, heat, and tracers
!! To properly incorporate GW effects (a-c) unresolved by DYCOREs we need GW physics
!! "Unified": a) all GW effects due to both dissipation/breaking; b) identical GW solvers for all GW sources; c) ability to replace solvers.
!! Unified Formalism:
!! 1. GW Sources: Stochastic and physics based mechanisms for GW-excitations in the lower atmosphere, calibrated by the high-res analyses/forecasts, and observations (3 types of GW sources: orography, convection, fronts/jets).
!! 2. GW Propagation: Unified solver for "propagation, dissipation and breaking" excited from all type of GW sources.
!! 3. GW Effects: Unified representation of GW impacts on the "resolved" flow for all sources (energy-balanced schemes for momentum, heat and mixing).
!! https://www.weather.gov/media/sti/nggps/Presentations%202017/02%20NGGPS_VYUDIN_2017_.pdf
!!
!! The ugwpv1_gsldrag scheme is activated by gwd_opt = 2 in the namelist.
!! The choice of schemes is activated at runtime by the following namelist options (boolean):
!! NA    do_ugwp_v0           -- activates V0 CIRES UGWP scheme - both orographic and non-stationary GWD is not active (NA)
!! NA    do_ugwp_v0_orog_only -- activates V0 CIRES UGWP scheme - orographic GWD only
!!       do_gsl_drag_ls_bl    -- activates RAP/HRRR (GSL) large-scale OGWD and blocking
!!       do_gsl_drag_ss       -- activates RAP/HRRR (GSL) small-scale OGWD
!!       do_gsl_drag_tofd     -- activates RAP/HRRR (GSL) turbulent orographic drag
!!       do_ugwp_v1           -- activates V1 CIRES UGWP scheme - both orographic and non-stationary GWD
!!       do_ugwp_v1_orog_only -- activates V1 CIRES UGWP scheme - orographic GWD only
!!       do_ugwp_v1_w_gsldrag -- activates V1 CIRES UGWP scheme with orographic drag of GSL
!! Note that only one "large-scale" scheme can be activated at a time.
!!

module ugwpv1_gsldrag

    use machine, only: kind_phys

    use cires_ugwpv1_triggers, only:  slat_geos5_2020, slat_geos5_tamp_v1
    use cires_ugwpv1_module,   only:  cires_ugwpv1_init, ngwflux_update, calendar_ugwp
    use cires_ugwpv1_module,   only:  knob_ugwp_version, cires_ugwp_dealloc, tamp_mpa
    use cires_ugwpv1_solv2,    only:  cires_ugwpv1_ngw_solv2
    use cires_ugwpv1_oro,      only:  orogw_v1

    use drag_suite,            only:  drag_suite_run

    implicit none

    private

    public ugwpv1_gsldrag_init, ugwpv1_gsldrag_run, ugwpv1_gsldrag_finalize

    logical :: is_initialized = .False.

contains

! ------------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme v0
! ------------------------------------------------------------------------
!>@brief The subroutine initializes the unified UGWP
!> \section arg_table_ugwpv1_gsldrag_init Argument Table
!! \htmlinclude ugwpv1_gsldrag_init.html
!!
! -----------------------------------------------------------------------
!
    subroutine ugwpv1_gsldrag_init  (                                          &
                me, master, nlunit, input_nml_file, logunit,                   &
                fn_nml2, jdat, lonr, latr, levs, ak, bk, dtp,                  &
                con_pi, con_rerth, con_p0,                                     &
                con_g, con_omega,  con_cp, con_rd, con_rv,con_fvirt,           &
                do_ugwp,do_ugwp_v0, do_ugwp_v0_orog_only, do_gsl_drag_ls_bl,   &
                do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1,                  &
                do_ugwp_v1_orog_only, do_ugwp_v1_w_gsldrag, errmsg, errflg)

    use ugwp_common

!----  initialization of unified_ugwp
    implicit none

    integer,              intent (in) :: me
    integer,              intent (in) :: master
    integer,              intent (in) :: nlunit
    character(len=*),     intent (in) :: input_nml_file(:)
    integer,              intent (in) :: logunit
    integer,              intent(in)  :: jdat(8)
    integer,              intent (in) :: lonr
    integer,              intent (in) :: levs
    integer,              intent (in) :: latr
    real(kind=kind_phys), intent (in) :: ak(levs+1), bk(levs+1)
    real(kind=kind_phys), intent (in) :: dtp

    real(kind=kind_phys), intent (in) :: con_p0, con_pi, con_rerth
    real(kind=kind_phys), intent(in)  :: con_g,  con_cp, con_rd, con_rv, con_omega, con_fvirt
    logical,              intent (in) :: do_ugwp

    logical,              intent (in) :: do_ugwp_v0, do_ugwp_v0_orog_only,  &
                                         do_gsl_drag_ls_bl, do_gsl_drag_ss, &
                                         do_gsl_drag_tofd, do_ugwp_v1,      &
                                         do_ugwp_v1_orog_only,do_ugwp_v1_w_gsldrag

    character(len=*), intent (in)  :: fn_nml2
    !character(len=*), parameter   :: fn_nml='input.nml'

    integer :: ios
    logical :: exists
    real    :: dxsg
    integer :: k

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
!============================================================================
!
!     gwd_opt => "1 and 2, 3, 22, 33' see previous GSL-commits
!                 related to GSL-oro drag suite
!    for use of the new-GSL/old-GFS/EMC inputs for sub-grid orography
!  see details inside /ufs-weather-model/FV3/io/FV3GFS_io.F90
!  FV3GFS_io.F90:    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
!  FV3GFS_io.F90:        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then
!  FV3GFS_io.F90:          if ( (Model%gwd_opt==3 .or. Model%gwd_opt==33) .or.    &
!  FV3GFS_io.F90:               ( (Model%gwd_opt==2 .or. Model%gwd_opt==22) .and. &
!
! gwd_opt=1 -current 14-element GFS-EMC subgrid-oro input
! gwd_opt=2 and 3    24-element -current 14-element GFS-EMC subgrid-oro input
! GSL uses the gwd_opt flag to control "extra" diagnostics (22 and 33)
! CCPP may use gwd_opt to determine 14 or 24 variables for the input
!      but at present you work with "nmtvr"
! GFS_GWD_generic.F90:      integer, intent(in) :: im, levs, nmtvr
!GFS_GWD_generic.F90:      real(kind=kind_phys), intent(in) :: mntvar(im,nmtvr)
!GFS_GWD_generic.F90:      if (nmtvr == 14) then  ! gwd_opt=1 current operational - as of 2014
!GFS_GWD_generic.F90:      elseif (nmtvr == 10) then ????
!GFS_GWD_generic.F90:      elseif (nmtvr == 6) then  ????
!GFS_GWD_generic.F90:      elseif (nmtvr == 24) then   ! GSD_drag_suite and unified_ugwp gwd_opt=2,3
!
! 1) gsldrag:   do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1
! 2) CIRES-v1:  do_ugwp_v1,        do_ugwp_v1_orog_only,  do_tofd,   ldiag_ugwp
!==============================================================================
    ! Test to make sure that at most only one large-scale/blocking
    ! orographic drag scheme is chosen
    if ( (do_ugwp_v0.and.(do_ugwp_v0_orog_only.or.do_gsl_drag_ls_bl.or.    &
                          do_ugwp_v1.or.do_ugwp_v1_orog_only))        .or. &
         (do_ugwp_v0_orog_only.and.(do_gsl_drag_ls_bl.or.do_ugwp_v1.or.    &
                                    do_ugwp_v1_orog_only))            .or. &
         (do_gsl_drag_ls_bl.and.do_ugwp_v1_orog_only)  ) then

       write(errmsg,'(*(a))') "Logic error: Only one large-scale&
          &/blocking scheme (do_ugwp_v0,do_ugwp_v0_orog_only,&
          &do_gsl_drag_ls_bl,do_ugwp_v1 or &
          &do_ugwp_v1_orog_only) can be chosen"
       errflg = 1
       return

    end if
!
    if ( do_ugwp_v0_orog_only .or. do_ugwp_v0) then
       print *,  ' ccpp do_ugwp_v0 active ', do_ugwp_v0
       print *,  ' ccpp do_ugwp_v1_orog_only active ', do_ugwp_v0_orog_only
       write(errmsg,'(*(a))') " the CIRES <ugwpv1_gsldrag> CCPP-suite does not &
         support <ugwp_v0> schemes "
       errflg = 1
       return
    endif
!
    if (do_ugwp_v1_w_gsldrag .and. do_ugwp_v1_orog_only ) then

       print *,  '  do_ugwp_v1_w_gsldrag ', do_ugwp_v1_w_gsldrag
       print *,  '  do_ugwp_v1_orog_only ', do_ugwp_v1_orog_only
       print *,  '  do_gsl_drag_ls_bl ',do_gsl_drag_ls_bl
       write(errmsg,'(*(a))') " the CIRES <ugwpv1_gsldrag> CCPP-suite intend to &
         support <ugwp_v1> with <gsldrag>  but  has Logic error"
       errflg = 1
       return
    endif
!==========================
!
! initialize ugwp_common
!   con_pi, con_rerth, con_p0,  con_g, con_omega,  con_cp, con_rd, con_rv,con_fvirt
!
!==========================

    pi    = con_pi
    arad  = con_rerth
    p0s   = con_p0
    grav  = con_g
    omega1= con_omega
    cpd   = con_cp
    rd    = con_rd
    rv    = con_rv
    fv    = con_fvirt

    grav2  = grav + grav; rgrav  = 1.0/grav ; rgrav2 = rgrav*rgrav
    rdi    = 1.0 / rd ; rcpd = 1./cpd;  rcpd2 = 0.5/cpd
    gor    = grav/rd
    gr2    = grav*gor
    grcp   = grav*rcpd
    gocp   = grcp
    rcpdl  = cpd*rgrav
    grav2cpd = grav*grcp

    pi2      = 2.*pi ;  pih = .5*pi
    rad_to_deg=180.0/pi
    deg_to_rad=pi/180.0

    bnv2min = (pi2/1800.)*(pi2/1800.)
    bnv2max = (pi2/30.)*(pi2/30.)
    dw2min  = 1.0
    velmin  = sqrt(dw2min)
    minvel  = 0.5

    omega2  = 2.*omega1
    omega3  = 3.*omega1

    hpscale = 7000. ; hpskm = hpscale*1.e-3
    rhp     = 1./hpscale
    rhp2 = 0.5*rhp; rh4 = 0.25*rhp
    rhp4 = rhp2 * rhp2
    khp  = rhp* rd/cpd
    mkzmin  = pi2/80.0e3
    mkz2min = mkzmin*mkzmin
    mkzmax  = pi2/500.
    mkz2max = mkzmax*mkzmax
    cdmin   = 2.e-2/mkzmax

    rcpdt  = rcpd/dtp

    if ( do_ugwp_v1 ) then
       call cires_ugwpv1_init (me, master, nlunit, logunit, jdat, con_pi,      &
                               con_rerth, fn_nml2, input_nml_file, lonr, latr, &
                               levs, ak, bk, con_p0, dtp, errmsg, errflg)
       if (errflg/=0) return
    end if

    if (me == master) then
       print *,  ' ccpp: ugwpv1_gsldrag_init   '

       print *,  ' ccpp do_ugwp_v1  flag ', do_ugwp_v1
       print *,  ' ccpp do_gsl_drag_ls_bl  flag ',    do_gsl_drag_ls_bl
       print *,  ' ccpp do_gsl_drag_ss  flag ' ,      do_gsl_drag_ss
       print *,  ' ccpp do_gsl_drag_tofd  flag ',     do_gsl_drag_tofd

       print *, ' ccpp: ugwpv1_gsldrag_init  '
    endif



    is_initialized = .true.


    end subroutine ugwpv1_gsldrag_init


! -----------------------------------------------------------------------
! finalize of ugwpv1_gsldrag   (_finalize)
! -----------------------------------------------------------------------

!>@brief The subroutine finalizes the CIRES UGWP

!> \section arg_table_ugwpv1_gsldrag_finalize Argument Table
!! \htmlinclude ugwpv1_gsldrag_finalize.html
!!

    subroutine ugwpv1_gsldrag_finalize(errmsg, errflg)

    implicit none
!
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

    call cires_ugwp_dealloc

    is_initialized = .false.

    end subroutine ugwpv1_gsldrag_finalize


! -----------------------------------------------------------------------
!    originally from ugwp_driver_v0.f
!    driver of cires_ugwp   (_driver)
! -----------------------------------------------------------------------
!   driver is called after pbl & before chem-parameterizations
! -----------------------------------------------------------------------
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
! -----------------------------------------------------------------------
!>@brief These subroutines and modules execute the CIRES UGWP Version 0
!>\defgroup ugwpv1_gsldrag_run Unified Gravity Wave Physics General Algorithm
!> @{
!! The physics of NGWs in the UGWP framework (Yudin et al. 2018 \cite yudin_et_al_2018) is represented by four GW-solvers, which is introduced in Lindzen (1981) \cite lindzen_1981, Hines (1997) \cite hines_1997, Alexander and Dunkerton (1999) \cite alexander_and_dunkerton_1999, and Scinocca (2003) \cite scinocca_2003. The major modification of these GW solvers is represented by the addition of the background dissipation of temperature and winds to the saturation criteria for wave breaking. This feature is important in the mesosphere and thermosphere for WAM applications and it considers appropriate scale-dependent dissipation of waves near the model top lid providing the momentum and energy conservation in the vertical column physics (Shaw and Shepherd 2009 \cite shaw_and_shepherd_2009). In the UGWP-v0, the modification of Scinocca (2003) \cite scinocca_2003 scheme for NGWs with non-hydrostatic and rotational effects for GW propagations and background dissipation is represented by the subroutine \ref fv3_ugwp_solv2_v0. In the next release of UGWP, additional GW-solvers will be implemented along with physics-based triggering of waves and stochastic approaches for selection of GW modes characterized by horizontal phase velocities, azimuthal directions and magnitude of the vertical momentum flux (VMF).
!!
!! In UGWP-v0, the specification for the VMF function is adopted from the GEOS-5 global atmosphere model of GMAO NASA/GSFC, as described in Molod et al. (2015) \cite molod_et_al_2015 and employed in the MERRRA-2 reanalysis (Gelaro et al., 2017 \cite gelaro_et_al_2017). The Fortran subroutine \ref slat_geos5_tamp describes the latitudinal shape of VMF-function as displayed in Figure 3 of Molod et al. (2015) \cite molod_et_al_2015. It shows that the enhanced values of VMF in the equatorial region gives opportunity to simulate the QBO-like oscillations in the equatorial zonal winds and lead to more realistic simulations of the equatorial dynamics in GEOS-5 operational and MERRA-2 reanalysis products. For the first vertically extended version of FV3GFS in the stratosphere and mesosphere, this simplified function of VMF allows us to tune the model climate and to evaluate multi-year simulations of FV3GFS with the MERRA-2 and ERA-5 reanalysis products, along with temperature, ozone, and water vapor observations of current satellite missions. After delivery of the UGWP-code, the EMC group developed and tested approach to modulate the zonal mean NGW forcing by 3D-distributions of the total precipitation as a proxy for the excitation of NGWs by convection and the vertically-integrated  (surface - tropopause) Turbulent Kinetic Energy (TKE). The verification scores with updated NGW forcing, as reported elsewhere by EMC researchers, display noticeable improvements in the forecast scores produced by FV3GFS configuration extended into the mesosphere.
!!
!> \section arg_table_ugwpv1_gsldrag_run Argument Table
!! \htmlinclude ugwpv1_gsldrag_run.html
!!
!> \section gen_ugwpv1_gsldrag CIRES UGWP Scheme General Algorithm
!! @{
     subroutine ugwpv1_gsldrag_run(me, master, im,  levs, ntrac, lonr, dtp, fhzero,kdt, &
          ldiag3d, lssav, flag_for_gwd_generic_tend, do_gsl_drag_ls_bl, do_gsl_drag_ss, &
          do_gsl_drag_tofd, do_ugwp_v1, do_ugwp_v1_orog_only, do_ugwp_v1_w_gsldrag,     &
          gwd_opt, do_tofd, ldiag_ugwp, cdmbgwd, jdat,                                  &
!          con_g, con_omega, con_pi, con_cp, con_rd, con_rv, con_rerth, con_fvirt,       &
          nmtvr, hprime, oc, theta, sigma, gamma, elvmax, clx, oa4,                     &
          varss,oc1ss,oa4ss,ol4ss, dx,  xlat, xlat_d, sinlat, coslat, area,             &
          rain, br1, hpbl, kpbl, slmsk,                                                 &
          ugrs, vgrs, tgrs, q1, prsi, prsl, prslk, phii, phil,  del, tau_amf,           &
          dudt_ogw, dvdt_ogw, du_ogwcol, dv_ogwcol,                                     &
          dudt_obl, dvdt_obl, du_oblcol, dv_oblcol,                                     &
          dudt_oss, dvdt_oss, du_osscol, dv_osscol,                                     &
          dudt_ofd, dvdt_ofd, du_ofdcol, dv_ofdcol,                                     &
          dudt_ngw, dvdt_ngw, dtdt_ngw, kdis_ngw, dudt_gw, dvdt_gw, dtdt_gw, kdis_gw,   &
      tau_ogw, tau_ngw,  tau_oss,                                                   &
          zogw,  zlwb,  zobl,  zngw,   dusfcg, dvsfcg,  dudt, dvdt, dtdt, rdxzb,        &
          ldu3dt_ogw, ldv3dt_ogw, ldt3dt_ogw, ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw,       &
      lprnt, ipr, errmsg, errflg)
!
!########################################################################
!  Attention New Arrays and Names must be ADDED inside
!
! a) /FV3/gfsphysics/GFS_layer/GFS_typedefs.meta
! b) /FV3/gfsphysics/GFS_layer/GFS_typedefs.F90
! c) /FV3/gfsphysics/GFS_layer/GFS_diagnostics.F90 "diag-cs is not tested"
!########################################################################

!

    use ugwp_common, only : con_pi => pi, con_g => grav,  con_rd   => rd,   &
                            con_rv => rv, con_cp => cpd,  con_fv   => fv,   &
                            con_rerth => arad, con_omega => omega1, rgrav

    implicit none

! Preference use    (im,levs) rather than (:,:) to avoid memory-leaks
!                    that found in Nov-Dec 2020
! order array-description control-logical
!                         other in-variables
!                         out-variables
!                         local-variables
!
!  unified GSL and CIRES diagnostics inside CCPP and GFS_typedefs.F90/GFS_diagnostics.F90
!
!
! interface variables
    logical,                 intent(in) :: ldiag3d, lssav
    logical,                 intent(in) :: flag_for_gwd_generic_tend
    logical,                 intent(in) :: lprnt

    integer,                 intent(in) :: ipr

! flags for choosing combination of GW drag schemes to run

    logical,  intent (in) :: do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd
    logical,  intent (in) :: do_ugwp_v1, do_ugwp_v1_orog_only, do_tofd, ldiag_ugwp
    logical,  intent (in) :: do_ugwp_v1_w_gsldrag                              ! combination of ORO and NGW schemes

    integer,                 intent(in) :: me, master, im, levs, ntrac,lonr
    real(kind=kind_phys),    intent(in) :: dtp, fhzero
    integer,                 intent(in) :: kdt, jdat(8)

! SSO parameters and variables
    integer,                 intent(in) :: gwd_opt                         !gwd_opt  and nmtvr are "redundant" controls
    integer,                 intent(in) :: nmtvr
    real(kind=kind_phys),    intent(in) :: cdmbgwd(4)                      ! for gsl_drag

    real(kind=kind_phys),    intent(in), dimension(im)       :: hprime, oc, theta, sigma, gamma

    real(kind=kind_phys),    intent(in), dimension(im)       :: elvmax
    real(kind=kind_phys),    intent(in), dimension(im, 4)    :: clx, oa4

    real(kind=kind_phys),    intent(in), dimension(im)       :: varss,oc1ss,dx
    real(kind=kind_phys),    intent(in), dimension(im, 4)    :: oa4ss,ol4ss

!=====
!ccpp-style passing constants, I prefer to take them out from the "call-subr" list
!=====
!    real(kind=kind_phys),    intent(in) :: con_g, con_omega, con_pi, con_cp, con_rd, &
!                                           con_rv, con_rerth, con_fvirt
! grids

    real(kind=kind_phys),    intent(in), dimension(im)         :: xlat, xlat_d, sinlat, coslat, area

! State vars + PBL/slmsk +rain

    real(kind=kind_phys),    intent(in), dimension(im, levs)   :: del, ugrs, vgrs, tgrs, prsl, prslk, phil
    real(kind=kind_phys),    intent(in), dimension(im, levs+1) :: prsi, phii
    real(kind=kind_phys),    intent(in), dimension(im, levs)   :: q1
    integer,                 intent(in), dimension(im)         :: kpbl

    real(kind=kind_phys),    intent(in), dimension(im) :: rain
    real(kind=kind_phys),    intent(in), dimension(im) :: br1, hpbl,  slmsk
!
! moved to GFS_phys_time_vary
!    real(kind=kind_phys),    intent(in), dimension(im) :: ddy_j1tau, ddy_j2tau
!    integer,                 intent(in), dimension(im) :: jindx1_tau, jindx2_tau
     real(kind=kind_phys),    intent(in), dimension(im) :: tau_amf

!Output (optional):

    real(kind=kind_phys), intent(out), dimension(im)  ::                  &
                            du_ogwcol,  dv_ogwcol,  du_oblcol, dv_oblcol, &
                            du_osscol,  dv_osscol,  du_ofdcol, dv_ofdcol
!
! we may add later but due to launch in the upper layes ~ mPa comparing to ORO Pa*(0.1)
!                            du_ngwcol, dv_ngwcol

    real(kind=kind_phys), intent(out), dimension(im)  :: dusfcg, dvsfcg
    real(kind=kind_phys), intent(out), dimension(im)  :: tau_ogw, tau_ngw, tau_oss

    real(kind=kind_phys), intent(out) , dimension(im, levs) ::    &
                          dudt_ogw, dvdt_ogw, dudt_obl, dvdt_obl, &
                          dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd

    real(kind=kind_phys), intent(out) , dimension(im, levs) :: dudt_ngw, dvdt_ngw, kdis_ngw
    real(kind=kind_phys), intent(out) , dimension(im, levs) :: dudt_gw,  dvdt_gw,  kdis_gw

    real(kind=kind_phys), intent(out) , dimension(im, levs) :: dtdt_ngw, dtdt_gw

    real(kind=kind_phys), intent(out) , dimension(im) ::  zogw,  zlwb,  zobl,  zngw
!
!
    real(kind=kind_phys), intent(inout), dimension(im, levs) :: dudt, dvdt, dtdt

!
! These arrays are only allocated if ldiag=.true.
!
! Version of COORDE updated by CCPP-dev for time-aver
!
    real(kind=kind_phys),    intent(inout), dimension(:,:)   :: ldu3dt_ogw, ldv3dt_ogw, ldt3dt_ogw
    real(kind=kind_phys),    intent(inout), dimension(:,:)   :: ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw



    real(kind=kind_phys),    intent(out), dimension(im)      :: rdxzb     ! for stoch phys. mtb-level

    character(len=*),        intent(out) :: errmsg
    integer,                 intent(out) :: errflg

! local variables
    integer :: i, k
    real(kind=kind_phys), dimension(im)       :: sgh30
    real(kind=kind_phys), dimension(im, levs) :: Pdvdt, Pdudt
    real(kind=kind_phys), dimension(im, levs) :: Pdtdt, Pkdis
!------------
!
! from ugwp_driver_v0.f -> cires_ugwp_initialize.F90 -> module ugwp_wmsdis_init
!  now in the namelist of cires_ugwp "knob_ugwp_tauamp" controls tamp_mpa
!
!        tamp_mpa =knob_ugwp_tauamp                         !amplitude for GEOS-5/MERRA-2
!------------
!    real(kind=kind_phys), parameter :: tamp_mpa_v0=30.e-3  ! large flux to help "GFS-ensembles" in July 2019

! switches that activate impact of OGWs and NGWs

!    integer :: nmtvr_temp

    real(kind=kind_phys), dimension(im, levs)   :: zmet  ! geopotential height at model Layer centers
    real(kind=kind_phys), dimension(im, levs+1) :: zmeti ! geopotential height at model layer interfaces


! ugwp_v1 local variables

    integer :: y4, month, day,  ddd_ugwp, curdate, curday

!  ugwp_v1 temporary (local) diagnostic variables from cires_ugwp_solv2_v1
!  diagnostics for wind and temp rms to compare with space-borne data and metrics
!   in the Middle atmosphere: 20-110 km ( not active in CCPP-style, oct 2020)
!    real(kind=kind_phys) :: tauabs(im,levs), wrms(im,levs), trms(im,levs)


    ! Initialize CCPP error handling variables

    errmsg = ''
    errflg = 0

! 1) ORO stationary GWs
!    ------------------
!
! for all oro-suites can uze geo-meters having "hpbl"
!
!
! All GW-schemes operate with Zmet =phil*inv_g, passing Zmet/Zmeti can be more robust
! + rho*dz = =delp *  inv_g   can be also pre-comp for all "GW-schemes"
!
       zmeti  = phii* rgrav
       zmet   = phil* rgrav

!===============================================================
! ORO-diag

      dudt_ogw(:,:)  = 0. ; dvdt_ogw(:,:)=0. ; dudt_obl(:,:)=0. ; dvdt_obl(:,:)=0.
      dudt_oss(:,:)  = 0. ; dvdt_oss(:,:)=0. ; dudt_ofd(:,:)=0. ; dvdt_ofd(:,:)=0.

      dusfcg (:)  = 0.  ;  dvsfcg(:) =0.

      du_ogwcol(:)=0. ; dv_ogwcol(:)=0. ; du_oblcol(:)=0. ; dv_oblcol(:)=0.
      du_osscol(:)=0. ; dv_osscol(:)=0. ;du_ofdcol(:)=0.  ; dv_ofdcol(:)=0.

!
       dudt_ngw(:,:)=0. ; dvdt_ngw(:,:)=0. ; dtdt_ngw(:,:)=0. ; kdis_ngw(:,:)=0.

! ngw+ogw - diag

       dudt_gw(:,:)=0. ;  dvdt_gw(:,:)=0.  ; dtdt_gw(:,:)=0.  ; kdis_gw(:,:)=0.
! source fluxes

      tau_ogw(:)=0. ; tau_ngw(:)=0. ;  tau_oss(:)=0.

! launch layers

      zlwb(:)= 0.  ; zogw(:)=0. ;  zobl(:)=0. ;  zngw(:)=0.
!===============================================================
!  diag tendencies due to all-SSO schemes (ORO-physics)
!  ogw + obl + oss + ofd ..... no explicit "lee wave trapping"
!===============================================================
     do k=1,levs
        do i=1,im
          Pdvdt(i,k) = 0.0
          Pdudt(i,k) = 0.0
          Pdtdt(i,k) = 0.0
          Pkdis(i,k) = 0.0
        enddo
      enddo
!
    ! Run the appropriate large-scale (large-scale GWD + blocking) scheme
    ! Note:  In case of GSL drag_suite, this includes ss and tofd

    if ( do_gsl_drag_ls_bl.or.do_gsl_drag_ss.or.do_gsl_drag_tofd) then
!
! to do: the zero diag and tendency values assigned inside "drag_suite_run" can be skipped :
!
! dudt_ogw, dvdt_ogw, dudt_obl, dvdt_obl,dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd
! du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol, du_osscol, dv_osscol, du_ofdcol dv_ofdcol
! dusfcg,  dvsfcg
!
!
       call drag_suite_run(im,levs, Pdvdt, Pdudt, Pdtdt,             &
                 ugrs,vgrs,tgrs,q1,                                  &
                 kpbl,prsi,del,prsl,prslk,phii,phil,dtp,             &
                 kdt,hprime,oc,oa4,clx,varss,oc1ss,oa4ss,            &
                 ol4ss,theta,sigma,gamma,elvmax,                     &
                  dudt_ogw, dvdt_ogw, dudt_obl, dvdt_obl,            &
                  dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd,            &
                  dusfcg,  dvsfcg,                                   &
                  du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol,        &
                  du_osscol, dv_osscol, du_ofdcol, dv_ofdcol,        &
                 slmsk,br1,hpbl, con_g,con_cp,con_rd,con_rv,         &
                 con_fv, con_pi, lonr,                               &
                 cdmbgwd(1:2),me,master,lprnt,ipr,rdxzb,dx,gwd_opt,  &
                 do_gsl_drag_ls_bl,do_gsl_drag_ss,do_gsl_drag_tofd,  &
                 errmsg,errflg)
!
! dusfcg = du_ogwcol + du_oblcol + du_osscol + du_ofdcol
!
!         if (kdt <= 2 .and. me == master) then
!      print *, ' unified drag_suite_run ', kdt
!      print *, ' GSL drag du/dt ', maxval(Pdudt)*86400, minval(Pdudt)*86400
!      print *, ' GSL drag dv/dt ', maxval(Pdvdt)*86400, minval(Pdvdt)*86400
!
! zero      print *, ' unified drag_GSL dT/dt ', maxval(Pdtdt)*86400, minval(Pdtdt)*86400
!
!      if (gwd_opt == 22 .or. gwd_opt == 33) then
!      print *, ' unified drag_GSL dUBL/dt ',  maxval(dudt_obl)*86400, minval(dudt_obl)*86400
!      print *, ' unified drag_GSL dVBL/dt ',  maxval(dvdt_obl)*86400, minval(dvdt_obl)*86400
!      print *, ' unified drag_GSL dUOGW/dt ', maxval(dudt_ogw)*86400, minval(dudt_ogw)*86400
!      print *, ' unified drag_GSL dVOGW/dt ', maxval(dvdt_ogw)*86400, minval(dvdt_ogw)*86400
!      print *, ' unified drag_GSL dUOss/dt ', maxval(dudt_oss)*86400, minval(dudt_oss)*86400
!      print *, ' unified drag_GSL dVOSS/dt ', maxval(dvdt_oss)*86400, minval(dvdt_oss)*86400
!      print *, ' unified drag_GSL dUOfd/dt ', maxval(dudt_ofd)*86400, minval(dudt_ofd)*86400
!      print *, ' unified drag_GSL dVOfd/dt ', maxval(dvdt_ofd)*86400, minval(dvdt_ofd)*86400
!      endif
!     endif

    else
!
! not gsldrag oro-scheme for example "do_ugwp_v1_orog_only"
!

    if ( do_ugwp_v1_orog_only ) then
!
! for TOFD we use now "varss" of GSL-drag  /not sgh30=abs(oro-oro_f)/
! only sum of integrated ORO+GW effects (dusfcg and dvsfcg) = sum(ogw + obl + oss*0 + ofd + ngw)
!
! OROGW_V1 introduce "orchestration" between OGW-effects and Mountain Blocking
!      it starts to examines options for the Scale-Aware (SA)formulation of SSO-effects
!      if ( me == master .and. kdt == 1) print *, ' bf orogw_v1 nmtvr=', nmtvr, ' do_tofd=', do_tofd

         if (gwd_opt ==1 )sgh30 = 0.15*hprime       ! portion of the mesoscale SSO (~[oro_unfilt -oro_filt)
         if (gwd_opt >1 ) sgh30 = varss             ! as in gsldrag: see drag_suite_run

       call orogw_v1 (im, levs,  lonr,  me, master,dtp, kdt, do_tofd,     &
                      xlat_d, sinlat, coslat, area,                       &
                      cdmbgwd(1:2), hprime, oc, oa4, clx, theta,          &
                      sigma, gamma, elvmax,  sgh30,  kpbl, ugrs,          &
                      vgrs, tgrs, q1, prsi,del,prsl,prslk, zmeti, zmet,   &
                      Pdvdt, Pdudt, Pdtdt, Pkdis, DUSFCg, DVSFCg,rdxzb,   &
                      zobl, zlwb, zogw, tau_ogw, dudt_ogw, dvdt_ogw,      &
                      dudt_obl, dvdt_obl,dudt_ofd, dvdt_ofd,              &
                      du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol,         &
                      du_ofdcol, dv_ofdcol, errmsg,errflg           )
!
! orogw_v1: dusfcg = du_ogwcol + du_oblcol  + du_ofdcol                           only 3 terms
!
!
!          if (kdt <= 2 .and. me == master) then
!
!       print *, ' unified_ugwp orogw_v1 ', kdt, me,  nmtvr
!       print *, ' unified_ugwp orogw_v1 du/dt ', maxval(Pdudt)*86400, minval(Pdudt)*86400
!       print *, ' unified_ugwp orogw_v1 dv/dt ', maxval(Pdvdt)*86400, minval(Pdvdt)*86400
!       print *, ' unified_ugwp orogw_v1 dT/dt ', maxval(Pdtdt)*86400, minval(Pdtdt)*86400
!       print *, ' unified_ugwp orogw_v1 dUBL/dt ', maxval(dudt_obl)*86400, minval(dudt_obl)*86400
!       print *, ' unified_ugwp orogw_v1 dVBL/dt ', maxval(dvdt_obl)*86400, minval(dvdt_obl)*86400
!      endif


    end if
!
!  for  old-fashioned GFS-style diag-cs like dt3dt(:.:, 1:14) collections
!
     if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
        do k=1,levs
          do i=1,im
             ldu3dt_ogw(i,k) = ldu3dt_ogw(i,k) + Pdudt(i,k)*dtp
             ldv3dt_ogw(i,k) = ldv3dt_ogw(i,k) + Pdvdt(i,k)*dtp
             ldt3dt_ogw(i,k) = ldt3dt_ogw(i,k) + Pdtdt(i,k)*dtp
          enddo
        enddo
      endif
   ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin non-stationary GW schemes
! ugwp_v1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (do_ugwp_v1) then

!==================================================================
!       call slat_geos5_tamp_v1(im, tamp_mpa, xlat_d, tau_ngw)
!
! 2020 updates of MERRA/GEOS tau_ngw for the C96-QBO FV3GFS-127L runs
!==================================================================

       call  slat_geos5_2020(im, tamp_mpa, xlat_d, tau_ngw)

       y4 = jdat(1); month = jdat(2); day = jdat(3)
!
! hour = jdat(5)
! fhour = float(hour)+float(jdat(6))/60. + float(jdat(7))/3600.
!       fhour = (kdt-1)*dtp/3600.
!       fhrday  = fhour/24.  - nint(fhour/24.)


       call calendar_ugwp(y4, month, day, ddd_ugwp)
       curdate = y4*1000 + ddd_ugwp
!
       call ngwflux_update(me, master, im, levs, kdt, ddd_ugwp,curdate, &
         tau_amf, xlat_d, sinlat,coslat, rain, tau_ngw)

       call cires_ugwpv1_ngw_solv2(me, master, im,   levs,  kdt, dtp,   &
                      tau_ngw, tgrs, ugrs,  vgrs,   q1, prsl, prsi,     &
                      zmet, zmeti,prslk,   xlat_d, sinlat, coslat,      &
                      dudt_ngw, dvdt_ngw, dtdt_ngw, kdis_ngw, zngw)
!
! =>  con_g, con_cp, con_rd, con_rv, con_omega,  con_pi, con_fvirt
!
!       if (me == master .and. kdt <= 2) then
!         print *
!         write(6,*)'FV3GFS finished fv3_ugwp_solv2_v1   '
!         write(6,*) ' non-stationary GWs with GMAO/MERRA GW-forcing '
!         print *
!
!      print *, ' ugwp_v1 ', kdt
!      print *, ' ugwp_v1 du/dt ', maxval(dudt_ngw)*86400, minval(dudt_ngw)*86400
!      print *, ' ugwp_v1 dv/dt ', maxval(dvdt_ngw)*86400, minval(dvdt_ngw)*86400
!      print *, ' ugwp_v1 dT/dt ', maxval(dtdt_ngw)*86400, minval(dtdt_ngw)*86400
!       endif


    end if   ! do_ugwp_v1

!
!  GFS-style diag dt3dt(:.:, 1:14)  time-averaged
!
      if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
        do k=1,levs
          do i=1,im
             ldu3dt_ngw(i,k) = ldu3dt_ngw(i,k) + dudt_ngw(i,k)*dtp
             ldv3dt_ngw(i,k) = ldv3dt_ngw(i,k) + dvdt_ngw(i,k)*dtp
             ldt3dt_ngw(i,k) = ldt3dt_ngw(i,k) + dtdt_ngw(i,k)*dtp
          enddo
        enddo
      endif

!
! get total sso-OGW + NGW
!
     dudt_gw =  Pdudt +dudt_ngw
     dvdt_gw =  Pdvdt +dvdt_ngw
     dtdt_gw =  Pdtdt +dtdt_ngw
     kdis_gw =  Pkdis +kdis_ngw
!
! accumulate "tendencies" as in the GFS-ipd (pbl + ugwp + zero-RF)
!
     dudt  = dudt  + dudt_ngw
     dvdt  = dvdt  + dvdt_ngw
     dtdt  = dtdt  + dtdt_ngw

    end subroutine ugwpv1_gsldrag_run
!! @}
!>@}
end module ugwpv1_gsldrag

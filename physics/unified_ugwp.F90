!>  \file unified_ugwp.F90
!! This file combines three gravity wave drag schemes under one ("unified_ugwp") suite:
!!      1) The "V0 CIRES UGWP" scheme (cires_ugwp.F90) as implemented in the FV3GFSv16 atmosphere model, which includes:
!!            a) the "traditional" EMC orograhic gravity wave drag and flow blocking scheme of gwdps.f
!!            b) the v0 cires ugwp non-stationary GWD scheme
!!      2) The GSL orographic drag suite (drag_suite.F90), as implmeneted in the RAP/HRRR, which includes:
!!            a) large-scale gravity wave drag and low-level flow blocking -- active at horizontal scales
!!               down to ~5km (Kim and Arakawa, 1995 \cite kim_and_arakawa_1995; Kim and Doyle, 2005 \cite kim_and_doyle_2005)
!!            b) small-scale gravity wave drag scheme -- active typically in stable PBL at horizontal grid resolutions down to ~1km
!!               (Steeneveld et al, 2008 \cite steeneveld_et_al_2008; Tsiringakis et al, 2017 \cite tsiringakis_et_al_2017)
!!            c) turbulent orographic form drag -- active at horizontal grid ersolutions down to ~1km
!!               (Beljaars et al, 2004 \cite beljaars_et_al_2004)
!!      3) The "V1 CIRES UGWP" scheme developed by Valery Yudin (University of Colorado, CIRES)
!! See Valery Yudin's presentation at 2017 NGGPS PI meeting:
!! Gravity waves (GWs): Mesoscale GWs transport momentum, energy (heat) , and create eddy mixing in the whole atmosphere domain; Breaking and dissipating GWs deposit: (a) momentum; (b) heat (energy); and create (c) turbulent mixing of momentum, heat, and tracers
!! To properly incorporate GW effects (a-c) unresolved by DYCOREs we need GW physics
!! "Unified": a) all GW effects due to both dissipation/breaking; b) identical GW solvers for all GW sources; c) ability to replace solvers.
!! Unified Formalism:
!! 1. GW Sources: Stochastic and physics based mechanisms for GW-excitations in the lower atmosphere, calibrated by the high-res analyses/forecasts, and observations (3 types of GW sources: orography, convection, fronts/jets).
!! 2. GW Propagation: Unified solver for "propagation, dissipation and breaking" excited from all type of GW sources.
!! 3. GW Effects: Unified representation of GW impacts on the "resolved" flow for all sources (energy-balanced schemes for momentum, heat and mixing).
!! https://www.weather.gov/media/sti/nggps/Presentations%202017/02%20NGGPS_VYUDIN_2017_.pdf
!!
!! The unified_ugwp scheme is activated by gwd_opt = 2 in the namelist.
!! The choice of schemes is activated at runtime by the following namelist options (boolean):
!!       do_ugwp_v0           -- activates V0 CIRES UGWP scheme - both orographic and non-stationary GWD
!!       do_ugwp_v0_orog_only -- activates V0 CIRES UGWP scheme - orographic GWD only
!!       do_gsl_drag_ls_bl    -- activates RAP/HRRR (GSL) large-scale GWD and blocking
!!       do_gsl_drag_ss       -- activates RAP/HRRR (GSL) small-scale GWD
!!       do_gsl_drag_tofd     -- activates RAP/HRRR (GSL) turbulent orographic drag
!!       do_ugwp_v1           -- activates V1 CIRES UGWP scheme - both orographic and non-stationary GWD
!!       do_ugwp_v1_orog_only -- activates V1 CIRES UGWP scheme - orographic GWD only
!! Note that only one "large-scale" scheme can be activated at a time.
!!

module unified_ugwp

    use machine, only: kind_phys

    use cires_ugwp_module, only: knob_ugwp_version, cires_ugwp_mod_init, cires_ugwp_mod_finalize

    use cires_ugwp_module_v1, only: cires_ugwp_init_v1, cires_ugwp_finalize, calendar_ugwp

    use gwdps, only: gwdps_run

    use drag_suite, only: drag_suite_run

    use cires_ugwp_orolm97_v1, only: gwdps_oro_v1

    use cires_ugwp_triggers_v1, only: slat_geos5_tamp_v1
    
    ! use cires_ugwp_ngw_utils, only: tau_limb_advance
    
    use cires_ugwp_solv2_v1_mod, only: cires_ugwp_solv2_v1

    implicit none

    private

    public unified_ugwp_init, unified_ugwp_run, unified_ugwp_finalize

    logical :: is_initialized = .False.

contains

! ------------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme v0
! ------------------------------------------------------------------------
!>@brief The subroutine initializes the unified UGWP
!> \section arg_table_unified_ugwp_init Argument Table
!! \htmlinclude unified_ugwp_init.html
!!
! -----------------------------------------------------------------------
!
    subroutine unified_ugwp_init (me, master, nlunit, input_nml_file, logunit, &
                fn_nml2, jdat, lonr, latr, levs, ak, bk, dtp, cdmbgwd, cgwf,   &
                con_pi, con_rerth, pa_rf_in, tau_rf_in, con_p0, do_ugwp,       &
                do_ugwp_v0, do_ugwp_v0_orog_only, do_gsl_drag_ls_bl,           &
                do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1,                  &
                do_ugwp_v1_orog_only, errmsg, errflg)

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
    real(kind=kind_phys), intent (in) :: ak(:), bk(:)
    real(kind=kind_phys), intent (in) :: dtp
    real(kind=kind_phys), intent (in) :: cdmbgwd(4), cgwf(2) ! "scaling" controls for "old" GFS-GW schemes
    real(kind=kind_phys), intent (in) :: pa_rf_in, tau_rf_in
    real(kind=kind_phys), intent (in) :: con_p0, con_pi, con_rerth
    logical,              intent (in) :: do_ugwp
    logical,              intent (in) :: do_ugwp_v0, do_ugwp_v0_orog_only,  &
                                         do_gsl_drag_ls_bl, do_gsl_drag_ss, &
                                         do_gsl_drag_tofd, do_ugwp_v1,      &
                                         do_ugwp_v1_orog_only

    character(len=*), intent (in) :: fn_nml2
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


    ! Test to make sure that at most only one large-scale/blocking
    ! orographic drag scheme is chosen
    if ( (do_ugwp_v0.and.(do_ugwp_v0_orog_only.or.do_gsl_drag_ls_bl.or.    &
                          do_ugwp_v1.or.do_ugwp_v1_orog_only))        .or. &
         (do_ugwp_v0_orog_only.and.(do_gsl_drag_ls_bl.or.do_ugwp_v1.or.    &
                                    do_ugwp_v1_orog_only))            .or. &
         (do_gsl_drag_ls_bl.and.(do_ugwp_v1.or.do_ugwp_v1_orog_only)) .or. &
         (do_ugwp_v1.and.do_ugwp_v1_orog_only) ) then

       write(errmsg,'(*(a))') "Logic error: Only one large-scale&
          &/blocking scheme (do_ugwp_v0,do_ugwp_v0_orog_only,&
          &do_gsl_drag_ls_bl,do_ugwp_v1 or &
          &do_ugwp_v1_orog_only) can be chosen"
       errflg = 1
       return

    end if


    if (is_initialized) return


    if ( do_ugwp_v0 ) then
       ! if (do_ugwp .or. cdmbgwd(3) > 0.0) then (deactivate effect of do_ugwp)
       if (cdmbgwd(3) > 0.0) then
         call cires_ugwp_mod_init (me, master, nlunit, input_nml_file, logunit, &
                                fn_nml2, lonr, latr, levs, ak, bk, con_p0, dtp, &
                                cdmbgwd(1:2), cgwf, pa_rf_in, tau_rf_in)
       else
         write(errmsg,'(*(a))') "Logic error: cires_ugwp_mod_init called but &
               &do_ugwp_v0 is true and cdmbgwd(3) <= 0"
         errflg = 1
         return
       end if
    end if


    if ( do_ugwp_v1 ) then
       call cires_ugwp_init_v1 (me, master, nlunit, logunit, jdat, con_pi,      &
                                con_rerth, fn_nml2, lonr, latr, levs, ak, bk,   &
                                con_p0, dtp, cdmbgwd(1:2), cgwf, pa_rf_in,      &
                                tau_rf_in, errmsg, errflg)
    end if

    is_initialized = .true.

    end subroutine unified_ugwp_init


! -----------------------------------------------------------------------
! finalize of unified_ugwp   (_finalize)
! -----------------------------------------------------------------------

!>@brief The subroutine finalizes the CIRES UGWP

!> \section arg_table_unified_ugwp_finalize Argument Table
!! \htmlinclude unified_ugwp_finalize.html
!!

    subroutine unified_ugwp_finalize(do_ugwp_v0,do_ugwp_v1,errmsg, errflg)

    implicit none
!
    logical,          intent (in) :: do_ugwp_v0, do_ugwp_v1
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

    if ( do_ugwp_v0 ) call cires_ugwp_mod_finalize()

    if ( do_ugwp_v1 ) call cires_ugwp_finalize()

    is_initialized = .false.

    end subroutine unified_ugwp_finalize


! -----------------------------------------------------------------------
!    originally from ugwp_driver_v0.f
!    driver of cires_ugwp   (_driver)
! -----------------------------------------------------------------------
!   driver is called after pbl & before chem-parameterizations
! -----------------------------------------------------------------------
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
! -----------------------------------------------------------------------
!>@brief These subroutines and modules execute the CIRES UGWP Version 0
!>\defgroup unified_ugwp_run Unified Gravity Wave Physics General Algorithm
!> @{
!! The physics of NGWs in the UGWP framework (Yudin et al. 2018 \cite yudin_et_al_2018) is represented by four GW-solvers, which is introduced in Lindzen (1981) \cite lindzen_1981, Hines (1997) \cite hines_1997, Alexander and Dunkerton (1999) \cite alexander_and_dunkerton_1999, and Scinocca (2003) \cite scinocca_2003. The major modification of these GW solvers is represented by the addition of the background dissipation of temperature and winds to the saturation criteria for wave breaking. This feature is important in the mesosphere and thermosphere for WAM applications and it considers appropriate scale-dependent dissipation of waves near the model top lid providing the momentum and energy conservation in the vertical column physics (Shaw and Shepherd 2009 \cite shaw_and_shepherd_2009). In the UGWP-v0, the modification of Scinocca (2003) \cite scinocca_2003 scheme for NGWs with non-hydrostatic and rotational effects for GW propagations and background dissipation is represented by the subroutine \ref fv3_ugwp_solv2_v0. In the next release of UGWP, additional GW-solvers will be implemented along with physics-based triggering of waves and stochastic approaches for selection of GW modes characterized by horizontal phase velocities, azimuthal directions and magnitude of the vertical momentum flux (VMF).
!!
!! In UGWP-v0, the specification for the VMF function is adopted from the GEOS-5 global atmosphere model of GMAO NASA/GSFC, as described in Molod et al. (2015) \cite molod_et_al_2015 and employed in the MERRRA-2 reanalysis (Gelaro et al., 2017 \cite gelaro_et_al_2017). The Fortran subroutine \ref slat_geos5_tamp describes the latitudinal shape of VMF-function as displayed in Figure 3 of Molod et al. (2015) \cite molod_et_al_2015. It shows that the enhanced values of VMF in the equatorial region gives opportunity to simulate the QBO-like oscillations in the equatorial zonal winds and lead to more realistic simulations of the equatorial dynamics in GEOS-5 operational and MERRA-2 reanalysis products. For the first vertically extended version of FV3GFS in the stratosphere and mesosphere, this simplified function of VMF allows us to tune the model climate and to evaluate multi-year simulations of FV3GFS with the MERRA-2 and ERA-5 reanalysis products, along with temperature, ozone, and water vapor observations of current satellite missions. After delivery of the UGWP-code, the EMC group developed and tested approach to modulate the zonal mean NGW forcing by 3D-distributions of the total precipitation as a proxy for the excitation of NGWs by convection and the vertically-integrated  (surface - tropopause) Turbulent Kinetic Energy (TKE). The verification scores with updated NGW forcing, as reported elsewhere by EMC researchers, display noticeable improvements in the forecast scores produced by FV3GFS configuration extended into the mesosphere.
!!
!> \section arg_table_unified_ugwp_run Argument Table
!! \htmlinclude unified_ugwp_run.html
!!
!> \section gen_unified_ugwp CIRES UGWP Scheme General Algorithm
!! @{
     subroutine unified_ugwp_run(me,  master, im,  levs, ntrac, dtp, fhzero, kdt,      &
         lonr, oro, oro_uf, hprime, nmtvr, oc, theta, sigma, gamma, elvmax, clx, oa4,  &
         varss,oc1ss,oa4ss,ol4ss,dx,dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,dusfc_ss,      &
         dvsfc_ss,dusfc_fd,dvsfc_fd,dtaux2d_ls,dtauy2d_ls,dtaux2d_bl,dtauy2d_bl,       &
         dtaux2d_ss,dtauy2d_ss,dtaux2d_fd,dtauy2d_fd,br1,hpbl,slmsk,                   &
         do_tofd, ldiag_ugwp, cdmbgwd, jdat, xlat, xlat_d, sinlat, coslat, area,       &
         ugrs, vgrs, tgrs, q1, prsi, prsl, prslk, phii, phil,                          &
         del, kpbl, dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,                &
         tau_tofd, tau_mtb, tau_ogw, tau_ngw, zmtb, zlwb, zogw,                        &
         dudt_mtb,dudt_ogw, dudt_tms, du3dt_mtb, du3dt_ogw, du3dt_tms,                 &
         dudt, dvdt, dtdt, rdxzb, con_g, con_omega, con_pi, con_cp, con_rd, con_rv,    &
         con_rerth, con_fvirt, rain, ntke, q_tke, dqdt_tke, lprnt, ipr,                &
         ldu3dt_ogw, ldv3dt_ogw, ldt3dt_ogw, ldu3dt_cgw, ldv3dt_cgw, ldt3dt_cgw,       &
         ldiag3d, lssav, flag_for_gwd_generic_tend, do_ugwp_v0, do_ugwp_v0_orog_only,  &
         do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1,              &
         do_ugwp_v1_orog_only, gwd_opt, errmsg, errflg)

    implicit none

    ! interface variables
    integer,                 intent(in) :: me, master, im, levs, ntrac, kdt, lonr, nmtvr
    integer,                 intent(in) :: gwd_opt
    integer,                 intent(in), dimension(im)       :: kpbl
    real(kind=kind_phys),    intent(in), dimension(im)       :: oro, oro_uf, hprime, oc, theta, sigma, gamma
    real(kind=kind_phys),    intent(in), dimension(im)       :: varss,oc1ss,oa4ss,ol4ss,dx
    logical,                 intent(in)                      :: flag_for_gwd_generic_tend
    ! elvmax is intent(in) for CIRES UGWP, but intent(inout) for GFS GWDPS
    real(kind=kind_phys),    intent(inout), dimension(im)    :: elvmax
    real(kind=kind_phys),    intent(in), dimension(im, 4)    :: clx, oa4
    real(kind=kind_phys),    intent(in), dimension(im)       :: xlat, xlat_d, sinlat, coslat, area
    real(kind=kind_phys),    intent(in), dimension(im, levs) :: del, ugrs, vgrs, tgrs, prsl, prslk, phil
    real(kind=kind_phys),    intent(in), dimension(im, levs+1) :: prsi, phii
    real(kind=kind_phys),    intent(in), dimension(im, levs) :: q1
    real(kind=kind_phys),    intent(in) :: dtp, fhzero, cdmbgwd(4)
    integer, intent(in) :: jdat(8)
    logical,                 intent(in) :: do_tofd, ldiag_ugwp

!Output (optional):
    real(kind=kind_phys), intent(out) ::                          &
      &                      dusfc_ls(:),dvsfc_ls(:),             &
      &                      dusfc_bl(:),dvsfc_bl(:),             &
      &                      dusfc_ss(:),dvsfc_ss(:),             &
      &                      dusfc_fd(:),dvsfc_fd(:)
    real(kind=kind_phys), intent(out) ::                          &
      &         dtaux2d_ls(:,:),dtauy2d_ls(:,:),                  &
      &         dtaux2d_bl(:,:),dtauy2d_bl(:,:),                  &
      &         dtaux2d_ss(:,:),dtauy2d_ss(:,:),                  &
      &         dtaux2d_fd(:,:),dtauy2d_fd(:,:)

    real(kind=kind_phys), intent(in) ::     br1(im),              &
      &                                     hpbl(im),             &
      &                                     slmsk(im)

    real(kind=kind_phys),    intent(out), dimension(im)      :: dusfcg, dvsfcg
    real(kind=kind_phys),    intent(out), dimension(im)      :: zmtb, zlwb, zogw, rdxzb
    real(kind=kind_phys),    intent(out), dimension(im)      :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
    real(kind=kind_phys),    intent(out), dimension(im, levs):: gw_dudt, gw_dvdt, gw_dtdt, gw_kdis
    real(kind=kind_phys),    intent(out), dimension(im, levs):: dudt_mtb, dudt_ogw, dudt_tms

    ! These arrays are only allocated if ldiag=.true.
    real(kind=kind_phys),    intent(inout), dimension(:,:)      :: ldu3dt_ogw, ldv3dt_ogw, ldt3dt_ogw
    real(kind=kind_phys),    intent(inout), dimension(:,:)      :: ldu3dt_cgw, ldv3dt_cgw, ldt3dt_cgw
    logical,                 intent(in)                         :: ldiag3d, lssav

    ! These arrays only allocated if ldiag_ugwp = .true.
    real(kind=kind_phys),    intent(out), dimension(:,:) :: du3dt_mtb, du3dt_ogw, du3dt_tms

    real(kind=kind_phys),    intent(inout), dimension(im, levs):: dudt, dvdt, dtdt

    real(kind=kind_phys),    intent(in) :: con_g, con_omega, con_pi, con_cp, con_rd, &
                                           con_rv, con_rerth, con_fvirt

    real(kind=kind_phys),    intent(in), dimension(im) :: rain

    integer,                 intent(in) :: ntke
    real(kind=kind_phys),    intent(in), dimension(:,:) :: q_tke, dqdt_tke

    logical, intent(in) :: lprnt
    integer, intent(in) :: ipr

    ! flags for choosing combination of GW drag schemes to run
    logical,              intent (in) :: do_ugwp_v0, do_ugwp_v0_orog_only,  &
                                         do_gsl_drag_ls_bl, do_gsl_drag_ss, &
                                         do_gsl_drag_tofd, do_ugwp_v1,      &
                                         do_ugwp_v1_orog_only

    character(len=*),        intent(out) :: errmsg
    integer,                 intent(out) :: errflg

    ! local variables
    integer :: i, k
    real(kind=kind_phys), dimension(im)       :: sgh30
    real(kind=kind_phys), dimension(im, levs) :: Pdvdt, Pdudt
    real(kind=kind_phys), dimension(im, levs) :: Pdtdt, Pkdis
    real(kind=kind_phys), dimension(im, levs) :: ed_dudt, ed_dvdt, ed_dtdt
    ! from ugwp_driver_v0.f -> cires_ugwp_initialize.F90 -> module ugwp_wmsdis_init
    real(kind=kind_phys), parameter :: tamp_mpa=30.e-3
    ! switches that activate impact of OGWs and NGWs (WL* how to deal with them? *WL)
    real(kind=kind_phys), parameter :: pogw=1., pngw=1., pked=1.
    real(kind=kind_phys), parameter :: fw1_tau=1.0
    integer :: nmtvr_temp

    real(kind=kind_phys), dimension(:,:), allocatable :: tke
    real(kind=kind_phys), dimension(:),   allocatable :: turb_fac, tem
    real(kind=kind_phys) :: rfac, tx1

    real(kind=kind_phys) :: inv_g
    real(kind=kind_phys), dimension(im, levs)   :: zmet  ! geopotential height at model Layer centers
    real(kind=kind_phys), dimension(im, levs+1) :: zmeti ! geopotential height at model layer interfaces


    ! ugwp_v1 local variables
    integer :: y4, month, day,  ddd_ugwp, curdate, curday
    integer :: hour
    real(kind=kind_phys) :: hcurdate, hcurday, fhour, fhrday
    integer :: kdtrest
    integer :: curday_ugwp
    integer :: curday_save=20150101
    logical :: first_qbo=.true.
    real    ::  hcurday_save =20150101.00
    save first_qbo, curday_save, hcurday_save


    ! ugwp_v1 temporary (local) diagnostic variables from cires_ugwp_solv2_v1
    real(kind=kind_phys) :: tauabs(im,levs), wrms(im,levs), trms(im,levs)


    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! 1) ORO stationary GWs
    !    ------------------

     zlwb(:)   = 0.

    ! Run the appropriate large-scale (large-scale GWD + blocking) scheme
    ! Note:  In case of GSL drag_suite, this includes ss and tofd

    if ( do_gsl_drag_ls_bl.or.do_gsl_drag_ss.or.do_gsl_drag_tofd ) then

       call drag_suite_run(im,levs,dvdt,dudt,dtdt,ugrs,vgrs,tgrs,q1, &
                 kpbl,prsi,del,prsl,prslk,phii,phil,dtp,             &
                 kdt,hprime,oc,oa4,clx,varss,oc1ss,oa4ss,            &
                 ol4ss,theta,sigma,gamma,elvmax,dtaux2d_ls,          &
                 dtauy2d_ls,dtaux2d_bl,dtauy2d_bl,dtaux2d_ss,        &
                 dtauy2d_ss,dtaux2d_fd,dtauy2d_fd,dusfcg,            &
                 dvsfcg,dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,         &
                 dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd,                &
                 slmsk,br1,hpbl,con_g,con_cp,con_rd,con_rv,          &
                 con_fvirt,con_pi,lonr,                              &
                 cdmbgwd(1:2),me,master,lprnt,ipr,rdxzb,dx,gwd_opt,  &
                 do_gsl_drag_ls_bl,do_gsl_drag_ss,do_gsl_drag_tofd,  &
                 errmsg,errflg)

    end if

    if ( do_ugwp_v1.or.do_ugwp_v1_orog_only ) then

       ! Valery's TOFD
       ! topo paras
       ! w/ orographic effects
       if(nmtvr == 14)then
         ! calculate sgh30 for TOFD
         sgh30 = abs(oro - oro_uf)
       ! w/o orographic effects
       else
         sgh30   = 0.
       endif

       inv_g = 1./con_g
       zmeti  = phii*inv_g
       zmet   = phil*inv_g

       call gwdps_oro_v1 (im, levs,  lonr,   do_tofd,                  &
                      Pdvdt, Pdudt, Pdtdt, Pkdis,                      &
                      ugrs , vgrs, tgrs, q1, KPBL, prsi,del,prsl,      &
                      prslk, zmeti, zmet, dtp, kdt, hprime, oc, oa4,   &
                      clx, theta, sigma, gamma, elvmax,                &
                      con_g, con_omega, con_rd, con_cp, con_rv,con_pi, &
                      con_rerth, con_fvirt, sgh30, DUSFCg, DVSFCg,     &
                      xlat_d, sinlat, coslat, area,cdmbgwd(1:2), me,   &
                      master, rdxzb, zmtb, zogw, tau_mtb, tau_ogw,     &
                      tau_tofd, du3dt_mtb, du3dt_ogw, du3dt_tms)

    end if

    if ( do_ugwp_v0.or.do_ugwp_v0_orog_only ) then

      do k=1,levs
        do i=1,im
          Pdvdt(i,k) = 0.0
          Pdudt(i,k) = 0.0
          Pdtdt(i,k) = 0.0
          Pkdis(i,k) = 0.0
        enddo
      enddo

      if (cdmbgwd(1) > 0.0 .or. cdmbgwd(2) > 0.0) then

        ! Override nmtvr with nmtvr_temp = 14 for passing into gwdps_run if necessary
        if ( nmtvr == 24 ) then  ! gwd_opt = 2, 22, 3, or 33
           nmtvr_temp = 14
        else
           nmtvr_temp = nmtvr
        end if

        call gwdps_run(im, levs, Pdvdt, Pdudt, Pdtdt,                  &
                   ugrs, vgrs, tgrs, q1,                               &
                   kpbl, prsi, del, prsl, prslk, phii, phil, dtp, kdt, &
                   hprime, oc, oa4, clx, theta, sigma, gamma,          &
                   elvmax, dusfcg, dvsfcg,                             &
                   con_g,  con_cp, con_rd, con_rv, lonr,               &
                   nmtvr_temp, cdmbgwd, me, lprnt, ipr, rdxzb,         &
                   errmsg, errflg)
        if (errflg/=0) return
      endif

      tau_mtb   = 0.0  ; tau_ogw   = 0.0 ;  tau_tofd = 0.0
      if (ldiag_ugwp) then
        du3dt_mtb = 0.0  ; du3dt_ogw = 0.0 ;  du3dt_tms= 0.0
      end if


      if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
        do k=1,levs
          do i=1,im
             ldu3dt_ogw(i,k) = ldu3dt_ogw(i,k) + Pdudt(i,k)*dtp
             ldv3dt_ogw(i,k) = ldv3dt_ogw(i,k) + Pdvdt(i,k)*dtp
             ldt3dt_ogw(i,k) = ldt3dt_ogw(i,k) + Pdtdt(i,k)*dtp
          enddo
        enddo
      endif
   
    end if 



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin non-stationary GW schemes
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! ugwp_v0 non-stationary GW drag
    !
    if (do_ugwp_v0) then

      if (cdmbgwd(3) > 0.0) then

        ! 2) non-stationary GW-scheme with GMAO/MERRA GW-forcing
        call slat_geos5_tamp(im, tamp_mpa, xlat_d, tau_ngw)

        if (abs(1.0-cdmbgwd(3)) > 1.0e-6) then
          if (cdmbgwd(4) > 0.0) then
            allocate(turb_fac(im))
            do i=1,im
              turb_fac(i) = 0.0
            enddo
            if (ntke > 0) then
              allocate(tke(im,levs))
              allocate(tem(im))
              tke(:,:) = q_tke(:,:) + dqdt_tke(:,:) * dtp
              tem(:)   = 0.0
              do k=1,(levs+levs)/3
                do i=1,im
                  turb_fac(i) = turb_fac(i) + del(i,k) * tke(i,k)
                  tem(i)      = tem(i)      + del(i,k)
                enddo
              enddo
              do i=1,im
                turb_fac(i) = turb_fac(i) / tem(i)
              enddo
              deallocate(tke)
              deallocate(tem)
            endif
            rfac = 86400000 / dtp
            do i=1,im
              tx1 = cdmbgwd(4)*min(10.0, max(turb_fac(i),rain(i)*rfac))
              tau_ngw(i) = tau_ngw(i) * max(0.1, min(5.0, tx1))
            enddo
            deallocate(turb_fac)
          endif
          do i=1,im
            tau_ngw(i) = tau_ngw(i) * cdmbgwd(3)
          enddo
        endif

        call fv3_ugwp_solv2_v0(im, levs, dtp, tgrs, ugrs, vgrs, q1,                        &
             prsl, prsi, phil, xlat_d, sinlat, coslat, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis, &
             tau_ngw, me, master, kdt)

        do k=1,levs
          do i=1,im
            gw_dtdt(i,k) = pngw*gw_dtdt(i,k)+ pogw*Pdtdt(i,k)
            gw_dudt(i,k) = pngw*gw_dudt(i,k)+ pogw*Pdudt(i,k)
            gw_dvdt(i,k) = pngw*gw_dvdt(i,k)+ pogw*Pdvdt(i,k)
            gw_kdis(i,k) = pngw*gw_kdis(i,k)+ pogw*Pkdis(i,k)
            ! accumulation of tendencies for CCPP to replicate EMC-physics updates (!! removed in latest code commit to VLAB)
            !dudt(i,k) = dudt(i,k) +gw_dudt(i,k)
            !dvdt(i,k) = dvdt(i,k) +gw_dvdt(i,k)
            !dtdt(i,k) = dtdt(i,k) +gw_dtdt(i,k)
          enddo
        enddo

      else  ! .not.(cdmbgwd(3) > 0.0)

        do k=1,levs
          do i=1,im
            gw_dtdt(i,k) = Pdtdt(i,k)
            gw_dudt(i,k) = Pdudt(i,k)
            gw_dvdt(i,k) = Pdvdt(i,k)
            gw_kdis(i,k) = Pkdis(i,k)
          enddo
        enddo

      endif  ! cdmbgwd(3) > 0.0

      if (pogw == 0.0) then
        tau_mtb  = 0. ; tau_ogw  = 0. ; tau_tofd = 0.
        dudt_mtb = 0. ; dudt_ogw = 0. ; dudt_tms = 0.
      endif

#if 0
    !=============================================================================
    ! make "ugwp eddy-diffusion" update for gw_dtdt/gw_dudt/gw_dvdt by solving
    ! vert diffusion equations & update "Statein%tgrs, Statein%ugrs, Statein%vgrs"
    !=============================================================================
    ! 3) application of "eddy"-diffusion to "smooth" UGWP-related tendencies
    !------------------------------------------------------------------------------
      do k=1,levs
        do i=1,im
          ed_dudt(i,k) = 0.0 ; ed_dvdt(i,k) = 0.0 ; ed_dtdt(i,k) = 0.0
        enddo
      enddo

      call edmix_ugwp_v0(im, levs, dtp, tgrs, ugrs, vgrs, q1,                &
           del, prsl, prsi, phil, prslk, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis, &
           ed_dudt, ed_dvdt, ed_dtdt, me, master, kdt)
      gw_dtdt = gw_dtdt*(1.-pked) +  ed_dtdt*pked
      gw_dvdt = gw_dvdt*(1.-pked) +  ed_dvdt*pked
      gw_dudt = gw_dudt*(1.-pked) +  ed_dudt*pked
#endif

      if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
        do k=1,levs
          do i=1,im
             ldu3dt_cgw(i,k) = ldu3dt_cgw(i,k) + (gw_dudt(i,k) - Pdudt(i,k))*dtp
             ldv3dt_cgw(i,k) = ldv3dt_cgw(i,k) + (gw_dvdt(i,k) - Pdvdt(i,k))*dtp
             ldt3dt_cgw(i,k) = ldt3dt_cgw(i,k) + (gw_dtdt(i,k) - Pdtdt(i,k))*dtp
          enddo
        enddo
      endif

    end if  ! do_ugwp_v0 


    !
    ! ugwp_v1 non-stationary GW drag
    !
    if (do_ugwp_v1) then

! --------    
! 2) non-stationary GWs with GEOS-5/MERRA GW-forcing
!    ----------------------------------------------
!--------       
! GMAO GEOS-5/MERRA GW-forcing  lat-dep
!--------
       call slat_geos5_tamp_v1(im, tamp_mpa, xlat_d, tau_ngw)
        
       y4 = jdat(1); month = jdat(2); day = jdat(3) ; hour = jdat(5)
        
       ! fhour = float(hour)+float(jdat(6))/60. + float(jdat(7))/3600.
       fhour = (kdt-1)*dtp/3600.
       fhrday  = fhour/24.  - nint(fhour/24.)
       fhour = fhrday*24.     
                
       call calendar_ugwp(y4, month, day, ddd_ugwp)
       curdate = y4*1000 + ddd_ugwp
       curday = y4*10000 + month*100 + day
       hcurdate = float(curdate) + fhrday
       hcurday  = float(curday)  + fhrday
!
       if (mod(fhour,fhzero) == 0 .or. first_qbo) then

         ! call tau_limb_advance(me, master, im, levs, ddd_ugwp, curdate,   &
         !    j1_tau, j2_tau, ddy_j1tau, ddy_j2tau, tau_sat, kdt )
        
         if (first_qbo) kdtrest = kdt        
         first_qbo = .false.
         curday_save = curday        
         hcurday_save= hcurday       
       endif
        
       ! tau_ngw = fw1_tau*tau_ngw + tau_sat*(1.-fw1_tau)
                
!       goto 111   
!       if (mod(fhour,fhzero) == 0 .or. first_qbo) then
                
!         call tau_qbo_advance(me, master, im, levs, ddd_ugwp, curdate, &
!            j1_tau, j2_tau, ddy_j1tau, ddy_j2tau,  j1_qbo, j2_qbo,      &
!            ddy_j1qbo, ddy_j2qbo, tau_sat, tau_qbo,  uqbo, ax_qbo, kdt     )
        
                
!         if (me == master) then
!           print *, ' curday_save first_qbo ', curday, curday_save, kdt   
!           print *, ' hcurdays ', hcurdate,  float(hour)/24.
!           print *, jdat(5), jdat(6),  jdat(7), (kdt-1)*dtp/3600., ' calendar '           
!!           print *, ' curday curday_ugwp first_qbo ', hcurday, first_qbo         
!!           print *, ' vay_tau-limb U' ,  maxval(uqbo), minval(uqbo) 
!!           print *, ' vay_tau-limb TS' , maxval(tau_sat), minval(tau_sat) 
!!           print *, ' vay_tau-limb TQ' , maxval(tau_qbo), minval(tau_qbo)                  
!         endif  


!         if (first_qbo) kdtrest = kdt        
!         first_qbo = .false.
!         curday_save = curday        
!         hcurday_save= hcurday       
!       endif
        

                

!       if (mod(kdt, 720) == 0 .and. me == master ) then       
!         print *, ' vay_qbo_U' , maxval(uqbo), minval(uqbo) , kdt                             
!       endif
        
!       wqbo = dtp/taurel
!       do k =1, levs
!!          sdexpz =  wqbo*vert_qbo(k)  
!          sdexpz =  0.25*vert_qbo(k)
!          do i=1, im
!!             if (dexpy(i) > 0.0) then
!             dforc = 0.25
!!             ugrs(i,k)  = ugrs(i,k)*(1.-dforc) + dforc*uqbo(i,levs+1-k)
!!             tgrs(i,k)  = tgrs(i,k)*(1.-dforc) + dforc*tqbo(i,levs+1-k)              
!!             endif 
!          enddo
!       enddo      
 
! 111      continue


       call cires_ugwp_solv2_v1(im,   levs,  dtp,                     &
                      tgrs, ugrs,  vgrs,   q1, prsl, prsi,            &
                      zmet, zmeti,prslk, xlat_d, sinlat, coslat,      &
                      con_g, con_cp, con_rd, con_rv, con_omega,       &
                      con_pi, con_fvirt,                              & 
                      gw_dudt, gw_dvdt, gw_dTdt, gw_kdis,             &
                      tauabs, wrms, trms,   tau_ngw, me, master, kdt)

       if (me == master .and. kdt < 2) then
         print *
         write(6,*)'FV3GFS finished fv3_ugwp_solv2_v1  in ugwp_driver_v0 '
         write(6,*) ' non-stationary GWs with GMAO/MERRA GW-forcing '
         print *
       endif

       do k=1,levs
         do i=1,im
           gw_dtdt(i,k) = pngw*gw_dtdt(i,k) + pogw*Pdtdt(i,k)
           gw_dudt(i,k) = pngw*gw_dudt(i,k) + pogw*Pdudt(i,k)
                !+(uqbo(i,levs+1-k)-ugrs(i,k))/21600.
           gw_dvdt(i,k) = pngw*gw_dvdt(i,k) + pogw*Pdvdt(i,k)
           gw_kdis(i,k) = pngw*gw_kdis(i,k)     ! + pogw*Pkdis(i,k)
         enddo
       enddo




       if (pogw == 0.0) then
!        zmtb = 0.;  zogw =0.
         tau_mtb   = 0.0  ; tau_ogw   = 0.0 ;  tau_tofd = 0.0
         du3dt_mtb = 0.0  ; du3dt_ogw = 0.0 ;  du3dt_tms= 0.0
       endif

!       return

!=============================================================================
! make "ugwp eddy-diffusion" update for gw_dtdt/gw_dudt/gw_dvdt by solving
! vert diffusion equations & update "Statein%tgrs, Statein%ugrs, Statein%vgrs"
!=============================================================================
!
! 3) application of "eddy"-diffusion to "smooth" UGWP-related tendencies
!------------------------------------------------------------------------------

!       ed_dudt(:,:) = 0.0 ; ed_dvdt(:,:) = 0.0 ; ed_dtdt(:,:) = 0.0



!       call edmix_ugwp_v1(im,   levs, dtp,                       &
!                        tgrs, ugrs, vgrs, q1, del,               &
!                         prsl, prsi, phil, prslk,                &
!                         gw_dudt, gw_dvdt, gw_dTdt, gw_kdis,     &
!                        ed_dudt, ed_dvdt, ed_dTdt,
!                         me, master, kdt )

!       do k=1,levs
!         do i=1,im
!           gw_dtdt(i,k) = gw_dtdt(i,k) + ed_dtdt(i,k)*pked
!           gw_dvdt(i,k) = gw_dvdt(i,k) + ed_dvdt(i,k)*pked
!           gw_dudt(i,k) = gw_dudt(i,k) + ed_dudt(i,k)*pked
!         enddo
!       enddo


    end if   ! do_ugwp_v1


    end subroutine unified_ugwp_run
!! @}
!>@}
end module unified_ugwp

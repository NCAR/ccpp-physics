!>  \file unified_ugwp.F90
!! This file combines three two orographic GW-schemes cires_ugwp.F90 and drag_suite.F90 under "unified_ugwp" suite:
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
!!       do_ugwp_v0_nst_only  -- activates V0 CIRES UGWP scheme - non-stationary GWD only
!!       do_gsl_drag_ls_bl    -- activates RAP/HRRR (GSL) large-scale GWD and blocking
!!       do_gsl_drag_ss       -- activates RAP/HRRR (GSL) small-scale GWD
!!       do_gsl_drag_tofd     -- activates RAP/HRRR (GSL) turbulent orographic drag
!! Note that only one "large-scale" scheme can be activated at a time.
!!

module unified_ugwp

    use machine, only: kind_phys

!    use cires_ugwp_module,   only: knob_ugwp_version, cires_ugwp_mod_init,   cires_ugwp_mod_finalize
    use cires_ugwpv0_module, only: knob_ugwp_version, cires_ugwpv0_mod_init, cires_ugwpv0_mod_finalize
    use gwdps, only: gwdps_run

    use drag_suite, only: drag_suite_run

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
                do_ugwp_v0, do_ugwp_v0_orog_only, do_ugwp_v0_nst_only,         &
                do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd,           &
                errmsg, errflg)

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
                                         do_ugwp_v0_nst_only,               &
                                         do_gsl_drag_ls_bl, do_gsl_drag_ss, &
                                         do_gsl_drag_tofd

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
    if ( (do_ugwp_v0.and.(do_ugwp_v0_orog_only.or.do_gsl_drag_ls_bl)) .or. &
         (do_ugwp_v0_orog_only.and.do_gsl_drag_ls_bl) ) then

       write(errmsg,'(*(a))') "Logic error: Only one large-scale&
          &/blocking scheme (do_ugwp_v0,do_ugwp_v0_orog_only,&
          &do_gsl_drag_ls_bl can be chosen"
       errflg = 1
       return

    end if


    if (is_initialized) return


    if ( do_ugwp_v0 .or. do_ugwp_v0_nst_only ) then
       ! if (do_ugwp .or. cdmbgwd(3) > 0.0) then (deactivate effect of do_ugwp)
       if (cdmbgwd(3) > 0.0) then
        call cires_ugwpv0_mod_init(me, master, nlunit, input_nml_file, logunit, &
                                fn_nml2, lonr, latr, levs, ak, bk, con_p0, dtp, &
                                cdmbgwd(1:2), cgwf, pa_rf_in, tau_rf_in)
       else
         write(errmsg,'(*(a))') "Logic error: cires_ugwp_mod_init called but &
               &do_ugwp_v0 or do_ugwp_v0_nst_only is true and cdmbgwd(3) <= 0"
         errflg = 1
         return
       end if
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

    subroutine unified_ugwp_finalize(do_ugwp_v0,do_ugwp_v0_nst_only,  &
                                     errmsg, errflg)

    implicit none
!
    logical,          intent (in) :: do_ugwp_v0, do_ugwp_v0_nst_only
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

    if ( do_ugwp_v0 .or. do_ugwp_v0_nst_only )  call cires_ugwpv0_mod_finalize()

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
!! The physics of NGWs in the UGWP framework (Yudin et al. 2018 \cite yudin_et_al_2018) is represented by four GW-solvers, which is introduced in Lindzen (1981) \cite lindzen_1981, Hines (1997) \cite hines_1997, Alexander and Dunkerton (1999) \cite alexander_and_dunkerton_1999, and Scinocca (2003) \cite scinocca_2003. The major modification of these GW solvers is represented by the addition of the background dissipation of temperature and winds to the saturation criteria for wave breaking. This feature is important in the mesosphere and thermosphere for WAM applications and it considers appropriate scale-dependent dissipation of waves near the model top lid providing the momentum and energy conservation in the vertical column physics (Shaw and Shepherd 2009 \cite shaw_and_shepherd_2009). In the UGWP-v0, the modification of Scinocca (2003) \cite scinocca_2003 scheme for NGWs with non-hydrostatic and rotational effects for GW propagations and backgroufnd dissipation is represented by the subroutine \ref fv3_ugwp_solv2_v0. In the next release of UGWP, additional GW-solvers will be implemented along with physics-based triggering of waves and stochastic approaches for selection of GW modes characterized by horizontal phase velocities, azimuthal directions and magnitude of the vertical momentum flux (VMF).
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
         dudt_mtb, dudt_tms, du3dt_mtb, du3dt_ogw, du3dt_tms,                          &
         dudt, dvdt, dtdt, rdxzb, con_g, con_omega, con_pi, con_cp, con_rd, con_rv,    &
         con_rerth, con_fvirt, rain, ntke, q_tke, dqdt_tke, lprnt, ipr,                &
         dtend, dtidx, index_of_temperature, index_of_x_wind, index_of_y_wind,         &
         index_of_process_orographic_gwd, index_of_process_nonorographic_gwd,          &
         ldiag3d, lssav, flag_for_gwd_generic_tend, do_ugwp_v0, do_ugwp_v0_orog_only,  &
         do_ugwp_v0_nst_only, do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd,     &
         gwd_opt, errmsg, errflg)

    implicit none

    ! interface variables
    integer,                 intent(in) :: me, master, im, levs, ntrac, kdt, lonr, nmtvr
    integer,                 intent(in) :: gwd_opt
    integer,                 intent(in), dimension(im)       :: kpbl
    real(kind=kind_phys),    intent(in), dimension(im)       :: oro, oro_uf, hprime, oc, theta, sigma, gamma
    real(kind=kind_phys),    intent(in), dimension(im)       :: varss,oc1ss, dx

!vay-nov 2020
    real(kind=kind_phys),    intent(in), dimension(im,4)     ::  oa4ss,ol4ss   
    
    logical,                 intent(in)                      :: flag_for_gwd_generic_tend
    
    ! elvmax is intent(in) for CIRES UGWPv1, but intent(inout) for GFS GWDPS
    
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
      &         dtaux2d_bl(:,:),dtauy2d_bl(:,:),                  &
      &         dtaux2d_ss(:,:),dtauy2d_ss(:,:),                  &
      &         dtaux2d_fd(:,:),dtauy2d_fd(:,:)

    real(kind=kind_phys), intent(in) ::     br1(im),              &
      &                                     hpbl(im),             &
      &                                     slmsk(im)

    real(kind=kind_phys),    intent(out), dimension(im)         :: dusfcg, dvsfcg
    real(kind=kind_phys),    intent(out), dimension(im)         :: zmtb, zlwb, zogw, rdxzb
    real(kind=kind_phys),    intent(out), dimension(im)         :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
    real(kind=kind_phys),    intent(out), dimension(im, levs)   :: gw_dudt, gw_dvdt, gw_dtdt, gw_kdis
    real(kind=kind_phys),    intent(out), dimension(:,:)        :: dudt_mtb, dudt_tms
    real(kind=kind_phys),    intent(out), dimension(:,:)        :: dtaux2d_ls, dtauy2d_ls

    ! The dtend array is are only allocated if ldiag=.true.
    real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
    integer, intent(in) :: dtidx(:,:), index_of_temperature, index_of_x_wind, &
         index_of_y_wind, index_of_process_nonorographic_gwd, &
         index_of_process_orographic_gwd
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
                                         do_ugwp_v0_nst_only,               &
                                         do_gsl_drag_ls_bl, do_gsl_drag_ss, &
                                         do_gsl_drag_tofd

    character(len=*),        intent(out) :: errmsg
    integer,                 intent(out) :: errflg

    ! local variables
    integer :: i, k
    real(kind=kind_phys), dimension(im)       :: sgh30
    real(kind=kind_phys), dimension(im, levs) :: Pdvdt, Pdudt
    real(kind=kind_phys), dimension(im, levs) :: Pdtdt, Pkdis
 
    real(kind=kind_phys), parameter :: tamp_mpa=30.e-3

    integer :: nmtvr_temp, idtend

    real(kind=kind_phys), dimension(:,:), allocatable :: tke
    real(kind=kind_phys), dimension(:),   allocatable :: turb_fac, tem
    real(kind=kind_phys) :: rfac, tx1

    real(kind=kind_phys) :: inv_g
    real(kind=kind_phys), dimension(im, levs)   :: zmet  ! geopotential height at model Layer centers
    real(kind=kind_phys), dimension(im, levs+1) :: zmeti ! geopotential height at model layer interfaces

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! 1) ORO stationary GWs
    !    ------------------

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
!
! put zeros due to xy GSL-drag style: dtaux2d_bl,dtauy2d_bl,dtaux2d_ss.......dusfc_ls,dvsfc_ls
!
        tau_mtb  = 0. ; tau_ogw  = 0. ; tau_tofd = 0.
        dudt_mtb = 0. ; dudt_tms = 0.
	
    end if

    if ( do_ugwp_v0.or.do_ugwp_v0_orog_only.or.do_ugwp_v0_nst_only ) then

      do k=1,levs
        do i=1,im
          Pdvdt(i,k) = 0.0
          Pdudt(i,k) = 0.0
          Pdtdt(i,k) = 0.0
          Pkdis(i,k) = 0.0
        enddo
      enddo

    end if

    if ( do_ugwp_v0.or.do_ugwp_v0_orog_only ) then

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
        idtend = dtidx(index_of_x_wind,index_of_process_orographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + Pdudt*dtp
        endif
        
        idtend = dtidx(index_of_y_wind,index_of_process_orographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + Pdvdt*dtp
        endif

        idtend = dtidx(index_of_temperature,index_of_process_orographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + Pdtdt*dtp
        endif
      endif
   
    end if 



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin non-stationary GW schemes
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! ugwp_v0 non-stationary GW drag
    !
    if (do_ugwp_v0.or.do_ugwp_v0_nst_only) then

      if (cdmbgwd(3) > 0.0) then

        ! 2) non-stationary GW-scheme with GMAO/MERRA GW-forcing
        call slat_geos5_tamp_v0(im, tamp_mpa, xlat_d, tau_ngw)

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
            gw_dtdt(i,k) = gw_dtdt(i,k)+ Pdtdt(i,k)
            gw_dudt(i,k) = gw_dudt(i,k)+ Pdudt(i,k)
            gw_dvdt(i,k) = gw_dvdt(i,k)+ Pdvdt(i,k)
            gw_kdis(i,k) = gw_kdis(i,k)+ Pkdis(i,k)
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
 
      if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
        idtend = dtidx(index_of_x_wind,index_of_process_nonorographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + Pdudt*dtp
        endif
        
        idtend = dtidx(index_of_y_wind,index_of_process_nonorographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + Pdvdt*dtp
        endif

        idtend = dtidx(index_of_temperature,index_of_process_nonorographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + Pdtdt*dtp
        endif
      endif

    end if  ! do_ugwp_v0.or.do_ugwp_v0_nst_only 


    end subroutine unified_ugwp_run
!! @}
!>@}
end module unified_ugwp

!>  \file cires_ugwp.F90
!! This file contains the Unified Gravity Wave Physics (UGWP) scheme by Valery Yudin (University of Colorado, CIRES)
!! See Valery Yudin's presentation at 2017 NGGPS PI meeting:
!! Gravity waves (GWs): Mesoscale GWs transport momentum, energy (heat) , and create eddy mixing in the whole atmosphere domain; Breaking and dissipating GWs deposit: (a) momentum; (b) heat (energy); and create (c) turbulent mixing of momentum, heat, and tracers
!! To properly incorporate GW effects (a-c) unresolved by DYCOREs we need GW physics
!! "Unified": a) all GW effects due to both dissipation/breaking; b) identical GW solvers for all GW sources; c) ability to replace solvers.
!! Unified Formalism:
!! 1. GW Sources: Stochastic and physics based mechanisms for GW-excitations in the lower atmosphere, calibrated by the high-res analyses/forecasts, and observations (3 types of GW sources: orography, convection, fronts/jets).
!! 2. GW Propagation: Unified solver for “propagation, dissipation and breaking” excited from all type of GW sources.
!! 3. GW Effects: Unified representation of GW impacts on the ‘resolved’ flow for all sources (energy-balanced schemes for momentum, heat and mixing).
!! https://www.weather.gov/media/sti/nggps/Presentations%202017/02%20NGGPS_VYUDIN_2017_.pdf

module cires_ugwp

    use machine, only: kind_phys

    use cires_ugwp_module, only: knob_ugwp_version, cires_ugwp_mod_init, cires_ugwp_mod_finalize

    implicit none

    private

    public cires_ugwp_init, cires_ugwp_run, cires_ugwp_finalize

    logical :: is_initialized = .False.

contains

! ------------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme v0
! ------------------------------------------------------------------------
!>@brief The subroutine initializes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_init Argument Table
!! \htmlinclude cires_ugwp_init.html
!!
#endif
! -----------------------------------------------------------------------
!
    subroutine cires_ugwp_init (me, master, nlunit, logunit, fn_nml2, &
                lonr, latr, levs, ak, bk, dtp, cdmvgwd, cgwf,         &
                pa_rf_in, tau_rf_in, con_p0, errmsg, errflg)

!----  initialization of cires_ugwp
    implicit none

    integer,              intent (in) :: me
    integer,              intent (in) :: master
    integer,              intent (in) :: nlunit
    integer,              intent (in) :: logunit
    integer,              intent (in) :: lonr
    integer,              intent (in) :: levs
    integer,              intent (in) :: latr
    real(kind=kind_phys), intent (in) :: ak(:), bk(:)
    real(kind=kind_phys), intent (in) :: dtp
    real(kind=kind_phys), intent (in) :: cdmvgwd(2), cgwf(2) ! "scaling" controls for "old" GFS-GW schemes
    real(kind=kind_phys), intent (in) :: pa_rf_in, tau_rf_in
    real(kind=kind_phys), intent (in) :: con_p0

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

    if (is_initialized) return

    call cires_ugwp_mod_init (me, master, nlunit, logunit, fn_nml2,  &
                              lonr, latr, levs, ak, bk, con_p0, dtp, &
                              cdmvgwd, cgwf, pa_rf_in, tau_rf_in)

    if (.not.knob_ugwp_version==0) then
       write(errmsg,'(*(a))') 'Logic error: CCPP only supports version zero of UGWP'
       errflg = 1
       return
    end if

    is_initialized = .true.

    end subroutine cires_ugwp_init


! -----------------------------------------------------------------------
! finalize of cires_ugwp   (_finalize)
! -----------------------------------------------------------------------

!>@brief The subroutine finalizes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_finalize Argument Table
!! \htmlinclude cires_ugwp_finalize.html
!!
#endif
    subroutine cires_ugwp_finalize(errmsg, errflg)

    implicit none
!
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

    call cires_ugwp_mod_finalize()

    is_initialized = .false.

    end subroutine cires_ugwp_finalize


! -----------------------------------------------------------------------
!    originally from ugwp_driver_v0.f
!    driver of cires_ugwp   (_driver)
! -----------------------------------------------------------------------
!   driver is called after pbl & before chem-parameterizations
! -----------------------------------------------------------------------
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
! -----------------------------------------------------------------------
!>@brief The subroutine executes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_run Argument Table
!! \htmlinclude cires_ugwp_run.html
!!
#endif

! subroutines original
     subroutine cires_ugwp_run(do_ugwp, me,  master, im,  levs, ntrac, dtp, kdt, lonr, &
         oro, oro_uf, hprime, nmtvr, oc, theta, sigma, gamma, elvmax, clx, oa4,        &
         do_tofd,  cdmbgwd, xlat, xlat_d, sinlat, coslat, area,                        &
         ugrs, vgrs, tgrs, qgrs, prsi, prsl, prslk, phii, phil,                        &
         del, kpbl, dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,                &
         tau_tofd, tau_mtb, tau_ogw, tau_ngw, zmtb, zlwb, zogw,                        &
         dudt_mtb,dudt_ogw, dudt_tms, dudt, dvdt, dtdt, rdxzb,                         &
         con_g, con_pi, con_cp, con_rd, con_rv, con_fvirt, errmsg, errflg)

    implicit none

    ! interface variables
    integer,                 intent(in) :: me, master, im, levs, ntrac, kdt, lonr, nmtvr
    integer,                 intent(in), dimension(im)       :: kpbl
    real(kind=kind_phys),    intent(in), dimension(im)       :: oro, oro_uf, hprime, oc, theta, sigma, gamma, elvmax
    real(kind=kind_phys),    intent(in), dimension(im, 4)    :: clx, oa4
    real(kind=kind_phys),    intent(in), dimension(im)       :: xlat, xlat_d, sinlat, coslat, area
    real(kind=kind_phys),    intent(in), dimension(im, levs) :: del, ugrs, vgrs, tgrs, prsl, prslk, phil
    real(kind=kind_phys),    intent(in), dimension(im, levs+1) :: prsi, phii
    real(kind=kind_phys),    intent(in), dimension(im, levs, ntrac):: qgrs
    real(kind=kind_phys),    intent(in) :: dtp, cdmbgwd(2)
    logical,                 intent(in) :: do_ugwp, do_tofd

    real(kind=kind_phys),    intent(out), dimension(im)      :: dusfcg, dvsfcg
    real(kind=kind_phys),    intent(out), dimension(im)      :: zmtb, zlwb, zogw, rdxzb
    real(kind=kind_phys),    intent(out), dimension(im)      :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
    real(kind=kind_phys),    intent(out), dimension(im, levs):: gw_dudt, gw_dvdt, gw_dtdt, gw_kdis
    real(kind=kind_phys),    intent(out), dimension(im, levs):: dudt_mtb, dudt_ogw, dudt_tms

    real(kind=kind_phys),    intent(inout), dimension(im, levs):: dudt, dvdt, dtdt

    real(kind=kind_phys),    intent(in) :: con_g, con_pi, con_cp, con_rd, con_rv, con_fvirt

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

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! wrap everything in a do_ugwp 'if test' in order not to break the namelist functionality
    if (do_ugwp) then

        ! topo paras
        ! w/ orographic effects
        if(nmtvr == 14)then
            ! calculate sgh30 for TOFD
            sgh30 = abs(oro - oro_uf)
        ! w/o orographic effects
        else
            sgh30   = 0.
        endif

        zlwb(:)   = 0.

        call GWDPS_V0(im, levs, lonr, do_tofd, Pdvdt, Pdudt, Pdtdt, Pkdis,          &
             ugrs, vgrs, tgrs, qgrs(:,:,1), kpbl, prsi,del,prsl, prslk, phii, phil, &
             dtp, kdt, sgh30, hprime, oc, oa4, clx, theta, sigma, gamma, elvmax,    &
             dusfcg, dvsfcg, xlat_d, sinlat, coslat, area, cdmbgwd,          &
             me, master, rdxzb, zmtb, zogw, tau_mtb, tau_ogw, tau_tofd, dudt_mtb, dudt_ogw, dudt_tms)


        ! 1) non-stationary GW-scheme with GMAO/MERRA GW-forcing
        call slat_geos5_tamp(im, tamp_mpa, xlat_d, tau_ngw)


        ! 2) non-stationary GW-scheme with GEOS-5/MERRA GW-forcing
        call fv3_ugwp_solv2_v0(im, levs, dtp, tgrs, ugrs, vgrs,qgrs(:,:,1), &
             prsl, prsi, phil, xlat_d, sinlat, coslat, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis, &
             tau_ngw, me, master, kdt)

        if(pogw /= 0.)then

            do k=1,levs
            do i=1,im
                gw_dtdt(i,k) = pngw*gw_dtdt(i,k)+ pogw*Pdtdt(i,k)
                gw_dudt(i,k) = pngw*gw_dudt(i,k)+ pogw*Pdudt(i,k)
                gw_dvdt(i,k) = pngw*gw_dvdt(i,k)+ pogw*Pdvdt(i,k)
                gw_kdis(i,k) = pngw*gw_kdis(i,k)+ pogw*Pkdis(i,k)

                ! accumulation of tendencies for CCPP to replicate EMC-physics updates (!! removed in latest code commit to VLAB)
                dudt(i,k) = dudt(i,k) +gw_dudt(i,k)
                dvdt(i,k) = dvdt(i,k) +gw_dvdt(i,k)
                dtdt(i,k) = dtdt(i,k) +gw_dtdt(i,k)
            enddo
            enddo

        else

            tau_mtb = 0.  ; tau_ogw =0.;     tau_tofd =0.
            dudt_mtb =0. ; dudt_ogw = 0.;  dudt_tms=0.

        endif

        return


        !=============================================================================
        ! make "ugwp eddy-diffusion" update for gw_dtdt/gw_dudt/gw_dvdt by solving
        ! vert diffusion equations & update "Statein%tgrs, Statein%ugrs, Statein%vgrs"
        !=============================================================================
        ! 3) application of "eddy"-diffusion to "smooth" UGWP-related tendencies
        !------------------------------------------------------------------------------
        ed_dudt(:,:) =0.; ed_dvdt(:,:) = 0. ; ed_dtdt(:,:) = 0.

        call edmix_ugwp_v0(im, levs, dtp, tgrs, ugrs, vgrs, qgrs(:,:,1), &
             del, prsl, prsi, phil, prslk, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis, &
             ed_dudt, ed_dvdt, ed_dtdt, me, master, kdt)
        gw_dtdt = gw_dtdt*(1.-pked) +  ed_dtdt*pked
        gw_dvdt = gw_dvdt*(1.-pked) +  ed_dvdt*pked
        gw_dudt = gw_dudt*(1.-pked) +  ed_dudt*pked



    endif ! do_ugwp

    end subroutine cires_ugwp_run


end module cires_ugwp

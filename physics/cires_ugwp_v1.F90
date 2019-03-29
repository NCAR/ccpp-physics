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

    use cires_ugwp_module, only: cires_ugwp_mod_init, cires_ugwp_mod_finalize, cires_ugwp_driver  


contains
! -----------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme
! -----------------------------------------------------------------------
!>@brief The subroutine initializes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_init Argument Table
!! | local_name       | standard_name                                                                 | long_name                                               | units  | rank | type      | kind      | intent | optional |
!! |------------------|-------------------------------------------------------------------------------|---------------------------------------------------------|--------|------|-----------|-----------|--------|----------|
!! | me               | mpi_rank                                                                      | MPI rank of current process                             | index  | 0    | integer   |           | in     | F        |
!! | master           | mpi_root                                                                      | MPI rank of master process                              | index  | 0    | integer   |           | in     | F        |
!! | nlunit           | iounit_namelist                                                               | fortran unit number for opening namelist file           | none   | 0    | integer   |           | in     | F        |
!! | logunit          | iounit_log                                                                    | fortran unit number for writing logfile                 | none   | 0    | integer   |           | in     | F        |
!! | fn_nml2          | namelist_filename                                                             | namelist filename for ugwp                              | none   | 0    | character | len=*     | in     | F        |
!! | lonr             | number_of_equatorial_longitude_points                                         | number of global points in x-dir (i) along the equator  | count  | 0    | integer   |           | in     | F        |
!! | latr             | number_of_laitude_points                                                      | number of global points in y-dir (j) along any meridian | count  | 0    | integer   |           | in     | F        |
!! | levs             | vertical_dimension                                                            | number of vertical levels                               | count  | 0    | integer   |           | in     | F        |
!! | ak               | a_parameter_of_the_hybrid_coordinate                                          | components of the hybrid coordinate                     | Pa     | 1    | real      | kind_phys | in     | F        |
!! | bk               | b_parameter_of_the_hybrid_coordinate                                          | components of the hybrid coordinate                     | none   | 1    | real      | kind_phys | in     | F        |
!! | pref             | n/a for pref                                                                  | reference pressure from radlw_datatb.f                  |?Pa/hPa?| 2    | real      | kind_phys | in     | F        |
!! | dtp              | time_step_for_physics                                                         | physics timestep                                        | s      | 0    | real      | kind_phys | in     | F        |
!! | cdmvgwd          | multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag | multiplication factors for cdmb and gwd                 | none   | 1    | real      | kind_phys | in     | F        |
!! | cgwf             | multiplication_factors_for_convective_gravity_wave_drag                       | multiplication factor for convective GWD                | none   | 1    | real      | kind_phys | in     | F        | 
!! | do_shoc          | flag_for_shoc                                                                 | flag to indicate use of SHOC                            | flag   | 0    | logical   |           | in     | F        |
!! | errmsg           | ccpp_error_message                                                            | error message for error handling in CCPP                | none   | 0    | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                               | error flag for error handling in CCPP                   | flag   | 0    | integer   |           | out    | F        |
!!
#endif
! -----------------------------------------------------------------------
!
! init  of cires_ugwp   (_init)  called from GFS_driver.F90
!
! -----------------------------------------------------------------------
    subroutine cires_ugwp_init (me, master, nlunit, logunit, fn_nml2, &
              lonr, latr, levs, ak, bk, pref, dtp, cdmvgwd, cgwf,  do_shoc, &
              errmsg, errflg)

    implicit none
    
    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit
    integer, intent (in) :: lonr
    integer, intent (in) :: levs
    integer, intent (in) :: latr   
    real,    intent (in) :: ak(levs+1), bk(levs+1), pref
    real,    intent (in) :: dtp   
    real,    intent (in) :: cdmvgwd(2), cgwf(2)             ! "scaling" controls for "old" GFS-GW schemes     
         
    character(len=*), intent (in) :: fn_nml2
    !character(len=*), parameter   :: fn_nml='input.nml' 

    logical,          intent( in) :: do_shoc       

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


    if (do_shoc) then
        write(errmsg,'(*(a))') 'SHOC is not currently compatible with CIRES UGWP'
        errflg = 1
        return
    endif

    call cires_ugwp_mod_init (me, master, nlunit, logunit, fn_nml2, &
                              lonr, latr, levs, ak, bk, pref, dtp, cdmvgwd, cgwf) 

    is_initialized = .true.


    end subroutine cires_ugwp_init




! -----------------------------------------------------------------------
! finalize of cires_ugwp   (_finalize) 
! -----------------------------------------------------------------------

!>@brief The subroutine finalizes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_finalize Argument Table
!! | local_name       | standard_name                                                                 | long_name                                               | units  | rank | type      | kind      | intent | optional |
!! |------------------|-----------|--------|------|-----------|-----------|--------|----------|-------|
!! | errmsg           | ccpp_error_message                                                            | error message for error handling in CCPP                | none   | 0    | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                               | error flag for error handling in CCPP                   | flag   | 0    | integer   |           | out    | F        |
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
! run of cires_ugwp(_run)  called from GFS_physics_driver.F90
! -----------------------------------------------------------------------
!    originally from cires_ugwp_module.F90
!    driver of cires_ugwp   (_driver)
! -----------------------------------------------------------------------
!   driver is called after pbl & before chem-parameterizations
! -----------------------------------------------------------------------
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
! -----------------------------------------------------------------------
!>@brief The subroutine executes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_run Argument Table
!! | local_name       | standard_name                                                                  | long_name                                                    | units     | rank | type      | kind      | intent | optional |
!! |------------------|--------------------------------------------------------------------------------|--------------------------------------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | me               | mpi_rank                                                                       | MPI rank of current process                                  | index     | 0    | integer   |           | in     | F        |
!! | master           | mpi_root                                                                       | MPI rank of master process                                   | index     | 0    | integer   |           | in     | F        |
!! | im               | horizontal_loop_extent                                                         | horizontal                                                   | count     | 0    | integer   |           | in     | F        |
!! | levs             | vertical_dimension                                                             | number of vertical levels                                    | count     | 0    | integer   |           | in     | F        |
!! | ntrac            | number_of_tracers                                                              | number of tracers                                            | count     | 0    | integer   |           | in     | F        |
!! | dtp              | time_step_for_physics                                                          | physics timestep                                             | s         | 0    | real      | kind_phys | in     | F        |
!! | kdt              | index_of_time_step                                                             | current forecast iteration                                   | index     | 0    | integer   |           | in     | F        |
!! | nvoro            | number_of_statistical_measures_of_subgrid_orography                            | number of topographic variables in GWD                       | count     | 0    | integer   |           | in     | F        |
!! | nmtvr            | number_of_statistical_measures_of_subgrid_orography                            | number of topographic variables in GWD                       | count     | 0    | integer   |           | in     | F        |
!! | lprnt            | flag_print                                                                     | control flag for diagnostic print out                        | flag      | 0    | logical   |           | in     | F        |
!! | lonr             | number_of_equatorial_longitude_points                                          | number of global points in x-dir (i) along the equator       | count     | 0    | integer   |           | in     | F        |
!! | pa_rf            | pressure_cutoff_for_rayleigh_damping                                           | pressure level from which Rayleigh Damping is applied        | Pa        | 0    | real      | kind_phys | in     | F        |
!! | tau_rf           | time_scale_for_rayleigh_damping                                                | time scale for Rayleigh damping in days                      | d         | 0    | real      | kind_phys | in     | F        |  
!! | cdmbgwd          | multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag  | multiplication factors for cdmb and gwd                      | none      | 1    | real      | kind_phys | in     | F        |
!! | xlat             | latitude                                                                       | grid latitude in radians                                     | radians   | 1    | real      | kind_phys | in     | F        | 
!! | sinlat           | sine_of_latitude                                                               | sine of the grid latitude                                    | none      | 1    | real      | kind_phys | in     | F        |
!! | coslat           | cosine_of_latitude                                                             | cosine of the grid latitude                                  | none      | 1    | real      | kind_phys | in     | F        |  
!! | Statein          | GFS_statein_type_instance                                                      | Fortran DDT containing prog state data in from dycore        | DDT       | 0    | GFS_statein_type |    | in     | F        |
!! | delp             | layer_pressure_thickness_for_radiation                                         | layer pressure thickness on radiation levels                 | hPa       | 2    | real      | kind_phys | in     | F        | 
!! | orostat          | statistical_measures_of_subgrid_orography                                      | orographic metrics                                           | various   |    2 | real      | kind_phys | in     | F        |
!! | kpbl             | vertical_index_at_top_of_atmosphere_boundary_layer                             | vertical index at top atmospheric boundary layer             | index     | 1    | integer   |           | in     | F        |
!! | dusfc            | instantaneous_x_stress_due_to_gravity_wave_drag                                | zonal surface stress due to orographic gravity wave drag     | Pa        | 1    | real      | kind_phys | in     | F        |
!! | dvsfc            | instantaneous_y_stress_due_to_gravity_wave_drag                                | meridional surface stress due to orographic gravity wave drag| Pa        | 1    | real      | kind_phys | in     | F        |
!! | errmsg           | ccpp_error_message                                                            | error message for error handling in CCPP                | none   | 0    | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                               | error flag for error handling in CCPP                   | flag   | 0    | integer   |           | out    | F        |
!!
#endif
   subroutine cires_ugwp_run                                     &
       (me, master, im, levs, ntrac,dtp, kdt, nvoro, nmtvr, lprnt, lonr,& 
        pa_rf, tau_rf, cdmbgwd,  xlat, sinlat,  coslat,        &
       Statein, delp,  orostat, kpbl,                                & 
       dusfc, dvsfc,                         &
! diagnostics (WL: do we need them?)
!      dudt,  dvdt, dtdt, kdis, axtot, axo, axc, axf,  aytot, ayo, ayc, ayf,                  &
!       eps_tot,  ekdis,  trig_okw, trig_fgf,                         &
!       dcheat, precip, cld_klevs, zmtb, scheat, dlength, cldf,       & 
! COORDE-2019 diagnostics    without 3d-fluxes:  tauz_ogw, tauz_ngw ....   (WL: do we need them?) 
!       taus_sso, taus_ogw, tauf_ogw, tauf_ngw,                       &
!       ugw_zmtb, ugw_zlwb, ugw_zogw, ugw_axmtb,ugw_axlwb, ugw_axtms, &
!
      ,errmsg, errflg ) 



               
    use machine,       only: kind_phys
    use physcons,      only: con_cp, con_fvirt, con_g, con_rd, pi
    use gfs_typedefs,  only: gfs_statein_type
    use ugwp_common,   only: omega2 
    use ugwp_okw_init,  only : &
        eff_okw, nstokw, nwokw, ch_okwp, nazokw, spf_okwp, xaz_okwp, yaz_okwp
    use ugwp_conv_init, only : &
        eff_con, nstcon, nwcon, ch_conv, nazcon, spf_conv, xaz_conv, yaz_conv
    use ugwp_fjet_init, only : &
        eff_fj,  nstfj,  nwfj,  ch_fjet, nazfj,  spf_fjet, xaz_fjet, yaz_fjet         

    implicit none
     
    ! interface variables     
    logical, intent(in)                 :: lprnt
    integer, intent(in)                 :: me, im, levs, ntrac, kdt, lonr, nvoro, nmtvr
    real(kind=kind_phys),  intent(in)   :: dtp
    real(kind=kind_phys),  intent(in)   :: pa_rf, tau_rf  
    real(kind=kind_phys),  intent(in)   :: cdmbgwd(2)
       
    type(GFS_statein_type),intent(in)   :: Statein

    integer, intent(in)   ::  kpbl(im)      
!    real,    intent(in)   ::  hpbl(im)     
    real(kind=kind_phys),  intent(in)   ::  orostat(im, 14)       
    real(kind=kind_phys),  intent(in)   ::  delp(im, levs)     

    real(kind=kind_phys),  intent(out),  dimension(im)       ::  xlat, xlatd, sinlat, coslat  
! diagnostics (WL: do we need them?)
!    real(kind=kind_phys),  intent(out),  dimension(im, levs) ::  trig_okw, trig_fgf     
!    real(kind=kind_phys),  intent(out),  dimension(im)       ::  precip              ! precip-n rates and   
!    integer                intent(out),  dimension(im, 3)    ::  cld_klevs           ! indices fo cloud top/bot/?
!    real(kind=kind_phys),  intent(out),  dimension(im, levs) ::  dcheat,  scheat     ! deep and shal conv heat tend.   
     
!    real(kind=kind_phys),  intent(out),  dimension(im)       ::  dlength           ! tail-grid  box scale in meters
!    real(kind=kind_phys),  intent(out),  dimension(im)       ::  cldf              ! "bizzard" old cgwd-tuning knobs dimensionless         
!===================
! tendency + kdis
!===================    
!    real(kind=kind_phys),  intent(out),  dimension(im, levs) ::  dudt,  dvdt, dtdt, kdis 
!    real(kind=kind_phys),  intent(out),  dimension(im, levs) ::  axtot, axo, axc, axf  
!    real(kind=kind_phys),  intent(out),  dimension(im, levs) ::  aytot, ayo, ayc, ayf  
!    real(kind=kind_phys),  intent(out),   dimension(im, levs) ::  eps_tot,  ekdis 
     
!    real(kind=kind_phys), dimension(im, levs) ::     eds_o, kdis_o
!    real(kind=kind_phys), dimension(im, levs) ::     eds_c, kdis_c
!    real(kind=kind_phys), dimension(im, levs) ::     eds_f, kdis_f                     
!    real(kind=kind_phys), dimension(im, levs) ::  ax_rf, ay_rf, eps_rf                     
!==================================================================================
! diagnostics for OGW & NGW  +  SSO effects axmtb, axlwb, axtms
!==================================================================================
!    real(kind=kind_phys),  intent(out), dimension(im)       ::   dusfc, dvsfc 
!    real(kind=kind_phys),  intent(out), dimension(im)       ::   taus_sso, taus_ogw, tauf_ogw, tauf_ngw
!    real(kind=kind_phys),  intent(out), dimension(im)       ::   ugw_zmtb, ugw_zlwb, ugw_zogw
!    real(kind=kind_phys),  intent(out), dimension(im, levs) ::   ugw_axmtb,ugw_axlwb, ugw_axtms
    !real(kind=kind_phys),  intent(out), dimension(im, levs) ::   tauz_ogw, tauz_ngw,  wtauz

! 
!    knob_ugwp_source=[ 1,    1,    1,       0 ]
!                       oro conv    nst   imbal-okw
!    local variables
!
!    integer :: i, j, k, istype, ido
!    
! internal diagnostics for oro-waves, lee waves, and mtb :  
!  
!    real(kind_phys), dimension(im) :: dusfc_mb, dvsfc_mb, dusfc_ogw, dvsfc_ogw
!    real(kind_phys), dimension(im) :: dusfc_lwb, dvsfc_lwb
!    real(kind_phys), dimension(im) :: zmtb, zlwb, zogw      ! GW-launch levels in "meters"   
!
!    real(kind_phys), dimension(im) ::  fcor,  c2f2
!
! three sources with different: a) spectra-content/azimuth; b) efficiency ;c) spectral shape
!     
!    real(kind_phys), dimension(im) ::  taub_con,  taub_fj, taub_okw
!    integer        , dimension(im) ::  klev_okw,  klev_fj, klev_con 
!    integer        , dimension(im) ::  if_okw,    if_con,  if_fj
!    integer                        ::  nf_okw,    nf_con,  nf_fj   

    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg


! local variables
    integer eps
    real del

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0


!    if (knob_ugwp_version == 1 ) then
!        if (kdt < 2  .and.  me == master) then
!            print *, ' VAY-attention UGWP-V1 cires_ugwp_driver '
!            print *, ' Only Test-mode by developers '
!            !stop ' cires_ugwp_driver Test-mode Jan 2019 '
!        endif  


    call cires_ugwp_driver                                        &
      (im, levs, ntrac,dtp, kdt, me, lprnt,   lonr,                 &
       pa_rf, tau_rf, cdmbgwd,  xlat, xlatd, sinlat,  coslat,        &
       Statein, delp,  orostat, kpbl,                                &
       dusfc, dvsfc, dudt,  dvdt, dtdt, kdis)
! diagnostics (WL: do we need them?)
!       axtot, axo, axc, axf,  aytot, ayo, ayc, ayf,                  &
!       eps_tot,  ekdis,  trig_okw, trig_fgf,                         &
!       dcheat, precip, cld_klevs, zmtb, scheat, dlength, cldf,       &
! COORDE-2019 diagnostics without 3d-fluxes:  tauz_ogw, tauz_ngw ....   (WL: do we need them?) 
!       taus_sso, taus_ogw, tauf_ogw, tauf_ngw,                       &
!       ugw_zmtb, ugw_zlwb, ugw_zogw, ugw_axmtb,ugw_axlwb, ugw_axtms, &
!       knob_ugwp_version)

! so far, not calculate tendencies for simplicity
!    do k=1,levs
!        do i=1,im
!            Pdtdt(i,k) = gw_dtdt(i,k)
!            Pdudt(i,k) = gw_dudt(i,k)
!            Pdvdt(i,k) = gw_dvdt(i,k)
!        enddo
!    enddo
        
!    else
!
!    if (kdt < 2 .and.  me == master) then
!        print *, ' VAY-attention UGWP-V0, Jan 2019 '
!    endif
!
!      eps = 1.0E-5    ! eps = the absolute accuracy requirment for the CDF
!      del = 2.0*eps
!    call cires_ugwp_driver_v0                                      &
!       (me, master,im, levs, ntrac, nvoro, dtp, kdt, lonr,    &
!       nmtvr,  do_tofd,cdmbgwd, xlat, xlatd, sinlat, coslat,           &
!       Statein, del, oro_stat, sgh30, kpbl,                         &
!       dusfc, dvsfc, dudt,  dvdt, dtdt, kdis,         &

!! WL: revisit from here 
!       tau_tms, tau_mtb,  tau_ogw,  tau_ngw,                        &
!       zm_mtb, zm_lwb, zm_ogw,  ax_mtb, ax_ogw, ax_tms,             &
!       Diag%zmtnblck )  
!cires_ugwp_driver_v0
!       (me,  master,im,  levs, ntrac, nvoro, dtp, kdt, imx,
!     &    nmtvr,  do_tofd,  cdmbgwd,  xlat, xlatd, sinlat, coslat,
!     &    Statein,  del, oro_stat, sgh30, kpbl,
!     &    dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,
!     &    tau_tofd, tau_mtb, tau_ogw, tau_ngw,
!     &    zmtb, zlwb, zogw, du3dt_mtb,du3dt_ogw, du3dt_tms,rdxzb )
!    if ( me == master)  print *, ' ugwp time-step=', kdt

    ! in case of errors, set errflg to a value != 0,
    ! create a meaningfull error message and return


    return
    end subroutine cires_ugwp_run




end module cires_ugwp

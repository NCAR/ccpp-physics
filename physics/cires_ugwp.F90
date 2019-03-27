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

    use cires_ugwp_module, only: cires_ugwp_mod_init, cires_ugwp_mod_finalize, cires_ugwp_driver, GWDPS_V0, slat_geos5_tamp, fv3_ugwp_solv2_v0 


contains
! ------------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme v0
! ------------------------------------------------------------------------
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
!! | dtp              | time_step_for_physics                                                         | physics timestep                                        | s      | 0    | real      | kind_phys | in     | F        |
!! | cdmvgwd          | multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag | multiplication factors for cdmb and gwd                 | none   | 1    | real      | kind_phys | in     | F        |
!! | cgwf             | multiplication_factors_for_convective_gravity_wave_drag                       | multiplication factor for convective GWD                | none   | 1    | real      | kind_phys | in     | F        | 
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
              lonr, latr, levs, ak, bk, dtp, cdmvgwd, cgwf,  &
              errmsg, errflg)

!----  initialization of cires_ugwp 
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
!! | local_name       | standard_name                                                                  | long_name                                                    | units     | rank | type      | kind      | intent | optional |
!! |------------------|--------------------------------------------------------------------------------|--------------------------------------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | Model            | GFS_control_type_instance                                                      | Fortran DDT containing FV3-GFS model control parameters      | DDT       | 0    | GFS_control_type  |   | in     | F        |
!! | me               | mpi_rank                                                                       | MPI rank of current process                                  | index     | 0    | integer   |           | in     | F        |
!! | master           | mpi_root                                                                       | MPI rank of master process                                   | index     | 0    | integer   |           | in     | F        |
!! | im               | horizontal_loop_extent                                                         | horizontal                                                   | count     | 0    | integer   |           | in     | F        |
!! | levs             | vertical_dimension                                                             | number of vertical levels                                    | count     | 0    | integer   |           | in     | F        |
!! | ntrac            | number_of_tracers                                                              | number of tracers                                            | count     | 0    | integer   |           | in     | F        |
!! | nvoro            | number_of_statistical_measures_of_subgrid_orography                            | number of topographic variables in GWD                       | count     | 0    | integer   |           | in     | F        |
!! | dtp              | time_step_for_physics                                                          | physics timestep                                             | s         | 0    | real      | kind_phys | in     | F        |
!! | kdt              | index_of_time_step                                                             | current forecast iteration                                   | index     | 0    | integer   |           | in     | F        |
!! | lonr             | number_of_equatorial_longitude_points                                          | number of global points in x-dir (i) along the equator       | count     | 0    | integer   |           | in     | F        |
!! | nmtvr            | number_of_statistical_measures_of_subgrid_orography                            | number of topographic variables in GWD                       | count     | 0    | integer   |           | in     | F        |
!! | do_tofd          | turb_oro_form_drag_flag                                                        | flag for turbulent orographic form drag                      | flag      | 0    | logical   |           | none   | F        | 
!! | cdmbgwd          | multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag  | multiplication factors for cdmb and gwd                      | none      | 1    | real      | kind_phys | in     | F        |
!! | xlat             | latitude                                                                       | grid latitude in radians                                     | radians   | 1    | real      | kind_phys | in     | F        | 
!! | xlat_d           | latitude_degree                                                                | latitude in degrees                                          | degree    | 1    | real      | kind_phys | in     | F        |
!! | sinlat           | sine_of_latitude                                                               | sine of the grid latitude                                    | none      | 1    | real      | kind_phys | in     | F        |
!! | coslat           | cosine_of_latitude                                                             | cosine of the grid latitude                                  | none      | 1    | real      | kind_phys | in     | F        |  
!! | Statein          | GFS_statein_type                                                               | definition of type GFS_statein_type                          | DDT       | 0    | GFS_statein_type|     | none   | F        |
!! | Sfcprop          | GFS_sfcprop_type                                                               | definition of type GFS_sfcprop_type                          | DDT       | 0    | GFS_sfcprop_type|     | none   | F        | 
!! | del              | air_pressure_difference_between_midlayers                                      | air pressure difference between midlayers                    | Pa        | 2    | real      | kind_phys | none   | F        |
!! | orostat          | statistical_measures_of_subgrid_orography                                      | orographic metrics                                           | various   | 2    | real      | kind_phys | in     | F        |
!! | kpbl             | vertical_index_at_top_of_atmosphere_boundary_layer                             | vertical index at top atmospheric boundary layer             | index     | 1    | integer   |           | in     | F        |
!! | dusfcg           | instantaneous_x_stress_due_to_gravity_wave_drag                                | zonal surface stress due to orographic gravity wave drag     | Pa        | 1    | real      | kind_phys | none   | F        |
!! | dvsfcg           | instantaneous_y_stress_due_to_gravity_wave_drag                                | meridional surface stress due to orographic gravity wave drag| Pa        | 1    | real      | kind_phys | none   | F        |
!! | gw_dudt          | tendency_of_x_wind_due_to_ugwp                                                 | zonal wind tendency due to UGWP                              | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | gw_dvdt          | tendency_of_y_wind_due_to_ugwp                                                 | meridional wind tendency due to UGWP                         | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | gw_dtdt          | tendency_of_air_temperature_due_to_ugwp                                        | air temperature tendency due to UGWP                         | K s-1     | 2    | real      | kind_phys | none   | F        |
!! | gw_kdis          | eddy_mixing_due_to_ugwp                                                        | eddy mixing due to UGWP                                      | m2 s-1    | 2    | real      | kind_phys | none   | F        |
!! | tau_tofd         | instantaneous_momentum_flux_due_to_turbulent_orographic_form_drag              | momentum flux or stress due to TOFD                          | Pa        | 2    | real      | kind_phys | none   | F        | 
!! | tau_mtb          | instantaneous_momentum_flux_due_to_mountain_blocking_drag                      | momentum flux or stress due to mountain blocking drag        | Pa        | 2    | real      | kind_phys | none   | F        | 
!! | tau_ogw          | instantaneous_momentum_flux_due_to_orographic_gravity_wave_drag                | momentum flux or stress due to orographic gravity wave drag  | Pa        | 2    | real      | kind_phys | none   | F        | 
!! | tau_ngw          | instantaneous_momentum_flux_due_to_nonstationary_gravity_wave                  | momentum flux or stress due to nonstationary gravity waves   | Pa        | 2    | real      | kind_phys | none   | F        | 
!! | zmtb             | height_of_mountain_blocking                                                    | height of mountain blocking drag                             | m         | 1    | real      | kind_phys | none   | F        |
!! | zlwb             | height_of_low_level_wave_breaking                                              | height of low level wave breaking                            | m         | 1    | real      | kind_phys | none   | F        |   
!! | zogw             | height_of_launch_level_of_orographic_gravity_wave                              | height of launch level of orographic gravity wave            | m         | 1    | real      | kind_phys | none   | F        | 
!! | du3dt_mtb        | instantaneous_change_in_x_wind_due_to_mountain_blocking_drag                   | instantaneous change in x wind due to mountain blocking drag | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | du3dt_ogw        | instantaneous_change_in_x_wind_due_to_orographic_gravity_wave_drag             | instantaneous change in x wind due to orographic gw drag     | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | du3dt_tms        | instantaneous_change_in_x_wind_due_to_turbulent_orographic_form_drag           | instantaneous change in x wind due to TOFD                   | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | rdxzb            | level_of_dividing_streamline                                                   | level of the dividing streamline                             | none      | 1    | real      | kind_phys | none   | F        |
!! | errmsg           | ccpp_error_message                                                             | error message for error handling in CCPP                     | none      | 0    | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                                | error flag for error handling in CCPP                        | flag      | 0    | integer   |           | out    | F        |
!!
#endif


! subroutines original
     subroutine cires_ugwp_run(Model, me,  master,
     &    im,  levs, ntrac, nvoro, dtp, kdt, lonr,
     &    nmtvr,  do_tofd,  cdmbgwd,  xlat, xlat_d, sinlat, coslat,
     &    Statein,  Sfcprop, del, orostat, kpbl,
     &    dusfcg, dvsfcg, 
! diag (tendencies due to ugwp) 
     &    gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,
! COORDE diag
     &    tau_tofd, tau_mtb, tau_ogw, tau_ngw,
     &    zmtb, zlwb, zogw, 
     &    du3dt_mtb,du3dt_ogw, du3dt_tms,
! 
     &    rdxzb,
     &    errmsg, errflg)

               
    use machine,       only: kind_phys
    use physcons,      only: con_g, con_pi, con_cp, con_rd, con_rv, con_fvirt
    use GFS_typedefs,  only: GFS_control_type, GFS_statein_type, GFS_sfcprop_type

    implicit none

    ! interface variables     
    integer,                 intent(in) :: me, master, im, levs, ntrac, nvoro, kdt, lonr, nmtvr
    integer,                 intent(in), dimension(im)       :: kpbl
    real(kind=kind_phys),    intent(in) :: dtp, cdmbgwd(2)
    real(kind=kind_phys),    intent(in), dimension(im)       :: xlat, xlat_d, sinlat, coslat
    real(kind=kind_phys),    intent(in), dimension(im, levs) :: del
    real(kind=kind_phys),    intent(in), dimension(im, nvoro):: orostat
    logical,                 intent(in) :: do_tofd

    type(GFS_control_type),  intent(in) :: Model
    type(GFS_statein_type),  intent(in) :: Statein
    type(GFS_Sfcprop_type),  intent(in) :: Sfcprop

    real(kind=kind_phys),    intent(out), dimension(im)      :: dusfcg, dvsfcg 
    real(kind=kind_phys),    intent(out), dimension(im)      :: zmtb, zlwb, zogw, rdxzb
    real(kind=kind_phys),    intent(out), dimension(im)      :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
    real(kind=kind_phys),    intent(out), dimension(im, levs):: gw_dudt, gw_dvdt, gw_dtdt, gw_kdis
    real(kind=kind_phys),    intent(out), dimension(im, levs):: du3dt_mtb, du3dt_ogw, du3dt_tms

    character(len=*),        intent(out) :: errmsg
    integer,                 intent(out) :: errflg

    ! local variables
    integer :: i, k
    real(kind=kind_phys) :: fdaily
    real(kind=kind_phys), dimension(im)       :: hprime, oc, theta, sigma, gamm, elvmax
    real(kind=kind_phys), dimension(im, 4)    :: clx, oa4
    real(kind=kind_phys), dimension(im, levs) :: Pdvdt, Pdudt
    real(kind=kind_phys), dimension(im, levs) :: Pdtdt, Pkdis
    ! ?? what is tamp_mpa? amplitude in a unit of mPa? need to confirm with VAY
    real(kind=kind_phys), parameter :: tamp_mpa=30.e-3
    ! switches that activate impact of OGWs and NGWs (WL* how to deal with them? *WL)
    real(kind=kind_phys), parameter :: pogw=1., pngw=1., pked=1.
    real(kind=kind_phys), parameter :: ftausec = 86400. 


    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0


    ! For COORDE averaging over fdaily. It's done before 3Diag fixes and averaging ingested using "fdaily"-factor
    fdaily  = dtp/ftausec
    if (Model%fhzero .ne. 0) then
        ftausec = Model%fhzero*3600
        fdaily  =  dtp/ftausec
    else
        print *, 'UGWP: Model%fhzero = 0., Bad Averaged-diagnostics '
    endif


    ! topo paras
    ! with orographic effects
    if(nmtvr == 14)then         
        do i=1,im
            ! calculate sgh30 for TOFD
            sgh30(i) = abs(Sfcprop%oro(i) - Sfcprop%oro_uf(i))
            hprime(i) = orostat(i,1)
            oc(i)     = orostat(i,2)
            oa4(i,1)  = orostat(i,3)
            oa4(i,2)  = orostat(i,4)
            oa4(i,3)  = orostat(i,5)
            oa4(i,4)  = orostat(i,6)
            clx(i,1)  = orostat(i,7)
            clx(i,2)  = orostat(i,8)
            clx(i,3)  = orostat(i,9)
            clx(i,4)  = orostat(i,10)
            theta(i)  = orostat(i,11)
            gamm(i)   = orostat(i,12)
            sigma(i)  = orostat(i,13)
            elvmax(i) = orostat(i,14)
        enddo
    ! no orographic effects
    else
        orostat = 0.   
        sgh30   = 0.
        oc      = 0. 
        oa4     = 0.
        clx     = 0.
        theta   = 0. 
        gamm    = 0. 
        sigma   = 0.
        elvmax  = 0.
        hprime  = 0.
    endif



    zlwb(:)   = 0.

    call GWDPS_V0(im, levs, lonr, do_tofd,
     &               Pdvdt, Pdudt, Pdtdt, Pkdis,
     &  Statein%ugrs, Statein%vgrs, Statein%tgrs,                        
     &  Statein%qgrs(:,:,1),kpbl, Statein%prsi,del,Statein%prsl,                          
     &  Statein%prslk, Statein%phii, Statein%phil, dtp, kdt,
     &  sgh30, hprime, oc, oa4, clx, theta, sigma, gamm, elvmax,
     &  dusfcg, dvsfcg,
     &  nmtvr, cdmbgwd, me, master, rdxzb,
     &  zmtb, zogw, tau_mtb, tau_ogw, tau_tofd,
     &  du3dt_mtb, du3dt_ogw, du3dt_tms)

    ! 1) non-stationary GW-scheme with GMAO/MERRA GW-forcing
    call slat_geos5_tamp(im, tamp_mpa, xlat_d, tau_ngw)

    ! 2) non-stationary GW-scheme with GEOS-5/MERRA GW-forcing
    call fv3_ugwp_solv2_v0(im, levs, dtp,
     &   Statein%tgrs, Statein%ugrs, Statein%vgrs,
     &   Statein%qgrs(:,:,1),
     &   Statein%prsl, Statein%prsi, Statein%phil, xlat_d,
     &   sinlat, coslat, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,
     &   tau_ngw, me, master, kdt)

    do k=1,levs
        do i=1,im
            gw_dtdt(i,k) = pngw*gw_dtdt(i,k)+ pogw*Pdtdt(i,k)
            gw_dudt(i,k) = pngw*gw_dudt(i,k)+ pogw*Pdudt(i,k)
            gw_dvdt(i,k) = pngw*gw_dvdt(i,k)+ pogw*Pdvdt(i,k)
            gw_kdis(i,k) = pngw*gw_kdis(i,k)+ pogw*Pkdis(i,k)
        enddo
    enddo  

    return



    end subroutine cires_ugwp_run




end module cires_ugwp

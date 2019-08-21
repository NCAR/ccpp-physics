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

    use cires_ugwp_module, only: knob_ugwp_version, cires_ugwp_mod_init, cires_ugwp_mod_finalize!, GWDPS_V0, slat_geos5_tamp, fv3_ugwp_solv2_v0, edmix_ugwp_v0, diff_1d_wtend, diff_1d_ptend
 

    implicit none

    logical :: is_initialized = .False.

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
!! | latr             | number_of_latitude_points                                                     | number of global points in y-dir (j) along the meridian | count  | 0    | integer   |           | none   | F        |
!! | levs             | vertical_dimension                                                            | number of vertical levels                               | count  | 0    | integer   |           | in     | F        |
!! | ak               | a_parameter_of_the_hybrid_coordinate                                          | a parameter for sigma pressure level calculations       | Pa     | 0    | real      | kind_phys | none   | F        |
!! | bk               | b_parameter_of_the_hybrid_coordinate                                          | b parameter for sigma pressure level calculations       | none   | 0    | real      | kind_phys | none   | F        |
!! | dtp              | time_step_for_physics                                                         | physics timestep                                        | s      | 0    | real      | kind_phys | in     | F        |
!! | cdmvgwd          | multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag | multiplication factors for cdmb and gwd                 | none   | 1    | real      | kind_phys | in     | F        |
!! | cgwf             | multiplication_factors_for_convective_gravity_wave_drag                       | multiplication factor for convective GWD                | none   | 1    | real      | kind_phys | in     | F        | 
!! | pa_rf_in         | pressure_cutoff_for_rayleigh_damping                                          | pressure level from which Rayleigh Damping is applied   | Pa     | 0    | real      | kind_phys | none   | F        |
!! | tau_rf_in        | time_scale_for_rayleigh_damping                                               | time scale for Rayleigh damping in days                 | d      | 0    | real      | kind_phys | none   | F        |
!! | errmsg           | ccpp_error_message                                                            | error message for error handling in CCPP                | none   | 0    | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                               | error flag for error handling in CCPP                   | flag   | 0    | integer   |           | out    | F        |
!!
#endif
! -----------------------------------------------------------------------
!
    subroutine cires_ugwp_init (me, master, nlunit, logunit, fn_nml2, &
                lonr, latr, levs, ak, bk, dtp, cdmvgwd, cgwf,  &
                pa_rf_in, tau_rf_in, errmsg, errflg)

    use physcons,      only: con_p0 ! pref in VAY's code

!----  initialization of cires_ugwp 
    implicit none
    
    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit
    integer, intent (in) :: lonr
    integer, intent (in) :: levs
    integer, intent (in) :: latr   
    real,    intent(in), dimension(levs+1):: ak, bk
    !real(kind=kind_phys),    intent (in) :: ak(levs+1), bk(levs+1)
    real,    intent (in) :: dtp   
    real,    intent (in) :: cdmvgwd(2), cgwf(2)             ! "scaling" controls for "old" GFS-GW schemes     
    real,    intent (in) :: pa_rf_in, tau_rf_in
      
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
                              lonr, latr, levs, ak, bk, con_p0, dtp, cdmvgwd, cgwf, &
                              pa_rf_in, tau_rf_in)

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
!! | local_name       | standard_name                                                                 | long_name                                               | units  | rank | type      | kind      | intent | optional |
!! |------------------|-------------------------------------------------------------------------------|---------------------------------------------------------|--------|------|-----------|-----------|--------|----------|
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
!! | do_ugwp          | do_ugwp                                                                        | flag to activate CIRES UGWP                                  | flag      | 0    | logical   |           | none   | F        |
!! | me               | mpi_rank                                                                       | MPI rank of current process                                  | index     | 0    | integer   |           | in     | F        |
!! | master           | mpi_root                                                                       | MPI rank of master process                                   | index     | 0    | integer   |           | in     | F        |
!! | im               | horizontal_loop_extent                                                         | horizontal                                                   | count     | 0    | integer   |           | in     | F        |
!! | levs             | vertical_dimension                                                             | number of vertical levels                                    | count     | 0    | integer   |           | in     | F        |
!! | ntrac            | number_of_tracers                                                              | number of tracers                                            | count     | 0    | integer   |           | in     | F        |
!! | dtp              | time_step_for_physics                                                          | physics timestep                                             | s         | 0    | real      | kind_phys | in     | F        |
!! | kdt              | index_of_time_step                                                             | current forecast iteration                                   | index     | 0    | integer   |           | in     | F        |
!! | lonr             | number_of_equatorial_longitude_points                                          | number of global points in x-dir (i) along the equator       | count     | 0    | integer   |           | in     | F        |
!! | oro              | orography                                                                      | orography                                                    | m         | 1    | real      | kind_phys | in     | F        |
!! | oro_uf           | orography_unfiltered                                                           | unfiltered orography                                         | m         | 1    | real      | kind_phys | in     | F        |
!! | hprime           | standard_deviation_of_subgrid_orography                                        | standard deviation of subgrid orography                      | m         | 1    | real      | kind_phys | out    | F        |
!! | nmtvr            | number_of_statistical_measures_of_subgrid_orography                            | number of topographic variables in GWD                       | count     | 0    | integer   |           | in     | F        |
!! | oc               | convexity_of_subgrid_orography                                                 | convexity of subgrid orography                               | none      | 1    | real      | kind_phys | out    | F        |
!! | oa4              | asymmetry_of_subgrid_orography                                                 | asymmetry of subgrid orography                               | none      | 2    | real      | kind_phys | out    | F        |
!! | clx              | fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height        | horizontal fraction of grid box covered by subgrid orography higher than critical height | frac    | 2 | real | kind_phys | out | F |
!! | theta            | angle_from_east_of_maximum_subgrid_orographic_variations                       | angle with_respect to east of maximum subgrid orographic variations                      | degrees | 1 | real | kind_phys | out | F |
!! | sigma            | slope_of_subgrid_orography                                                     | slope of subgrid orography                                   | none      | 1    | real      | kind_phys | out    | F        |
!! | gamma            | anisotropy_of_subgrid_orography                                                | anisotropy of subgrid orography                              | none      | 1    | real      | kind_phys | out    | F        |
!! | elvmax           | maximum_subgrid_orography                                                      | maximum of subgrid orography                                 | m         | 1    | real      | kind_phys | out    | F        |
!! | do_tofd          | turb_oro_form_drag_flag                                                        | flag for turbulent orographic form drag                      | flag      | 0    | logical   |           | none   | F        | 
!! | cdmbgwd          | multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag  | multiplication factors for cdmb and gwd                      | none      | 1    | real      | kind_phys | in     | F        |
!! | xlat             | latitude                                                                       | grid latitude in radians                                     | radians   | 1    | real      | kind_phys | in     | F        | 
!! | xlat_d           | latitude_degree                                                                | latitude in degrees                                          | degree    | 1    | real      | kind_phys | in     | F        |
!! | sinlat           | sine_of_latitude                                                               | sine of the grid latitude                                    | none      | 1    | real      | kind_phys | in     | F        |
!! | coslat           | cosine_of_latitude                                                             | cosine of the grid latitude                                  | none      | 1    | real      | kind_phys | in     | F        |  
!! | area             | cell_area                                                                      | area of the grid cell                                        | m2        | 1    | real      | kind_phys | none   | F        |
!! | ugrs             | x_wind                                                                         | zonal wind                                                   | m s-1     | 2    | real      | kind_phys | in     | F        |
!! | vgrs             | y_wind                                                                         | meridional wind                                              | m s-1     | 2    | real      | kind_phys | in     | F        |
!! | tgrs             | air_temperature                                                                | model layer mean temperature                                 | K         | 2    | real      | kind_phys | in     | F        |
!! | qgrs             | tracer_concentration                                                           | model layer mean tracer concentration                        | kg kg-1   | 3    | real      | kind_phys | in     | F        |
!! | prsi             | air_pressure_at_interface                                                      | air pressure at model layer interfaces                       | Pa        | 2    | real      | kind_phys | in     | F        |
!! | prsl             | air_pressure                                                                   | mean layer pressure                                          | Pa        | 2    | real      | kind_phys | in     | F        |
!! | prslk            | dimensionless_exner_function_at_model_layers                                   | dimensionless Exner function at model layer centers          | none      | 2    | real      | kind_phys | in     | F        |
!! | phii             | geopotential_at_interface                                                      | geopotential at model layer interfaces                       | m2 s-2    | 2    | real      | kind_phys | in     | F        |
!! | phil             | geopotential                                                                   | geopotential at model layer centers                          | m2 s-2    | 2    | real      | kind_phys | in     | F        |
!! | del              | air_pressure_difference_between_midlayers                                      | air pressure difference between midlayers                    | Pa        | 2    | real      | kind_phys | none   | F        |
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
!! | dudt_mtb         | instantaneous_change_in_x_wind_due_to_mountain_blocking_drag                   | instantaneous change in x wind due to mountain blocking drag | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | dudt_ogw         | instantaneous_change_in_x_wind_due_to_orographic_gravity_wave_drag             | instantaneous change in x wind due to orographic gw drag     | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | dudt_tms         | instantaneous_change_in_x_wind_due_to_turbulent_orographic_form_drag           | instantaneous change in x wind due to TOFD                   | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | dudt             | tendency_of_x_wind_due_to_model_physics                                        | zonal wind tendency due to model physics                     | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | dvdt             | tendency_of_y_wind_due_to_model_physics                                        | meridional wind tendency due to model physics                | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | dtdt             | tendency_of_air_temperature_due_to_model_physics                               | air temperature tendency due to model physics                | K s-1     | 2    | real      | kind_phys | none   | F        |
!! | rdxzb            | level_of_dividing_streamline                                                   | level of the dividing streamline                             | none      | 1    | real      | kind_phys | none   | F        |
!! | errmsg           | ccpp_error_message                                                             | error message for error handling in CCPP                     | none      | 0    | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                                | error flag for error handling in CCPP                        | flag      | 0    | integer   |           | out    | F        |
!!
#endif

! subroutines original
     subroutine cires_ugwp_run(do_ugwp, me,  master, im,  levs, ntrac, dtp, kdt, lonr,  &
         oro, oro_uf, hprime, nmtvr, oc, theta, sigma, gamma, elvmax, clx, oa4,& 
         do_tofd,  cdmbgwd,  xlat, xlat_d, sinlat, coslat, area,               &
         ugrs, vgrs, tgrs, qgrs, prsi, prsl, prslk, phii, phil,                &
         del, kpbl, dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,        &
         tau_tofd, tau_mtb, tau_ogw, tau_ngw, zmtb, zlwb, zogw,                & 
         dudt_mtb,dudt_ogw, dudt_tms, dudt, dvdt, dtdt, rdxzb,                 &
         errmsg, errflg)

               
    use machine,       only: kind_phys
    use physcons,      only: con_g, con_pi, con_cp, con_rd, con_rv, con_fvirt

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
    IF (do_ugwp) THEN

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



    ENDIF ! do_ugwp


    end subroutine cires_ugwp_run




end module cires_ugwp

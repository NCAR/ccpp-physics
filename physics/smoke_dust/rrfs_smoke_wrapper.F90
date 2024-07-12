!>\file rrfs_smoke_wrapper.F90
!! This file is CCPP driver of RRFS Smoke and Dust
!! Haiqin.Li@noaa.gov 02/2021

 module rrfs_smoke_wrapper

   use mpi_f08
   use machine ,              only : kind_phys
   use rrfs_smoke_config,     only : kemit, dust_opt, seas_opt, do_plumerise,           &
                                     addsmoke_flag, plumerisefire_frq, wetdep_ls_opt,   &
                                     drydep_opt, pm_settling, aero_ind_fdb, ebb_dcycle, &
                                     dbg_opt,smoke_forecast,wetdep_ls_alpha,do_rrfs_sd, &
                                     ebb_dcycle, extended_sd_diags,add_fire_heat_flux,  &
                                     num_moist, num_chem, num_emis_seas, num_emis_dust, &
                                     p_qv, p_atm_shum, p_atm_cldq,                      &
                                     p_smoke, p_dust_1, p_coarse_pm, epsilc, n_dbg_lines
   use dust_data_mod,         only : dust_alpha, dust_gamma, dust_moist_opt,            &
                                     dust_moist_correction, dust_drylimit_factor
   use seas_mod,              only : gocart_seasalt_driver
   use dust_fengsha_mod,      only : gocart_dust_fengsha_driver
   use dep_dry_simple_mod,    only : dry_dep_driver_simple
   use dep_dry_emerson_mod,   only : dry_dep_driver_emerson
   use module_wetdep_ls,      only : wetdep_ls
   use module_plumerise,      only : ebu_driver
   use module_add_emiss_burn, only : add_emis_burn
   use coarsepm_settling_mod, only : coarsepm_settling_driver

   implicit none

   private

   public :: rrfs_smoke_wrapper_run, rrfs_smoke_wrapper_init

   integer :: plume_wind_eff

contains

!>\defgroup rrfs_smoke_wrapper rrfs-sd emission driver Module
!> \ingroup gsd_chem_group
!! This is the rrfs-sd emission driver Module

!> \section arg_table_rrfs_smoke_wrapper_init Argument Table
!! \htmlinclude rrfs_smoke_wrapper_init.html
!!
  subroutine rrfs_smoke_wrapper_init( seas_opt_in,                                & ! sea salt namelist 
                              drydep_opt_in, pm_settling_in,                      & ! Dry Dep namelist
                              wetdep_ls_opt_in,wetdep_ls_alpha_in,                & ! Wet dep namelist
                              rrfs_sd, do_plumerise_in, plumerisefire_frq_in,     & ! smoke namelist 
                              plume_wind_eff_in,add_fire_heat_flux_in,            & ! smoke namelist
                              addsmoke_flag_in, ebb_dcycle_in, smoke_forecast_in, & ! Smoke namelist
                              dust_opt_in, dust_alpha_in, dust_gamma_in,          & ! Dust namelist
                              dust_moist_opt_in,                                  & ! Dust namelist
                              dust_moist_correction_in, dust_drylimit_factor_in,  & ! Dust namelist                        
                              aero_ind_fdb_in,                                    & ! Feedback namelist
                              extended_sd_diags_in,dbg_opt_in,                    & ! Other namelist
                              errmsg, errflg, n_dbg_lines_in                      )
                             

!>-- Namelist
  real(kind_phys), intent(in) :: dust_alpha_in, dust_gamma_in, wetdep_ls_alpha_in
  real(kind_phys), intent(in) :: dust_moist_correction_in
  real(kind_phys), intent(in) :: dust_drylimit_factor_in
  integer,         intent(in) :: dust_opt_in,dust_moist_opt_in, wetdep_ls_opt_in, pm_settling_in, seas_opt_in
  integer,         intent(in) :: drydep_opt_in
  logical,         intent(in) :: aero_ind_fdb_in,dbg_opt_in, extended_sd_diags_in, add_fire_heat_flux_in
  integer,         intent(in) :: smoke_forecast_in, plume_wind_eff_in, plumerisefire_frq_in, n_dbg_lines_in
  integer,         intent(in) :: addsmoke_flag_in, ebb_dcycle_in
  logical,         intent(in) :: do_plumerise_in, rrfs_sd
  character(len=*),intent(out):: errmsg
  integer,        intent(out) :: errflg

     errmsg = ''
     errflg = 0

!>-- Assign namelist values 
  !>-Dust
     dust_alpha            = dust_alpha_in
     dust_gamma            = dust_gamma_in
     dust_moist_opt        = dust_moist_opt_in 
     dust_moist_correction = dust_moist_correction_in
     dust_drylimit_factor  = dust_drylimit_factor_in
     dust_opt              = dust_opt_in
  !>-Sea Salt
     seas_opt              = seas_opt_in
  !>-Dry and wet deposition
     drydep_opt            = drydep_opt_in
     pm_settling           = pm_settling_in
     wetdep_ls_opt         = wetdep_ls_opt_in
     wetdep_ls_alpha       = wetdep_ls_alpha_in
  !>-Smoke
     do_rrfs_sd            = rrfs_sd
     ebb_dcycle            = ebb_dcycle_in
     do_plumerise          = do_plumerise_in
     plumerisefire_frq     = plumerisefire_frq_in
     addsmoke_flag         = addsmoke_flag_in
     smoke_forecast        = smoke_forecast_in
     plume_wind_eff        = plume_wind_eff_in
     add_fire_heat_flux    = add_fire_heat_flux_in
  !>-Feedback
     aero_ind_fdb          = aero_ind_fdb_in
  !>-Other
     extended_sd_diags     = extended_sd_diags_in
     dbg_opt               = dbg_opt_in
     n_dbg_lines           = n_dbg_lines_in

  end subroutine rrfs_smoke_wrapper_init

!! \section arg_table_rrfs_smoke_wrapper_run Argument Table
!! \htmlinclude rrfs_smoke_wrapper_run.html
!!
!>\section rrfs_smoke_wrapper rrfs-sd Scheme General Algorithm
!> @{
    subroutine rrfs_smoke_wrapper_run(im, kte, kme, ktau, dt, garea, land, jdate,          &
                   u10m, v10m, ustar, rlat, rlon, tskin, pb2d, t2m, dpt2m,                 &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum, w,                  &
                   nsoil, smc, tslb, vegtype_dom, vegtype_frac, soiltyp, nlcat,            &
                   dswsfc, zorl, snow, julian,recmol,                                      &
                   idat, rain_cpl, rainc_cpl, hf2d, g, pi, con_cp, con_rd, con_fv,         &
                   dust12m_in, emi_ant_in, smoke_RRFS, smoke2d_RRFS,                       &
                   ntrac, qgrs, gq0, chem3d, tile_num,                                     &
                   ntsmoke, ntdust, ntcoarsepm, imp_physics, imp_physics_thompson,         &
                   nwfa, nifa, emanoc, emdust, emseas, drydep_flux_out, wetdpr,            &
                   ebb_smoke_in, frp_output, coef_bb, fire_type_out,                       &
                   ebu_smoke,fhist,min_fplume,                                             &
                   max_fplume, hwp, hwp_ave, wetness, ndvel, ddvel_inout,                  &
                   peak_hr_out,lu_nofire_out,lu_qfire_out,                                 &
                   fire_heat_flux_out, frac_grid_burned_out, kpbl,oro,                     &
                   uspdavg, hpbl_thetav, mpicomm, mpirank, mpiroot, errmsg,errflg          )
        
    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,nsoil,tile_num,jdate(8),idat(8)
    integer,        intent(in) :: ntrac, ntsmoke, ntdust, ntcoarsepm, ndvel, nlcat
    real(kind_phys),intent(in) :: dt, julian, g, pi, con_cp, con_rd, con_fv

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer,         dimension(:),     intent(in)    :: land, vegtype_dom, soiltyp
    real(kind_phys), dimension(:,:),   intent(in), optional :: smc, tslb
    real(kind_phys), dimension(:,:,:), intent(in), optional :: dust12m_in
    real(kind_phys), dimension(:,:,:), intent(in), optional :: smoke_RRFS
    real(kind_phys), dimension(:,:),   intent(in), optional :: smoke2d_RRFS
    real(kind_phys), dimension(:,:),   intent(in), optional :: emi_ant_in
    real(kind_phys), dimension(:),     intent(in)    :: u10m, v10m, ustar, dswsfc,         &
                           recmol, garea, rlat,rlon, tskin, pb2d, zorl, snow,              &
                           rain_cpl, rainc_cpl, hf2d, t2m, dpt2m 
    real(kind_phys), dimension(:,:),   intent(in)    :: vegtype_frac
    real(kind_phys), dimension(:,:),   intent(in)    :: ph3d, pr3d
    real(kind_phys), dimension(:,:),   intent(in)    :: phl3d, prl3d, tk3d, us3d, vs3d, spechum, w
    real(kind_phys), dimension(:,:,:), intent(inout) :: qgrs, gq0
    real(kind_phys), dimension(:,:,:), intent(inout), optional :: chem3d
    real(kind_phys), dimension(:),     intent(inout), optional :: emdust, emseas, emanoc
    real(kind_phys), dimension(:),     intent(inout), optional :: ebb_smoke_in,coef_bb, frp_output, fhist
    real(kind_phys), dimension(:,:),   intent(inout), optional :: ebu_smoke
    real(kind_phys), dimension(:),     intent(out  ), optional :: fire_heat_flux_out, frac_grid_burned_out
    real(kind_phys), dimension(:),     intent(inout), optional :: max_fplume, min_fplume, uspdavg, hpbl_thetav
    real(kind_phys), dimension(:),     intent(inout), optional :: hwp, peak_hr_out
    real(kind_phys), dimension(:),     intent(inout), optional :: hwp_ave
    real(kind_phys), dimension(:,:),   intent(inout), optional :: nwfa, nifa
    real(kind_phys), dimension(:,:),   intent(inout), optional :: ddvel_inout
    real(kind_phys), dimension(:,:),   intent(inout), optional :: drydep_flux_out
    real(kind_phys), dimension(:,:),   intent(inout), optional :: wetdpr
    real(kind_phys), dimension(:),     intent(in),    optional :: wetness
    real(kind_phys), dimension(:),     intent(out),   optional :: lu_nofire_out,lu_qfire_out
    integer,         dimension(:),     intent(out),   optional :: fire_type_out
    integer,                           intent(in)    :: imp_physics, imp_physics_thompson
    integer,         dimension(:),     intent(in)    :: kpbl
    real(kind_phys), dimension(:),     intent(in)    :: oro
    character(len=*),                  intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

!>-- Local Variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: ebu
    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,  &
                     p_phy,pi_phy,wind_phy,theta_phy,z_at_w, dz8w, p8w, t8w,      &
                     rho_phy, vvel, zmid
    real(kind_phys), dimension(ims:im, jms:jme) :: frp_inst, u10, v10, ust, tsk,  &
                     xland, xlat, xlong, dxy, pbl, hfx, rnav, hwp_local,          &
                     wetdpr_smoke_local, wetdpr_dust_local, wetdpr_coarsepm_local
!>- sea salt & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, 1, jms:jme, 1:num_emis_seas  )  :: emis_seas
    real(kind_phys), dimension(ims:im, jms:jme) :: seashelp
!>-- indexes, time
    integer :: ide, ime, ite, kde, julday
    real(kind_phys) :: gmt
!>- dust & chemistry variables
    real(kind_phys), dimension(ims:im, jms:jme) :: ssm, rdrag, uthr, snowh  ! fengsha dust
    real(kind_phys), dimension(ims:im, jms:jme) :: rmol, swdown, znt, clayf, sandf
    real(kind_phys), dimension(ims:im, nlcat, jms:jme) :: vegfrac
    real(kind_phys), dimension(ims:im, nsoil, jms:jme) :: smois, stemp
    real(kind_phys), dimension(ims:im, 1:1, jms:jme, 1:num_emis_dust) :: emis_dust
    integer,         dimension(ims:im, jms:jme) :: isltyp, ivgtyp
!>- plume variables
    ! -- buffers
    real(kind_phys), dimension(ims:im, jms:jme )  :: coef_bb_dc, flam_frac, frp_in,        &
                                          fire_hist, peak_hr, lu_nofire, lu_qfire, ebu_in, &
                                                     fire_end_hr, hwp_day_avg, kpbl_thetav,&
                                                     uspdavg2, hpbl_thetav2
    integer,         dimension(ims:im, jms:jme )  :: min_fplume2, max_fplume2, fire_type
    logical :: call_plume, reset_hwp_ave, avg_hwp_ave
!>- optical variables
    real(kind_phys), dimension(ims:im, jms:jme, ndvel) :: ddvel, settling_flux, drydep_flux_local
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, ndvel) :: vgrav
!>-- anthropogentic variables
    real(kind_phys), dimension(ims:im) :: emis_anoc
    real(kind_phys), dimension(ims:im, jms:jme, 1) :: sedim
!> -- parameter to caluclate wfa&ifa (m)
    real(kind_phys), parameter :: mean_diameter1= 4.E-8, sigma1=1.8
    real(kind_phys), parameter :: mean_diameter2= 1.E-6, sigma2=1.8
    real(kind_phys) :: fact_wfa, fact_ifa
!> -- aerosol density (kg/m3)
    real(kind_phys), parameter :: density_dust= 2.6e+3, density_sulfate=1.8e+3
    real(kind_phys), parameter :: density_oc  = 1.4e+3, density_seasalt=2.2e+3
    real(kind_phys), dimension(im) :: daero_emis_wfa, daero_emis_ifa
!> -- other
    real(kind_phys), dimension(im) :: wdgust, snoweq
    integer :: current_month, current_hour, hour_int
    real(kind_phys) :: curr_secs
    real(kind_phys) :: factor, factor2, factor3
    integer :: nbegin, nv
    integer :: i, j, k, kp, n
! MPI variables
    integer :: mpiid
    type(MPI_comm), intent(in) :: mpicomm
    integer, intent(in) :: mpirank
    integer, intent(in) :: mpiroot

    mpiid = mpirank

    errmsg = ''
    errflg = 0

    if (.not. do_rrfs_sd) return

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    min_fplume2 = 0
    max_fplume2 = 0
    uspdavg2 = 0.
    hpbl_thetav2 = 0.
    emis_seas   = 0.
    emis_dust   = 0.
    peak_hr     = 0.
    fire_type   = 0
    lu_qfire    = 0.
    lu_nofire   = 0.
    flam_frac   = 0.
    daero_emis_wfa = 0.
    daero_emis_ifa = 0.

    rnav = 0.

    curr_secs = ktau * dt
    current_month=jdate(2)      ! needed for the dust input data
    current_hour =jdate(5)+1    ! =1 at 00Z
    hour_int=floor(ktau*dt/3600.)      ! hours since the simulation start
    gmt    = real(mod(idat(5)+hour_int,24))
    julday = int(julian)

    do nv=1,ndvel
    do i=its,ite
      ddvel(i,1,nv)=ddvel_inout(i,nv)
    enddo
    enddo

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rnav(i,1)=max((rain_cpl(i)-rainc_cpl(i))*1000., 0.) ! meter to mm
! coef_bb initializes as clear_val (from GFS_typedefs.F90)
! at ktau = 1, coef_bb_dc is set = 1.0
     coef_bb_dc(i,1) = coef_bb(i)
! fhist   initializes as 1.        (from GFS_typedefs.F90)
     fire_hist (i,1) = fhist  (i)
     peak_hr   (i,1) = peak_hr_out(i)
    enddo

    ! Is this a reset timestep (00:00 + dt)?
    reset_hwp_ave = mod(int(curr_secs-dt),3600) == 0
    avg_hwp_ave   = mod(int(curr_secs),3600) == 0

    ! plumerise frequency in minutes set up by the namelist input
    call_plume       = (do_plumerise .and. (plumerisefire_frq > 0))
    if (call_plume) call_plume = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0) .or. (ktau == 2)
    
!>- get ready for chemistry run
    call rrfs_smoke_prep(                                               &
        ktau,current_month, current_hour, gmt, con_rd, con_fv, con_cp,  &
        u10m,v10m,ustar,land,garea,rlat,rlon,tskin,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                 &
        nsoil,smc,tslb,vegtype_dom,soiltyp,                             &
        nlcat,vegtype_frac,dswsfc,zorl,                                 &
        snow,dust12m_in,emi_ant_in,smoke_RRFS,smoke2d_RRFS,coef_bb_dc,  &
        hf2d, pb2d, g, pi, hour_int, peak_hr,                           &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,pi_phy,wind_phy,theta_phy,          &
        rho_phy,dz8w,p8w,t8w,recmol,                                    &
        z_at_w,vvel,zmid,                                               &
        ntrac,gq0,                                                      &
        num_chem,num_moist,                                             &
        ntsmoke, ntdust,ntcoarsepm,                                     &
        moist,chem,ebu_in,kpbl_thetav,ebb_smoke_in,                     &
        fire_hist,frp_in, hwp_day_avg, fire_end_hr,                     &
        emis_anoc,smois,stemp,ivgtyp,isltyp,vegfrac,rmol,swdown,znt,    &
        hfx,pbl,snowh,clayf,rdrag,sandf,ssm,uthr,oro, hwp_local,        &
        t2m,dpt2m,wetness,kpbl,                                         &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte                                       )

    IF (ktau==1) THEN
      ebu = 0.
      do j=jts,jte
      do i=its,ite
        ebu(i,kts,j)= ebu_in(i,j)
        do k=kts+1,kte
         ebu(i,k,j)= 0.
        enddo
      enddo
      enddo
    ELSE
      do k=kts,kte
      do i=its,ite
      ! ebu is divided by coef_bb_dc since it is applied in the output
         ebu(i,k,1)=ebu_smoke(i,k) / coef_bb_dc(i,1)
      enddo
      enddo
    ENDIF

!RAR: change this to the fractional LU type; fire_type: 0- no fires, 1- Ag
! or urban fires, 2- prescribed fires in wooded area, 3- wildfires 
    if (ebb_dcycle==2) then
    do j=jts,jte
      do i=its,ite
        if (ebu_in(i,j)<0.01) then
           fire_type(i,j) = 0
           lu_nofire(i,j) = 1.0
        else
          ! Permanent wetlands, snow/ice, water, barren tundra 
          lu_nofire(i,j) = vegfrac(i,11,j) + vegfrac(i,15,j) + vegfrac(i,17,j) + vegfrac(i,20,j)
          ! cropland, urban, cropland/natural mosaic, barren and sparsely vegetated
          lu_qfire(i,j)  = vegfrac(i,12,j) + vegfrac(i,13,j) + vegfrac(i,14,j) + vegfrac(i,16,j)
          if (lu_nofire(i,j)>0.95) then
             fire_type(i,j) = 0
          else if (lu_qfire(i,j)>0.95) then
             fire_type(i,j) = 1
          else
             fire_type(i,j) = 3  ! RAR: need to add another criteria for fire_type=2, i.e. prescribed fires 
          end if
        end if
      end do
    end do
    endif

!>- compute sea-salt
    ! -- compute sea salt (opt=1)
    if (seas_opt == 1) then
    call gocart_seasalt_driver(dt,rri,t_phy,                            &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,ust,p8w,tsk,              &
        xland,xlat,xlong,dxy,g,emis_seas,pi,                            &
        seashelp,num_emis_seas,num_chem,seas_opt,                       &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)
    endif

    !-- compute dust (opt=5)
    if (dust_opt==1) then
       call gocart_dust_fengsha_driver(dt,chem,rho_phy,                 &
            smois,stemp,p8w,ssm,                                        &
            isltyp,snowh,xland,dxy,g,emis_dust,ust,znt,                 &
            clayf,sandf,rdrag,uthr,                                     &
            num_emis_dust,num_chem,nsoil,                               &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte)
    end if

    ! compute wild-fire plumes
    !-- to add a namelist option to turn on/off plume raising
    !--- replace plumerise_driver with HRRR-smoke 05/10/2021
    !-- /scratch2/BMC/ap-fc/Ravan/rapid-refresh/WRFV3.9/smoke
    ! Every hour (per namelist) the ebu_driver is called to calculate ebu, but
    ! the plumerise is controlled by the namelist option of plumerise_flag
    if (add_fire_heat_flux) then
     WRITE(1000+mpiid,*) 'Entered add_fire_heat_flux at timestep:',ktau
     do i = its,ite
       if ( coef_bb_dc(i,1)*frp_in(i,1) .ge. 1.E7 ) then
          fire_heat_flux_out(i) = min(max(0.,0.88*coef_bb_dc(i,1)*frp_in(i,1) / &
                                  0.55/dxy(i,1)) ,5000.) ! JLS - W m-2 [0 - 10,000]
          frac_grid_burned_out(i) = min(max(0., 1.3*0.0006*coef_bb_dc(i,1)*frp_in(i,1)/dxy(i,1) ),1.)
       else
          fire_heat_flux_out(i)   = 0.0
          frac_grid_burned_out(i) = 0.0
       endif
     enddo
    endif
    if (call_plume) then
        ! Apply the diurnal cycle coefficient to frp_inst ()
        do j=jts,jte
        do i=its,ite
         frp_inst(i,j) = frp_in(i,j)*coef_bb_dc(i,j)
        enddo
        enddo

        call ebu_driver (                                              &
                   flam_frac,ebu_in,ebu,                               &
                   theta_phy,moist(:,:,:,p_qv),                        &
                   rho_phy,vvel,u_phy,v_phy,pi_phy,wind_phy,           &
                   z_at_w,zmid,g,con_cp,con_rd,                        &
                   frp_inst, min_fplume2, max_fplume2,                 &
                   plume_wind_eff,                                     &
                   kpbl_thetav,                                        &
                   ids,ide, jds,jde, kds,kde,                          &
                   ims,ime, jms,jme, kms,kme,                          &
                   its,ite, jts,jte, kts,kte, errmsg, errflg, curr_secs, &
                   xlat, xlong, uspdavg2, hpbl_thetav2, mpiid  )
        if(errflg/=0) return
    end if

    ! -- add biomass burning emissions at every timestep
    if (addsmoke_flag == 1) then
    call add_emis_burn(dt,dz8w,rho_phy,pi,                           &
                       chem,julday,gmt,xlat,xlong,                   &
                       fire_end_hr, peak_hr,curr_secs,               &
                       coef_bb_dc,fire_hist,hwp_local,hwp_day_avg,   &
                       swdown,ebb_dcycle,ebu_in,ebu,fire_type,       &
                       ids,ide, jds,jde, kds,kde,                    &
                       ims,ime, jms,jme, kms,kme,                    &
                       its,ite, jts,jte, kts,kte , mpiid             )
    endif

    !>-- compute coarsepm setting if using simple dry dep option and
    !    pm_settling is on. This is necessary becasue the simple scheme
    !    does not have an explicty settling routine, Emersion (opt=1) does.
    if (drydep_opt == 2 .and. pm_settling == 1) then
    call coarsepm_settling_driver(dt,t_phy,                          &
                                  chem(:,:,:,p_coarse_pm),           &
                                  rho_phy,dz8w,p8w,p_phy,sedim,      &
                                  dxy,g,1,                           &
                                  ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte          )
    endif
    !>-- compute dry deposition, based on Emerson et al., (2020)
    if (drydep_opt == 1) then
    call dry_dep_driver_emerson(rmol,ust,znt,ndvel,ddvel,            &
       vgrav,chem,dz8w,snowh,t_phy,p_phy,rho_phy,ivgtyp,g,dt,        & 
       pm_settling,drydep_flux_local,settling_flux,dbg_opt,          &
       ids,ide, jds,jde, kds,kde,                                    &
       ims,ime, jms,jme, kms,kme,                                    &
       its,ite, jts,jte, kts,kte, curr_secs, mpiid, xlat, xlong      )
       do nv=1,ndvel
          do i=its,ite
             ddvel_inout(i,nv)=ddvel(i,1,nv)
          enddo
       enddo
    !>-- compute dry deposition based on simple parameterization (HRRR-Smoke)
    elseif (drydep_opt == 2) then
    call dry_dep_driver_simple(rmol,ust,ndvel,ddvel,                 &
       ids,ide, jds,jde, kds,kde,                                    &
       ims,ime, jms,jme, kms,kme,                                    &
       its,ite, jts,jte, kts,kte                                     )
       do nv=1,ndvel
          do i=its,ite
             ddvel_inout(i,nv)=ddvel(i,1,nv)
          enddo
       enddo
    else
       ddvel_inout(:,:)=0.
    endif

!>- large-scale wet deposition
    if (wetdep_ls_opt == 1) then
       call  wetdep_ls(dt,chem,rnav,moist,                       &
                     rho_phy,num_chem,num_moist,ndvel, dz8w,vvel,&
                     wetdpr_smoke_local, wetdpr_dust_local,      &
                     wetdpr_coarsepm_local,                      &
                     ids,ide, jds,jde, kds,kde,                  &
                     ims,ime, jms,jme, kms,kme,                  &
                     its,ite, jts,jte, kts,kte                   )
       if ( extended_sd_diags .or. dbg_opt) then
         do i = its, ite
            wetdpr(i,1) = wetdpr(i,1) + wetdpr_smoke_local   (i,1)
            wetdpr(i,2) = wetdpr(i,2) + wetdpr_dust_local    (i,1)
            wetdpr(i,3) = wetdpr(i,3) + wetdpr_coarsepm_local(i,1)
         enddo
       endif
    endif

! Smoke emisisons diagnostic, RAR: let's multiply by coef_bb_dc before output
! Since ebu_smoke includes coef_bb_dc, we need to divide by coef_bb_dc when it
! comes back into the wrapper.
    do k=kts,kte
     do i=its,ite
       ebu_smoke(i,k)=ebu(i,k,1) * coef_bb_dc(i,1)
     enddo
    enddo

!---- diagnostic output of hourly wildfire potential (07/2021)
    if (ktau == 1 .or. reset_hwp_ave) then
       hwp_ave = 0.
    endif
    hwp = 0.
    do i=its,ite
      hwp(i)=hwp_local(i,1)
      hwp_ave(i) = hwp_ave(i) + hwp(i)*dt
      if ( ktau == 1) then
         hwp_ave(i) = hwp_ave(i) / dt
      elseif ( avg_hwp_ave ) then
         hwp_ave(i) = hwp_ave(i) / 3600._kind_phys
      endif
    enddo


!---- diagnostic output of dry deposition & gravitational settling fluxes
    if ( drydep_opt == 1 .and. (extended_sd_diags .or. dbg_opt) ) then
      do nv = 1, ndvel
       do i=its,ite
          drydep_flux_out(i,nv) = drydep_flux_out(i,nv)     +  &
                                  drydep_flux_local(i,1,nv) !+  &
                                  !settling_flux(i,1,nv)
       enddo
      enddo
    endif
!-------------------------------------
!---- put smoke stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntsmoke )  = min(5000.,max(epsilc,chem(i,k,1,p_smoke ))) 
       gq0(i,k,ntdust  )  = min(200.,max(epsilc,chem(i,k,1,p_dust_1)))
       gq0(i,k,ntcoarsepm)= min(5000.,max(epsilc,chem(i,k,1,p_coarse_pm)))
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntsmoke   )= gq0(i,k,ntsmoke )
       qgrs(i,k,ntdust    )= gq0(i,k,ntdust  )
       qgrs(i,k,ntcoarsepm)= gq0(i,k,ntcoarsepm)
       chem3d(i,k,1       )= gq0(i,k,ntsmoke )
       chem3d(i,k,2       )= gq0(i,k,ntdust  )
       chem3d(i,k,3       )= gq0(i,k,ntcoarsepm)
     enddo
    enddo
!-------------------------------------
!-- to output for diagnostics
    do i = 1, im
! RAR: let's remove the seas and ant. OC
     emseas     (i) = emis_seas(i,1,1,1)*1.e+9   ! size bin 1 sea salt emission: ug/m2/s
     emanoc     (i) = emis_anoc (i)              ! anthropogenic organic carbon: ug/m2/s
     emdust     (i) = emis_dust(i,1,1,1) + emis_dust(i,1,1,2) +   &
                      emis_dust(i,1,1,3) + emis_dust(i,1,1,4) ! dust emission: ug/m2/s
     coef_bb    (i) = coef_bb_dc(i,1)
     frp_output (i) = coef_bb_dc(i,1)*frp_in(i,1)
     fhist      (i) = fire_hist (i,1)
     min_fplume (i) = real(min_fplume2(i,1))
     max_fplume (i) = real(max_fplume2(i,1))
     fire_type_out(i)=fire_type(i,1)
     lu_nofire_out(i)=lu_nofire(i,1)
     lu_qfire_out (i)=lu_qfire(i,1)
    enddo

    do i = 1, im
      peak_hr_out(i) =  peak_hr(i,1)
    enddo

!-- to provide real aerosol emission for Thompson MP
    if (imp_physics == imp_physics_thompson .and. aero_ind_fdb) then
      fact_wfa = 1.e-9*6.0/pi*exp(4.5*log(sigma1)**2)/mean_diameter1**3
      fact_ifa = 1.e-9*6.0/pi*exp(4.5*log(sigma2)**2)/mean_diameter2**3

      do i = 1, im
         daero_emis_wfa(i) =(emanoc(i)+ebu_smoke(i,kemit))/density_oc + emseas(i)/density_seasalt
         daero_emis_ifa(i) = emdust(i)/density_dust

         daero_emis_wfa(i) = daero_emis_wfa(i)*fact_wfa*rri(i,kemit,1)/dz8w(i,kemit,1)
         daero_emis_ifa(i) = daero_emis_ifa(i)*fact_ifa*rri(i,kemit,1)/dz8w(i,kemit,1)

         nwfa(i,kemit)     = nwfa(i,kemit) + daero_emis_wfa(i)*dt
         nifa(i,kemit)     = nifa(i,kemit) + daero_emis_ifa(i)*dt

        if(land(i).eq.1)then
         nwfa(i,kemit)    = nwfa(i,kemit)*(1. - 0.10*dt/86400.)  !-- mimicking dry deposition
         nifa(i,kemit)    = nifa(i,kemit)*(1. - 0.10*dt/86400.)  !-- mimicking dry deposition
        else
         nwfa(i,kemit)    = nwfa(i,kemit)*(1. - 0.05*dt/86400.)  !-- mimicking dry deposition
         nifa(i,kemit)    = nifa(i,kemit)*(1. - 0.05*dt/86400.)  !-- mimicking dry deposition
        endif
         nwfa(i,kemit)    = MIN(2.E10,nwfa(i,kemit))
         nifa(i,kemit)    = MIN(9999.E6,nifa(i,kemit))
      enddo
    endif

 end subroutine rrfs_smoke_wrapper_run

 subroutine rrfs_smoke_prep(                                               &
        ktau,current_month,current_hour,gmt,con_rd,con_fv,con_cp,          &
        u10m,v10m,ustar,land,garea,rlat,rlon,ts2d,                         &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                    &
        nsoil,smc,tslb,vegtype_dom,soiltyp,nlcat,vegtype_frac,dswsfc,zorl, &
        snow_cpl,dust12m_in,emi_ant_in,smoke_RRFS,smoke2d_RRFS,coef_bb_dc, &
        hf2d, pb2d, g, pi, hour_int, peak_hr,                              &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                              &
        rri,t_phy,u_phy,v_phy,p_phy,pi_phy,wind_phy,theta_phy,             &
        rho_phy,dz8w,p8w,t8w,recmol,                                       &
        z_at_w,vvel,zmid,                                                  &
        ntrac,gq0,                                                         &
        num_chem, num_moist,                                               &
        ntsmoke, ntdust, ntcoarsepm,                                       &
        moist,chem,ebu_in,kpbl_thetav,ebb_smoke_in,                        &
        fire_hist,frp_in, hwp_day_avg, fire_end_hr,                        &
        emis_anoc,smois,stemp,ivgtyp,isltyp,vegfrac,rmol,swdown,           &
        znt,hfx,pbl,snowh,clayf,rdrag,sandf,ssm,uthr,oro,hwp_local,        &
        t2m,dpt2m,wetness,kpbl,                                            &
        ids,ide, jds,jde, kds,kde,                                         &
        ims,ime, jms,jme, kms,kme,                                         &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    integer, intent(in) :: current_month, current_hour, hour_int, nlcat

    !FV3 input variables
    integer, intent(in) :: nsoil, ktau
    integer, dimension(ims:ime), intent(in) :: land, vegtype_dom, soiltyp, kpbl
    integer, intent(in) :: ntrac
    real(kind=kind_phys), intent(in) :: g, pi, gmt, con_rd, con_fv, con_cp
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         u10m, v10m, ustar, garea, rlat, rlon, ts2d, dswsfc,       &
         zorl, snow_cpl, pb2d, hf2d, oro, t2m, dpt2m, wetness, recmol
    real(kind=kind_phys), dimension(ims:ime, nlcat),   intent(in) :: vegtype_frac
    real(kind=kind_phys), dimension(ims:ime, nsoil),   intent(in) :: smc,tslb
    real(kind=kind_phys), dimension(ims:ime, 12, 5),   intent(in) :: dust12m_in
    real(kind=kind_phys), dimension(ims:ime, 24, 2),   intent(in) :: smoke_RRFS
! This is a place holder for ebb_dcycle == 2, currently set to hold a single
! value, which is the previous day's average of hwp, frp, ebb, fire_end
    real(kind=kind_phys), dimension(ims:ime,     4),   intent(in) :: smoke2d_RRFS
    real(kind=kind_phys), dimension(ims:ime,     1),   intent(in) :: emi_ant_in
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !rrfs-sd variables
    integer,intent(in) ::  num_chem, num_moist, ntsmoke, ntdust, ntcoarsepm
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte


    real(kind_phys), dimension(ims:ime, jms:jme),intent(out) :: ebu_in
    
    integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel,               &
         zmid, pi_phy, theta_phy, wind_phy
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         u10, v10, ust, tsk, xland, xlat, xlong, dxy, rmol, swdown, znt,      &
         pbl, hfx, snowh, clayf, rdrag, sandf, ssm, uthr, hwp_local
    real(kind_phys), dimension(ims:ime, nlcat, jms:jme), intent(out) :: vegfrac
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, nsoil, jms:jme), intent(out) :: smois,stemp
    real(kind_phys), dimension(ims:ime,jms:jme), intent(inout) :: frp_in, fire_end_hr, fire_hist, coef_bb_dc
    real(kind_phys), dimension(ims:ime,jms:jme), intent(inout) :: hwp_day_avg, peak_hr
    real(kind_phys), dimension(ims:ime), intent(inout) :: emis_anoc,ebb_smoke_in
    real(kind_phys), parameter :: conv_frp = 1.e+06_kind_phys  ! FRP conversion factor, MW to W
    real(kind_phys), parameter :: frpc  = 1._kind_phys         ! FRP conversion factor (Regional)

    ! -- local variables
    integer i,ip,j,k,k1,kp,kk,kkp,nv,l,ll,n,nl
    real(kind_phys) :: SFCWIND,WIND,DELWIND,DZ,wdgust,snoweq,THETA
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: THETAV
    real(kind_phys), dimension(ims:ime, jms:jme) :: windgustpot
    real(kind_phys), dimension(ims:ime, jms:jme),  intent(out)   :: kpbl_thetav
    real(kind_phys), parameter :: delta_theta4gust = 0.5
    real(kind=kind_phys),parameter :: p1000mb = 100000.

    ! -- initialize fire emissions
    ebu_in         = 0._kind_phys
    ebb_smoke_in   = 0._kind_phys
    emis_anoc      = 0._kind_phys
    frp_in         = 0._kind_phys
    hwp_day_avg    = 0._kind_phys
    fire_end_hr    = 0._kind_phys

    ! -- initialize output arrays
    isltyp         = 0._kind_phys
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    theta_phy      = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    wind_phy       = 0._kind_phys
    p_phy          = 0._kind_phys
    pi_phy         = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    vvel           = 0._kind_phys
    zmid           = 0._kind_phys
    u10            = 0._kind_phys
    v10            = 0._kind_phys
    ust            = 0._kind_phys
    tsk            = 0._kind_phys
    xland          = 0._kind_phys
    xlat           = 0._kind_phys
    xlong          = 0._kind_phys
    dxy            = 0._kind_phys
    vegfrac        = 0._kind_phys
    rmol           = 0._kind_phys
    swdown         = 0._kind_phys
    znt            = 0._kind_phys
    hfx            = 0._kind_phys
    pbl            = 0._kind_phys
    snowh          = 0._kind_phys
    clayf          = 0._kind_phys
    rdrag          = 0._kind_phys
    sandf          = 0._kind_phys 
    ssm            = 0._kind_phys
    uthr           = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys

    do i=its,ite
     u10  (i,1)=u10m (i)
     v10  (i,1)=v10m (i)
     tsk  (i,1)=ts2d (i)
     ust  (i,1)=ustar(i)
     dxy  (i,1)=garea(i)
     xland(i,1)=real(land(i))
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     swdown(i,1)=dswsfc(i)
     znt  (i,1)=zorl(i)*0.01
     hfx  (i,1)=hf2d(i)
     pbl  (i,1)=pb2d(i)
     snowh(i,1)=snow_cpl(i)*0.001
     clayf(i,1)=dust12m_in(i,current_month,1)
     rdrag(i,1)=dust12m_in(i,current_month,2)
     sandf(i,1)=dust12m_in(i,current_month,3)
     ssm  (i,1)=dust12m_in(i,current_month,4)
     uthr (i,1)=dust12m_in(i,current_month,5)
     ivgtyp (i,1)=vegtype_dom (i)
     isltyp (i,1)=soiltyp(i)
     do nl = 1,nlcat
       vegfrac(i,nl,1)=vegtype_frac (i,nl)
     enddo
     rmol (i,1)=recmol (i)
    enddo
   

    do k=1,nsoil
     do j=jts,jte
      do i=its,ite
       smois(i,k,j)=smc(i,k)
       stemp(i,k,j)=tslb(i,k)
      enddo
     enddo
    enddo

    do j=jts,jte
     do i=its,ite
       z_at_w(i,kts,j)=max(0.,ph3d(i,1)/g)
     enddo
    enddo

    do j=jts,jte
     do k=kts,kte
      do i=its,ite
        dz8w(i,k,j)=abs(ph3d(i,k+1)-ph3d(i,k))/g
        z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
      enddo
     enddo
    enddo

    do j=jts,jte
     do k=kts,kte+1
      do i=its,ite
        p8w(i,k,j)=pr3d(i,k)
      enddo
     enddo
    enddo

    do j=jts,jte
      do k=kts,kte+1
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(i,kkp)
          p_phy(i,k,j)=prl3d(i,kkp)
          u_phy(i,k,j)=us3d(i,kkp)
          v_phy(i,k,j)=vs3d(i,kkp)
          pi_phy(i,k,j) = con_cp*(p_phy(i,k,j)/p1000mb)**(con_rd/con_cp)
          theta_phy(i,k,j) = t_phy(i,k,j)/pi_phy(i,k,j)*con_cp
          wind_phy(i,k,j) = sqrt(u_phy(i,k,j)**2 + v_phy(i,k,j)**2)
          ! from mp_thompson.F90 ; rho = con_eps*prsl/(con_rd*tgrs*(qv+con_eps))
          ! from mynnd
          rho_phy(i,k,j)=p_phy(i,k,j)/(con_rd*t_phy(i,k,j)) !*(1.+con_fv*spechum(i,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(i,kkp)*rri(i,k,j)/g 
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(i,kkp,1)
          moist(i,k,j,2)=gq0(i,kkp,2)
          if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          !--
          zmid(i,k,j)=phl3d(i,kkp)/g
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=2,kte
        do i=its,ite
          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
        enddo
      enddo
    enddo

    do j=jts,jte
      do i=its,ite
        t8w(i,1,j)=t_phy(i,1,j)
        t8w(i,kte+1,j)=t_phy(i,kte,j)
      enddo
    enddo

    ! -- anthropogenic organic carbon
    do i=its,ite
      emis_anoc(i) = emi_ant_in(i,1)
    enddo

!---- Calculate PBLH and K-PBL based on virtual potential temperature profile
!---- First calculate THETAV
    do j = jts,jte
    do i = its,ite
    do k = kts,kte
       THETA = t_phy(i,k,j) * (1.E5/p_phy(i,k,j))**0.286
       THETAV(i,k,j) = THETA * (1. + 0.61 * (moist(i,k,j,p_qv)))
    enddo
    enddo
    enddo
!---- Now use the UPP code to deterimine the height and level
    do i = its, ite
    do j = jts, jte
       if ( THETAV(i,kts+1,j) .lt. ( THETAV(i,kts,j) + delta_theta4gust) ) then
          do k = kts+1, kte
             k1 = k
!--- give theta-v at the sfc a 0.5K boost in the PBLH definition
             if ( THETAV(i,kts+k-1,j) .gt. ( THETAV(i,kts,j) + delta_theta4gust) ) then
                exit
             endif
          enddo
          kpbl_thetav(i,j) = k1
       else
          kpbl_thetav(i,j) = kts + 1
       endif
   enddo
   enddo

!---- Calculate wind gust potential and HWP
    do i = its,ite
       SFCWIND          = sqrt(u10m(i)**2+v10m(i)**2)
       windgustpot(i,1) = SFCWIND
       if (kpbl_thetav(i,1)+1 .ge. kts+1 ) then
          do k=kts+1,int(kpbl_thetav(i,1))+1
             WIND = sqrt(us3d(i,k)**2+vs3d(i,k)**2)
             DELWIND = WIND - SFCWIND
             DZ = zmid(i,k,1) - oro(i)
             DELWIND = DELWIND*(1.0-MIN(0.5,DZ/2000.))
             windgustpot(i,1) = max(windgustpot(i,1),SFCWIND+DELWIND)
          enddo
       endif
    enddo
    hwp_local = 0.
    do i=its,ite
      wdgust=max(windgustpot(i,1),3.)
      snoweq=max((25.-snow_cpl(i))/25.,0.)
      hwp_local(i,1)=0.177*wdgust**0.97*max(t2m(i)-dpt2m(i),15.)**1.03*((1.-wetness(i))**0.4)*snoweq   ! Eric update   11/2023
    enddo
! Set paramters for ebb_dcycle option
    if (ebb_dcycle == 1 ) then
      if (hour_int .le. 24) then
          do j=jts,jte
           do i=its,ite
            ebu_in     (i,j) = smoke_RRFS(i,hour_int+1,1) ! smoke
            frp_in     (i,j) = smoke_RRFS(i,hour_int+1,2)*conv_frp ! frp
          ! These 2 arrays aren't needed for this option
          !  fire_end_hr(i,j) = 0.0
          !  hwp_day_avg(i,j) = 0.0
            ebb_smoke_in (i) = ebu_in(i,j)
           enddo
          enddo
      endif
    endif
    ! RAR: here we need to initialize various arrays in order to apply HWP to
    ! diurnal cycle
    ! if ebb_dcycle/=2 then those arrays=0, we need to read in temporal 
    if (ebb_dcycle == 2) then
      do i=its, ite
       do j=jts, jte 
         ebu_in      (i,j) = smoke2d_RRFS(i,1)!/86400.
         frp_in      (i,j) = smoke2d_RRFS(i,2)*conv_frp
         fire_end_hr (i,j) = smoke2d_RRFS(i,3)
         hwp_day_avg (i,j) = smoke2d_RRFS(i,4)
         ebb_smoke_in(i  ) = ebu_in(i,j)
       enddo
      enddo
    end if

    if (ktau==1) then
     do j=jts,jte
       do i=its,ite
    ! GFS_typedefs.F90 initializes this = 1, but should be OK to duplicate, RAR??
          fire_hist   (i,j) = 1.
          coef_bb_dc  (i,j) = 1.
          if (xlong(i,j)<230.) then
              peak_hr(i,j)= 0.0* 3600.    ! peak at 24 UTC, fires in Alaska
          elseif(xlong(i,j)<245.) then
              peak_hr(i,j)= 23.0* 3600.
          elseif (xlong(i,j)<260.) then
              peak_hr(i,j)= 22.0* 3600.    ! peak at 22 UTC, fires in the western US
          elseif (xlong(i,j)<275.) then
              peak_hr(i,j)= 21.0* 3600.
          elseif (xlong(i,j)<290.) then         ! peak at 20 UTC, fires in the eastern US
              peak_hr(i,j)= 20.0* 3600.
          else
              peak_hr(i,j)= 19.0* 3600.
          endif
       enddo
     enddo
    endif

    ! We will add a namelist variable, real :: flam_frac_global, RAR??
    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_smoke    )=max(epsilc,gq0(i,k,ntsmoke   ))
       chem(i,k,jts,p_dust_1   )=max(epsilc,gq0(i,k,ntdust    ))
       chem(i,k,jts,p_coarse_pm)=max(epsilc,gq0(i,k,ntcoarsepm))
     enddo
    enddo
 
  end subroutine rrfs_smoke_prep
  

!> @}
  end module rrfs_smoke_wrapper

!>\file rrfs_smoke_wrapper.F90
!! This file is CCPP driver of RRFS Smoke and Dust
!! Haiqin.Li@noaa.gov 02/2021

 module rrfs_smoke_wrapper

   use machine ,        only : kind_phys
   use rrfs_smoke_config
   use dust_data_mod
   use seas_mod,        only : gocart_seasalt_driver
   use dust_fengsha_mod,only : gocart_dust_fengsha_driver
   use plume_data_mod
   use module_plumerise1 !plume_rise_mod
   use module_add_emiss_burn
   use coarsepm_settling_mod
   use dep_dry_mod
   use module_wetdep_ls

   implicit none

   private

   public :: rrfs_smoke_wrapper_run

contains

!>\defgroup rrfs_smoke_wrapper rrfs-sd emission driver Module
!> \ingroup gsd_chem_group
!! This is the rrfs-sd emission driver Module
!! \section arg_table_rrfs_smoke_wrapper_run Argument Table
!! \htmlinclude rrfs_smoke_wrapper_run.html
!!
!>\section rrfs_smoke_wrapper rrfs-sd Scheme General Algorithm
!> @{
    subroutine rrfs_smoke_wrapper_run(im, kte, kme, ktau, dt, garea, land, jdate,          &
                   u10m, v10m, ustar, rlat, rlon, tskin, pb2d, t2m, dpt2m,                 &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum, w,                  &
                   nsoil, smc, vegtype, soiltyp, sigmaf, dswsfc, zorl,snow,                &
                   julian, idat, rain_cpl, rainc_cpl, exch, hf2d, g, pi, con_cp, con_rd,   &
                   dust12m_in, emi_in, smoke_RRFS, ntrac, qgrs, gq0, chem3d, tile_num,     &
                   ntsmoke, ntdust, ntcoarsepm, imp_physics, imp_physics_thompson,         &
                   nwfa, nifa, emanoc, emdust, emseas,                                     &
                   ebb_smoke_hr, frp_hr, frp_std_hr,                                       &
                   coef_bb, ebu_smoke,fhist, min_fplume, max_fplume, hwp, wetness,         &
                   smoke_ext, dust_ext, ndvel, ddvel_inout,rrfs_sd,                        &
                   dust_alpha_in, dust_gamma_in, fire_in,                                  &
                   seas_opt_in, dust_opt_in, drydep_opt_in, coarsepm_settling_in,          &
                   do_plumerise_in, plumerisefire_frq_in, addsmoke_flag_in,                &
                   wetdep_ls_opt_in,wetdep_ls_alpha_in,                                    &
                   smoke_forecast_in, aero_ind_fdb_in,dbg_opt_in,errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,nsoil,tile_num,jdate(8),idat(8)
    integer,        intent(in) :: ntrac, ntsmoke, ntdust, ntcoarsepm, ndvel
    real(kind_phys),intent(in) :: dt, julian, g, pi, con_cp, con_rd
    logical,        intent(in) :: aero_ind_fdb_in,dbg_opt_in
    integer,        intent(in) :: smoke_forecast_in

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(:), intent(in) :: land, vegtype, soiltyp
    real(kind_phys), dimension(:,:), intent(in) :: smc
    real(kind_phys), dimension(:,:,:), intent(in) :: dust12m_in
    real(kind_phys), dimension(:,:,:), intent(in) :: smoke_RRFS
    real(kind_phys), dimension(:,:), intent(in) :: emi_in
    real(kind_phys), dimension(:), intent(in) :: u10m, v10m, ustar, dswsfc,    &
                garea, rlat,rlon, tskin, pb2d, sigmaf, zorl, snow,             &
                rain_cpl, rainc_cpl, hf2d, t2m, dpt2m
    real(kind_phys), dimension(:,:), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(:,:), intent(in) :: phl3d, prl3d, tk3d,         &
                us3d, vs3d, spechum, exch, w
    real(kind_phys), dimension(:,:,:), intent(inout) :: qgrs, gq0
    real(kind_phys), dimension(:,:,:), intent(inout) :: chem3d
    real(kind_phys), dimension(:), intent(inout) :: emdust, emseas, emanoc
    real(kind_phys), dimension(:), intent(inout) :: ebb_smoke_hr, frp_hr, frp_std_hr
    real(kind_phys), dimension(:), intent(inout) :: coef_bb, fhist
    real(kind_phys), dimension(:,:), intent(inout) :: ebu_smoke
    real(kind_phys), dimension(:,:), intent(inout) :: fire_in
    real(kind_phys), dimension(:), intent(inout) :: max_fplume, min_fplume       
    real(kind_phys), dimension(:), intent(  out) :: hwp
    real(kind_phys), dimension(:,:), intent(out) :: smoke_ext, dust_ext
    real(kind_phys), dimension(:,:), intent(inout) :: nwfa, nifa
    real(kind_phys), dimension(:,:), intent(inout) :: ddvel_inout
    real (kind=kind_phys), dimension(:), intent(in) :: wetness
    integer, intent(in   ) :: imp_physics, imp_physics_thompson
    real (kind=kind_phys), intent(in) :: dust_alpha_in, dust_gamma_in, wetdep_ls_alpha_in
    integer,        intent(in) :: seas_opt_in, dust_opt_in, drydep_opt_in,        &
                                  coarsepm_settling_in, plumerisefire_frq_in,     &
                                  addsmoke_flag_in, wetdep_ls_opt_in
    logical, intent(in   ) :: do_plumerise_in, rrfs_sd
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: ebu
    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,  &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, vvel, zmid, exch_h

    real(kind_phys), dimension(ims:im, jms:jme) :: u10, v10, ust, tsk,            &
                     xland, xlat, xlong, dxy, pbl, hfx, rcav, rnav

!>- sea salt & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, 1, jms:jme, 1:num_emis_seas  )  :: emis_seas
    real(kind_phys), dimension(ims:im, jms:jme) :: seashelp

    integer :: ide, ime, ite, kde, julday

!>- dust & chemistry variables
    real(kind_phys), dimension(ims:im, jms:jme) :: ssm, rdrag, uthr, snowh  ! fengsha dust
    real(kind_phys), dimension(ims:im, jms:jme) :: vegfrac, rmol, swdown, znt, clayf, sandf
    real(kind_phys), dimension(ims:im, nsoil, jms:jme) :: smois
    real(kind_phys), dimension(ims:im, 1:1, jms:jme, 1:num_emis_dust) :: emis_dust
    integer,         dimension(ims:im, jms:jme) :: isltyp, ivgtyp

!>- plume variables
    ! -- buffers
    real(kind_phys), dimension(ims:im, jms:jme) :: ebu_in
    real(kind_phys), dimension(ims:im, jms:jme, num_frp_plume ) :: plume_frp
    real(kind_phys), dimension(ims:im, jms:jme )  :: coef_bb_dc, flam_frac,             &
                     fire_hist, peak_hr
    real(kind_phys), dimension(ims:im,kms:kme,jms:jme ) :: ext3d_smoke, ext3d_dust
    integer,         dimension(ims:im, jms:jme )  :: min_fplume2, max_fplume2
    logical :: call_fire
!>- optical variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: rel_hum
    real(kind_phys), dimension(ims:im, jms:jme, ndvel) :: ddvel

!>-- anthropogentic variables
    real(kind_phys), dimension(ims:im) :: emis_anoc
    real(kind_phys), dimension(ims:im, jms:jme, 1) :: sedim

    real(kind_phys) :: gmt

!> -- parameter to caluclate wfa&ifa (m)
    real(kind_phys), parameter :: mean_diameter1= 4.E-8, sigma1=1.8
    real(kind_phys), parameter :: mean_diameter2= 1.E-6, sigma2=1.8
    real(kind_phys), parameter :: kappa_oc      = 0.2
    real(kind_phys), parameter :: kappa_dust    = 0.04
    real(kind_phys) :: fact_wfa, fact_ifa
!> -- aerosol density (kg/m3)
    real(kind_phys), parameter :: density_dust= 2.6e+3, density_sulfate=1.8e+3
    real(kind_phys), parameter :: density_oc  = 1.4e+3, density_seasalt=2.2e+3

    real(kind_phys), dimension(im) :: daero_emis_wfa, daero_emis_ifa
!>-- local variables
    real(kind_phys), dimension(im) :: wdgust, snoweq
    integer :: current_month, current_hour, hour_int
    real(kind_phys) :: curr_secs
    real(kind_phys) :: factor, factor2, factor3
    integer :: nbegin, nv
    integer :: i, j, k, kp, n

    errmsg = ''
    errflg = 0

    if (.not. rrfs_sd) return

!>-- options to turn on/off sea-salt, dust, plume-rising
    seas_opt          = seas_opt_in
    dust_opt          = dust_opt_in
    drydep_opt        = drydep_opt_in
    do_plumerise      = do_plumerise_in
    plumerisefire_frq = plumerisefire_frq_in
    addsmoke_flag     = addsmoke_flag_in
    smoke_forecast    = smoke_forecast_in
    aero_ind_fdb      = aero_ind_fdb_in
    dbg_opt           = dbg_opt_in
    wetdep_ls_opt     = wetdep_ls_opt_in
    wetdep_ls_alpha   = wetdep_ls_alpha_in
    coarsepm_settling = coarsepm_settling_in

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    min_fplume2 = 0
    max_fplume2 = 0
    emis_seas   = 0.
    emis_dust   = 0.
    peak_hr     = 0.
    flam_frac   = 0.
    ext3d_smoke = 0.
    ext3d_dust  = 0.
    daero_emis_wfa = 0.
    daero_emis_ifa = 0.

    rcav = 0.
    rnav = 0.

    curr_secs = ktau * dt
    current_month=jdate(2)      ! needed for the dust input data
    current_hour =jdate(5)+1    ! =1 at 00Z
    hour_int=ktau*dt/3600.    ! hours since the simulation start
    gmt    = real(mod(idat(5)+hour_int,24))
    julday = int(julian)

    do nv=1,ndvel
    do i=its,ite
      ddvel(i,1,nv)=ddvel_inout(i,nv)
    enddo
    enddo

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
     rnav(i,1)=max((rain_cpl(i)-rainc_cpl(i))*1000., 0.) ! meter to mm
     coef_bb_dc(i,1) = coef_bb(i)
     fire_hist (i,1) = fhist  (i)
    enddo


    ! plumerise frequency in minutes set up by the namelist input
    call_fire       = (do_plumerise .and. (plumerisefire_frq > 0))
    if (call_fire) call_fire = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0) .or. (ktau == 2)
    
!>- get ready for chemistry run
    call rrfs_smoke_prep(                                               &
        current_month, current_hour, gmt,                               &
        u10m,v10m,ustar,land,garea,rlat,rlon,tskin,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,exch,w,            &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow,dust12m_in,emi_in,smoke_RRFS,                              &
        hf2d, pb2d, g, pi, hour_int,                                    &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,exch_h,                                                     &
        z_at_w,vvel,zmid,                                               &
        ntrac,gq0,                                                      &
        num_chem,num_moist,                                             &
        ntsmoke, ntdust,ntcoarsepm,                                     &
        moist,chem,plume_frp,ebu_in,                                    &
        ebb_smoke_hr, frp_hr, frp_std_hr, emis_anoc,                    &
        smois,ivgtyp,isltyp,vegfrac,rmol,swdown,znt,hfx,pbl,            &
        snowh,clayf,rdrag,sandf,ssm,uthr,rel_hum,                       &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

! Make this global, calculate at 1st time step only
!>-- for plumerise --

     do j=jts,jte
       do i=its,ite
         peak_hr(i,j)= fire_in(i,10)
       enddo
     enddo

     IF (ktau==1) THEN
     do j=jts,jte
       do i=its,ite
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
     ENDIF

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
       ebu(i,k,1)=ebu_smoke(i,k)
      enddo
      enddo
     ENDIF


!>- compute sea-salt
    ! -- compute sea salt (opt=2)
    if (seas_opt == 2) then
    call gocart_seasalt_driver(dt,rri,t_phy,                            &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,ust,p8w,tsk,              &
        xland,xlat,xlong,dxy,g,emis_seas,pi,                            &
        seashelp,num_emis_seas,num_chem,seas_opt,                       &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)
    endif

    !-- compute dust (opt=5)
    if (dust_opt==DUST_OPT_FENGSHA) then
       ! Set at compile time in dust_data_mod:
       dust_alpha = dust_alpha_in
       dust_gamma = dust_gamma_in
       call gocart_dust_fengsha_driver(dt,chem,rho_phy,smois,p8w,ssm,   &
            isltyp,vegfrac,snowh,xland,dxy,g,emis_dust,ust,znt,         &
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
    if (call_fire) then
        call ebu_driver (                                              &
                   flam_frac,ebu_in,ebu,                          &
                   t_phy,moist(:,:,:,p_qv),                            &
                   rho_phy,vvel,u_phy,v_phy,p_phy,                     &
                   z_at_w,zmid,g,con_cp,con_rd,                        &
                   plume_frp, min_fplume2, max_fplume2,                &   ! new approach
                   ids,ide, jds,jde, kds,kde,                          &
                   ims,ime, jms,jme, kms,kme,                          &
                   its,ite, jts,jte, kts,kte, errmsg, errflg           )
        if(errflg/=0) return
    end if

    ! -- add biomass burning emissions at every timestep
    if (addsmoke_flag == 1) then
    call add_emis_burn(dt,dz8w,rho_phy,rel_hum,chem,                 &
                       julday,gmt,xlat,xlong,                        &
                       ivgtyp, vegfrac, peak_hr,                     &   ! RAR
                       curr_secs,ebu,                                &
                       coef_bb_dc,fire_hist,ext3d_smoke,ext3d_dust,  &
                       rcav, rnav,swdown,smoke_forecast,             &
                       ids,ide, jds,jde, kds,kde,                    &
                       ims,ime, jms,jme, kms,kme,                    &
                       its,ite, jts,jte, kts,kte                     )
    endif

    !>-- compute coarsepm setting
    if (coarsepm_settling == 1) then
    call coarsepm_settling_driver(dt,t_phy,rel_hum,                  &
                                  chem(:,:,:,p_coarse_pm),           &
                                  rho_phy,dz8w,p8w,p_phy,sedim,      &
                                  dxy,g,1,                           &
                                  ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte          )
    endif
    !>-- compute dry deposition
    if (drydep_opt == 1) then

    call dry_dep_driver(rmol,ust,ndvel,ddvel,rel_hum,                &
       ids,ide, jds,jde, kds,kde,                                    &
       ims,ime, jms,jme, kms,kme,                                    &
       its,ite, jts,jte, kts,kte)

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
       call  wetdep_ls(dt,chem,rnav,moist,                      &
                     rho_phy,num_chem,num_moist,dz8w,vvel,      &
                     ids,ide, jds,jde, kds,kde,                 &
                     ims,ime, jms,jme, kms,kme,                 &
                     its,ite, jts,jte, kts,kte)
    endif

    do k=kts,kte
     do i=its,ite
       ebu_smoke(i,k)=ebu(i,k,1)
     enddo
    enddo


!---- diagnostic output of hourly wildfire potential (07/2021)
    hwp = 0.
    do i=its,ite
      wdgust(i)=max(1.68*sqrt(us3d(i,1)**2+vs3d(i,1)**2),3.)
      snoweq(i)=max((25.-snow(i))/25.,0.)
      hwp(i)=0.237*wdgust(i)**1.11*max(t2m(i)-dpt2m(i),15.)**0.92*((1.-wetness(i))**6.95)*snoweq(i) ! Eric 08/2022
    enddo
    
!---- diagnostic output of smoke & dust optical extinction (12/2021)
    do k=kts,kte
     do i=its,ite
       smoke_ext(i,k) = ext3d_smoke(i,k,1)
       dust_ext (i,k) = ext3d_dust (i,k,1)
     enddo
    enddo
!-------------------------------------
!---- put smoke stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntsmoke )  = min(5000.,max(epsilc,chem(i,k,1,p_smoke ))) 
       gq0(i,k,ntdust  )  = min(100.,max(epsilc,chem(i,k,1,p_dust_1)))
       gq0(i,k,ntcoarsepm)= min(1000.,max(epsilc,chem(i,k,1,p_coarse_pm)))
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
     emseas     (i) = emis_seas(i,1,1,1)*1.e+9   ! size bin 1 sea salt emission: ug/m2/s
     emdust     (i) = emis_dust(i,1,1,1) + emis_dust(i,1,1,2) +   &
                      emis_dust(i,1,1,3) + emis_dust(i,1,1,4) ! dust emission: ug/m2/s
     emanoc     (i) = emis_anoc (i)              ! anthropogenic organic carbon: ug/m2/s
     coef_bb    (i) = coef_bb_dc(i,1)
     fhist      (i) = fire_hist (i,1)
     min_fplume (i) = real(min_fplume2(i,1))
     max_fplume (i) = real(max_fplume2(i,1))
     emseas     (i) = sandf(i,1) ! sand for dust
     emanoc     (i) = uthr (i,1) ! u threshold for dust
    enddo

    do i = 1, im
      fire_in(i,10) =  peak_hr(i,1)
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

 subroutine rrfs_smoke_prep(                                            &
        current_month,current_hour,gmt,                                 &
        u10m,v10m,ustar,land,garea,rlat,rlon,ts2d,                      &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,exch,w,            &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow_cpl,dust12m_in,emi_in,smoke_RRFS,                          &
        hf2d, pb2d, g, pi, hour_int,                                    &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,exch_h,                                                     &
        z_at_w,vvel,zmid,                                               &
        ntrac,gq0,                                                      &
        num_chem, num_moist,                                            &
        ntsmoke, ntdust, ntcoarsepm,                                    &
        moist,chem,plume_frp,ebu_in,                                    &
        ebb_smoke_hr, frp_hr, frp_std_hr, emis_anoc,                    &
        smois,ivgtyp,isltyp,vegfrac,rmol,swdown,znt,hfx,pbl,            &
        snowh,clayf,rdrag,sandf,ssm,uthr,rel_hum,                       &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    integer, intent(in) :: current_month, current_hour, hour_int

    !FV3 input variables
    integer, intent(in) :: nsoil
    integer, dimension(ims:ime), intent(in) :: land, vegtype, soiltyp
    integer, intent(in) :: ntrac
    real(kind=kind_phys), intent(in) :: g, pi, gmt
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         u10m, v10m, ustar, garea, rlat, rlon, ts2d, sigmaf, dswsfc,       &
         zorl, snow_cpl, pb2d, hf2d
    real(kind=kind_phys), dimension(ims:ime, nsoil),   intent(in) :: smc 
    real(kind=kind_phys), dimension(ims:ime, 12, 5),   intent(in) :: dust12m_in
    real(kind=kind_phys), dimension(ims:ime, 24, 3),   intent(in) :: smoke_RRFS
    real(kind=kind_phys), dimension(ims:ime,     1),   intent(in) :: emi_in
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,exch,w
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !rrfs-sd variables
    integer,intent(in) ::  num_chem, num_moist, ntsmoke, ntdust, ntcoarsepm
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte


    real(kind_phys), dimension(ims:ime, jms:jme),intent(out) :: ebu_in
    real(kind_phys), dimension(ims:ime, jms:jme, num_frp_plume), intent(out) :: plume_frp
    
    integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel,               &
         zmid, exch_h, rel_hum
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         u10, v10, ust, tsk, xland, xlat, xlong, dxy, vegfrac, rmol, swdown, znt,      &
         pbl, hfx, snowh, clayf, rdrag, sandf, ssm, uthr
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, nsoil, jms:jme), intent(out) :: smois
    real(kind_phys), dimension(ims:ime), intent(inout) :: ebb_smoke_hr, frp_hr, frp_std_hr
    real(kind_phys), dimension(ims:ime), intent(inout) :: emis_anoc
   !real(kind_phys), dimension(ims:ime, jms:jme, num_plume_data) :: plume
    real(kind_phys), parameter :: conv_frp = 1.e+06_kind_phys  ! FRP conversion factor, MW to W
    real(kind_phys), parameter :: frpc  = 1._kind_phys          ! FRP conversion factor (Regional)

    ! -- local variables
    integer i,ip,j,k,kp,kk,kkp,nv,l,ll,n

    ! -- initialize fire emissions
    !plume          = 0._kind_phys
    plume_frp      = 0._kind_phys
    ebu_in         = 0._kind_phys
    ebb_smoke_hr   = 0._kind_phys
    emis_anoc      = 0._kind_phys
    frp_hr         = 0._kind_phys
    frp_std_hr     = 0._kind_phys

    ! -- initialize output arrays
    isltyp         = 0._kind_phys
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    vvel           = 0._kind_phys
    zmid           = 0._kind_phys
    exch_h         = 0._kind_phys
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
    rel_hum        = 0._kind_phys

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
     ivgtyp (i,1)=vegtype(i)
     isltyp (i,1)=soiltyp(i)
     vegfrac(i,1)=sigmaf (i)
    enddo
   
    rmol=0.

    do k=1,nsoil
     do j=jts,jte
      do i=its,ite
       smois(i,k,j)=smc(i,k)
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
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(i,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(i,kkp)*rri(i,k,j)/g 
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(i,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=gq0(i,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=gq0(i,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          !rel_hum(i,k,j) = min(0.95,spechum(i,kkp))
          rel_hum(i,k,j) = min(0.95, moist(i,k,j,1) / &
            (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
            (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          rel_hum(i,k,j) = max(0.1,rel_hum(i,k,j))
          !--
          zmid(i,k,j)=phl3d(i,kkp)/g
        enddo
      enddo
    enddo

    ! -- the imported atmospheric heat diffusivity is only available up to kte-1
     do k=kts,kte-1
       do i=its,ite
          exch_h(i,k,1)=exch(i,k)
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
      emis_anoc(i) = emi_in(i,1)
    enddo

    if (hour_int<24) then
        do j=jts,jte
         do i=its,ite
          ebb_smoke_hr(i)  = smoke_RRFS(i,hour_int+1,1) ! smoke
          frp_hr      (i)  = smoke_RRFS(i,hour_int+1,2) ! frp
          frp_std_hr  (i)  = smoke_RRFS(i,hour_int+1,3) ! std frp
          ebu_in    (i,j)  = ebb_smoke_hr(i)
          plume_frp(i,j,p_frp_hr ) = conv_frp* frp_hr      (i)
          plume_frp(i,j,p_frp_std) = conv_frp* frp_std_hr  (i)
         enddo
        enddo
    endif

    ! We will add a namelist variable, real :: flam_frac_global

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

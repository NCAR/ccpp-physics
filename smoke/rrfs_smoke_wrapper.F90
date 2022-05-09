!>\file rrfs_smoke_wrapper.F90
!! This file is CCPP RRFS smoke driver
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
   use dep_dry_mod
   use rrfs_smoke_data

   implicit none

   private

   public :: rrfs_smoke_wrapper_run

contains

!>\defgroup rrfs_smoke_wrapper GSD Chem emission driver Module  
!> \ingroup gsd_chem_group
!! This is the GSD Chem emission driver Module
!! \section arg_table_rrfs_smoke_wrapper_run Argument Table
!! \htmlinclude rrfs_smoke_wrapper_run.html
!!
!>\section rrfs_smoke_wrapper GSD Chemistry Scheme General Algorithm
!> @{
    subroutine rrfs_smoke_wrapper_run(im, kte, kme, ktau, dt, garea, land, jdate,      &
                   u10m, v10m, ustar, rlat, rlon, tskin, pb2d, t2m, dpt2m,             &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum, w,              &
                   nsoil, smc, vegtype, soiltyp, sigmaf, dswsfc, zorl,snow,            &
                   julian, idat, rain_cpl, rainc_cpl, exch, hf2d, g, pi, con_cp, con_rd,   &
                   dust12m_in, emi_in, smoke_GBBEPx, ntrac, qgrs, gq0, chem3d, tile_num,   &
                   ntsmoke, ntdust, imp_physics, imp_physics_thompson,                 &
                   nwfa, nifa, emanoc,                                                 &
                   emdust, emseas, ebb_smoke_hr, frp_hr, frp_std_hr,                   &
                   coef_bb, ebu_smoke,fhist, min_fplume, max_fplume, hwp,              &
                   smoke_ext, dust_ext,                                                &
                   seas_opt_in, dust_opt_in, biomass_burn_opt_in, drydep_opt_in,       &
                   do_plumerise_in, plumerisefire_frq_in, addsmoke_flag_in,            &
                   smoke_forecast_in, aero_ind_fdb_in,dbg_opt_in,errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,nsoil,tile_num,jdate(8),idat(8)
    integer,        intent(in) :: ntrac, ntsmoke, ntdust
    real(kind_phys),intent(in) :: dt, julian, g, pi, con_cp, con_rd
    logical,        intent(in) :: smoke_forecast_in,aero_ind_fdb_in,dbg_opt_in

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(:), intent(in) :: land, vegtype, soiltyp
    real(kind_phys), dimension(:,:), intent(in) :: smc
    real(kind_phys), dimension(:,:,:), intent(in) :: dust12m_in
    real(kind_phys), dimension(:,:,:), intent(in) :: smoke_GBBEPx
    real(kind_phys), dimension(:,:), intent(in) :: emi_in
    real(kind_phys), dimension(:), intent(in) :: u10m, v10m, ustar, dswsfc,      &
                garea, rlat,rlon, tskin, pb2d, sigmaf, zorl, snow,                &
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
    real(kind_phys), dimension(:), intent(inout) :: max_fplume, min_fplume       
    real(kind_phys), dimension(:), intent(  out) :: hwp
    real(kind_phys), dimension(:,:), intent(out) :: smoke_ext, dust_ext
    real(kind_phys), dimension(:,:), intent(inout) :: nwfa, nifa
    integer, intent(in   ) :: imp_physics, imp_physics_thompson
    integer,        intent(in) :: seas_opt_in, dust_opt_in, biomass_burn_opt_in,  &
                                  drydep_opt_in, plumerisefire_frq_in, addsmoke_flag_in
    logical, intent(in   ) :: do_plumerise_in
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
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  :: dry_fall
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
    real(kind_phys), dimension(ims:im,kms:kme,jms:jme ) :: aod3d_smoke, aod3d_dust
    integer,         dimension(ims:im, jms:jme )  :: min_fplume2, max_fplume2
    real(kind_phys) :: dtstep
    logical :: call_plume, scale_fire_emiss
!>- optical variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: rel_hum

!>-- anthropogentic variables
!   real(kind_phys), dimension(ims:im, kms:kemit, jms:jme, 1:num_emis_ant) :: emis_ant
    real(kind_phys), dimension(ims:im) :: emis_anoc

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: ac3, ahno3, anh3, asulf, cor3, h2oai, h2oaj, nu3
    real(kind_phys), dimension(ims:im, jms:jme) :: dep_vel_o3, e_co

    real(kind_phys) :: gmt
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

!> -- parameter to caluclate wfa&ifa (m)
    real(kind_phys), parameter :: mean_diameter1= 4.E-8, sigma1=1.8
    real(kind_phys), parameter :: mean_diameter2= 1.E-6, sigma2=1.8
    real(kind_phys), parameter :: kappa_oc      = 0.2
    real(kind_phys), parameter :: kappa_dust    = 0.04
    real(kind_phys) :: fact_wfa, fact_ifa
!> -- aerosol density (kg/m3)
    real(kind_phys), parameter :: density_dust= 2.6e+3, density_sulfate=1.8e+3
    real(kind_phys), parameter :: density_oc  = 1.4e+3, density_seasalt=2.2e+3

    real(kind_phys) :: daero_emis_wfa, daero_emis_ifa
!>-- local variables
    real(kind_phys), dimension(im) :: wdgust, snoweq
    integer :: current_month, current_hour
    real(kind_phys) :: curr_secs
    real(kind_phys) :: factor, factor2, factor3
    integer :: nbegin, nv, nvv
    integer :: i, j, jp, k, kp, n

    type(smoke_data), pointer :: data

    data => get_thread_smoke_data()

    errmsg = ''
    errflg = 0

!>-- options to turn on/off sea-salt, dust, plume-rising
    seas_opt          = seas_opt_in
    dust_opt          = dust_opt_in
    biomass_burn_opt  = biomass_burn_opt_in
    drydep_opt        = drydep_opt_in
    do_plumerise      = do_plumerise_in
    plumerisefire_frq = plumerisefire_frq_in
    addsmoke_flag     = addsmoke_flag_in
    smoke_forecast    = smoke_forecast_in
    aero_ind_fdb      = aero_ind_fdb_in
    dbg_opt           = dbg_opt_in

    !print*,'hli ktau',ktau
    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    h2oai = 0.
    h2oaj = 0.
    nu3   = 0.
    ac3   = 0.
    cor3  = 0.
    asulf = 0.
    ahno3 = 0.
    anh3  = 0.
    e_co  = 0.
    dep_vel_o3 = 0.

    min_fplume2 = 0
    max_fplume2 = 0
    emis_seas   = 0.
    emis_dust   = 0.
    peak_hr     = 0.
    flam_frac   = 0.
    aod3d_smoke = 0.
    aod3d_dust  = 0.

    rcav = 0.
    rnav = 0.

    curr_secs = ktau * dt
    current_month=jdate(2)
    current_hour =jdate(5)+1
    gmt    = real(idat(5))
    julday = int(julian)

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
     rnav(i,1)=max((rain_cpl(i)-rainc_cpl(i))*1000., 0.) ! meter to mm
     coef_bb_dc(i,1) = coef_bb(i)
     fire_hist (i,1) = fhist  (i)
    enddo


    ! plumerise frequency in minutes set up by the namelist input
    call_plume       = (biomass_burn_opt == BURN_OPT_ENABLE) .and. (plumerisefire_frq > 0)
    if (call_plume) &
       call_plume    = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0)         &
                        .or. (ktau == 2)
    
    !scale_fire_emiss = .false.

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

!>- get ready for chemistry run
    call rrfs_smoke_prep(                                               &
        ktau, current_month, current_hour,                              &
        u10m,v10m,ustar,land,garea,rlat,rlon,tskin,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,exch,w,            &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow,dust12m_in,emi_in,smoke_GBBEPx,                            &
        hf2d, pb2d, g, pi,                                              &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,exch_h,                                                     &
        z_at_w,vvel,zmid,                                               &
        ntrac,gq0,                                                      &
        num_chem, num_moist, ppm2ugkg,                                  &
        ntsmoke, ntdust,                                                &
        moist,chem,plume_frp,ebu_in,                                    &
        ebb_smoke_hr, frp_hr, frp_std_hr, emis_anoc,                    &
        smois,ivgtyp,isltyp,vegfrac,rmol,swdown,znt,hfx,pbl,            &
        snowh,clayf,rdrag,sandf,ssm,uthr,rel_hum,                       &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

! Make this global, calculate at 1st time step only
!>-- for plumerise --
!IF (ktau==1) THEN
     do j=jts,jte
       do i=its,ite
          if (xlong(i,j)<-130.) then
              peak_hr(i,j)= 0.0* 3600.    ! peak at 24 UTC, fires in Alaska
          elseif(xlong(i,j)<-115.) then
              peak_hr(i,j)= 23.0* 3600.
          elseif (xlong(i,j)<-100.) then
              peak_hr(i,j)= 22.0* 3600.    ! peak at 22 UTC, fires in the western US
          elseif (xlong(i,j)<-85.) then
              peak_hr(i,j)= 21.0* 3600.
          elseif (xlong(i,j)<-70.) then         ! peak at 20 UTC, fires in the eastern US
              peak_hr(i,j)= 20.0* 3600.
          else
              peak_hr(i,j)= 19.0* 3600.
          endif
       enddo
     enddo
!END IF

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
    ! -- compute sea salt
    if (seas_opt >= SEAS_OPT_DEFAULT) then
    call gocart_seasalt_driver(ktau,dt,rri,t_phy,moist,                 &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,ust,p8w,tsk,              &
        xland,xlat,xlong,dxy,g,emis_seas,pi,                            &
        seashelp,num_emis_seas,num_moist,num_chem,seas_opt,             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)
    endif

    !-- compute dust
    select case (dust_opt)
      case (DUST_OPT_FENGSHA)
       ! Set at compile time in dust_data_mod:
       call gocart_dust_fengsha_driver(data,dt,chem,rho_phy,smois,p8w,ssm,   &
            isltyp,vegfrac,snowh,xland,dxy,g,emis_dust,ust,znt,         &
            clayf,sandf,rdrag,uthr,                                     &
            num_emis_dust,num_moist,num_chem,nsoil,                     &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte)
    end select

    ! compute wild-fire plumes
    !-- to add a namelist option to turn on/off plume raising
    !--- replace plumerise_driver with HRRR-smoke 05/10/2021
    !-- /scratch2/BMC/ap-fc/Ravan/rapid-refresh/WRFV3.9/smoke
    ! Every hour (per namelist) the ebu_driver is called to calculate ebu, but
    ! the plumerise is controlled by the namelist option of plumerise_flag
    if (call_plume) then
!       WRITE(*,*) 'plumerise is called at ktau= ',ktau
        call ebu_driver (                                              &
                   data,flam_frac,ebu_in,ebu,                          &
                   t_phy,moist(:,:,:,p_qv),                            &
                   rho_phy,vvel,u_phy,v_phy,p_phy,                     &
                   z_at_w,zmid,ktau,g,con_cp,con_rd,                   &
                   plume_frp, min_fplume2, max_fplume2,                &   ! new approach
                   ids,ide, jds,jde, kds,kde,                          &
                   ims,ime, jms,jme, kms,kme,                          &
                   its,ite, jts,jte, kts,kte, errmsg, errflg           )
        if(errflg/=0) return
    end if

    ! -- add biomass burning emissions at every timestep
    if (addsmoke_flag == 1) then
    call add_emis_burn(data,dtstep,ktau,dz8w,rho_phy,rel_hum,chem,        &
                       julday,gmt,xlat,xlong,                        &
                       ivgtyp, vegfrac, peak_hr,                     &   ! RAR
                       curr_secs,ebu,                                &
                       coef_bb_dc,fire_hist,aod3d_smoke,aod3d_dust,  &
        !              scalar(ims,kms,jms,P_QNWFA),scalar(ims,kms,jms,P_QNIFA), ! &
                       rcav, rnav,swdown,smoke_forecast,             &
                       ids,ide, jds,jde, kds,kde,                    &
                       ims,ime, jms,jme, kms,kme,                    &
                       its,ite, jts,jte, kts,kte                     )
    endif
!       WRITE(*,*) 'after add_emis_burn  at ktau= ',ktau

    !>-- compute dry deposition
    if (drydep_opt == 1) then
    call dry_dep_driver(data,ktau,dt,julday,current_month,t_phy,p_phy,    &
       moist,p8w,rmol,rri,gmt,t8w,rcav,                              &
       chem,rho_phy,dz8w,exch_h,hfx,                                 &
       ivgtyp,tsk,swdown,vegfrac,pbl,ust,znt,zmid,z_at_w,            &
       xland,xlat,xlong,h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,        &
       anh3,dry_fall,dep_vel_o3,g,                                   &
       e_co,kemit,snowh,numgas,                                      &
       num_chem,num_moist,                                           &
       ids,ide, jds,jde, kds,kde,                                    &
       ims,ime, jms,jme, kms,kme,                                    &
       its,ite, jts,jte, kts,kte)
    endif
!    WRITE(*,*) 'dry depostion is done at ktau= ',ktau

    do k=kts,kte
     do i=its,ite
       ebu_smoke(i,k)=ebu(i,k,1)
     enddo
    enddo


!---- diagnostic output of hourly wildfire potential (07/2021)
    hwp = 0.
    do i=its,ite
      wdgust(i)=1.68*sqrt(us3d(i,1)**2+vs3d(i,1)**2)
      snoweq(i)=max((25.-snow(i)*1000.)/25.,0.)
      !hwp(i)=44.09*wdgust(i)**1.82*max(0.,t2m(i)-dpt2m(i))**0.61*max(0.,1.-smc(i,1))**14.0*snoweq(i)*sigmaf(i)
      hwp(i)=44.09*wdgust(i)**1.82*(t2m(i)-dpt2m(i))**0.61*(1.-smc(i,1))**14.0*snoweq(i)*sigmaf(i)
    enddo
    
!---- diagnostic output of smoke & dust optical extinction (12/2021)
    do k=kts,kte
     do i=its,ite
       smoke_ext(i,k) = aod3d_smoke(i,k,1) 
       dust_ext (i,k) = aod3d_dust (i,k,1) 
     enddo
    enddo
!-------------------------------------
!---- put smoke stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntsmoke )=ppm2ugkg(p_smoke ) * max(epsilc,chem(i,k,1,p_smoke)) !
       gq0(i,k,ntdust  )=ppm2ugkg(p_dust_1) * max(epsilc,chem(i,k,1,p_dust_1))
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntsmoke )= gq0(i,k,ntsmoke ) 
       qgrs(i,k,ntdust  )= gq0(i,k,ntdust  ) 
       chem3d(i,k,1     )= gq0(i,k,ntsmoke )
       chem3d(i,k,2     )= gq0(i,k,ntdust  )
     enddo
    enddo
!-------------------------------------
!-- to output for diagnostics
!    WRITE(*,*) 'rrfs nwfa/nifa 1 at ktau= ',ktau
    do i = 1, im
     emseas     (i) = emis_seas  (i,1,1,1)*1.e+9   ! size bin 1 sea salt emission: ug/m2/s
     emdust     (i) = emis_dust  (i,1,1,1)         ! size bin 1 dust emission    : ug/m2/s
     emanoc     (i) = emis_anoc   (i)              ! anthropogenic organic carbon: ug/m2/s
     coef_bb    (i) = coef_bb_dc (i,1)
     fhist      (i) = fire_hist  (i,1)
     min_fplume (i) = real(min_fplume2(i,1))
     max_fplume (i) = real(max_fplume2(i,1))
    enddo

!    WRITE(*,*) 'rrfs nwfa/nifa 2 at ktau= ',ktau
!-- to provide real aerosol emission for Thompson MP
    if (imp_physics == imp_physics_thompson .and. aero_ind_fdb) then
      fact_wfa = 1.e-9*6.0/pi*exp(4.5*log(sigma1)**2)/mean_diameter1**3
      fact_ifa = 1.e-9*6.0/pi*exp(4.5*log(sigma2)**2)/mean_diameter2**3

      do i = its, ite
       do k = kts, kte
        if (k==1)then
         daero_emis_wfa =(emanoc(i)+ebu_smoke(i,k))/density_oc + emseas(i)/density_seasalt
        else
         daero_emis_wfa = ebu_smoke(i,k)/density_oc 
        endif
         daero_emis_wfa = kappa_oc* daero_emis_wfa*fact_wfa*rri(i,k,1)/dz8w(i,k,1) ! consider  using dust tracer

         nwfa(i,k)      = nwfa(i,k) + daero_emis_wfa*dt
         nifa(i,k)      = gq0(i,k,ntdust)/density_dust*fact_ifa*kappa_dust  ! Check the formula

        if(land(i).eq.1)then
         nwfa(i,k)      = nwfa(i,k)*(1 - 0.10*dt/86400.)  !-- mimicking dry deposition
        else
         nwfa(i,k)      = nwfa(i,k)*(1 - 0.05*dt/86400.)  !-- mimicking dry deposition
        endif
       enddo
      enddo
    endif
!    WRITE(*,*) 'rrfs smoke wrapper is done at ktau= ',ktau

 end subroutine rrfs_smoke_wrapper_run

 subroutine rrfs_smoke_prep(                                            &
        ktau,current_month,current_hour,                                &
        u10m,v10m,ustar,land,garea,rlat,rlon,ts2d,                      &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,exch,w,            &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow_cpl,dust12m_in,emi_in,smoke_GBBEPx,                        &
        hf2d, pb2d, g, pi,                                              &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,exch_h,                                                     &
        z_at_w,vvel,zmid,                                               &
        ntrac,gq0,                                                      &
        num_chem, num_moist, ppm2ugkg,                                  &
        ntsmoke, ntdust,                                                &
       !num_emis_ant,                                                   &
       !emis_ant,                                                       &
        moist,chem,plume_frp,ebu_in,                                    &
        ebb_smoke_hr, frp_hr, frp_std_hr, emis_anoc,                    &
        smois,ivgtyp,isltyp,vegfrac,rmol,swdown,znt,hfx,pbl,            &
        snowh,clayf,rdrag,sandf,ssm,uthr,rel_hum,                       &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    integer, intent(in) :: ktau, current_month, current_hour

    !FV3 input variables
    integer, intent(in) :: nsoil
    integer, dimension(ims:ime), intent(in) :: land, vegtype, soiltyp
    integer, intent(in) :: ntrac
    real(kind=kind_phys), intent(in) :: g, pi
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         u10m, v10m, ustar, garea, rlat, rlon, ts2d, sigmaf, dswsfc,       &
         zorl, snow_cpl, pb2d, hf2d
    real(kind=kind_phys), dimension(ims:ime, nsoil),   intent(in) :: smc 
    real(kind=kind_phys), dimension(ims:ime, 12, 5),   intent(in) :: dust12m_in
    real(kind=kind_phys), dimension(ims:ime, 24, 3),   intent(in) :: smoke_GBBEPx
    real(kind=kind_phys), dimension(ims:ime,     1),   intent(in) :: emi_in
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,exch,w
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
   !integer,intent(in) ::  num_emis_ant
    integer,intent(in) ::  num_chem, num_moist, ntsmoke, ntdust
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte


   !real(kind_phys), dimension(ims:ime, kms:kemit, jms:jme, num_emis_ant), intent(inout) :: emis_ant
    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg
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
    integer i,ip,j,jp,k,kp,kk,kkp,nv,l,ll,n

    ! -- initialize fire emissions
    !plume          = 0._kind_phys
    plume_frp      = 0._kind_phys
    ebu_in         = 0._kind_phys
    ebb_smoke_hr   = 0._kind_phys
    emis_anoc      = 0._kind_phys

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

    !if (ktau <= 1) then
    !  emis_ant = 0.
    ! !emis_vol = 0.
    !end if

    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,1)/g)
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=abs(ph3d(ip,kp+1)-ph3d(ip,kp))/g
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          p8w(i,k,j)=pr3d(ip,kp)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(ip,kkp)
          p_phy(i,k,j)=prl3d(ip,kkp)
          u_phy(i,k,j)=us3d(ip,kkp)
          v_phy(i,k,j)=vs3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(ip,kkp)*rri(i,k,j)/g 
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          rel_hum(i,k,j) = .95
          rel_hum(i,k,j) = MIN( .95, moist(i,k,j,1) / &
            (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
            (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          rel_hum(i,k,j)=max(0.1,rel_hum(i,k,j))
          !--
          zmid(i,k,j)=phl3d(ip,kkp)/g
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

    ! -- only used in phtolysis....
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

 !   select case (plumerise_flag)
 !     case (FIRE_OPT_GBBEPx)
        do j=jts,jte
         do i=its,ite
          ebb_smoke_hr(i)  = smoke_GBBEPx(i,current_hour,1) ! smoke
          frp_hr      (i)  = smoke_GBBEPx(i,current_hour,2) ! frp
          frp_std_hr  (i)  = smoke_GBBEPx(i,current_hour,3) ! std frp
          ebu_in    (i,j)  = ebb_smoke_hr(i)
          plume_frp(i,j,p_frp_hr     ) = conv_frp* frp_hr  (i)
          plume_frp(i,j,p_frp_std    ) = conv_frp* frp_std_hr  (i)
         enddo
        enddo
 !     case default
 !   end select

    ! We will add a namelist variable, real :: flam_frac_global

    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_smoke )=max(epsilc,gq0(i,k,ntsmoke )/ppm2ugkg(p_smoke))
       chem(i,k,jts,p_dust_1)=max(epsilc,gq0(i,k,ntdust  )/ppm2ugkg(p_dust_1))
     enddo
    enddo
 


  end subroutine rrfs_smoke_prep

!> @}
  end module rrfs_smoke_wrapper

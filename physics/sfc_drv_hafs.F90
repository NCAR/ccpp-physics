!>  \file sfc_drv_hafs.f
!!  This file contains the Noah land surface scheme driver.

!> This module contains the CCPP-compliant Noah land surface scheme driver for the hurricane application.
      module lsm_noah_hafs

      implicit none

      private

      public :: lsm_noah_hafs_init, lsm_noah_hafs_run, lsm_noah_hafs_finalize

      contains

!>\ingroup Noah_LSM_hafs
!! \section arg_table_lsm_noah_hafs_init Argument Table
!! \htmlinclude lsm_noah_hafs_init.html
!!
      subroutine lsm_noah_hafs_init(lsm, lsm_noah_hafs, restart, veg_data_choice, soil_data_choice, ialb, ncol, nsoil, vtype, snoalb, isurban, sthick, errmsg, errflg)

      use machine, only : kind_phys
      use module_sf_noahlsm,  only: maxalb
      
      implicit none
      
      integer,              intent(in)  :: lsm, lsm_noah_hafs, veg_data_choice, soil_data_choice, ialb, ncol, nsoil
      real(kind=kind_phys), dimension(ncol), intent(in)  :: vtype
      logical,              intent(in)  :: restart
      
      integer, intent(inout) :: isurban
      real(kind=kind_phys), dimension(ncol), intent(inout) :: snoalb
      real(kind=kind_phys), dimension(nsoil), intent(inout) :: sthick
      
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg
      
      ! Local variables
      
      character(len=256)                  :: mminlu, mminsl
      
      integer :: i, k
      !integer :: ite, ide, itf
      !integer               :: ids,ide, jds,jde, kds,kde, &
      !                         ims,ime, jms,jme, kms,kme, &
      !                         its,ite, jts,jte, kts,kte
      
      real(kind=kind_phys), parameter, dimension(4) :: zsoil = (/ -0.1,-0.4,-1.0,-2.0/) !what if nsoil /= 4?  
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      if (lsm/=lsm_noah_hafs) then
        write(errmsg,'(*(a))') "Logic error: namelist choice of LSM is different from NOAH HAFS"
        errflg = 1
        return
      end if
      
      select case (veg_data_choice)
       case (0)
         mminlu = 'USGS'
         isurban = 1
       case (1)
         mminlu = 'MODIFIED_IGBP_MODIS_NOAH'
         isurban = 13
       case (3)
         mminlu = 'NLCD40'
         isurban = 13
       case (4)
         mminlu = 'USGS-RUC'
         isurban = 1
       case (5)
         mminlu = 'MODI-RUC'
         isurban = 13
       case default
         errmsg = 'The value of the ivegsrc physics namelist parameter is incompatible with this version of NOAH LSM'
         errflg = 1
         return
      end select
      
      select case (soil_data_choice)
       case (1)
         mminsl = 'STAS'
       case (2)
         mminsl = 'STAS-RUC'
       case default
         errmsg = 'The value of the isot physics namelist parameter is incompatible with this version of NOAH LSM'
         errflg = 1
         return
      end select
      
      call soil_veg_gen_parm(trim(mminlu), trim(mminsl), errmsg, errflg)
      
      ! Set internal dimensions
      ! ids = 1
      ! ims = 1
      ! its = 1
      ! ide = ncol
      ! ime = ncol
      ! ite = ncol
      ! jds = 1
      ! jms = 1
      ! jts = 1
      ! jde = 1
      ! jme = 1
      ! jte = 1
      ! kds = 1
      ! kms = 1
      ! kts = 1
      ! kde = nlev
      ! kme = nlev
      ! kte = nlev
      
      if (.not. restart) then
        do i = 1, ncol      
          if(ialb == 0) then
             snoalb(i) = maxalb(int(0.5 + vtype(i)))*0.01
          endif
        end do
      end if
      
      sthick(1) = - zsoil(1)
      do k = 2, nsoil
        sthick(k) = zsoil(k-1) - zsoil(k)
      enddo
      
      end subroutine lsm_noah_hafs_init


!! \section arg_table_lsm_noah_hafs_finalize Argument Table
!! \htmlinclude lsm_noah_hafs_finalize.html
!!
      subroutine lsm_noah_hafs_finalize(errmsg, errflg)

      implicit none

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      end subroutine lsm_noah_hafs_finalize


!>\defgroup Noah_LSM_hafs Noah LSM Model for the hurricane application
!! \section arg_table_lsm_noah_hafs_run Argument Table
!! \htmlinclude lsm_noah_hafs_run.html
!!
!> \section general_noah_hafs_drv GFS sfc_drv General Algorithm
!>  @{
      subroutine lsm_noah_hafs_run                                      &
     &     (im, land, flag_iter, srflag, isurban, opt_thcnd, dt, rhowater, &
            eps, epsm1, cp, rd, rvrdm1, sigma, cph2o, cpice, lsubf, zlvl, nsoil, &
            sthick, lwdn, soldn, solnet, sfcprs, tprcp, sfctmp, q1, prslki, &
            vegtyp, soiltyp, slopetyp, shmin, shmax, snoalb, tbot, wind, shdfac, &
            albbrd, z0brd, z0k, emissi, canopy, t1, stc, smc, swc, snwdph, sneqv, &
            ch, ribb, eta_kinematic, shflx, embrd, ec, edir, ett, esnow, etp, &
            ssoil, snohf, sncovr, snowc, runoff, drain, stm, qsurf, smcwlt, smcref, errmsg, errflg &
     &     )
!
      use machine , only : kind_phys
      use module_sf_noahlsm, only: sflx, lutype, sltype
      use funcphys, only : fpvs

      implicit none
      
      integer,                             intent(in) :: im, isurban, opt_thcnd, nsoil
      real(kind=kind_phys),                intent(in) :: dt, rhowater, eps, epsm1, cp, rd, rvrdm1, sigma, cph2o, cpice, lsubf

      integer, dimension(im),              intent(in) :: vegtyp, soiltyp, slopetyp
      logical, dimension(im),              intent(in) :: flag_iter, land
      real(kind=kind_phys), dimension(im), intent(in) :: srflag, zlvl, lwdn, soldn, solnet, sfcprs, tprcp, sfctmp, q1, prslki, shmin, shmax, snoalb, tbot, wind
      real(kind=kind_phys), dimension(nsoil), intent(in) :: sthick

      real(kind=kind_phys), dimension(im), intent(inout) :: shdfac, albbrd, z0brd, z0k, emissi, canopy, t1, snwdph, sneqv, ch, ribb
      real(kind=kind_phys), dimension(im,nsoil), intent(inout) :: stc, smc, swc
      
      real(kind=kind_phys), dimension(im), intent(out) :: embrd, eta_kinematic, shflx, ec, edir, ett, esnow, etp, ssoil, snohf, sncovr, snowc, runoff, drain, stm, qsurf, smcwlt, smcref

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!     local Variables
      integer :: i, k
      logical, parameter :: local = .false.
      logical, parameter :: ua_phys = .false.
      logical, parameter :: rdlai2d = .false.
      logical, parameter :: usemonalb = .true. !Biswas says true for HWRF
      real(kind=kind_phys), parameter :: aoasis = 1.0 !hard-coded to 1 in module_sf_noahdrv or set to value from urban module?
      integer, parameter :: fasdas = 0 ! = 1 if using "flux-adjusting surface data assimilation system"
      real(kind=kind_phys) :: prcp, rho, qs1, q2k, th2, dqsdt2, dummy, cmc, snowhk, chk, sheat, flx1, flx2, flx3, runoff1, runoff2, soilm, smcmax 
      
      
      integer :: nroot
      real(kind=kind_phys) :: albedok, eta, fdown, drip, dew, beta, snomlt, runoff3, rc, pc, rsmin, xlai, rcs, rct, rcq, rcsoil, soilw, snotime1, smcdry
      real (kind=kind_phys), dimension(nsoil) :: et, smav
      real(kind=kind_phys) :: sfcheadrt, infxsrt, etpnd1 !doesn't appear to be used unless WRF_HYDRO preprocessor directive is defined and no documentation
      real(kind=kind_phys) :: xsda_qfx, hfx_phy, qfx_phy, xqnorm, hcpct_fasdas !only used if fasdas = 1 ("flux-adjusting surface data assimilation system")
      
      !GJF: 
      ! albedok is an output from sflx (but is not sent outside of the scheme in GFS) and is defined as:
      !     SURFACE ALBEDO INCLUDING SNOW EFFECT (UNITLESS FRACTION)
      !                =SNOW-FREE ALBEDO (ALB) WHEN SNEQV=0, OR
      !                =FCT(MSNOALB,ALB,VEGTYP,SHDFAC,SHDMIN) WHEN SNEQV>0
      !     if needed by other schemes or diagnostics, one needs to add it to the host model and CCPP metadata; could also just pass in dummy argument
      ! eta is an output from sflx, but eta_kinematic is what is passed out
      ! fdown is an output from sflx (but is not sent outside of the scheme in GFS) and is defined as:
      !     Radiation forcing at the surface (W m-2) = SOLDN*(1-alb)+LWDN
      ! et is an output from sflx (but is not sent outside of the scheme in GFS) and is defined as:
      !     plant transpiration from a particular root (soil) layer (W m-2)
      ! drip is an output from sflx (but is not sent outside of the scheme in GFS) and is defined as:
      !     through-fall of precip and/or dew in excess of canopy water-holding capacity (m)
      ! dew is an output from sflx (but is not sent outside of the scheme in GFS) and is defined as:
      !     dewfall (or frostfall for t<273.15) (m)
      ! beta is an output from sflx (but is not sent outside of the scheme in GFS) and is defined as:
      !     ratio of actual/potential evap (dimensionless)
      ! snomlt is an output from sflx (but is not sent outside of the scheme in GFS) and is defined as:
      !     snow melt (m) (water equivalent)
      ! rc is an output from sflx (but is not sent outside of the scheme in GFS) and is defined as:
      !     canopy resistance (s m-1)
      ! pc : plant coefficient (unitless fraction, 0-1) where pc*etp = actual transp
      ! rsmin : minimum canopy resistance (s m-1)
      ! xlai: leaf area index (dimensionless)
      ! rcs: incoming solar rc factor (dimensionless)
      ! rct: air temperature rc factor (dimensionless)
      ! rcq: atmos vapor pressure deficit rc factor (dimensionless)
      ! rcsoil: soil moisture rc factor (dimensionless)
      ! soilw: available soil moisture in root zone (unitless fraction between smcwlt and smcmax)
      ! smav: soil moisture availability for each layer, as a fraction between smcwlt and smcmax.
      ! snotime1: no documentation in module_sf_noahlsm.F, but described as "initial number of timesteps since last snowfall" in module_sf_noahdrv.F; related to CCPP nondimensional_snow_age for NoahMP? Since inout, need to initialize here?
      ! smcdry: dry soil moisture threshold where direct evap frm top layer ends (volumetric)
      ! nroot: number of root layers, a function of veg type, determined in subroutine redprm.
      
      !variables associated with UA_PHYS (not used for now)
      real(kind=kind_phys) :: flx4, fvb, fbur, fgsn
      
      REAL, PARAMETER  :: A2=17.67,A3=273.15,A4=29.65,   &
                          A23M4=A2*(A3-A4)
      
!> - Initialize CCPP error handling variables

      errmsg = ''
      errflg = 0
      
      do i=1, im
        if (flag_iter(i) .and. land(i)) then
          !GJF: module_sf_noahdrv.F from WRF has hardcoded slopetyp = 1; why? replicate here?
          !GJF: shdfac is zeroed out for particular combinations of vegetation table source and vegetation types; replicate here?
          
          !GJF: could potentially pass in pre-calculated rates instead of calculating here
          prcp = rhowater * tprcp(i) / dt
          
          !GJF: The GFS version of NOAH prepares the specific humidity in sfc_drv.f as follows:
          q2k   = max(q1(i), 1.e-8)
          rho   = sfcprs(i) / (rd*t1(i)*(1.0+rvrdm1*q2k))
          qs1   = fpvs( sfctmp(i) )
          qs1   = max(eps*qs1 / (sfcprs(i)+epsm1*qs1), 1.e-8)
          q2k   = min(qs1, q2k)
          
          !GJF: could potentially pass in pre-calcualted potential temperature if other schemes also need it (to avoid redundant calculation)
          th2 = t1(i) * prslki(i)
          
          !GJF: module_sf_noahdrv.F from WRF modifies dqsdt2 if the surface has snow.
          dqsdt2=qs1*a23m4/(sfctmp(i)-a4)**2
          
          !GJF: convert canopy moisture from kg m-2 to m
          canopy(i) = max(canopy(i), 0.0) !check for positive values in sfc_drv.f
          cmc = canopy(i)/rhowater
          
          !GJF: snow depth passed in to NOAH is conditionally modified differently in GFS and WRF:
          if ( (sneqv(i) /= 0.0 .and. snwdph(i) == 0.) .or. (snwdph(i) <= sneqv(i)) ) then
            snowhk = 5.*sneqv(i)
          endif
          !GJF: GFS version:
          ! if (sneqv(i) /= 0.0 .and. snwdph(i) == 0.0) then
          !   snowhk = 10.0 * sneqv(i)
          ! endif
          
          !GJF: calculate conductance from surface exchange coefficient
          chk = ch(i)  * wind(i)
          
          call sflx (i, 1, srflag(i), &
                   isurban, dt, zlvl(i), nsoil, sthick,              &    !c
                   local,                                            &    !L
                   lutype, sltype,                                   &    !CL
                   lwdn(i),soldn(i),solnet(i),sfcprs(i),prcp,sfctmp(i),q2k,dummy,&    !F
                   dummy,dummy, dummy,                               &    !F
                   th2,qs1,dqsdt2,                                   &    !I
                   vegtyp(i),soiltyp(i),slopetyp(i),shdfac(i),shmin(i),shmax(i),       &    !I
                   albbrd(i), snoalb(i), tbot(i), z0brd(i), z0k(i), emissi(i), embrd(i),  &    !S
                   cmc,t1(i),stc(i,:),smc(i,:),swc(i,:),snowhk,sneqv(i),albedok,chk,dummy,&    !H
                   cp, rd, sigma, cph2o, cpice, lsubf,&
                   eta,sheat,eta_kinematic(i),fdown,                   &    !O
                   ec(i),edir(i),et,ett(i),esnow(i),drip,dew,                    &    !O
                   beta,etp(i),ssoil(i),                                   &    !O
                   flx1,flx2,flx3,                                   &    !O
                   flx4,fvb,fbur,fgsn,ua_phys,                       &    !UA 
                   snomlt,sncovr(i),                                    &    !O
                   runoff1,runoff2,runoff3,                          &    !O
                   rc,pc,rsmin,xlai,rcs,rct,rcq,rcsoil,              &    !O
                   soilw,soilm,qsurf(i),smav,                              &    !D
                   rdlai2d,usemonalb,                                &
                   snotime1,                                         &
                   ribb(i),                                             &
                   smcwlt(i),smcdry,smcref(i),smcmax,nroot,                &
                   sfcheadrt,                                   &    !I
                   infxsrt,etpnd1,opt_thcnd,aoasis,              &    !O
                   xsda_qfx, hfx_phy, qfx_phy, xqnorm, fasdas, hcpct_fasdas,     &   ! fasdas
                   errflg, errmsg)
          if (errflg > 0) return
          !set fasdas = 0; all other vars can be dummy, I think
          
          canopy(i) = cmc*rhowater
          snwdph(i) = snowhk
           
          shflx(i) = sheat / (cp*rho)
          
          !aggregating several outputs into one like GFS sfc_drv.F
          snohf(i) = flx1 + flx2 + flx3
          
          snowc(i) = sncovr(i) !GJF: redundant?
          
          !convert from m s-1 to kg m-2 s-1 by multiplying by rhowater
          runoff(i) = runoff1 * rhowater
          drain(i) = runoff2 * rhowater
          stm(i) = soilm * rhowater
          
          !wet1(i) = smc(i,1) / smcmax !Sarah Lu added 09/09/2010 (for GOCART)
          
        endif
      end do
      
      
      end subroutine lsm_noah_hafs_run
      
      subroutine soil_veg_gen_parm( mminlu, mminsl, errmsg, errflg)
        !use namelist_soilveg_hafs
        use module_sf_noahlsm,  only: shdtbl, nrotbl, rstbl, rgltbl, hstbl, snuptbl, & ! begin land use / vegetation variables
                                      maxalb, laimintbl, laimaxtbl, z0mintbl, z0maxtbl, &
                                      albedomintbl, albedomaxtbl, ztopvtbl,zbotvtbl, &
                                      emissmintbl, emissmaxtbl, topt_data, cmcmax_data, &
                                      cfactr_data, rsmax_data, bare, natural, &
                                      low_density_residential, high_density_residential, &
                                      high_intensity_industrial, lucats, lutype, &  !end land use / vegetation variables
                                      bb,drysmc,f11,                           & ! begin soil variables
                                      maxsmc, refsmc,satpsi,satdk,satdw, wltsmc,qtz,&
                                      slcats, sltype, &                         ! end soil variables
                                      slope_data, sbeta_data,fxexp_data,csoil_data,salp_data,refdk_data,           & ! begin NOAH "general" variables
                                      refkdt_data,frzk_data,zbot_data,  smlow_data,smhigh_data,        &
                                      czil_data, lvcoef_data, slpcats ! end NOAH "general" variables
        implicit none

        character(len=*), intent(in) :: mminlu, mminsl
        character(len=*), intent(inout) :: errmsg
        integer,          intent(inout) :: errflg
        
        integer :: lumatch, iindex, lc, num_slope, iunit_noah
        integer :: ierr
        integer , parameter :: open_ok = 0
        logical :: opened

        character*128 :: mess , message
        character*256 :: a_string
        integer , parameter :: loop_max   = 10
        integer             :: loop_count, i
        
!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!                 Z0: Roughness length (m)
!             SHDFAC: Green vegetation fraction (in percentage)
!  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
!          the monthly green vegetation data
!             CMXTBL: MAX CNPY Capacity (m)
!             NROTBL: Rooting depth (layer)
!              RSMIN: Mimimum stomatal resistance (s m-1)
!              RSMAX: Max. stomatal resistance (s m-1)
!                RGL: Parameters used in radiation stress function
!                 HS: Parameter used in vapor pressure deficit functio
!               TOPT: Optimum transpiration air temperature. (K)
!             CMCMAX: Maximum canopy water capacity
!             CFACTR: Parameter used in the canopy inteception calculati
!               SNUP: Threshold snow depth (in water equivalent m) that
!                     implies 100% snow cover
!                LAI: Leaf area index (dimensionless)
!             MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
!
      iunit_noah = -1
      do i = 20,99
        inquire ( i , opened = opened )
        if ( .not. opened ) then
          iunit_noah = i
          exit
        endif
      enddo

      if ( iunit_noah < 0 ) then
        errflg = 1
        errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: '//   &
                 'can not find unused fortran unit to read.'
        return
      endif
      
      open(iunit_noah, file='VEGPARM.TBL',form='formatted',status='old',iostat=ierr)
      if(ierr .ne. open_ok ) then
        errflg = 1
        errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: failure opening VEGPARM.TBL'
        return
      end if
      
      lumatch=0

      loop_count = 0
      read (iunit_noah,fmt='(a)',end=2002) a_string
      find_lutype : do while (lumatch == 0)
         read (iunit_noah,*,end=2002)lutype
         read (iunit_noah,*)lucats,iindex
         if(lutype.eq.mminlu)then
            !write( mess , * ) 'landuse type = ' // trim ( lutype ) // ' found', lucats,' categories'
            !call wrf_message( mess )
            lumatch=1
         else
            loop_count = loop_count+1
            !call wrf_message ( "skipping over lutype = " // trim ( lutype ) )
            find_vegetation_parameter_flag : do
               read (iunit_noah,fmt='(a)', end=2002) a_string
               if ( a_string(1:21) .eq. 'Vegetation Parameters' ) then
                  exit find_vegetation_parameter_flag
               else if ( loop_count .ge. loop_max ) then
                  errflg = 1
                  errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: too many loops in VEGPARM.TBL'
                  return
               endif
            enddo find_vegetation_parameter_flag
         endif
      enddo find_lutype
      
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
      if ( size(shdtbl)       < lucats .or. &
           size(nrotbl)       < lucats .or. &
           size(rstbl)        < lucats .or. &
           size(rgltbl)       < lucats .or. &
           size(hstbl)        < lucats .or. &
           size(snuptbl)      < lucats .or. &
           size(maxalb)       < lucats .or. &
           size(laimintbl)    < lucats .or. &
           size(laimaxtbl)    < lucats .or. &
           size(z0mintbl)     < lucats .or. &
           size(z0maxtbl)     < lucats .or. &
           size(albedomintbl) < lucats .or. &
           size(albedomaxtbl) < lucats .or. &
           size(ztopvtbl) < lucats .or. &
           size(zbotvtbl) < lucats .or. &
           size(emissmintbl ) < lucats .or. &
           size(emissmaxtbl ) < lucats ) then
         errflg = 1
         errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: table sizes too small for value of lucats'
         return
      endif

      if(lutype.eq.mminlu)then
        do lc=1,lucats
          read (iunit_noah,*)iindex,shdtbl(lc),                        &
                            nrotbl(lc),rstbl(lc),rgltbl(lc),hstbl(lc), &
                            snuptbl(lc),maxalb(lc), laimintbl(lc),     &
                            laimaxtbl(lc),emissmintbl(lc),             &
                            emissmaxtbl(lc), albedomintbl(lc),         &
                            albedomaxtbl(lc), z0mintbl(lc), z0maxtbl(lc),&
                            ztopvtbl(lc), zbotvtbl(lc)
        enddo
        
        read (iunit_noah,*)
        read (iunit_noah,*)topt_data
        read (iunit_noah,*)
        read (iunit_noah,*)cmcmax_data
        read (iunit_noah,*)
        read (iunit_noah,*)cfactr_data
        read (iunit_noah,*)
        read (iunit_noah,*)rsmax_data
        read (iunit_noah,*)
        read (iunit_noah,*)bare
        read (iunit_noah,*)
        read (iunit_noah,*)natural
        read (iunit_noah,*)
        read (iunit_noah,*)
        read (iunit_noah,fmt='(a)') a_string
        if ( a_string(1:21) .eq. 'Vegetation Parameters' ) then
           errflg = 1
           errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: expected low and high density residential, and high density industrial information in VEGPARM.TBL'
           return
        endif
        read (iunit_noah,*)low_density_residential
        read (iunit_noah,*)
        read (iunit_noah,*)high_density_residential
        read (iunit_noah,*)
        read (iunit_noah,*)high_intensity_industrial
      endif

2002   continue

      close (iunit_noah)
      if (lumatch == 0) then
         errflg = 1
         errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: land use dataset '//mminlu//' not found in VEGPARM.TBL.'
         return
      endif
      
      
      !CALL wrf_dm_bcast_string  ( LUTYPE  , 4 )
      !CALL wrf_dm_bcast_integer ( LUCATS  , 1 )
      !CALL wrf_dm_bcast_integer ( IINDEX  , 1 )
      !CALL wrf_dm_bcast_integer ( LUMATCH , 1 )
      !CALL wrf_dm_bcast_real    ( SHDTBL  , NLUS )
      !CALL wrf_dm_bcast_real    ( NROTBL  , NLUS )
      !CALL wrf_dm_bcast_real    ( RSTBL   , NLUS )
      !CALL wrf_dm_bcast_real    ( RGLTBL  , NLUS )
      !CALL wrf_dm_bcast_real    ( HSTBL   , NLUS )
      !CALL wrf_dm_bcast_real    ( SNUPTBL , NLUS )
      !CALL wrf_dm_bcast_real    ( LAIMINTBL    , NLUS )
      !CALL wrf_dm_bcast_real    ( LAIMAXTBL    , NLUS )
      !CALL wrf_dm_bcast_real    ( Z0MINTBL     , NLUS )
      !CALL wrf_dm_bcast_real    ( Z0MAXTBL     , NLUS )
      !CALL wrf_dm_bcast_real    ( EMISSMINTBL  , NLUS )
      !CALL wrf_dm_bcast_real    ( EMISSMAXTBL  , NLUS )
      !CALL wrf_dm_bcast_real    ( ALBEDOMINTBL , NLUS )
      !CALL wrf_dm_bcast_real    ( ALBEDOMAXTBL , NLUS )
      !CALL wrf_dm_bcast_real    ( ZTOPVTBL , NLUS )
      !CALL wrf_dm_bcast_real    ( ZBOTVTBL , NLUS )
      !CALL wrf_dm_bcast_real    ( MAXALB  , NLUS )
      !CALL wrf_dm_bcast_real    ( TOPT_DATA    , 1 )
      !CALL wrf_dm_bcast_real    ( CMCMAX_DATA  , 1 )
      !CALL wrf_dm_bcast_real    ( CFACTR_DATA  , 1 )
      !CALL wrf_dm_bcast_real    ( RSMAX_DATA  , 1 )
      !CALL wrf_dm_bcast_integer ( BARE    , 1 )
      !CALL wrf_dm_bcast_integer ( NATURAL , 1 )
      !CALL wrf_dm_bcast_integer ( LOW_DENSITY_RESIDENTIAL , 1 )
      !CALL wrf_dm_bcast_integer ( HIGH_DENSITY_RESIDENTIAL , 1 )
      !CALL wrf_dm_bcast_integer ( HIGH_INTENSITY_INDUSTRIAL , 1 )
      
!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
            
      open(iunit_noah, file='SOILPARM.TBL',form='formatted',status='old',iostat=ierr)
      if(ierr .ne. open_ok ) then
        errflg = 1
        errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: failure opening SOILPARM.TBL'
        return
      end if

      !write(mess,*) 'input soil texture classification = ', trim ( mminsl )
      !call wrf_message( mess )

      lumatch=0

      read (iunit_noah,*)
      read (iunit_noah,2000,end=2003)sltype
2000   format (a4)
      read (iunit_noah,*)slcats,iindex
      if(sltype.eq.mminsl)then
          !write( mess , * ) 'soil texture classification = ', trim ( sltype ) , ' found', &
          !      slcats,' categories'
          !call wrf_message ( mess )
        lumatch=1
      endif
! prevent possible array overwrite, bill bovermann, ibm, may 6, 2008
      if ( size(bb    ) < slcats .or. &
           size(drysmc) < slcats .or. &
           size(f11   ) < slcats .or. &
           size(maxsmc) < slcats .or. &
           size(refsmc) < slcats .or. &
           size(satpsi) < slcats .or. &
           size(satdk ) < slcats .or. &
           size(satdw ) < slcats .or. &
           size(wltsmc) < slcats .or. &
           size(qtz   ) < slcats  ) then
         errflg = 1
         errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: table sizes too small for value of slcats'
         return
      endif
      if(sltype.eq.mminsl)then
        do lc=1,slcats
            read (iunit_noah,*) iindex,bb(lc),drysmc(lc),f11(lc),maxsmc(lc),&
                      refsmc(lc),satpsi(lc),satdk(lc), satdw(lc),   &
                      wltsmc(lc), qtz(lc)
        enddo
      endif

2003   continue

      close (iunit_noah)
            

      ! CALL wrf_dm_bcast_integer ( LUMATCH , 1 )
      ! CALL wrf_dm_bcast_string  ( SLTYPE  , 4 )
      ! CALL wrf_dm_bcast_string  ( MMINSL  , 4 )  ! since this is reset above, see oct2 ^
      ! CALL wrf_dm_bcast_integer ( SLCATS  , 1 )
      ! CALL wrf_dm_bcast_integer ( IINDEX  , 1 )
      ! CALL wrf_dm_bcast_real    ( BB      , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( DRYSMC  , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( F11     , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( MAXSMC  , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( REFSMC  , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( SATPSI  , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( SATDK   , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( SATDW   , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( WLTSMC  , NSLTYPE )
      ! CALL wrf_dm_bcast_real    ( QTZ     , NSLTYPE )

      if(lumatch.eq.0)then
          errflg = 1
          errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: soil texture dataset '//mminsl//' not found in SOILPARM.TBL.'
          return
      endif

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
            
      open(iunit_noah, file='GENPARM.TBL',form='formatted',status='old',iostat=ierr)
      if(ierr .ne. open_ok ) then
        errflg = 1
        errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: failure opening GENPARM.TBL'
        return
      end if

      read (iunit_noah,*)
      read (iunit_noah,*)
      read (iunit_noah,*) num_slope

      slpcats=num_slope
! prevent possible array overwrite, bill bovermann, ibm, may 6, 2008
      if ( size(slope_data) < num_slope ) then
        errflg = 1
        errmsg = 'module_lsm_noah_hafs: set_soil_veg_parm: num_slope too large for slope_data array'
        return
      endif

      do lc=1,slpcats
          read (iunit_noah,*)slope_data(lc)
      enddo

      read (iunit_noah,*)
      read (iunit_noah,*)sbeta_data
      read (iunit_noah,*)
      read (iunit_noah,*)fxexp_data
      read (iunit_noah,*)
      read (iunit_noah,*)csoil_data
      read (iunit_noah,*)
      read (iunit_noah,*)salp_data
      read (iunit_noah,*)
      read (iunit_noah,*)refdk_data
      read (iunit_noah,*)
      read (iunit_noah,*)refkdt_data
      read (iunit_noah,*)
      read (iunit_noah,*)frzk_data
      read (iunit_noah,*)
      read (iunit_noah,*)zbot_data
      read (iunit_noah,*)
      read (iunit_noah,*)czil_data
      read (iunit_noah,*)
      read (iunit_noah,*)smlow_data
      read (iunit_noah,*)
      read (iunit_noah,*)smhigh_data
      read (iunit_noah,*)
      read (iunit_noah,*)lvcoef_data
      close (iunit_noah)
    

    ! call wrf_dm_bcast_integer ( num_slope    ,  1 )
    ! call wrf_dm_bcast_integer ( slpcats      ,  1 )
    ! call wrf_dm_bcast_real    ( slope_data   ,  nslope )
    ! call wrf_dm_bcast_real    ( sbeta_data   ,  1 )
    ! call wrf_dm_bcast_real    ( fxexp_data   ,  1 )
    ! call wrf_dm_bcast_real    ( csoil_data   ,  1 )
    ! call wrf_dm_bcast_real    ( salp_data    ,  1 )
    ! call wrf_dm_bcast_real    ( refdk_data   ,  1 )
    ! call wrf_dm_bcast_real    ( refkdt_data  ,  1 )
    ! call wrf_dm_bcast_real    ( frzk_data    ,  1 )
    ! call wrf_dm_bcast_real    ( zbot_data    ,  1 )
    ! call wrf_dm_bcast_real    ( czil_data    ,  1 )
    ! call wrf_dm_bcast_real    ( smlow_data   ,  1 )
    ! call wrf_dm_bcast_real    ( smhigh_data  ,  1 )
    ! call wrf_dm_bcast_real    ( lvcoef_data  ,  1 )

      end subroutine soil_veg_gen_parm
!-----------------------------
!> @}

end module lsm_noah_hafs

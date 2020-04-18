!>  \file sfc_noah_wrfv4_interstitial.F90
!!  This file contains data preparation for the WRFv4 version of Noah LSM as part of a GFS-based suite.

!> This module contains the CCPP-compliant data preparation for the WRFv4 version of Noah LSM.
      module sfc_noah_wrfv4_pre

      implicit none

      public :: sfc_noah_wrfv4_pre_init, sfc_noah_wrfv4_pre_run, sfc_noah_wrfv4_pre_finalize
      
      private

      logical :: is_initialized = .false.
      
      contains

!> \ingroup NOAH_LSM_WRFv4
!! \section arg_table_sfc_noah_wrfv4_pre_init Argument Table
!! \htmlinclude sfc_noah_wrfv4_pre_init.html
!!
      subroutine sfc_noah_wrfv4_pre_init(lsm, lsm_noah_wrfv4, veg_data_choice, &
          soil_data_choice, isurban, isice, iswater, errmsg, errflg)

      use machine, only : kind_phys
      
      implicit none
      
      integer,          intent(in)  :: lsm, lsm_noah_wrfv4, &
                                       veg_data_choice, soil_data_choice
      
      integer,          intent(inout) :: isurban, isice, iswater
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      ! Local variables
      
      character(len=256) :: mminlu, mminsl
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (is_initialized) return
      
      if (lsm/=lsm_noah_wrfv4) then
        write(errmsg,'(*(a))') "Logic error: namelist choice of LSM is different from NOAH WRFv4"
        errflg = 1
        return
      end if
      
      select case (veg_data_choice)
       case (0)
         mminlu = 'USGS'
         isurban = 1
         isice = 24
         iswater = 16
       case (1)
         mminlu = 'MODIFIED_IGBP_MODIS_NOAH'
         isurban = 13
         isice = 15
         iswater = 17
       case (3)
         mminlu = 'NLCD40'
         isurban = 13
         isice = 15 !or 22?
         iswater = 17 !or 21?
       case (4)
         mminlu = 'USGS-RUC'
         isurban = 1
         isice = 24
         iswater = 16
       case (5)
         mminlu = 'MODI-RUC'
         isurban = 13
         isice = 15
         iswater = 17
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
      
      is_initialized = .true.
      
      end subroutine sfc_noah_wrfv4_pre_init


!! \section arg_table_sfc_noah_wrfv4_pre_finalize Argument Table
!! \htmlinclude sfc_noah_wrfv4_pre_finalize.html
!!
      subroutine sfc_noah_wrfv4_pre_finalize(errmsg, errflg)

      implicit none

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      end subroutine sfc_noah_wrfv4_pre_finalize


!> \ingroup NOAH_LSM_WRFv4 Noah LSM from WRFv4 pre-scheme data preparation
!! \section arg_table_sfc_noah_wrfv4_pre_run Argument Table
!! \htmlinclude sfc_noah_wrfv4_pre_run.html
!!
!> \section general_noah_wrfv4_pre NOAH LSM WRFv4 pre-scheme data preparation General Algorithm
!>  @{
      subroutine sfc_noah_wrfv4_pre_run (im, nsoil, ialb, isice, land,         &
        flag_guess, flag_iter, restart, first_time_step, flag_lsm,             &
        flag_lsm_glacier, dt, rhowater, rd, rvrdm1, eps, epsm1, sfcprs, tprcp, &
        sfctmp, q1, prslki, wind, snwdph, cm, ch, weasd, tsfc, vtype, smc,     &
        stc, slc, snoalb, prcp, q2k, rho1, qs1, th1, dqsdt2, canopy, cmc,      &
        snowhk, chk, cmm, chh, weasd_save, snwdph_save, tsfc_save, canopy_save,&
        smc_save, stc_save, slc_save, ep, evap, hflx, gflux, drain, evbs, evcw,&
        trans, sbsno, snowc, snohf, sthick, errmsg, errflg)

      use machine , only : kind_phys
      use funcphys, only : fpvs
      use module_sf_noahlsm,  only: maxalb

      implicit none

      !GJF: Data preparation and output preparation from SFLX follows the GFS physics code (sfc_drv.F)
      ! rather than the WRF code (module_sf_noahdrv.F) in order to "fit in" with other GFS physics-based
      ! suites. Another version of this scheme (and the associated post) could potentially be
      ! created from the WRF version. No attempt was made to test sensitivities to either approach. 
      ! Note that the version of NOAH LSM expected here is "generic" - there are no urban, fasdas, or
      ! or University of Arizona(?) additions.
      
      integer,                             intent(in) :: im, nsoil, ialb, isice
      logical,                             intent(in) :: restart, first_time_step
      real(kind=kind_phys),                intent(in) :: dt, rhowater, rd, rvrdm1, eps, epsm1

      logical, dimension(im),              intent(in) :: flag_guess, flag_iter, land
      real(kind=kind_phys), dimension(im), intent(in) :: sfcprs, tprcp, sfctmp, q1, prslki, wind, cm, ch, snwdph
      real(kind=kind_phys), dimension(im), intent(in) :: weasd, tsfc, vtype
      real(kind=kind_phys), dimension(im,nsoil), intent(in) :: smc, stc, slc

      logical, dimension(im), intent(inout) :: flag_lsm, flag_lsm_glacier
      real(kind=kind_phys), dimension(im), intent(inout) :: snoalb, prcp, q2k, rho1, qs1, th1, dqsdt2, canopy, cmc, snowhk, chk, cmm, chh
      real(kind=kind_phys), dimension(im), intent(inout) :: weasd_save, snwdph_save, tsfc_save, canopy_save
      real(kind=kind_phys), dimension(im,nsoil), intent(inout) :: smc_save, stc_save, slc_save
      real(kind=kind_phys), dimension(im), intent(inout) :: ep, evap, hflx, gflux, drain, evbs, evcw, trans, sbsno, snowc, snohf
      real(kind=kind_phys), dimension(nsoil), intent(inout) :: sthick

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!     local Variables
      integer :: i, k
      real(kind=kind_phys) :: sneqv
      
      REAL, PARAMETER  :: A2=17.67,A3=273.15,A4=29.65,   &
                          A23M4=A2*(A3-A4)
      real(kind=kind_phys), parameter, dimension(4) :: zsoil = (/ -0.1,-0.4,-1.0,-2.0/) !what if nsoil /= 4?  
      
!> - Initialize CCPP error handling variables

      errmsg = ''
      errflg = 0
      
      !from module_sf_noahdrv.F/lsminit
      if (.not. restart .and. first_time_step .and. ialb == 0) then
        do i = 1, im      
             snoalb(i) = maxalb(int(0.5 + vtype(i)))*0.01
        end do
      end if
      
      do i=1, im
        if (land(i) .and. flag_guess(i)) then
          weasd_save(i) = weasd(i)
          snwdph_save(i) = snwdph(i)
          tsfc_save(i) = tsfc(i)
          canopy_save(i) = canopy(i)
        
          do k=1,nsoil
            smc_save(i,k) = smc(i,k)
            stc_save(i,k) = stc(i,k)
            slc_save(i,k) = slc(i,k)
          end do
        end if
      end do
      
      sthick(1) = - zsoil(1)
      do k = 2, nsoil
        sthick(k) = zsoil(k-1) - zsoil(k)
      enddo
      
      flag_lsm(:) = .false.
      flag_lsm_glacier(:) = .false.
      do i=1, im
        if (flag_iter(i) .and. land(i)) then
          if (vtype(i) == isice) then
            flag_lsm_glacier(i) = .true.
          else
            flag_lsm(i) = .true.
          end if
          !GJF: module_sf_noahdrv.F from WRF has hardcoded slopetyp = 1; why? replicate here?
          !GJF: shdfac is zeroed out for particular combinations of vegetation table source and vegetation types; replicate here?
          
          ep(i)     = 0.0
          evap (i)  = 0.0
          hflx (i)  = 0.0
          gflux(i)  = 0.0
          drain(i)  = 0.0

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0
          snowc(i)  = 0.0
          snohf(i)  = 0.0
          
          !GJF: could potentially pass in pre-calculated rates instead of calculating here
          prcp(i) = rhowater * tprcp(i) / dt
          
          !GJF: The GFS version of NOAH prepares the specific humidity in sfc_drv.f as follows:
          q2k(i)   = max(q1(i), 1.e-8)
          rho1(i)  = sfcprs(i) / (rd*sfctmp(i)*(1.0+rvrdm1*q2k(i)))
          
          qs1(i)   = fpvs( sfctmp(i) )
          qs1(i)   = max(eps*qs1(i) / (sfcprs(i)+epsm1*qs1(i)), 1.e-8)
          q2k(i)   = min(qs1(i), q2k(i))
          
          !GJF: could potentially pass in pre-calcualted potential temperature if other schemes also need it (to avoid redundant calculation)
          th1(i) = sfctmp(i) * prslki(i)
          
          !GJF: module_sf_noahdrv.F from WRF modifies dqsdt2 if the surface has snow.
          dqsdt2(i)=qs1(i)*a23m4/(sfctmp(i)-a4)**2
          
          !GJF: convert canopy moisture from kg m-2 to m
          canopy(i) = max(canopy(i), 0.0) !check for positive values in sfc_drv.f
          cmc(i) = canopy(i)/rhowater
          
          !GJF: snow depth passed in to NOAH is conditionally modified differently in GFS and WRF:
          sneqv = weasd(i) * 0.001
          snowhk(i) = snwdph(i) * 0.001
          if ( (sneqv /= 0.0 .and. snowhk(i) == 0.) .or. (snowhk(i) <= sneqv) ) then
            snowhk(i) = 5.*sneqv
          end if
          !GJF: GFS version:
          ! if (sneqv(i) /= 0.0 .and. snwdph(i) == 0.0) then
          !   snowhk(i) = 10.0 * sneqv(i)
          ! endif
          
          !GJF: calculate conductance from surface exchange coefficient
          chk(i) = ch(i)  * wind(i)
          
          chh(i) = chk(i) * rho1(i)
          cmm(i) = cm(i) * wind(i)
          

!GJF: If the perturbations of vegetation fraction is desired, one could uncomment this code
! and add appropriate arguments to make this work. This is from the GFS version of NOAH LSM
! in sfc_drv.f.
          
!>  - Call surface_perturbation::ppfbet() to perturb vegetation fraction that goes into gsflx().
!  perturb vegetation fraction that goes into sflx, use the same
!  perturbation strategy as for albedo (percentile matching)
!! Following Gehne et al. (2018) \cite gehne_et_al_2018, a perturbation of vegetation
!! fraction is added to account for the uncertainty. A percentile matching technique
!! is applied to guarantee the perturbed vegetation fraction is bounded between 0 and
!! 1. The standard deviation of the perturbations is 0.25 for vegetation fraction of
!! 0.5 and the perturbations go to zero as vegetation fraction  approaches its upper
!! or lower bound.
          ! vegfp  = vegfpert(i)                    ! sfc-perts, mgehne
          ! if (pertvegf(1)>0.0) then
          !         ! compute beta distribution parameters for vegetation fraction
          !         mv = shdfac
          !         sv = pertvegf(1)*mv*(1.-mv)
          !         alphav = mv*mv*(1.0-mv)/(sv*sv)-mv
          !         betav  = alphav*(1.0-mv)/mv
          !         ! compute beta distribution value corresponding
          !         ! to the given percentile albPpert to use as new albedo
          !         call ppfbet(vegfp,alphav,betav,iflag,vegftmp)
          !         shdfac = vegftmp
          ! endif
! *** sfc-perts, mgehne
        endif
      end do
      
      
      end subroutine sfc_noah_wrfv4_pre_run
      
      subroutine soil_veg_gen_parm( mminlu, mminsl, errmsg, errflg)
        !this routine is mostly taken from module_sf_noahdrv.F in WRF
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
        errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: '//   &
                 'can not find unused fortran unit to read.'
        return
      endif
      
      open(iunit_noah, file='VEGPARM.TBL',form='formatted',status='old',iostat=ierr)
      if(ierr .ne. open_ok ) then
        errflg = 1
        errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: failure opening VEGPARM.TBL'
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
                  errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: too many loops in VEGPARM.TBL'
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
         errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: table sizes too small for value of lucats'
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
           errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: expected low and high density residential, and high density industrial information in VEGPARM.TBL'
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
         errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: land use dataset '//mminlu//' not found in VEGPARM.TBL.'
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
        errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: failure opening SOILPARM.TBL'
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
         errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: table sizes too small for value of slcats'
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
          errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: soil texture dataset '//mminsl//' not found in SOILPARM.TBL.'
          return
      endif

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
            
      open(iunit_noah, file='GENPARM.TBL',form='formatted',status='old',iostat=ierr)
      if(ierr .ne. open_ok ) then
        errflg = 1
        errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: failure opening GENPARM.TBL'
        return
      end if

      read (iunit_noah,*)
      read (iunit_noah,*)
      read (iunit_noah,*) num_slope

      slpcats=num_slope
! prevent possible array overwrite, bill bovermann, ibm, may 6, 2008
      if ( size(slope_data) < num_slope ) then
        errflg = 1
        errmsg = 'sfc_noah_wrfv4_interstitial: set_soil_veg_parm: num_slope too large for slope_data array'
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

      end module sfc_noah_wrfv4_pre

      module sfc_noah_wrfv4_post
        
      implicit none

      private

      public :: sfc_noah_wrfv4_post_init, sfc_noah_wrfv4_post_run, sfc_noah_wrfv4_post_finalize

      contains
        
      subroutine sfc_noah_wrfv4_post_init ()
      end subroutine sfc_noah_wrfv4_post_init
      
      subroutine sfc_noah_wrfv4_post_finalize ()
      end subroutine sfc_noah_wrfv4_post_finalize
      
!! \section arg_table_sfc_noah_wrfv4_post_run Argument Table
!! \htmlinclude sfc_noah_wrfv4_post_run.html
!!
      subroutine sfc_noah_wrfv4_post_run (im, nsoil, land, flag_guess, flag_lsm, &
        rhowater, cp, hvap, cmc, rho1, sheat, eta, flx1, flx2, flx3, sncovr, runoff1,&
        runoff2, soilm, snowhk, weasd_save, snwdph_save, tsfc_save, tsurf,     &
        canopy_save, smc_save, stc_save, slc_save, smcmax, canopy, shflx,      &
        lhflx, snohf, snowc, runoff, drain, stm, weasd, snwdph, tsfc, smc, stc,& 
        slc, wet1, errmsg, errflg)
      
      use machine, only : kind_phys
      
      implicit none
      
      integer, intent(in) :: im, nsoil
      logical, dimension(im), intent(in) :: land, flag_guess, flag_lsm
      real(kind=kind_phys),   intent(in) :: rhowater, cp, hvap
      real(kind=kind_phys), dimension(im), intent(in) :: cmc, rho1, sheat, eta, &
        flx1, flx2, flx3, sncovr, runoff1, runoff2, soilm, snowhk
      real(kind=kind_phys), dimension(im), intent(in) :: weasd_save, snwdph_save, tsfc_save, tsurf, canopy_save, smcmax
      real(kind=kind_phys), dimension(im,nsoil), intent(in) :: smc_save, stc_save, slc_save
      
      real(kind=kind_phys), dimension(im), intent(inout) :: canopy, shflx, lhflx, &
        snohf, snowc, runoff, drain, stm, wet1
      real(kind=kind_phys), dimension(im), intent(inout) :: weasd, snwdph, tsfc
      real(kind=kind_phys), dimension(im, nsoil), intent(inout) :: smc, stc, slc
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      !local variables
      integer :: i, k
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      do i=1, im
        if (flag_lsm(i)) then
          canopy(i) = cmc(i)*rhowater
          snwdph(i) = 1000.0*snowhk(i)
          
          shflx(i) = sheat(i) / (cp*rho1(i))
          lhflx(i) = eta(i) / (hvap*rho1(i))
          
          !aggregating several outputs into one like GFS sfc_drv.F
          snohf(i) = flx1(i) + flx2(i) + flx3(i)
           
          snowc(i) = sncovr(i) !GJF: redundant?
       
          !convert from m s-1 to kg m-2 s-1 by multiplying by rhowater
          runoff(i) = runoff1(i) * rhowater
          drain(i) = runoff2(i) * rhowater
          
          stm(i) = soilm(i) * rhowater
      
          wet1(i) = smc(i,1) / smcmax(i) !Sarah Lu added 09/09/2010 (for GOCART)
        end if
      end do
      
      do i=1, im
        if (land(i)) then 
          if (flag_guess(i)) then
            weasd(i) = weasd_save(i)
            snwdph(i) = snwdph_save(i)
            tsfc(i) = tsfc_save(i)
            canopy(i) = canopy_save(i)
            
            do k=1,nsoil
              smc(i,k) = smc_save(i,k)
              stc(i,k) = stc_save(i,k)
              slc(i,k) = slc_save(i,k)
            end do
            
          else
            tsfc(i) = tsurf(i)
          end if
        end if
      end do

      end subroutine sfc_noah_wrfv4_post_run
  
      end module sfc_noah_wrfv4_post

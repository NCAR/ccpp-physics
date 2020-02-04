!>\file HWRF_rrtmg_driver.F

      MODULE HWRF_rrtmg_driver

      USE module_model_constants
      USE module_ra_rrtmg_lw   , ONLY : rrtmg_lwrad
      USE module_ra_rrtmg_sw   , ONLY : rrtmg_swrad

      IMPLICIT NONE

!MZ:  HWRF namelist
!      integer, parameter, private  :: aer_opt        = 1
!      integer, parameter, private  :: o3input        = 2
!      integer, parameter, private  :: swint_opt      = 0
!      integer, parameter, private  :: ra_call_offset = -1
!      integer, parameter, private  :: ICLOUD         = 3
!      integer, parameter, private  :: cldovrlp       = 4
      real,    parameter, private  :: cam_abs_freq_s = 21600. !default CAM clearsky longwave absorption calculation frequency
!      INTEGER, parameter, private  :: aer_type       = 1      ! aerosol type: 1 is SF79 rural, 2 is SF79 urban 
      LOGICAL, parameter, private  :: is_CAMMGMP_used  = .false. ! CAM
      LOGICAL, parameter, private  :: warm_rain  = .false.
!      INTEGER, parameter, private  :: sf_surface_physics = 2
!      INTEGER, parameter, private  :: calc_clean_atm_diag = 0

      CONTAINS

      SUBROUTINE radiation_driver (ALBEDO                               &
               ,DT                                                      &
              ,DZ8W   ,EMISS  ,GLW     ,GMT    ,GSW                     & 
              ,ITIMESTEP,JULDAY, JULIAN,JULYR                           &
              ,NPHS,  O3RAD                                             &
              ,o3input, aer_opt, swint_opt, ra_call_offset              &
              ,icloud, cldovrlp                                         &
              ,sf_surface_physics                                       &
              ,P8W  ,P ,PI  ,        RADT                               & 
              ,RHO    ,RLWTOA  ,RTHRATEN                                &
              ,RTHRATENLW    ,RTHRATENSW   ,HRSWPD, HRLWPD              &
              ,SNOW   ,STEPRA ,SWDOWN  ,SWDOWNC                         &
              ,T8W     ,T ,                  TSK                        &
              ,XICE ,XLAND  ,XLAT ,XLONG ,YR                            &
              ,sinlat, coslat, solhr                                    &
              ,COSZEN,SOLCON                                            &
!mz              ,Z                                                        &
              ,ALEVSIZ, no_src_types                                    &
              ,LEVSIZ, N_OZMIXM,   N_AEROSOLC                           &
              ,PAERLEV                                                  & 
              ,XTIME                                                    &
!mz              ,CURR_SECS                                                &
              ,its,ite,jts,jte                                          &
              ,IDS,IDE, JDS,JDE, KDS,KDE                                &
              ,IMS,IME, JMS,JME, KMS,KME                                &
              ,kts, kte                                                 &
              , CLDFRA                                                  &
#if (EM_CORE == 1)
              , lradius,iradius                                         &
#endif
              , re_cloud, re_ice, re_snow                               & ! G. Thompson
              , has_reqc, has_reqi, has_reqs                            & ! G. Thompson
!mz              , PB                                                      &
              , F_ICE_PHY,F_RAIN_PHY                                    &
              , QV                                                      &!, F_QV  &
              , QC                                                      &!, F_QC&
              , QR                                                      &!, F_QR&
              , QI                                                      &!, F_QI&
              , QS                                                      &!, F_QS&
              , QG                                                      &!, F_QG&
!              , QNDROP, F_QNDROP                                        &
              ,ACSWUPT   ,ACSWUPTC            &
              ,ACSWDNT   ,ACSWDNTC            &
              ,ACSWUPB   ,ACSWUPBC            &
              ,ACSWDNB   ,ACSWDNBC            &
              ,ACLWUPT   ,ACLWUPTC            &
              ,ACLWDNT   ,ACLWDNTC            &
              ,ACLWUPB   ,ACLWUPBC            &
              ,ACLWDNB   ,ACLWDNBC            &
              ,SWUPT ,SWUPTC, SWUPTCLN        &
              ,SWDNT ,SWDNTC, SWDNTCLN        &
              ,SWUPB ,SWUPBC, SWUPBCLN        &
              ,SWDNB ,SWDNBC, SWDNBCLN        &
              ,LWUPT ,LWUPTC, LWUPTCLN        &
              ,LWDNT ,LWDNTC, LWDNTCLN        &
              ,LWUPB ,LWUPBC, LWUPBCLN        &
              ,LWDNB ,LWDNBC, LWDNBCLN        &
              ,LWCF                           &
              ,SWCF                           &
              ,dx,dy                          &
!              , PINA, aodtot           &
              ,OZMIXM, PIN                    &
              ,CALC_CLEAN_ATM_DIAG            &
#if ( WRF_CHEM == 1 )
              ,AER_RA_FEEDBACK                &
              ,TAUAER300, TAUAER400 & ! jcb
              ,TAUAER600, TAUAER999 & ! jcb
              ,GAER300, GAER400, GAER600, GAER999 & ! jcb
              ,WAER300, WAER400, WAER600, WAER999 & ! jcb
              ,TAUAERlw1,  TAUAERlw2  &
              ,TAUAERlw3,  TAUAERlw4  &
              ,TAUAERlw5,  TAUAERlw6  &
              ,TAUAERlw7,  TAUAERlw8  &
              ,TAUAERlw9,  TAUAERlw10   &
              ,TAUAERlw11, TAUAERlw12   &
              ,TAUAERlw13, TAUAERlw14   &
              ,TAUAERlw15, TAUAERlw16  &
              ,progn                                            &
#endif
              ,SWUPFLX,SWUPFLXC,SWDNFLX,SWDNFLXC                          & ! Optional
              ,LWUPFLX,LWUPFLXC,LWDNFLX,LWDNFLXC                          & ! Optional
              ,ALSWVISDIR, ALSWVISDIF, ALSWNIRDIR, ALSWNIRDIF             & !fds ssib alb comp (06/2010)
              ,SWVISDIR, SWVISDIF, SWNIRDIR, SWNIRDIF                     & !fds ssib swr comp (06/2010)
              ,swddir,swddni,swddif                                       & ! jararias 2013/08
              ,swddirc,swddnic                                            & ! PAJ: clear-sky direct irradiance
!              ,swdown_ref,swddir_ref,coszen_ref                           & !,Gx,gg,Bx,bb               &
              ,mp_physics, imp_physics_fer_hires                          &
!              ,EFCG,EFCS,EFIG,EFIS,EFSG,EFSS                                   &!mzaercu_opt   
              ,mpirank, mpiroot,errflg, errmsg ) 
!-------------------------------------------------------------------------
   !  This driver calls subroutines for the radiation parameterizations.
   !
   !  short wave radiation choices:
   !  1. swrad (19??)
   !  4. rrtmg_sw - Added November 2008, MJIacono, AER, Inc.
   !
   !  long wave radiation choices:
   !  1. rrtmlwrad
   !  4. rrtmg_lw - Added November 2008, MJIacono, AER, Inc.
   !
!----------------------------------------------------------------------
      !mz: use GFS utilities
      use module_radiation_astronomy,only: calc_coszmn
        IMPLICIT NONE
!======================================================================
! Grid structure in physics part of WRF
! 
!-------------------------------------
! The horizontal velocities used in the physics are unstaggered 
! relative to temperature/moisture variables. All predicted 
! variables are carried at half levels except w, which is at full 
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!==================================================================
! Definitions
!-----------
! Theta      potential temperature (K)
! Qv         water vapor mixing ratio (kg/kg)
! Qc         cloud water mixing ratio (kg/kg)
! Qr         rain water mixing ratio (kg/kg)
! Qi         cloud ice mixing ratio (kg/kg)
! Qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- RTHRATEN      Theta tendency 
!                 due to radiation (K/s)
!-- RTHRATENLW    Theta tendency 
!                 due to long wave radiation (K/s)
!-- RTHRATENSW    Theta temperature tendency 
!                 due to short wave radiation (K/s)
!-- dt            time step (s)
!-- itimestep     number of time steps
!-- GLW           downward long wave flux at ground surface (W/m^2)
!-- GSW           net short wave flux at ground surface (W/m^2)
!-- SWDOWN        downward short wave flux at ground surface (W/m^2)
!-- SWDOWNC       clear-sky downward short wave flux at ground surface (W/m^2; optional; for AQ)
!-- RLWTOA        upward long wave at top of atmosphere (w/m2)
!-- RSWTOA        upward short wave at top of atmosphere (w/m2)
!-- XLAT          latitude, south is negative (degree)
!-- XLONG         longitude, west is negative (degree)
!-- ALBEDO        albedo (between 0 and 1)
!-- CLDFRA        cloud fraction (between 0 and 1)
!-- CLDFRA_DP     cloud fraction from deep cloud in a cumulus scheme
!-- CLDFRA_SH     cloud fraction from shallow cloud in a cumulus scheme
!-- CLDFRA_MP_ALL cloud fraction from CAMMGMP microphysics scheme
!-- EMISS         surface emissivity (between 0 and 1)
!-- rho_phy       density (kg/m^3)
!-- rr            dry air density (kg/m^3)
!-- moist         moisture array (4D - last index is species) (kg/kg)
!-- n_moist       number of moisture species
!-- qndrop        Cloud droplet number (#/kg)
!-- p8w           pressure at full levels (Pa)
!-- p_phy         pressure (Pa)
!-- Pb            base-state pressure (Pa)
!-- pi_phy        exner function (dimensionless)
!-- dz8w          dz between full levels (m)
!-- t_phy         temperature (K)
!-- t8w           temperature at full levels (K)
!-- GMT           Greenwich Mean Time Hour of model start (hour)
!-- JULDAY        the initial day (Julian day)
!-- RADT          time for calling radiation (min)
!-- ra_call_offset -1 (old) means usually just before output, 0 after
!-- DEGRAD        conversion factor for 
!                 degrees to radians (pi/180.) (rad/deg)
!-- DPD           degrees per day for earth's 
!                 orbital position (deg/day)
!-- R_d           gas constant for dry air (J/kg/K)
!-- CP            heat capacity at constant pressure for dry air (J/kg/K)
!-- G             acceleration due to gravity (m/s^2)
!-- rvovrd        R_v divided by R_d (dimensionless)
!-- XTIME         time since simulation start (min)
!-- DECLIN        solar declination angle (rad)
!-- SOLCON        solar constant (W/m^2)
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!
!==================================================================
!


   INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde,       &
                                       ims,ime, jms,jme, kms,kme,       &
                                       its,ite, jts,jte, kts,kte  
   INTEGER,      INTENT(IN   )    ::   mpirank, mpiroot
   character(len=*),          intent(  out) :: errmsg
   integer,                   intent(  out) :: errflg

   INTEGER,      INTENT(IN   )    ::   mp_physics, imp_physics_fer_hires 
   INTEGER,      INTENT(IN   )    ::   STEPRA     
   INTEGER,      INTENT(IN   )    ::   AER_OPT
   INTEGER,      INTENT(IN   )    ::   O3INPUT
   INTEGER,      INTENT(IN   )    ::   swint_opt
   INTEGER,      INTENT(IN   )    ::   ra_call_offset
   INTEGER,      INTENT(IN   )    ::   icloud
   INTEGER,      INTENT(IN   )    ::   cldovrlp
!   REAL,         INTENT(IN   )    ::   cam_abs_freq_s
   INTEGER,      INTENT(IN   )    ::   sf_surface_physics
   INTEGER,      INTENT(IN   )    ::   alevsiz, no_src_types
   INTEGER,      INTENT(IN   )    ::   levsiz, n_ozmixm
   INTEGER,      INTENT(IN   )    ::   paerlev, n_aerosolc      
   REAL,         INTENT(IN   )       ::   RADT
   REAL, DIMENSION( ims:ime, jms:jme ),                                 &
         INTENT(IN   )  ::                                 XLAND,       &
                                                            XICE,       &
                                                             TSK,       &
                                                          !VEGFRA, &
                                                            SNOW 
   REAL,  DIMENSION( ims:ime, levsiz, jms:jme, n_ozmixm ),  OPTIONAL,    &
          INTENT(IN   ) ::                                  OZMIXM
!   REAL,  DIMENSION( ims:ime, jms:jme ), OPTIONAL,                &
!          INTENT(INOUT) ::                                  AODTOT

   REAL,  DIMENSION( ims:ime, levsiz, jms:jme, n_ozmixm ) :: OZFLG

   REAL,  DIMENSION(levsiz), OPTIONAL, INTENT(IN )  ::     PIN
!   REAL,  DIMENSION(alevsiz), OPTIONAL, INTENT(IN )  ::     PINA

   INTEGER, INTENT(IN   )  ::   julyr
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN ) ::                                     dz8w, &
!                                                               z, &
                                                             p8w, &
                                                               p, &
                                                              pi, &
                                                               t, &
                                                             t8w, &
                                                             rho

!  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL,        &
!        INTENT(IN)     ::                            EFCG,       & 
!                                                     EFCS,       &
!                                                     EFIG,       &
!                                                     EFIS,       &
!                                                     EFSG,       &
!                                                     EFSS

#if ( WRF_CHEM == 1 )
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL ,       &
         INTENT(IN ) ::  tauaer300,tauaer400,tauaer600,tauaer999, & ! jcb
                                 gaer300,gaer400,gaer600,gaer999, & ! jcb
                                 waer300,waer400,waer600,waer999


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL ,       &
         INTENT(IN ) ::  tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4, & ! czhao 
                         tauaerlw5,tauaerlw6,tauaerlw7,tauaerlw8, & ! czhao 
                         tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12, & ! czhao 
                         tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16

   INTEGER, OPTIONAL, INTENT(IN   )    :: progn
    INTEGER, INTENT(IN   ), OPTIONAL  ::   aer_ra_feedback
#endif
   INTEGER, INTENT(IN   )  ::   calc_clean_atm_diag
   

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT)  ::                              RTHRATEN, &
                                                      RTHRATENLW, &
                                                      RTHRATENSW


   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  & 
         INTENT(OUT)  ::                                  HRLWPD, &
                                                          HRSWPD

   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                      ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC,          &
                      ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,          &
                      ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC,          &
                      ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC

! TOA and surface, upward and downward, total, clear (no cloud), and clean (no aerosol) fluxes
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
              SWUPT,  SWUPTC, SWUPTCLN,  SWDNT,  SWDNTC, SWDNTCLN,&
              SWUPB,  SWUPBC, SWUPBCLN,  SWDNB,  SWDNBC, SWDNBCLN,&
              LWUPT,  LWUPTC, LWUPTCLN,  LWDNT,  LWDNTC, LWDNTCLN,&
              LWUPB,  LWUPBC, LWUPBCLN,  LWDNB,  LWDNBC, LWDNBCLN


! Upward and downward, total and clear sky layer fluxes (W m-2)
   REAL, DIMENSION( ims:ime, kms:kme+2, jms:jme ),                &
         OPTIONAL, INTENT(INOUT) ::                               &
                               SWUPFLX,SWUPFLXC,SWDNFLX,SWDNFLXC, &
                               LWUPFLX,LWUPFLXC,LWDNFLX,LWDNFLXC

   REAL, DIMENSION( ims:ime, jms:jme ),          OPTIONAL ,       &
         INTENT(INOUT)  ::                                  SWCF, &
                                                            LWCF
   REAL, DIMENSION( ims:ime, jms:jme ),          OPTIONAL ,       &
         INTENT(IN   )  ::                            ALSWVISDIR, &
                                                      ALSWVISDIF, &
                                                      ALSWNIRDIR, &
                                                      ALSWNIRDIF
   REAL, DIMENSION( ims:ime, jms:jme ),          OPTIONAL ,       &
         INTENT(OUT  )  ::                              SWVISDIR, &
                                                        SWVISDIF, &
                                                        SWNIRDIR, &
                                                        SWNIRDIF
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(IN   )  ::                                  XLAT, &
                                                           XLONG, &
                                                          ALBEDO, &
                                                           EMISS
!
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(INOUT)  ::                                   GSW, &
                                                             GLW

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)  ::   SWDOWN

! ------------------------------------------------------------------------------ jararias 2013/08/10 -----------
   REAL, DIMENSION( ims:ime, jms:jme ),  INTENT(OUT) :: swddir, & ! All-sky SW broadband surface direct irradiance
                                                        swddni, & ! All-sky SW broadband surface direct normal irradiance
                                                        swddif    ! All-sky SW broadband surface diffuse irradiance
!   REAL, DIMENSION( ims:ime, jms:jme ),  INTENT(INOUT) ::              & !Gx,Bx,gg,bb, & ! For SW sza-interpolation
!                                                          swdown_ref,  &
!                                                          swddir_ref,  &
!                                                          coszen_ref
! ------------------------------------------------------------------------------ jararias 2013/11    -----------
!
   REAL, INTENT(IN  )   ::                                GMT,dt, &
                                                   julian, xtime
   INTEGER, INTENT(IN  ),OPTIONAL ::                          YR
!
   INTEGER, INTENT(IN  ) ::                    JULDAY, itimestep
!   REAL, INTENT(IN ),OPTIONAL     ::                    CURR_SECS
!mz: dyn_em specific
!   LOGICAL, INTENT(IN ),OPTIONAL  ::              ADAPT_STEP_FLAG

   INTEGER,INTENT(IN)                                       :: NPHS
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(INOUT)  ::                   RLWTOA 
!mz                                                      RSWTOA,     
!mz                                                      ACFRST,     & 
!mz                                                      ACFRCV       

!mz   INTEGER,DIMENSION( ims:ime, jms:jme ),INTENT(INOUT)        ::  &
!mz                                                          NCFRST, & 
!mz                                                          NCFRCV   


!
! Optional 
!
!MZ:IKJ
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         OPTIONAL,                                                &
         INTENT(INOUT) ::                                 CLDFRA   


!MZ:IKJ
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT):: re_cloud, re_ice, re_snow
   INTEGER, INTENT(INOUT):: has_reqc, has_reqi, has_reqs

!MZ:IKJ
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                     &
         OPTIONAL,                                                   &
         INTENT(IN   ) ::                                            &
                                                          F_ICE_PHY, &
                                                         F_RAIN_PHY
                                                !      CLDFRA_MP_ALL     !mz EM  specific

#if (EM_CORE == 1)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                     &
         OPTIONAL,                                                   &
         INTENT(IN   ) ::                                            &
                                                            LRADIUS, &
                                                            IRADIUS
#endif

   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         OPTIONAL,                                                &
         INTENT(OUT) ::                 SWDOWNC, SWDDIRC, SWDDNIC
!MZ:IKJ
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         OPTIONAL,                                                &
         INTENT(INOUT ) ::                                        &
!mz                                                            pb &
                                        qv,qc,qr,qi,qs,qg    !,qndrop

!mz     LOGICAL, OPTIONAL ::     f_qv,f_qc,f_qr,f_qi,f_qs,f_qg !,f_qndrop
     LOGICAL, PARAMETER :: f_qv = .true.
     LOGICAL, PARAMETER :: f_qc = .true.
     LOGICAL, PARAMETER :: f_qr = .true.
     LOGICAL, PARAMETER :: f_qi = .false.
     LOGICAL, PARAMETER :: f_qs = .true.
     LOGICAL, PARAMETER :: f_qg = .true.



     REAL, OPTIONAL, INTENT(IN) :: dx,dy
     REAL, INTENT(IN) :: SOLCON
!mz     REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(IN)  :: ht

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL ,       &
         INTENT(INOUT)  ::                       o3rad


! LOCAL  VAR
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) ::    CEMISS
   REAL, DIMENSION( ims:ime, jms:jme ) ::             coszr
   REAL, DIMENSION( ims:ime, levsiz, jms:jme )  ::    ozmixt
   REAL, DIMENSION( ims:ime, jms:jme )  :: xlatd, xlond

   REAL    ::    DECLIN   !,SOLCON 
   INTEGER ::    i,j,k,ij
!   INTEGER ::    STEPABS
   LOGICAL ::     use_aer3d !,doabsems,
   INTEGER ::    s

   REAL    ::    OBECL,SINOB,SXLONG,ARG,DECDEG,                  &
                 DJUL,RJUL,ECCFAC
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) :: qi_temp,qc_temp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) :: qi_save,qc_save
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) :: qs_save


   REAL    ::    gridkm

   REAL    ::    next_rad_time, DTaccum
!mz   LOGICAL ::    run_param , doing_adapt_dt , decided
   REAL    ::    cldji,cldlji
!------------------------------------------------------------------
! solar related variables are added to declaration
!mz*: use radiation_astronomy
!-------------------------------------------------
!mz   REAL, OPTIONAL, INTENT(OUT) :: DECLINX,SOLCONX
   REAL,   INTENT(IN) ::  SOLHR
   REAL,  DIMENSION( ims:ime, jms:jme), INTENT(IN) ::  SINLAT, COSLAT
   REAL,  DIMENSION( ims:ime, jms:jme), INTENT(IN) :: COSZEN
!local
!   REAL,  DIMENSION( ims:ime, jms:jme) :: coszdg
!------------------------------------------------------------------
   REAL :: DEGRAD=3.1415926/180.     !mz
   real :: ioh,kt,airmass,kd
   real, dimension(ims:ime,jms:jme) :: coszen_loc,hrang_loc

!  the following three arrays may be dimensioned by (ims,ime,kms,kme,jms,jme,aerosol_vars)
   real, dimension(:,:,:,:), pointer :: tauaer_sw=>null(), ssaaer_sw=>null(), asyaer_sw=>null()


         gridkm = SQRT(dx*0.001*dx*0.001 + dy*0.001*dy*0.001)

      if (itimestep .LE. 100 .and. mpirank == mpiroot) then
         WRITE ( 0 , * ) 'Grid spacing in km ', dx, dy, gridkm
      endif


!MZ
!   if(swint_opt.eq.1) then
!         CALL radconst(XTIME,DECLIN,SOLCON,JULIAN,               &
!                       DEGRAD,DPD                                )
!         call calc_coszen(ims,ime,jms,jme,its,ite,jts,jte, &
!                          julian,xtime,gmt,declin,degrad,  &
!                          xlong,xlat,coszen_loc,hrang_loc)
!   end if


! CAM-specific additional radiation frequency - cam_abs_freq_s (=21600s by default)
!     STEPABS = nint(cam_abs_freq_s/(dt*STEPRA))*STEPRA
!     IF (itimestep .eq. 1 .or. mod(itimestep,STEPABS) .eq. 1 + ra_call_offset &
!                                        .or. STEPABS .eq. 1 ) THEN
!       doabsems = .true.
!     ELSE
!       doabsems = .false.
!     ENDIF


   use_aer3d=.false.

   ! Stubs to ensure these variables are not passed in unallocated:
   if(.not.associated(tauaer_sw)) allocate(tauaer_sw(1, 1, 1, 1))
   if(.not.associated(ssaaer_sw)) allocate(ssaaer_sw(1, 1, 1, 1))
   if(.not.associated(asyaer_sw)) allocate(asyaer_sw(1, 1, 1, 1))

!---------------
!> - Calculate constant for short wave radiation
! moved up and out of OMP loop because it only needs to be computed once
! and because it is not entirely thread-safe (XT24, TOLOCTM and XXLAT need
! their thread-privacy)  JM 20100217
!       CALL radconst(XTIME,DECLIN,SOLCON,JULIAN,               &
!     &              DEGRAD,DPD                                )

!mz
!     IF(PRESENT(declinx).AND.PRESENT(solconx))THEN
! saved to state arrays used in surface driver
!       declinx=declin
!       solconx=solcon
!     ENDIF

!   outputs: coszen, hrang
!mz HWRF call
!     call calc_coszen(ims,ime,jms,jme,its,ite,jts,jte,  &
!                      julian,xtime+radt*0.5,gmt, &
!                      declin,degrad,xlong,xlat,coszen,hrang)



!     if(mpirank == mpiroot) then
!       write(0,*)'coszmn: max/min(xlong) = ',maxval(xlong),minval(xlong)
!       write(0,*)'coszmn: max/min(sinlat) = ',maxval(sinlat),minval(sinlat)
!       write(0,*)'coszmn: max/min(coslat) = ',maxval(coslat),minval(coslat)
!       write(0,*)'coszmn: solhr = ', solhr
!     endif

!       call calc_coszmn( ims,ime,jms,jme,its,ite,jts,jte,               &
!     &                   xlong,sinlat,coslat,solhr,                     &     ! ---  inputs
!     &                   coszen, coszdg     )                           ! ---  outputs

!     if(mpirank == mpiroot) then
!       write(0,*)'coszmn: max/min(coszen) = ',maxval(coszen),minval(coszen)
!       write(0,*)'coszmn: max/min(coszdg) = ',maxval(coszdg),minval(coszdg)
!     endif



! initialize data

!mz     if ((itimestep.eq.1).and.(swint_opt.eq.1)) then
!mz        do j=jts,jte
!mz           do i=its,ite
!mz              Bx(i,j)=0.
!mz              bb(i,j)=0.
!mz              Gx(i,j)=0.
!mz              gg(i,j)=0.
!           end do
!        end do
!     end if

!mz move to above
     DO j=jts,jte
     DO i=its,ite
        GSW(I,J)=0.
        GLW(I,J)=0.
        SWDOWN(I,J)=0.
        swddir(i,j)=0.  
        swddni(i,j)=0.  
        swddif(i,j)=0. 
        XLATD(I,J)=XLAT(I,J)/DEGRAD        !radians-> degree
        XLOND(I,J)=XLONG(I,J)/DEGRAD
     ENDDO
     ENDDO

     DO j=jts,jte
     DO k=kts,kte+1
     DO i=its,ite
        RTHRATEN(I,K,J)=0.
        RTHRATENLW(I,K,J)=0.
        RTHRATENSW(I,K,J)=0.
        HRLWPD(I,K,J)=0. 
        HRSWPD(I,K,J)=0. 
        CEMISS(I,K,J)=0.0
     ENDDO
     ENDDO
     ENDDO

     IF ( PRESENT( SWUPFLX ) ) THEN
        DO j=jts,jte
        DO k=kts,kte+2
        DO i=its,ite
           SWUPFLX(I,K,J) = 0.0
           SWDNFLX(I,K,J) = 0.0
           SWUPFLXC(I,K,J) = 0.0
           SWDNFLXC(I,K,J) = 0.0
           LWUPFLX(I,K,J) = 0.0
           LWDNFLX(I,K,J) = 0.0
           LWUPFLXC(I,K,J) = 0.0
           LWDNFLXC(I,K,J) = 0.0
        ENDDO
        ENDDO
        ENDDO
     ENDIF

! Fill temporary water variable
        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
           qc_temp(I,K,J)=qc(I,K,J)
        ENDDO
        ENDDO
        ENDDO

! Choose how to compute cloud fraction 
     DO j=jts,jte
     DO k=kts,kte
     DO i=its,ite
        CLDFRA(i,k,j) = 0.
     END DO
     END DO
     END DO


           DO j = jts,jte
           DO k = kts,kte
           DO i = its,ite
              qc_save(i,k,j) = qc(i,k,j)
              qi_save(i,k,j) = qi(i,k,j)
           ENDDO
           ENDDO
           ENDDO
           
!mz           IF (PRESENT(F_QS)) THEN
              DO j = jts,jte
              DO k = kts,kte
              DO i = its,ite
                 qs_save(i,k,j) = qs(i,k,j)
              ENDDO
             ENDDO
              ENDDO
!mz           ELSE
!              DO j = jts,jte
!              DO k = kts,kte
!              DO i = its,ite
!                 qs_save(i,k,j) = 0.0
!              ENDDO
!              ENDDO
!              ENDDO
!           ENDIF

           
           if(mpirank == mpiroot) then
              write(0,*)"in HWRF rad driver: call &
     & cldfra3 to use gthompson cloud fraction scheme"
             ! write(0,*)'cldfra3: max/min(p)     = ', maxval(p), minval(p)
             ! write(0,*)'cldfra3: max/min(t)     = ', maxval(t), minval(t)
             ! write(0,*)'cldfra3: max/min(qv)    = ', maxval(qv), minval(qv)
             ! write(0,*)'cldfra3: max/min(qc)    = ', maxval(qc), minval(qc)
             ! write(0,*)'cldfra3: max/min(qi)    = ', maxval(qi), minval(qi)
             ! write(0,*)'cldfra3: max/min(qs)    = ', maxval(qs), minval(qs)
             ! write(0,*)'cldfra3: max/min(rho)   = ', maxval(rho), minval(rho)
             ! write(0,*)'cldfra3: max/min(xland) = ', maxval(xland), minval(xland)
           endif
           CALL cal_cldfra3(CLDFRA, qv, qc, qi, qs,                     &
     &                 p,t,rho, XLAND, gridkm,                          &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
     &                 its,ite, jts,jte, kts,kte)
           if(mpirank == mpiroot) then
            write(0,*)'cal_cldfra3::max/min(cldfra) =', maxval(cldfra), &
     &                 minval(cldfra)
           endif

!     Interpolating climatological ozone and aerosol to model time and levels
!     Adapted from camrad code
     IF ( o3input .EQ. 2 ) THEN
! this is the CAM version interpolate to model time-step
        if(mpirank == mpiroot) then
             write(0,*)'Have o3rad'
             write(0,*)'julday, julyr, julian= ',julday, julyr, julian
             write(0,*)'max/min(ozmixm) = ', maxval(ozmixm), minval(ozmixm)
        endif
! interpolate to model time-step
        call ozn_time_int(julday,julian,ozmixm,ozmixt,levsiz,n_ozmixm,    &
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )
! interpolate to model levels
        call ozn_p_int(p ,pin, levsiz, ozmixt, o3rad, &
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )
        if(mpirank == mpiroot) then
             write(0,*)'max/min(o3rad) = ', maxval(o3rad), minval(o3rad)
        endif

     ENDIF


             if(mpirank == mpiroot) then 
               write(0,*)"call rrtmg_lw"
               write(0,*)"max/min(p8w) = ", maxval(p8w),minval(p8w)
               write(0,*)"max/min(p) = ", maxval(p),minval(p)
                write(0,*)"max/min(pi) = ", maxval(pi),minval(pi)
                write(0,*)"max/min(dz8w) = ", maxval(dz8w),minval(dz8w)
             endif
             CALL RRTMG_LWRAD(                                      &
                  RTHRATENLW=RTHRATEN,                              &
                  HRLWPD=HRLWPD,                                    & ! J. Henderson AER
                  LWUPT=LWUPT,LWUPTC=LWUPTC,LWUPTCLN=LWUPTCLN,      &
                  LWDNT=LWDNT,LWDNTC=LWDNTC,LWDNTCLN=LWDNTCLN,      &
                  LWUPB=LWUPB,LWUPBC=LWUPBC,LWUPBCLN=LWUPBCLN,      &
                  LWDNB=LWDNB,LWDNBC=LWDNBC,LWDNBCLN=LWDNBCLN,      &
                  GLW=GLW,OLR=RLWTOA,LWCF=LWCF,                     &
                  EMISS=EMISS,                                      &
                  P8W=p8w,P3D=p,PI3D=pi,DZ8W=dz8w,TSK=tsk,T3D=t,    &
                  T8W=t8w,RHO3D=rho,R=R_d,G=G,                      &
                  ICLOUD=icloud,WARM_RAIN=warm_rain,CLDFRA3D=CLDFRA,&
                  cldovrlp=cldovrlp,                                & ! J. Henderson AER: cldovrlp namelist value
#if (EM_CORE == 1)
                  LRADIUS=lradius, IRADIUS=iradius,                 &
#endif
                  IS_CAMMGMP_USED=is_cammgmp_used,                  &
                  F_ICE_PHY=F_ICE_PHY,F_RAIN_PHY=F_RAIN_PHY,        &
                  XLAND=XLAND,XICE=XICE,SNOW=SNOW,                  &
                  QV3D=QV,QC3D=QC,QR3D=QR,                          &
                  QI3D=QI,QS3D=QS,QG3D=QG,                          &
                  O3INPUT=O3INPUT,O33D=O3RAD,                       &
                  F_QV=F_QV,F_QC=F_QC,F_QR=F_QR,                    &
                  F_QI=F_QI,F_QS=F_QS,F_QG=F_QG,                    &
                  RE_CLOUD=re_cloud,RE_ICE=re_ice,RE_SNOW=re_snow,  & ! G. Thompson
                  has_reqc=has_reqc,has_reqi=has_reqi,has_reqs=has_reqs, & ! G. Thompson
#if ( WRF_CHEM == 1 )
                  TAUAERLW1=tauaerlw1,TAUAERLW2=tauaerlw2,          & ! jcb
                  TAUAERLW3=tauaerlw3,TAUAERLW4=tauaerlw4,          & ! jcb
                  TAUAERLW5=tauaerlw5,TAUAERLW6=tauaerlw6,          & ! jcb
                  TAUAERLW7=tauaerlw7,TAUAERLW8=tauaerlw8,          & ! jcb
                  TAUAERLW9=tauaerlw9,TAUAERLW10=tauaerlw10,        & ! jcb
                  TAUAERLW11=tauaerlw11,TAUAERLW12=tauaerlw12,      & ! jcb
                  TAUAERLW13=tauaerlw13,TAUAERLW14=tauaerlw14,      & ! jcb
                  TAUAERLW15=tauaerlw15,TAUAERLW16=tauaerlw16,      & ! jcb
                  aer_ra_feedback=aer_ra_feedback,                  &
                  progn=progn,                                      &
#endif
                  calc_clean_atm_diag=calc_clean_atm_diag,          &
!mz                  QNDROP3D=qndrop,F_QNDROP=f_qndrop,                &
!ccc Added for time-varying trace gases.
                  YR=YR,JULIAN=JULIAN,                              &
!ccc
                  mp_physics=mp_physics,                            &
                  imp_physics_fer_hires= imp_physics_fer_hires,     & 
                  IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde,&
                  IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme,&
                  ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte,&
                  LWUPFLX=LWUPFLX,LWUPFLXC=LWUPFLXC,                &
                  LWDNFLX=LWDNFLX,LWDNFLXC=LWDNFLXC,                &
                  mpirank=mpirank,mpiroot=mpiroot,                  &  
                  errmsg=errmsg, errflg=errflg)

        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
           RTHRATENLW(I,K,J)=RTHRATEN(I,K,J)
! OLR ALSO WILL CONTAIN OUTGOING LONGWAVE FOR RRTM (NMM HAS NO OLR ARRAY)
!mz           IF(PRESENT(OLR) .AND. K .EQ. 1)OLR(I,J)=RLWTOA(I,J)
        ENDDO
        ENDDO
        ENDDO
           
!mz

             if(mpirank==mpiroot) then
         write(0,*)'mz: max/min(RTHRATENLW) = ',                       &
     &    maxval(RTHRATENLW),minval(RTHRATENLW)


                   write(0,*) 'CALL HWRF_rrtmg_sw'
                   write(0,*)'mz*rrtmg_swrad: solcon = ',solcon
                   
             endif
                 
             CALL RRTMG_SWRAD(                                         &
                     RTHRATENSW=RTHRATENSW,                            &
                     HRSWPD=HRSWPD,                                    & ! J. Henderson AER
                     SWUPT=SWUPT,SWUPTC=SWUPTC,SWUPTCLN=SWUPTCLN,      &
                     SWDNT=SWDNT,SWDNTC=SWDNTC,SWDNTCLN=SWDNTCLN,      &
                     SWUPB=SWUPB,SWUPBC=SWUPBC,SWUPBCLN=SWUPBCLN,      &
                     SWDNB=SWDNB,SWDNBC=SWDNBC,SWDNBCLN=SWDNBCLN,      &
                     SWCF=SWCF,GSW=GSW,                                &
                     XLAT=XLATD,XLONG=XLOND,                           &
                     RADT=RADT,DEGRAD=DEGRAD,                          &
                     COSZR=COSZR,JULDAY=JULDAY,SOLCON=SOLCON,          &
                     ALBEDO=ALBEDO,t3d=t,t8w=t8w,TSK=TSK,              &
                     p3d=p,p8w=p8w,pi3d=pi,rho3d=rho,                  &
                     dz8w=dz8w,CLDFRA3D=CLDFRA,                        &
#if (EM_CORE == 1)
                     LRADIUS=lradius, IRADIUS=iradius,                 &
#endif
                     IS_CAMMGMP_USED=is_cammgmp_used,                  &
                     R=R_D,G=G,              &
                     ICLOUD=icloud,WARM_RAIN=warm_rain,                &
                     cldovrlp=cldovrlp,                                & ! J. Henderson AER: cldovrlp namelist value
                     F_ICE_PHY=F_ICE_PHY,F_RAIN_PHY=F_RAIN_PHY,        &
                     XLAND=XLAND,XICE=XICE,SNOW=SNOW,                  &
                     QV3D=qv,QC3D=qc,QR3D=qr,                          &
                     QI3D=qi,QS3D=qs,QG3D=qg,                          &
                     O3INPUT=O3INPUT,O33D=O3RAD,                       &
                     AER_OPT=AER_OPT, no_src=no_src_types,             & !mz  aerod=aerod
                     ALSWVISDIR=alswvisdir ,ALSWVISDIF=alswvisdif,     &  !Zhenxin ssib alb comp (06/2010)
                     ALSWNIRDIR=alswnirdir ,ALSWNIRDIF=alswnirdif,     &  !Zhenxin ssib alb comp (06/2010)
                     SWVISDIR=swvisdir ,SWVISDIF=swvisdif,             &  !Zhenxin ssib swr comp (06/2010)
                     SWNIRDIR=swnirdir ,SWNIRDIF=swnirdif,             &  !Zhenxin ssib swr comp (06/2010)
                     SF_SURFACE_PHYSICS=sf_surface_physics,            &  !Zhenxin ssib sw_phy   (06/2010)
                     F_QV=f_qv,F_QC=f_qc,F_QR=f_qr,                    &
                     F_QI=f_qi,F_QS=f_qs,F_QG=f_qg,                    &
                     RE_CLOUD=re_cloud,RE_ICE=re_ice,RE_SNOW=re_snow,  & ! G. Thompson
                     has_reqc=has_reqc,has_reqi=has_reqi,has_reqs=has_reqs, & ! G. Thompson
#if ( WRF_CHEM == 1 )
                     TAUAER300=tauaer300,TAUAER400=tauaer400,          & ! jcb
                     TAUAER600=tauaer600,TAUAER999=tauaer999,          & ! jcb
                     GAER300=gaer300,GAER400=gaer400,                  & ! jcb
                     GAER600=gaer600,GAER999=gaer999,                  & ! jcb
                     WAER300=waer300,WAER400=waer400,                  & ! jcb
                     WAER600=waer600,WAER999=waer999,                  & ! jcb
                     aer_ra_feedback=aer_ra_feedback,                  &
!jdfcz               progn=progn,prescribe=prescribe,                  &
                     progn=progn,                                      &
#endif
                     calc_clean_atm_diag=calc_clean_atm_diag,          &
!                     QNDROP3D=qndrop,F_QNDROP=f_qndrop,                &
                     mp_physics=mp_physics,                            &
                     imp_physics_fer_hires=imp_physics_fer_hires,      &
                     IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde,&
                     IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme,&
                     ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte,&
                     SWUPFLX=SWUPFLX,SWUPFLXC=SWUPFLXC,                &
                     SWDNFLX=SWDNFLX,SWDNFLXC=SWDNFLXC,                &
                     tauaer3d_sw=tauaer_sw,                             & ! jararias 2013/11
                     ssaaer3d_sw=ssaaer_sw,                             & ! jararias 2013/11
                     asyaer3d_sw=asyaer_sw,                             & ! jararias 2013/11
                     use_aer3d=use_aer3d,                               & ! bug fix, SGT 2015/03
                     swddir=swddir,swddni=swddni,swddif=swddif,         & ! jararias 2013/08/10
                     swdownc=swdownc, swddnic=swddnic, swddirc=swddirc, & ! PAJ
                     xcoszen=coszen,julian=julian,                      &
                     errflg=errflg, errmsg=errmsg               ) ! jararias 2013/08/14

             DO j=jts,jte
             DO k=kts,kte
             DO i=its,ite
                RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+RTHRATENSW(I,K,J)
             ENDDO
             ENDDO
             ENDDO


        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
           RTHRATENSW(I,K,J)=RTHRATEN(I,K,J)-RTHRATENLW(I,K,J)
        ENDDO
        ENDDO
        ENDDO

        DO j=jts,jte
        DO i=its,ite
           SWDOWN(I,J)=GSW(I,J)/(1.-ALBEDO(I,J))
        ENDDO
        ENDDO

             if(mpirank==mpiroot) then
            write(0,*)'mz: max/min(RTHRATENSW) = ',                     &
     &    maxval(RTHRATENSW),minval(RTHRATENSW)
          write(0,*)'mz: max/min(SWDOWN) = ',                           &
     &    maxval(SWDOWN),minval(SWDOWN)
             endif




!        IF (PRESENT(CLDFRA) .AND.                                       &
!            PRESENT(F_QC) .AND. PRESENT ( F_QI ) ) THEN
!           CALL wrf_debug (150, 'DEBUG-GT, back to micro-only Qc and Qi')
!           write(0,*)'DEBUG-HWRF-RRTMG, back to micro-only Qc and Qi'
           DO j = jts,jte
           DO k = kts,kte
           DO i = its,ite
              qc(i,k,j) = qc_save(i,k,j)
              qi(i,k,j) = qi_save(i,k,j)
           ENDDO
           ENDDO
           ENDDO
!mz           IF (PRESENT(F_QS)) THEN
              DO j = jts,jte
              DO k = kts,kte
              DO i = its,ite
                 qs(i,k,j) = qs_save(i,k,j)
              ENDDO
              ENDDO
              ENDDO
!mz           ENDIF

!mz        IF (PRESENT(CLDFRA) ) THEN
!           CALL wrf_debug (150, 'DEBUG-GT, back to micro-only Qc and Qi')
!           IF (PRESENT(F_QC)) THEN
              DO j = jts,jte
              DO k = kts,kte
              DO i = its,ite
                 qc(i,k,j) = qc_save(i,k,j)
              ENDDO
              ENDDO
              ENDDO
!           ENDIF
!           IF (PRESENT(F_QI)) THEN
              DO j = jts,jte
              DO k = kts,kte
              DO i = its,ite
                 qi(i,k,j) = qi_save(i,k,j)
              ENDDO
              ENDDO
              ENDDO
!           ENDIF
!           IF (PRESENT(F_QS)) THEN
              DO j = jts,jte
              DO k = kts,kte
              DO i = its,ite
                 qs(i,k,j) = qs_save(i,k,j)
              ENDDO
              ENDDO
              ENDDO
!           ENDIF
!mz        ENDIF

      ! parameters update for SW surface fluxes interpolation
!mz      IF (swint_opt.EQ.1) THEN
!mz         ! interpolation applies on all-sky fluxes (swddir, swdown)
!mz         CALL update_swinterp_parameters(ims,ime,jms,jme,its,ite,jts,jte,   &
!                                         coszen,coszen_loc,swddir,swdown,   &
!                                         swddir_ref,bb,Bx,swdown_ref,gg,Gx, &
!                                         coszen_ref                         )
!      ENDIF

!mz   ENDDO
!mz   !$OMP END PARALLEL DO

   IF ( associated(tauaer_sw) ) deallocate(tauaer_sw)
   IF ( associated(ssaaer_sw) ) deallocate(ssaaer_sw)
   IF ( associated(asyaer_sw) ) deallocate(asyaer_sw)

!mz   ENDIF Radiation_step

 ! jararias, aug 2013
 ! SW surface fluxes interpolation (meaningful when not in a Radiation_step)
!mz
! if (swint_opt .eq. 1) then
!    call wrf_debug(100,'SW surface irradiance interpolation')
!
!      call interp_sw_radiation(ims,ime,jms,jme,its,ite,jts,jte,  &
!                               coszen_ref,coszen_loc,swddir_ref, &
!                               bb,Bx,swdown_ref,gg,Gx,albedo,    &
!                               swdown,swddir,swddni,swddif,gsw   )
! end if



   IF(PRESENT(LWUPTC))THEN
!  NMM calls the driver every RADT time steps, EM calls every DT
#if (EM_CORE == 1)
   DTaccum = DT
#else
   DTaccum = RADT*60
#endif
!mz   !$OMP PARALLEL DO   &
!mz   !$OMP PRIVATE ( ij ,i,j,k,its,ite,jts,jte)

!mz   DO ij = 1 , num_tiles
!mz     its = i_start(ij)
!mz     ite = i_end(ij)
!mz     jts = j_start(ij)
!mz     jte = j_end(ij)

        DO j=jts,jte
        DO i=its,ite
           ACLWUPT(I,J) = ACLWUPT(I,J) + LWUPT(I,J)*DTaccum
           ACLWUPTC(I,J) = ACLWUPTC(I,J) + LWUPTC(I,J)*DTaccum
           ACLWDNT(I,J) = ACLWDNT(I,J) + LWDNT(I,J)*DTaccum
           ACLWDNTC(I,J) = ACLWDNTC(I,J) + LWDNTC(I,J)*DTaccum
           ACLWUPB(I,J) = ACLWUPB(I,J) + LWUPB(I,J)*DTaccum
           ACLWUPBC(I,J) = ACLWUPBC(I,J) + LWUPBC(I,J)*DTaccum
           ACLWDNB(I,J) = ACLWDNB(I,J) + LWDNB(I,J)*DTaccum
           ACLWDNBC(I,J) = ACLWDNBC(I,J) + LWDNBC(I,J)*DTaccum
        ENDDO
        ENDDO
!mz   ENDDO
   ENDIF

   IF(PRESENT(SWUPTC))THEN
!  NMM calls the driver every RADT time steps, EM calls every DT
#if (EM_CORE == 1)
   DTaccum = DT
#else
   DTaccum = RADT*60
#endif
!mz   !$OMP PARALLEL DO   &
!mz   !$OMP PRIVATE ( ij ,i,j,k,its,ite,jts,jte)

!mz   DO ij = 1 , num_tiles
!mz     its = i_start(ij)
!mz     ite = i_end(ij)
!mz     jts = j_start(ij)
!mz     jte = j_end(ij)

        DO j=jts,jte
        DO i=its,ite
           ACSWUPT(I,J) = ACSWUPT(I,J) + SWUPT(I,J)*DTaccum
           ACSWUPTC(I,J) = ACSWUPTC(I,J) + SWUPTC(I,J)*DTaccum
           ACSWDNT(I,J) = ACSWDNT(I,J) + SWDNT(I,J)*DTaccum
           ACSWDNTC(I,J) = ACSWDNTC(I,J) + SWDNTC(I,J)*DTaccum
           ACSWUPB(I,J) = ACSWUPB(I,J) + SWUPB(I,J)*DTaccum
           ACSWUPBC(I,J) = ACSWUPBC(I,J) + SWUPBC(I,J)*DTaccum
           ACSWDNB(I,J) = ACSWDNB(I,J) + SWDNB(I,J)*DTaccum
           ACSWDNBC(I,J) = ACSWDNBC(I,J) + SWDNBC(I,J)*DTaccum
        ENDDO
        ENDDO
!mz   ENDDO
!mz   !$OMP END PARALLEL DO
   ENDIF


!mz:cldt is EM specific
! compute cloud diagnosis (random overlapping)
! IF ( PRESENT ( CLDFRA ) .AND. PRESENT ( CLDT ) .AND.        &
!      PRESENT ( F_QC ) .AND. PRESENT ( F_QI ) ) THEN

!mz   DO ij = 1 , num_tiles
!mz     its = i_start(ij)
!mz     ite = i_end(ij)
!mz     jts = j_start(ij)
!mz     jte = j_end(ij)

!        DO j=jts,jte
!        DO i=its,ite
!          cldji=1.0
!          do k=kte-1,kts,-1
!            cldji=cldji*(1.0-cldfra(i,k,j))
!          enddo
!          cldt(i,j)=1.0-cldji
!        END DO
!        END DO
!mz    END DO
! END IF

   END SUBROUTINE radiation_driver


!---------------------------------------------------------------------
   SUBROUTINE radconst(XTIME,DECLIN,SOLCON,JULIAN,                   &
                       DEGRAD,DPD                                    )
!---------------------------------------------------------------------
!mz   USE module_wrf_error
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   REAL, INTENT(IN   )      ::       DEGRAD,DPD,XTIME,JULIAN
   REAL, INTENT(OUT  )      ::       DECLIN,SOLCON
   REAL                     ::       OBECL,SINOB,SXLONG,ARG,  &
                                     DECDEG,DJUL,RJUL,ECCFAC
!
! !DESCRIPTION:
! Compute terms used in radiation physics 
!EOP

! for short wave radiation

   DECLIN=0.
   SOLCON=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.
        
   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)
        
!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:
        
   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)
   SXLONG=SXLONG*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)
   DECDEG=DECLIN/DEGRAD
!----SOLAR CONSTANT ECCENTRICITY FACTOR (PALTRIDGE AND PLATT 1976)
   DJUL=JULIAN*360./365.
   RJUL=DJUL*DEGRAD
   ECCFAC=1.000110+0.034221*COS(RJUL)+0.001280*SIN(RJUL)+0.000719*  &
          COS(2*RJUL)+0.000077*SIN(2*RJUL)
   SOLCON=1370.*ECCFAC
   
   END SUBROUTINE radconst

   SUBROUTINE calc_coszen(ims,ime,jms,jme,its,ite,jts,jte,  &
                          julian,xtime,gmt, &
                          declin,degrad,xlon,xlat,coszen,hrang)
       ! Added Equation of Time correction : jararias, 2013/08/10
       implicit none
       integer, intent(in) :: ims,ime,jms,jme,its,ite,jts,jte
       real, intent(in)    :: julian,declin,xtime,gmt,degrad
       real, dimension(ims:ime,jms:jme), intent(in)    :: xlat,xlon
       real, dimension(ims:ime,jms:jme), intent(inout) :: coszen,hrang

       integer :: i,j
       real    :: da,eot,xt24,tloctm,xxlat

       da=6.2831853071795862*(julian-1)/365.
       eot=(0.000075+0.001868*cos(da)-0.032077*sin(da) &
            -0.014615*cos(2*da)-0.04089*sin(2*da))*(229.18)
       xt24=mod(xtime,1440.)+eot
       do j=jts,jte
          do i=its,ite
             tloctm=gmt+xt24/60.+xlon(i,j)/15.
             hrang(i,j)=15.*(tloctm-12.)*degrad
             xxlat=xlat(i,j)*degrad
             coszen(i,j)=sin(xxlat)*sin(declin) &
                        +cos(xxlat)*cos(declin) *cos(hrang(i,j))
          enddo
       enddo
   END SUBROUTINE calc_coszen

   subroutine update_swinterp_parameters(ims,ime,jms,jme,its,ite,jts,jte, &
                                         coszen,coszen_loc,swddir,swdown, &
                                         swddir_ref,bb,Bx,                &
                                         swdown_ref,gg,Gx,                &
                                         coszen_ref                       )
      ! Author: jararias 2013/11
      implicit None
      integer, intent(in) :: ims,ime,jms,jme,its,ite,jts,jte
      real, dimension(ims:ime,jms:jme), intent(in)    :: coszen,coszen_loc,swddir,swdown
      real, dimension(ims:ime,jms:jme), intent(inout) :: swddir_ref,bb,Bx, &
                                                         swdown_ref,gg,Gx, &
                                                         coszen_ref

      integer :: i,j
      real :: swddir_0,swdown_0,coszen_0
      real, parameter :: coszen_min=1e-4

      do j=jts,jte
         do i=its,ite
            if ((coszen(i,j).gt.coszen_min) .and. (coszen_loc(i,j).gt.coszen_min)) then
               ! parameters update for DIR
               if (Bx(i,j).le.0) then
                  swddir_0 =(coszen_loc(i,j)/coszen(i,j))*swddir(i,j) ! linear first guess estimation
                  coszen_0 =coszen_loc(i,j)
               else
                  swddir_0 =swddir_ref(i,j)
                  coszen_0 =coszen_ref(i,j)
               end if
               if ((coszen(i,j)/coszen_0).lt.1.) then
                  bb(i,j) =log(max(1.,swddir(i,j))/max(1.,swddir_0)) / log(min(1.-1e-4,coszen(i,j)/coszen_0))
               elseif ((coszen(i,j)/coszen_0).gt.1) then
                  bb(i,j) =log(max(1.,swddir(i,j))/max(1.,swddir_0)) / log(max(1.+1e-4,coszen(i,j)/coszen_0))
               else
                  bb(i,j) =0.
               end if
               bb(i,j) =max(-.5,min(2.5,bb(i,j)))
               Bx(i,j) =swddir(i,j)/(coszen(i,j)**bb(i,j))

               !write(wrf_err_message,*) 'XXX I=',i,' J=',j,'  Bx=',Bx(i,j),'  bb=',bb(i,j),'  swddir=',swddir(i,j), &
               !                         '  swddir_0=',swddir_0,'  coszen=',coszen(i,j),'  coszen_0=',coszen_0
               !call wrf_debug(1,wrf_err_message)

               ! parameters update for GHI
               if (Gx(i,j).le.0) then
                  swdown_0 =(coszen_loc(i,j)/coszen(i,j))*swdown(i,j) ! linear first guess estimation
                  coszen_0 =coszen_loc(i,j)
               else
                  swdown_0 =swdown_ref(i,j)
                  coszen_0 =coszen_ref(i,j)
               end if
               if ((coszen(i,j)/coszen_0).lt.1.) then
                  gg(i,j) =log(max(1.,swdown(i,j))/max(1.,swdown_0)) / log(min(1.-1e-4,coszen(i,j)/coszen_0))
               elseif ((coszen(i,j)/coszen_0).gt.1) then
                  gg(i,j) =log(max(1.,swdown(i,j))/max(1.,swdown_0)) / log(max(1.+1e-4,coszen(i,j)/coszen_0))
               else
                  gg(i,j) =0.
               end if
               gg(i,j) =max(-.5,min(2.5,gg(i,j)))
               Gx(i,j) =swdown(i,j)/(coszen(i,j)**gg(i,j))
            else
               Bx(i,j) =0.
               bb(i,j) =0.
               Gx(i,j) =0.
               gg(i,j) =0.
            end if

            ! saving last SW run in state variables
            coszen_ref(i,j) =coszen(i,j)
            swdown_ref(i,j) =swdown(i,j)
            swddir_ref(i,j) =swddir(i,j)

            !if ((i.eq.20).and.(j.eq.20)) then
            !   write(wrf_err_message,'("   RADSTEP : tn=",I4," csz_0=",F9.6," csz=",F9.6," csz_1=",F9.6," Gx=",F14.2," gg=",F9.5,  &
            !                           " Bx=",F14.2," bb=",F9.5)') itimestep,coszen_0,coszen_loc(i,j),coszen(i,j),Gx(i,j),gg(i,j), &
            !                           Bx(i,j),bb(i,j)
            !   call wrf_debug(1,wrf_err_message)
            !end if

         end do
      end do

   end subroutine update_swinterp_parameters

   subroutine interp_sw_radiation(ims,ime,jms,jme,its,ite,jts,jte,  &
                                  coszen_ref,coszen_loc,swddir_ref, &
                                  bb,Bx,swdown_ref,gg,Gx,albedo,    &
                                  swdown,swddir,swddni,swddif,gsw   )
      ! Author: jararias 2013/11
      implicit None
      integer, intent(in) :: ims,ime,jms,jme,its,ite,jts,jte
      real, dimension(ims:ime,jms:jme), intent(in) :: coszen_ref,coszen_loc, &
                                                      swddir_ref,Bx,bb,      &
                                                      swdown_ref,Gx,gg,      &
                                                      albedo

      real, dimension(ims:ime,jms:jme), intent(inout) :: swddir,swdown, &
                                                         swddif,swddni, gsw

      integer :: i,j
      real, parameter :: coszen_min=1e-4

      do j=jts,jte
         do i=its,ite
            ! sza interpolation of surface fluxes
            if ((coszen_ref(i,j).gt.coszen_min) .and. (coszen_loc(i,j).gt.coszen_min)) then
               if ((bb(i,j).eq.-0.5).or.(bb(i,j).eq.2.5).or.(bb(i,j).eq.0.0)) then
                  swddir(i,j) =(coszen_loc(i,j)/coszen_ref(i,j))*swddir_ref(i,j)
               else
                  swddir(i,j) =Bx(i,j)*(coszen_loc(i,j)**bb(i,j))
               end if
               if ((gg(i,j).eq.-0.5).or.(gg(i,j).eq.2.5).or.(gg(i,j).eq.0.0)) then
                  swdown(i,j) =(coszen_loc(i,j)/coszen_ref(i,j))*swdown_ref(i,j)
               else
                  swdown(i,j) =Gx(i,j)*(coszen_loc(i,j)**gg(i,j))
               end if
               swddif(i,j) =swdown(i,j)-swddir(i,j)
               swddni(i,j) =swddir(i,j)/coszen_loc(i,j)
               gsw(i,j)    =swdown(i,j)*(1.-albedo(i,j))
            else
               swddir(i,j) =0.
               swdown(i,j) =0.
               swddif(i,j) =0.
               swddni(i,j) =0.
               gsw(i,j)    =0.
            end if
         end do
      end do
   end subroutine interp_sw_radiation


!+---+-----------------------------------------------------------------+
!..Cloud fraction scheme by G. Thompson (NCAR-RAL), not intended for
!.. combining with any cumulus or shallow cumulus parameterization
!.. scheme cloud fractions.  This is intended as a stand-alone for
!.. cloud fraction and is relatively good at getting widespread stratus
!.. and stratoCu without caring whether any deep/shallow Cu param schemes
!.. is making sub-grid-spacing clouds/precip.  Under the hood, this
!.. scheme follows Mocko and Cotton (1995) in applicaiton of the
!.. Sundqvist et al (1989) scheme but using a grid-scale dependent
!.. RH threshold, one each for land v. ocean points based on
!.. experiences with HWRF testing.
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+

      SUBROUTINE cal_cldfra3(CLDFRA, qv, qc, qi, qs,                    &
     &                 p,t,rho, XLAND, gridkm,                          &
!    &                 rand_perturb_on, kme_stoch, rand_pert,           &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
     &                 its,ite, jts,jte, kts,kte)
!
      USE module_mp_thompson   , ONLY : rsif, rslf
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: ids,ide, jds,jde, kds,kde,                  &
     &                      ims,ime, jms,jme, kms,kme,                  &
!    &                      kme_stoch,                                  &
     &                      its,ite, jts,jte, kts,kte

!     INTEGER, INTENT(IN):: rand_perturb_on
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN):: qv,p,t,rho
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT):: qc,qi,qs
!     REAL, DIMENSION(ims:ime,kms:kme_stoch,jms:jme), INTENT(IN):: rand_pert
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN):: XLAND

      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT):: cldfra
      REAL, INTENT(IN):: gridkm

!..Local vars.
      REAL:: RH_00L, RH_00O, RH_00, RHI_max, entrmnt
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme):: qvsat
      INTEGER:: i,j,k
      REAL:: TK, TC, qvsi, qvsw, RHUM, xx, yy
      REAL, DIMENSION(kts:kte):: qvs1d, cfr1d, T1d,                     &
     &                           P1d, R1d, qc1d, qi1d, qs1d

      character*512 dbg_msg
      LOGICAL:: debug_flag

!+---+

!..First cut scale-aware. Higher resolution should require closer to
!.. saturated grid box for higher cloud fraction.  Simple functions
!.. chosen based on Mocko and Cotton (1995) starting point and desire
!.. to get near 100% RH as grid spacing moves toward 1.0km, but higher
!.. RH over ocean required as compared to over land.

      RH_00L = 0.7  + SQRT(1./(25.0+gridkm*gridkm*gridkm))
      RH_00O = 0.81 + SQRT(1./(50.0+gridkm*gridkm*gridkm))

      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
         RHI_max = 0.0
         CLDFRA(I,K,J) = 0.0

         if (qc(i,k,j).gt.1.E-6 .or. qi(i,k,j).ge.1.E-7 .or. qs(i,k,j).gt.1.E-5) then
            CLDFRA(I,K,J) = 1.0
            qvsat(i,k,j) = qv(i,k,j)
         else
            TK   = t(i,k,j)
            TC   = TK - 273.16

            qvsw = rslf(P(i,k,j), TK)
            qvsi = rsif(P(i,k,j), TK)

            if (tc .ge. -12.0) then
               qvsat(i,k,j) = qvsw
            elseif (tc .lt. -20.0) then
               qvsat(i,k,j) = qvsi
            else
               qvsat(i,k,j) = qvsw - (qvsw-qvsi)*(-12.0-tc)/(-12.0+20.)
            endif
            RHUM = MAX(0.01, MIN(qv(i,k,j)/qvsat(i,k,j), 0.9999))

            IF ((XLAND(I,J)-1.5).GT.0.) THEN                             !--- Ocean
               RH_00 = RH_00O
            ELSE                                                         !--- Land
               RH_00 = RH_00L
            ENDIF

            if (tc .ge. -12.0) then
               RHUM = MIN(0.999, RHUM)
               CLDFRA(I,K,J) = MAX(0.0, 1.0-SQRT((1.0-RHUM)/(1.-RH_00)))
            elseif (tc.lt.-12..and.tc.gt.-70. .and. RHUM.gt.RH_00L) then
               RHUM = MAX(0.01, MIN(qv(i,k,j)/qvsat(i,k,j), 1.0 - 1.E-6))
               CLDFRA(I,K,J) = MAX(0., 1.0-SQRT((1.0-RHUM)/(1.0-RH_00L)))
            endif
            CLDFRA(I,K,J) = MIN(0.90, CLDFRA(I,K,J))

         endif
      ENDDO
      ENDDO
      ENDDO

!..Prepare for a 1-D column to find various cloud layers.

      DO j = jts,jte
      DO i = its,ite
!        if (i.gt.10.and.i.le.20 .and. j.gt.10.and.j.le.20) then
!          debug_flag = .true.
!        else
!           debug_flag = .false.
!        endif

!        if (rand_perturb_on .eq. 1) then
!           entrmnt = MAX(0.01, MIN(0.99, 0.5 + rand_pert(i,1,j)*0.5))
!        else
            entrmnt = 0.5
!        endif

         DO k = kts,kte
            qvs1d(k) = qvsat(i,k,j)
            cfr1d(k) = cldfra(i,k,j)
            T1d(k) = t(i,k,j)
            P1d(k) = p(i,k,j)
            R1d(k) = rho(i,k,j)
            qc1d(k) = qc(i,k,j)
            qi1d(k) = qi(i,k,j)
            qs1d(k) = qs(i,k,j)
         ENDDO

!     if (debug_flag) then
!       WRITE (dbg_msg,*) 'DEBUG-GT: finding cloud layers at point  (', i, ', ', j, ')'
!       CALL wrf_debug (150, dbg_msg)
!     endif
         call find_cloudLayers(qvs1d, cfr1d, T1d, P1d, R1d, entrmnt,    &
     &                         debug_flag, qc1d, qi1d, qs1d, kts,kte)

         DO k = kts,kte
            cldfra(i,k,j) = cfr1d(k)
            qc(i,k,j) = qc1d(k)
            qi(i,k,j) = qi1d(k)
         ENDDO
      ENDDO
      ENDDO


      END SUBROUTINE cal_cldfra3

!+---+-----------------------------------------------------------------+
!..From cloud fraction array, find clouds of multi-level depth and compute
!.. a reasonable value of LWP or IWP that might be contained in that depth,
!.. unless existing LWC/IWC is already there.

      SUBROUTINE find_cloudLayers(qvs1d, cfr1d, T1d, P1d, R1d, entrmnt, &
     &                            debugfl, qc1d, qi1d, qs1d, kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts, kte
      LOGICAL, INTENT(IN):: debugfl
      REAL, INTENT(IN):: entrmnt
      REAL, DIMENSION(kts:kte), INTENT(IN):: qvs1d,T1d,P1d,R1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: cfr1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qc1d, qi1d, qs1d

!..Local vars.
      REAL, DIMENSION(kts:kte):: theta, dz
      REAL:: Z1, Z2, theta1, theta2, ht1, ht2
      INTEGER:: k, k2, k_tropo, k_m12C, k_m40C, k_cldb, k_cldt, kbot
      LOGICAL:: in_cloud
      character*512 dbg_msg

!+---+

      k_m12C = 0
      k_m40C = 0
      DO k = kte, kts, -1
         theta(k) = T1d(k)*((100000.0/P1d(k))**(287.05/1004.))
         if (T1d(k)-273.16 .gt. -40.0 .and. P1d(k).gt.7000.0) k_m40C = MAX(k_m40C, k)
         if (T1d(k)-273.16 .gt. -12.0 .and. P1d(k).gt.10000.0) k_m12C = MAX(k_m12C, k)
      ENDDO
      if (k_m40C .le. kts) k_m40C = kts
      if (k_m12C .le. kts) k_m12C = kts

      Z2 = 44307.692 * (1.0 - (P1d(kte)/101325.)**0.190)
      DO k = kte-1, kts, -1
         Z1 = 44307.692 * (1.0 - (P1d(k)/101325.)**0.190)
         dz(k+1) = Z2 - Z1
         Z2 = Z1
      ENDDO
      dz(kts) = dz(kts+1)

!..Find tropopause height, best surrogate, because we would not really
!.. wish to put fake clouds into the stratosphere.  The 10/1500 ratio
!.. d(Theta)/d(Z) approximates a vertical line on typical SkewT chart
!.. near typical (mid-latitude) tropopause height.  Since messy data
!.. could give us a false signal of such a transition, do the check over 
!.. three K-level change, not just a level-to-level check.  This method
!.. has potential failure in arctic-like conditions with extremely low
!.. tropopause height, as would any other diagnostic, so ensure resulting
!.. k_tropo level is above 4km.

      DO k = kte-3, kts, -1
         theta1 = theta(k)
         theta2 = theta(k+2)
         ht1 = 44307.692 * (1.0 - (P1d(k)/101325.)**0.190)
         ht2 = 44307.692 * (1.0 - (P1d(k+2)/101325.)**0.190)
         if ( (((theta2-theta1)/(ht2-ht1)) .lt. 10./1500. ) .AND.       &
     &                       (ht1.lt.19000.) .and. (ht1.gt.4000.) ) then 
            goto 86
         endif
      ENDDO
 86   continue
      k_tropo = MAX(kts+2, k+2)

!     if (debugfl) then
!     print*, ' FOUND TROPOPAUSE ', k_tropo, ' near ', ht2, ' m'
!       WRITE (dbg_msg,*) 'DEBUG-GT: FOUND TROPOPAUSE ', k_tropo, ' near ', ht2, ' m'
!       CALL wrf_debug (150, dbg_msg)
!     endif

!..Eliminate possible fractional clouds above supposed tropopause.
      DO k = k_tropo+1, kte
         if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.0.999) then
            cfr1d(k) = 0.
         endif
      ENDDO

!..We would like to prevent fractional clouds below LCL in idealized
!.. situation with deep well-mixed convective PBL, that otherwise is
!.. likely to get clouds in more realistic capping inversion layer.

      kbot = kts+2
      DO k = kbot, k_m12C
         if ( (theta(k)-theta(k-1)) .gt. 0.05E-3*dz(k)) EXIT
      ENDDO
      kbot = MAX(kts+1, k-2)
      DO k = kts, kbot
         if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.0.999) cfr1d(k) = 0.
      ENDDO


!..Starting below tropo height, if cloud fraction greater than 1 percent,
!.. compute an approximate total layer depth of cloud, determine a total 
!.. liquid water/ice path (LWP/IWP), then reduce that amount with tuning 
!.. parameter to represent entrainment factor, then divide up LWP/IWP
!.. into delta-Z weighted amounts for individual levels per cloud layer. 

      k_cldb = k_tropo
      in_cloud = .false.
      k = k_tropo
      DO WHILE (.not. in_cloud .AND. k.gt.k_m12C)
         k_cldt = 0
         if (cfr1d(k).ge.0.01) then
            in_cloud = .true.
            k_cldt = MAX(k_cldt, k)
         endif
         if (in_cloud) then
            DO k2 = k_cldt-1, k_m12C, -1
               if (cfr1d(k2).lt.0.01 .or. k2.eq.k_m12C) then
                  k_cldb = k2+1
                  goto 87
               endif
            ENDDO
 87         continue
            in_cloud = .false.
         endif
         if ((k_cldt - k_cldb + 1) .ge. 2) then
!     if (debugfl) then
!           print*, 'An ice cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       WRITE (dbg_msg,*) 'DEBUG-GT: An ice cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       CALL wrf_debug (150, dbg_msg)
!     endif
            call adjust_cloudIce(cfr1d, qi1d, qs1d, qvs1d, T1d,R1d,dz,  &
     &                           entrmnt, k_cldb,k_cldt,kts,kte)
            k = k_cldb
         else
            if (cfr1d(k_cldb).gt.0.and.qi1d(k_cldb).lt.1.E-6)           &
     &               qi1d(k_cldb)=1.E-5*cfr1d(k_cldb)
         endif


         k = k - 1
      ENDDO


      k_cldb = k_tropo
      in_cloud = .false.
      k = k_m12C + 2
      DO WHILE (.not. in_cloud .AND. k.gt.kbot)
         k_cldt = 0
         if (cfr1d(k).ge.0.01) then
            in_cloud = .true.
            k_cldt = MAX(k_cldt, k)
         endif
         if (in_cloud) then
            DO k2 = k_cldt-1, kbot, -1
               if (cfr1d(k2).lt.0.01 .or. k2.eq.kbot) then
                  k_cldb = k2+1
                  goto 88
               endif
            ENDDO
 88         continue
            in_cloud = .false.
         endif
         if ((k_cldt - k_cldb + 1) .ge. 2) then
!     if (debugfl) then
!           print*, 'A water cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       WRITE (dbg_msg,*) 'DEBUG-GT: A water cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       CALL wrf_debug (150, dbg_msg)
!     endif
            call adjust_cloudH2O(cfr1d, qc1d, qvs1d, T1d,R1d,dz,        &
     &                           entrmnt, k_cldb,k_cldt,kts,kte)
            k = k_cldb
         else
            if (cfr1d(k_cldb).gt.0.and.qc1d(k_cldb).lt.1.E-6)           &
     &               qc1d(k_cldb)=1.E-5*cfr1d(k_cldb)
         endif
         k = k - 1
      ENDDO

!..Do a final total column adjustment since we may have added more than 1mm
!.. LWP/IWP for multiple cloud decks.

      call adjust_cloudFinal(cfr1d, qc1d, qi1d, R1d,dz, kts,kte,k_tropo)

!     if (debugfl) then
!     print*, ' Made-up fake profile of clouds'
!     do k = kte, kts, -1
!        write(*,'(i3, 2x, f8.2, 2x, f9.2, 2x, f6.2, 2x,  f15.7, 2x, f15.7)') &
!    &        K, T1d(k)-273.15, P1d(k)*0.01, cfr1d(k)*100., qc1d(k)*1000.,qi1d(k)*1000.
!     enddo
!       WRITE (dbg_msg,*) 'DEBUG-GT:  Made-up fake profile of clouds'
!       CALL wrf_debug (150, dbg_msg)
!       do k = kte, kts, -1
!          write(dbg_msg,'(f8.2, 2x, f9.2, 2x, f6.2, 2x,  f15.7, 2x, f15.7)') &
!    &          T1d(k)-273.15, P1d(k)*0.01, cfr1d(k)*100., qc1d(k)*1000.,qi1d(k)*1000.
!          CALL wrf_debug (150, dbg_msg)
!       enddo
!     endif


      END SUBROUTINE find_cloudLayers

!+---+-----------------------------------------------------------------+

      SUBROUTINE adjust_cloudIce(cfr,qi,qs,qvs, T,Rho,dz, entr, k1,k2,kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: k1,k2, kts,kte
      REAL, INTENT(IN):: entr
      REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, qvs, T, Rho, dz
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qi, qs
      REAL:: iwc, max_iwc, tdz, this_iwc, this_dz, iwp_exists
      INTEGER:: k, kmid

      tdz = 0.
      do k = k1, k2
         tdz = tdz + dz(k)
      enddo
      kmid = NINT(0.5*(k1+k2))
      max_iwc = ABS(qvs(k2-1)-qvs(k1))
!     print*, ' max_iwc = ', max_iwc, ' over DZ=',tdz

      iwp_exists = 0.
      do k = k1, k2
         iwp_exists = iwp_exists + (qi(k)+qs(k))*Rho(k)*dz(k)
      enddo
      if (iwp_exists .gt. 1.0) RETURN

      this_dz = 0.0
      do k = k1, k2
         if (k.eq.k1) then
            this_dz = this_dz + 0.5*dz(k)
         else
            this_dz = this_dz + dz(k)
         endif
         this_iwc = max_iwc*this_dz/tdz
         iwc = MAX(1.E-6, this_iwc*(1.-entr))
         if (cfr(k).gt.0.01.and.cfr(k).lt.0.99.and.T(k).ge.203.16) then
            qi(k) = qi(k) + 0.1*cfr(k)*iwc
         elseif (qi(k).lt.1.E-5.and.cfr(k).ge.0.99.and.T(k).ge.203.16) then
            qi(k) = qi(k) + 0.01*iwc
         endif
      enddo

      END SUBROUTINE adjust_cloudIce

!+---+-----------------------------------------------------------------+

      SUBROUTINE adjust_cloudH2O(cfr, qc, qvs, T,Rho,dz, entr, k1,k2,kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: k1,k2, kts,kte
      REAL, INTENT(IN):: entr
      REAL, DIMENSION(kts:kte):: cfr, qc, qvs, T, Rho, dz
      REAL:: lwc, max_lwc, tdz, this_lwc, this_dz, lwp_exists
      INTEGER:: k, kmid

      tdz = 0.
      do k = k1, k2
         tdz = tdz + dz(k)
      enddo
      kmid = NINT(0.5*(k1+k2))
      max_lwc = ABS(qvs(k2-1)-qvs(k1))
!     print*, ' max_lwc = ', max_lwc, ' over DZ=',tdz

      lwp_exists = 0.
      do k = k1, k2
         lwp_exists = lwp_exists + qc(k)*Rho(k)*dz(k)
      enddo
      if (lwp_exists .gt. 1.0) RETURN

      this_dz = 0.0
      do k = k1, k2
         if (k.eq.k1) then
            this_dz = this_dz + 0.5*dz(k)
         else
            this_dz = this_dz + dz(k)
         endif
         this_lwc = max_lwc*this_dz/tdz
         lwc = MAX(1.E-6, this_lwc*(1.-entr))
         if (cfr(k).gt.0.01.and.cfr(k).lt.0.99.and.T(k).lt.298.16.and.T(k).ge.253.16) then
            qc(k) = qc(k) + cfr(k)*cfr(k)*lwc
         elseif (cfr(k).ge.0.99.and.qc(k).lt.1.E-5.and.T(k).lt.298.16.and.T(k).ge.253.16) then
            qc(k) = qc(k) + 0.1*lwc
         endif
      enddo

      END SUBROUTINE adjust_cloudH2O

!+---+-----------------------------------------------------------------+

!..Do not alter any grid-explicitly resolved hydrometeors, rather only
!.. the supposed amounts due to the cloud fraction scheme.

      SUBROUTINE adjust_cloudFinal(cfr, qc, qi, Rho,dz, kts,kte,k_tropo)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts,kte,k_tropo
      REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, Rho, dz
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qc, qi
      REAL:: lwp, iwp, xfac
      INTEGER:: k

      lwp = 0.
      do k = kts, k_tropo
         if (cfr(k).gt.0.0) then
            lwp = lwp + qc(k)*Rho(k)*dz(k)
         endif
      enddo

      iwp = 0.
      do k = kts, k_tropo
         if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
            iwp = iwp + qi(k)*Rho(k)*dz(k)
         endif
      enddo

      if (lwp .gt. 1.5) then
         xfac = 1./lwp
         do k = kts, k_tropo
            if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
               qc(k) = qc(k)*xfac
            endif
         enddo
      endif

      if (iwp .gt. 1.5) then
         xfac = 1./iwp
         do k = kts, k_tropo
            if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
               qi(k) = qi(k)*xfac
            endif
         enddo
      endif

      END SUBROUTINE adjust_cloudFinal

!+---+-----------------------------------------------------------------+


SUBROUTINE ozn_time_int(julday,julian,ozmixm,ozmixt,levsiz,num_months,  &
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )

! adapted from oznint from CAM module
!  input: ozmixm - read from physics_init
! output: ozmixt - time interpolated

!  USE module_ra_cam_support, ONLY : getfactors

   IMPLICIT NONE

   INTEGER,    INTENT(IN) ::           ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   INTEGER,      INTENT(IN   )    ::   levsiz, num_months

   REAL,  DIMENSION( ims:ime, levsiz, jms:jme, num_months ),      &
          INTENT(IN   ) ::                                  ozmixm

   INTEGER, INTENT(IN )      ::        JULDAY
   REAL,    INTENT(IN )      ::        JULIAN

   REAL,  DIMENSION( ims:ime, levsiz, jms:jme ),      &
          INTENT(OUT  ) ::                                  ozmixt

   !Local
   REAL      :: intJULIAN
   integer   :: np1,np,nm,m,k,i,j
   integer   :: IJUL
   integer, dimension(12) ::  date_oz
   data date_oz/16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350/
   real, parameter :: daysperyear = 365.  ! number of days in a year
   real      :: cdayozp, cdayozm
   real      :: fact1, fact2, deltat
   logical   :: finddate
   logical   :: ozncyc
   CHARACTER(LEN=256) :: msgstr

   ozncyc = .true.
   ! JULIAN starts from 0.0 at 0Z on 1 Jan.
   intJULIAN = JULIAN + 1.0       ! offset by one day
! jan 1st 00z is julian=1.0 here
   IJUL=INT(intJULIAN)
!  Note that following will drift.
!    Need to use actual month/day info to compute julian.
   intJULIAN=intJULIAN-FLOAT(IJUL)
   IJUL=MOD(IJUL,365)
   IF(IJUL.EQ.0)IJUL=365
   intJULIAN=intJULIAN+IJUL
   np1=1
   finddate=.false.

!  do m=1,num_months
   do m=1,12
      if(date_oz(m).gt.intjulian.and..not.finddate) then
        np1=m
        finddate=.true.
      endif
   enddo
   cdayozp=date_oz(np1)

   if(np1.gt.1) then
      cdayozm=date_oz(np1-1)
      np=np1
      nm=np-1
   else
      cdayozm=date_oz(12)
      np=np1
      nm=12
   endif

!  call getfactors(ozncyc,np1, cdayozm, cdayozp,intjulian, &
!                   fact1, fact2)
!
! Determine time interpolation factors.  Account for December-January
! interpolation if dataset is being cycled yearly.
!
   if (ozncyc .and. np1 == 1) then                      ! Dec-Jan interpolation
      deltat = cdayozp + daysperyear - cdayozm
      if (intjulian > cdayozp) then                     ! We are in December
         fact1 = (cdayozp + daysperyear - intjulian)/deltat
         fact2 = (intjulian - cdayozm)/deltat
      else                                              ! We are in January
         fact1 = (cdayozp - intjulian)/deltat
         fact2 = (intjulian + daysperyear - cdayozm)/deltat
      end if
   else
      deltat = cdayozp - cdayozm
      fact1 = (cdayozp - intjulian)/deltat
      fact2 = (intjulian - cdayozm)/deltat
   end if
!
! Time interpolation.
!
      do j=jts,jte
      do k=1,levsiz
      do i=its,ite
            ozmixt(i,k,j) = ozmixm(i,k,j,nm)*fact1 + ozmixm(i,k,j,np)*fact2
      end do
      end do
      end do

END SUBROUTINE ozn_time_int

SUBROUTINE ozn_p_int(p ,pin, levsiz, ozmixt, o3vmr, &
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )

!-----------------------------------------------------------------------
!
! Purpose: Interpolate ozone from current time-interpolated values to model levels
!
! Method: Use pressure values to determine interpolation levels
!
! Author: Bruce Briegleb
! WW: Adapted for general use
!
!--------------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------------
!
! Arguments
!
   INTEGER,    INTENT(IN) ::           ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   integer, intent(in) :: levsiz              ! number of ozone layers

   real, intent(in) :: p(ims:ime,kms:kme,jms:jme)   ! level pressures (mks, bottom-up)
   real, intent(in) :: pin(levsiz)        ! ozone data level pressures (mks, top-down)
   real, intent(in) :: ozmixt(ims:ime,levsiz,jms:jme) ! ozone mixing ratio

   real, intent(out) :: o3vmr(ims:ime,kms:kme,jms:jme) ! ozone volume mixing ratio
!
! local storage
!
   real    pmid(its:ite,kts:kte)
   integer i,j                 ! longitude index
   integer k, kk, kkstart, kout! level indices
   integer kupper(its:ite)     ! Level indices for interpolation
   integer kount               ! Counter
   integer ncol, pver

   real    dpu                 ! upper level pressure difference
   real    dpl                 ! lower level pressure difference

   ncol = ite - its + 1
   pver = kte - kts + 1

   do j=jts,jte
!
! Initialize index array
!
!  do i=1, ncol
   do i=its, ite
      kupper(i) = 1
   end do
!
! Reverse the pressure array, and pin is in Pa, the same as model pmid
!
      do k = kts,kte
         kk = kte - k + kts
      do i = its,ite
         pmid(i,kk) = p(i,k,j)
      enddo
      enddo

   do k=1,pver

      kout = pver - k + 1
!     kout = k
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = levsiz
!     do i=1,ncol
      do i=its,ite
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
      do kk=kkstart,levsiz-1
!        do i=1,ncol
         do i=its,ite
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
         if (kount.eq.ncol) then
!           do i=1,ncol
            do i=its,ite
               dpu = pmid(i,k) - pin(kupper(i))
               dpl = pin(kupper(i)+1) - pmid(i,k)
               o3vmr(i,kout,j) = (ozmixt(i,kupper(i),j)*dpl + &
                             ozmixt(i,kupper(i)+1,j)*dpu)/(dpl + dpu)
            end do
            goto 35
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top ozone data level for at least some
! of the longitude points.
!
!     do i=1,ncol
      do i=its,ite
         if (pmid(i,k) .lt. pin(1)) then
            o3vmr(i,kout,j) = ozmixt(i,1,j)*pmid(i,k)/pin(1)
         else if (pmid(i,k) .gt. pin(levsiz)) then
            o3vmr(i,kout,j) = ozmixt(i,levsiz,j)
         else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            o3vmr(i,kout,j) = (ozmixt(i,kupper(i),j)*dpl + &
                          ozmixt(i,kupper(i)+1,j)*dpu)/(dpl + dpu)
         end if
      end do

      if (kount.gt.ncol) then
!        call endrun ('OZN_P_INT: Bad ozone data: non-monotonicity suspected')
!mz         call wrf_error_fatal ('OZN_P_INT: Bad ozone data: non-monotonicity suspected')
      end if
35    continue

   end do
   end do

   return
END SUBROUTINE ozn_p_int

!mz to-do work: move to radiation_gases
SUBROUTINE aer_time_int(julday,julian,aerodm,aerodt,levsiz,num_months,no_src,  &
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )

! adapted from oznint from CAM module
!  input: aerodm - read from physics_init
! output: aerodt - time interpolated

!  USE module_ra_cam_support, ONLY : getfactors

   IMPLICIT NONE

   INTEGER,    INTENT(IN) ::           ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   INTEGER,      INTENT(IN   )    ::   levsiz, num_months, no_src

   REAL,  DIMENSION( ims:ime, levsiz, jms:jme, num_months, no_src ),      &
          INTENT(IN   ) ::                                  aerodm

   INTEGER, INTENT(IN )      ::        JULDAY
   REAL,    INTENT(IN )      ::        JULIAN

   REAL,  DIMENSION( ims:ime, levsiz, jms:jme, no_src ),      &
          INTENT(OUT  ) ::                                  aerodt

   !Local
   REAL      :: intJULIAN
   integer   :: np1,np,nm,m,k,i,j,s
   integer   :: IJUL
   integer, dimension(12) ::  date_oz
   data date_oz/16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350/
   real, parameter :: daysperyear = 365.  ! number of days in a year
   real      :: cdayozp, cdayozm
   real      :: fact1, fact2, deltat
   logical   :: finddate
   logical   :: ozncyc
   CHARACTER(LEN=256) :: msgstr

   ozncyc = .true.
   ! JULIAN starts from 0.0 at 0Z on 1 Jan.
   intJULIAN = JULIAN + 1.0       ! offset by one day
! jan 1st 00z is julian=1.0 here
   IJUL=INT(intJULIAN)
!  Note that following will drift.
!    Need to use actual month/day info to compute julian.
   intJULIAN=intJULIAN-FLOAT(IJUL)
   IJUL=MOD(IJUL,365)
   IF(IJUL.EQ.0)IJUL=365
   intJULIAN=intJULIAN+IJUL
   np1=1
   finddate=.false.

!  do m=1,num_months
   do m=1,12
      if(date_oz(m).gt.intjulian.and..not.finddate) then
        np1=m
        finddate=.true.
      endif
   enddo
   cdayozp=date_oz(np1)

   if(np1.gt.1) then
      cdayozm=date_oz(np1-1)
      np=np1
      nm=np-1
   else
      cdayozm=date_oz(12)
      np=np1
      nm=12
   endif

!  call getfactors(ozncyc,np1, cdayozm, cdayozp,intjulian, &
!                   fact1, fact2)
!
! Determine time interpolation factors.  Account for December-January
! interpolation if dataset is being cycled yearly.
!
   if (ozncyc .and. np1 == 1) then                      ! Dec-Jan interpolation
      deltat = cdayozp + daysperyear - cdayozm
      if (intjulian > cdayozp) then                     ! We are in December
         fact1 = (cdayozp + daysperyear - intjulian)/deltat
         fact2 = (intjulian - cdayozm)/deltat
      else                                              ! We are in January
         fact1 = (cdayozp - intjulian)/deltat
         fact2 = (intjulian + daysperyear - cdayozm)/deltat
      end if
   else
      deltat = cdayozp - cdayozm
      fact1 = (cdayozp - intjulian)/deltat
      fact2 = (intjulian - cdayozm)/deltat
   end if
!
! Time interpolation.
!
      do s=1, no_src
      do j=jts,jte
      do k=1,levsiz
      do i=its,ite
            aerodt(i,k,j,s) = aerodm(i,k,j,nm,s)*fact1 + aerodm(i,k,j,np,s)*fact2
      end do
      end do
      end do
      end do

END SUBROUTINE aer_time_int

!mz to-do work move to radiation_gases
SUBROUTINE aer_p_int(p ,pin, levsiz, aerodt, aerod, no_src, pf, totaod,   &
                     ids , ide , jds , jde , kds , kde ,     &
                     ims , ime , jms , jme , kms , kme ,     &
                     its , ite , jts , jte , kts , kte )

!-----------------------------------------------------------------------
!
! Purpose: Interpolate aerosol from current time-interpolated values to model levels
!
! Method: Use pressure values to determine interpolation levels
!
! Author: Bruce Briegleb
! WW: Adapted for general use
!
!   p:  model level pressure at half levels (Pa, bottom-up)
!   pf: model level pressure at full levles (Pa, bottom-up)
!
!--------------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------------
!
! Arguments
!
   INTEGER,    INTENT(IN) ::           ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   integer, intent(in) :: levsiz              ! number of aerosol layers
   integer, intent(in) :: no_src              ! types of aerosol 

   real, intent(in) :: p(ims:ime,kms:kme,jms:jme)
   real, intent(in) :: pf(ims:ime,kms:kme,jms:jme)
   real, intent(in) :: pin(levsiz)        ! aerosol data level pressures (mks, top-down)
   real, intent(in) :: aerodt(ims:ime,levsiz,jms:jme,1:no_src) ! aerosol optical depth

   real, intent(out) :: aerod(ims:ime,kms:kme,jms:jme,1:no_src) ! aerosol optical depth
   real, intent(out) :: totaod(ims:ime,jms:jme)                 ! total aerosol optical depth
!
! local storage
!
   real    pmid(its:ite,kts:kte)
   integer i,j                 ! longitude index
   integer k, kk, kkstart, kout! level indices
   integer kupper(its:ite)     ! Level indices for interpolation
   integer kount               ! Counter
   integer ncol, pver, s

   real    dpu                 ! upper level pressure difference
   real    dpl                 ! lower level pressure difference
   real    dpm                 ! pressure difference in a model layer surrounding half p

   ncol = ite - its + 1
   pver = kte - kts + 1

   do s=1,no_src
   do j=jts,jte
!
! Initialize index array
!
   do i=its, ite
      kupper(i) = 1
   end do
!
! The pressure from incoming data is in hPa and top-down, 
!     while model pressure is in Pa and bottom-up
!
      do k = kts,kte
         kk = kte - k + kts
      do i = its,ite
         pmid(i,kk) = p(i,k,j)*0.01
      enddo
      enddo

   do k=1,pver

      kout = pver - k + 1
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = levsiz
      do i=its,ite
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
      do kk=kkstart,levsiz-1
         do i=its,ite
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
         if (kount.eq.ncol) then
            do i=its,ite
               dpu = pmid(i,k) - pin(kupper(i))
               dpl = pin(kupper(i)+1) - pmid(i,k)
               dpm = pf(i,kout,j) - pf(i,kout+1,j)
               aerod(i,kout,j,s) = (aerodt(i,kupper(i),j,s)*dpl + &
                             aerodt(i,kupper(i)+1,j,s)*dpu)/(dpl + dpu)
               aerod(i,kout,j,s) = aerod(i,kout,j,s)*dpm
            end do
            goto 35
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top aerosol data level for at least some
! of the longitude points.
!
      do i=its,ite
         if (pmid(i,k) .lt. pin(1)) then
            dpm = pf(i,kout,j) - pf(i,kout+1,j)
            aerod(i,kout,j,s) = aerodt(i,1,j,s)*pmid(i,k)/pin(1)
            aerod(i,kout,j,s) = aerod(i,kout,j,s)*dpm
         else if (pmid(i,k) .gt. pin(levsiz)) then
            dpm = pf(i,kout,j) - pf(i,kout+1,j)
            aerod(i,kout,j,s) = aerodt(i,levsiz,j,s)
            aerod(i,kout,j,s) = aerod(i,kout,j,s)*dpm
         else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            dpm = pf(i,kout,j) - pf(i,kout+1,j)
            aerod(i,kout,j,s) = (aerodt(i,kupper(i),j,s)*dpl + &
                          aerodt(i,kupper(i)+1,j,s)*dpu)/(dpl + dpu)
            aerod(i,kout,j,s) = aerod(i,kout,j,s)*dpm
         end if
      end do

      if (kount.gt.ncol) then
!         call wrf_error_fatal ('AER_P_INT: Bad aerosol data: non-monotonicity suspected')
      end if
35    continue

   end do
   end do
   end do

   do j=jts,jte
   do i=its,ite
      totaod(i,j) = 0.
   end do
   end do

   do s=1,no_src
   do j=jts,jte
   do k=1,pver
   do i=its,ite
      totaod(i,j) = totaod(i,j) + aerod(i,k,j,s)
   end do
   end do
   end do
   end do

   return
END SUBROUTINE aer_p_int

END MODULE HWRF_rrtmg_driver

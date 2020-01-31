!>\file module_HWRF_radiation.F
!! This is the HWRF RRTMG driver wrapper module.

      MODULE HWRF_radiation
!
      use machine, only : kind_phys
!-----------------------------------------------------------------------
      USE MODULE_MODEL_CONSTANTS
      !USE module_radiation_astronomy ,ONLY : CAL_MON_DAY,ZENITH
      USE HWRF_rrtmg_driver
!-----------------------------------------------------------------------
      implicit none
!
      public :: HWRF_radiation_init, HWRF_radiation_run,                &
     &          HWRF_radiation_finalize
      
      private
      LOGICAL :: acswalloc = .false.
      LOGICAL :: aclwalloc = .false.


      CONTAINS
     
!! \section arg_table_HWRF_radiation_init Argument Table
!! \htmlinclude HWRF_radiation_init.html
!!
      SUBROUTINE HWRF_radiation_init(ncol, nlev, KDT,                   &
     &                RTHRATEN,RTHRATENLW,                              &
     &                RTHRATENSW,CLDFRA,                                &
     &                levsiz,XLAT,XLONG,num_ozmixm,                     &
     &                alevsiz,no_src_types,                             &
     &                has_reqc,has_reqi,has_reqs,                       &
     &                ozmixm,pin,                                       &
     &                aerodm,pina,                                      &
     &                mpicomm, mpiroot,mpirank,                         &
     &                errmsg, errflg )
!---------------------------------------------------------------------
      USE module_ra_rrtmg_lw  , ONLY : rrtmg_lwinit
      USE module_ra_rrtmg_sw  , ONLY : rrtmg_swinit
      USE module_radiation_gases, ONLY : AEROSOL_IN, oznini
!---------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER,               INTENT(IN   )    :: levsiz, num_ozmixm
      INTEGER,               INTENT(IN   )    :: NCOL, nlev, KDT
      INTEGER,               INTENT(IN   )    :: alevsiz, no_src_types
      INTEGER,               INTENT(INOUT)    :: has_reqc, has_reqi,    &
     &                                           has_reqs
      INTEGER,               INTENT(IN   )    :: mpirank
      INTEGER,               INTENT(IN   )    :: mpiroot
      INTEGER,               INTENT(IN   )    :: mpicomm
      REAL(KIND_PHYS),  DIMENSION(1:ncol) , INTENT(IN) ::  XLAT, XLONG
      REAL(KIND_PHYS),  DIMENSION(1:ncol, levsiz, num_ozmixm ),         &
     &                                              OPTIONAL,           &
     &                                          INTENT(INOUT) :: OZMIXM

      REAL(KIND_PHYS),                                                  &
     &      DIMENSION(1:ncol, alevsiz, num_ozmixm-1, no_src_types ),    &
     &                                                   OPTIONAL,      &
     &     INTENT(INOUT) ::                                  aerodm
      REAL(KIND_PHYS), DIMENSION(levsiz), OPTIONAL, INTENT(INOUT) ::PIN
      REAL(KIND_PHYS), DIMENSION(alevsiz), OPTIONAL, INTENT(INOUT)::PINA
      REAL(KIND_PHYS), DIMENSION(1:ncol, 1:nlev) , INTENT(INOUT) ::     &
     &                                                        RTHRATEN, &
     &                                                      RTHRATENLW, &
     &                                                      RTHRATENSW, &
     &                                                      CLDFRA
      CHARACTER(LEN=*),          INTENT(  OUT) :: errmsg
      INTEGER,                   INTENT(  OUT) :: errflg

      !local
      LOGICAL, PARAMETER             :: allowed_to_read =.true.
      INTEGER                        :: ids, ide, jds, jde, kds, kde,   &
     &                                  ims, ime, jms, jme, kms, kme,   &
     &                                  its, ite, jts, jte, kts, kte

      INTEGER :: i, j, k
      REAL(KIND_PHYS) :: DEGRAD=3.1415926/180. 
      !mz in degree
      REAL(KIND_PHYS),  DIMENSION(1:ncol) :: xlatd, xlond


      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (aclwalloc .and. acswalloc) return

      if (mpirank==mpiroot) then
        write(0,*) '---------------------------------------------------'
        write(0,*) '-                   WARNING                       -'
        write(0,*) '-  the CCPP HWRF RRTMG LW/SW scheme is currently  -'
        write(0,*) '-  under development, use at your own risk        -'
        write(0,*) '-                   WARNING                       -'
        write(0,*) '---------------------------------------------------'
      end if

      !Determin if we will compute and pass radiative effective radii of 
      !cloud water, ice, and snow. Currently not supported in
      !Ferrier-Aligo with RRTMG-LW/SW.
      has_reqc = 0
      has_reqi = 0
      has_reqs = 0


      ! Set internal dimensions
      ids = 1
      ims = 1
      its = 1
      ide = ncol
      ime = ncol
      ite = ncol
      jds = 1
      jms = 1
      jts = 1
      jde = 1
      jme = 1
      jte = 1
      kds = 1
      kms = 1
      kts = 1
      kde = nlev
      kme = nlev
      kte = nlev


!---------------------------------------------------------------------

!-- calculate radiation time step

!      STEPRA = nint(RADT*60./DT)
!      STEPRA = max(STEPRA,1)

!-- initialization
!MZ: rad_reset 
      IF(KDT .LE. 1) THEN
      DO k=1,nlev
      DO i=1,ncol
        RTHRATEN(i,k)=0.
        RTHRATENLW(i,k)=0.
        RTHRATENSW(i,k)=0.
        CLDFRA(i,k)=0.
      ENDDO
      ENDDO
      ENDIF
!mz convert xlat/xlong to degree
     
      DO i=1,ncol       
        xlatd(i)= xlat(i)/degrad
        xlond(i) = xlong(i)/degrad
      enddo

!-- use CAM ozone and some aerosol profiles in HWRF rad schemes
!   n_ozmixm: no of months; levsiz: = 59, vertical dim
!   Read in CAM ozone data, and interpolate data to model grid
!   Interpolation of aerosols is done on domain 1 only

      CALL oznini(ozmixm,pin,levsiz,num_ozmixm,XLATD,                   &
     &               ids, ide, jds, jde, kds, kde,                      &
     &               ims, ime, jms, jme, kms, kme,                      &
     &               its, ite, jts, jte, kts, kte,                      &
     &               mpicomm, mpirank, mpiroot, errflg, errmsg)
                   
      if (errflg /= 0) return
     

      CALL aerosol_in(aerodm,pina,alevsiz,num_ozmixm-1,                 &
     &               no_src_types,XLATD,XLOND,                          &
     &               ids, ide, jds, jde, kds, kde,                      &
     &               ims, ime, jms, jme, kms, kme,                      &
     &               its, ite, jts, jte, kts, kte)
      if (errflg /= 0) return

      CALL rrtmg_lwinit(                                                &
     &               allowed_to_read,                                   &
     &               ids, ide, jds, jde, kds, kde,                      &
     &               ims, ime, jms, jme, kms, kme,                      &
     &               its, ite, jts, jte, kts, kte,                      &
     &               mpirank, mpiroot, mpicomm,                         &
     &               errmsg, errflg       )
     
      if (errflg /= 0) return
      if(mpirank ==mpiroot)write (0,*)'HWRF rrtmg_lwinit completed...'

      aclwalloc = .true.

      CALL rrtmg_swinit(                                                &
     &               allowed_to_read ,                                  &
     &               ids, ide, jds, jde, kds, kde,                      &
     &               ims, ime, jms, jme, kms, kme,                      &
     &               its, ite, jts, jte, kts, kte,                      &
     &               mpirank, mpiroot, mpicomm,                         &
     &               errmsg, errflg )
      if (errflg /= 0) return
      if (mpirank == mpiroot)write(0,*)'HWRF rrtmg_swinit completed' 

      acswalloc = .true.

      END SUBROUTINE HWRF_radiation_init

      SUBROUTINE HWRF_radiation_finalize()
      END SUBROUTINE HWRF_radiation_finalize
!

!> \section arg_table_HWRF_radiation_run Argument Table
!! \htmlinclude HWRF_radiation_run.html
!!
       SUBROUTINE HWRF_radiation_run (NCOL,NLEV, NTSD,DT                &
     &                    ,JULDAY,JULYR,JULIAN                          &
     &                    ,IHRST,NPHS                                   &
     &                    ,NRADS,NRADL, DX ,p8w,prsl,tsfc,T,Q           &
     &                    ,QV,QC,QI,QR,QS,QG                            & !MOIST: dry mixing ratio
!     &                    ,EPSR                                        & !Radtend%semis
     &                    ,F_ICE,F_RAIN,slmsk,CLDFRA                       &
     &                    ,RLWTT,RSWTT,RLWIN,RSWIN,RSWINC,RSWOUT,GSW    &
     &                    ,RLWTOA,CZMEAN                                &
     &                    ,SNOW,SICE, NUM_OZMIXM,OZMIXM,PIN,LEVSIZ      &
     &                    ,RTHRATEN,RTHRATENLW,RTHRATENSW               &
     &                    ,re_cloud,re_ice,re_snow                      &
     &                    ,has_reqc,has_reqi,has_reqs                   &
     &                    ,swddir, swddni,swddif, swddirc,swddnic       &
     &                    ,sinlat,coslat,solhr,solcon, radtend          &
     &                    ,xlat,xlong                                   &
     &                    ,paerlev, ALEVSIZ,no_src_types                &
     &                    ,icloud,cldovrlp                              &
     &                    ,RA_CALL_OFFSET, o3input,aer_opt,o3rad        &
     &                    ,sf_surface_physics                           &
     &                    ,swint_opt                                    &
     &                    ,calc_clean_atm_diag                          &
     &                    ,mp_physics, imp_physics_fer_hires            &
     &                    ,hrswpd,hrlwpd                                &
     &                    ,SWUPT,SWUPTC,SWUPTCLN,SWDNT,SWDNTC, SWDNTCLN &
     &                    ,SWUPB,SWUPBC,SWUPBCLN,SWDNB,SWDNBC,SWDNBCLN  &
     &                    ,LWUPT,LWUPTC,LWUPTCLN,LWDNT,LWDNTC,LWDNTCLN  &
     &                    ,LWUPB,LWUPBC,LWUPBCLN,LWDNB,LWDNBC,LWDNBCLN  &
     &                    ,ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC            &
     &                    ,ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC            &
     &                    ,ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC            &
     &                    ,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC            &
!mz-new
!     &                    ,SWUPFLX,SWUPFLXC,SWDNFLX,SWDNFLXC            &
!     &                    ,LWUPFLX,LWUPFLXC,LWDNFLX,LWDNFLXC            &
!     &                    ,ALSWVISDIR, ALSWVISDIF                       &
!     &                    ,ALSWNIRDIR, ALSWNIRDIF                       &
!mz-new
     &                    ,SWVISDIR ,SWVISDIF                           &
     &                    ,SWNIRDIR, SWNIRDIF                           &
     &                    ,MPIROOT, MPIRANK,MPICOMM                     &
     &                    ,ERRMSG, ERRFLG    )

!***********************************************************************
      USE module_radiation_astronomy ,ONLY : CAL_MON_DAY,ZENITH
      USE GFS_typedefs,               ONLY : GFS_radtend_type

      IMPLICIT NONE
     
      type(GFS_radtend_type),              intent(in) :: Radtend

      INTEGER,         INTENT(IN) :: IHRST,JULDAY,JULYR                 &
     &                              ,NPHS,NRADL,NRADS,NTSD              &
     &                              ,NUM_OZMIXM

      INTEGER,         INTENT(IN) :: NCOL,NLEV 
      REAL(KIND_PHYS), INTENT(IN) :: DX(1:ncol)
      integer,         intent(in)    :: mp_physics
      integer,         intent(in)    :: imp_physics_fer_hires

      REAL(KIND_PHYS), INTENT(IN) :: OZMIXM(1:NCOL,LEVSIZ,NUM_OZMIXM)   &
     &                               ,PIN(LEVSIZ)
      integer,         intent(in) :: swint_opt, calc_clean_atm_diag
      INTEGER,         INTENT(IN) :: o3input, aer_opt
      INTEGER,         INTENT(IN) :: ICLOUD,ra_call_offset
      INTEGER,         INTENT(IN) :: cldovrlp        
      INTEGER,         INTENT(IN) :: alevsiz, no_src_types
      INTEGER,         INTENT(IN) :: levsiz !, n_ozmixm
      INTEGER,         INTENT(IN) :: paerlev!, n_aerosolc
!      REAL(KIND_PHYS), INTENT(IN) :: cam_abs_freq_s
!      LOGICAL,         INTENT(IN) :: warm_rain
      INTEGER,    OPTIONAL, INTENT(IN   )    ::  sf_surface_physics

      REAL(KIND_PHYS), INTENT(IN) :: DT,JULIAN
      REAL(KIND_PHYS), DIMENSION(1:ncol,1:nlev), OPTIONAL ,             &
     &                 INTENT(INOUT)  ::                       o3rad


      REAL(KIND_PHYS), DIMENSION(1:ncol),      INTENT(IN)    ::         & !ALBEDO  &
!     &                                              EPSR                &
     &                                             SICE,slmsk          &
     &                                             ,SNOW
!MZ      REAL(KIND_PHYS),DIMENSION(1:ncol),       INTENT(INOUT) :: CUPPT
      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlev),INTENT(IN)    :: Q,T,    &
     &                                                           prsl
      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlev+1),INTENT(IN)  :: p8w 
      REAL(KIND_PHYS),DIMENSION(1:ncol),       INTENT(IN)  :: tsfc

      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlev),INTENT(IN)    :: F_ICE   &
     &                                                         ,F_RAIN
      REAL(KIND_PHYS),DIMENSION(1:NCOL),       INTENT(INOUT) :: SWDDIR  &
     &                                                         ,SWDDNI  &
     &                                                         ,SWDDIF  &
     &                                                         ,SWDDNIC &
     &                                                         ,SWDDIRC 
      REAL(KIND_PHYS), DIMENSION(1:NCOL), INTENT(in) :: sinlat, coslat
      REAL(KIND_PHYS), INTENT(in) :: solhr,solcon
!      REAL(KIND_PHYS), DIMENSION(1:NCOL), INTENT(OUT) :: COSZEN
!mz      REAL(KIND_PHYS), OPTIONAL, DIMENSION(1:NCOL), INTENT(OUT) :: HRANG
      REAL(KIND_PHYS),DIMENSION(1:ncol) , INTENT(IN) ::  XLAT, XLONG
      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlev),INTENT(INOUT) ::RTHRATEN &
     &                                                      ,RTHRATENSW &
     &                                                       ,RTHRATENLW
      REAL(KIND_PHYS),DIMENSION(1:ncol,1:nlev),INTENT(INOUT) ::    QV   &
     &                                                            ,QC   &
     &                                                            ,QI   &
     &                                                            ,QR   &
     &                                                            ,QS   &
     &                                                            ,QG
      REAL(kind_phys),DIMENSION(1:ncol),INTENT(INOUT) ::                &
     &                                                 RLWIN,RLWTOA     &
     &                                                ,RSWIN,RSWOUT     &
     &                                                ,GSW              &
     &                                                ,RSWINC
      REAL(kind_phys),DIMENSION(1:ncol,1:nlev),INTENT(INOUT) ::         &
     &                                                         RLWTT    &
     &                                                        ,RSWTT
      REAL(kind_phys),DIMENSION(1:ncol),INTENT(INOUT) ::                &
     &                                                        CZMEAN    
      REAL(kind_phys),DIMENSION(1:ncol,1:nlev),INTENT(INOUT) :: CLDFRA
      REAL(kind_phys), DIMENSION(1:ncol), INTENT(INOUT) ::              &
     &                      ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC,          &
     &                      ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,          &
     &                      ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC,          &
     &                      ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC

      ! TOA and surface, upward and downward, total and clear fluxes
      REAL(kind_phys), DIMENSION(1:ncol), OPTIONAL, INTENT(INOUT) ::    &
     &         SWUPT,  SWUPTC, SWUPTCLN,  SWDNT,  SWDNTC, SWDNTCLN,     &
     &         SWUPB,  SWUPBC, SWUPBCLN,  SWDNB,  SWDNBC, SWDNBCLN,     &
     &         LWUPT,  LWUPTC, LWUPTCLN,  LWDNT,  LWDNTC, LWDNTCLN,     &
     &         LWUPB,  LWUPBC, LWUPBCLN,  LWDNB,  LWDNBC, LWDNBCLN

      !mz-new added  ----
      ! Upward and downward, total and clear sky layer fluxes (W m-2)
!      REAL(KIND_PHYS), DIMENSION( 1:ncol, 1:nlev+2 ),                   &
!     &   OPTIONAL, INTENT(INOUT) ::                                     &
!     &                         SWUPFLX,SWUPFLXC,SWDNFLX,SWDNFLXC,       &
!     &                         LWUPFLX,LWUPFLXC,LWDNFLX,LWDNFLXC
      ! ---- fds (06/2010) ssib alb components ------------
!      REAL(KIND_PHYS), DIMENSION( 1:ncol ),            OPTIONAL ,       &
!     &   INTENT(IN   )  ::                                  ALSWVISDIR, &
!     &                                                      ALSWVISDIF, &
!     &                                                      ALSWNIRDIR, &
!     &                                                       ALSWNIRDIF
      ! ---- fds (06/2010) ssib swr components ------------
      REAL(kind_phys), DIMENSION(1:ncol),  OPTIONAL,                    &
     &                                    INTENT(OUT  )  ::   SWVISDIR, &
     &                                                        SWVISDIF, &
     &                                                        SWNIRDIR, &
     &                                                        SWNIRDIF
      !..Additions for coupling cloud physics effective radii and radiation.  
      REAL(kind_phys),DIMENSION(1:ncol,1:nlev),INTENT(INOUT):: re_cloud,&
     &                                                      re_ice,     &
     &                                                      re_snow
      INTEGER, INTENT(INOUT):: has_reqc, has_reqi, has_reqs

      !..Output daily longwave and shortwave heating rates: 
      REAL(kind_phys),DIMENSION(1:ncol,1:nlev),INTENT(OUT):: hrlwpd,    &
     &                                                       hrswpd
      integer,                   intent(in)    :: mpirank
      integer,                   intent(in)    :: mpiroot
      INTEGER,               INTENT(IN   )    :: mpicomm
      character(len=*),          intent(  out) :: errmsg
      integer,                   intent(  out) :: errflg

!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER :: I,IENDX,II,ISTAT,J,JDAY,JMONTH,K,KMNTH,N,NRAD,         &
     &           NUM_AEROSOLC
      INTEGER,DIMENSION(3) :: IDAT
      INTEGER,DIMENSION(12) :: MONTH=(/31,28,31,30,31,30,31,31          &
     &                                ,30,31,30,31/)
!
      REAL(kind_phys) :: CAPA,DAYI,FICE,FRAIN,GMT,HOUR,PLYR,            &
     &                   QW,RADT,TIMES,WC,TDUM,XTIME
!
      REAL(kind_phys),DIMENSION(1:nlev-1) :: QL,TL
!
      REAL(kind_phys),DIMENSION(1:ncol) :: CZEN                         &
     &                                  ,SWCF, LWCF                     &
     &                                  ,REXNSFC,SWNETDN                &
     &                                  ,TOT,TOTLWDN,TOTSWDN,TOTSWDNC   &
     &                                  ,XLAND
      REAL(kind_phys),DIMENSION(1:ncol,1:nlev) :: CLFR,DZ               &
     &                                          ,P_PHY,PI_PHY           &
     &                                          ,RR,T8W                 &
     &                                          ,THRATENLW,THRATENSW    &
     &                                          ,TH_PHY,T_PHY
      REAL(kind_phys) :: DXKM, DYKM
      REAL(KIND_PHYS), DIMENSION(1:ncol) :: sm
      integer :: ntsd_rad

      INTEGER            :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,myis,myis1,myie,myie1                       &
     &                     ,myjs,myjs2,myje,myje2

      LOGICAL, PARAMETER             :: allowed_to_read =.true.

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      !mz 
      if(mpirank == mpiroot) write(0,*)'start HWRF_radiation_run '
      if(mpirank == mpiroot) write(0,*)'ntsd, nrads,nradl = ',          &
     &   ntsd, nrads, nradl       
      ! Check initialization state
      if (.not.acswalloc .or. .not. aclwalloc) then
          write(errmsg, fmt='((a))') 'HWRF_radiation_run called before  &
     &                                HWRF_radiation_init'
          errflg = 1
          return
      end if

      ! Set internal dimensions
      ids    =    1
      ims    =    1
      its    =    1
      myis   =    1
      myis1  =    1
      ide    =    ncol
      ime    =    ncol
      ite    =    ncol
      myie   =    ncol
      myie1  =    ncol
      jds    =    1
      jms    =    1
      jts    =    1
      myjs   =    1
      myjs2  =    1
      jde    =    1
      jme    =    1
      jte    =    1
      myje   =    1
      myje2  =    1
      kds    =    1
      kms    =    1
      kts    =    1
      kde    =    nlev
      kme    =    nlev
      kte    =    nlev

!      CALL rrtmg_lwinit(                                                &
!     &               allowed_to_read,                                   &
!     &               ids, ide, jds, jde, kds, kde,                      &
!     &               ims, ime, jms, jme, kms, kme,                      &
!     &               its, ite, jts, jte, kts, kte,                      &
!     &               mpirank, mpiroot, mpicomm,                         &
!     &               errmsg, errflg       )

!      if (errflg /= 0) return
!      if(mpirank ==mpiroot)write (0,*)'HWRF rrtmg_lwinit completed...'

!      aclwalloc = .true.

!      CALL rrtmg_swinit(                                                &
!     &               allowed_to_read ,                                  &
!     &               ids, ide, jds, jde, kds, kde,                      &
!     &               ims, ime, jms, jme, kms, kme,                      &
!     &               its, ite, jts, jte, kts, kte,                      &
!     &               mpirank, mpiroot, mpicomm,                         &
!     &               errmsg, errflg )     
!      if (errflg /= 0) return                         
!      if (mpirank == mpiroot)write(0,*)'HWRF rrtmg_swinit completed'   
!      acswalloc = .true.                 
!----------------------------------------------------------------------
!***  RADIATION
!----------------------------------------------------------------------
!
!***  When allocating CAM radiation 4d arrays (ozmixm, aerosolc), 
!***  the following two scalars are not needed.
!
      NUM_AEROSOLC=1
!
!MZ      IF(grid%ntsd<=0)THEN
!MZ        NTSD_rad=grid%ntsd
!MZ      ELSE
!
!***  Call radiation just BEFORE the top of the hour
!***  so that updated fields are written to history files.
!
!mz
      IF (NTSD > 0) then
        NTSD_rad=ntsd+1
      ENDIF

!mz*:minute_since_simulation_start
      XTIME = ntsd*dt/60.

      if(mpirank == mpiroot) then
         write(0,*)'ntsd,ntsd_rad, nrads,nradl,xtime = ',               &
     &   ntsd,ntsd_rad, nrads, nradl, xtime
      endif

!mz      IF(MOD(NTSD_rad,NRADS)==0.OR. MOD(NTSD_rad,NRADL)==0)THEN

!
!-----------------------------------------------------------------------
!***** NOTE: THIS IS HARDWIRED FOR CALLS TO LONGWAVE AND SHORTWAVE
!*****       AT EQUAL INTERVALS
!-----------------------------------------------------------------------
!
      NRAD=NRADS
      RADT=DT*NRADS/60.
!
!mz
      DO I = 1,NCOL
         if (slmsk(i) == 0.) then !sea
           sm(i) = 1.
          else
           sm(i) =0.   !land or ice (mz: is it true?)
          endif
      ENDDO

!-----------------------------------------------------------------------
!
      CAPA=R_D/CP
!mz: wrf-nmm
!      DXKM=grid%dlmd*0.01745329*6371200. ! numbers from module_initialize_real.F
!      DYKM=grid%dy_nmm
! In FV3, dx in meter
!mz: use the dx of the 1st i point to get a value of dxkm and dykm for
!cloud fraction calculation
       DXKM = dx(1)*0.001
       DYKM = dx(1)*0.001
!
!-----------------------------------------------------------------------
!
!mz      DO J=MYJS2,MYJE2
!mz      DO I=MYIS1,MYIE1

      DO I=1,NCOL
!
!mz        PDSL(I,J)=PD(I,J)*RES(I,J)
!mz        P8W(I,KTE+1,J)=PT
!        XLAND(I,J)=SM(I,J)+1.
         XLAND(I)=SM(I)+1.
!mz        PSFC=PD(I,J)+PDTOP+PT
!mz        REXNSFC(I,J)=(PSFC*1.E-5)**CAPA
!mz        TSFC(I,J)=THS(I,J)*REXNSFC(I,J)
!        T8W(I,KTS,J)=TSFC(I,J)
        T8W(I,1) = TSFC(I)
!MZ:eta1-interface sigma value in pressure domain 
!   eta2-interface sigma value in sigma domain
!  pdtop-mass at i,j in pressure domain, Pa
!  pdsl- sigma-domain pressure at sigma=1
!MZ        P8W(I,KTS,J)=ETA1(KTS)*PDTOP+ETA2(KTS)*PDSL(I,J)+PT
!prsi: air pressure at interface in Pa
!        P8W(I,KTS,J)=prsi(i,kts)
     
!        Z_PHY(I,KTS,J)=Z(I,J,KTS)
!        HT(I,J)=Z(I,J,KTS)
      ENDDO
!      ENDDO
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
      !DO J=MYJS2,MYJE2
      DO I = 1,NCOL      !MYIS1,MYIE1
        DO K =1,NLEV     !KTS,KTE
          !mz DPL=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
          !! QL(K)=MAX(Q(I,J,K),EPSQ)
          QL(K) = MAX(Q(I,K),EPSQ)
          !mz PLYR=AETA1(K)*PDTOP+AETA2(K)*PDSL(I,J)+PT
          !!PLYR = PRSL(I,J,K)
          PLYR = PRSL(I,K)
          !!TL(K)=T(I,J,K)
          TL(K)=T(I,K)
!
          !!RR(I,K,J)=PLYR/(R_D*TL(K)*(1.+P608*QL(K)))
          RR(I,K)=PLYR/(R_D*TL(K)*(1.+P608*QL(K)))
          T_PHY(I,K)=TL(K)
          TH_PHY(I,K)=TL(K)*(1.E5/PLYR)**CAPA
          !mz P8W(I,K+1,J)=ETA1(K+1)*PDTOP+ETA2(K+1)*PDSL(I,J)+PT
          P_PHY(I,K)=PLYR
          PI_PHY(I,K)=(PLYR*1.E-5)**CAPA
          DZ(I,K)=TL(K)*(P608*QL(K)+1.)*R_D                             &
     &                 *(P8W(I,K)-P8W(I,K+1))                           &
     &                 /(P_PHY(I,K)*G)
!
          RTHRATEN(I,K)=0.
          THRATENLW(I,K)=0.
          THRATENSW(I,K)=0.

        ENDDO
!
        DO K=2, NLEV    !KTS+1,KTE
          T8W(I,K)=0.5*(TL(K-1)+TL(K))
        ENDDO
!        T8W(I,KTE+1,J)=-1.E20 
! For RRTM 
        T8W(I,KTE+1)=T8W(I,KTE) + 0.5*(T8W(I,KTE)-T8W(I,KTE-1))
!
      ENDDO
!      ENDDO
!
      GMT=REAL(IHRST)
!
      DO K=1,NLEV 
      !  DO J=JMS,JME
      DO I=1,NCOL 
          CLDFRA(I,K)=0.
      ENDDO
      ENDDO
!
        DO I=1, ncol 
!          CFRACH(I,J)=0.
!          CFRACL(I,J)=0.
!          CFRACM(I,J)=0.
          CZMEAN(I)=0.
!          SIGT4(I,J)=0.
          TOTSWDN(I)=0.   ! TOTAL (clear+cloudy sky) shortwave down at the surface
          TOTSWDNC(I)=0.  ! CLEAR SKY shortwave down at the surface
          SWNETDN(I)=0.   ! Net (down - up) total (clear+cloudy sky) shortwave at the surface
          TOTLWDN(I)=0.   ! Total longwave down at the surface
!mz          CUPPTR(I,J)=CUPPT(I,J)   ! Temporary array set to zero in radiation
!mz          HTOPR(I,J) =0.
!mz          HBOTR(I,J) = REAL(KTE+1)
!          HBOTR(I,J) =0.
          SWCF(I) =0.
          LWCF(I) =0.
!
!
        ENDDO
!
!-------------------------------------------------------------------
!
!***  CALL THE INNER DRIVER.
!
!-----------------------------------------------------------------------

      if (mpirank == mpiroot) then
          write(0,*)'mz: NRAD, RADT =', NRAD,RADT
     !     write(0,*)'mz: max/min(dx) =', maxval(dx), minval(dx)
     !     write(0,*)'mz: max/min(sm) =', maxval(sm), minval(sm)
     !     write(0,*)'mz: max/min(tsfc) =', maxval(tsfc), minval(tsfc)
          write(0,*)'mz: max/min(Radtend%sfalb) =',                     &
     &               maxval(Radtend%sfalb), minval(Radtend%sfalb)
          write(0,*)'mz: max/min(Radtend%coszen) =',                    &
     &               maxval(Radtend%coszen), minval(Radtend%coszen)
          write(0,*)'mz: max/min(Radtend%semis) =',                     &
     &               maxval(Radtend%semis), minval(Radtend%semis)
      endif


      CALL RADIATION_DRIVER (ALBEDO=Radtend%sfalb                       &
     &                 ,CZMEAN=CZMEAN ,DT=DT                            &
     &                 ,DZ8W=DZ,EMISS=Radtend%semis,GLW=TOTLWDN, GMT=GMT&
     &                 ,GSW=SWNETDN                                     &
     &                 ,ITIMESTEP=NTSD_rad ,JULDAY=JULDAY               &
     &                 ,JULIAN=JULIAN,JULYR=JULYR                       &
     &                 ,NPHS=NPHS, O3RAD=o3rad                          &
     &                 ,O3INPUT=O3INPUT,AER_OPT=AER_OPT                 &
     &                 ,swint_opt=swint_opt                             &
     &                 ,RA_CALL_OFFSET=RA_CALL_OFFSET                   &
     &                 ,ICLOUD=ICLOUD,cldovrlp=cldovrlp                 &
     &                 ,SF_SURFACE_PHYSICS=SF_SURFACE_PHYSICS           &
     &                 ,P8W=P8W,P=P_PHY,PI=PI_PHY,RADT=RADT             &
     &                 ,RHO=RR, RLWTOA=RLWTOA ,RTHRATEN=RTHRATEN        &
     &                 ,RTHRATENLW=THRATENLW,RTHRATENSW=THRATENSW       &
     &                 ,HRSWPD=HRSWPD, HRLWPD=HRLWPD                    &
     &                 ,SNOW=SNOW ,STEPRA=NRAD ,SWDOWN=TOTSWDN          &
     &                 ,SWDOWNC=TOTSWDNC                                &
     &                 ,T8W=T8W ,T=T_PHY,TSK=TSFC                       &
     &                 ,XICE=SICE,XLAND=XLAND                           &
     &                 ,XLAT=XLAT,XLONG=XLONG                           &
     &                 ,YR=JULYR                                        &
     &                 ,sinlat=sinlat                                   &
     &                 ,coslat=coslat                                   &
     &                 ,solhr= solhr                                    &
     &                 ,coszen=Radtend%coszen                           &
     &                 ,solcon = solcon                                 &
!mz     &                 ,hrang=hrang                                     &
     &                 ,ALEVSIZ=ALEVSIZ,no_src_types=no_src_types       &
     &                 ,LEVSIZ=LEVSIZ,N_OZMIXM=NUM_OZMIXM               &
     &                 ,N_AEROSOLC=NUM_AEROSOLC,PAERLEV=PAERLEV         &
     &                 ,XTIME=XTIME                                     &
     &                 ,ITS=ITS,ITE=ITE                                 &
     &                 ,JTS=JTS,JTE=JTE                                 &
     &                 ,IDS=IDS,IDE=IDE,JDS=JDS                         &
     &                 ,JDE=JDE,KDS=KDS,KDE=KDE                         &
     &                 ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=KMS,KME=KME &
     &                 ,KTS=KTS,KTE=KTE                                 &
     &                 ,CLDFRA=CLFR                                     &
     & ,re_cloud=re_cloud ,re_ice=re_ice ,re_snow=re_snow               &
     & ,has_reqc=has_reqc ,has_reqi=has_reqi, has_reqs=has_reqs         &
     &                 ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN               &
     &          ,QV=QV ,QC=QC ,QR=QR ,QI=QI ,QS=QS ,QG=QG               &
     &          ,ACSWUPT=ACSWUPT ,ACSWUPTC=ACSWUPTC                     &
     &          ,ACSWDNT=ACSWDNT ,ACSWDNTC=ACSWDNTC                     &
     &          ,ACSWUPB=ACSWUPB ,ACSWUPBC=ACSWUPBC                     &
     &          ,ACSWDNB=ACSWDNB ,ACSWDNBC=ACSWDNBC                     &
     &          ,ACLWUPT=ACLWUPT ,ACLWUPTC=ACLWUPTC                     &
     &          ,ACLWDNT=ACLWDNT ,ACLWDNTC=ACLWDNTC                     &
     &          ,ACLWUPB=ACLWUPB ,ACLWUPBC=ACLWUPBC                     &
     &          ,ACLWDNB=ACLWDNB ,ACLWDNBC=ACLWDNBC                     &
     &          ,SWUPT=SWUPT,SWUPTC=SWUPTC,SWUPTCLN=SWUPTCLN            &
     &          ,SWDNT=SWDNT,SWDNTC=SWDNTC,SWDNTCLN=SWDNTCLN            &
     &          ,SWUPB=SWUPB,SWUPBC=SWUPBC,SWUPBCLN=SWUPBCLN            &
     &          ,SWDNB=SWDNB,SWDNBC=SWDNBC,SWDNBCLN=SWDNBCLN            &
     &          ,LWUPT=LWUPT,LWUPTC=LWUPTC,LWUPTCLN=LWUPTCLN            &
     &          ,LWDNT=LWDNT,LWDNTC=LWDNTC,LWDNTCLN=LWDNTCLN            &
     &          ,LWUPB=LWUPB,LWUPBC=LWUPBC,LWUPBCLN=LWUPBCLN            &
     &          ,LWDNB=LWDNB,LWDNBC=LWDNBC,LWDNBCLN=LWDNBCLN            &
     &                 ,LWCF=LWCF,SWCF=SWCF                             &
     &                 ,DX=DXKM,DY=DYKM                                 &
     &                 ,OZMIXM=OZMIXM,PIN=PIN                           &
     &         ,CALC_CLEAN_ATM_DIAG=CALC_CLEAN_ATM_DIAG                 &
!mz new
!     &         ,SWUPFLX=swupflx,SWUPFLXC=swupflxc                       &
!     &         ,SWDNFLX=swdnflx,SWDNFLXC=swdnflxc                       & ! Optional
!     &         ,LWUPFLX=lwupflx,LWUPFLXC=lwupflxc                       &
!     &         ,LWDNFLX=lwdnflx,LWDNFLXC=lwdnflxc                       & ! Optional
!     &         ,ALSWVISDIR=alswvisdir, ALSWVISDIF=alswvisdif            &
!     &         ,ALSWNIRDIR=alswnirdir, ALSWNIRDIF=alswnirdif            & !fds ssib alb comp (06/2010)
!mz new
     &        ,SWVISDIR=swvisdir ,SWVISDIF=swvisdif                     &
     &        ,SWNIRDIR=swnirdir ,SWNIRDIF=swnirdif                     &
     &                 ,SWDDIR=swddir,SWDDNI=swddni,SWDDIF=swddif       &
     &                 ,SWDDIRC=swddirc,SWDDNIC=swddnic                 &
     &                 ,MP_PHYSICS=MP_PHYSICS                           &
     &                 ,imp_physics_fer_hires = imp_physics_fer_hires   &
     &        ,MPIRANK=mpirank, MPIROOT=mpiroot                         &
     &        ,ERRFLG=errflg, ERRMSG=errmsg )


!
!-----------------------------------------------------------------------
!
!***  UPDATE FLUXES AND TEMPERATURE TENDENCIES.
!
!-----------------------------------------------------------------------
!***  SHORTWAVE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!mz      nrads_block: IF(MOD(NTSD_rad,NRADS)==0)THEN
!-----------------------------------------------------------------------
!***  COMPUTE CZMEAN FOR NON-GFDL SHORTWAVE
!-----------------------------------------------------------------------
!mz          DO I=1,ncol
         
!            CZMEAN(I)=0.
!            TOT(I)=0.
!          ENDDO
!
!          CALL CAL_MON_DAY(JULDAY,JULYR,JMONTH,JDAY)
!          IDAT(1)=JMONTH
!          IDAT(2)=JDAY
!          IDAT(3)=JULYR
!
!          DO II=0,NRADS,NPHS
!            TIMES=NTSD_rad*DT+II*DT
!            CALL ZENITH(TIMES,DAYI,HOUR,IDAT,IHRST,XLONG,XLAT,CZEN      &
!     &                 ,MYIS,MYIE,MYJS,MYJE                             &
!     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                         &
!     &                 ,IMS,IME,JMS,JME,KMS,KME                         &
!     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)
!
!mz!$omp parallel do                                                       &
!mz!$omp& private(i,j)
           ! DO J=MYJS,MYJE
!            DO I= 1,ncol !MYIS,MYIE
!              IF(CZEN(I)>0.)THEN
!                CZMEAN(I)=CZMEAN(I)+CZEN(I)
!                TOT(I)=TOT(I)+1.
!              ENDIF
!            ENDDO
!
!          ENDDO
!
!mz!$omp parallel do                                                       &
!mz!$omp& private(i,j)
          !DO J=MYJS,MYJE
!          DO I=1,ncol  !MYIS,MYIE
!            IF(TOT(I)>0.)CZMEAN(I)=CZMEAN(I)/TOT(I)
!          ENDDO
          !ENDDO
!
!-----------------------------------------------------------------------
!***  COMPUTE TOTAL SFC SHORTWAVE DOWN FOR NON-GFDL SCHEMES
!-----------------------------------------------------------------------
!
!mz!$omp parallel do                                                       &
!mz!$omp& private(i,j)
          !DO J=MYJS2,MYJE2
          DO I=1,ncol !YIS1,MYIE1
!
!MZ            IF(HBM2(I,J)>0.5)THEN
!              TOTSWDN(I)=SWNETDN(I)/(1.-ALBEDO(I))
              TOTSWDN(I)=SWNETDN(I)/(1.-Radtend%sfalb(I))
!
!--- No value currently available for clear-sky solar fluxes from
!    non GFDL schemes, though it's needed for air quality forecasts.
!    For the time being, set to the total downward solar fluxes.
!
              TOTSWDNC(I)=TOTSWDN(I)
!MZ            ENDIF
!
          ENDDO
          !ENDDO
!
!        ENDIF   !End non-GFDL block
!-----------------------------------------------------------------------
!
!mz!$omp parallel do                                                       &
!mz!$omp& private(i,iendx,j)
        !DO J=MYJS2,MYJE2
        !  IENDX=MYIE1
        !  IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=1,ncol   !MYIS1,IENDX
!
            RSWIN(I)=TOTSWDN(I)
            RSWINC(I)=TOTSWDNC(I)
            RSWOUT(I)=TOTSWDN(I)-SWNETDN(I)
!
          ENDDO
        !ENDDO
!
!mz!$omp parallel do                                                       &
!mz!$omp& private(i,iendx,j,k)
        !DO J=MYJS2,MYJE2
        !  IENDX=MYIE1
        !  IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=1,ncol !MYIS1,IENDX
            DO K=1, nlev !KTS,KTE
              RSWTT(I,K)=THRATENSW(I,K)*PI_PHY(I,K)
            ENDDO
!
          ENDDO
!        ENDDO
!
!mz      ENDIF nrads_block
!
!-----------------------------------------------------------------------
!***  LONGWAVE
!-----------------------------------------------------------------------
!
!mz      nradl_block: IF(MOD(NTSD_rad,NRADL)==0)THEN
!
!$omp parallel do                                                       &
!$omp& private(i,iendx,j)
        !DO J=MYJS2,MYJE2
        !  IENDX=MYIE1
        !  IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=1,ncol !MYIS1,IENDX
!
!MZ: hbm2=1. in global FV3
!            IF(HBM2(I,J)>0.5)THEN
              TDUM=T(I,1)
!              SIGT4(I,J)=STBOLT*TDUM*TDUM*TDUM*TDUM
              RLWIN(I)=TOTLWDN(I)
!            ENDIF
!
          ENDDO
        !ENDDO
!
!mz$omp parallel do                                                       &
!mz$omp& private(i,iendx,j,k)
        !DO J=MYJS2,MYJE2
        !  IENDX=MYIE1
        !  IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
!
          DO K=1,nlev !KTS,KTE
          DO I=1, ncol !MYIS1,IENDX
!MZ            IF(HBM2(I,J)>0.5)THEN
                RLWTT(I,K)=THRATENLW(I,K)*PI_PHY(I,K)
!MZ            ENDIF
          ENDDO
          ENDDO
!
!        ENDDO
!
!      ENDIF nradl_block
!
!-----------------------------------------------------------------------
!***  STORE 3D CLOUD FRACTIONS.
!-----------------------------------------------------------------------
!
!mz!$omp parallel do                                                       &
!mz!$omp& private(i,iendx,j,k)
      DO K=1, nlev !KTS,KTE
        !DO J=MYJS2,MYJE2
        !  IENDX=MYIE1
        !  IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=1,ncol !MYIS1,IENDX
            CLDFRA(I,K)=CLFR(I,K)
          ENDDO
        !ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  RESET THE DIAGNOSTIC CONVECTIVE CLOUD TOPS/BOTTOMS AFTER
!***  EACH RADIATION CALL.
!-----------------------------------------------------------------------
!
!mz!$omp parallel do                                                       &
!mz!!$omp& private(i,iendx,j)
      !DO J=MYJS2,MYJE2
      !  IENDX=MYIE1
      !  IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
!        DO I=1, ncol !MYIS1,IENDX
!          HBOT(I)=HBOTR(I)
!          HTOP(I,J)=HTOPR(I,J)
!          CUPPT(I,J)=CUPPTR(I,J)
!        ENDDO
      !ENDDO
!
!-----------------------------------------------------------------------
!***  ZERO OUT BOUNDARY ROWS.
!-----------------------------------------------------------------------
!
!MZ
!      DO J=JTS,JTE
!      DO I=ITS,ITE
!        IF(HBM2(I,J)<0.5)THEN
!          ACFRST(I,J)=0.
!          ACFRCV(I,J)=0.
!!          CFRACL(I,J)=0.
!!          CFRACM(I,J)=0.
!!          CFRACH(I,J)=0.
!          RSWTOA(I,J)=0.
!          RLWTOA(I,J)=0.
!        ENDIF
!      ENDDO
!      ENDDO

!
      !DO J=jts,min(jde-1,jte)
      !DO I=its,min(ide-1,ite)
      do i=1,ncol
          gsw(I)=rswin(I)-rswout(I)
      ENDDO
      !ENDDO
     
!mz      ENDIF    !radiation loop

!
!-----------------------------------------------------------------------
!
      END SUBROUTINE HWRF_radiation_run
      END MODULE HWRF_radiation

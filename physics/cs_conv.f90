module cs_conv
!---------------------------------------------------------------------------------
! Purpose:
!
!>---------------------------------------------------------------------------------
! Purpose:
!
! Interface for Chikira-Sugiyama convection scheme 
!
! Author: Minoru Chikira
! History:
!  June 26 2014  D. Dazlich - Modified for GFS
!  Apr 10 2015 : S. Moorthi - check for allocatable arrays and fix argument for cbmfx
!  Oct    2015 : D. Dazlich - Add computation of updraft area fraction (sigma) for
!                             diagnostic purposes.
!  Aug    2016 : D. Dazlich - Create flux form of tendencies and multiply by
!                             Arakawa-Wu functions of sigma
!  Sep    2016 : S. Moorthi - found two bugs - cleanup and some optimization
!  Oct    2016 : S. Moorthi - added sigma affects on tracers and CUMFLX and CUMDET
!                             made many cosmetic changes
!  Nov    2016 : S. Moorthi - further optimization and cleanup and several bug fixes
!
!  Arakawa-Wu implemtation: for background, consult An Introduction to the 
!      General Circulation of the Atmosphere, Randall, chapter six.
!    Traditional parameterizations compute tendencies like those in eq 103, 105 and 106.
!    Because Arakawa-Wu applies different functions to different components to the
!    terms within these equations, it requires the terms used in alternate eqns 91 - 93.
!    The code required to compute these terms is added within, and the appropriate
!    functions of updraft area fraction (sigma) are applied. Thus, AW requires three
!    steps:
!       computation of the updraft area fraction
!       alternative representation of the tendency terms
!       application of functions of sigma to the alternative tendency terms
!       here, and in gbphys to the large-scale microphysics tendencies.
!
!  The bulk of AW is implemented within subroutine CS_CUMLUS, and the routines it calls.
!
!---------------------------------------------------------------------------------
!
  use machine ,   only : r8    => kind_phys
  use physcons,   only : cp    => con_cp,   grav   => con_g,                   &
     &                   rair  => con_rd,   rvap   => con_rv,                  &
     &                   cliq  => con_cliq, cvap   => con_cvap,                &
     &                   epsv  => con_eps,  epsvm1 => con_epsm1,               &
     &                   epsvt => con_fvirt,                                   &
     &                   el    => con_hvap, emelt  => con_hfus, t0c => con_t0c
  use funcphys, only : fpvs ! this is saturation vapor pressure in funcphys.f

  
  implicit none

  private                ! Make default type private to the module

   real(r8), parameter :: zero=0.0d0, one=1.0d0, half=0.5d0
   real(r8), parameter :: cpoel=cp/el, cpoesub=cp/(el+emelt), esubocp=1.0/cpoesub, &
                          elocp=el/cp, oneocp=one/cp, gocp=grav/cp, gravi=one/grav,&
                          emeltocp=emelt/cp, cpoemelt=cp/emelt
   real(r8), parameter :: fact1=(cvap-cliq)/rvap, fact2=el/rvap-fact1*t0c
   logical,  parameter :: adjustp=.true.

! Tuning parameters set from namelist
!
!  real(r8), save, public :: CLMD = 0.6,    & ! entrainment efficiency
   real(r8), save, public :: CLMD = 0.7,    & ! entrainment efficiency
                             PA=0.15,       & ! factor for buoyancy to affect updraft velocity
                             CPRES = 0.55,  & ! pressure factor for momentum transport
                             ALP0 = 8.0e7     ! alpha parameter in prognostic closure

!DD next two parameters control partitioning of water between detrainment
!DD   and precipitation. Decrease for more precip
!M REAL(r8), public, save ::  PRECZ0 = 1.5e3_r8   ! default = 1.5e3
!M REAL(r8), public, save ::  PRECZ0 = 1.5e3_r8   ! default = 1.5e3
!  REAL(r8), public, save ::  PRECZ0 = 1.5e3_r8   ! default = 1.5e3
!  REAL(r8), public, save ::  PRECZH = 4.e3_r8    ! default = 4.e3

!  REAL(r8), public, save ::  PRECZ0 = 1.0e3_r8   ! default = 1.5e3
!  REAL(r8), public, save ::  PRECZH = 3.e3_r8    ! default = 4.e3

!  REAL(r8), public, save ::  PRECZ0 = 0.5e3_r8   ! default = 1.5e3
!  REAL(r8), public, save ::  PRECZH = 2.e3_r8    ! default = 4.e3

   real(r8), public       ::  precz0, preczh
!
! Private data
!
  real(r8), parameter         :: unset_r8 = -999._r8   ! missing value
!
  integer :: iulog ! unit to write debugging and diagnostic output
                    !DD Note - see if I can find corresponding variable in a GFS module
!
! Shared variables
!
  integer, parameter :: ITI = 2, ITL = 3  ! index of ice and liquid water

! logical            :: outputflag(100)

!DD  integer, save :: ICHNK        ! chunk identifier

!   [INTERNAL PARM]   !DD moved to module scope and allocatable

! logical, save, dimension(50) :: OTSPT1, OTSPT2
  integer, save, dimension(50) :: IMFXR   ! 0: mass fixer is not applied
                                          !    tracers which may become negative
                                          !    values e.g. subgrid-PDFs
                                          ! 1: mass fixer is applied, total mass
                                          !    may change through cumulus scheme
                                          !    e.g. moisture, liquid cloud, ice
                                          !    cloud, aerosols
                                          ! 2: mass fixer is applied, total mass
                                          !    never change through cumulus scheme
                                          !    e.g. CO2

! LOGICAL,   SAVE, ALLOCATABLE, DIMENSION(:) :: OTSPT1   ! tracer transport by updraft, downdraft on/off
                                                         ! should not include subgrid PDF and turbulence
!  LOGICAL,  SAVE, ALLOCATABLE, DIMENSION(:) :: OTSPT2   ! tracer transport by subsidence on/off
                                                         ! should include subgrid PDF and turbulence
!  INTEGER,  SAVE, ALLOCATABLE, DIMENSION(:) :: IMFXR 
!  REAL(r8), SAVE, ALLOCATABLE, DIMENSION(:) :: FSCAV    !DD    split declaration and initialization
!  REAL(r8), SAVE, ALLOCATABLE, DIMENSION(:) :: FSWTR    !DD    split declaration and initialization
!
!
!

!  PUBLIC: interfaces
!
   public cs_convr        ! CS scheme main driver
  
   contains

!---------------------------------------------------------------------------------
! use GFS functions
   function FQSAT( T, P )   ! calculate saturation water vapor 

   implicit none
  
   real(r8)             :: FQSAT  ! saturation water vapor
   real(r8), intent(in) :: T      ! temperature [K]
   real(r8), intent(in) :: P      ! pressure [Pa]

!  real(r8), parameter  :: one_m10=1.0d-10,   &
!                          ES0   = 611._r8,   &   ! saturation e at 0 deg C (Pa)
!                          TQICE = 273.15_r8, &   ! T threshold for ice QSAT
!                          TMELT = 273.15_r8      ! melting point of water
  
!DD  FQSAT = EPSV * ES0 / P &
!DD        * EXP( (EL+EMELT/2._r8*(1._r8-SIGN(1._r8,T-TQICE))) &
!DD               /RVAP *( 1._r8/TMELT - 1._r8/T )           )

   FQSAT = min(p,fpvs(T))      !DD this is saturation vapor pressure

!  FQSAT = EPSV * FQSAT / P    !DD This is saturation mixing ratio
!  FQSAT = EPSV * FQSAT / (max(p+epsvm1*fqsat,ONE_M10))  !DD&Moo This is saturation specific humidity
   FQSAT = min(EPSV*FQSAT/max(p+epsvm1*fqsat,1.0e-10), 1.0)  !DD&Moo This is saturation specific humidity

   end function FQSAT
!---------------------------------------------------------------------------------
!  following GFS
   function FDQSAT( T, QS )   ! calculate d(qs)/dT

   implicit none
  
   real(r8)             :: FDQSAT ! d(QSAT)/d(T)
   real(r8), intent(in) :: T      ! temperature [K]
   real(r8), intent(in) :: QS     ! saturation water vapor [kg/kg]
   real(r8)             :: wrk
  
   real(r8), parameter  :: fact1=(cvap-cliq)/rvap,fact2=el/rvap-fact1*t0c
  
!DD  FDQSAT = (EL+EMELT/2._r8*(1._r8-SIGN(1._r8,T-TMELT))) &
!DD         * QS / ( RVAP * T*T )

   wrk    = 1.0 / t
   FDQSAT = qs * wrk * (fact1 + fact2*wrk)
!  FDQSAT = qs * (fact1 / t + fact2 / (t**2))


   end function FDQSAT
!---------------------------------------------------------------------------------
   subroutine cs_convr(IM     , IJSDIM ,  KMAX     , NTR     , nctp,     & !DD dimensions
                       otspt  , lat    ,kdt    ,                         &
                       t      , q      ,  prec     , clw     ,           &
                       zm     , zi     ,  pap      , paph    ,           &
                       delta  , delti  ,  ud_mf    , dd_mf   , dt_mf,    &
                       u      , v      ,  fscav    , fswtr,              &
                       cbmfx  , mype   ,  wcbmaxm  , precz0in, preczhin, &
                       sigmai , sigma  ,  vverti   , do_aw, do_awdd,     &
                       lprnt, ipr,                                       &
! for coupling to Morrison microphysics
                       QLCN, QICN, w_upi, cf_upi, CNV_MFD, CNV_PRC3,     &
                       CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,ncld)

!---------------------------------------------------------------------------------
! Purpose:
!
! Main driver for Chikira-Sugiyama convective scheme 
!
! Author: Minoru Chikira
!
!---------------------------------------------------------------------------------

   implicit none
!
! input arguments
!
   INTEGER, INTENT(IN)     :: IM,IJSDIM, KMAX, NTR, mype, nctp, ncld, kdt,lat !! DD, for GFS, pass in
   logical, intent(in)     :: otspt(ntr,2)


   real(r8), intent(inout) :: t(IM,KMAX)          ! temperature at mid-layer (K)
   real(r8), intent(inout) :: q(IM,KMAX)          ! water vapor array including moisture (kg/kg)
   real(r8), intent(inout) :: clw(IM,KMAX,ntr-1)  ! tracer array including cloud condensate (kg/kg)
   real(r8), intent(in)    :: pap(IM,KMAX)        ! pressure at mid-layer (Pa)
   real(r8), intent(in)    :: paph(IM,KMAX+1)     ! pressure at boundaries (Pa)
   real(r8), intent(in)    :: zm(IM,KMAX)         ! geopotential at mid-layer (m)
   real(r8), intent(in)    :: zi(IM,KMAX+1)       ! geopotential at boundaries (m)
   real(r8), intent(in)    :: fscav(ntr), fswtr(ntr), wcbmaxm(ijsdim)
   real(r8), intent(in)    :: precz0in, preczhin
! added for cs_convr
   real(r8), intent(inout) :: u(IM,KMAX)          ! zonal wind at mid-layer (m/s)
   real(r8), intent(inout) :: v(IM,KMAX)          ! meridional wind at mid-layer (m/s)
   
   real(r8), intent(in)    :: DELTA               ! physics time step
   real(r8), intent(in)    :: DELTI               ! dynamics time step (model time increment in seconds)
   logical,  intent(in)    :: do_aw, do_awdd
!
! modified arguments
!
   real(r8), intent(inout) :: CBMFX(IM,nctp)      ! cloud base mass flux (kg/m2/s)
!
! output arguments
!
!  updraft, downdraft, and detrainment mass flux (kg/m2/s)
   real(r8), intent(inout), dimension(IJSDIM,KMAX) :: ud_mf, dd_mf, dt_mf
   
   real(r8), intent(out)   :: prec(IJSDIM)        ! precipitation at surface (including snowfall) (kg/m2/s)
   real(r8), intent(out), dimension(ijsdim,kmax) :: qlcn, qicn, w_upi,cnv_mfd, cnv_prc3,&
                                                    cnv_dqldt, clcn, cnv_fice,          &
                                                    cnv_ndrop, cnv_nice, cf_upi

!DDsigma - output added for AW sigma diagnostics
!  interface sigma and vertical velocity by cloud type (1=sfc) 
   real(r8), intent(out), dimension(IM,KMAX,nctp)  :: sigmai, vverti
   real(r8), intent(out), dimension(IM,KMAX)       :: sigma  ! sigma  sigma totaled over cloud type - on interfaces (1=sfc)
!   sigma  terms in eq 91 and 92
   real(r8), dimension(IM,KMAX)                    :: sfluxterm, qvfluxterm, condterm
!DDsigma
!
! output arguments of CS_CUMLUS
!
   real(r8) GTT(IJSDIM,KMAX)           ! temperature tendency [K/s]
   real(r8) GTQ(IJSDIM,KMAX,NTR)       ! tracer tendency [kg/kg/s]
   real(r8) GTU(IJSDIM,KMAX)           ! zonal velocity tendency [m/s2]
   real(r8) GTV(IJSDIM,KMAX)           ! meridional velocity tendency [m/s2]
   real(r8) GTPRP(IJSDIM,KMAX)         ! precipitation (including snowfall) flux at interfaces [kg/m2/s]
   real(r8) GSNWP(IJSDIM,KMAX)         ! snowfall flux at interfaces [kg/m2/s]

!  real(r8) CMDET(IJSDIM,KMAX)         ! detrainment mass flux [kg/m2/s]
!  real(r8) GTLDET(IJSDIM,KMAX)        ! cloud liquid tendency by detrainment [1/s]
!  real(r8) GTIDET(IJSDIM,KMAX)        ! cloud ice tendency by detrainment [1/s]

!DD removed as output arguments
!  real(r8) :: jctop(IJSDIM)           ! o row of top-of-deep-convection indices passed out.
!  real(r8) :: jcbot(IJSDIM)           ! o row of base of cloud indices passed out.

!   The following commented by moorthi to save memory for now - oct 2016
!  real(r8) :: dlf(IJSDIM,KMAX)        ! scattered version of the detraining cld h2o tend (kg/kg/s)
!  real(r8) :: pflx(IJSDIM,KMAX+1)     ! scattered precip flux at each level
!  real(r8) :: cme(IJSDIM,KMAX)        ! condensation - evaporation
!  real(r8) :: rliq(IJSDIM)            ! reserved liquid (not yet in cldliq) for energy integrals (m/s)
!  real(r8) :: flxprec(IJSDIM,KMAX+1)  ! precipitation flux (including snowfall) at interfaces (kg/m2/s)
!  real(r8) :: flxsnow(IJSDIM,KMAX+1)  ! snowfall flux at interfaces (kg/m2/s)

   integer  KT(IJSDIM,nctp)            ! cloud top index for each cloud type

   real(r8) :: cape(IJSDIM)            ! convective available potential energy (J/kg)
   real(r8) :: snow(IJSDIM)            ! snowfall at surface (kg/m2/s)

!
! input arguments of CS_CUMLUS
!
   real(r8) GDT(IJSDIM,KMAX)           ! temperature [K]
   real(r8) GDQ(IJSDIM,KMAX,NTR)       ! tracers including moisture [kg/kg]  !DDsigmadiag
   real(r8) GDU(IJSDIM,KMAX)           ! zonal wind [m/s]
   real(r8) GDV(IJSDIM,KMAX)           ! meridional wind [m/s]
   real(r8) GDTM(IJSDIM,KMAX+1)        ! temperature at boundaries of layers [K]
   real(r8) GDP(IJSDIM,KMAX)           ! pressure [Pa]
   real(r8) GDPM(IJSDIM,KMAX+1)        ! pressure at boundaries of layers [Pa]
   real(r8) GDZ(IJSDIM,KMAX)           ! altitude [m]
   real(r8) GDZM(IJSDIM,KMAX+1)        ! altitude at boundaries of layers [m]
   real(r8) delp(IJSDIM,KMAX)        ! altitude at boundaries of layers [m]
!
! local variables
!
!DD   real(r8) :: zs(IJSDIM)           ! surface height [m]

   integer KTMAX(IJSDIM)               ! max of KT
   real(r8)    :: ftintm, wrk, wrk1, tem
   integer i, k, n, ISTS, IENS, kp1, ipr
!  integer i, k, n, iunit

!DD borrowed from RAS to go form total condensate to ice/water separately
!  parameter (tf=130.16, tcr=160.16, tcrf=1.0/(tcr-tf),tcl=2.0)
!  parameter (tf=230.16, tcr=260.16, tcrf=1.0/(tcr-tf))
   real(r8), parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf), tcl=2.0
   logical, save       :: first=.true.
   logical lprnt


   precz0 = precz0in
   preczh = preczhin
!
!  lprnt = lat == 15 .and. kdt <= 2
   if (first) then
!    write(1000+mype,*)' precz0=',precz0,' preczh=',preczh,' nctp=',nctp
     do i=1,ntr
       IMFXR(i)  = 0
     enddo
!    IMFXR(1)   = 1
!    IMFXR(ITL) = 1
!    IMFXR(ITI) = 1
     first      = .false.
   endif
!
   ISTS = 1
   IENS = IJSDIM

   do k=1,KMAX+1
     do i=1,IJSDIM
       GDZM(i,k) = zi(i,k) * gravi
       GDPM(i,k) = paph(i,k)
     enddo
   enddo

   do k=1,KMAX
     do i=1,IJSDIM
       GDT(i,k)   = t(i,k)
       GDU(i,k)   = u(i,k)
       GDV(i,k)   = v(i,k)
       GDZ(i,k)   = zm(i,k) * gravi
       GDP(i,k)   = pap(i,k)
       GDQ(i,k,1) = q(i,k)
       delp(i,k)  = paph(i,k) - paph(i,k+1)
     enddo
   enddo

!DD following adapted from ras
   if (clw(1,1,2) <= -999.0) then  ! input ice/water are together
     do k=1,kmax
       do i=1,IJSDIM
         tem = clw(i,k,1) * MAX(ZERO, MIN(ONE, (TCR-t(i,k))*TCRF))
         clw(i,k,2) = clw(i,k,1) - tem
         clw(i,k,1) = tem
       enddo
     enddo
   endif
!DD end ras adaptation
   do k=1,kmax
     do i=1,ijsdim
       tem = min(clw(i,k,1), 0.0)
       wrk = min(clw(i,k,2), 0.0)
       clw(i,k,1) = clw(i,k,1) - tem
       clw(i,k,2) = clw(i,k,2) - wrk
       gdq(i,k,1) = gdq(i,k,1) + tem + wrk
     enddo
   enddo

   do n=2,NTR
     do k=1,KMAX
       do i=1,IJSDIM
         GDQ(i,k,n) = clw(i,k,n-1) !DDsigmadiag
       enddo
     enddo
   enddo
!***************************************************************************************
!  iunit = 400 + mype
!  write(iunit,*)kmax,'kmax',delta,'delta',im,'im',ijsdim,'ijsdim',iens,'iens',ists,'ists' !DDdebug
!  write(iunit,*),i  !DDdebug
!  do i = 1, 1                  !DDdebug
!    write(iunit,*)'gdt'        !DDdebug
!    write(iunit,*)gdt(I,:)     !DDdebug
!    write(iunit,*)'gdu'        !DDdebug
!    write(iunit,*)gdu(I,:)     !DDdebug
!    write(iunit,*)'gdv'        !DDdebug
!    write(iunit,*)gdv(I,:)     !DDdebug
!    do k = 1,ntr  !DDdebug
!      write(iunit,*)'gdq',k    !DDdebug
!      write(iunit,*)gdq(I,:,k) !DDdebug
!    enddo  !DDdebug
!    write(iunit,*)'gdz'        !DDdebug
!    write(iunit,*)gdz(I,:)     !DDdebug
!    write(iunit,*)'gdp'        !DDdebug
!    write(iunit,*)gdp(I,:)     !DDdebug
!    write(iunit,*)'gdzm'       !DDdebug
!    write(iunit,*)gdzm(I,:)    !DDdebug
!    write(iunit,*)'gdpm'       !DDdebug
!    write(iunit,*)gdpm(I,:)    !DDdebug
!    write(iunit,*)'cbmfx'      !DDdebug
!    write(iunit,*)cbmfx(I,:)   !DDdebug
!  enddo  !DDdebug
!***************************************************************************************
!
! calculate temperature at interfaces
!
!  call TINTP( IJSDIM, KMAX  ,     & !DD dimensions
!              GDTM,               & ! output
!              GDT, GDP, GDPM,     & ! input
!              ISTS, IENS      )     ! active array size

   DO K=2,KMAX
     DO I=ISTS,IENS
       wrk  = one / GDP(I,K)
       wrk1 = one / LOG(GDP(I,K-1)*wrk)
       FTINTM    = wrk1 * LOG(GDPM(I,K)*wrk)
       GDTM(I,K) = FTINTM*GDT(I,K-1) + (one-FTINTM)*GDT(I,K)
      ENDDO
   ENDDO

   DO I=ISTS,IENS
     GDTM(I,KMAX+1) = GDT(I,KMAX)
     GDTM(I,1)      = GDT(I,1)        ! Is this a good approximation ? - Moorthi
   ENDDO

!DDsigma - initialize the sigma diagnostics
   do n=1,nctp
     do k=1,kmax
       do i=ists,iens
         sigmai(i,k,n) = zero !DDsigma
         vverti(i,k,n) = zero !DDsigma
       enddo
     enddo
   enddo
   do k=1,kmax
     do i=ists,iens
       sigma(i,k)  = zero !DDsigma
     enddo
   enddo
!
! call main routine
!
!***************************************************************************************
   call CS_CUMLUS (im    , IJSDIM, KMAX  , NTR   ,    &  !DD dimensions
                   otspt(1,1), otspt(1,2), lprnt, ipr,&
                   GTT   , GTQ   , GTU   , GTV   ,    & ! output
!                  CMDET , GTLDET, GTIDET,            & ! output
!                  GTPRP , GSNWP , GMFX0 ,            & ! output
!                  GMFX1 , cape  , KT    ,            & ! output
!                  dt_mf , GTLDET, GTIDET,            & ! output
                   dt_mf ,                            & ! output
                   GTPRP , GSNWP , ud_mf ,            & ! output
                   dd_mf , cape  , KT    ,            & ! output
                   CBMFX ,                            & ! modified
                   GDT   , GDQ   , GDU   , GDV   ,    & ! input
                   GDTM  ,                            & ! input
                   GDP   , GDPM  , GDZ   , GDZM  ,    & ! input
                   DELTA , DELTI , ISTS  , IENS, mype,& ! input
                   fscav,  fswtr,  wcbmaxm, nctp,     &
                   sigmai, sigma,  vverti,            & ! input/output !DDsigma
                   sfluxterm, qvfluxterm, do_aw, do_awdd)!DDsigmadiag, output
!
!
!DD detrainment has to be added in for GFS
!
!  if (lprnt) write(0,*)' aft cs_cum gtqi=',gtq(ipr,:,2)
!  if (lprnt) write(0,*)' aft cs_cum gtql=',gtq(ipr,:,3)
   do n=2,NTR
     do k=1,KMAX
       do i=1,IJSDIM
         clw(i,k,n-1) = GDQ(i,k,n) + GTQ(i,k,n) * delta
!        clw(i,k,1)   = GDQ(i,k,2) + (gtq(i,k,2) + gtidet(i,k)) * delta
!        clw(i,k,2)   = GDQ(i,k,3) + (gtq(i,k,3) + gtldet(i,k)) * delta
       enddo
     enddo
   enddo

!  if (ntr > 3) then    ! update tracers
!    do n=4,ntr
!      do k=1,kmax
!        do i=1,ijsdim
!          clw(i,k,n-1) = gdq(i,k,n) + gtq(i,k,n) * delta
!        enddo
!      enddo
!    enddo
!  endif
!
   do k=1,KMAX
     do i=1,IJSDIM
!DD    heat(i,KMAX-k+1) = CP*GTT(i,k) - EMELT*GTIDET(i,k)
!DD    dlf (i,k)     = GTLDET(i,k) + GTIDET(i,k)
!DD    rliq(i)       = (GTLDET(i,k)+GTIDET(i,k))*(GDPM(i,k+1)-GDPM(i,k))/GRAV

       q(i,k)        = GDQ(i,k,1)  + GTQ(i,k,1) * delta
       t(i,k)        = GDT(i,k)    + GTT(i,k)   * delta
       u(i,k)        = GDU(i,k)    + GTU(i,k)   * delta
       v(i,k)        = GDV(i,k)    + GTV(i,k)   * delta
!
!      not used for now - moorthi
!       flxprec(i,k) = GTPRP(i,k)
!       flxsnow(i,k) = GSNWP(i,k)

! Set the mass fluxes.
!       ud_mf  (i,k) = GMFX0(i,k)
!       dd_mf  (i,k) = GMFX1(i,k)
!       dt_mf  (i,k) = CMDET(i,k)
!      if (lprnt .and. i == ipr) write(0,*)' k=',k,'in cs_conv qv=',q(ipr,k)&
!    ,       ' GDQ=',gdq(ipr,k,1),' gtq=',GTQ(ipr,k,1)*delta,' kdt=',kdt
     enddo
   enddo
!  if (lprnt) write(0,*)' in cs_conv qv=',q(ipr,1:35)

   if (ncld == 2) then  ! for 2M microphysics, always output these variables
     if (do_aw) then
       do k=1,KMAX
         kp1 = min(k+1,kmax)
         do i=1,IJSDIM
           qicn(i,k)      = max(0.0, clw(i,k,1)-gdq(i,k,2))
           qlcn(i,k)      = max(0.0, clw(i,k,2)-gdq(i,k,3))

!!         qicn(i,k)      = max(0.0, (gtq(i,k,2)+gtidet(i,k)) * delta)
!!         qlcn(i,k)      = max(0.0, (gtq(i,k,3)+gtldet(i,k)) * delta)
           cnv_fice(i,k)  = qicn(i,k) / max(1.0e-10,qicn(i,k)+qlcn(i,k))
!
           CNV_MFD(i,k)   = dt_mf(i,k) * (1.0/delta)
!!         CNV_DQLDT(i,k) = dt_mf(i,k) * max(0.0,gtidet(i,k)+gtldet(i,k))
           CNV_DQLDT(i,k) = (qicn(i,k)+qlcn(i,k)) / delta
           CNV_PRC3(i,k)  = 0.0
           CNV_NDROP(i,k) = 0.0
           CNV_NICE(i,k)  = 0.0
           cf_upi(i,k)    = max(0.0, min(0.5, 0.5*(sigma(i,k)+sigma(i,kp1))))
           CLCN(i,k)      = cf_upi(i,k)                     !downdraft is below updraft
!!         clcn(i,k)      = max(0.0,min(0.01*log(1.0+500*ud_mf(i,k)/delta),0.25))

           w_upi(i,k)     = 0.0
!!         w_upi(i,k)     = ud_mf(i,k)*(t(i,k)+epsvt*gdq(i,k,1)) * rair &
!!                       / (delta*max(cf_upi(i,k),1.e-12)*gdp(i,k))
         enddo
       enddo
!!     do n=1,nctp
         do k=1,kmax
           do i=1,ijsdim
!!           w_upi(i,k) = w_upi(i,k) + 0.25*(sigmai(i,k,n)+sigmai(i,k+1,n))     &
!!                                   *      (vverti(i,k,n)+vverti(i,k+1,n))
             tem = 0.0
             do n=1,nctp
               tem = tem + sigmai(i,k,n)
               w_upi(i,k) = w_upi(i,k) + sigmai(i,k,n) * vverti(i,k,n)
             enddo
             w_upi(i,k) = w_upi(i,k) / max (1.0e-10,tem)

!!         cf_upi(i,k)    = max(0.0,min(0.01*log(1.0+500*ud_mf(i,k)/delta),0.25))
!!   &                                               500*ud_mf(i,k)/delta),0.60))
!!         CLCN(i,k)      = cf_upi(i,k)                     !downdraft is below updraft

!!         w_upi(i,k)     = ud_mf(i,k)*(t(i,k)+epsvt*gdq(i,k,1)) * rair &
!!                       / (delta*max(cf_upi(i,k),1.e-12)*gdp(i,k))

           enddo
         enddo
!!     enddo
!!     do k=1,kmax
!!       do i=1,ijsdim
!!         w_upi(i,k) = w_upi(i,k) / max(1.0e-9, 0.5*(sigma(i,k)+sigma(i,k+1)))
!!       enddo
!!     enddo
     else
       do k=1,KMAX
         do i=1,IJSDIM
           qicn(i,k)      = max(0.0, clw(i,k,1)-gdq(i,k,2))
           qlcn(i,k)      = max(0.0, clw(i,k,2)-gdq(i,k,3))
!          qicn(i,k)      = max(0.0, (gtq(i,k,2)+gtidet(i,k)) * delta)
!          qlcn(i,k)      = max(0.0, (gtq(i,k,3)+gtldet(i,k)) * delta)
           cnv_fice(i,k)  = qicn(i,k) / max(1.0e-10,qicn(i,k)+qlcn(i,k))
! 
           CNV_MFD(i,k)   = dt_mf(i,k) * (1/delta)
!          CNV_DQLDT(i,k) = max(0.0,gtidet(i,k)+gtldet(i,k))
           CNV_DQLDT(i,k) = (qicn(i,k)+qlcn(i,k)) / delta
           CNV_PRC3(i,k)  = 0.0
           CNV_NDROP(i,k) = 0.0
           CNV_NICE(i,k)  = 0.0
           cf_upi(i,k)    = max(0.0,min(0.01*log(1.0+500*ud_mf(i,k)/delta),0.25))
!    &                                               500*ud_mf(i,k)/delta),0.60))
           CLCN(i,k)      = cf_upi(i,k)                     !downdraft is below updraft
           
           w_upi(i,k)     = ud_mf(i,k)*(t(i,k)+epsvt*gdq(i,k,1)) * rair &
                          / (delta*max(cf_upi(i,k),1.e-12)*gdp(i,k))
         enddo
       enddo
     endif
   endif

!****************************************************************************
!  do i=1,1                     !DDdebug 
!    write(iunit,*)'gtt'        !DDdebug
!    write(iunit,*)gtt(I,:)     !DDdebug
!    do k = 1,ntr               !DDdebug
!      write(iunit,*)'gtq',k    !DDdebug
!      write(iunit,*)gtq(I,:,k) !DDdebug
!    enddo                      !DDdebug
!    write(iunit,*)'gtu'        !DDdebug
!    write(iunit,*)gtu(I,:)     !DDdebug
!    write(iunit,*)'gtv'        !DDdebug
!    write(iunit,*)gtv(I,:)     !DDdebug
!    write(iunit,*)'gtprp'      !DDdebug
!    write(iunit,*)gtprp(I,:)   !DDdebug
!    write(iunit,*)'gsnwp'      !DDdebug
!    write(iunit,*)gsnwp(I,:)   !DDdebug
!    write(iunit,*)'gmfx0'      !DDdebug
!    write(iunit,*)gmfx0(I,:)   !DDdebug
!    write(iunit,*)'gmfx1'      !DDdebug
!    write(iunit,*)gmfx1(I,:)   !DDdebug
!    write(iunit,*)'cmdet'      !DDdebug
!    write(iunit,*)cmdet(I,:)   !DDdebug
!    write(iunit,*)'cbmfx'      !DDdebug
!    write(iunit,*)cbmfx(I,:)   !DDdebug
!    write(iunit,*)'kt'         !DDdebug
!    write(iunit,*)kt(I,:)      !DDdebug
!    write(iunit,*)'cape'       !DDdebug
!    write(iunit,*)cape(I)      !DDdebug
!    write(iunit,*)'gtldet'     !DDdebug
!    write(iunit,*)gtldet(I,:)  !DDdebug
!    write(iunit,*)'gtidet'     !DDdebug
!    write(iunit,*)gtldet(I,:)  !DDdebug
!  enddo                        !DDdebug
!****************************************************************************
!
   KTMAX = 1
   do n=1,nctp
     do i=1,IJSDIM
        KTMAX(i) = max(KTMAX(i), KT(i,n))
     enddo
   enddo
!
   do i=1,IJSDIM
!    jctop(i) = KTMAX(i)
     prec(i)  = GTPRP(i,1)
     snow(i)  = GSNWP(i,1)
!    rliq(i)  = rliq(i)/1000._r8      ! kg/m2/s => m/s
   enddo
!  if (lprnt) then
!    write(0,*)' aft cs_cum prec=',prec(ipr),'GTPRP=',GTPRP(ipr,1)
!  endif
  
!  cme   = zero    ! temporarily set to be zero
!  pflx  = zero    ! temporarily set to be zero
!  jcbot = 1       ! set to be the lowest layer
!    if (lprnt) then
!    write(2000+mype,*)' gdq=',gdq(13,:,1)
!    write(2000+mype,*)' q=',q(13,:)
!    endif

!    if (do_aw) then
!    call moist_bud(ijsdim,ijsdim,im,kmax,mype,kdt,grav,delta,delp,prec &
!    ,              gdq(1,1,1), gdq(1,1,2), gdq(1,1,3)                  &
!    ,              q,clw(1,1,1),clw(1,1,2),'cs_conv_aw')
!    endif

   end subroutine cs_convr


!************************************************************************
!* Original source code in MIROC5
!*
!* PACKAGE PCUMC  !!  physics: cumulus parameterization with
!*                             state-dependent entrainment rate
!*                             developed by Minoru Chikira
!* [Note]
!* -This routine works as the prognostic Arakawa-Schubert scheme
!*  if OPT_ASMODE is specified.
!* -Specify OPT_NS02 to use entrainment rate of Neggers et al. (2002)
!* -Specify OPT_CUMBGT to check water and energy budget.
!* -Specify OPT_CUMCHK to check range of output values.
!*
!*   [HIS] 08/09/19(chikira)   MIROC4.1
!*         08/10/30(hiro)      CMT modified
!*         08/11/11(chikira)   Neggers et al. (2002)
!*         08/12/3 (chikira)   downdraft detrainment modified
!*         08/12/3 (chikira)   COSP output
!*         09/02/24(chikira)   fix convective inhibition
!*         09/04/16(hiro)      CMIP5 output (cbasep,ctopp)
!*         09/09/03(yokohata)  COSP
!*         10/11/19(toshi)     small bug fix
!*         14/02/07(chikira)   CUMDWN bug fix, CMT modified
!************************************************************************
! cumulus main routine
! --------------------
   SUBROUTINE CS_CUMLUS (im    , IJSDIM, KMAX  , NTR   ,   & !DD dimensions
                         otspt1, otspt2, lprnt, ipr,       &
                         GTT   , GTQ   , GTU   , GTV   ,   & ! output
!                        CMDET , GTLDET, GTIDET,           & ! output
                         CMDET ,                           & ! output
                         GTPRP , GSNWP , GMFX0 ,           & ! output
                         GMFX1 , CAPE  , KT    ,           & ! output
!                        CUMCLW, CUMFRC,
                         CBMFX ,                            & ! modified
                         GDT   , GDQ   , GDU   , GDV   ,    & ! input
                         GDTM  ,                            & ! input
                         GDP   , GDPM  , GDZ   , GDZM  ,    & ! input
!                        GDCFRC,
                         DELTA , DELTI , ISTS  , IENS, mype,& ! input
                         fscav,  fswtr,  wcbmaxm, nctp,     & !
                         sigmai, sigma,  vverti,            & ! input/output !DDsigma
                         sfluxterm, qvfluxterm, do_aw, do_awdd ) ! output !DDsigmadiag
!
   IMPLICIT NONE
      
   INTEGER, INTENT(IN)   :: im, IJSDIM, KMAX, NTR, mype, nctp, ipr !! DD, for GFS, pass in
   logical, intent(in)   :: do_aw, do_awdd  ! switch to apply Arakawa-Wu to the tendencies
   logical, intent(in)   :: otspt1(ntr), otspt2(ntr), lprnt
!
! [OUTPUT]
   REAL(r8), INTENT(OUT) :: GTT   (IJSDIM, KMAX     ) !! heating rate
   REAL(r8), INTENT(OUT) :: GTQ   (IJSDIM, KMAX, NTR) !! change in q
   REAL(r8), INTENT(OUT) :: GTU   (IJSDIM, KMAX     ) !! tendency of u
   REAL(r8), INTENT(OUT) :: GTV   (IJSDIM, KMAX     ) !! tendency of v
   REAL(r8), INTENT(OUT) :: CMDET (IJSDIM, KMAX     ) !! detrainment mass flux

!  REAL(r8), INTENT(OUT) :: GTLDET(IJSDIM, KMAX     ) !! cloud liquid tendency by detrainment
!  REAL(r8), INTENT(OUT) :: GTIDET(IJSDIM, KMAX     ) !! cloud ice tendency by detrainment

! assuming there is no flux  at the top of the atmospherea - Moorthi
   REAL(r8), INTENT(OUT) :: GTPRP (IJSDIM, KMAX     ) !! rain+snow flux
   REAL(r8), INTENT(OUT) :: GSNWP (IJSDIM, KMAX     ) !! snowfall flux
   REAL(r8), INTENT(OUT) :: GMFX0 (IJSDIM, KMAX     ) !! updraft mass flux
   REAL(r8), INTENT(OUT) :: GMFX1 (IJSDIM, KMAX     ) !! downdraft mass flux

   REAL(r8), INTENT(OUT) :: CAPE  (IJSDIM           )
   INTEGER , INTENT(OUT) :: KT    (IJSDIM, NCTP     ) !! cloud top
!
!  [MODIFIED]
   REAL(r8), INTENT(INOUT) :: CBMFX (IM, NCTP)        !! cloud base mass flux

!DDsigma - output added for AW sigma diagnostics
! sigma and vert. velocity as a function of cloud type (1==sfc)
   real(r8), intent(out), dimension(IM,KMAX,nctp)   :: sigmai, vverti
   real(r8), intent(out), dimension(IM,KMAX)        :: sigma       !DDsigma sigma totaled over cloud type - on interfaces (1=sfc)

! for computing AW flux form of tendencies
! The tendencies are summed over all cloud types
   real(r8), intent(out), dimension(IM,KMAX) ::    &  !DDsigmadiag
                                    sfluxterm, qvfluxterm    ! tendencies of DSE and water vapor due to eddy mass flux
   real(r8), dimension(IM,KMAX) ::  qlfluxterm, qifluxterm   ! tendencies of cloud water and cloud ice due to eddy mass flux

! The fluxes are for an individual cloud type and reused.
!  condtermt, condtermq  are eddy flux of temperature and water vapor
   real(r8), dimension(IM,KMAX) :: condtermt, condtermq, frzterm, &
                                   prectermq, prectermfrz
!
!  [INPUT]
   REAL(r8), INTENT(IN) :: GDT   (IJSDIM, KMAX     ) ! temperature T
   REAL(r8), INTENT(IN) :: GDQ   (IJSDIM, KMAX, NTR) ! humidity, tracer  !DDsigmadiag
   REAL(r8), INTENT(IN) :: GDU   (IJSDIM, KMAX     ) ! westerly u
   REAL(r8), INTENT(IN) :: GDV   (IJSDIM, KMAX     ) ! southern wind v
   REAL(r8), INTENT(IN) :: GDTM  (IJSDIM, KMAX+1   ) ! temperature T
   REAL(r8), INTENT(IN) :: GDP   (IJSDIM, KMAX     ) ! pressure P
   REAL(r8), INTENT(IN) :: GDPM  (IJSDIM, KMAX+1   ) ! pressure (half lev)
   REAL(r8), INTENT(IN) :: GDZ   (IJSDIM, KMAX     ) ! altitude
   REAL(r8), INTENT(IN) :: GDZM  (IJSDIM, KMAX+1   ) ! altitude
   REAL(r8), INTENT(IN) :: DELTA                     ! delta(t) (dynamics)
   REAL(r8), INTENT(IN) :: DELTI                     ! delta(t) (internal variable)
   INTEGER, INTENT(IN)  :: ISTS, IENS                ! array range

   real(r8), intent(in) :: fscav(ntr), fswtr(ntr), wcbmaxm(ijsdim)
!
!  [INTERNAL WORK]
   REAL(r8)     GPRCC (IJSDIM, NTR)       ! rainfall
   REAL(r8)     GSNWC (IJSDIM)            ! snowfall
   REAL(r8)     CUMCLW(IJSDIM, KMAX)      ! cloud water in cumulus
   REAL(r8)     CUMFRC(IJSDIM)            ! cumulus cloud fraction
!COSP
!  REAL(r8)     QLIQC (IJSDIM, KMAX)      ! cumulus cloud liquid water [kg/kg]
!  REAL(r8)     QICEC (IJSDIM, KMAX)      ! cumulus cloud ice [kg/kg]
!  REAL(r8)     GPRCPF(IJSDIM, KMAX)      ! rainfall flux at full level
!  REAL(r8)     GSNWPF(IJSDIM, KMAX)      ! snowfall flux at full level
!
   REAL(r8)     GTCFRC(IJSDIM, KMAX)      ! change in cloud fraction
   REAL(r8)     FLIQC (IJSDIM, KMAX)      ! liquid ratio in cumulus
!
   REAL(r8)     GDCFRC(IJSDIM, KMAX)      ! cloud fraction
!
!  REAL(r8)     GDQI  (IJSDIM, KMAX)      ! cloud ice
!  REAL(r8)     GTQI  (IJSDIM, KMAX)      ! tendency of cloud ice
!  REAL(r8)     GTQL  (IJSDIM, KMAX)      ! tendency of cloud liquid
!
   REAL(r8)     GDW   (IJSDIM, KMAX)      ! total water
   REAL(r8)     DELP  (IJSDIM, KMAX)
   REAL(r8)     GDQS  (IJSDIM, KMAX)      ! saturate moisture
   REAL(r8)     FDQS  (IJSDIM, KMAX)
   REAL(r8)     GAM   (IJSDIM, KMAX)
   REAL(r8)     GDS   (IJSDIM, KMAX)      ! dry static energy
   REAL(r8)     GDH   (IJSDIM, KMAX)      ! moist static energy
   REAL(r8)     GDHS  (IJSDIM, KMAX)      ! saturate MSE
!
   REAL(r8)     GCYM  (IJSDIM, KMAX)      ! norm. mass flux (half lev)
   REAL(r8)     GCHB  (IJSDIM)            ! cloud base MSE-Li*Qi
   REAL(r8)     GCWB  (IJSDIM)            ! cloud base total water
   REAL(r8)     GCUB  (IJSDIM)            ! cloud base U
   REAL(r8)     GCVB  (IJSDIM)            ! cloud base V
   REAL(r8)     GCIB  (IJSDIM)            ! cloud base ice
   REAL(r8)     ELAM  (IJSDIM, KMAX, NCTP)! entrainment (rate*massflux)
   REAL(r8)     GCYT  (IJSDIM, NCTP)      ! norm. mass flux @top
   REAL(r8)     GCHT  (IJSDIM, NCTP)      ! cloud top MSE
   REAL(r8)     GCQT  (IJSDIM, NCTP)      ! cloud top q
   REAL(r8)     GCwT  (IJSDIM)            ! cloud top total water
   REAL(r8)     GCUT  (IJSDIM, NCTP)      ! cloud top U
   REAL(r8)     GCVT  (IJSDIM, NCTP)      ! cloud top V
   REAL(r8)     GCLT  (IJSDIM, NCTP)      ! cloud top cloud water
   REAL(r8)     GCIT  (IJSDIM, NCTP)      ! cloud top cloud ice
   REAL(r8)     GTPRT (IJSDIM, NCTP)      ! precipitation/M
   REAL(r8)     GCLZ  (IJSDIM, KMAX)      ! cloud liquid for each CTP
   REAL(r8)     GCIZ  (IJSDIM, KMAX)      ! cloud ice for each CTP

   REAL(r8)     ACWF  (IJSDIM, NCTP)      ! cloud work function
   REAL(r8)     GPRCIZ(IJSDIM, KMAX)      ! precipitation
   REAL(r8)     GSNWIZ(IJSDIM, KMAX)      ! snowfall
   REAL(r8)     GTPRC0(IJSDIM)            ! precip. before evap.

   REAL(r8)     GMFLX (IJSDIM, KMAX)      ! mass flux (updraft+downdraft)
   REAL(r8)     QLIQ  (IJSDIM, KMAX)      ! total cloud liquid
   REAL(r8)     QICE  (IJSDIM, KMAX)      ! total cloud ice
   REAL(r8)     GPRCI (IJSDIM, KMAX)      ! rainfall generation
   REAL(r8)     GSNWI (IJSDIM, KMAX)      ! snowfall generation

   REAL(r8)     GPRCP (IJSDIM, KMAX)      ! rainfall flux
!
   REAL(r8)     GTEVP (IJSDIM, KMAX)      ! evaporation+sublimation
   REAL(r8)     GMDD  (IJSDIM, KMAX)      ! downdraft mass flux

!  REAL(r8)     CUMHGT(IJSDIM, NCTP)      ! cloud top height
!  REAL(r8)     CTOPP (IJSDIM)            ! cloud top pressure

   REAL(r8)     GDZTR (IJSDIM)            ! tropopause height
   REAL(r8)     FLIQOU(IJSDIM, KMAX)      ! liquid ratio in cumulus
   INTEGER      KB    (IJSDIM)
   INTEGER      KSTRT (IJSDIM)            ! tropopause level
   REAL(r8)     GAMX
   REAL(r8)     CIN   (IJSDIM)
   INTEGER      JBUOY (IJSDIM)
   REAL(r8)     DELZ, BUOY, DELWC, DELER
   REAL(r8)     WCBX (IJSDIM)
!  REAL(r8)     ERMR  (NCTP)              ! entrainment rate (ASMODE)
!  SAVE         ERMR
   INTEGER      KTMX  (NCTP)              ! max of cloud top
   INTEGER      KTMXT                     ! max of cloud top
   REAL(r8)     TIMED
   REAL(r8)     GDCLDX, GDMU2X, GDMU3X
!
   LOGICAL      OOUT1, OOUT2

   REAL(r8)     HBGT (IJSDIM)             ! imbalance in column heat
   REAL(r8)     WBGT (IJSDIM)             ! imbalance in column water
   
!DDsigma begin local work variables - all on model interfaces (sfc=1)
   REAL(r8)     lamdai                    ! lamda for cloud type ctp
   REAL(r8)     gdqm,      gdlm, gdim     !  water vaper
   REAL(r8)     gdtrm                     !  water vaper tracer
   character(len=4) :: cproc  !DDsigmadiag

! the following are new arguments to cumup to get them out for AW
   REAL(r8)   wcv   (IJSDIM, KMAX)        ! in-cloud vertical velocity
   REAL(r8)   GCTM  (IJSDIM, KMAX)        ! cloud T (half lev)   !DDsigmadiag make output
   REAL(r8)   GCQM  (IJSDIM, KMAX)        ! cloud q (half lev)   !DDsigmadiag make output
   REAL(r8)   GCwM  (IJSDIM, KMAX)        ! cloud q (half lev)   !DDsigmadiag make output
   REAL(r8)   GCiM  (IJSDIM, KMAX)        ! cloud q (half lev)   !DDsigmadiag make output
   REAL(r8)   GClM  (IJSDIM, KMAX)        ! cloud q (half lev)   !DDsigmadiag make output
   REAL(r8)   GChM  (IJSDIM, KMAX)        ! cloud q (half lev)   !DDsigmadiag make output

! eddy flux profiles for dse, water vapor, cloud water, cloud ice
   REAL(r8), dimension(Kmax+1)        :: sfluxtem, qvfluxtem, qlfluxtem, qifluxtem

! tendency profiles - condensation heating, condensation moistening, heating due to
!                     freezing, total precip production, frozen precip production
   REAL(r8), dimension(ijsdim,Kmax)   :: dtcondtem, dqcondtem, dtfrztem, dqprectem,& ! Moorthi
                                         dfrzprectem, lamdaprod ! product of (1+lamda) through cloud type ctp
   REAL(r8), dimension(ijsdim,Kmax)   :: dtevap, dqevap, dtmelt, dtsubl

! factor to modify precip rate to force conservation of water. With bug fixes it's
!    not doing anything now.
   REAL(r8), dimension(ijsdim)        :: moistening_aw
   real(r8), dimension(ijsdim,kmax)   :: gctbl, gcqbl,gcwbl, gcqlbl, gcqibl, & !DDsigmadiag updraft profiles below cloud Base
                                         sigmad           ! downdraft area fraction
! rhs_q, rhs_h are residuals of condensed water, MSE budgets to compute condensation,
!                   and heating due to freezing
   real(r8)                           :: rhs_q, rhs_h, fsigma, delpinv
!  real(r8)                           :: rhs_q, rhs_h, sftem, qftem, qlftem, qiftem, &
!                                        fsigma  ! factor to reduce mass flux terms (1-sigma**2) for AW
!DDsigma end local work variables
!
! profiles of heating due to precip evaporation, melting and sublimation, and the
!     evap, melting and sublimation rates.

   REAL(r8)    dtdwn  (IJSDIM, KMAX) ! t  tendency downdraft detrainment
   REAL(r8)    dqvdwn (IJSDIM, KMAX) ! qv tendency downdraft detrainment
   REAL(r8)    dqldwn (IJSDIM, KMAX) ! ql tendency downdraft detrainment
   REAL(r8)    dqidwn (IJSDIM, KMAX) ! qi tendency downdraft detrainment

!DDsigma end local work variables
!
!  [INTERNAL PARM]
   REAL(r8), parameter :: WCBMIN = zero       ! min. of updraft velocity at cloud base
!M REAL(r8) :: WCBMAX = 1.4_r8      ! max. of updraft velocity at cloud base
!M wcbas commented by Moorthi since it is not used
!M REAL(r8) :: WCBAS  = 2._r8       ! updraft velocity**2 at cloud base (ASMODE)
!M REAL(r8) :: ERAMIN = 1.e-5_r8    ! min. of entrainment rate
                                    ! used only in OPT_ASMODE
!M REAL(r8) :: ERAMAX = 2.e-3_r8    ! max. of entrainment rate
                                    ! used only in OPT_ASMODE
   LOGICAL  :: OINICB = .false.     ! set 0.d0 to CBMFX when .true.

!  REAL(r8) :: VARMIN = 1.e-13_r8   ! minimum of PDF variance
!  REAL(r8) :: VARMAX = 5.e-7_r8    ! maximum of PDF variance
!  REAL(r8) :: SKWMAX = 0.566_r8    ! maximum of PDF skewness

   REAL(r8) :: PSTRMX = 400.e2_r8   ! max P of tropopause
   REAL(r8) :: PSTRMN = 50.e2_r8    ! min P of tropopause
   REAL(r8) :: GCRSTR = 1.e-4_r8    ! crit. dT/dz tropopause

   real(kind=r8)             :: tem, esat, mflx_e, cbmfl, tem1, tem2, tem3
   INTEGER                   :: KBMX, I, K, CTP, ierr, n, kp1, km1, kk, kbi
!
   LOGICAL, SAVE :: OFIRST = .TRUE.   ! called first time?
!
!  [ONCE]
   IF (OFIRST) THEN

     OFIRST = .FALSE.

!    fscav = 0._r8
!    fswtr = 0._r8
!    write(0,*)' NTR in cs_conv=',ntr,' mype=',mype
!    do n=1,ntr
!      FSCAV(n) = 0._r8       !DD    split declaration and initialization
!      FSWTR(n) = 0._r8       !DD    split declaration and initialization
!    enddo

     IF (OINICB) THEN
       CBMFX = zero
     ENDIF
   ENDIF                  ! ofirst if
!
   do n=1,ntr
     do k=1,kmax
       do i=1,ijsdim
         gtq(i,k,n) = zero
       enddo
     enddo
   enddo

   do k=1,kmax
     do i=1,ijsdim
       gtt(i,k)         = zero
       gtu(i,k)         = zero
       gtv(i,k)         = zero
!      gtqi(i,k)        = zero
!      gtql(i,k)        = zero
       gmflx(i,k)       = zero
       gmfx0(i,k)       = zero
       gprci(i,k)       = zero
       gsnwi(i,k)       = zero
       qliq(i,k)        = zero
       qice(i,k)        = zero
       gtcfrc(i,k)      = zero
       cumclw(i,k)      = zero
       fliqc(i,k)       = zero
       fliqou(i,k)      = zero
!      gprcpf(i,k)      = zero
!      gsnwpf(i,k)      = zero
       sfluxterm(i,k)   = zero
       qvfluxterm(i,k)  = zero
       qlfluxterm(i,k)  = zero
       qifluxterm(i,k)  = zero
       condtermt(i,k)   = zero
       condtermq(i,k)   = zero
       frzterm(i,k)     = zero
       prectermq(i,k)   = zero
       prectermfrz(i,k) = zero
       dtdwn(i,k)       = zero
       dqvdwn(i,k)      = zero
       dqidwn(i,k)      = zero
       dqvdwn(i,k)      = zero
     enddo
   enddo
   do i=1,ijsdim
     gprcc(i,:) = zero
     gtprc0(i)  = zero
     hbgt(i)    = zero
     wbgt(i)    = zero
     gdztr(i)   = zero
     kstrt(i)   = kmax
   enddo

   do k=1,kmax
     do i=1,ijsdim
!      GDQI(i,k) = GDQ(i,k,ITI)
       GDW(i,k)  = GDQ(i,k,1) + GDQ(i,k,ITL) + GDQ(i,k,iti)
     enddo
   enddo
!
   DO K=1,KMAX
     DO I=ISTS,IENS
       DELP(I,K) = GDPM(I,K) - GDPM(I,K+1)
       esat      = min(gdp(i,k), fpvs(gdt(i,k)))
       GDQS(I,K) = min(EPSV*esat/max(gdp(i,k)+epsvm1*esat, 1.0e-10), 1.0)
!      FDQS(I,K) = FDQSAT(GDT(I,K), GDQS(I,K))
       tem       = one / GDT(I,K)
       FDQS(I,K) = GDQS(I,K) * tem * (fact1 + fact2*tem)
       GAM (I,K) = ELOCP*FDQS(I,K)
       GDS (I,K) = CP*GDT(I,K) + GRAV*GDZ(I,K) ! layer dry static energy
       GDH (I,K) = GDS(I,K) + EL*GDQ(I,K,1)    ! layer moist static energy
       GDHS(I,K) = GDS(I,K) + EL*GDQS(I,K)     ! layer sat. moist static energy
     ENDDO
   ENDDO
!
!        < tropopause >
!
   DO K=1,KMAX
     kp1 = k + 1
     DO I=ISTS,IENS
       GAMX = (GDTM(I,KP1)-GDTM(I,K)) / (GDZM(I,KP1)-GDZM(I,K))
       IF ((GDP(I,K) < PSTRMX .AND. GAMX > GCRSTR) .OR. GDP(I,K) < PSTRMN) THEN
          KSTRT(I) = MIN(K, KSTRT(I))
       ENDIF
     ENDDO
   ENDDO
   DO I=ISTS,IENS
     K = KSTRT(I)
     GDZTR(I) = GDZM(I,K)
   ENDDO
!
!DDsigma - arguments added to get subcloud profiles in updraft
!          so AW eddy flux tendencies can be computed

!! Cloud Base properties
   CALL CUMBAS(IJSDIM, KMAX  ,                           & !DD dimensions
               KB    , GCYM  , KBMX  ,                   & ! output
               GCHB  , GCWB  , GCUB  , GCVB  ,           & ! output
               GCIB  ,                                   & ! output
               GDH   , GDW   , GDHS  , GDQS  ,           & ! input
               GDQ(1,1,iti)  , GDU   , GDV   , GDZM  ,   & ! input
               GDPM  , FDQS  , GAM   ,                   & ! input
               ISTS  , IENS                  ,           & !)   ! input
               gctbl, gcqbl,gdq(1,1,1),gcwbl, gcqlbl, gcqibl) ! sub cloud tendencies
!

!DDsigma some initialization  before summing over cloud type
   do k=1,kmax    ! Moorthi
     do i=1,ijsdim
       lamdaprod(i,k)   = one
       dqcondtem(i,k)   = zero
       dqprectem(i,k)   = zero
       dfrzprectem(i,k) = zero
       dtfrztem(i,k)    = zero
       dtcondtem(i,k)   = zero
     enddo
   enddo

   DO CTP=1,NCTP                   ! loop over cloud types

     tem = ctp / DBLE(NCTP)
     do i=1,ijsdim
       DELWC   = tem *  (WCBMAXm(i) - WCBMIN)
       WCBX(I) = DELWC * DELWC
     enddo

! getting more incloud profiles of variables to compute eddy flux tendencies
!    and condensation rates

!! CUMUP computes In-cloud Properties

     CALL CUMUP(IJSDIM, KMAX, NTR   ,                               & !DD dimensions
                ACWF(1,CTP) , ELAM(1,1,CTP),                        & ! output
                GCLZ        , GCIZ        , GPRCIZ      , GSNWIZ,   & ! output
                GCYT(1,CTP) , GCHT(1,CTP) , GCQT (1,CTP),           & ! output
                GCLT(1,CTP) , GCIT(1,CTP) , GTPRT(1,CTP),           & ! output
                GCUT(1,CTP) , GCVT(1,CTP) ,                         & ! output
                KT  (1,CTP) , KTMX(CTP)   ,                         & ! output
                GCYM  ,                                             & ! modified
                wcv   ,                                             & ! !DD-sigma new output
                GCHB  , GCWB  , GCUB  , GCVB  ,                     & ! input  !DDsigmadiag
                GCIB  ,                                             & ! input
                GDU   , GDV   , GDH   , GDW   ,                     & ! input
                GDHS  , GDQS  , GDT   , GDTM  ,                     & ! input
                GDQ   , GDQ(1,1,iti)  , GDZ   , GDZM  ,             & ! input
                GDPM  , FDQS  , GAM   , GDZTR ,                     & ! input
                CPRES , WCBX  ,                                     & ! input
!               CPRES , WCBX  , ERMR(CTP),                          & ! input
                KB    , CTP   , ISTS  , IENS  ,                     & ! input
                gctm, gcqm, gcwm, gchm, gcwt, gclm, gcim,           & ! additional incloud profiles and cloud top total water
                cbmfx(1,ctp), dtcondtem, dqcondtem, dtfrztem  ) !DDsigmadiag
!
!! CUMBMX computes Cloud Base Mass Flux

     CALL CUMBMX(IJSDIM, KMAX  ,                                     & !DD dimensions
                 CBMFX(1,CTP),                                       & ! modified
                 ACWF (1,CTP), GCYT(1,CTP), GDZM     ,               & ! input
                 GDW         , GDQS       , DELP     ,               & ! input
                 KT   (1,CTP), KTMX(CTP)  , KB       ,               & ! input
                 DELTI       , ISTS       , IENS       )
                 
!DDsigma -  begin sigma computation
! At this point cbmfx is updated and we have everything we need to compute sigma

     if (do_aw) then
       do i=ISTS,IENS
         do k=1,kmax+1     ! initialize eddy fluxes for this cloud time
           sfluxtem(k)  = zero
           qvfluxtem(k) = zero
           qlfluxtem(k) = zero
           qifluxtem(k) = zero
         enddo

         cbmfl = cbmfx(i,ctp)
         kk    = kt(i,ctp)      ! cloud top index

         if(cbmfl > zero) then  ! this should avoid zero wcv in the denominator
           kbi = kb(i)          ! cloud top index
           do k=2,kbi           ! compute eddy fluxes below cloud base
             tem = - gcym(i,k) * cbmfl

! first get environment variables at layer interface
             GDQM  = half * (GDQ(I,K,1) + GDQ(I,K-1,1))
             GDlM  = half * (GDQ(I,K,3) + GDQ(I,K-1,3))
             GDiM  = half * (GDQ(I,K,2) + GDQ(I,K-1,2))
!            GDwM  = half * (GDw(I,K)   + GDw(I,K-1))

! flux = mass flux * (updraft variable minus environment variable)
!centered differences
!            sfluxtem(k)  = tem * (gdtm(i,k)-gctbl(i,k))
!            qvfluxtem(k) = tem * (gdqm-gcqbl(i,k))
!            qlfluxtem(k) = tem * (gdlm-gcqlbl(i,k))
!            qifluxtem(k) = tem * (gdim-gcqibl(i,k))

!upstream - This better matches what the original CS tendencies do
             sfluxtem(k)  = tem * (gdt(i,k)+gocp*(gdz(i,k)-gdzm(i,k))-gctbl(i,k))
             qvfluxtem(k) = tem * (gdq(i,k,1)-gcqbl(i,k))
             qlfluxtem(k) = tem * (gdq(i,k,3)-gcqlbl(i,k))
             qifluxtem(k) = tem * (gdq(i,k,2)-gcqibl(i,k))

           enddo
           do k=kbi,kk             ! loop from cloud base to cloud top
             km1 = k - 1
             rhs_h = zero
             rhs_q = zero
! get environment variables interpolated to layer interface
             GDQM   = half * (GDQ(I,K,1) + GDQ(I,KM1,1))  ! as computed in cumup
!            GDwM   = half * (GDw(I,K)   + GDw(I,KM1 ))
             GDlM   = half * (GDQ(I,K,3) + GDQ(I,KM1,3))
             GDiM   = half * (GDQ(I,K,2) + GDQ(I,KM1,2))
             mflx_e = gcym(i,k) * cbmfl          ! mass flux at level k for cloud ctp

! this is the computation of lamda for a cloud type, and then updraft area fraction
! (sigmai for a single cloud type)
!            gdtvm  = gdtm(i,k) * (1 + epsvt * gdqm)
!            gdrhom = gdpm(i,k) / (rair * gdtvm)        ! gas law
!            gdrhom = gdpm(i,k) / (rair * gdtm(i,k)*(one+epsvt*gdqm))        ! gas law
!            lamdai = mflx_e / (gdrhom*wcv(i,k))

             lamdai = mflx_e * rair * gdtm(i,k)*(one+epsvt*gdqm)         &
                    / (gdpm(i,k)*wcv(i,k))
             lamdaprod(i,k)  = lamdaprod(i,k) * (one+lamdai)
             vverti(i,k,ctp) = wcv(i,k)
             sigmai(i,k,ctp) = lamdai / lamdaprod(i,k)
             sigma(i,k) = sigma(i,k) + sigmai(i,k,ctp)
!            fsigma     = 1.0   ! no aw effect, comment following lines to undo AW
!            fsigma     = (one - sigmai(i,k,ctp)*sigmai(i,k,ctp))
             fsigma     = one - sigma(i,k)
!            fsigma     = (one - sigmai(i,k,ctp)) * (one - sigmai(i,k,ctp))

! compute tendencies based on mass flux, and tendencies based on condensation
! fsigma is the AW reduction of flux tendencies

             if(k > kbi) then   ! uncomment for test
! flux = mass flux * (updraft variable minus environment variable)

               tem = - fsigma * mflx_e
!centered
!              sfluxtem(k)  = tem * (gdtm(i,k)+gocp*gdzm(i,k)-gctm(i,k))
!              qvfluxtem(k) = tem * (gdqm-gcqm(i,k))
!              qlfluxtem(k) = tem * (gdlm-gclm(i,k))
!              qifluxtem(k) = tem * (gdim-gcim(i,k))

!upstream  - This better matches what the original CS tendencies do
               if(k < kk) then
                 sfluxtem(k)  = tem * (gdt(i,k)+gocp*gdz(i,k)-gctm(i,k))
                 qvfluxtem(k) = tem * (gdq(i,k,1)-gcqm(i,k))
                 qlfluxtem(k) = tem * (gdq(i,k,3)-gclm(i,k))
                 qifluxtem(k) = tem * (gdq(i,k,2)-gcim(i,k))
               else
! centered at top of cloud
                 sfluxtem(k)  = tem * (gdtm(i,k)+gocp*gdzm(i,k)-gctm(i,k))
                 qvfluxtem(k) = tem * (gdqm-gcqm(i,k))
                 qlfluxtem(k) = tem * (gdlm-gclm(i,k))
                 qifluxtem(k) = tem * (gdim-gcim(i,k))
               endif

! the condensation terms - these come from the MSE and condensed water budgets for
!   an entraining updraft
!              if(k > kb(i)) then  ! comment for test
!              if(k <= kk) then             ! Moorthi
!              if(k < kt(i,ctp)) then
!                rhs_h = cbmfl*(gcym(i,k)*gchm(i,k) - (gcym(i,km1)*gchm(i,km1) &
!                                    + GDH(I,Km1 )*(gcym(i,k)-gcym(i,km1))) )
!                rhs_q = cbmfl*(gcym(i,k)*(gcwm(i,k)-gcqm(i,k))                &
!                                 - (gcym(i,km1)*(gcwm(i,km1)-gcqm(i,km1))          &
!                                 + (GDw( I,Km1 )-gdq(i,km1,1))*(gcym(i,k)-gcym(i,km1))) )
                 tem   = cbmfl * (one - sigma(i,k))
                 tem1  = gcym(i,k)   * (one - sigma(i,k))
                 tem2  = gcym(i,km1) * (one - sigma(i,km1))
                 rhs_h = cbmfl * (tem1*gchm(i,k) - (tem2*gchm(i,km1) &
                                                 + GDH(I,Km1)*(tem1-tem2)) )
                 rhs_q = cbmfl * (tem1*(gcwm(i,k)-gcqm(i,k))         &
                               - (tem2*(gcwm(i,km1)-gcqm(i,km1))     &
                               + (GDw(I,Km1)-gdq(i,km1,1))*(tem1-tem2)) )

!              ELSE
!                rhs_h = cbmfl*(gcht(i,ctp) - (gcym(i,k-1)*gchm(i,k-1) + GDH( I,K-1 )*(gcyt(i,ctp)-gcym(i,k-1))) )
!                rhs_q = cbmfl*((gcwt(i)-gcqt(i,ctp)) - (gcym(i,k-1)*(gcwm(i,k-1)-gcqm(i,k-1)) + (GDw( I,K-1 )-gdq(i,k-1,1))*(gcyt(i,ctp)-gcym(i,k-1))) )
!              endif

!
               dqcondtem(i,km1)   = -rhs_q                             ! condensation
!              dqprectem(i,km1)   = cbmfl * (GPRCIZ(i,k) + GSNWIZ(i,k))
               dqprectem(i,km1)   = tem * (GPRCIZ(i,k) + GSNWIZ(i,k))  ! total precip production
!              dfrzprectem(i,km1) = cbmfl * GSNWIZ(i,k)
               dfrzprectem(i,km1) = tem * GSNWIZ(i,k)                  ! production of frozen precip
               dtfrztem(i,km1)    = rhs_h*oneocp                       ! heating due to freezing
! total temperature tendency due to in cloud microphysics
               dtcondtem(i,km1)   = - elocp * dqcondtem(i,km1) + dtfrztem(i,km1)

             endif ! if(k > kbi) then
           enddo   ! end of k=kbi,kk loop

         endif     ! end of if(cbmfl > zero)
    
    
! get tendencies by difference of fluxes, sum over cloud type

         do k = 1,kk
           delpinv          = grav / delp(I,k)
! cloud microphysical tendencies for single cloud type
           dtcondtem(i,k)   = dtcondtem(i,k)   * delpinv
           dqcondtem(i,k)   = dqcondtem(i,k)   * delpinv
           dqprectem(i,k)   = dqprectem(i,k)   * delpinv
           dtfrztem(i,k)    = dtfrztem(i,k)    * delpinv
! sum cloud microphysical tendencies over all cloud types
           condtermt(i,k)   = condtermt(i,k)   + dtcondtem(i,k) 
           condtermq(i,k)   = condtermq(i,k)   + dqcondtem(i,k) 
           prectermq(i,k)   = prectermq(i,k)   + dqprectem(i,k) 
           prectermfrz(i,k) = prectermfrz(i,k) + dfrzprectem(i,k) 
           frzterm(i,k)     = frzterm(i,k)     + dtfrztem(i,k)

! flux tendencies - compute the vertical flux divergence
           sfluxterm(i,k)  = sfluxterm(i,k)  - (sfluxtem(k+1)  - sfluxtem(k))  * delpinv
           qvfluxterm(i,k) = qvfluxterm(i,k) - (qvfluxtem(k+1) - qvfluxtem(k)) * delpinv
           qlfluxterm(i,k) = qlfluxterm(i,k) - (qlfluxtem(k+1) - qlfluxtem(k)) * delpinv
           qifluxterm(i,k) = qifluxterm(i,k) - (qifluxtem(k+1) - qifluxtem(k)) * delpinv
           
         enddo
        
       enddo           ! end of i loop

     endif             ! end of do_aw if !DDsigma -  end sigma computation for AW

!
! Cloud Mass Flux & Precip.
     CALL CUMFLX(IM    , IJSDIM, KMAX  ,                             & !DD dimensions
                 GMFX0 , GPRCI , GSNWI ,                             & ! output
                 QLIQ  , QICE  , GTPRC0,                             & ! output
                 CBMFX(1,CTP)  , GCYM        , GPRCIZ    , GSNWIZ ,  & ! input
                 GTPRT(1,CTP)  , GCLZ        , GCIZ      ,           & ! input
                 KB            , KT(1,CTP)   , KTMX(CTP) ,           & ! input
                 ISTS          , IENS, sigma                        )    ! input

   ENDDO          ! end of cloud type ctp loop
   
!
   do k=1,kmax
     do i=ists,iens
       GMFLX(I,k) = GMFX0(I,k) ! contains net updraft mass flux for all clouds
     enddo
   enddo
   KTMXT = 3
   DO CTP=1,NCTP
     IF (KTMX(CTP) > KTMXT) KTMXT = KTMX(CTP)
   ENDDO
   DO K=1,KTMXT
     DO I=ISTS,IENS
       CUMCLW(I,K) = QLIQ(I,K) + QICE(I,K)
       IF (CUMCLW(I,K) > zero) THEN
            FLIQC(I,K)  = QLIQ(I,K) / CUMCLW(I,K)
            FLIQOU(I,K) = FLIQC(I,K)
       ENDIF
     ENDDO
   ENDDO
!
! Cumulus Cloudiness
   CALL CUMCLD(IJSDIM, KMAX  ,                                & !DD dimensions
               CUMCLW, QLIQ  , QICE  , FLIQC  ,               & ! modified
               CUMFRC,                                        & ! output
               GMFLX , KTMXT , ISTS  , IENS    )                ! input
!
! Cloud Detrainment Heating
   CALL CUMDET(im    , IJSDIM, KMAX  , NTR   ,                 & !DD dimensions
               CMDET ,                                         & ! output
!              CMDET , GTLDET, GTIDET,                         & ! output
               GTT   , GTQ   , GTCFRC, GTU   , GTV   ,         & ! modified
!              GTQI  ,                                         & ! modified
               GDH   , GDQ   , GDCFRC, GDU   , GDV   ,         & ! input
               CBMFX , GCYT  , DELP  , GCHT  , GCQT  ,         & ! input
               GCLT  , GCIT  , GCUT  , GCVT  , GDQ(1,1,iti),   & ! input
               KT    , ISTS  , IENS, nctp, sigmai      )   ! input

!   if (lprnt) write(0,*)' after cumdet gtqi=',gtq(ipr,:,2)

!for now area fraction of the downdraft is zero, it will be computed
!  within cumdwn and applied there
! Get AW downdraft eddy flux and microphysical tendencies out of downdraft code.
   do k=1,kmax
     do i=ists,iens
       sigmad(i,k)  = zero
     enddo
   enddo

! cumulus downdraft - Melt & Freeze & Evaporation
   CALL CUMDWN(IM    , IJSDIM, KMAX  , NTR   ,           & ! DD dimensions
               GTT   , GTQ   , GTU   , GTV   ,           & ! modified
                       GMFLX ,                           & ! modified updraft+downdraft flux
!              GTQI  , GMFLX ,                           & ! modified
               GPRCP , GSNWP , GTEVP , GMDD  ,           & ! output
               GPRCI , GSNWI ,                           & ! input
               GDH   , GDW   , GDQ   , GDQ(1,1,iti) ,    & ! input
               GDQS  , GDS   , GDHS  , GDT   ,           & ! input
               GDU   , GDV   , GDZ   ,                   & ! input
               GDZM  , GCYM  , FDQS  , DELP  ,           & ! input
               sigmad, do_aw , do_awdd,                  & ! DDsigma input
               dtmelt, dtevap, dtsubl,                   & ! DDsigma input
               dtdwn , dqvdwn, dqldwn, dqidwn,           & ! DDsigma input
               KB    , KTMXT , ISTS  , IENS    )           ! input

!   if (lprnt) write(0,*)' after cumdwn gtqi=',gtq(ipr,:,2)
! here we substitute the AW tendencies into tendencies to be passed out
!  if (do_aw) then
     do k=1,kmax
       do i=ists,iens
         sigma(i,k)  = sigma(i,k) + sigmad(i,k)
       enddo
     enddo

! AW lump all heating together, compute qv term
     do k=1,kmax
       do i=ists,iens
         dqevap(i,k) = - dtevap(i,k)*cpoel - dtsubl(i,k)*cpoesub
         dtevap(i,k) =   dtevap(i,k) + dtsubl(i,k)
         dtsubl(i,k) = zero
       enddo
     enddo

!    do i=1,ijsdim
!      moistening_aw(i) = zero
!    enddo
!    DO K = 1, KMAX
!      DO I = ISTS, IENS
!        tem              = frzterm(i,k)*cpoEMELT - prectermfrz(i,k)
!        gtt(i,k)         = dtdwn(i,k)  + sfluxterm(i,k)  + condtermt(i,k) &
!                         + dtmelt(i,k) + dtevap(i,k)
!        gtq(i,k,1)       = dqvdwn(i,k) + qvfluxterm(i,k) + condtermq(i,k) &
!                         + dqevap(i,k)
!        gtq(i,k,itl)     = dqldwn(i,k) + qlfluxterm(i,k) - condtermq(i,k) &
!                         - prectermq(i,k) - tem
!        gtq(i,k,iti)     = dqidwn(i,k) + qifluxterm(i,k) + tem
! detrainment terms get zeroed
!        gtldet(i,k)      = zero
!        gtidet(i,k)      = zero
! column-integrated total water tendency - used to impose water conservation
!        moistening_aw(i) = moistening_aw(i)                               &
!                         + (gtq(i,k,1)+gtq(i,k,itl)+gtq(i,k,iti)) * delp(i,k)*gravi
!      ENDDO
!    ENDDO
!
! This code ensures conservation of water. In fact, no adjustment of the precip
!   is occuring now which is a good sign! DD
!    DO I=ISTS,IENS
!      if(gprcp(i,1)+gsnwp(i,1) > 1.e-12_r8) then
!        moistening_aw(i) = - moistening_aw(i) / (gprcp(i,1)+gsnwp(i,1))
!      endif
!      if (abs(1.0-moistening_aw(i)) > 0.3 .and. gprcp(i,1) > 0.0) &
!         write(1000+mype,*)' moistening_aw=', &
!         moistening_aw(i),' i=',i,' xlon=',xlon(i),' xlat=',xlat(i),' kdt=',kdt&
!        , ' gprcp=',gprcp(i,1:5)
!    ENDDO
!     write(1000+mype,*)' moistening_aw=',moistening_aw
!    do k=1,kmax
!      DO I = ISTS, IENS
!        gprcp(i,k) = gprcp(i,k) * moistening_aw(i)
!        gsnwp(i,k) = gsnwp(i,k) * moistening_aw(i)
!      ENDDO
!    enddo

!  else

! Cloud Subsidence Heating
     CALL CUMSBH(IM    , IJSDIM, KMAX  , NTR   ,           & !DD dimensions
                 GTT   , GTQ   ,                           & ! modified
!                GTT   , GTQ   , GTQI  ,                   & ! modified
                 GTU   , GTV   ,                           & ! modified
                 GDH   , GDQ   , GDQ(1,1,iti)  ,           & ! input
                 GDU   , GDV   ,                           & ! input
                 DELP  , GMFLX , GMFX0 ,                   & ! input
                 KTMXT , CPRES , ISTS  , IENS )   ! input

!   if (lprnt) write(0,*)' after cumsbh gtqi=',gtq(ipr,:,2)
!  endif
!
! for now the following routines appear to be of no consequence to AW - DD
!
! Tracer Updraft
   CALL CUMUPR(im    , IJSDIM, KMAX  , NTR   ,           & !DD dimensions
               GTQ   , GPRCC ,                           & ! modified
               GDQ   , CBMFX , ELAM  , GDZ   , GDZM  ,   & ! input
               GCYM  , GCYT  , GCQT  , GCLT  , GCIT  ,   & ! input
               GTPRT , GTEVP , GTPRC0,                   & ! input
               KB    , KBMX  , KT    , KTMX  , KTMXT ,   & ! input
               DELP  , OTSPT1, ISTS  , IENS,             & ! input
               fscav, fswtr, nctp)
!
! Tracer Downdraft
   CALL CUMDNR(im    ,IJSDIM , KMAX  , NTR   ,           & !DD dimensions
               GTQ   ,                                   & ! modified
               GDQ   , GMDD  , DELP  ,                   & ! input
               KTMXT , OTSPT1, ISTS  , IENS )              ! input
!
! Tracer Subsidence
   CALL CUMSBR(im    , IJSDIM,KMAX  , NTR   ,            & !DD dimensions
               GTQ   ,                                   & ! modified
               GDQ   , DELP  ,                           & ! input
               GMFLX , KTMXT , OTSPT2,                   & ! input
               ISTS  , IENS            )                   ! input

!   if (lprnt) write(0,*)' after cumsbr gtqi=',gtq(ipr,:,2)
!
!  do k=1,kmax
!    do i=ISTS,IENS
!      GTQ(I,k,ITI) = GTQI(I,k)
!    enddo
!  enddo
!
! Tracer mass fixer without detrainment
   CALL CUMFXR(IM    , IJSDIM, KMAX  , NTR   ,           & !DD dimensions
               GTQ   ,                                   & ! modified
               GDQ   , DELP  , DELTA , KTMXT , IMFXR,    & ! input
               ISTS  , IENS                            )   ! input
!
!  do k=1,kmax
!    do i=ISTS,IENS
!      GTQL(I,k) = GTQ(I,k,ITL) + GTLDET(I,k) + GTIDET(I,k)
!    enddo
!  enddo
!
! Tracer mass fixer with detrainment
!  CALL CUMFXR1(IM  , IJSDIM, KMAX  ,                           & !DD dimensions
!               GTQL        ,                                   & ! modified
!               GDQ(1,1,ITL), DELP, DELTA, KTMXT, IMFXR(ITL),   & ! input
!               ISTS        , IENS                            )   ! input
!
   DO K=1,KMAX
     DO I=ISTS, IENS
!      GTLDET(I,k) = GTQL(I,k) - GTQ(I,k,ITL) - GTIDET(I,k)

! tendencies of subgrid PDF (turned off)
!      GDCLDX = GDCFRC( I,K ) + GTCFRC( I,K )*DELTA
!      GDCLDX = MIN( MAX( GDCLDX, 0.D0 ), one )
!      GTCFRC( I,K ) = ( GDCLDX - GDCFRC( I,K ) )/DELTA
!
!      GDMU2X = GDQ( I,K,IMU2 ) + GTQ( I,K,IMU2 )*DELTA
!      GDMU2X = MIN( MAX( GDMU2X,VARMIN ),VARMAX )
!      GDMU3X = GDQ( I,K,IMU3 ) + GTQ( I,K,IMU3 )*DELTA
!      GDMU3X = MIN( MAX( GDMU3X,-SKWMAX ),SKWMAX )
!      GTQ( I,K,IMU2 ) = ( GDMU2X - GDQ( I,K,IMU2 ))/DELTA
!      GTQ( I,K,IMU3 ) = ( GDMU3X - GDQ( I,K,IMU3 ))/DELTA
!
       tem = DELP(I,K)*GRAVI
       HBGT(I) = HBGT(I) + (CP*GTT(I,K) + EL*GTQ(I,K,1)                         &
                           - EMELT*GTQ(I,K,ITI)) * tem
!                          - EMELT*(GTQ(I,K,ITI)+GTIDET(I,K))) * tem
       WBGT(I) = WBGT(I) + (GTQ(I,K,1)   + GTQ(I,K,ITL) + GTQ(I,K,ITI)) * tem 
!                                        + GTLDET(I,K)  + GTIDET(I,K)) * tem
     ENDDO
   ENDDO
  
 
! here we substitute the AW tendencies into tendencies to be passed out
   if(do_aw) then
     do i=1,ijsdim
       moistening_aw(i) = zero
     enddo
     tem2 = one / delta
     DO K=1,KMAX
       DO I=ISTS,IENS
         tem = frzterm(i,k)*cpoEMELT - prectermfrz(i,k)
         gtt(i,k)         = dtdwn(i,k)  + sfluxterm(i,k)  + condtermt(i,k)         &
                          + dtmelt(i,k) + dtevap(i,k)
         gtq(i,k,1)       = dqvdwn(i,k) + qvfluxterm(i,k) + condtermq(i,k)         &
                          + dqevap(i,k)
         gtq(i,k,itl)     = dqldwn(i,k) + qlfluxterm(i,k) - condtermq(i,k)         &
                          - prectermq(i,k) - tem
         gtq(i,k,iti)     = dqidwn(i,k) + qifluxterm(i,k) + tem
! detrainment terms get zeroed
!        gtldet(i,k)      = zero
!        gtidet(i,k)      = zero

         tem1             = - gdq(i,k,itl)*tem2
         if (gtq(i,k,itl) < tem1) then
           tem3           = gtq(i,k,itl) - tem1
           gtq(i,k,1)     = gtq(i,k,1) + tem3
           gtq(i,k,itl)   = tem1
           gtt(i,k)       = gtt(i,k) - elocp*tem3
         endif
         tem1             = - gdq(i,k,iti)*tem2
         if (gtq(i,k,iti) < tem1) then
           tem3           = gtq(i,k,iti) - tem1
           gtq(i,k,1)     = gtq(i,k,1) + tem3
           gtq(i,k,iti)   = tem1
           gtt(i,k)       = gtt(i,k) - esubocp*tem3
         endif
!        tem1             = - gdq(i,k,1)*tem2
!        if (gtq(i,k,1) < tem1) then
!          gtt(i,k)       = gtt(i,k) + elocp*(gtq(i,k,1)-tem1)
!          gtq(i,k,1)     = tem1
!        endif
           
! column-integrated total water tendency - to be used to impose water conservation
         moistening_aw(i) = moistening_aw(i)                                       &
                          + (gtq(i,k,1)+gtq(i,k,itl)+gtq(i,k,iti)) * delp(i,k)/grav
       ENDDO
     ENDDO
   endif

!   if (lprnt) then
!     write(0,*)' after doaw dqvdwn=',dqvdwn(ipr,:)
!     write(0,*)' after doaw qvfluxterm=',qvfluxterm(ipr,:)
!     write(0,*)' after doaw dqevap=',dqevap(ipr,:)
!     write(0,*)' after doaw condtermq=',condtermq(ipr,:)
!     write(0,*)' after doaw dqidwn=',dqidwn(ipr,:)
!     write(0,*)' after doaw qifluxterm=',qifluxterm(ipr,:)
!     write(0,*)' after doaw prectermfrz=',prectermfrz(ipr,:)
!     write(0,*)' after doaw frzterm=',frzterm(ipr,:)
!     write(0,*)' after doaw gtqv=',gtq(ipr,:,1)
!     write(0,*)' after doaw gtqi=',gtq(ipr,:,2)
!     write(0,*)' after doaw gtql=',gtq(ipr,:,3)
!   endif
!
   DO I=ISTS,IENS
     HBGT(I)  = HBGT(I) - EMELT*GSNWC(I)
     WBGT(I)  = WBGT(I) + GPRCC(I,1) + GSNWC(I)
!    CTOPP(I) = 1.D6
   ENDDO
!
!  The following commented out because they are unused
!  DO CTP=1,NCTP
!    DO I=ISTS, IENS
!      kk = kt(i,ctp)
!      IF (KK > KB(I) ) THEN
!        CUMHGT(I,CTP) = GDZ(I,KK)
!        CTOPP(I)      = MIN(CTOPP(I), GDP(I,KK))
!      ELSE
!        CUMHGT (I,CTP) = -999.D0
!      ENDIF
!    ENDDO
!  ENDDO
!  DO I=ISTS,IENS
!    IF(CTOPP(I) >= 1.D6) THEN
!     CTOPP(I) = -999.D0
!    ENDIF
!  ENDDO
!
! This code ensures conservation of water. In fact, no adjustment of the precip
!   is occuring now which is a good sign! DD
   if(do_aw .and. adjustp) then
     DO I = ISTS, IENS
       if(gprcp(i,1)+gsnwp(i,1) > 1.e-12_r8) then
         moistening_aw(i) = - moistening_aw(i) / (gprcp(i,1)+gsnwp(i,1))
       endif
     ENDDO
     do k=1,kmax
       DO I = ISTS, IENS
         gprcp(i,k) = gprcp(i,k) * moistening_aw(i)
         gsnwp(i,k) = gsnwp(i,k) * moistening_aw(i)
       ENDDO
     enddo

!    if (lprnt) then
!      write(1000+mype,*)' moistening_aw=',moistening_aw(1:ijsdim)
!      write(1000+mype,*)' gprcp=',gprcp(1:ijsdim,1)
!    endif

   endif


! commnting out becasue these are not used
!  DO K=1,KMAX
!    kp1 = min(k+1,kmax)
!    DO I=ISTS,IENS
!      GPRCPF(I,K) = half * (GPRCP(I,K) + GPRCP(I,KP1))
!      GSNWPF(I,K) = half * (GSNWP(I,K) + GSNWP(I,KP1))
!    ENDDO
!  ENDDO
!
   do i=ISTS,IENS
      GPRCC(I,1) = GPRCP(I,1)
      GSNWC(I  ) = GSNWP(I,1)
   enddo
   do k=1,kmax
     do i=ISTS,IENS
       GTPRP(I,k) = GPRCP(I,k) + GSNWP(I,k)
     enddo
   enddo
!

!COSP
!necessary?
!  DO K=1,KMAX
!    DO I=ISTS,IENS
!      QLIQC(I,K) = QLIQ(I,K)
!      QICEC(I,K) = QICE(I,K)
!    ENDDO
!  ENDDO
!
!  IF ( OOUT1 .OR. OOUT2 ) THEN
     DO I=ISTS,IENS
       CAPE(i)  = zero
       CIN(i)   = zero
       JBUOY(i) = 0
     enddo
     DO K=2,KMAX
       DO I=ISTS,IENS
         IF (K >= KB(I)) THEN
            BUOY = (GDH(I,1)-GDHS(I,K)) / ((one+ELOCP*FDQS(I,K)) * CP*GDT(I,K))
         ELSE
            BUOY = (GDS(I,1)-GDS(I,K)) / (CP*GDT(I,K))
         END IF
            DELZ = GDZM(I,K+1) - GDZM(I,K)
         IF (BUOY > zero .AND. JBUOY(I) /=  0) THEN
            CAPE(I) = CAPE(I) + BUOY*GRAV*DELZ
              JBUOY(I) = 2
         ELSEIF (BUOY < zero .AND. JBUOY(I) /= 2) THEN
            CIN(I) = CIN(I) - BUOY*GRAV*DELZ
            JBUOY(I) = 1
         ENDIF
       ENDDO
     ENDDO
     DO I=ISTS,IENS
       IF (JBUOY(I) /= 2) CIN(I) = -999.D0
     ENDDO

!DD provide GFS with a separate downdraft mass flux
     DO K=1,KMAX
       DO I=ISTS,IENS
         GMFX1(I,K) = GMFX0(I,K) - GMFLX(I,K)
       ENDDO
     ENDDO
     
!
      END SUBROUTINE CS_CUMLUS
!***********************************************************************
      SUBROUTINE CUMBAS                            & !! cloud base
               ( IJSDIM, KMAX  ,                   & !DD dimensions
                 KB    , GCYM  , KBMX  ,           & ! output
                 GCHB  , GCWB  , GCUB  , GCVB  ,   & ! output
                 GCIB  ,                           & ! output
                 GDH   , GDW   , GDHS  , GDQS  ,   & ! input
                 GDQI  , GDU   , GDV   , GDZM  ,   & ! input
                 GDPM  , FDQS  , GAM   ,           & ! input
                 ISTS  , IENS , gctbl, gcqbl ,gdq, gcwbl, gcqlbl, gcqibl   )   ! input  !DDsigmadiag add updraft profiles below cloud base
!
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJSDIM, KMAX ! DD, for GFS, pass in
!
!   [OUTPUT]
      INTEGER    KB    (IJSDIM)         ! cloud base
      REAL(r8)   GCYM  (IJSDIM, KMAX)   ! norm. mass flux (half lev)
      INTEGER    KBMX
      REAL(r8)   GCHB  (IJSDIM)         ! cloud base MSE
      REAL(r8)   GCWB  (IJSDIM)         ! cloud base total water
      REAL(r8)   GCUB  (IJSDIM)         ! cloud base U
      REAL(r8)   GCVB  (IJSDIM)         ! cloud base V
      REAL(r8)   GCIB  (IJSDIM)         ! cloud base ice

!DDsigma added to arglist for AW, subcloud updraft profiles: temperature, water vapor
!                               total water, cloud water, and cloud ice respectively
      REAL(r8), dimension(ijsdim,kmax) :: gctbl, gcqbl, gcwbl, gcqlbl, gcqibl   !DDsigmadiag
!
!   [INPUT]
      REAL(r8)   GDH   (IJSDIM, KMAX)        ! moist static energy
      REAL(r8)   GDW   (IJSDIM, KMAX)        ! total water
      REAL(r8)   GDq   (IJSDIM, KMAX)        ! water vapor !DDsigmadiag
      REAL(r8)   GDHS  (IJSDIM, KMAX)        ! saturate MSE
      REAL(r8)   GDQS  (IJSDIM, KMAX)        ! saturate humidity
      REAL(r8)   GDQI  (IJSDIM, KMAX)        ! cloud ice
      REAL(r8)   GDU   (IJSDIM, KMAX)        ! u-velocity
      REAL(r8)   GDV   (IJSDIM, KMAX)        ! v-velocity
      REAL(r8)   GDZM  (IJSDIM, KMAX+1)      ! Altitude (half lev)
      REAL(r8)   GDPM  (IJSDIM, KMAX+1)      ! pressure (half lev)
      REAL(r8)   FDQS  (IJSDIM, KMAX)
      REAL(r8)   GAM   (IJSDIM, KMAX)
      INTEGER    ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)   CBASE (IJSDIM)              ! one over cloud base height
!     REAL(r8)   CBASEP(IJSDIM)              ! cloud base pressure
      REAL(r8)   DELZ, QSL, GAMX, wrk
      REAL(r8), dimension(ijsdim,kmax) :: gchbl   !DDsigmadiag
      real(r8), dimension(ijsdim)      :: gcqb
      INTEGER    I, K, kp1
!
!   [INTERNAL PARM]
      INTEGER :: KMAXM1
      INTEGER :: KLCLB                       !! LCL base level
      INTEGER :: KCB                         !! fix cloud bottom
      INTEGER :: KBMAX                       !! cloud base max
      INTEGER :: KBOFS                       !! cloud base offset

      KMAXM1 = KMAX-1
      KLCLB  = 1                   ! LCL base level
      KCB    = 0                   ! fix cloud bottom
      KBMAX  = KMAXM1              ! cloud base max
      KBOFS  = 0                   ! cloud base offset
!
      do k=1,kmax
        do i=ists,iens
          GCYM(I,k) = zero
        enddo
      enddo
!
      IF (KCB > 0) THEN
        DO I=ISTS,IENS
          KB(I) = KCB
        ENDDO
      ELSE
        DO I=ISTS,IENS
          KB(I) = KBMAX
        ENDDO
        DO K=KBMAX-1,KLCLB+1,-1
          DO I=ISTS,IENS
            GAMX = FDQS(I,K) / (one+GAM(I,K)) * oneocp
            QSL  = GDQS(I,K) + GAMX * (GDH(I,KLCLB)-GDHS(I,K))
            IF (GDW(I,KLCLB) >= QSL) THEN
              KB(I) = K + KBOFS
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
      KBMX = 1
      DO I=ISTS,IENS
        KBMX = MAX(KBMX, KB(I))
        CBASE (I) = one / (GDZM(I,KB(I)) - GDZM(I,1))
!       CBASEP(I) = GDPM(I,KB(I))
      ENDDO
!
      DO K=1,KBMX
        DO I=ISTS,IENS
          IF (K <= KB(I)) THEN
            GCYM(I,K) = sqrt((GDZM(I,K)-GDZM(I,1)) *  CBASE(i))
          ENDIF
        ENDDO
      ENDDO
!
      DO I=ISTS,IENS
        GCHB(I) = zero
        GCWB(I) = zero
        GCUB(I) = zero
        GCVB(I) = zero
        GCIB(I) = zero
        GCQB(I) = zero
      ENDDO
      do k=1,kmax
        do i=ists,iens
          GChbl(i,k)  = zero
          gcqbl(i,k)  = zero
          gcqlbl(i,k) = zero
          gcqibl(i,k) = zero
          gctbl(i,k)  = zero
          gcwbl(i,k)  = zero
        enddo
      enddo
!
      DO K=1,KBMX
        kp1 = min(k+1, kmax)
        DO I=ISTS,IENS
          IF (K < KB(I)) THEN
            DELZ    = GCYM(I,Kp1) - GCYM(I,K)
            GCHB(I) = GCHB(I) + DELZ * GDH (I,K)
            GCWB(I) = GCWB(I) + DELZ * GDW (I,K)
            GCUB(I) = GCUB(I) + DELZ * GDU (I,K)
            GCVB(I) = GCVB(I) + DELZ * GDV (I,K)
            GCIB(I) = GCIB(I) + DELZ * GDQI(I,K)
            GCqB(I) = GCqB(I) + DELZ * GDQ (I,K)
! get subcloud profiles to pass out and do AW eddy flux tendencies
!   removing the normalized mass flux weighting
            wrk       = one /  gcym(i,kp1)
            gchbl(i,kp1)  = gchb(i) * wrk
            gcqbl(i,kp1)  = gcqb(i) * wrk
            gcqibl(i,kp1) = gcib(i) * wrk
            gcwbl(i,kp1)  = gcwb(i) * wrk
            gcqlbl(i,kp1) = gcwbl(i,kp1) - (gcqibl(i,kp1)+gcqbl(i,kp1))
            gctbl(i,kp1)  = (gchbl(i,kp1) - grav*gdzm(i,kp1) - el*gcqbl(i,kp1)) * oneocp
          ENDIF
        ENDDO
      ENDDO
!
      END SUBROUTINE CUMBAS
!***********************************************************************
      SUBROUTINE CUMUP                             & !! in-cloud properties
               ( IJSDIM, KMAX  , NTR   ,           & !DD dimensions
                 ACWF  , ELAM  ,                   & ! output
                 GCLZ  , GCIZ  , GPRCIZ, GSNWIZ,   & ! output
                 GCYT  , GCHT  , GCQT  ,           & ! output
                 GCLT  , GCIT  , GTPRT ,           & ! output
                 GCUT  , GCVT  ,                   & ! output
                 KT    , KTMX  ,                   & ! output
                 GCYM  ,                           & ! modified
                 wcv   ,                           & ! !DDsigma new output
                 GCHB  , GCWB  , GCUB  , GCVB  ,   & ! input   !DDsigmadiag
                 GCIB  ,                           & ! input
                 GDU   , GDV   , GDH   , GDW   ,   & ! input
                 GDHS  , GDQS  , GDT   , GDTM  ,   & ! input
                 GDQ   , GDQI  , GDZ   , GDZM  ,   & ! input
                 GDPM  , FDQS  , GAM   , GDZTR ,   & ! input
                 CPRES , WCB   ,                   & ! input
!                CPRES , WCB   , ERMR  ,           & ! input
                 KB    , CTP   , ISTS  , IENS,     &  ! input
                 gctm  , gcqm, gcwm, gchm, gcwt, gclm, gcim, &
                 cbmfx , dtcond, dqcond, dtfrz    )   !DDsigmadiag
!
!DD AW the above line of arguments were previously local, and often scalars.
!  Dimensions were added to them to save profiles for each grid point.
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, KMAX, NTR     ! DD, for GFS, pass in
!
!   [OUTPUT]
      REAL(r8)   ACWF  (IJSDIM)           ! cloud work function
      REAL(r8)   ELAM  (IJSDIM, KMAX)     ! entrainment (rate*massflux)
      REAL(r8)   GCLZ  (IJSDIM, KMAX)     ! cloud liquid water*eta
      REAL(r8)   GCIZ  (IJSDIM, KMAX)     ! cloud ice*eta
      REAL(r8)   GPRCIZ(IJSDIM, KMAX)     ! rain generation*eta
      REAL(r8)   GSNWIZ(IJSDIM, KMAX)     ! snow generation*eta
      REAL(r8)   GCYT  (IJSDIM)           ! norm. mass flux @top
      REAL(r8)   GCHT  (IJSDIM)           ! cloud top MSE*eta
      REAL(r8)   GCQT  (IJSDIM)           ! cloud top moisture*eta
      REAL(r8)   GCLT  (IJSDIM)           ! cloud top liquid water*eta
      REAL(r8)   GCIT  (IJSDIM)           ! cloud top ice*eta
      REAL(r8)   GTPRT (IJSDIM)           ! cloud top (rain+snow)*eta
      REAL(r8)   GCUT  (IJSDIM)           ! cloud top u*eta
      REAL(r8)   GCVT  (IJSDIM)           ! cloud top v*eta
      REAL(r8)   GCwT  (IJSDIM)           ! cloud top v*eta
      INTEGER    KT    (IJSDIM)           ! cloud top
      INTEGER    KTMX                     ! max of cloud top
      REAL(r8)   WCV   (IJSDIM, KMAX)     ! updraft velocity (half lev) !DD sigma make output
!
!   [MODIFIED]
      REAL(r8)   GCYM  (IJSDIM, KMAX)     ! norm. mass flux
!
!   [INPUT]
      REAL(r8)   GCHB  (IJSDIM)           ! MSE at cloud base
      REAL(r8)   GCWB  (IJSDIM)           ! total water @cloud base
      REAL(r8)   GCUB  (IJSDIM)           ! U at cloud base
      REAL(r8)   GCVB  (IJSDIM)           ! V at cloud base
      REAL(r8)   GCIB  (IJSDIM)           ! cloud ice at cloud base
      REAL(r8)   GDU   (IJSDIM, KMAX)     ! U
      REAL(r8)   GDV   (IJSDIM, KMAX)     ! V
      REAL(r8)   GDH   (IJSDIM, KMAX)     ! moist static energy
      REAL(r8)   GDW   (IJSDIM, KMAX)     ! total water
      REAL(r8)   GDHS  (IJSDIM, KMAX)     ! saturation MSE
      REAL(r8)   GDQS  (IJSDIM, KMAX)     ! saturation q
      REAL(r8)   GDT   (IJSDIM, KMAX)     ! T
      REAL(r8)   GDTM  (IJSDIM, KMAX+1)   ! T (half lev)
      REAL(r8)   GDQ   (IJSDIM, KMAX, NTR)! q  !!DDsigmadiag
      REAL(r8)   GDQI  (IJSDIM, KMAX)     ! cloud ice
      REAL(r8)   GDZ   (IJSDIM, KMAX)     ! z
      REAL(r8)   GDZM  (IJSDIM, KMAX+1)   ! z (half lev)
      REAL(r8)   GDPM  (IJSDIM, KMAX+1)   ! p (half lev)
      REAL(r8)   FDQS  (IJSDIM, KMAX)
      REAL(r8)   GAM   (IJSDIM, KMAX)
      REAL(r8)   GDZTR (IJSDIM)           ! tropopause height
      REAL(r8)   CPRES                    ! pres. fac. for cum. fric.
      REAL(r8)   WCB(ijsdim)              ! updraft velocity**2 @base
!     REAL(r8)   ERMR                     ! entrainment rate (ASMODE)
      INTEGER    KB    (IJSDIM)
      INTEGER    CTP, ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     myGCHt                 ! cloud top h *eta (half lev)
      REAL(r8)     GCHMZ (IJSDIM, KMAX)   ! cloud h *eta (half lev)
      REAL(r8)     GCWMZ (IJSDIM, KMAX)   ! cloud Qt*eta (half lev)
      REAL(r8)     GCqMZ (IJSDIM, KMAX)   ! cloud qv*eta (half lev)
      REAL(r8)     GCUMZ (IJSDIM, KMAX)   ! cloud U *eta (half lev)
      REAL(r8)     GCVMZ (IJSDIM, KMAX)   ! cloud V *eta (half lev)
      REAL(r8)     GCIMZ (IJSDIM, KMAX)   ! cloud Qi*eta (half lev)
      REAL(r8)     GTPRMZ(IJSDIM, KMAX)   ! rain+snow *eta (half lev)
!
      REAL(r8)     BUOY  (IJSDIM, KMAX)     ! buoyancy
      REAL(r8)     BUOYM (IJSDIM, KMAX)   ! buoyancy (half lev)
      REAL(r8)     WCM   (IJSDIM, KMAX)   ! updraft velocity**2 (half lev)
!DD sigma make output     REAL(r8)     WCV   ( IJSDIM, KMAX+1 )   !! updraft velocity (half lev)
      REAL(r8)     GCY   (IJSDIM, KMAX)     ! norm. mass flux
      REAL(r8)     ELAR  (IJSDIM, KMAX)     ! entrainment rate
!
      REAL(r8)     GCHM  (IJSDIM, KMAX)   ! cloud MSE (half lev)
      REAL(r8)     GCWM  (IJSDIM, KMAX)   ! cloud Qt  (half lev)  !DDsigmadiag
      REAL(r8)     GCTM  (IJSDIM, KMAX)   ! cloud T (half lev)   !DDsigmadiag make output
      REAL(r8)     GCQM  (IJSDIM, KMAX)   ! cloud q (half lev)   !DDsigmadiag make output
      REAL(r8)     dtcond(IJSDIM, KMAX)   ! in cloud condensation heating DDsigmadiag make output
      REAL(r8)     dqcond(IJSDIM, KMAX)   ! in cloud condensation water vapor tendency   !DDsigmadiag make output
      REAL(r8)     dtfrz (IJSDIM, KMAX)   ! in cloud temperature tendency due to freezing  !DDsigmadiag make output
      REAL(r8)     cbmfx (IJSDIM)         ! cloud base mass flux   !DDsigmadiag make output
      REAL(r8)     GCLM  (IJSDIM, KMAX)   ! cloud liquid ( half lev)
      REAL(r8)     GCIM  (IJSDIM, KMAX)   ! cloud ice (half lev)
      REAL(r8)     GCUM  (IJSDIM, KMAX)   ! cloud U (half lev)
      REAL(r8)     GCVM  (IJSDIM, KMAX)   ! cloud V (half lev)
!
      REAL(r8), dimension(IJSDIM) :: WCM_, ELARM1, GDZMKB
      REAL(r8)     GDQSM, GDHSM, GDQM, GDSM, GDCM, FDQSM, GCCM, gdtrm,  &
                   DELZ, ELADZ, DCTM , CPGMI, DELC, FICE, ELARM2,GCCMZ, &
                   PRECR, GTPRIZ, DELZL, GCCT, DCT, WCVX, PRCZH, wrk
      INTEGER      K, I, kk, km1, kp1
      CHARACTER    CTNUM*2
!
!DD#ifdef OPT_CUMBGT
!DD   REAL(r8)     HBGT  (IJSDIM)           ! heat budget
!DD   REAL(r8)     WBGT  (IJSDIM)           ! water budget
!DD   REAL(r8)     PBGT  (IJSDIM)           ! precipitation budget
!DD   REAL(r8)     MBGT  (IJSDIM)           ! mass budget
!DD   REAL(r8)     GTPRX (IJSDIM)           ! (rain+snow)*eta at top
!DD   REAL(r8)     GSNWT (IJSDIM)           ! cloud top snow*eta
!DD   REAL(r8)     HBMX, WBMX, PBMX, MBMX
!DD   SAVE       HBMX, WBMX, PBMX, MBMX
!DD#endif
!
!   [INTERNAL PARAM]

      REAL(r8), SAVE :: CLMP
!DD   REAL(r8) ::  PRECZ0 = 1.5e3_r8   ! move to module scope for tuning
!DD   REAL(r8) ::  PRECZH = 4.e3_r8    ! move to module scope for tuning
      REAL(r8) ::  ZTREF  = 1._r8
      REAL(r8) ::  PB     = 1._r8
!m    REAL(r8) ::  TAUZ   = 5.e3_r8
      REAL(r8) ::  TAUZ   = 1.e4_r8
      REAL(r8) ::  ELMD   = 2.4e-3     ! for Neggers and Siebesma (2002)
      REAL(r8) ::  ELAMIN = zero       ! min. of entrainment rate
      REAL(r8) ::  ELAMAX = 4.e-3      ! max. of entrainment rate
!m    REAL(r8) ::  ELAMAX = 5.e-3      ! max. of entrainment rate
      REAL(r8) ::  WCCRT  = zero
!m    REAL(r8) ::  WCCRT  = 0.01
      REAL(r8) ::  TSICE  = 268.15_r8  ! compatible with macrop_driver
      REAL(r8) ::  TWICE  = 238.15_r8  ! compatible with macrop_driver
      REAL(r8) ::  EPSln  = 1.e-10
      
!     REAL(r8) ::  esat, tem
      REAL(r8) ::  esat, tem, rhs_h, rhs_q

      LOGICAL, SAVE :: OFIRST = .TRUE.
!
!   [INTERNAL FUNC]
      REAL(r8)     FPREC   ! precipitation ratio in condensate
      REAL(r8)     FRICE   ! ice ratio in cloud water
      REAL(r8)     Z       ! altitude
      REAL(r8)     ZH      ! scale height
      REAL(r8)     T       ! temperature
!
      FPREC(Z,ZH) = MIN(MAX(one-EXP(-(Z-PRECZ0)/ZH), zero), one)
      FRICE(T)    = MIN(MAX((TSICE-T)/(TSICE-TWICE), zero), one)
!
! Note: iteration is not made to diagnose cloud ice for simplicity
!
      IF (OFIRST) THEN
        CLMP   = 2.D0*(one-CLMD)*PA
        OFIRST = .FALSE.
      ENDIF
   
      do i=ists,iens
        ACWF (I) = zero
        GCYT (I) = zero
        GCHT (I) = zero
        GCQT (I) = zero
        GCLT (I) = zero
        GCIT (I) = zero
        GTPRT(I) = zero
        GCUT (I) = zero
        GCVT (I) = zero
        GCwT (I) = zero
      enddo
      do k=1,kmax
        do i=ists,iens
          ELAM  (I,k) = unset_r8
          GCLZ  (I,k) = zero
          GCIZ  (I,k) = zero
          GPRCIZ(I,k) = zero
          GSNWIZ(I,k) = zero
!
          GCHMZ (I,k) = zero
          GCWMZ (I,k) = zero
          GCqMZ (I,k) = zero
          GCIMZ (I,k) = zero
          GCUMZ (I,k) = zero
          GCVMZ (I,k) = zero
          GTPRMZ(I,k) = zero

          dtcond(i,k) = zero
          dqcond(i,k) = zero
          dtfrz(i,k)  = zero
!
          BUOY  (I,k) = unset_r8
          BUOYM (I,k) = unset_r8
          WCM   (I,k) = unset_r8
          WCV   (I,k) = unset_r8
          GCY   (I,k) = unset_r8
          ELAR  (I,k) = unset_r8
!
          GCHM  (I,k) = unset_r8
          GCWM  (I,k) = unset_r8
          GCTM  (I,k) = unset_r8
          GCQM  (I,k) = unset_r8
          GCLM  (I,k) = unset_r8
          GCIM  (I,k) = unset_r8
          GCUM  (I,k) = unset_r8
          GCVM  (I,k) = unset_r8
        enddo
      enddo

!#ifdef SYS_SX
      DO K=1,KMAX
        DO I=ISTS, IENS
          IF (K > KB(I)) THEN
            GCYM(I,K) = zero
          ENDIF
        ENDDO
      ENDDO
!#else
!      DO I=ISTS,IENS
!        GCYM(I,KB(I)+1:KMAX) = zero
!      ENDDO
!#endif
      DO I=ISTS,IENS
        GDZMKB(I) = GDZM(I,KB(I))     ! cloud base height
     ENDDO
!
!     < cloud base properties >
!
      DO I=ISTS,IENS
        K = KB(I)
        GCHM(I,K) = GCHB(I)
        GCWM(I,K) = GCWB(I)
        WCM (I,K) = WCB(i)
        GCUM(I,K) = GCUB(I)
        GCVM(I,K) = GCVB(I)
!
        esat  = min(gdpm(i,k), fpvs(gdtm(i,k)))
        GDQSM = min(EPSV*esat/max(gdpm(i,k)+epsvm1*esat, 1.0e-10), 1.0)
        gdsm  = CP*GDTM(I,K) + GRAV*GDZMKB(I)        ! dse
        GDHSM = gdsm + EL*GDQSM                      ! saturated mse
!       FDQSM = FDQSAT(GDTM(I,K), GDQSM)
        tem   = one / GDTM(I,K)
        FDQSM = GDQSM * tem * (fact1 + fact2*tem)
!
        tem       = one / (CP+EL*FDQSM)
        DCTM      = (GCHB(I) - GDHSM) * tem
        GCQM(I,K) = min(GDQSM + FDQSM*DCTM, GCWM(I,K))
        GCCM      = MAX(GCWM(I,K)-GCQM(I,K), zero)
!       GCTM(I,K) = GDT(I,K) + DCTM                  ! old
        GCTM(I,K) = (GCHB(I) - gdsm - el*gcqm(i,k)) * oneocp + dctm  ! new
!
        GCIM(I,K) = FRICE(GCTM(I,K)) * GCCM          ! cloud base ice
        GCLM(I,K) = MAX(GCCM-GCIM(I,K), zero)        ! cloud base liquid
        GCHM(I,K) = GCHM(I,K) + EMELT * (GCIM(I,K)-GCIB(I))
        DCTM      = (GCHM(I,K) - GDHSM) * tem
!       GCTM(I,K) = dctm + GDT(I,K) + gocp*gdzm(i,k) !DD old AW convert to DSE
        GCTM(I,K) = dctm + (GCHB(I) - el*gcqm(i,k)) * oneocp ! new, make dse
!
        GDQM  = half * (GDQ(I,K,1)     + GDQ(I,K-1,1))
        GDCM  = half * (GDQ(I,K,ITL)   + GDQI(I,K)                &
                     +  GDQ(I,K-1,ITL) + GDQI(I,K-1))
                       
!
        BUOYM(I,K) = (DCTM/GDTM(I,K) + EPSVT*(GCQM(I,K)-GDQM) - GCCM + GDCM )*GRAV
!
!DD#ifdef OPT_ASMODE
!DD     ELARM1(I) = ERMR
!DD#elif defined OPT_NS02
!DD     ELARM1(I) = ELMD / SQRT(WCM(I,K))
!DD#else
        ELARM1(I) = CLMD*PA*BUOYM(I,K)/WCM(I,K)
!DD#endif
        ELARM1(I) = MIN(MAX(ELARM1(I), ELAMIN), ELAMAX)
!
        GCHMZ (I,K)   = GCHM(I,K)
        GCWMZ (I,K)   = GCWM(I,K)
        GCqMZ (I,K)   = GCqM(I,K)
        GCUMZ (I,K)   = GCUM(I,K)
        GCVMZ (I,K)   = GCVM(I,K)
        GCIMZ (I,K)   = GCIM(I,K)
        WCM_(I)       = WCM(I,K)
      ENDDO
!
!     < in-cloud properties >
!
      DO K=3,KMAX
        km1 = k - 1
        DO I=ISTS,IENS
          IF (K > KB(I) .AND. WCM_(I) > WCCRT) THEN
            WCV(I,KM1) = SQRT(MAX(WCM_(I), zero))
            DELZ       = GDZM(I,K) - GDZM(I,KM1)
            GCYM(I,K)  = GCYM(I,KM1) * EXP(ELARM1(I)*DELZ)
            ELADZ      = GCYM(I,K) - GCYM(I,KM1)
!
            GCHMZ(I,K) = GCHMZ(I,KM1) + GDH(I,KM1)*ELADZ
            GCWMZ(I,K) = GCWMZ(I,KM1) + GDW(I,KM1)*ELADZ
!
            esat  = min(gdpm(i,k), fpvs(gdtm(i,k)))
            GDQSM = min(EPSV*esat/max(gdpm(i,k)+epsvm1*esat, 1.0e-10), 1.0)
            GDHSM = CP*GDTM(I,K ) + GRAV*GDZM(I,K) + EL*GDQSM
!           FDQSM = FDQSAT(GDTM(I,K), GDQSM)
            tem   = one / GDTM(I,K)
            FDQSM = GDQSM * tem * (fact1 + fact2*tem)
            CPGMI = one / (CP + EL*FDQSM)

            PRCZH = PRECZH * MIN(GDZTR(I)/ZTREF, one)
            PRECR = FPREC(GDZM(I,K)-GDZMKB(I), PRCZH )
!
            wrk   = one / GCYM(I,K)
            DCTM        = (GCHMZ(I,K)*wrk - GDHSM) * CPGMI
            GCQMZ(i,k)  = (GDQSM+FDQSM*DCTM) * GCYM(I,K)
            GCQMZ(i,k)  = MIN(GCQMZ(i,k), GCWMZ(I,K))
            GTPRMZ(I,K) = PRECR * (GCWMZ(I,K)-GCQMZ(i,k))
            GTPRMZ(I,K) = MAX(GTPRMZ(I,K), GTPRMZ(I,KM1))
            GCCMZ       = GCWMZ(I,K) - GCQMZ(i,k) - GTPRMZ(I,K )
            DELC        = MIN(GCCMZ, zero)
            GCCMZ       = GCCMZ - DELC
            GCQMZ(i,k)  = GCQMZ(i,k) + DELC
!
            FICE          = FRICE(GDTM(I,K)+DCTM )
            GCIMZ(I,K)    = FICE*GCCMZ
            GSNWIZ(I,KM1) = FICE * (GTPRMZ(I,K)-GTPRMZ(I,KM1))
            GCHMZ(I,K)    = GCHMZ(I,K) + EMELT * (GCIMZ(I,K) + GSNWIZ(I,KM1) &
                                       - GCIMZ(I,KM1) - GDQI(I,KM1)*ELADZ)
            DCTM          = (GCHMZ(I,K)*wrk - GDHSM) * CPGMI
!
            GDQM          = half * (GDQ(I,K,1)     + GDQ(I,KM1,1))
            GDCM          = half * (GDQ(I,K,ITL)   + GDQI(I,K)          &
                                 +  GDQ(I,KM1,ITL) + GDQI(I,KM1))
            GCQM(I,K)     = GCQMZ(i,k)*wrk
            GCCM          = GCCMZ*wrk
!
            BUOYM(I,K)    = (DCTM/GDTM(I,K)                              &
                          + EPSVT*(GCQM(I,K)-GDQM )-GCCM+GDCM) * GRAV
            BUOY(I,KM1)   = half * (BUOYM(I,K)+BUOYM(I,KM1))
!
!DD#ifdef OPT_ASMODE
!DD         WCM(I,K ) &
!DD          = (WCM_(I) + 2.D0*PA*DELZ*BUOY(I,KM1) ) &
!DD               / (one + 2.D0*PB*DELZ*ERMR)
!DD#elif OPT_NS02
!DD         WCM(I,K ) = WCM_(I ) &
!DD           + 2.D0*DELZ*(PA*BUOYM(I,KM1)-ELMD*WCV(I,KM1))
!DD         WCM(I,K ) = MAX(WCM(I,K ), zero )
!DD         WCVX = SQRT(half*(WCM(I,K )+WCM_(I)))
!DD         WCM(I,K) = WCM_(I)  + 2.D0*DELZ*(PA*BUOY(I,KM1)-ELMD*WCVX)
!DD#else
            IF (BUOY(I,KM1) > zero) THEN
              WCM(I,K) = (WCM_(I) + CLMP*DELZ*BUOY(I,KM1)) &
                       / (one + DELZ/TAUZ)
            ELSE
              WCM(I,K) = (WCM_(I) + PA*(DELZ+DELZ)*BUOY(I,KM1) ) &
                       / (one + DELZ/TAUZ + (DELZ+DELZ)*ELAMIN )
            ENDIF
!DD#endif
!
!DD#ifdef OPT_ASMODE
!DD         ELARM2 = ERMR
!DD#elif OPT_NS02
!DD         ELARM2 = ELMD/SQRT(MAX(WCM(I,K), EPSln))
!DD#else
            ELARM2        = CLMD*PA*BUOYM(I,K) / MAX(WCM(I,K), EPSln)
!DD#endif
            ELARM2        = MIN(MAX(ELARM2, ELAMIN), ELAMAX)
            ELAR(I,KM1)   = half * (ELARM1(I) + ELARM2)
            GCYM(I,K)     = GCYM(I,KM1) * EXP(ELAR(I,KM1)*DELZ)
            ELADZ         = GCYM(I,K) - GCYM(I,KM1)
            ELAM(I,KM1)   = ELADZ / DELZ
!
            GCHMZ(I,K)    = GCHMZ(I,KM1) + GDH(I,KM1)*ELADZ
            GCWMZ(I,K)    = GCWMZ(I,KM1) + GDW(I,KM1)*ELADZ
            GCUMZ(I,K)    = GCUMZ(I,KM1) + GDU(I,KM1)*ELADZ
            GCVMZ(I,K)    = GCVMZ(I,KM1) + GDV(I,KM1)*ELADZ
!
            wrk           = one / GCYM(I,K)
            DCTM          = (GCHMZ(I,K)*wrk - GDHSM) * CPGMI
            GCQMZ(i,k)    = (GDQSM+FDQSM*DCTM) * GCYM(I,K)
            GCQMZ(i,k)    = MIN(GCQMZ(i,k), GCWMZ(I,K))
            GTPRMZ(I,K)   = PRECR * (GCWMZ(I,K)-GCQMZ(i,k))
            GTPRMZ(I,K)   = MAX(GTPRMZ(I,K), GTPRMZ(I,KM1))
            GCCMZ         = GCWMZ(I,K) - GCQMZ(i,k) - GTPRMZ(I,K)
            DELC          = MIN(GCCMZ, zero)
            GCCMZ         = GCCMZ - DELC
            GCQMZ(i,k)    = GCQMZ(i,k) + DELC
            GCCM          = GCCMZ*wrk
            GCQM(I,K)     = GCQMZ(i,k)*wrk
!
            FICE          = FRICE(GDTM(I,K)+DCTM )
            GCIMZ(I,K)    = FICE*GCCMZ
            GCIM(I,K)     = GCIMZ(I,K)*wrk
            GCLM(I,K)     = MAX(GCCM-GCIM(I,K), zero)
            GTPRIZ        = GTPRMZ(I,K) - GTPRMZ(I,KM1)
            GSNWIZ(I,KM1) = FICE*GTPRIZ
            GPRCIZ(I,KM1) = (one-FICE )*GTPRIZ
            GCHMZ(I,K)    = GCHMZ(I,K) + EMELT*(GCIMZ(I,K) + GSNWIZ(I,KM1) &
                                        - GCIMZ(I,KM1) - GDQI(I,KM1)*ELADZ )
            GCHM(I,K)     = GCHMZ(I,K)*wrk
            DCTM          = (GCHM(I,K)-GDHSM) * CPGMI
!           GCTM(I,K)     = dctm + GDTM(I,K) + gocp*gdzm(i,k)          ! old, make dse
            GCTM(I,K)     = dctm + (GCHM(I,K) - el*gcqm(i,k)) * oneocp ! new, make dse
!
            GCWM(I,K)     = GCWMZ(I,K)*wrk
            GCUM(I,K)     = GCUMZ(I,K)*wrk
            GCVM(I,K)     = GCVMZ(I,K)*wrk
            DELZL         = GDZ(I,KM1)-GDZM(I,KM1)
            GCY (I,KM1)   = GCYM(I,KM1) * EXP(ELAR(I,KM1)*DELZL)
            GCLZ(I,KM1)   = half * (GCLM(I,K) + GCLM(I,KM1)) * GCY(I,KM1)
            GCIZ(I,KM1)   = half * (GCIM(I,K) + GCIM(I,KM1)) * GCY(I,KM1)
            IF (BUOY(I,KM1) > zero) THEN
              ACWF(I)     = ACWF(I) + BUOY(I,KM1)*GCY(I,KM1)*DELZ
            ENDIF
!
            ELARM1(I)     = ELARM2
            WCM_(I)       = WCM(I,K)

            rhs_h = cbmfx(i)*(gchmz(i,k) - (gchmz(i,km1) + GDH(I,KM1)*ELADZ))
            rhs_q = cbmfx(i)*(gcwmz(i,k)-gcqmz(i,k)                   &
                           - (gcwmz(i,km1)-gcqmz(i,km1)               &
                           + (GDw(I,KM1)-gdq(i,km1,1))*ELADZ))
            dqcond(i,km1) = -rhs_q
            dtfrz(i,km1)  =  rhs_h*oneocp
            dtcond(i,km1) = -ELocp*DQCOND(i,km1)

          ENDIF ! IF (K > KB(I) .AND. WCM_(I) > WCCRT) THEN
        ENDDO
      ENDDO
!
!     < find cloud top >
!
      DO I=ISTS,IENS
        KT(I) = -1
      enddo
      DO K=KMAX,2,-1
        DO I=ISTS,IENS
          IF (K > KB(I) .AND. KT(I)  == -1                            &
              .AND. BUOYM(I,K) > zero .AND. WCM(I,K) > WCCRT) THEN
            KT(I) = K
          ENDIF
        ENDDO
      ENDDO
!
      KTMX = 2
      DO I=ISTS,IENS
        kt(i) = min(kt(i), kmax-1)
        KTMX  = max(ktmx, KT(I))
      ENDDO
!
      DO I=ISTS,IENS
        kk = kt(i)
        IF (KK > 0 ) then
          do k=kk+1,kmax
            GCYM(I,K) = zero
          enddo
          do k=kk,kmax
            GCLZ  (I,K) = zero 
            GCIZ  (I,K) = zero
            GPRCIZ(I,K) = zero
            GSNWIZ(I,K) = zero
            dtcond(i,k) = 0.0
            dqcond(i,k) = 0.0
            dtfrz(i,k)  = 0.0
          enddo
        ELSE
          do k=kb(i)+1,kmax
            GCYM(I,K) = zero
          enddo
          do k=1,kmax
            GCLZ  (I,k) = zero
            GCIZ  (I,k) = zero
            GPRCIZ(I,k) = zero
            GSNWIZ(I,k) = zero
            dtcond(i,k) = 0.0
            dqcond(i,k) = 0.0
            dtfrz(i,k)  = 0.0
          enddo
        ENDIF
      ENDDO
!
!     < cloud top properties >
!
      DO I=ISTS,IENS
        IF (KT(I) > 0) THEN
          K   = KT(I)
          kp1 = k + 1
          GCYT(I)     = GCY(I,K)
          ELADZ       = GCYT(I) - GCYM(I,K)
          ELAM(I,K)   = ELADZ / (GDZ(I,K)-GDZM(I,K))
!
          GCHT(I)     = GCHMZ(I,K) + GDH(I,K)*ELADZ
          GCWT(i)     = GCWMZ(I,K) + GDW(I,K)*ELADZ
          GCUT(I)     = GCUMZ(I,K) + GDU(I,K)*ELADZ
          GCVT(I)     = GCVMZ(I,K) + GDV(I,K)*ELADZ
!
          DCT         = (GCHT(I)/GCYT(I) - GDHS(I,K)) &
                      / (CP*(one + GAM(I,K)))
          GCQT(I)     = (GDQS(I,K) + FDQS(I,K)*DCT) * GCYT(I)
          GCQT(I)     = MIN(GCQT(I), GCWT(i))
          PRCZH       = PRECZH * MIN(GDZTR(I)/ZTREF, one)
          GTPRT(I)    = FPREC(GDZ(I,K)-GDZMKB(I), PRCZH) * (GCWT(i)-GCQT(I))
          GTPRT(I)    = MAX(GTPRT(I), GTPRMZ(I,K))
          GCCT        = GCWT(i) - GCQT(I) - GTPRT(I)
          DELC        = MIN(GCCT, zero)
          GCCT        = GCCT - DELC
          GCQT(I)     = GCQT(I) + DELC
!
          FICE        = FRICE(GDT(I,K)+DCT)
          GCIT(I)     = FICE*GCCT
          GCLT(I)     = (one-FICE) * GCCT
          GTPRIZ      = GTPRT(I) - GTPRMZ(I,K)
          GPRCIZ(I,K) = (one-FICE) * GTPRIZ
          GSNWIZ(I,K) = FICE * GTPRIZ
          GCHT(I)     = GCHT(I)                                          &
                      + EMELT * (GCIT(I) + GSNWIZ(I,K) - GCIMZ(I,K) - GDQI(I,K)*ELADZ)
!
          GCUT(I)     = GCUT(I)*(one-CPRES) + GCY(I,K)*GDU(I,K)*CPRES
          GCVT(I)     = GCVT(I)*(one-CPRES) + GCY(I,K)*GDV(I,K)*CPRES
          GCLZ(I,K)   = GCLT(I)
          GCIZ(I,K)   = GCIT(I)

!DD AW get  the cloud top values denormalized and put into profile
          mygcht      = gcht(I) - el*(gcwt(i) - gcqt(i))

          wrk         = one / gcyt(i)
          gctm(i,kp1) = wrk * (mygcht - el*gcqt(i)) * oneocp
!Moorthi  gcqm(i,kp1) = gcqt(i)
          gcqm(i,kp1) = gcqt(i)*wrk     ! check this - oct17 2016
          gcim(i,kp1) = gcit(i)*wrk
          gclm(i,kp1) = gclt(i)*wrk
!
          rhs_q = cbmfx(i)*( gcwt(i)-gcqt(i) - (gcwmz(i,k)-gcqmz(i,k) &
                                      + (GDw(I,K)-gdq(i,k,1))*ELADZ) )
          dqcond(i,k) = -rhs_q
          rhs_h       = cbmfx(i)*(gcht(i) - (gchmz(i,k) + GDH(I,K)*ELADZ))

          dtfrz(i,k)  = rhs_h * oneocp
          dtcond(i,k) = -ELocp*DQCOND(i,k)
        ENDIF
      ENDDO
!
!DD#ifdef OPT_CUMBGT   /* budget check */
!DD   HBGT ( ISTS:IENS ) = 0.D0
!DD   WBGT ( ISTS:IENS ) = 0.D0
!DD   PBGT ( ISTS:IENS ) = 0.D0
!DD   MBGT ( ISTS:IENS ) = 0.D0
!DD   GTPRX( ISTS:IENS ) = 0.D0
!DD   GSNWT( ISTS:IENS ) = 0.D0
!DD!
!DD   IF ( CTP .EQ. 1 ) THEN
!DD     HBMX = 0.D0
!DD     WBMX = 0.D0
!DD     PBMX = 0.D0
!DD     MBMX = 0.D0
!DD   END IF
!DD!
!DD   DO K = 2, KMAX
!DD     DO I = ISTS, IENS
!DD       IF ( K .GE. KB( I ) .AND. K .LT. KT( I ) ) THEN
!DD         ELADZ = GCYM( I,K+1 ) - GCYM( I,K )
!DD         DELZ  = GDZM( I,K+1 ) - GDZM( I,K )
!DD         HBGT( I ) = HBGT( I ) + ( GDH( I,K )-EMELT*GDQI( I,K ) )*ELADZ
!DD         WBGT( I ) = WBGT( I ) + GDW( I,K )*ELADZ
!DD         MBGT( I ) = MBGT( I ) + ELAM( I,K )*DELZ
!DD         GTPRX( I ) = GTPRX( I ) + GPRCIZ( I,K ) + GSNWIZ( I,K )
!DD         GSNWT( I ) = GSNWT( I ) + GSNWIZ( I,K )
!DD       END IF
!DD     END DO
!DD   END DO
!DD!
!DD   DO I = ISTS, IENS
!DD     IF ( KT( I ) .GT. KB( I ) ) THEN
!DD       ELADZ = GCYT( I ) - GCYM( I,KT(I) )
!DD       DELZ  = GDZ( I,KT(I) )-GDZM( I,KT(I) )
!DD       GTPRX( I ) = GTPRX( I ) + GPRCIZ( I,KT(I) ) + GSNWIZ( I,KT(I) )
!DD       GSNWT( I ) = GSNWT( I ) + GSNWIZ( I,KT(I) )
!DD       HBGT( I )  = HBGT( I )  + GCHB( I ) - EMELT*GCIB( I )        &
!DD                  + ( GDH( I,KT(I) )-EMELT*GDQI( I,KT(I) ) ) *ELADZ &
!DD                  - ( GCHT(I)-EMELT*( GCIT(I)+GSNWT(I) ) )
!DD       WBGT( I ) = WBGT( I ) &
!DD                 + GCWB( I ) + GDW( I,KT(I) )*ELADZ &
!DD                 - GCQT( I ) - GCLT( I ) - GCIT( I ) &
!DD                 - GTPRT( I )
!DD       MBGT( I ) = MBGT( I ) + one + ELAM( I,KT(I) )*DELZ &
!DD                 - GCYT( I )
!DD       PBGT( I ) = GTPRT( I ) - GTPRX( I )
!DD!
!DD       IF ( ABS( HBGT(I) ) .GT. ABS( HBMX ) ) HBMX = HBGT(I)
!DD       IF ( ABS( WBGT(I) ) .GT. ABS( WBMX ) ) WBMX = WBGT(I)
!DD       IF ( ABS( PBGT(I) ) .GT. ABS( PBMX ) ) PBMX = PBGT(I)
!DD       IF ( ABS( MBGT(I) ) .GT. ABS( MBMX ) ) MBMX = MBGT(I)
!DD     END IF
!DD   END DO
!DD!
!DD   IF ( CTP .EQ. NCTP ) THEN
!DD     WRITE( iulog,* ) '### CUMUP(rank=',irank,'): energy imbalance =', HBMX
!DD     WRITE( iulog,* ) '### CUMUP(rank=',irank,'): water imbalance =', WBMX
!DD     WRITE( iulog,* ) '### CUMUP(rank=',irank,'): precipitation imbalance =', PBMX
!DD     WRITE( iulog,* ) '### CUMUP(rank=',irank,'): mass imbalance =', MBMX
!DD   END IF
!DD#endif
!
!      WRITE( CTNUM, '(I2.2)' ) CTP
!
      END SUBROUTINE CUMUP
!***********************************************************************
      SUBROUTINE CUMBMX   & !! cloud base mass flux
               ( IJSDIM, KMAX  ,           & !DD dimensions
                 CBMFX ,                   & ! modified
                 ACWF  , GCYT  , GDZM  ,   & ! input
                 GDW   , GDQS  , DELP  ,   & ! input
                 KT    , KTMX  , KB    ,   & ! input
                 DELT  , ISTS  , IENS    )   ! input
!
!
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: IJSDIM, KMAX  ! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(r8)     CBMFX (IJSDIM)          ! cloud base mass flux
!
!   [INPUT]
      REAL(r8)     ACWF  (IJSDIM)          ! cloud work function
      REAL(r8)     GCYT  (IJSDIM)          ! norm mass flux @top
      REAL(r8)     GDZM  (IJSDIM, KMAX+1)  ! height
      REAL(r8)     GDW   (IJSDIM, KMAX)    ! total water
      REAL(r8)     GDQS  (IJSDIM, KMAX)    ! saturate humidity
      REAL(r8)     DELP  (IJSDIM, KMAX)    ! delt pressure
      INTEGER      KT    (IJSDIM)          ! cloud top
      INTEGER      KTMX                    ! max. of cloud top
      INTEGER      KB    (IJSDIM)          ! cloud base
      REAL(r8)     DELT                    ! time step
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8), dimension(ijsdim) :: QX, QSX, RHM
      INTEGER      I, K
      REAL(r8)     ALP, FMAX1, wrk
!
!   [INTERNAL PARAM]
      REAL(r8) :: FMAX   = 1.5e-2_r8         ! maximum flux
      REAL(r8) :: RHMCRT = zero              ! critical val. of RH@ all could
!     REAL(r8) :: RHMCRT = 0.5_r8            ! critical val. of RH@ all could
      REAL(r8) :: ALP1   = zero
      REAL(r8) :: TAUD   = 1.e3_r8
      REAL(r8) :: ZFMAX  = 3.5e3_r8
      REAL(r8) :: ZDFMAX = 5.e2_r8
!     REAL(r8) :: FMAXP  = 2._r8
      REAL(r8) :: EPSln  = 1.e-10_r8
!
      do i=ists,iens
        qx(i)  = zero
        qsx(i) = zero
      enddo
!
      DO K=1,KTMX
        DO I=ISTS,IENS
          IF (K >= KB(I) .AND. K <= KT(I)) THEN
            QX (I) = QX (I) + GDW (I,K) * DELP(I,K)
            QSX(I) = QSX(I) + GDQS(I,K) * DELP(I,K)
          ENDIF
        ENDDO
      ENDDO
      DO I=ISTS,IENS
        RHM(I) = min(one, max(zero, QX(I)/MAX(QSX(I),EPSln)))
      ENDDO
!
      wrk = one + delt/(taud+taud)
      DO I=ISTS,IENS
        IF (KT(I) > KB(I) .AND. RHM(I) >= RHMCRT) THEN
          ALP      = ALP0 + ALP1 * (GDZM(I,KT(I))-GDZM(I,KB(I)))
          FMAX1    = (one - TANH((GDZM(I,1)-ZFMAX)/ZDFMAX)) * half
!         FMAX1    = FMAX * FMAX1**FMAXP
          FMAX1    = FMAX * FMAX1 * FMAX1
!         CBMFX(I) = CBMFX(I) + MAX(ACWF(I), zero)/(ALP+ALP)*DELT
!         CBMFX(I) = CBMFX(I) / (one + DELT/(TAUD+TAUD))
          CBMFX(I) = (CBMFX(I) + MAX(ACWF(I), zero)/(ALP+ALP)*DELT) * wrk
          CBMFX(I) = MIN(max(CBMFX(I), zero), FMAX1/GCYT(I))
        ELSE
          CBMFX(I) = zero
        ENDIF
      ENDDO
!
      END SUBROUTINE CUMBMX
!***********************************************************************
      SUBROUTINE CUMFLX                                   & !! cloud mass flux
                      ( IM    , IJSDIM, KMAX  ,           & !DD dimensions
                        GMFLX , GPRCI , GSNWI ,           & ! output
                        QLIQ  , QICE  , GTPRC0,           & ! output
                        CBMFX , GCYM  , GPRCIZ, GSNWIZ,   & ! input
                        GTPRT , GCLZ  , GCIZ  ,           & ! input
                        KB    , KT    , KTMX  ,           & ! input
                        ISTS  , IENS  , sigma           )   ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, KMAX, IM            !! DD, for GFS, pass in
!
!   [OUTPUT]
      REAL(r8)     GMFLX (IJSDIM, KMAX)     !! mass flux
      REAL(r8)     GPRCI (IJSDIM, KMAX)     !! rainfall generation
      REAL(r8)     GSNWI (IJSDIM, KMAX)     !! snowfall generation
      REAL(r8)     QLIQ  (IJSDIM, KMAX)     !! cloud liquid
      REAL(r8)     QICE  (IJSDIM, KMAX)     !! cloud ice
      REAL(r8)     GTPRC0(IJSDIM)           !! precip. before evap.
!
!   [INPUT]
      REAL(r8)     CBMFX (IJSDIM)           !! cloud base mass flux
      REAL(r8)     GCYM  (IJSDIM, KMAX)     !! normalized mass flux
      REAL(r8)     sigma (IJSDIM, KMAX)     !! AW sigma
      REAL(r8)     GPRCIZ(IJSDIM, KMAX)     !! precipitation/M
      REAL(r8)     GSNWIZ(IJSDIM, KMAX)     !! snowfall/M
      REAL(r8)     GTPRT (IJSDIM)           !! rain+snow @top
      REAL(r8)     GCLZ  (IJSDIM, KMAX)     !! cloud liquid/M
      REAL(r8)     GCIZ  (IJSDIM, KMAX)     !! cloud ice/M
      real(r8)     tem
      INTEGER      KB    (IJSDIM)           !! cloud base
      INTEGER      KT    (IJSDIM)           !! cloud top
      INTEGER      KTMX                     !! max of cloud top
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER    I, K
!
!M    DO K=1,KTMX
!M      DO I=ISTS,IENS
!M        GMFLX(I,K) = GMFLX(I,K) + GCYM(I,K)*CBMFX(I)
!M      ENDDO
!M    ENDDO
!
      DO K=1,KTMX
        DO I=ISTS,IENS
          tem = CBMFX(I) * (one - sigma(i,k))
          GMFLX(I,K) = GMFLX(I,K) + tem * GCYM(I,K)
          GPRCI(I,K) = GPRCI(I,K) + tem * GPRCIZ(I,K)
          GSNWI(I,K) = GSNWI(I,K) + tem * GSNWIZ(I,K)
          QLIQ(I,K)  = QLIQ (I,K) + tem * GCLZ(I,K)
          QICE(I,K)  = QICE (I,K) + tem * GCIZ(I,K)

!         GMFLX(I,K) = GMFLX(I,K) + GCYM(I,K)   * CBMFX(I)
!         GPRCI(I,K) = GPRCI(I,K) + GPRCIZ(I,K) * CBMFX(I)
!         GSNWI(I,K) = GSNWI(I,K) + GSNWIZ(I,K) * CBMFX(I)
!         QLIQ(I,K)  = QLIQ (I,K) + GCLZ(I,K)   * CBMFX(I)
!         QICE(I,K)  = QICE (I,K) + GCIZ(I,K)   * CBMFX(I)
        ENDDO
      ENDDO
!
      DO I= ISTS,IENS
        GTPRC0(I) = GTPRC0(I) + GTPRT(I) * CBMFX(I)
      ENDDO
!
!M    DO K = 1, KTMX
!M      DO I = ISTS, IENS
!M        QLIQ(I,K) = QLIQ(I,K) + GCLZ(I,K)*CBMFX(I)
!M        QICE(I,K) = QICE(I,K) + GCIZ(I,K)*CBMFX(I)
!M      ENDDO
!M    ENDDO
!
      END SUBROUTINE CUMFLX
!***********************************************************************
      SUBROUTINE CUMDET                                    & !! detrainment
               ( im    , IJSDIM, KMAX  , NTR   ,           & !DD dimensions
                 CMDET ,                                   & ! output
!                CMDET , GTLDET, GTIDET,                   & ! output
                 GTT   , GTQ   , GTCFRC, GTU   , GTV   ,   & ! modified
!                GTQI  ,                                   & ! modified
                 GDH   , GDQ   , GDCFRC, GDU   , GDV   ,   & ! input
                 CBMFX , GCYT  , DELP  , GCHT  , GCQT  ,   & ! input
                 GCLT  , GCIT  , GCUT  , GCVT  , GDQI  ,   & ! input
                 KT    , ISTS  , IENS  , nctp  , sigi    )   ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: im, IJSDIM, KMAX, NTR, nctp      !! DD, for GFS, pass in
!
!   [OUTPUT]
      REAL(r8)     CMDET (IJSDIM, KMAX)   !! detrainment mass flux
!     REAL(r8)     GTLDET(IJSDIM, KMAX)   !! cloud liquid tendency by detrainment
!     REAL(r8)     GTIDET(IJSDIM, KMAX)   !! cloud ice tendency by detrainment
!
!   [MODIFY]
      REAL(r8)     GTT   (IJSDIM, KMAX)   !! temperature tendency
      REAL(r8)     GTQ   (IJSDIM, KMAX, NTR)   !! moisture tendency
      REAL(r8)     GTCFRC(IJSDIM, KMAX)   !! cloud fraction tendency
      REAL(r8)     GTU   (IJSDIM, KMAX)   !! u tendency
      REAL(r8)     GTV   (IJSDIM, KMAX)   !! v tendency
!     REAL(r8)     GTQI  (IJSDIM, KMAX)   !! cloud ice tendency
!
!   [INPUT]
      REAL(r8)     GDH   (IJSDIM, KMAX)   !! moist static energy
      REAL(r8)     GDQ   (IJSDIM, KMAX, NTR) !! humidity qv
      REAL(r8)     GDCFRC(IJSDIM, KMAX)   !! cloud fraction
      REAL(r8)     GDU   (IJSDIM, KMAX)
      REAL(r8)     GDV   (IJSDIM, KMAX)
      REAL(r8)     DELP  (IJSDIM, KMAX)
      REAL(r8)     CBMFX (IM,     NCTP)   !! cloud base mass flux
      REAL(r8)     GCYT  (IJSDIM, NCTP)   !! detraining mass flux
      REAL(r8)     GCHT  (IJSDIM, NCTP)   !! detraining MSE
      REAL(r8)     GCQT  (IJSDIM, NCTP)   !! detraining qv
      REAL(r8)     GCLT  (IJSDIM, NCTP)   !! detraining ql
      REAL(r8)     GCIT  (IJSDIM, NCTP)   !! detraining qi
      REAL(r8)     GCUT  (IJSDIM, NCTP)   !! detraining u
      REAL(r8)     GCVT  (IJSDIM, NCTP)   !! detraining v
      REAL(r8)     GDQI  (IJSDIM, KMAX)   !! cloud ice
      REAL(r8)     sigi  (IJSDIM, KMAX,nctp) !! cloud fraction
      INTEGER      KT    (IJSDIM, NCTP)   !! cloud top
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     sigma(ijsdim)
      REAL(r8)     GTHCI, GTQVCI, GTQLCI, GTQICI, GTXCI, tem
!M    REAL(r8)     GTCCI
!M    REAL(r8)     GTUCI, GTVCI
      INTEGER      I, K, CTP, kk
!
!
!PARALLEL_FORBID

      do k=1,kmax
        DO I=ISTS,IENS
          CMDET (I,k) = zero
!         GTLDET(I,k) = zero
!         GTIDET(I,k) = zero
        enddo
      enddo
      do i=ists,iens
        sigma(i) = zero
      enddo

!PARALLEL_FORBID
      DO CTP=1,NCTP
        DO I=ISTS,IENS
          K = KT(I,CTP)
          IF (K > 0) THEN
            sigma(i) = sigma(i) + sigi(i,k,ctp)
            tem    = CBMFX(I,CTP) * (one - sigma(i))
            GTXCI  = GRAV/DELP(I,K)*tem

            GTHCI  = GTXCI * (GCHT(I,CTP) - GCYT(I,CTP)*GDH(I,K))
            GTQVCI = GTXCI * (GCQT(I,CTP) - GCYT(I,CTP)*GDQ(I,K,1))
            GTQLCI = GTXCI * (GCLT(I,CTP) - GCYT(I,CTP)*GDQ(I,K,ITL))
            GTQICI = GTXCI * (GCIT(I,CTP) - GCYT(I,CTP)*GDQI(I,K))
!
            GTQ(I,K,1)   = GTQ(I,K,1) + GTQVCI
            GTT(I,K)     = GTT(I,K) + (GTHCI - EL*GTQVCI) * oneocp
! ql tendency by detrainment is treated by stratiform scheme
            GTQ(I,K,ITL) = GTQ(I,K,ITL) + GTQLCI
!           GTLDET(I,K)  = GTLDET(I,K) + GTQLCI
! qi tendency by detrainment is treated by stratiform scheme
!           GTQI  (I,K)  = GTQI(I,K)   + GTQICI
!           GTIDET(I,K)  = GTIDET(I,K) + GTQICI
            GTQ(I,K,ITI) = GTQ(I,K,ITI) + GTQICI

            GTCFRC(I,K)  = GTCFRC(I,K) + GTXCI * (GCYT(I,CTP) - GCYT(I,CTP)*GDCFRC(I,K))
            GTU(I,K)     = GTU(I,K)    + GTXCI * (GCUT(I,CTP) - GCYT(I,CTP)*GDU(I,K))
            GTV(I,K)     = GTV(I,K)    + GTXCI * (GCVT(I,CTP) - GCYT(I,CTP)*GDV(I,K))
!
            CMDET(I,K )  = CMDET(I,K) + GCYT(I,CTP) * tem
          ENDIF
        ENDDO
      ENDDO
!
      END SUBROUTINE CUMDET
!***********************************************************************
      SUBROUTINE CUMSBH                         & !! adiabat. descent
               ( IM    , IJSDIM, KMAX  , NTR   ,& !DD dimensions
                 GTT   , GTQ   ,                & ! modified
!                GTT   , GTQ   , GTQI  ,        & ! modified
                 GTU   , GTV   ,                & ! modified
                 GDH   , GDQ   , GDQI  ,        & ! input
                 GDU   , GDV   ,                & ! input
                 DELP  , GMFLX , GMFX0 ,        & ! input
                 KTMX  , CPRES , ISTS  , IENS )   ! input
!
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, IM, KMAX, NTR             !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(r8)     GTT   (IJSDIM, KMAX)      !! Temperature tendency
      REAL(r8)     GTQ   (IJSDIM, KMAX, NTR) !! Moisture etc tendency
!     REAL(r8)     GTQI  (IJSDIM, KMAX)
      REAL(r8)     GTU   (IJSDIM, KMAX)      !! u tendency
      REAL(r8)     GTV   (IJSDIM, KMAX)      !! v tendency
!
!   [INPUT]
      REAL(r8)     GDH   (IJSDIM, KMAX)
      REAL(r8)     GDQ   (IJSDIM, KMAX, NTR) !! humidity etc
      REAL(r8)     GDQI  (IJSDIM, KMAX)
      REAL(r8)     GDU   (IJSDIM, KMAX)
      REAL(r8)     GDV   (IJSDIM, KMAX)
      REAL(r8)     DELP  (IJSDIM, KMAX)
      REAL(r8)     GMFLX (IJSDIM, KMAX)      !! mass flux (updraft+downdraft)
      REAL(r8)     GMFX0 (IJSDIM, KMAX)      !! mass flux (updraft only)
      INTEGER      KTMX
      REAL(r8)     CPRES                     !! pressure factor for cumulus friction
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER      I, K, KM, KP
      REAL(r8)     SBH0, SBQ0, SBL0, SBI0, SBC0,  SBS0,        &
                   SBH1, SBQ1, SBL1, SBI1, SBC1,  SBS1,   FX1, &
                   SBU0, SBV0, SBU1, SBV1, GTHCI, GTQVCI,      &
                   GTQLCI, GTQICI, GTM2CI, GTM3CI, wrk, wrk1
!M    REAL(r8)     GTUCI, GTVCI, wrk, wrk1
      REAL(r8)     FX(ISTS:IENS)

      REAL(r8), dimension(IJSDIM, KMAX)  :: GTLSBH, GTISBH
!
!
      FX     = zero
      do k=1,kmax
        do i=ists,iens
          GTLSBH(i,k) = zero
          GTISBH(i,k) = zero
        enddo
      enddo
!
      DO K=KTMX,1,-1
        KM = MAX(K-1, 1)
        KP = MIN(K+1, KMAX)
        DO I=ISTS,IENS
          SBH0 = GMFLX(I,KP) * (GDH(I,KP)-GDH(I,K))
          SBQ0 = GMFLX(I,KP) * (GDQ(I,KP,1)-GDQ(I,K,1))
          SBL0 = GMFLX(I,KP) * (GDQ(I,KP,ITL )-GDQ(I,K,ITL))
          SBI0 = GMFLX(I,KP) * (GDQI(I,KP)-GDQI(I,K))
          SBU0 = GMFLX(I,KP) * (GDU(I,KP)-GDU(I,K))           &
               - GMFX0(I,KP) * (GDU(I,KP)-GDU(I,K))*CPRES
          SBV0 = GMFLX(I,KP) * (GDV(I,KP)-GDV(I,K))           &
               - GMFX0(I,KP) * (GDV(I,KP)-GDV(I,K))*CPRES
!
          SBH1 = GMFLX(I,K) * (GDH(I,K)-GDH(I,KM))
          SBQ1 = GMFLX(I,K) * (GDQ(I,K,1)-GDQ(I,KM,1))
          SBL1 = GMFLX(I,K) * (GDQ(I,K,ITL)-GDQ(I,KM,ITL))
          SBI1 = GMFLX(I,K) * (GDQI(I,K)-GDQI(I,KM))
          SBU1 = GMFLX(I,K) * (GDU(I,K)-GDU(I,KM))            &
               - GMFX0(I,K) * (GDU(I,K)-GDU(I,KM))*CPRES
          SBV1 = GMFLX(I,K) * (GDV(I,K)-GDV(I,KM))            &
               - GMFX0(I,K) * (GDV(I,K)-GDV(I,KM))*CPRES
!
!#ifndef SYS_SX   /* original */
          IF (GMFLX(I,K) > GMFLX(I,KP)) THEN
             FX1 = half
          ELSE
             FX1 = zero
          ENDIF
!#else            /* optimized for NEC SX series */
!         FX1 = 0.25D0 - SIGN(0.25D0,GMFLX(I,K+1)-GMFLX(I,K)) !! 0.5 or 0.
!#endif
!
          wrk    = GRAV / DELP(I,K)
          wrk1   = one - FX(I)
          GTHCI  = wrk * (wrk1*SBH0 + FX1 *SBH1)
          GTQVCI = wrk * (wrk1*SBQ0 + FX1 *SBQ1)
          GTQLCI = wrk * (wrk1*SBL0 + FX1 *SBL1)
          GTQICI = wrk * (wrk1*SBI0 + FX1 *SBI1)
!M        GTUCI  = wrk * (wrk1*SBU0 + FX1 *SBU1)
!M        GTVCI  = wrk * (wrk1*SBV0 + FX1 *SBV1)
!
          GTT (I,K    ) = GTT(I,K)     + (GTHCI-EL*GTQVCI)*oneocp
          GTQ (I,K,1  ) = GTQ(I,K,1)   +  GTQVCI
          GTQ (I,K,ITL) = GTQ(I,K,ITL) +  GTQLCI
          GTQ (I,K,ITI) = GTQ(I,K,ITI) +  GTQICI
!         GTQI(I,K)     = GTQI(I,K)    +  GTQICI
!M        GTU (I,K)     = GTU(I,K)     +  GTUCI
!M        GTV (I,K)     = GTV(I,K)     +  GTVCI
          GTU (I,K)     = GTU(I,K)     +  wrk * (wrk1*SBU0 + FX1*SBU1)
          GTV (I,K)     = GTV(I,K)     +  wrk * (wrk1*SBV0 + FX1*SBV1)

          GTLSBH(I,K)   = GTQLCI
          GTISBH(I,K)   = GTQICI
!
!         SBC0 = GMFLX(I,K+1) * (GDQ(I,KP,IMU2)-GDQ(I,K,IMU2))
!         SBS0 = GMFLX(I,K+1) * (GDQ(I,KP,IMU3)-GDQ(I,K,IMU3))
!         SBC1 = GMFLX(I,K)   * (GDQ(I,K,IMU2)-GDQ(I,KM,IMU2))
!         SBS1 = GMFLX(I,K)   * (GDQ(I,K,IMU3)-GDQ(I,KM,IMU3))
!         GTM2CI = GRAV/DELP(I,K)
!     &          *(( one-FX(I))*SBC0 + FX1 *SBC1)
!         GTM3CI = GRAV/DELP(I,K)
!     &          *((one-FX(I))*SBS0 + FX1 *SBS1)
!         GTQ(I,K,IMU2) = GTQ(I,K,IMU2) + GTM2CI
!         GTQ(I,K,IMU3) = GTQ(I,K,IMU3) + GTM3CI
!
          FX(I) = FX1
        enddo
      enddo
!
      END SUBROUTINE CUMSBH
!***********************************************************************
!
      SUBROUTINE CUMDWN                            & ! Freeze & Melt & Evaporation
               ( IM    , IJSDIM, KMAX  , NTR   ,   & !DD dimensions
                 GTT   , GTQ   , GTU   , GTV   ,   & ! modified
                         GMFLX ,                   & ! modified
!                GTQI  , GMFLX ,                   & ! modified
                 GPRCP , GSNWP , GTEVP , GMDD  ,   & ! output
                 GPRCI , GSNWI ,                   & ! input
                 GDH   , GDW   , GDQ   , GDQI  ,   & ! input
                 GDQS  , GDS   , GDHS  , GDT   ,   & ! input
                 GDU   , GDV   , GDZ   ,           & ! input
                 GDZM  , GCYM  , FDQS  , DELP  ,   & ! input
                 sigmad, do_aw , do_awdd       ,   & !DDsigma input
                 gtmelt, gtevap, gtsubl,           & !DDsigma input
                 dtdwn , dqvdwn, dqldwn, dqidwn,   & !DDsigma input
                 KB    , KTMX  , ISTS  , IENS    )   ! input
!
! DD AW : modify to get eddy fluxes and microphysical tendencies for AW
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX, NTR   ! DD, for GFS, pass in
      logical, intent(in) :: do_aw, do_awdd
!
!   [MODIFY]
      REAL(r8)     GTT   (IJSDIM, KMAX)       ! Temperature tendency
      REAL(r8)     GTQ   (IJSDIM, KMAX, NTR)  ! Moisture etc tendency
      REAL(r8)     GTU   (IJSDIM, KMAX)       ! u tendency
      REAL(r8)     GTV   (IJSDIM, KMAX)       ! v tendency
!     REAL(r8)     GTQI  (IJSDIM, KMAX)       ! cloud ice tendency
      REAL(r8)     GMFLX (IJSDIM, KMAX)       ! mass flux
!
!   [OUTPUT]
      REAL(r8)     GPRCP (IJSDIM, KMAX)       ! rainfall flux
      REAL(r8)     GSNWP (IJSDIM, KMAX)       ! snowfall flux
      REAL(r8)     GTEVP (IJSDIM, KMAX)       ! evaporation+sublimation
      REAL(r8)     GMDD  (IJSDIM, KMAX)       ! downdraft mass flux

!AW microphysical tendencies
      REAL(r8)     gtmelt (IJSDIM, KMAX)      ! t tendency ice-liq
      REAL(r8)     gtevap (IJSDIM, KMAX)      ! t tendency liq-vapor
      REAL(r8)     gtsubl (IJSDIM, KMAX)      ! t tendency ice-vapor
!AW eddy flux tendencies
      REAL(r8)     dtdwn  (IJSDIM, KMAX)      ! t tendency downdraft detrainment
      REAL(r8)     dqvdwn (IJSDIM, KMAX)      ! qv tendency downdraft detrainment
      REAL(r8)     dqldwn (IJSDIM, KMAX)      ! ql tendency downdraft detrainment
      REAL(r8)     dqidwn (IJSDIM, KMAX)      ! qi tendency downdraft detrainment
! AW downdraft area fraction (assumed zero for now)
      REAL(r8)     sigmad (IM,KMAX)           !DDsigma cloud downdraft area fraction

!   [INPUT]
      REAL(r8)     GPRCI (IJSDIM, KMAX)     ! rainfall generation
      REAL(r8)     GSNWI (IJSDIM, KMAX)     ! snowfall generation
      REAL(r8)     GDH   (IJSDIM, KMAX)     ! moist static energy
      REAL(r8)     GDW   (IJSDIM, KMAX)     ! total water
      REAL(r8)     GDQ   (IJSDIM, KMAX, NTR)! humidity etc
      REAL(r8)     GDQI  (IJSDIM, KMAX)     ! cloud ice
      REAL(r8)     GDQS  (IJSDIM, KMAX)     ! saturate humidity
      REAL(r8)     GDS   (IJSDIM, KMAX)     ! dry static energy
      REAL(r8)     GDHS  (IJSDIM, KMAX)     ! saturate moist static energy
      REAL(r8)     GDT   (IJSDIM, KMAX)     ! air temperature T
      REAL(r8)     GDU   (IJSDIM, KMAX)     ! u-velocity
      REAL(r8)     GDV   (IJSDIM, KMAX)     ! v-velocity
      REAL(r8)     GDZ   (IJSDIM, KMAX)     ! altitude
      REAL(r8)     GDZM  (IJSDIM, KMAX+1)   ! altitude (half lev)
      REAL(r8)     GCYM  (IJSDIM, KMAX)     ! norm. mass flux
      REAL(r8)     FDQS  (IJSDIM, KMAX)
      REAL(r8)     DELP  (IJSDIM, KMAX)
      INTEGER      KB    (IJSDIM)
      INTEGER      KTMX, ISTS, IENS
!
!   [INTERNAL WORK]
! Note: Some variables have 3-dimensions for the purpose of budget check.
      REAL(r8)     EVAPD (IJSDIM, KMAX)      ! evap. in downdraft
      REAL(r8)     SUBLD (IJSDIM, KMAX)      ! subl. in downdraft
      REAL(r8)     EVAPE (IJSDIM, KMAX)      ! evap. in environment
      REAL(r8)     SUBLE (IJSDIM, KMAX)      ! subl. in environment
      REAL(r8)     EVAPX (IJSDIM, KMAX)      ! evap. env. to DD
      REAL(r8)     SUBLX (IJSDIM, KMAX)      ! subl. env. to DD
      REAL(r8)     GMDDE (IJSDIM, KMAX)      ! downdraft entrainment
      REAL(r8)     SNMLT (IJSDIM, KMAX)      ! melt - freeze
      REAL(r8)     GCHDD (IJSDIM, KMAX)      ! MSE detrainment
      REAL(r8)     GCWDD (IJSDIM, KMAX)      ! water detrainment
      REAL(r8)     GTTEV (IJSDIM, KMAX)      ! T tendency by evaporation
      REAL(r8)     GTQEV (IJSDIM, KMAX)      ! q tendency by evaporation
      REAL(r8)     GCHD  (ISTS:IENS)         ! downdraft MSE
      REAL(r8)     GCWD  (ISTS:IENS)         ! downdraft q
! profiles of downdraft variables for AW flux tendencies
      REAL(r8)     GCdseD(ISTS:IENS, KMAX)   ! downdraft dse
      REAL(r8)     GCqvD (ISTS:IENS, KMAX)   ! downdraft qv
      REAL(r8)     GCqlD (ISTS:IENS, KMAX)   ! downdraft ql
      REAL(r8)     GCqiD (ISTS:IENS, KMAX)   ! downdraft qi

      REAL(r8)     GCUD  (ISTS:IENS)         ! downdraft u
      REAL(r8)     GCVD  (ISTS:IENS)         ! downdraft v
      REAL(r8)     FSNOW (ISTS:IENS)
      REAL(r8)     GMDDD (ISTS:IENS)

      REAL(r8)     GDTW,   GCHX,   GCTX,  GCQSX, GTPRP, EVSU,  GTEVE, LVIC,   &
                   DQW,    DTW,    GDQW,  DZ,    GCSD,  FDET,  GDHI,  GMDDX,  &
                   GMDDMX, GCHDX,  GCWDX, GCUDD, GCVDD, GTHCI, GTQVCI,        &
                   GTQLCI, GTQICI, wrk, wrk1, wrk2, wrk3, wrk4,               &
                   WMX, HMX, DDWMX, DDHMX, dp_above, dp_below, fsigma,        &
                   fmelt, fevp

!M    REAL(r8)     GTHCI, GTQVCI, GTQLCI, GTQICI, GTUCI, GTVCI
!DD#ifdef OPT_CUMBGT
! Water, energy, downdraft water and downdraft energy budgets
      REAL(r8), dimension(ISTS:IENS) :: WBGT, HBGT, DDWBGT, DDHBGT, tx1
      integer      ij, i, k, kp1
!DD#endif
!
!   [INTERNAL PARM]
      REAL(r8), parameter :: TWSNOW = 273.15_r8   ! wet-bulb temp. rain/snow
      REAL(r8), parameter :: FTMLT  = 4._r8       ! temp. factor for melt
      REAL(r8), parameter :: GMFLXC = 5.e-2_r8    ! critical mass flux
      REAL(r8), parameter :: VTERMS = 2._r8       ! terminal velocity of snowflake
      REAL(r8), parameter :: MELTAU = 10._r8      ! melting timescale
!
      REAL(r8), parameter :: EVAPR  = 0.3_r8      ! evaporation factor
!     REAL(r8), parameter :: EVAPR  = 0._r8       ! evaporation factor
      REAL(r8), parameter :: REVPDD = 1._r8       ! max rate of DD to evapolation
      REAL(r8), parameter :: RDDR   = 5.e-4_r8    ! DD rate (T0 R0 W0)^-1
!     REAL(r8), parameter :: RDDR   = 0._r8       ! DD rate (T0 R0 W0)^-1
      REAL(r8), parameter :: RDDMX  = 0.5_r8      ! norm. flux of downdraft
      REAL(r8), parameter :: VTERM  = 5._r8       ! term. vel. of precip.
      REAL(r8), parameter :: EVATAU = 2._r8       ! evaporation/sublimation timescale
      REAL(r8), parameter :: ZDMIN  = 5.e2_r8     ! min altitude of downdraft detrainment
      real(r8), parameter :: evapovtrm=EVAPR/VTERM

!NOTE
! downdraft area ffraction still needs to be computed for AW, assumed zero for now,
!   as passed in.

!
! Note: It is assumed that condensate evaporates in downdraft air.
!
      do k=1,kmax
        do i=ists,iens
          GPRCP (I,k) = zero
          GSNWP (I,k) = zero
          GMDD  (I,k) = zero
          GTEVP (I,k) = zero
          EVAPD (I,k) = zero
          SUBLD (I,k) = zero 
          EVAPE (I,k) = zero 
          SUBLE (I,k) = zero 
          EVAPX (I,k) = zero 
          SUBLX (I,k) = zero 
          GMDDE (I,k) = zero 
          SNMLT (I,k) = zero 
          GCHDD (I,k) = zero 
          GCWDD (I,k) = zero 
          GTTEV (I,k) = zero 
          GTQEV (I,k) = zero 
          GCdseD(I,k) = zero 
          GCqvD (I,k) = zero 
          GCqlD (I,k) = zero 
          GCqiD (I,k) = zero 
          gtevap(I,k) = zero 
          gtmelt(I,k) = zero 
          gtsubl(I,k) = zero 
        enddo
      enddo
!  testing on oct 17 2016
      if (do_aw) then
        if (.not. do_awdd) then
          do k=1,kmax
            do i=ists,iens
              dtdwn (i,k) = gtt(i,k)
              dqvdwn(i,k) = gtq(i,k,1)
              dqldwn(i,k) = gtq(i,k,itl)
              dqidwn(i,k) = gtq(i,k,iti)
            enddo
          enddo
        else
          do k=1,kmax
            do i=ists,iens
              dtdwn (I,k) = zero
              dqvdwn(I,k) = zero
              dqldwn(I,k) = zero
              dqidwn(I,k) = zero
            enddo
          enddo
        endif
      endif
!
      do i=ists,iens
        GCHD(I) = zero
        GCWD(I) = zero
        GCUD(I) = zero
        GCVD(I) = zero
      enddo
!
      DO K=KTMX,1,-1   ! loop A
        kp1 = min(k+1,kmax)
!
!     < precipitation melt & freeze >
!
        DO I=ISTS,IENS
          GTPRP = GPRCP(I,KP1) + GSNWP(I,KP1)
          IF (GTPRP > zero) THEN
             FSNOW(I) = GSNWP(I,KP1) / GTPRP
          ELSE
             FSNOW(I) = zero
          ENDIF
          LVIC  = ELocp + EMELTocp*FSNOW(I)
          GDTW  = GDT(I,K) - LVIC*(GDQS(I,K) - GDQ(I,K,1)) &
                           / (one + LVIC*FDQS(I,K))
          IF (GDTW  < TWSNOW) THEN
            GSNWP(I,K) = GSNWP(I,KP1) + GPRCI(I,K) + GSNWI(I,K)
            GTTEV(I,K) = EMELToCP*GPRCI(I,K) * GRAV/DELP(I,K)
            SNMLT(I,K) = -GPRCI(I,K)
          ELSE
            DZ   = GDZM(I,KP1) - GDZM(I,K)
            FMELT      = (one + FTMLT*(GDTW - TWSNOW))     &
                       * (one - TANH(GMFLX(I,KP1)/GMFLXC)) &
                       * (one - TANH(VTERMS*MELTAU/DZ))
            SNMLT(I,K) = GSNWP(I,KP1)*max(min(FMELT, one), zero)
            GSNWP(I,K) = GSNWP(I,KP1)+GSNWI(I,K) - SNMLT(I,K)
            GPRCP(I,K) = GPRCP(I,KP1)+GPRCI(I,K) + SNMLT(I,K)
            GTTEV(I,K) = -EMELToCP*SNMLT(I,K) * GRAV/DELP(I,K)
          ENDIF
!DD heating rate due to precip melting for AW
          gtmelt(i,k) = gtmelt(i,k) + GTTEV(I,K)
        ENDDO
!
!     < downdraft >
!
        DO I=ISTS,IENS   ! loop B
          wrk  = grav / delp(i,k)
          wrk1 = oneocp * wrk
          DZ   = GDZM(I,KP1) - GDZM(I,K)
          FEVP = (one - TANH(EVATAU*VTERM/DZ))
          IF (GMDD(I,KP1) > zero) THEN
            GCHX  = GCHD(I) / GMDD(I,KP1)
            GCTX  = GDT(I,K)  + (GCHX-GDHS(I,K)) / (CP+EL*FDQS(I,K))
            GCQSX = GDQS(I,K) + FDQS(I,K) * (GCTX - GDT(I,K))
            GCQSX = GCQSX*GMDD(I,KP1)
            EVSU  = MAX(GCQSX-GCWD(I), zero) * FEVP
            GTPRP = GPRCP(I,K) + GSNWP(I,K)
            IF (GTPRP > zero)  THEN
              FSNOW(I) = GSNWP(I,K) / GTPRP
            ELSE
              FSNOW(I) = zero
            ENDIF
            EVAPD(I,K)  = min(EVSU*(one-FSNOW(I)), GPRCP(I,K))
            SUBLD(I,K)  = min(EVSU*FSNOW(I), GSNWP(I,K))
            GPRCP(I,K)  = GPRCP(I,K) - EVAPD(I,K)
            GSNWP(I,K)  = GSNWP(I,K) - SUBLD(I,K)
! temperature tendencies due to evaporation and sublimation of precip
!  This is within downdraft
            gtevap(i,k) = gtevap(i,k) - elocp   * evapd(i,k) * wrk
            gtsubl(i,k) = gtsubl(i,k) - esubocp * subld(i,k) * wrk
            GCWD(I)     = GCWD(I) + EVAPD(I,K) + SUBLD(I,K)
            GCHD(I)     = GCHD(I) - EMELT*SUBLD(I,K)
          ENDIF

          GMDD(I,K) = GMDD(I,KP1)
!
          LVIC = ELocp + EMELTocp*FSNOW(I)
          DQW  = (GDQS(I,K) - GDW(I,K)) / (one + LVIC*FDQS(I,K))
          DQW  = MAX(DQW, zero)
          DTW  = LVIC*DQW
          GDQW = GDW(I,K) + DQW*FEVP
!
          EVSU = min(one, EVAPOVTRM*DQW*DZ*FEVP)
          EVAPE(I,K) = EVSU*GPRCP(I,K)
          SUBLE(I,K) = EVSU*GSNWP(I,K)
          GTEVP(I,K) = EVAPD(I,K) + SUBLD(I,K) + EVAPE(I,K) + SUBLE(I,K)
!
          GTPRP      = GPRCP(I,K) + GSNWP(I,K)
          GPRCP(I,K) = GPRCP(I,K) - EVAPE(I,K)
          GSNWP(I,K) = GSNWP(I,K) - SUBLE(I,K)
! additional temperature tendencies due to evaporation and sublimation of precip
!    This is outside of downdraft
          gtevap(i,k) = gtevap(i,k) - el*evape(i,k) * wrk1
          gtsubl(i,k) = gtsubl(i,k) - (el+emelt)*suble(i,k) * wrk1
!
          GMDDD(I) = zero
          IF (GDZ(I,K)-GDZM(I,1) > ZDMIN) THEN
            GTEVE      = EVAPE(I,K) + SUBLE(I,K)
            GMDDMX     = REVPDD*GTEVE/MAX(DQW, 1.D-10)
            GMDDE(I,K) = RDDR * (DTW*GTPRP*DELP(I,K))
            GMDDE(I,K) = MAX(MIN(GMDDE(I,K), GMDDMX), zero)
            GMDDX      = GMDD(I,KP1) + GMDDE(I,K)
            EVSU       = GMDDE(I,K)*DQW*FEVP
            IF (GTEVE > zero) THEN
              FSNOW(I) = SUBLE(I,K) / GTEVE
            ELSE
              FSNOW(I) = zero
            END IF
            EVAPX(I,K) = (one-FSNOW(I)) * EVSU
            SUBLX(I,K) = FSNOW(I) * EVSU
!
            IF (GMDDX > zero) THEN
              GDHI  = GDH(I,K) - EMELT*GDQI(I,K)
              GCHDX = GCHD(I) + GDHI*GMDDE(I,K) - EMELT*SUBLX(I,K)
              GCWDX = GCWD(I) + GDQW*GMDDE(I,K)
              GCSD  = (GCHDX - EL*GCWDX) / GMDDX
              IF (GCSD < GDS(I,K)) THEN
                GCHD(I)    = GCHDX
                GCWD(I)    = GCWDX
                GCUD(I)    = GCUD(I) + GDU(I,K)*GMDDE(I,K)
                GCVD(I)    = GCVD(I) + GDV(I,K)*GMDDE(I,K)
                GMDD(I,K)  = GMDDX
                EVAPE(I,K) = EVAPE(I,K) - EVAPX(I,K)
                SUBLE(I,K) = SUBLE(I,K) - SUBLX(I,K)
                EVAPD(I,K) = EVAPD(I,K) + EVAPX(I,K)
                SUBLD(I,K) = SUBLD(I,K) + SUBLX(I,K)
                GMDDD(I)   = zero
              ELSE
                GMDDE(I,K) = zero
                GMDDD(I)   = GMDD(I,KP1)
              ENDIF
            ENDIF
          ELSE
            GMDDD(I) = DZ / (GDZM(I,KP1)-GDZM(I,1)) * GMDD(I,KP1)
          ENDIF
!
          GMDDD(I) = MAX(GMDDD(I), GMDD(I,K)-RDDMX*GMFLX(I,K))
!
          IF (GMDDD(I) > zero) THEN
            FDET       = GMDDD(I)/GMDD(I,K)
            GCHDD(I,K) = FDET*GCHD(I)
            GCWDD(I,K) = FDET*GCWD(I)
            GCUDD      = FDET*GCUD(I)
            GCVDD      = FDET*GCVD(I)
!
            GTHCI  =  wrk * (GCHDD(I,K) - GMDDD(I)*GDH(I,K))
            GTQVCI =  wrk * (GCWDD(I,K) - GMDDD(I)*GDQ(I,K,1))
            GTQLCI = -wrk * GMDDD(I)*GDQ(I,K,ITL)
            GTQICI = -wrk * GMDDD(I)*GDQI(I,K)
!
            GTT (I,K)     = GTT(I,K)     + (GTHCI - EL*GTQVCI)*oneoCP
            GTQ (I,K,1)   = GTQ(I,K,1)   + GTQVCI
            GTQ (I,K,ITL) = GTQ(I,K,ITL) + GTQLCI
            GTQ (I,K,ITI) = GTQ(I,K,ITI) + GTQICI
!           GTQI(I,K)     = GTQI(I,K)    + GTQICI

            GTU (I,K) = GTU(I,K) + wrk * (GCUDD - GMDDD(I)*GDU(I,K))
            GTV (I,K) = GTV(I,K) + wrk * (GCVDD - GMDDD(I)*GDV(I,K))
!
            GCHD(I)   = GCHD(I)   - GCHDD(I,K)
            GCWD(I)   = GCWD(I)   - GCWDD(I,K)
            GCUD(I)   = GCUD(I)   - GCUDD
            GCVD(I)   = GCVD(I)   - GCVDD
            GMDD(I,K) = GMDD(I,K) - GMDDD(I)
          ENDIF
          GCdseD(I,K) = GCHD(I) - el*GCWD(I)
          GCqvD (I,K) = GCWD(I)
        ENDDO   ! loop B
!
      ENDDO   ! loop A
!
      do i=ists,iens
        tx1(i) = GRAV / DELP(I,1)
      enddo
      DO K=1,KTMX
        kp1 = min(k+1,kmax)
        DO I=ISTS,IENS
          wrk    = tx1(i)
          tx1(i) = GRAV / DELP(I,kp1)
            
          GTTEV(I,K) = GTTEV(I,K) - wrk                              &
                     * (ELocp*EVAPE(I,K)+(ELocp+EMELTocp)*SUBLE(I,K))
          GTT(I,K)   = GTT(I,K) + GTTEV(I,K)
!
          GTQEV(I,K) = GTQEV(I,K) + (EVAPE(I,K)+SUBLE(I,K)) * wrk
          GTQ(I,K,1) = GTQ(I,K,1) + GTQEV(I,K)
!
          GMFLX(I,K) = GMFLX(I,K) - GMDD(I,K)

! AW tendencies due to vertical divergence of eddy fluxes
          if (do_awdd .and. k > 1) then
            fsigma        = one - sigmad(i,kp1)
            dp_below      = wrk    * (one - sigmad(i,k))
            dp_above      = tx1(i) * (one - sigmad(i,kp1))

            wrk1          = gmdd(i,kp1) * (gdt(i,k)+gocp*gdz(i,k)) - gcdsed(i,kp1)*oneocp
            wrk2          = gmdd(i,kp1) * gdq(i,k,1) - gcqvd(i,kp1)
            wrk3          = gmdd(i,kp1) * gdq(i,k,itl)
            wrk4          = gmdd(i,kp1) * gdqi(i,k)

            dtdwn(i,k)    = dtdwn(i,k)    + dp_below * wrk1
            dqvdwn(i,k)   = dqvdwn(i,k)   + dp_below * wrk2
            dqldwn(i,k)   = dqldwn(i,k)   + dp_below * wrk3 ! gcqld=0   - gcqld(i,k))
            dqidwn(i,k)   = dqidwn(i,k)   + dp_below * wrk4 ! gcqid=0   - gcqid(i,k))

            dtdwn(i,kp1)  = dtdwn(i,kp1)  - dp_above * wrk1
            dqvdwn(i,kp1) = dqvdwn(i,kp1) - dp_above * wrk2
            dqldwn(i,kp1) = dqldwn(i,kp1) - dp_above * wrk3 ! gcqld=0   - gcqld(i,k))
            dqidwn(i,kp1) = dqidwn(i,kp1) - dp_above * wrk4 ! gcqid=0   - gcqid(i,k))
          endif

        ENDDO   ! end of i loop
      ENDDO     ! end of k loop
!
      if (.not. do_awdd) then
        do k=1,kmax
          do i=ists,iens
            dtdwn(i,k)  = gtt(i,k)     - dtdwn(i,k)
            dqvdwn(i,k) = gtq(i,k,1)   - dqvdwn(i,k)
            dqldwn(i,k) = gtq(i,k,itl) - dqldwn(i,k)
            dqidwn(i,k) = gtq(i,k,iti) - dqidwn(i,k)
!!          dqidwn(i,k) = gtqi(i,k)    - dqidwn(i,k)
          enddo
        enddo
      endif
!
      END SUBROUTINE CUMDWN
!***********************************************************************
      SUBROUTINE CUMCLD                             & !! cloudiness
               ( IJSDIM, KMAX  ,                    & !DD dimensions
                 CUMCLW, QLIQ  , QICE  , FLIQC  ,   & ! modified
                 CUMFRC,                            & ! output
                 GMFLX , KTMX  , ISTS, IENS   )       ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, KMAX  ! DD, for GFS, pass in
!
!   [OUTPUT]
      REAL(r8)     CUMFRC(IJSDIM)          ! cumulus cloud fraction
!
!   [MODIFY]
      REAL(r8)     CUMCLW(IJSDIM, KMAX)    ! cloud water in cumulus
      REAL(r8)     QLIQ  (IJSDIM, KMAX)    ! cloud liquid
      REAL(r8)     QICE  (IJSDIM, KMAX)    ! cloud ice
      REAL(r8)     FLIQC (IJSDIM, KMAX)    ! liquid ratio in cumulus
!
!   [INPUT]
      REAL(r8)     GMFLX (IJSDIM, KMAX)   ! cumulus mass flux
      INTEGER      KTMX
      INTEGER      ISTS, IENS
!
!   [WORK]
      INTEGER      I, K
      REAL(r8)     CUMF, QC, wrk
      LOGICAL, SAVE :: OFIRST = .TRUE.
!
!   [INTERNAL PARAM]
      REAL(r8) :: FACLW  = 0.1_r8     ! Mc->CLW
      REAL(r8) :: CMFMIN = 2.e-3_r8   ! Mc->cloudiness
      REAL(r8) :: CMFMAX = 3.e-1_r8   ! Mc->cloudiness
      REAL(r8) :: CLMIN  = 1.e-3_r8   ! cloudiness Min.
      REAL(r8) :: CLMAX  = 0.1_r8     ! cloudiness Max.
      REAL(r8), SAVE :: FACLF
!
      IF (OFIRST) THEN
         FACLF = (CLMAX-CLMIN) / LOG(CMFMAX/CMFMIN)
         OFIRST = .FALSE.
      END IF
!
      CUMFRC(ISTS:IENS) = zero
      DO K=1,KTMX
        DO I=ISTS,IENS
          CUMFRC(I) = MAX(CUMFRC(I), GMFLX(I,K))
        ENDDO
      ENDDO
      DO I=ISTS,IENS
        IF (CUMFRC(I) > zero) THEN
          CUMF      = LOG(MAX(CUMFRC(I), CMFMIN)/CMFMIN)
          CUMFRC(I) = MIN(FACLF*CUMF+CLMIN, CLMAX)
        ENDIF
      ENDDO
!
      DO K=1,KTMX
        DO I=ISTS,IENS
          IF (GMFLX(I,K) > zero) THEN
            wrk = FACLW / GMFLX(I,K) * CUMFRC(I)
            QLIQ  (I,K) = wrk * QLIQ(I,K)
            QICE  (I,K) = wrk * QICE(I,K)
            CUMCLW(I,K) = wrk * CUMCLW(I,K)
            QC          = QLIQ(I,K) + QICE(I,K)
            IF (QC > zero) THEN
              FLIQC(I,K) = QLIQ(I,K) / QC
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
      END SUBROUTINE CUMCLD
!***********************************************************************
      SUBROUTINE CUMUPR                                    & !! Tracer Updraft
               ( im    , IJSDIM, KMAX  , NTR   ,           & !DD dimensions
                 GTR   , GPRCC ,                           & ! modified
                 GDR   , CBMFX , ELAM  , GDZ   , GDZM  ,   & ! input
                 GCYM  , GCYT  , GCQT  , GCLT  , GCIT  ,   & ! input
                 GTPRT , GTEVP , GTPRC0,                   & ! input
                 KB    , KBMX  , KT    , KTMX  , KTMXT ,   & ! input
                 DELP  , OTSPT , ISTS  , IENS,             & ! input
                 fscav, fswtr, nctp)
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: im, IJSDIM, KMAX, NTR, nctp             !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(r8)     GTR   (IJSDIM, KMAX, NTR)
      REAL(r8)     GPRCC (IJSDIM, NTR)
!
!   [INPUT]
      REAL(r8)     GDR   (IJSDIM, KMAX, NTR)
      REAL(r8)     CBMFX (IM, NCTP)
      REAL(r8)     ELAM  (IJSDIM, KMAX, NCTP)
      REAL(r8)     GDZ   (IJSDIM, KMAX)
      REAL(r8)     GDZM  (IJSDIM, KMAX+1)
      REAL(r8)     GCYM  (IJSDIM, KMAX)
      REAL(r8)     GCYT  (IJSDIM, NCTP)
      REAL(r8)     GCQT  (IJSDIM, NCTP)
      REAL(r8)     GCLT  (IJSDIM, NCTP)
      REAL(r8)     GCIT  (IJSDIM, NCTP)
      REAL(r8)     GTPRT (IJSDIM, NCTP)
      REAL(r8)     GTEVP (IJSDIM, KMAX)
      REAL(r8)     GTPRC0(IJSDIM)   !! precip. before evap.
      real(r8)     fscav(ntr), fswtr(ntr)
      INTEGER      KB    (IJSDIM )
      INTEGER      KBMX
      INTEGER      KT    (IJSDIM, NCTP)
      INTEGER      KTMX  (NCTP)
      INTEGER      KTMXT
      REAL(r8)     DELP  (IJSDIM, KMAX)
      LOGICAL      OTSPT (NTR)              !! transport with this routine?
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER      I, K, LT, TP, CTP
      REAL(r8)                            :: GCRTD, SCAV, GCWT, GPRCR, evpf, cbmfxl
      REAL(r8), dimension(ists:iens)      :: GCRB, GCRT,  DR,   gtprc0i
!     REAL(r8), dimension(ists:iens,kmax) :: DGCB, DZ,    RDZM,  EVPF
      REAL(r8), dimension(ists:iens,kmax) ::       DZ,    RDZM
!     REAL(r8), dimension(ists:iens,nctp) :: DZT,  RGCWT, MASK1, MASK2
      REAL(r8), dimension(ists:iens,nctp) :: DZT,  RGCWT, MASK1
!
!     DO K=1,KBMX
!       DO I=ISTS,IENS
!         DGCB(I,K) = GCYM(I,K+1) - GCYM(I,K)
!       ENDDO
!     ENDDO
      do i=ists,iens
        if (gtprc0(i) > zero) then
          gtprc0i(i) = one / gtprc0(i)
        else
          gtprc0i(i) = zero
        endif
      enddo
      DO K=1,KTMXT
        DO I=ISTS,IENS
          DZ  (I,K) = GDZM(I,K+1) - GDZM(I,K)
          RDZM(I,K) = GRAV / DELP(I,K)
!         EVPF(I,K) = zero
!         IF (GTPRC0(I) > zero) THEN
!           EVPF(I,K) = GTEVP(I,K) / GTPRC0(I)
!         ENDIF
        ENDDO
      ENDDO
      DO CTP=1,NCTP
        DO I=ISTS,IENS
          K = KT(I,CTP)
!
          GCWT = GCQT(I,CTP) + GCLT(I,CTP) + GCIT(I,CTP)
          RGCWT(I,CTP) = zero
          IF (GCWT > zero) THEN
            RGCWT(I,CTP) = one / GCWT
          ENDIF
!
          MASK1(I,CTP) = zero
          DZT  (I,CTP) = zero
          IF (K > KB(I)) THEN
            MASK1(I,CTP) = one
            DZT  (I,CTP) = GDZ(I,K) - GDZM(I,K)
          ENDIF
!         MASK2(I,CTP) = zero
!         IF (CBMFX(I,CTP) > zero) then
!           MASK2(I,CTP) = one
!         ENDIF
        ENDDO
      ENDDO
!
      DO LT=1,NTR   ! outermost tracer LT loop
!
        IF (OTSPT(LT)) THEN
          GCRB = zero
          DO K=1,KBMX
            DO I=ISTS,IENS
              IF (K < KB(I)) THEN
!               GCRB(I) = GCRB(I) + DGCB(I,K) * GDR(I,K,LT)
                GCRB(I) = GCRB(I) + (GCYM(I,K+1)-GCYM(I,K))* GDR(I,K,LT)
              ENDIF
            ENDDO
          ENDDO
!
          DO CTP=1,NCTP
            DR = zero
            DO K=2,KTMX(CTP)
              DO I=ISTS,IENS
                IF (K >= KB(I) .AND.  K < KT(I,CTP)) THEN
                  DR(I) = DR(I) + DZ(I,K) * ELAM(I,K,CTP) * GDR(I,K,LT)
                ENDIF
              ENDDO
            ENDDO
!
            DO I=ISTS,IENS
              K = MAX(KT(I,CTP),1)
              DR(I)   = DR(I) + DZT(I,CTP) * ELAM(I,K,CTP) * GDR (I,K,LT) &
                                           * MASK1(I,CTP)
              GCRT(I) = (GCRB(I) + DR(I))  * MASK1(I,CTP)
!
              SCAV    = FSCAV(LT)*GTPRT(I,CTP) + FSWTR(LT)*GTPRT(I,CTP)*RGCWT(I,CTP)
              SCAV    = MIN(SCAV, one)
              GCRTD   = GCRT(I) * (one - SCAV)
              cbmfxl  = max(zero, CBMFX(I,CTP))
              GPRCR   = SCAV * GCRT(I) * CBMFXl

              GTR(I,K,LT) = GTR(I,K,LT) + RDZM(I,K) * CBMFXl          &
                              * (GCRTD - GCYT(I,CTP) * GDR(I,K,LT))
              GPRCC(I,LT) = GPRCC(I,LT) + GPRCR

!             GPRCR   = SCAV * GCRT(I) * CBMFX(I,CTP)
!             GTR(I,K,LT) = GTR(I,K,LT) + RDZM(I,K) * CBMFX(I,CTP) &
!                             * (GCRTD - GCYT(I,CTP) * GDR(I,K,LT)) * MASK2(I,CTP)
!             GPRCC(I,LT) = GPRCC(I,LT) + GPRCR * MASK2(I,CTP)
            ENDDO
          ENDDO
!
          DO K=KTMXT,1,-1
            DO I=ISTS,IENS
              evpf = GTEVP(i,k) * gtprc0i(i)
              GTR(I,K,LT) = GTR(I,K,LT) + RDZM(I,K) * GPRCC(I,LT) * EVPF
              GPRCC(I,LT) = GPRCC(I,LT) * (one - EVPF)
!             GTR(I,K,LT) = GTR(I,K,LT) + RDZM(I,K) * GPRCC(I,LT) * EVPF(I,K)
!             GPRCC(I,LT) = GPRCC(I,LT) * (one - EVPF(I,K))
            ENDDO
          ENDDO
!
        ENDIF
!
      ENDDO   ! outermost tracer LT loop
!
      END SUBROUTINE CUMUPR
!***********************************************************************
      SUBROUTINE CUMDNR                                 & !! Tracer Downdraft
                      ( IM    , IJSDIM, KMAX  , NTR   , & !DD dimensions
                        GTR   ,                         & ! modified
                        GDR   , GMDD  , DELP  ,         & ! input
                        KTMX  , OTSPT , ISTS  , IENS )    ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX, NTR             !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(r8)     GTR   (IJSDIM, KMAX, NTR)   ! Temperature tendency
!
!   [INPUT]
      REAL(r8)     GDR   (IJSDIM, KMAX, NTR)
      REAL(r8)     GMDD  (IJSDIM, KMAX)        ! downdraft mass flux
      REAL(r8)     DELP  (IJSDIM, KMAX  )
      LOGICAL      OTSPT (NTR)
      INTEGER      KTMX, ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     GCRD  (ISTS:IENS)           ! downdraft q
      REAL(r8)     GMDDE, GMDDD, GCRDD
      INTEGER      I, K, LT, kp1
!
!
      DO LT=1,NTR
        IF (OTSPT(LT)) THEN
          GCRD = zero
          DO K=KTMX,1,-1
            kp1 = min(k+1,kmax)
            DO I=ISTS,IENS
              GMDDE = GMDD(I,K) - GMDD(I,KP1)
              IF (GMDDE >= zero) THEN
                GCRD(I) = GCRD(I) + GDR(I,K,LT)*GMDDE
              ELSEIF (GMDD(I,KP1) > zero) THEN
                GMDDD = - GMDDE
                GCRDD = GMDDD/GMDD(I,KP1) * GCRD(I)
                GTR(I,K,LT) = GTR(I,K,LT) + GRAV/DELP(I,K) &
                                          * (GCRDD - GMDDD*GDR(I,K,LT))
                GCRD(I)     = GCRD(I) - GCRDD
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
      END SUBROUTINE CUMDNR
!***********************************************************************
      SUBROUTINE CUMSBR                                      & !! Tracer Subsidence
                      ( IM    , IJSDIM, KMAX  , NTR   ,      & !DD dimensions
                        GTR   ,                              & ! modified
                        GDR   , DELP  ,                      & ! input
                        GMFLX , KTMX  , OTSPT ,              & ! input
                        ISTS, IENS )                           ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX, NTR             !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(r8)     GTR   (IJSDIM, KMAX, NTR)   !! tracer tendency
!
!   [INPUT]
      REAL(r8)     GDR   (IJSDIM, KMAX, NTR)   !! tracer
      REAL(r8)     DELP  (IJSDIM, KMAX)
      REAL(r8)     GMFLX (IJSDIM, KMAX)        !! mass flux
      INTEGER      KTMX
      LOGICAL      OTSPT (NTR)                 !! tracer transport on/off
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER      I, K, KM, KP, LT
      REAL(r8)     SBR0, SBR1, FX1
      REAL(r8)     FX(ISTS:IENS)
!
      DO LT=1,NTR
        IF (OTSPT(LT)) THEN
          DO I=ISTS,IENS
            FX(I) = zero
          enddo
          DO K=KTMX,1,-1
            KM = MAX(K-1, 1)
            KP = MIN(K+1, KMAX)
            DO I=ISTS,IENS
              SBR0 = GMFLX(I,KP) * (GDR(I,KP,LT) - GDR(I,K,LT))
              SBR1 = GMFLX(I,K)  * (GDR(I,K,LT)  - GDR(I,KM,LT))
              IF (GMFLX(I,K) > GMFLX(I,KP)) THEN
                FX1 = half
              ELSE
                FX1 = zero
              END IF
              GTR(I,K,LT) = GTR(I,K,LT) + GRAV/DELP(I,K)              &
                                        * ((one-FX(I))*SBR0 + FX1*SBR1)
              FX(I) = FX1
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
      END SUBROUTINE CUMSBR
!*********************************************************************
      SUBROUTINE CUMFXR                                           & ! Tracer mass fixer
                      ( IM    , IJSDIM, KMAX  , NTR   ,           & !DD dimensions
                        GTR   ,                                   & ! modified
                        GDR   , DELP  , DELTA , KTMX  , IMFXR ,   & ! input
                        ISTS  , IENS                            )   ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX, NTR             !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(r8)     GTR   (IJSDIM, KMAX, NTR)   ! tracer tendency
!
!   [INPUT]
      REAL(r8)     GDR   (IJSDIM, KMAX, NTR)   ! tracer
      REAL(r8)     DELP  (IJSDIM, KMAX)
      REAL(r8)     DELTA                       ! time step
      INTEGER      KTMX
      INTEGER      IMFXR (NTR)
        ! 0: mass fixer is not applied
        !    tracers which may become negative values
        !    e.g. subgrid-PDFs
        ! 1: mass fixer is applied, total mass may change through cumulus scheme
        !    e.g. moisture, liquid cloud, ice cloud, aerosols
        ! 2: mass fixer is applied, total mass never change through cumulus scheme
        !    e.g. CO2
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     GDR1
      REAL(r8)     GDR2  (ISTS:IENS, KMAX)
      REAL(r8), dimension(ISTS:IENS) :: TOT0, TOT1, TRAT
      REAL(r8)     FWAT
      INTEGER      I, K, LT
!
! Attention: tracers are forced to be positive unless IMFXR=0.
!
      DO LT=1,NTR
        SELECT CASE (IMFXR(LT))
          CASE (0)
            CYCLE
          CASE (1)
            FWAT = one
          CASE (2)
            FWAT = zero
          CASE DEFAULT
            EXIT
        END SELECT
!
        DO I=ISTS,IENS
          TOT0(I) = zero
          TOT1(I) = zero
        enddo
!
        DO K=KTMX,1,-1
          DO I=ISTS,IENS
            IF (GTR(I,K,LT) /= zero) THEN
              GDR1      = GDR(I,K,LT) + DELTA*GTR(I,K,LT)
              GDR2(I,K) = MAX(GDR1, zero)
              GDR1      = GDR1 * FWAT + GDR(I,K,LT)*(one - FWAT)
              TOT0(I)   = TOT0(I) + GDR1 *(DELP(I,K)*GRAVI)
              TOT1(I)   = TOT1(I) + GDR2(I,K)*(DELP(I,K)*GRAVI)
            ENDIF
          ENDDO
        ENDDO
!
        DO I=ISTS,IENS
          IF (TOT1(I) > zero ) THEN
            TRAT(I) = MAX(TOT0(I), zero) / TOT1(I)
          ELSE
            TRAT(I) = one
          ENDIF
        ENDDO
!
        DO K=KTMX,1,-1
          DO I=ISTS,IENS
            IF (GTR(I,K,LT) /= zero ) THEN
              GDR2(I,K   ) = GDR2(I,K)*TRAT(I)
              GTR (I,K,LT) = (GDR2(I,K)-GDR(I,K,LT)) / DELTA
            ENDIF
          ENDDO
        ENDDO
!
      ENDDO   ! LT-loop
!
      END SUBROUTINE CUMFXR
!*********************************************************************
      SUBROUTINE CUMFXR1                                   & ! Tracer mass fixer
               ( IM    , IJSDIM, KMAX  ,                   & !DD dimensions
                 GTR   ,                                   & ! modified
                 GDR   , DELP  , DELTA , KTMX  , IMFXR ,   & ! input
                 ISTS  , IENS                            )   ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX  ! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(r8)     GTR   (IJSDIM, KMAX)      ! tracer tendency
!
!   [INPUT]
      REAL(r8)     GDR   (IJSDIM, KMAX)      ! tracer
      REAL(r8)     DELP  (IJSDIM, KMAX)
      REAL(r8)     DELTA                     ! time step
      INTEGER      KTMX
      INTEGER      IMFXR
        ! 0: mass fixer is not applied
        !    tracers which may become negative values
        !    e.g. subgrid-PDFs
        ! 1: mass fixer is applied, total mass may change through cumulus scheme
        !    e.g. moisture, liquid cloud, ice cloud, aerosols
        ! 2: mass fixer is applied, total mass never change through cumulus scheme
        !    e.g. CO2
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     GDR1
      REAL(r8)     GDR2  (ISTS:IENS, KMAX)
      REAL(r8), dimension(ISTS:IENS) :: TOT0, TOT1, TRAT
      REAL(r8)     FWAT
      INTEGER      I, K
!
! Attention: tracers are forced to be positive unless IMFXR=0.
!
      SELECT CASE (IMFXR)
        CASE (0)
          RETURN
        CASE (1)
          FWAT = one
        CASE (2)
          FWAT = zero
        CASE DEFAULT
          RETURN
      END SELECT
!
      DO I=ISTS,IENS
        TOT0(I) = zero
        TOT1(I) = zero
      enddo
!
      DO K=KTMX,1,-1
        DO I=ISTS,IENS
          IF (GTR(I,K) /= zero) THEN
            GDR1      = GDR(I,K) + DELTA*GTR(I,K)
            GDR2(I,K) = MAX(GDR1, zero)
            GDR1      = GDR1*FWAT + GDR(I,K)*(one - FWAT)
            TOT0(I)   = TOT0(I) + GDR1 *(DELP(I,K)*GRAVI)
            TOT1(I)   = TOT1(I) + GDR2(I,K)*(DELP(I,K)*GRAVI)
          ENDIF
        ENDDO
      ENDDO
!
      DO I=ISTS,IENS
        IF (TOT1(I) > zero) THEN
          TRAT(I) = MAX(TOT0(I), zero) / TOT1(I)
        ELSE
          TRAT(I) = one
        ENDIF
      ENDDO
!
      DO K=KTMX,1,-1
        DO I=ISTS,IENS
          IF (GTR(I,K) /= zero) THEN
            GDR2(I,K) = GDR2(I,K)*TRAT(I)
            GTR (I,K) = (GDR2(I,K)-GDR(I,K)) / DELTA
          ENDIF
        ENDDO
      ENDDO
!
      END SUBROUTINE CUMFXR1
!*********************************************************************
      SUBROUTINE CUMCHK                                   & ! check range of output values
                      ( IJSDIM, KMAX  , NTR   ,           & !DD dimensions
                        GTT   , GTQ   , GTU   , GTV   ,   & ! input
                        GTCFRC, GPRCC , GSNWC , CUMCLW,   & ! input
                        CUMFRC, FLIQC , GTPRP ,           & ! input
                        ISTS  , IENS                    )   ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, KMAX, NTR ! DD, for GFS, pass in
!
!   [INPUT]
      REAL(r8)     GTT   (IJSDIM, KMAX)      ! heating rate
      REAL(r8)     GTQ   (IJSDIM, KMAX, NTR) ! change in q
      REAL(r8)     GTU   (IJSDIM, KMAX)      ! tendency of u
      REAL(r8)     GTV   (IJSDIM, KMAX)      ! tendency of v
      REAL(r8)     GPRCC (IJSDIM, NTR )      ! rainfall
      REAL(r8)     GSNWC (IJSDIM)            ! snowfall
      REAL(r8)     CUMCLW(IJSDIM, KMAX)      ! cloud water in cumulus
      REAL(r8)     CUMFRC(IJSDIM)            ! cumulus cloud fraction
      REAL(r8)     GTCFRC(IJSDIM, KMAX)      ! change in cloud fraction
      REAL(r8)     FLIQC (IJSDIM, KMAX)      ! liquid ratio in cumulus
      REAL(r8)     GTPRP (IJSDIM, KMAX)      ! rain+snow flux
!
      INTEGER    ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER    I, K
!
!   [INTERNAL PARM]
      REAL(r8) :: GTTMAX  = 1.e-2_r8
      REAL(r8) :: GTQVMAX = 1.e-4_r8
      REAL(r8) :: GTQLMAX = 1.e-5_r8
      REAL(r8) :: GTUMAX  = 1.e-2_r8
      REAL(r8) :: GTVMAX  = 1.e-2_r8
      REAL(r8) :: GTCFMAX = 1.e-3_r8
      REAL(r8) :: PRCCMAX = 1.e-2_r8
      REAL(r8) :: SNWCMAX = 1.e-2_r8
      REAL(r8) :: CLWMAX  = 1.e-3_r8
      REAL(r8) :: TPRPMAX = 1.e-2_r8
      REAL(r8) :: GTQIMAX = 1.e-5_r8
      REAL(r8) :: GTM2MAX = 1._r8
      REAL(r8) :: GTM3MAX = 1._r8
!
      DO K=1,KMAX
        DO I=ISTS, IENS
          IF (ABS(GTT(I,K)) > GTTMAX) THEN
            WRITE(iulog,*) '### CUMCHK: GTT(',I,',',K,')=',GTT(I,K)
          ENDIF
          IF (ABS(GTQ(I,K,1) ) > GTQVMAX) THEN
            WRITE(iulog,*) '### CUMCHK: GTQ(',I,',',K,',1 )=', GTQ(I,K,1)
          ENDIF
          IF (ABS(GTQ(I,K,ITL)) > GTQLMAX) THEN
            WRITE(iulog,*) '### CUMCHK: GTQ(',I,',',K,',ITL )=', GTQ(I,K,ITL)
          ENDIF
          IF (ABS(GTU(I,K)) > GTUMAX) THEN
            WRITE(iulog,*) '### CUMCHK: GTU(',I,',',K,')=',GTU(I,K)
          END IF
          IF (ABS(GTV(I,K)) > GTVMAX) THEN
            WRITE(iulog,*) '### CUMCHK: GTV(',I,',',K,')=',GTV(I,K)
          ENDIF
          IF (ABS(GTCFRC(I,K)) > GTCFMAX) THEN
            WRITE(iulog,*) '### CUMCHK: GTCFRC(',I,',',K,')=', GTCFRC(I,K)
          ENDIF
          IF (CUMCLW(I,K) > CLWMAX .OR. CUMCLW(I,K) < zero) THEN
            WRITE(iulog,*) '### CUMCHK: CUMCLW(',I,',',K,')=', CUMCLW(I,K)
          ENDIF
          IF (FLIQC(I,K) > one .OR.  FLIQC(I,K) < zero) THEN
            WRITE(iulog,*) '### CUMCHK: FLIQC(',I,',',K,')=', FLIQC(I,K)
          ENDIF
          IF (GTPRP(I,K) > TPRPMAX .OR.  GTPRP(I,K) < zero) THEN
            WRITE(iulog,*) '### CUMCHK: GTPRP(',I,',',K,')=', GTPRP(I,K)
          ENDIF
          IF (ABS(GTQ(I,K,ITI)) > GTQIMAX) THEN
            WRITE(iulog,*) '### CUMCHK: GTQ(',I,',',K,',ITI )=', GTQ(I,K,ITI)
          ENDIF
!         IF (ABS(GTQ(I,K,IMU2) ) > GTM2MAX) THEN
!           WRITE(iulog,*) '### CUMCHK: GTQ(',I,',',K,',IMU2 )=', GTQ(I,K,IMU2)
!         ENDIF
!         IF (ABS(GTQ(I,K,IMU3)) > GTM3MAX) THEN
!           WRITE(iulog,*) '### CUMCHK: GTQ(',I,',',K,',IMU3 )=', GTQ(I,K,IMU3)
!         ENDIF
        ENDDO
      ENDDO
!
      DO I=ISTS,IENS
        IF (GPRCC(I,1) > PRCCMAX .OR. GPRCC(I,1) < zero) THEN
          WRITE(iulog,*) '### CUMCHK: GPRCC(',I,')=',GPRCC(I,1)
        END IF
        IF (GSNWC(I) > SNWCMAX .OR. GSNWC(I) < zero) THEN
          WRITE(iulog,*) '### CUMCHK: GSNWC(',I,')=',GSNWC(I)
        END IF
        IF (CUMFRC(I) > one .OR.  CUMFRC(I) < zero) THEN
          WRITE(iulog,*) '### CUMCHK: CUMFRC(',I,')=',CUMFRC(I)
        ENDIF
      ENDDO
!
      END SUBROUTINE CUMCHK
!***********************************************************************
      SUBROUTINE TINTP                          & ! vertical interpolation of temperature
                     ( IJSDIM, KMAX  ,          & !DD dimensions
                       GDTM  ,                  & ! output
                       GDT   , GDP   , GDPM ,   & ! input
                       ISTS  , IENS           )   ! input

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJSDIM, KMAX       ! DD, for GFS, pass in
!*
!*   [OUTPUT]
      REAL(r8)     GDTM  (IJSDIM, KMAX+1)     ! temperature (half lev)
!*
!*   [INPUT]
      REAL(r8)     GDT   (IJSDIM, KMAX)       ! temperature (full lev)
      REAL(r8)     GDP   (IJSDIM, KMAX)       ! pressure (full lev)
      REAL(r8)     GDPM  (IJSDIM, KMAX+1)     ! pressure (half lev)
      INTEGER      ISTS, IENS                   ! range of active grids
!*
!*   [INTERNAL WORK]
!     REAL(r8)     FTINT ( KMAX )               ! intrp. coef.
!     REAL(r8)     FTINTM( KMAX )               ! intrp. coef.
      real (r8)  :: wrk, wrk1, ftintm

      INTEGER    I, K
!*
!*          < interp. temp. >
!*
      DO K=2,KMAX
        DO I=ISTS,IENS
          wrk  = one / GDP(I,K)
          wrk1 = one / LOG(GDP(I,K-1)*wrk)
          FTINTM      = wrk1 * LOG(GDPM(I,K)*wrk)
          GDTM(I,K) = FTINTM *GDT(I,K-1) + (1.0-FTINTM)*GDT(I,K)
!         FTINTM( K ) = wrk1 * LOG(GDPM(I,K)*wrk)
!         FTINT ( K ) = wrk1 * LOG(GDP(I,K-1)/GDPM(I,K))
!         GDTM( I,K ) = FTINTM(K)*GDT(I,K-1) + FTINT(K)*GDT(I,K)
         ENDDO
      ENDDO

      DO I = ISTS, IENS
        GDTM(I,KMAX+1) = GDT(I,KMAX)
        GDTM(I,1     ) = GDT(I,1)
      ENDDO

      RETURN
      END SUBROUTINE TINTP
!***********************************************************************

end module cs_conv

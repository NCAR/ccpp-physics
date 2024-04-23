!>  \file cs_conv.F90
!!  This file contains the Chikira-Sugiyama Convection scheme.

module cs_conv
!---------------------------------------------------------------------------------
! Purpose:
!
!>---------------------------------------------------------------------------------
! Purpose:
!
!> Interface for Chikira-Sugiyama convection scheme 
!!
!! Author: Minoru Chikira
!---------------------------------------------------------------------------------
!
  use machine ,   only : kind_phys
  use physcons,   only : cp    => con_cp,   grav   => con_g,                   &
     &                   rair  => con_rd,   rvap   => con_rv,                  &
     &                   cliq  => con_cliq, cvap   => con_cvap,                &
     &                   epsv  => con_eps,  epsvm1 => con_epsm1,               &
     &                   epsvt => con_fvirt,                                   &
     &                   el    => con_hvap, emelt  => con_hfus, t0c => con_t0c
  use funcphys,   only : fpvs ! this is saturation vapor pressure in funcphys.f

  
  implicit none

  private                ! Make default type private to the module

   real(kind_phys), parameter :: zero=0.0d0,  one=1.0d0, half=0.5d0
   real(kind_phys), parameter :: cpoel=cp/el, cpoesub=cp/(el+emelt), esubocp=1.0/cpoesub, &
                          elocp=el/cp, oneocp=one/cp, gocp=grav/cp, gravi=one/grav,&
                          emeltocp=emelt/cp, cpoemelt=cp/emelt, epsln=1.e-10_kind_phys

   real(kind_phys), parameter :: fact1=(cvap-cliq)/rvap, fact2=el/rvap-fact1*t0c !< to calculate d(qs)/dT

   logical,  parameter :: adjustp=.true.
!  logical,  parameter :: adjustp=.false.

! Tuning parameters set from namelist
!
!  real(kind_phys), parameter, public :: CLMD = 0.60,   & !< entrainment efficiency (now thru argument)
   real(kind_phys), parameter, public ::                &
                                  PA=0.15,       & !< factor for buoyancy to affect updraft velocity
                                  CPRES = 0.55,  & !< pressure factor for momentum transport
                                  ALP0 = 5.0e7,  & !< alpha parameter in prognostic closure
!                                 ALP0 = 8.0e7,  & !< alpha parameter in prognostic closure
!                                 CLMP = (one-CLMD)*(PA+PA), &
!                                 CLMDPA = CLMD*PA,          &
                                  spblmin=0.05, &  !< minimum cloudbase height in p/ps
                                  spblmax=0.30, &  !< maximum cloudbase height in p/ps
!                                 spblcrit=0.03, & !< minimum cloudbase height in p/ps
!                                 spblcrit=0.035,& !< minimum cloudbase height in p/ps
!                                 spblcrit=0.025,& !< minimum cloudbase height in p/ps
                                  cincrit= -150.0
!                                 cincrit= -120.0
!                                 cincrit= -100.0

!DD precz0 and  preczh control partitioning of water between detrainment
!DD   and precipitation. Decrease for more precip

   real(kind_phys), public       ::  precz0, preczh, clmd, clmp, clmdpa
   real(kind_phys), public, parameter :: c0t=0.002, d0t=0.002
!
! Private data
!
  real(kind_phys), parameter     :: unset_kind_phys = -999._kind_phys   ! missing value
!
  integer :: iulog !< unit to write debugging and diagnostic output
                   !DD Note - see if I can find corresponding variable in a GFS module
!
! Shared variables
!
  integer, parameter :: ITI = 2, ITL = 3  !< index of ice and liquid water

  integer, save, dimension(50) :: IMFXR   !< 0: mass fixer is not applied
                                          !!    tracers which may become negative
                                          !!    values e.g. subgrid-PDFs
                                          !! 1: mass fixer is applied, total mass
                                          !!    may change through cumulus scheme
                                          !!    e.g. moisture, liquid cloud, ice
                                          !!    cloud, aerosols
                                          !! 2: mass fixer is applied, total mass
                                          !!    never change through cumulus scheme
                                          !!    e.g. CO2

!  PUBLIC: interfaces
!
   public  cs_conv_run         ! CS scheme main driver
  
   contains

!>\defgroup cs_scheme Chikira-Sugiyama Cumulus Scheme Module
!> \brief The subroutine contains the main driver for Chikira-Sugiyama convective scheme.
!!
!! \author Minoru Chikira
!!
!! History:
!! - Jun 26 2014 : D. Dazlich - Modified for GFS
!! - Apr 10 2015 : S. Moorthi - check for allocatable arrays and fix argument for cbmfx
!! - Oct    2015 : D. Dazlich - Add computation of updraft area fraction (sigma) for
!!                             diagnostic purposes.
!! - Aug    2016 : D. Dazlich - Create flux form of tendencies and multiply by
!!                             Arakawa-Wu functions of sigma
!! - Sep    2016 : S. Moorthi - found two bugs - cleanup and some optimization
!! - Oct    2016 : S. Moorthi - added sigma affects on tracers and CUMFLX and CUMDET
!!                             made many cosmetic changes
!! - Nov    2016 : S. Moorthi - further optimization and cleanup and several bug fixes
!! - April  2017 : S. moorthi - many changes including removing elam and making gcym
!!                             a function of cloud type.  This makes it possible for
!!                             AW affect propagate to other routines such as CUMUPR
!! - Apr 12, 2017 : S. Moorthi - Added flx_form logical and relevant code to compute AW
!!                             without flux form when false.
!! - May 17, 2017 : S. Moorthi - Added routine CUMSBW for just momentum change
!!                             in advective form
!! - Sep 08, 2017 : D. Dazlich - tracers in flux form for AW
!! - Nov     2017 : S. Moorthi - fix some bugs and fix fluxform for tracers
!! - Nov 22  2017 : S. Moorthi - add kcnv array to identify points where deep convection
!!                             operates - 0 - no convection 1 - with convection
!! - Jan 30 2018  : S. Moorthi - fixed sigmad dimension error in CUMDWN and an error when adjustp=.true.
!! - May    2018  : S. Moorthi - modified cumup to compute total workfunction (positive plus negative)
!!                             and negative part only and to let a particular ensemble exist only if
!!                             the ratio of negative to total is less than some prescribed percent.
!!                             Also, added an extra iteration in this k loop. Reduced some memory.
!! - June   2018  : S. Moorthi - the output mass fluxes ud_mf, dd_mf and dt_mf are over time step delta
!!
!! \b Arakawa-Wu \b implemtation: 
!! for background, consult An Introduction to the
!! General Circulation of the Atmosphere, Randall, chapter six.
!! Traditional parameterizations compute tendencies like those in eq 103, 105 and 106.
!! Because Arakawa-Wu applies different functions to different components to the
!! terms within these equations, it requires the terms used in alternate eqns 91 - 93.
!! The code required to compute these terms is added within, and the appropriate
!! functions of updraft area fraction (sigma) are applied. Thus, AW requires three
!! steps:
!! -# computation of the updraft area fraction
!! -# alternative representation of the tendency terms
!! -# application of functions of sigma to the alternative tendency terms
!!  here, and in gbphys to the large-scale microphysics tendencies.
!!
!!  The bulk of AW is implemented within subroutine CS_CUMLUS(), and the routines it calls.
!!
!!
!! JLS NOTE:  The convective mass fluxes (dt_mf, dd_mf and ud_mf) passed in and out of cs_conv have not been multiplied by
!!            the timestep (kg/m2/sec) as they are in all other convective schemes.  EMC is aware of this problem, 
!!            and in the future will be fixing this discrepancy.  In the meantime, CCPP will use the same mass flux standard_name
!!            and long_name as the other convective schemes, where the units are in kg/m2. (Aug 2018)
!!
!! \section arg_table_cs_conv_run Argument Table
!! \htmlinclude cs_conv_run.html
!!
!!  \section general_cs_conv CS Convection Scheme General Algorithm
!> @{
   subroutine cs_conv_run(         IJSDIM ,  KMAX     , ntracp1 , NN,       &
                          NTR    , nctp   ,                                 & !DD dimensions
                          otspt  , lat    ,  kdt      ,                     &
                          t      , q      ,  rain1    , clw     ,           &
                          zm     , zi     ,  pap      , paph    ,           &
                          delta  , delti  ,  ud_mf    , dd_mf   , dt_mf,    &
                          u      , v      ,  fscav    , fswtr,              &
                          cbmfx  , mype   ,  wcbmaxm  , precz0in, preczhin, &
                          clmdin , sigma  , do_aw     , do_awdd , flx_form, &
                          lprnt  , ipr, kcnv,                               &
                          QLCN, QICN, w_upi, cf_upi, CNV_MFD,               & ! for coupling to MG microphysics
                          CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,       &
                          mp_phys,errmsg,errflg)


   implicit none
!
! input arguments
!
   INTEGER, INTENT(IN)     :: IJSDIM, KMAX, ntracp1, nn, NTR, mype, nctp, mp_phys, kdt, lat !! DD, for GFS, pass in
   logical, intent(in)     :: otspt(:,:)          ! otspt(:,1) - on/off switch for tracer transport by updraft and
                                                  !              downdraft. should not include subgrid PDF and turbulence
                                                  ! otspt(:,2) - on/off switch for tracer transport by subsidence
                                                  !              should include subgrid PDF and turbulence

   real(kind_phys), intent(inout) :: t(:,:)          ! temperature at mid-layer (K)
   real(kind_phys), intent(inout) :: q(:,:)          ! water vapor array including moisture (kg/kg)
   real(kind_phys), intent(inout) :: clw(:,:,:)      ! tracer array including cloud condensate (kg/kg)
   real(kind_phys), intent(in)    :: pap(:,:)        ! pressure at mid-layer (Pa)
   real(kind_phys), intent(in)    :: paph(:,:)       ! pressure at boundaries (Pa)
   real(kind_phys), intent(in)    :: zm(:,:)         ! geopotential at mid-layer (m)
   real(kind_phys), intent(in)    :: zi(:,:)         ! geopotential at boundaries (m)
   real(kind_phys), intent(in)    :: fscav(:), fswtr(:), wcbmaxm(:)
   real(kind_phys), intent(in)    :: precz0in, preczhin, clmdin
! added for cs_convr
   real(kind_phys), intent(inout) :: u(:,:)          ! zonal wind at mid-layer (m/s)
   real(kind_phys), intent(inout) :: v(:,:)          ! meridional wind at mid-layer (m/s)

   real(kind_phys), intent(in)    :: DELTA           ! physics time step
   real(kind_phys), intent(in)    :: DELTI           ! dynamics time step (model time increment in seconds)
   logical,  intent(in)    :: do_aw, do_awdd, flx_form
!
! modified arguments
!
   real(kind_phys), intent(inout), optional :: CBMFX(:,:)      ! cloud base mass flux (kg/m2/s)
!
! output arguments
!
!  updraft, downdraft, and detrainment mass flux (kg/m2/s)
   real(kind_phys), intent(inout), dimension(:,:), optional :: ud_mf
   real(kind_phys), intent(inout), dimension(:,:) :: dd_mf, dt_mf
   
   real(kind_phys), intent(out)   :: rain1(:)        ! lwe thickness of deep convective precipitation amount (m)
! GJF* These variables are conditionally allocated depending on whether the
!     Morrison-Gettelman microphysics is used, so they must be declared 
!     using assumed shape.
   real(kind_phys), intent(out), dimension(:,:), optional :: qlcn, qicn, w_upi,cnv_mfd, &
                                                   cnv_dqldt, clcn, cnv_fice, &
                                                   cnv_ndrop, cnv_nice, cf_upi
! *GJF
   logical, intent(in)    :: lprnt
   integer, intent(in)    :: ipr
   integer, intent(inout) :: kcnv(:)          ! zero if no deep convection and 1 otherwise
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!DDsigma - output added for AW sigma diagnostics
!  interface sigma and vertical velocity by cloud type (1=sfc) 
!  real(kind_phys), intent(out), dimension(:,:,:)  :: sigmai, vverti
   real(kind_phys), intent(out), dimension(:,:)    :: sigma  ! sigma  sigma totaled over cloud type - on interfaces (1=sfc)
!   sigma  terms in eq 91 and 92
!  real(kind_phys), dimension(IJSDIM,KMAX)                    :: sfluxterm, qvfluxterm, condterm
!DDsigma
!
! output arguments of CS_CUMLUS
!
   real(kind_phys), dimension(IJSDIM,KMAX+1,nctp)  :: vverti, sigmai

   real(kind_phys) GTT(IJSDIM,KMAX)           !< temperature tendency [K/s]
   real(kind_phys) GTQ(IJSDIM,KMAX,NTR)       !< tracer tendency [kg/kg/s]
   real(kind_phys) GTU(IJSDIM,KMAX)           !< zonal velocity tendency [m/s2]
   real(kind_phys) GTV(IJSDIM,KMAX)           !< meridional velocity tendency [m/s2]
   real(kind_phys) CMDET(IJSDIM,KMAX)         !< detrainment mass flux [kg/m2/s]
   real(kind_phys) GTPRP(IJSDIM,KMAX+1)       !< precipitation (including snowfall) flux at interfaces [kg/m2/s]
   real(kind_phys) GSNWP(IJSDIM,KMAX+1)       !< snowfall flux at interfaces [kg/m2/s]
   real(kind_phys) GMFX0(IJSDIM,KMAX+1)       !< updraft mass flux [kg/m2/s]
   real(kind_phys) GMFX1(IJSDIM,KMAX+1)       !< downdraft mass flux [kg/m2/s]
   integer  KT(IJSDIM,nctp)            !< cloud top index for each cloud type

   real(kind_phys) :: cape(IJSDIM)            !< convective available potential energy (J/kg)
   real(kind_phys) :: prec(IJSDIM)            !< precipitation at surface (including snowfall) (kg/m2/s)
   real(kind_phys) :: snow(IJSDIM)            !< snowfall at surface (kg/m2/s)
!
! input arguments of CS_CUMLUS
!
   real(kind_phys) GDT(IJSDIM,KMAX)           !< temperature [K]
   real(kind_phys) GDQ(IJSDIM,KMAX,NTR)       !< tracers including moisture [kg/kg]  !DDsigmadiag
   real(kind_phys) GDU(IJSDIM,KMAX)           !< zonal wind [m/s]
   real(kind_phys) GDV(IJSDIM,KMAX)           !< meridional wind [m/s]
   real(kind_phys) GDTM(IJSDIM,KMAX+1)        !< temperature at boundaries of layers [K]
   real(kind_phys) GDP(IJSDIM,KMAX)           !< pressure [Pa]
   real(kind_phys) GDPM(IJSDIM,KMAX+1)        !< pressure at boundaries of layers [Pa]
   real(kind_phys) GDZ(IJSDIM,KMAX)           !< altitude [m]
   real(kind_phys) GDZM(IJSDIM,KMAX+1)        !< altitude at boundaries of layers [m]
   real(kind_phys) delp(IJSDIM,KMAX)          !< pressure difference between layers [Pa]
   real(kind_phys) delpi(IJSDIM,KMAX)         !< grav/delp
!
! local variables
!
!DD   real(kind_phys) :: zs(IJSDIM)           !< surface height [m]

   integer KTMAX(IJSDIM)               !< max of KT
   real(kind_phys)    :: ftintm, wrk, wrk1, tem
   integer i, k, n, ISTS, IENS, kp1

!DD borrowed from RAS to go form total condensate to ice/water separately
!  parameter (tf=130.16, tcr=160.16, tcrf=1.0/(tcr-tf),tcl=2.0)
!  parameter (tf=230.16, tcr=260.16, tcrf=1.0/(tcr-tf))
   real(kind_phys), parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf), tcl=2.0
   logical, save       :: first=.true.

   ! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

!  lprnt = kdt == 1 .and. mype == 38
!  ipr = 43

   precz0 = precz0in
   preczh = preczhin
   clmd   = clmdin
   CLMP   = (one-CLMD)*(PA+PA)
   CLMDPA = CLMD*PA
!
   if (first) then
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
       delpi(i,k) = grav / delp(i,k)
     enddo
   enddo

!DD following adapted from ras
!> -# Following the Relaxed Arakawa Schubert Scheme (RAS;  
!! Moorthi and Suarez 1992 \cite moorthi_and_suarez_1992 ), 
!! separate total condensate between ice and water.
!! The ratio of cloud ice to cloud water is determined by a linear function
!! of temperature:
!!\f[
!! F_i(T)= (T_2-T)/(T_2-T_1)
!!\f]
!! where T is temperature, and\f$T_1\f$ and \f$T_2\f$ are set as tcf=263.16
!! and tf= 233.16 
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
!  if (lprnt) write(0,*)'in cs clw1b=',clw(ipr,:,1),' kdt=',kdt
!  if (lprnt) write(0,*)'in cs clw2b=',clw(ipr,:,2),' kdt=',kdt

   do n=2,NTR
     do k=1,KMAX
       do i=1,IJSDIM
         GDQ(i,k,n) = clw(i,k,n-1)
       enddo
     enddo
   enddo
!  if (lprnt) write(0,*)' incs tke=',gdq(ipr,1:25,ntr)
!
!***************************************************************************************
!
!> -# Calculate temperature at interfaces
!

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

!> -# Initialize the sigma diagnostics
   do n=1,nctp
     do k=1,kmax+1
       do i=ists,iens
         vverti(i,k,n) = zero
         sigmai(i,k,n) = zero
       enddo
     enddo
   enddo
   do k=1,kmax+1
     do i=ists,iens
       sigma(i,k)  = zero
     enddo
   enddo
!
!> -# Call cs_cumlus() for the main CS cumulus parameterization
   call CS_CUMLUS (IJSDIM, IJSDIM, KMAX  , NTR   ,    &  !DD dimensions
                   otspt(1:ntr,1), otspt(1:ntr,2),    &
                   lprnt , ipr   ,                    &
                   GTT   , GTQ   , GTU   , GTV   ,    & ! output
                   CMDET ,                            & ! output
                   GTPRP , GSNWP , GMFX0 ,            & ! output
                   GMFX1 , cape  , KT    ,            & ! output
                   CBMFX ,                            & ! modified
                   GDT   , GDQ   , GDU   , GDV   ,    & ! input
                   GDTM  ,                            & ! input
                   GDP   , GDPM  , GDZ   , GDZM  ,    & ! input
                   delp  , delpi ,                    &
                   DELTA , DELTI , ISTS  , IENS, mype,& ! input
                   fscav,  fswtr,  wcbmaxm, nctp,     &
                   sigmai, sigma,  vverti,            & ! input/output !DDsigma
                   do_aw, do_awdd, flx_form)
!
!
!DD detrainment has to be added in for GFS
!
!  if (lprnt) write(0,*)' aft cs_cum gtqi=',gtq(ipr,:,2)
!  if (lprnt) write(0,*)' aft cs_cum gtql=',gtq(ipr,:,3)


   do n=2,NTR
     do k=1,KMAX
       do i=1,IJSDIM
         clw(i,k,n-1) = max(zero, GDQ(i,k,n) + GTQ(i,k,n) * delta)
       enddo
     enddo
   enddo
!  if (lprnt) write(0,*)' aftcs_cum tkein=',gdq(ipr,1:25,ntr),' delta=',delta
!  if (lprnt) write(0,*)' aftcs_cum tke=',clw(ipr,1:25,ntr-1)
!  if (lprnt) write(0,*)'in cs clw1a=',clw(ipr,:,1),' kdt=',kdt
!  if (lprnt) write(0,*)'in cs clw2a=',clw(ipr,:,2),' kdt=',kdt
!
   do k=1,KMAX
     do i=1,IJSDIM
       q(i,k)        = max(zero, GDQ(i,k,1)  + GTQ(i,k,1) * delta)
       t(i,k)        = GDT(i,k)    + GTT(i,k)   * delta
       u(i,k)        = GDU(i,k)    + GTU(i,k)   * delta
       v(i,k)        = GDV(i,k)    + GTV(i,k)   * delta
! Set the mass fluxes.
        ud_mf  (i,k) = GMFX0(i,k)
        dd_mf  (i,k) = GMFX1(i,k)
        dt_mf  (i,k) = CMDET(i,k)
     enddo
   enddo


   if (mp_phys == 10) then  ! for 2M microphysics, always output these variables
     if (do_aw) then
       do k=1,KMAX
         kp1 = min(k+1,kmax)
         do i=1,IJSDIM
           qicn(i,k)      = max(0.0, clw(i,k,1)-gdq(i,k,2))
           qlcn(i,k)      = max(0.0, clw(i,k,2)-gdq(i,k,3))

           
           wrk = qicn(i,k) + qlcn(i,k)
           if (wrk > 1.0e-12) then
             cnv_fice(i,k)  = qicn(i,k) / wrk
           else
             cnv_fice(i,k)  = 0.0
           endif
!
!          CNV_MFD(i,k)   = dt_mf(i,k) * (1.0/delta)
           CNV_MFD(i,k)   = dt_mf(i,k)
           CNV_DQLDT(i,k) = wrk / delta
!          CNV_PRC3(i,k)  = 0.0
           CNV_NDROP(i,k) = 0.0
           CNV_NICE(i,k)  = 0.0
           cf_upi(i,k)    = max(0.0,min(0.01*log(1.0+500*ud_mf(i,k)),0.1))
!          CLCN(i,k)      = cf_upi(i,k)                     !downdraft is below updraft
!!         clcn(i,k)      = max(0.0,min(0.01*log(1.0+500*ud_mf(i,k)/delta),0.25))

           w_upi(i,k)     = 0.0
!!         w_upi(i,k)     = ud_mf(i,k)*(t(i,k)+epsvt*gdq(i,k,1)) * rair &
!!                       / (delta*max(cf_upi(i,k),1.e-12)*gdp(i,k))
         enddo
       enddo
       do k=1,KMAX
         do i=1,ijsdim
           do n=1,nctp
             w_upi(i,k) = w_upi(i,k) + vverti(i,k,n)
           enddo
           if (sigma(i,k) > 1.0e-10) then
             w_upi(i,k) = w_upi(i,k) / sigma(i,k)
           else
             w_upi(i,k) = 0.0
           endif
         enddo
       enddo
     else
       do k=1,KMAX
         do i=1,IJSDIM
           qicn(i,k)      = max(0.0, clw(i,k,1)-gdq(i,k,2))
           qlcn(i,k)      = max(0.0, clw(i,k,2)-gdq(i,k,3))
           cnv_fice(i,k)  = qicn(i,k) / max(1.0e-10,qicn(i,k)+qlcn(i,k))
! 
!          CNV_MFD(i,k)   = dt_mf(i,k) * (1/delta)
           CNV_MFD(i,k)   = dt_mf(i,k)
           CNV_DQLDT(i,k) = (qicn(i,k)+qlcn(i,k)) / delta
!          CNV_PRC3(i,k)  = 0.0
           CNV_NDROP(i,k) = 0.0
           CNV_NICE(i,k)  = 0.0
           cf_upi(i,k)    = max(0.0,min(0.01*log(1.0+500*ud_mf(i,k)),0.1))
!    &                                               500*ud_mf(i,k)),0.60))
!          CLCN(i,k)      = cf_upi(i,k)                     !downdraft is below updraft
           
           w_upi(i,k)     = ud_mf(i,k)*(t(i,k)+epsvt*gdq(i,k,1)) * rair &
                          / (max(cf_upi(i,k),1.e-12)*gdp(i,k))
         enddo
       enddo
     endif
   endif

!****************************************************************************
 
   KTMAX = 1
   do n=1,nctp
     do i=1,IJSDIM
        KTMAX(i) = max(KTMAX(i), KT(i,n))
     enddo
   enddo
!
   do i=1,IJSDIM
     prec(i)  = GTPRP(i,1)
     snow(i)  = GSNWP(i,1)
     if (prec(i)+snow(i) > 0.0) then
       kcnv(i) = 1
     else
       kcnv(i) = 0
     endif
   enddo

!> -# Multiply mass fluxes by the time step

   do k=1,kmax
     do i=1,ijsdim
       ud_mf(i,k) = ud_mf(i,k) * delta
       dd_mf(i,k) = dd_mf(i,k) * delta
       dt_mf(i,k) = dt_mf(i,k) * delta
     enddo
   enddo

!   rain1(:) = prec(:) * (delta*0.001)      ! Convert prec flux kg/m2/sec to rain1 in m
    do i= 1, IJSDIM
     rain1(i) = prec(i) * (delta*0.001)
    enddo

!  if (lprnt) then
!    write(0,*)' aft cs_cum prec=',prec(ipr),'GTPRP=',GTPRP(ipr,1)
!  endif
  

!    if (do_aw) then
!    call moist_bud(ijsdim,ijsdim,im,kmax,mype,kdt,grav,delta,delp,prec &
!    ,              gdq(1,1,1), gdq(1,1,2), gdq(1,1,3)                  &
!    ,              q,clw(1,1,1),clw(1,1,2),'cs_conv_aw')
!    endif

   end subroutine cs_conv_run
!> @}


!************************************************************************
!>\ingroup cs_scheme
!! Main subroutine for the cumulus parameterization with
!! state-dependent entrainment rate developed by Minoru Chikira.
!!
!! - This routine works as the prognostic Arakawa-Schubert scheme
!!   if OPT_ASMODE is specified.
!! - Specify OPT_NS02 to use entrainment rate of Neggers et al. (2002)
!! - Specify OPT_CUMBGT to check water and energy budget.
!! - Specify OPT_CUMCHK to check range of output values.
!!
!! History(yy/mm/dd):
!! - 08/09/19(chikira)   MIROC4.1
!! - 08/10/30(hiro)      CMT modified
!! - 08/11/11(chikira)   Neggers et al. (2002)
!! - 08/12/3 (chikira)   downdraft detrainment modified
!! - 08/12/3 (chikira)   COSP output
!! - 09/02/24(chikira)   fix convective inhibition
!! - 09/04/16(hiro)      CMIP5 output (cbasep,ctopp)
!! - 09/09/03(yokohata)  COSP
!! - 10/11/19(toshi)     small bug fix
!! - 14/02/07(chikira)   CUMDWN bug fix, CMT modified
!!\section gen_cs_cumlus CSAW cs_cumlus General Algorithm
!> @{
   SUBROUTINE CS_CUMLUS (im    , IJSDIM, KMAX  , NTR   ,    & !DD dimensions
                         otspt1, otspt2, lprnt , ipr   ,    &
                         GTT   , GTQ   , GTU   , GTV   ,    & ! output
                         CMDET ,                            & ! output
                         GTPRP , GSNWP , GMFX0 ,            & ! output
                         GMFX1 , CAPE  , KT    ,            & ! output
                         CBMFX ,                            & ! modified
                         GDT   , GDQ   , GDU   , GDV   ,    & ! input
                         GDTM  ,                            & ! input
                         GDP   , GDPM  , GDZ   , GDZM  ,    & ! input
                         delp  , delpinv ,                  &
                         DELTA , DELTI , ISTS  , IENS, mype,& ! input
                         fscav,  fswtr,  wcbmaxm, nctp,     & !
                         sigmai, sigma,  vverti,            & ! input/output !DDsigma
                         do_aw, do_awdd, flx_form)
!
   IMPLICIT NONE
      
   Integer, parameter    :: ntrq=4                    ! starting index for tracers
   INTEGER, INTENT(IN)   :: im, IJSDIM, KMAX, NTR, mype, nctp, ipr !! DD, for GFS, pass in
   logical, intent(in)   :: do_aw, do_awdd, flx_form  ! switch to apply Arakawa-Wu to the tendencies
   logical, intent(in)   :: otspt1(ntr), otspt2(ntr), lprnt
   REAL(kind_phys),intent(in)   :: DELP  (IJSDIM, KMAX)
   REAL(kind_phys),intent(in)   :: DELPINV (IJSDIM, KMAX)
!
! [OUTPUT]
   REAL(kind_phys), INTENT(OUT) :: GTT   (IJSDIM, KMAX     ) ! heating rate
   REAL(kind_phys), INTENT(OUT) :: GTQ   (IJSDIM, KMAX, NTR) ! change in q
   REAL(kind_phys), INTENT(OUT) :: GTU   (IJSDIM, KMAX     ) ! tendency of u
   REAL(kind_phys), INTENT(OUT) :: GTV   (IJSDIM, KMAX     ) ! tendency of v
   REAL(kind_phys), INTENT(OUT) :: CMDET (IJSDIM, KMAX     ) ! detrainment mass flux
   REAL(kind_phys) :: GTLDET( IJSDIM, KMAX      ) ! cloud liquid tendency by detrainment
   REAL(kind_phys) :: GTIDET( IJSDIM, KMAX      ) ! cloud ice tendency by detrainment
! assuming there is no flux  at the top of the atmospherea - Moorthi
   REAL(kind_phys), INTENT(OUT) :: GTPRP (IJSDIM, KMAX+1   ) ! rain+snow flux
   REAL(kind_phys), INTENT(OUT) :: GSNWP (IJSDIM, KMAX+1   ) ! snowfall flux
   REAL(kind_phys), INTENT(OUT) :: GMFX0 (IJSDIM, KMAX+1   ) ! updraft mass flux
   REAL(kind_phys), INTENT(OUT) :: GMFX1 (IJSDIM, KMAX+1   ) ! downdraft mass flux

   REAL(kind_phys), INTENT(OUT) :: CAPE  (IJSDIM           )
   INTEGER , INTENT(OUT) :: KT    (IJSDIM, NCTP     ) ! cloud top
!
!  [MODIFIED]
   REAL(kind_phys), INTENT(INOUT) :: CBMFX ( IM, NCTP        ) !! cloud base mass flux

   !DDsigma - output added for AW sigma diagnostics
   real(kind_phys), intent(out)   :: sigmai(IM,KMAX+1,nctp)  !DDsigma  sigma by cloud type - on interfaces (1=sfc)
   real(kind_phys), intent(out)   :: vverti(IM,KMAX+1,nctp)  !DDsigma  vert. vel. by cloud type - on interfaces (1=sfc)
   real(kind_phys), intent(out)   :: sigma(IM,KMAX+1)        !DDsigma  sigma totaled over cloud type - on interfaces (1=sfc)
   
! for computing AW flux form of tendencies
!  real(kind_phys), dimension(IM,KMAX) ::    &  !DDsigmadiag
!      sfluxterm, qvfluxterm
!  real(kind_phys), dimension(IM,KMAX) ::    &  !DDsigmadiag
!      qlfluxterm, qifluxterm
!  real(kind_phys),  dimension(ijsdim,kmax,ntrq:ntr) :: trfluxterm ! tendencies of tracers due to eddy mass flux
   real(kind_phys), dimension(IM,KMAX) ::    &  !DDsigmadiag
       condtermt, condtermq, frzterm, prectermq, prectermfrz
   !DDsigma

!
!  [INPUT]
   REAL(kind_phys), INTENT(IN) :: GDT   (IJSDIM, KMAX     ) ! temperature T
   REAL(kind_phys), INTENT(IN) :: GDQ   (IJSDIM, KMAX, NTR) ! humidity, tracer  !DDsigmadiag
   REAL(kind_phys), INTENT(IN) :: GDU   (IJSDIM, KMAX     ) ! westerly u
   REAL(kind_phys), INTENT(IN) :: GDV   (IJSDIM, KMAX     ) ! southern wind v
   REAL(kind_phys), INTENT(IN) :: GDTM  (IJSDIM, KMAX+1   ) ! temperature T
   REAL(kind_phys), INTENT(IN) :: GDP   (IJSDIM, KMAX     ) ! pressure P
   REAL(kind_phys), INTENT(IN) :: GDPM  (IJSDIM, KMAX+1   ) ! pressure (half lev)
   REAL(kind_phys), INTENT(IN) :: GDZ   (IJSDIM, KMAX     ) ! altitude
   REAL(kind_phys), INTENT(IN) :: GDZM  (IJSDIM, KMAX+1   ) ! altitude
   REAL(kind_phys), INTENT(IN) :: DELTA                     ! delta(t) (dynamics)
   REAL(kind_phys), INTENT(IN) :: DELTI                     ! delta(t) (internal variable)
   INTEGER,  INTENT(IN) :: ISTS, IENS                ! array range

   real(kind_phys), intent(in) :: fscav(ntr), fswtr(ntr), wcbmaxm(ijsdim)
!
!  [INTERNAL WORK]
   REAL(kind_phys), allocatable :: GPRCC (:, :)  ! rainfall
   REAL(kind_phys)     GSNWC ( IJSDIM            ) ! snowfall
   REAL(kind_phys)     CUMCLW( IJSDIM, KMAX      ) ! cloud water in cumulus
   REAL(kind_phys)     CUMFRC( IJSDIM            ) ! cumulus cloud fraction
!COSP
   REAL(kind_phys)     QLIQC ( IJSDIM, KMAX   )    ! cumulus cloud liquid water [kg/kg]
   REAL(kind_phys)     QICEC ( IJSDIM, KMAX   )    ! cumulus cloud ice [kg/kg]
   REAL(kind_phys)     GPRCPF( IJSDIM, KMAX   )    ! rainfall flux at full level
   REAL(kind_phys)     GSNWPF( IJSDIM, KMAX   )    ! snowfall flux at full level
!
   REAL(kind_phys)     GTCFRC( IJSDIM, KMAX      ) ! change in cloud fraction
   REAL(kind_phys)     FLIQC ( IJSDIM, KMAX      ) ! liquid ratio in cumulus
!
!#ifdef OPT_CHASER
!      REAL(kind_phys)     RFXC  ( IJSDIM, KMAX+1    ) ! precipi. flx [kg/m2/s]
!      REAL(kind_phys)     SFXC  ( IJSDIM, KMAX+1    ) ! ice/snow flx [kg/m2/s]
!      INTEGER      LEVCUM( IJSDIM, KMAX      ) ! flag for cum. cloud top
!      REAL(kind_phys)     LNFRC ( IJSDIM, KMAX      ) ! areal rates of clouds
!      REAL(kind_phys)     REVC  ( IJSDIM, KMAX      ) ! evaporation rates
!#endif
!
   REAL(kind_phys)     GDCFRC( IJSDIM, KMAX      ) ! cloud fraction
!
!   REAL(kind_phys)     GTQL  ( IJSDIM, KMAX )      ! tendency of cloud liquid
!
   REAL(kind_phys)     GDW   ( IJSDIM, KMAX )      ! total water
   REAL(kind_phys)     GDQS  ( IJSDIM, KMAX )      ! saturate moisture
   REAL(kind_phys)     FDQS  ( IJSDIM, KMAX )
   REAL(kind_phys)     GAM   ( IJSDIM, KMAX )
   REAL(kind_phys)     GDS   ( IJSDIM, KMAX )      ! dry static energy
   REAL(kind_phys)     GDH   ( IJSDIM, KMAX )      ! moist static energy
   REAL(kind_phys)     GDHS  ( IJSDIM, KMAX )      ! saturate MSE
!
   REAL(kind_phys)     GCYM  ( IJSDIM, KMAX, NCTP )      ! norm. mass flux (half lev)
   REAL(kind_phys)     GCHB  ( IJSDIM )            ! cloud base MSE-Li*Qi
   REAL(kind_phys)     GCWB  ( IJSDIM )            ! cloud base total water
   REAL(kind_phys)     GCtrB  ( IJSDIM, ntrq:ntr )            ! cloud base water vapor tracer
   REAL(kind_phys)     GCUB  ( IJSDIM )            ! cloud base U
   REAL(kind_phys)     GCVB  ( IJSDIM )            ! cloud base V
   REAL(kind_phys)     GCIB  ( IJSDIM )            ! cloud base ice
   REAL(kind_phys)     ELAM  ( IJSDIM, KMAX, NCTP )   ! entrainment (rate*massflux)
   REAL(kind_phys)     GCYT  ( IJSDIM, NCTP )      ! norm. mass flux @top
   REAL(kind_phys)     GCHT  ( IJSDIM, NCTP )      ! cloud top MSE
   REAL(kind_phys)     GCQT  ( IJSDIM, NCTP )      ! cloud top q
   REAL(kind_phys)     GCwT  ( IJSDIM )      ! cloud top total water
   REAL(kind_phys)     GCUT  ( IJSDIM, NCTP )      ! cloud top U
   REAL(kind_phys)     GCVT  ( IJSDIM, NCTP )      ! cloud top V
   REAL(kind_phys)     GCLT  ( IJSDIM, NCTP )      ! cloud top cloud water
   REAL(kind_phys)     GCIT  ( IJSDIM, NCTP )      ! cloud top cloud ice
   REAL(kind_phys)     GCtrT (IJSDIM, ntrq:ntr, NCTP) ! cloud top tracer
   REAL(kind_phys)     GTPRT ( IJSDIM, NCTP )      ! precipitation/M
   REAL(kind_phys)     GCLZ  ( IJSDIM, KMAX )      ! cloud liquid for each CTP
   REAL(kind_phys)     GCIZ  ( IJSDIM, KMAX )      ! cloud ice for each CTP

   REAL(kind_phys)     ACWF  ( IJSDIM       )      ! cloud work function
   REAL(kind_phys)     GPRCIZ( IJSDIM, KMAX+1, NCTP )    ! precipitation
   REAL(kind_phys)     GSNWIZ( IJSDIM, KMAX+1, NCTP )    ! snowfall
   REAL(kind_phys)     GTPRC0( IJSDIM       )      ! precip. before evap.

   REAL(kind_phys)     GMFLX ( IJSDIM, KMAX+1 )    ! mass flux (updraft+downdraft)
   REAL(kind_phys)     QLIQ  ( IJSDIM, KMAX   )    ! total cloud liquid
   REAL(kind_phys)     QICE  ( IJSDIM, KMAX   )    ! total cloud ice
   REAL(kind_phys)     GPRCI ( IJSDIM, KMAX   )    ! rainfall generation
   REAL(kind_phys)     GSNWI ( IJSDIM, KMAX   )    ! snowfall generation

   REAL(kind_phys)     GPRCP ( IJSDIM, KMAX+1 )    ! rainfall flux
!
   REAL(kind_phys)     GTEVP ( IJSDIM, KMAX   )    ! evaporation+sublimation
   REAL(kind_phys)     GMDD  ( IJSDIM, KMAX+1 )    ! downdraft mass flux

   REAL(kind_phys)     CUMHGT( IJSDIM, NCTP   )    ! cloud top height
   REAL(kind_phys)     CTOPP ( IJSDIM         )    ! cloud top pressure

   REAL(kind_phys)     GDZTR ( IJSDIM         )   ! tropopause height
   REAL(kind_phys)     FLIQOU( IJSDIM, KMAX   )   ! liquid ratio in cumulus
!#ifdef OPT_CHASER
!      REAL(kind_phys)     TOPFLX( IJSDIM, NCTP   )    !! flux at each cloud top
!#endif
   INTEGER    KB    ( IJSDIM )
   INTEGER    KSTRT ( IJSDIM ) ! tropopause level
   REAL(kind_phys)   GAMX
   REAL(kind_phys)   CIN   ( IJSDIM )
   INTEGER    JBUOY ( IJSDIM )
   REAL(kind_phys)   DELZ, BUOY, DELWC, DELER
!M REAL(kind_phys)   WCB   ( NCTP )                ! updraft velocity**2 @base
!M SAVE       WCB
   REAL(kind_phys)   WCBX (IJSDIM)
!   REAL(kind_phys)   ERMR  ( NCTP )                ! entrainment rate (ASMODE)
!   SAVE       ERMR
   INTEGER    KTMX  ( NCTP )                ! max of cloud top
   INTEGER    KTMXT                         ! max of cloud top
   REAL(kind_phys)   TIMED
   REAL(kind_phys)   GDCLDX, GDMU2X, GDMU3X
!
   LOGICAL    OOUT1, OOUT2
   INTEGER    KBMX, I, K, CTP, ierr, n, kp1, l, l1, kk, kbi, kmi, km1
   real(kind_phys) tem1, tem2, tem3, cbmfl, mflx_e, teme, tems

   REAL(kind_phys)     HBGT ( IJSDIM )     ! imbalance in column heat
   REAL(kind_phys)     WBGT ( IJSDIM )     ! imbalance in column water
   
   !DDsigma begin local work variables - all on model interfaces (sfc=1)
   REAL(kind_phys)     lamdai( IJSDIM, KMAX+1, nctp )         ! lamda for cloud type ctp
   REAL(kind_phys)     lamdaprod( IJSDIM, KMAX+1   )   ! product of (1+lamda) through cloud type ctp
   REAL(kind_phys)     gdrhom         !  density
   REAL(kind_phys)     gdtvm          !  virtual temperature
   REAL(kind_phys)     gdqm, gdwm,gdlm, gdim   !  water vaper
   REAL(kind_phys)     gdtrm(ntrq:ntr)           ! tracer
   character(len=4) :: cproc  !DDsigmadiag

   ! the following are new arguments to cumup to get them out 
   REAL(kind_phys)     wcv( IJSDIM, KMAX+1, nctp)        ! in-cloud vertical velocity
   REAL(kind_phys)     GCTM  ( IJSDIM, KMAX+1 )   ! cloud T (half lev)   !DDsigmadiag make output
   REAL(kind_phys)     GCQM  ( IJSDIM, KMAX+1, nctp )   ! cloud q (half lev)   !DDsigmadiag make output
   REAL(kind_phys)     GCwM  ( IJSDIM, KMAX+1, nctp )   ! cloud q (half lev)   !DDsigmadiag make output
   REAL(kind_phys)     GCiM  ( IJSDIM, KMAX+1 )   ! cloud q (half lev)   !DDsigmadiag make output
   REAL(kind_phys)     GClM  ( IJSDIM, KMAX+1 )   ! cloud q (half lev)   !DDsigmadiag make output
   REAL(kind_phys)     GChM  ( IJSDIM, KMAX+1, nctp )   ! cloud q (half lev)   !DDsigmadiag make output
   REAL(kind_phys)   GCtrM (IJSDIM, KMAX, ntrq:ntr) ! cloud tracer (half lev) !DDsigmadiag make output
      
! these are the fluxes at the interfaces - AW will operate on them
   REAL(kind_phys), dimension(ijsdim,Kmax+1,nctp) :: sfluxtem, qvfluxtem, qlfluxtem, qifluxtem
   REAL(kind_phys), dimension(ijsdim,Kmax+1,ntrq:ntr,nctp) :: trfluxtem  ! tracer
      
   REAL(kind_phys), dimension(ijsdim,Kmax+1) :: dtcondtem, dqcondtem, dtfrztem, dqprectem,dfrzprectem
   REAL(kind_phys), dimension(ijsdim,Kmax) :: dtevap, dqevap, dtmelt, dtsubl
   REAL(kind_phys), dimension(ijsdim) :: moistening_aw
   real(kind_phys) rhs_q, rhs_h, sftem, qftem, qlftem, qiftem
   real(kind_phys), dimension(ijsdim,kmax+1) :: gctbl, gcqbl,gcwbl, gcqlbl, gcqibl !DDsigmadiag updraft profiles below cloud Base
   real(kind_phys), dimension(ijsdim,kmax,ntrq:ntr) :: gctrbl    !DDsigmadiag tracer updraft profiles below cloud Base
   real(kind_phys), dimension(ijsdim,kmax+1) :: sigmad
   real(kind_phys) :: fsigma( IJSDIM, KMAX+1 )  ! factor to reduce mass flux terms (1-sigma**2) for AW
   real(kind_phys) :: lamdamax  ! for sorting lamda values
   integer loclamdamax
   real(kind_phys) :: pr_tot, pr_ice, pr_liq
!DDsigma end local work variables
!
!  [INTERNAL PARM]
   REAL(kind_phys) :: WCBMIN = 0._kind_phys       ! min. of updraft velocity at cloud base
!M REAL(kind_phys) :: WCBMAX = 1.4_kind_phys      ! max. of updraft velocity at cloud base
!M wcbas commented by Moorthi since it is not used
!M REAL(kind_phys) :: WCBAS  = 2._kind_phys       ! updraft velocity**2 at cloud base (ASMODE)
!M REAL(kind_phys) :: ERAMIN = 1.e-5_kind_phys    ! min. of entrainment rate
                                    ! used only in OPT_ASMODE
!M REAL(kind_phys) :: ERAMAX = 2.e-3_kind_phys    ! max. of entrainment rate
                                    ! used only in OPT_ASMODE
! downdraft mass flux terms now slot nctp+1 in the *fluxterm arrays
    REAL(kind_phys)     dtdwn ( IJSDIM, KMAX   ) ! t tendency downdraft detrainment
    REAL(kind_phys)     dqvdwn ( IJSDIM, KMAX   ) ! qv tendency downdraft detrainment
    REAL(kind_phys)     dqldwn ( IJSDIM, KMAX   ) ! ql tendency downdraft detrainment
    REAL(kind_phys)     dqidwn ( IJSDIM, KMAX   ) ! qi tendency downdraft detrainment
    REAL(kind_phys), dimension(ijsdim,kmax,ntrq:ntr) :: dtrdwn    ! tracer tendency downdraft detrainment

    LOGICAL  :: OINICB = .false.     ! set 0.d0 to CBMFX

    REAL(kind_phys) :: VARMIN = 1.e-13_kind_phys   ! minimum of PDF variance
    REAL(kind_phys) :: VARMAX = 5.e-7_kind_phys    ! maximum of PDF variance
    REAL(kind_phys) :: SKWMAX = 0.566_kind_phys    ! maximum of PDF skewness

    REAL(kind_phys) :: PSTRMX = 400.e2_kind_phys   ! max P of tropopause
    REAL(kind_phys) :: PSTRMN = 50.e2_kind_phys    ! min P of tropopause
    REAL(kind_phys) :: GCRSTR = 1.e-4_kind_phys    ! crit. dT/dz tropopause

        ! 0: mass fixer is not applied
        !    tracers which may become negative values
        !    e.g. subgrid-PDFs
        ! 1: mass fixer is applied, total mass may change through cumulus scheme
        !    e.g. moisture, liquid cloud, ice cloud, aerosols
        ! 2: mass fixer is applied, total mass never change through cumulus scheme
        !    e.g. CO2
    real(kind=kind_phys), parameter  :: zero=0.0, one=1.0
    real(kind=kind_phys)             :: tem, esat
!
    LOGICAL, SAVE :: OFIRST = .TRUE.   ! called first time?
!
   IF (OFIRST) THEN
     OFIRST = .FALSE.
     IF (OINICB) THEN
       CBMFX = zero
     ENDIF
   ENDIF
!

   kp1 = kmax + 1
   do n=1,ntr
     do k=1,kmax
       do i=1,ijsdim
         gtq(i,k,n) = zero
       enddo
     enddo
   enddo
   do k=1,kmax+1
     do i=1,ijsdim
       gmflx(i,k)   = zero
       gmfx0(i,k)   = zero
     enddo
   enddo
   do k=1,kmax
     do i=1,ijsdim
       gtt(i,k)     = zero
       gtu(i,k)     = zero
       gtv(i,k)     = zero
       gprci(i,k)   = zero
       gsnwi(i,k)   = zero
       qliq(i,k)    = zero
       qice(i,k)    = zero
!      gtcfrc(i,k)  = zero
!      cumclw(i,k)  = zero
!      fliqc(i,k)   = zero
       fliqou(i,k)  = zero
       gprcpf(i,k)  = zero
       gsnwpf(i,k)  = zero
       cmdet(i,k)   = zero
     enddo
   enddo
   if (flx_form) then
     do ctp = 1,nctp
       do k=1,kp1
         do i=1,ijsdim
           sfluxtem(i,k,ctp)   = zero
           qvfluxtem(i,k,ctp)  = zero
           qlfluxtem(i,k,ctp)  = zero
           qifluxtem(i,k,ctp)  = zero
         enddo
       enddo
       do n = ntrq,ntr
         do k=1,kp1
           do i=1,ijsdim
            trfluxtem(i,k,n,ctp) = zero
           enddo
         enddo
       enddo
     enddo
       do k=1,kmax
         do i=1,ijsdim
           condtermt(i,k)   = zero
           condtermq(i,k)   = zero
           frzterm(i,k)     = zero
           prectermq(i,k)   = zero
           prectermfrz(i,k) = zero
         enddo
       enddo
      do k=1,kmax
        do i=1,ijsdim
          dtdwn(i,k)       = zero
          dqvdwn(i,k)      = zero
          dqldwn(i,k)      = zero
          dqidwn(i,k)      = zero
        enddo
      enddo
      do n = ntrq,ntr
        do k=1,kmax
          do i=1,ijsdim
            dtrdwn(i,k,n)     = zero
          enddo
        enddo
      enddo
   endif
   do i=1,ijsdim
!    gprcc(i,:)   = zero
!    gmflx(i,kp1) = zero
     gmfx0(i,kp1) = zero
     gtprc0(i)    = zero
!    hbgt(i)      = zero
!    wbgt(i)      = zero
     gdztr(i)     = zero
     kstrt(i)     = kmax
   enddo

   do k=1,kmax
     do i=1,ijsdim
       GDW(i,k)  = GDQ(i,k,1) + GDQ(i,k,ITL) + GDQ(i,k,iti)
     enddo
   enddo
!> -# Compute layer saturate moisture \f$Q_i\f$(GDQS) and 
!! saturate moist static energy (GDHS; see Appendix B in
!! Chikira and Sugiyama (2010) \cite Chikira_2010)
   DO K=1,KMAX
     DO I=ISTS,IENS
       esat      = min(gdp(i,k), fpvs(gdt(i,k)))
       GDQS(I,K) = min(EPSV*esat/max(gdp(i,k)+epsvm1*esat, 1.0e-10), 0.1)
       tem       = one / GDT(I,K)
       FDQS(I,K) = GDQS(I,K) * tem * (fact1 + fact2*tem) ! calculate d(qs)/dT
       GAM (I,K) = ELOCP*FDQS(I,K)
       GDS (I,K) = CP*GDT(I,K) + GRAV*GDZ(I,K) ! layer dry static energy
       GDH (I,K) = GDS(I,K) + EL*GDQ(I,K,1)    ! layer moist static energy
       GDHS(I,K) = GDS(I,K) + EL*GDQS(I,K)     ! layer sat. moist static energy
     ENDDO
   ENDDO
!
!        < tropopause >
!
!> -# Compute tropopause height (GDZTR)
   DO K=1,KMAX
     DO I=ISTS,IENS
       GAMX = (GDTM(I,K+1)-GDTM(I,K)) / (GDZM(I,K+1)-GDZM(I,K))
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

!> -# Call cumbas() to compute cloud base properties
   CALL CUMBAS(IJSDIM, KMAX  ,                           & !DD dimensions
               KB    , GCYM(:,:,1)  , KBMX  ,            & ! output
               ntr   , ntrq  ,                           &
               GCHB  , GCWB  , GCUB  , GCVB  ,           & ! output
               GCIB  , gctrb,                            & ! output
               GDH   , GDW   , GDHS  , GDQS  ,           & ! input
               GDQ(:,:,iti)  , GDU   , GDV   , GDZM  ,   & ! input
               GDPM  , FDQS  , GAM   ,                   & ! input
               lprnt,  ipr,                              &
               ISTS  , IENS                  ,           & !)   ! input
               gctbl, gcqbl,gdq,gcwbl, gcqlbl, gcqibl, gctrbl) ! sub cloud tendencies
!
!> -# Compute CAPE and CIN
!
     DO I=ISTS,IENS
       CAPE(i)  = zero
       CIN(i)   = zero
       JBUOY(i) = 0
     enddo
     DO K=2,KMAX
       DO I=ISTS,IENS
         if (kb(i) > 0) then
           IF (K >= KB(I)) THEN
              BUOY = (GDH(I,1)-GDHS(I,K)) / ((one+ELOCP*FDQS(I,K)) * CP*GDT(I,K))
           ELSE
              BUOY = (GDS(I,1)-GDS(I,K)) / (CP*GDT(I,K))
           END IF
           IF (BUOY > zero .AND. JBUOY(I) >=  -1) THEN
              CAPE(I) = CAPE(I) + BUOY * GRAV * (GDZM(I,K+1) - GDZM(I,K))
              JBUOY(I) = 2
           ELSEIF (BUOY < zero .AND. JBUOY(I) /= 2) THEN
              CIN(I) = CIN(I) + BUOY * GRAV * (GDZM(I,K+1) - GDZM(I,K))
              JBUOY(I) = -1
           ENDIF
         endif
       ENDDO
     ENDDO
     DO I=ISTS,IENS
       IF (JBUOY(I) /= 2) CIN(I) = -999.D0
       if (cin(i) < cincrit) kb(i) = -1
     ENDDO

!DDsigma some initialization  before summing over cloud type
!> -# Initialize variables before summing over cloud types
   if(flx_form) then
   do k=1,kp1    ! Moorthi
     do i=1,ijsdim
       lamdaprod(i,k)   = one
       sigma(i,k) = 0.0
     enddo
   enddo

   do ctp=1,nctp 
     do k=1,kp1
       do i=1,ijsdim
         lamdai(i,k,ctp) = zero
         sigmai(i,k,ctp) = zero
         vverti(i,k,ctp) = zero
       enddo
     enddo
   enddo
   endif

   do ctp=2,nctp
     do k=1,kmax
       do i=1,ijsdim
         gcym(i,k,ctp) = gcym(i,k,1)
       enddo
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


! DH* GNU crashes - check all arguments to CUMUP for their dimensions
! before and after CUMUP (i.e. here), and inside the routine, in
! particular: gctm, gcqm, gcwm, gchm, gcwt, gclm, gcim,gctrm
! also, inside, check that no reads/writes out of bounds occur *DH
!> -# Call cumup() to compute in-cloud properties
     CALL CUMUP(IJSDIM, KMAX, NTR,   ntrq,                          & !DD dimensions
                ACWF        ,                                       & ! output
                GCLZ        , GCIZ        , GPRCIZ(:,:,CTP), GSNWIZ(:,:,CTP),   & ! output
                GCYT(:,CTP) , GCHT(:,CTP) , GCQT (:,CTP),           & ! output
                GCLT(:,CTP) , GCIT(:,CTP) , GTPRT(:,CTP),           & ! output
                GCUT(:,CTP) , GCVT(:,CTP) , gctrt(:,ntrq:ntr,ctp),  & ! output
                KT  (:,CTP) , KTMX(CTP)   ,                         & ! output
                GCYM(:,:,CTP) ,                                     & ! modified
                wcv(:,:,CTP)  ,                                     & ! !DD-sigma new output
                GCHB  , GCWB  , GCUB  , GCVB  ,                     & ! input  !DDsigmadiag
                GCIB  , gctrb ,                                     & ! input
                GDU   , GDV   , GDH   , GDW   ,                     & ! input
                GDHS  , GDQS  , GDT   , GDTM  ,                     & ! input
                GDQ   , GDQ(:,:,iti)  , GDZ   , GDZM  ,             & ! input
                GDPM  , FDQS  , GAM   , GDZTR ,                     & ! input
                CPRES , WCBX  ,                                     & ! input
                KB    , CTP   , ISTS  , IENS  ,                     & ! input
                gctm  , gcqm(:,:,CTP), gcwm(:,:,CTP), gchm(:,:,CTP),&
                gcwt, gclm, gcim, gctrm,                            & ! additional incloud profiles and cloud top total water
                lprnt , ipr )
!
!> -# Call cumbmx() to compute cloud base mass flux
     CALL CUMBMX(IJSDIM, KMAX,                                      & !DD dimensions
                 CBMFX(:,CTP),                                      & ! modified
                 ACWF        , GCYT(:,CTP), GDZM     ,              & ! input
                 GDW         , GDQS       , DELP     ,              & ! input
                 KT   (:,CTP), KTMX(CTP)  , KB       ,              & ! input
                 DELTI       , ISTS       , IENS       )
                 
!DDsigma -  begin sigma computation
! At this point cbmfx is updated and we have everything we need to compute sigma

     if (flx_form) then
       do k=1,kmax + 1    ! Moorthi
         do i=1,ijsdim
           dqcondtem(i,k)   = zero
           dqprectem(i,k)   = zero
           dfrzprectem(i,k) = zero
           dtfrztem(i,k)    = zero
           dtcondtem(i,k)   = zero
         enddo
       enddo

       do i=ISTS,IENS
       cbmfl = cbmfx(i,ctp)
       kk    = kt(i,ctp)      ! cloud top index

       if(cbmfl > zero) then  ! this should avoid zero wcv in the denominator
         kbi = kb(i)          ! cloud base index
         do k=kbi,kk          ! loop from cloud base to cloud top
           km1 = k - 1
! get environment variables interpolated to layer interface
           GDQM   = half * (GDQ(I,K,1) + GDQ(I,KM1,1))  ! as computed in cumup
!          GDwM   = half * (GDw(I,K)   + GDw(I,KM1 ))
           GDlM   = half * (GDQ(I,K,itl) + GDQ(I,KM1,itl))
           GDiM   = half * (GDQ(I,K,iti) + GDQ(I,KM1,iti))
           do n = ntrq,NTR
             GDtrM(n)   = half * (GDQ(I,K,n) + GDQ(I,KM1,n))  ! as computed in cumup
           enddo
           mflx_e = gcym(i,k,ctp) * cbmfl          ! mass flux at level k for cloud ctp


!> -# Compute lamda for a cloud type and then updraft area fraction
!! (sigmai) following Equations 23 and 12 of 
!! Arakawa and Wu (2013) \cite arakawa_and_wu_2013 , respectively

             lamdai(i,k,ctp) = mflx_e * rair * gdtm(i,k)*(one+epsvt*gdqm)         &
                    / (gdpm(i,k)*wcv(i,k,ctp))
                    
! just compute lamdai here, we will compute sigma, sigmai, and vverti outside
!    the cloud type loop after we can sort lamdai
!             lamdaprod(i,k)  = lamdaprod(i,k) * (one+lamdai(i,k,ctp))
!
!!            vverti(i,k,ctp) = wcv(i,k)
!!            sigmai(i,k,ctp) = lamdai / lamdaprod(i,k)
!!            sigma(i,k) = max(zero, min(one, sigma(i,k) + sigmai(i,k,ctp)))
!
!             sigmai(i,k,ctp)          = lamdai(i,k,ctp) / lamdaprod(i,k)
!             sigma(i,k)      = max(zero, min(one, sigma(i,k) + sigmai(i,k,ctp)))
!             vverti(i,k,ctp) = sigmai(i,k,ctp) * wcv(i,k,ctp)


! sigma effect won't be applied until later, when lamda is sorted
!            fsigma     = 1.0   ! no aw effect, comment following lines to undo AW
!             fsigma     = one - sigma(i,k)

!> -# Compute tendencies based on mass flux and condensation
! fsigma is the AW reduction of flux tendencies

             if(k == kbi) then
               do l=2,kbi           ! compute eddy fluxes below cloud base
!                 tem = - fsigma * gcym(i,l,ctp) * cbmfl
                 tem = - gcym(i,l,ctp) * cbmfl

! first get environment variables at layer interface
                 l1 = l - 1
                 GDQM = half * (GDQ(I,l,1) + GDQ(I,l1,1))
                 GDlM = half * (GDQ(I,l,itl) + GDQ(I,l1,itl))
                 GDiM = half * (GDQ(I,l,iti) + GDQ(I,l1,iti))
!!               GDwM  = half * (GDw(I,l)   + GDw(I,l1))
                 do n = ntrq,NTR
                   GDtrM(n) = half * (GDQ(I,l,n) + GDQ(I,l1,n))  ! as computed in cumup
                 enddo

! flux = mass flux * (updraft variable minus environment variable)
!centered differences
                 sfluxtem(i,l,ctp)  = tem * (gdtm(i,l)-gctbl(i,l))
                 qvfluxtem(i,l,ctp) = tem * (gdqm-gcqbl(i,l))
                 qlfluxtem(i,l,ctp) = tem * (gdlm-gcqlbl(i,l))
                 qifluxtem(i,l,ctp) = tem * (gdim-gcqibl(i,l))
                 do n = ntrq,NTR
                   trfluxtem(i,l,n,ctp) = tem * (gdtrm(n)-gctrbl(i,l,n))
                 enddo
!     if(lprnt .and. i == ipr) write(0,*)' l=',l,' kbi=',kbi,' tem =', tem,' trfluxtem=',trfluxtem(l,ntr),&
!     ' gdtrm=',gdtrm(ntr),' gctrbl=',gctrbl(i,l,ntr),' gq=',GDQ(I,l,ntr),GDQ(I,l1,ntr),' l1=',l1,' ctp=',ctp,&
!    ' fsigma=',fsigma,' gcym=',gcym(i,l,ctp),' cbmfl=',cbmfl,' sigma=',sigma(i,k)

!  The following commented out by Moorthi on April 13, 2018 because tke below
!  cloud base becomes too large otherwise when shoc is used

!upstream - This better matches what the original CS tendencies do
!                sfluxtem(l)  = tem * (gdt(i,l)+gocp*(gdz(i,l)-gdzm(i,l))-gctbl(i,l))
!                qvfluxtem(l) = tem * (gdq(i,l,1)-gcqbl(i,l))
!                qlfluxtem(l) = tem * (gdq(i,l,3)-gcqlbl(i,l))
!                qifluxtem(l) = tem * (gdq(i,l,2)-gcqibl(i,l))
!                do n = ntrq,NTR
!                  trfluxtem(l,n)  = tem * (gdq(i,l,n)-gctrbl(i,l,n))
!                enddo

               enddo
             else
! flux = mass flux * (updraft variable minus environment variable)

!              tem = - fsigma * mflx_e
               tem = - mflx_e
!centered
               sfluxtem(i,k,ctp)  = tem * (gdtm(i,k)+gocp*gdzm(i,k)-gctm(i,k))
               qvfluxtem(i,k,ctp) = tem * (gdqm-gcqm(i,k,ctp))
               qlfluxtem(i,k,ctp) = tem * (gdlm-gclm(i,k))
               qifluxtem(i,k,ctp) = tem * (gdim-gcim(i,k))
               do n = ntrq,NTR
                 trfluxtem(i,k,n,ctp) = tem * (gdtrm(n)-gctrm(i,k,n))
               enddo

!upstream  - This better matches what the original CS tendencies do
!              if(k < kk) then
!                sfluxtem(k)  = tem * (gdt(i,k)+gocp*gdz(i,k)-gctm(i,k))
!                qvfluxtem(k) = tem * (gdq(i,k,1)-gcqm(i,k))
!                qlfluxtem(k) = tem * (gdq(i,k,3)-gclm(i,k))
!                qifluxtem(k) = tem * (gdq(i,k,2)-gcim(i,k))
!                do n = ntrq,NTR
!                  trfluxtem(k,n)  = tem * (gdq(i,k,n)-gctrm(i,k,n))
!                enddo
!    if(lprnt .and. i == ipr) write(0,*)' k=',k,' kbi=',kbi,' tem =', tem,' kk=',kk,&
!     ' gctrm=',gctrm(i,k,ntr),' gdq=',gdq(I,k,ntr),' gctrm=',gctrm(I,k,ntr),' ctp=',ctp,&
!    ' fsigma=',fsigma,' mflx_e=',mflx_e,' trfluxtemk=',trfluxtem(k,ntr),' sigma=',sigma(i,k)

!              else
! centered at top of cloud
!                sfluxtem(k)  = tem * (gdtm(i,k)+gocp*gdzm(i,k)-gctm(i,k))
!                qvfluxtem(k) = tem * (gdqm-gcqm(i,k))
!                qlfluxtem(k) = tem * (gdlm-gclm(i,k))
!                qifluxtem(k) = tem * (gdim-gcim(i,k))
!                do n = ntrq,NTR
!                  trfluxtem(k,n)  = tem * (gdtrm(n)-gctrm(i,k,n))
!                enddo
!              endif

!    if(lprnt .and. i == ipr) write(0,*)' k=',k,' kbi=',kbi,' tem =', tem,' kk=',kk,&
!     ' gctrm=',gctrm(i,k,ntr),' gdtrm=',gdtrm(ntr),' gctrm=',gctrm(I,k,ntr),' ctp=',ctp,&
!    ' fsigma=',fsigma,' mflx_e=',mflx_e,' trfluxtemk=',trfluxtem(k,ntr),' sigma=',sigma(i,k)


             endif ! if(k > kbi) then
         enddo     ! end of k=kbi,kk loop

       endif       ! end of if(cbmfl > zero)
    
    
        
     enddo           ! end of i loop
    endif         ! if (flx_form)
!
! we don't reduce these values in AW, just the tendencies due to fluxes
!     do i=ists,iens
!       if (cbmfx(i,ctp) > zero) then
!         tem = one - sigma(i,kt(i,ctp))
!         gcyt(i,ctp)  = tem * gcyt(i,ctp)
!         gtprt(i,ctp) = tem * gtprt(i,ctp)
!         gclt(i,ctp)  = tem * gclt(i,ctp)
!         gcht(i,ctp)  = tem * gcht(i,ctp)
!         gcqt(i,ctp)  = tem * gcqt(i,ctp)
!         gcit(i,ctp)  = tem * gcit(i,ctp)
!         do n = ntrq,ntr
!           gctrt(i,n,ctp)  = tem * gctrt(i,n,ctp)
!         enddo
!         gcut(i,ctp)  = tem * gcut(i,ctp)
!         gcvt(i,ctp)  = tem * gcvt(i,ctp)
!         do k=1,kmax
!           kk = kb(i)         
!           if (k < kk) then
!             tem  = one - sigma(i,kk)
!             tem1 = tem
!           else
!             tem = one - sigma(i,k)
!             tem1 = one - 0.5*(sigma(i,k)+sigma(i,k-1))
!           endif
!           gcym(i,k,ctp) = tem  * gcym(i,k,ctp)
!           gprciz(i,k)   = tem1 * gprciz(i,k)
!           gsnwiz(i,k)   = tem1 * gsnwiz(i,k)
!           gclz(i,k)     = tem1 * gclz(i,k)
!           gciz(i,k)     = tem1 * gciz(i,k)
!         enddo
!       endif
!     enddo

!
!> -# Call cumflx() to compute cloud mass flux and precipitation
     CALL CUMFLX(IM    , IJSDIM, KMAX  ,                               & !DD dimensions
                 GMFX0 , GPRCI , GSNWI , CMDET,                        & ! output
                 QLIQ  , QICE  , GTPRC0,                               & ! output
                 CBMFX(:,CTP)  , GCYM(:,:,ctp), GPRCIZ(:,:,CTP), GSNWIZ(:,:,CTP) ,    & ! input
                 GTPRT(:,CTP)  , GCLZ         , GCIZ     , GCYT(:,ctp),& ! input
                 KB            , KT(:,CTP)    , KTMX(CTP) ,            & ! input
                 ISTS          , IENS                               )    ! input

   ENDDO      ! end of cloud type ctp loop
   
!> -# Compute net updraft mass flux for all clouds
   do k=1,kmax
     do i=ists,iens
       GMFLX(I,k) = GMFX0(I,k) ! contains net updraft mass flux for all clouds
     enddo
   enddo
   KTMXT = 3
   DO CTP=1,NCTP
     IF (KTMX(CTP) > KTMXT) KTMXT = KTMX(CTP)
   ENDDO

!  DO K=1,KTMXT
!    DO I=ISTS,IENS
!      CUMCLW(I,K) = QLIQ(I,K) + QICE(I,K)
!      IF (CUMCLW(I,K) > zero) THEN
!           FLIQC(I,K)  = QLIQ(I,K) / CUMCLW(I,K)
!      ENDIF
!    ENDDO
!  ENDDO
!
! Cumulus Cloudiness
!  CALL CUMCLD(IJSDIM, KMAX  ,                                & !DD dimensions
!              CUMCLW, QLIQ  , QICE  , FLIQC  ,               & ! modified
!              CUMFRC,                                        & ! output
!              GMFLX , KTMXT , ISTS  , IENS    )                ! input
!
!  - Call cumdet() to compute cloud detrainment heating.
   if (.not. flx_form) then
     CALL CUMDET(im    , IJSDIM, KMAX  , NTR   , ntrq  ,      & !DD dimensions
                 GTT   , GTQ   ,         GTU   , GTV   ,      & ! modified
                 GDH   , GDQ   ,         GDU   , GDV   ,      & ! input
!                GTT   , GTQ   , GTCFRC, GTU   , GTV   ,      & ! modified
!                GDH   , GDQ   , GDCFRC, GDU   , GDV   ,      & ! input
                 CBMFX , GCYT  , DELPInv , GCHT  , GCQT  ,      & ! input
                 GCLT  , GCIT  , GCUT  , GCVT  , GDQ(:,:,iti),& ! input
                 gctrt ,                                      &
                 KT    , ISTS  , IENS, nctp              )      ! input
   endif

!for now area fraction of the downdraft is zero, it will be computed
!  within cumdwn and applied there. So we will get the total sigma now before calling it,
!  and apply to the diabatic terms from the updrafts.

!    if (do_aw.and.flx_form) then
    if (flx_form) then
      do k=1,kp1
        do i=ists,iens
           lamdamax = maxval(lamdai(i,k,:))
           do while (lamdamax > zero)
              loclamdamax = maxloc(lamdai(i,k,:),dim=1)
              lamdaprod(i,k)  = lamdaprod(i,k) * (one+lamdai(i,k,loclamdamax))
              sigmai(i,k,loclamdamax)          = lamdai(i,k,loclamdamax) / lamdaprod(i,k)
              sigma(i,k)      = max(zero, min(one, sigma(i,k) + sigmai(i,k,loclamdamax)))
              vverti(i,k,loclamdamax) = sigmai(i,k,loclamdamax) * wcv(i,k,loclamdamax)
              
              ! make this lamdai negative so it won't be counted again
              lamdai(i,k,loclamdamax) = -lamdai(i,k,loclamdamax)
              ! get new lamdamax
              lamdamax =  maxval(lamdai(i,k,:))
           enddo
           ! restore original values of lamdai
           lamdai(i,k,:) = abs(lamdai(i,k,:))
!          write(6,'(i2,14f7.4)') k,sigmai(i,k,:)
        enddo
      enddo
    endif

! the condensation terms - these come from the MSE and condensed water budgets for
!   an entraining updraft
  if(flx_form) then
   DO CTP=1,NCTP                   ! loop over cloud types
       dtcondtem(:,:) = zero
       dqcondtem(:,:) = zero
       dqprectem(:,:) = zero
       dfrzprectem(:,:) = zero
       dtfrztem(:,:) = zero
       do i=ISTS,IENS
       cbmfl = cbmfx(i,ctp)
       kk    = kt(i,ctp)      ! cloud top index
       if(cbmfl > zero) then  ! this should avoid zero wcv in the denominator
         kbi = kb(i)          ! cloud base index
         do k=kbi,kk          ! loop from cloud base to cloud top
           km1 = k - 1
           rhs_h = zero
           rhs_q = zero
             if(k > kbi) then
!                tem   = cbmfl * (one - sigma(i,k))
                 tem   = cbmfl * (one - 0.5*(sigma(i,k)+sigma(i,km1)))
                 tem1  = gcym(i,k,ctp) * (one - sigma(i,k))
                 tem2  = gcym(i,km1,ctp) * (one - sigma(i,km1))
                 rhs_h = cbmfl * (tem1*gchm(i,k,ctp) - (tem2*gchm(i,km1,ctp) &
                                               + GDH(I,Km1)*(tem1-tem2)) )
                 rhs_q = cbmfl * (tem1*(gcwm(i,k,ctp)-gcqm(i,k,ctp))         &
                               - (tem2*(gcwm(i,km1,ctp)-gcqm(i,km1,ctp))     &
                               + (GDw(I,Km1)-gdq(i,km1,1))*(tem1-tem2)) )
!
                 dqcondtem(i,km1)   = -rhs_q                             ! condensation
                 dqprectem(i,km1)   = tem * (GPRCIZ(i,k,ctp) + GSNWIZ(i,k,ctp))  ! total precip production
                 dfrzprectem(i,km1) = tem * GSNWIZ(i,k,ctp)                  ! production of frozen precip
                 dtfrztem(i,km1)    = rhs_h*oneocp                       ! heating due to freezing
! total temperature tendency due to in cloud microphysics
                 dtcondtem(i,km1)   = - elocp * dqcondtem(i,km1) + dtfrztem(i,km1)

             endif ! if(k > kbi) then
         enddo     ! end of k=kbi,kk loop

       endif       ! end of if(cbmfl > zero)
    
    
! get tendencies by difference of fluxes, sum over cloud type

         do k = 1,kk
! sum single cloud microphysical tendencies over all cloud types
           condtermt(i,k)   = condtermt(i,k)   + dtcondtem(i,k)   * delpinv(i,k)
           condtermq(i,k)   = condtermq(i,k)   + dqcondtem(i,k)   * delpinv(i,k)
           prectermq(i,k)   = prectermq(i,k)   + dqprectem(i,k)   * delpinv(i,k)
           prectermfrz(i,k) = prectermfrz(i,k) + dfrzprectem(i,k) * delpinv(i,k)
           frzterm(i,k)     = frzterm(i,k)     + dtfrztem(i,k)    * delpinv(i,k)

!     if (lprnt .and. i == ipr) write(0,*)' k=',k,' trfluxtem=',trfluxtem(k+1,ntr),trfluxtem(k,ntr),&
!       ' ctp=',ctp,' trfluxterm=',trfluxterm(i,k,ntr)
         enddo
        
     enddo           ! end of i loop
   enddo  ! end of nctp loop
  endif
!downdraft sigma and mass-flux tendency terms are now put into 
! the nctp+1 slot of the cloud-type dimensiond variables

    do k=1,kmax
      do i=ists,iens
        sigmad(i,k)  = zero
      enddo
    enddo

!> -# Call cumdwn() to compute cumulus downdraft and assocated melt, freeze 
!! and evaporation
   CALL CUMDWN(IM, IJSDIM, KMAX, NTR, ntrq, nctp,        & ! DD dimensions
               GTT   , GTQ   , GTU   , GTV   ,           & ! modified
                       GMFLX ,                           & ! modified updraft+downdraft flux
               GPRCP , GSNWP , GTEVP , GMDD  ,           & ! output
               GPRCI , GSNWI ,                           & ! input
               GDH   , GDW   , GDQ   , GDQ(:,:,iti) ,    & ! input
               GDQS  , GDS   , GDHS  , GDT   ,           & ! input
               GDU   , GDV   , GDZ   ,                   & ! input
               GDZM  ,         FDQS  , DELP  , DELPInv ,   & ! input
               sigmad, do_aw , do_awdd, flx_form,        & ! DDsigma input
               dtmelt, dtevap, dtsubl,                   & ! DDsigma input
               dtdwn , dqvdwn, dqldwn, dqidwn,           & ! DDsigma input
               dtrdwn,                                   &
               KB    , KTMXT , ISTS  , IENS    )           ! input


!  sigma = sigma + sigmad

!> -# Call cumsbw() to compute cloud subsidence heating
   if (.not. flx_form) then
!  Cloud Subsidence Heating
!  -----------------------=
     CALL CUMSBH(IM    , IJSDIM, KMAX  , NTR   , ntrq  ,   & !DD dimensions
                 GTT   , GTQ   ,                           & ! modified
                 GTU   , GTV   ,                           & ! modified
                 GDH   , GDQ   , GDQ(:,:,iti)  ,           & ! input
                 GDU   , GDV   ,                           & ! input
                 DELPINV , GMFLX , GMFX0 ,                   & ! input
                 KTMXT , CPRES , kb, ISTS  , IENS )   ! input
   else
     CALL CUMSBW(IM    , IJSDIM, KMAX  ,                   & !DD dimensions
                 GTU   , GTV   ,                           & ! modified
                 GDU   , GDV   ,                           & ! input
                 DELPINV , GMFLX , GMFX0 ,                   & ! input
                 KTMXT , CPRES , kb, ISTS  , IENS )          ! input

   endif
!
! for now the following routines appear to be of no consequence - DD
!
   if (.not. flx_form) then
! Tracer Updraft properties
!  -------------
     allocate (gprcc(ijsdim,ntr))
     do n=1,ntr
       do i=1,ijsdim
         gprcc(i,n) = zero
       enddo
     enddo
     CALL CUMUPR(im    , IJSDIM, KMAX  , NTR   ,           & !DD dimensions
                 GTQ   , GPRCC ,                           & ! modified
                 GDQ   , CBMFX ,                           & ! input
                 GCYM  , GCYT  , GCQT  , GCLT  , GCIT  ,   & ! input
                 GTPRT , GTEVP , GTPRC0,                   & ! input
                 KB    , KBMX  , KT    , KTMX  , KTMXT ,   & ! input
                 DELPInv , OTSPT1, ISTS  , IENS,             & ! input
                 fscav , fswtr, nctp)
!
! Tracer Change due to Downdraft
!  ---------------
     CALL CUMDNR(im    ,IJSDIM , KMAX  , NTR   ,           & !DD dimensions
                 GTQ   ,                                   & ! modified
                 GDQ   , GMDD  , DELPInv ,                   & ! input
                 KTMXT , OTSPT1, ISTS  , IENS )              ! input
!!
!! Tracer change due to Subsidence
!! ---------------
!! This will be done by cumsbh, now DD 20170907
!    CALL CUMSBR(im    , IJSDIM, KMAX  , NTR  ,NCTP,       & !DD dimensions
!                GTQ   ,                                   & ! modified
!                GDQ   , DELPI ,                           & ! input
!                GMFLX , KTMXT , OTSPT2,                   & ! input
!                ISTS  , IENS            )                   ! input

   endif   

! if this tracer not advected zero it out
   DO n = ntrq,NTR
     if (.not. OTSPT2(n)) then
       DO K=1,KMAX
         DO I=ISTS,IENS
           gtq(i,k,n) = 0.0
         ENDDO
       ENDDO
     endif
   ENDDO
     
!  if(do_aw .and. flx_form) then ! compute AW tendencies
!> -# Compute AW tendencies of T, ql and qi
   if(flx_form) then ! compute AW tendencies
                                 ! AW lump all heating together, compute qv term
                                 
! sigma interpolated to the layer for condensation, etc. terms, precipitation
     if(do_aw) then
       do k=1,kmax
         kp1 = k+1
         do i=1,ijsdim
           fsigma(i,k)     = one - half*(sigma(i,k)+sigma(i,kp1))
         enddo
       enddo
     else
       do k=1,kmax+1
         do i=1,ijsdim
           fsigma(i,k)     = one
         enddo
       enddo
     endif

! AW adjustment of precip fluxes from downdraft model
     if(do_aw) then
       kp1 = kmax+1
       DO I=ISTS,IENS
         GSNWP( I,kp1 ) = zero
         GPRCP( I,kp1 ) = zero
       ENDDO
       tem1 = cpoemelt/grav
       tem2 = cpoel/grav
       tem3 = cpoesub/grav
       DO K=KMAX,1,-1
         kp1 = k+1
         DO I=ISTS,IENS
           tem = -dtmelt(i,k) * delp(i,k) * tem1
           teme = -dtevap(i,k) * delp(i,k) * tem2
           tems = -dtsubl(i,k) * delp(i,k) * tem3
           GSNWP(I,k) = GSNWP(I,kp1) + fsigma(i,k) * (GSNWI(i,k) - tem - tems)
           GPRCP(I,k) = GPRCP(I,kp1) + fsigma(i,k) * (GPRCI(i,k) + tem - teme) 
         ENDDO
       ENDDO
     endif


! some of the above routines have set the tendencies and they need to be 
!    reinitialized, gtt not needed, but gtq needed Anning 5/25/2020
     do n=1,ntr
       do k=1,kmax
         do i=1,ijsdim
           gtq(i,k,n) = zero
         enddo
       enddo
     enddo
!    do k=1,kmax
!      do i=1,ijsdim
!        gtt(i,k)     = zero
!      enddo
!    enddo
     do k=1,kmax
       do i=ists,iens
         dqevap(i,k) = - dtevap(i,k)*cpoel - dtsubl(i,k)*cpoesub
         dtevap(i,k) =   dtevap(i,k) + dtsubl(i,k)
         dtsubl(i,k) = zero
       enddo
     enddo


! diabatic terms from updraft and downdraft models          
     DO K=1,KMAX
       DO I=ISTS,IENS
         tem = frzterm(i,k)*cpoEMELT - prectermfrz(i,k)
!        gtt(i,k)         = gtt(i,k) + fsigma(i,k)*(dtmelt(i,k) + dtevap(i,k)) + condtermt(i,k)
!        gtq(i,k,1)       = gtq(i,k,1) + fsigma(i,k)*dqevap(i,k) + condtermq(i,k)
!        gtq(i,k,itl)     = gtq(i,k,itl) -  (condtermq(i,k)  + prectermq(i,k) + tem) 
!        gtq(i,k,iti)     = gtq(i,k,iti) + tem 
         gtt(i,k)         = dtdwn(i,k)  + condtermt(i,k)         &
                          + fsigma(i,k)*(dtmelt(i,k) + dtevap(i,k))
         gtq(i,k,1)       = dqvdwn(i,k) + condtermq(i,k)         &
                          + fsigma(i,k) * dqevap(i,k)
         gtq(i,k,itl)     = dqldwn(i,k)  - condtermq(i,k)         &
                          - prectermq(i,k) - tem
         gtq(i,k,iti)     = dqidwn(i,k) + tem


! detrainment terms get zeroed
!        gtldet(i,k)      = zero
!        gtidet(i,k)      = zero
       ENDDO
     ENDDO
!! flux tendencies - compute the vertical flux divergence
     DO ctp =1,nctp
       DO I=ISTS,IENS
         cbmfl = cbmfx(i,ctp)
         kk    = kt(i,ctp)      ! cloud top index
         if(cbmfl > zero) then  ! this should avoid zero wcv in the denominator
           DO K=1,kk
             kp1 = k+1
             gtt(i,k) = gtt(i,k) - (fsigma(i,kp1)*sfluxtem(i,kp1,ctp)   &
                                        - fsigma(i,k)*sfluxtem(i,k,ctp))  * delpinv(i,k)      
             gtq(i,k,1) = gtq(i,k,1) - (fsigma(i,kp1)*qvfluxtem(i,kp1,ctp)   &
                                        - fsigma(i,k)*qvfluxtem(i,k,ctp))  * delpinv(i,k)         
             gtq(i,k,itl) = gtq(i,k,itl) - (fsigma(i,kp1)*qlfluxtem(i,kp1,ctp)   &
                                        - fsigma(i,k)*qlfluxtem(i,k,ctp))  * delpinv(i,k) 
             gtq(i,k,iti) = gtq(i,k,iti) - (fsigma(i,kp1)*qifluxtem(i,kp1,ctp)   &
                                        - fsigma(i,k)*qifluxtem(i,k,ctp))  * delpinv(i,k)
           ENDDO
! replace tracer tendency only if to be advected.
           DO n = ntrq,NTR
             if (OTSPT2(n)) then
               DO K=1,kk
                 kp1 = k+1
                 gtq(i,k,n) = - (fsigma(i,kp1)*trfluxtem(i,kp1,n,ctp)   &
                                        - fsigma(i,k)*trfluxtem(i,k,n,ctp))  * delpinv(i,k)
               ENDDO
             endif
           ENDDO
         end if
       ENDDO
     ENDDO

!      if(kdt>4) stop 1000
     DO I=ISTS,IENS
       moistening_aw(i) = zero
     enddo

! adjust tendencies that will lead to negative water mixing ratios
     tem2 = one / delta
     DO K=1,KMAX
       DO I=ISTS,IENS
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
         tem1             = - gdq(i,k,1)*tem2
         if (gtq(i,k,1) < tem1) then
           gtt(i,k)       = gtt(i,k) + elocp*(gtq(i,k,1)-tem1)
           gtq(i,k,1)     = tem1
         endif
           
! column-integrated total water tendency - to be used to impose water conservation
         moistening_aw(i) = moistening_aw(i)                                       &
                          + (gtq(i,k,1)+gtq(i,k,itl)+gtq(i,k,iti)) * delp(i,k) * gravi
       ENDDO
     ENDDO

! replace tracer tendency only if to be advected.
     DO n = ntrq,NTR
       if (OTSPT2(n)) then
         DO K=1,KMAX
           DO I=ISTS,IENS
             gtq(i,k,n) = gtq(i,k,n) + dtrdwn(i,k,n)
           ENDDO
         ENDDO
       endif
     ENDDO
!    if (lprnt) write(0,*)' endcs_cum gtq=',gtq(ipr,1:25,ntr)
!    if (lprnt) write(0,*)' endcs_cum trfluxterm=',trfluxterm(ipr,1:25,ntr)

   endif        ! if (flx_form)

!!!! this section may need adjustment for cloud ice and water with flux_form
!
!  do k=1,kmax
!    do i=ISTS,IENS
!      GTQ(I,k,ITI) = GTQI(I,k)
!    enddo
!  enddo
!
!> -# Call cumfxr() for tracer mass fixer without detrainment
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
!!!!! end fixer section

!  DO K=1,KMAX
!    DO I=ISTS, IENS
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
!      tem = DELP(I,K)*GRAVI
!      HBGT(I) = HBGT(I) + (CP*GTT(I,K) + EL*GTQ(I,K,1)                         &
!                          - EMELT*GTQ(I,K,ITI)) * tem
!                          - EMELT*(GTQ(I,K,ITI)+GTIDET(I,K))) * tem
!      WBGT(I) = WBGT(I) + (GTQ(I,K,1)   + GTQ(I,K,ITL) + GTQ(I,K,ITI)) * tem 
!                                        + GTLDET(I,K)  + GTIDET(I,K)) * tem
!    ENDDO
!  ENDDO
  

!
!  DO I=ISTS,IENS
!    HBGT(I)  = HBGT(I) - EMELT*GSNWC(I)
!    WBGT(I)  = WBGT(I) + GPRCC(I,1) + GSNWC(I)
!    CTOPP(I) = 1.D6
!  ENDDO
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
!> -# Ensures conservation of water. 
!In fact, no adjustment of the precip
!   is occuring now which is a good sign! DD
   if(flx_form) then
     DO I = ISTS, IENS
       if(gprcp(i,1)+gsnwp(i,1) > 1.e-12_kind_phys) then
         moistening_aw(i) = -moistening_aw(i) / (gprcp(i,1)+gsnwp(i,1))
!        print*,'moistening_aw',moistening_aw(i)
         gprcp(i,:) = gprcp(i,:) * moistening_aw(i)
         gsnwp(i,:) = gsnwp(i,:) * moistening_aw(i)
       endif
     END DO
   endif
   
! second method of determining sfc precip only
! if(flx_form) then
!    DO I = ISTS, IENS
!      pr_tot = zero
!      pr_liq = zero
!      pr_ice = zero
!      do k = 1,kmax
!        pr_tot = pr_tot - (gtq(i,k,1)+gtq(i,k,itl)+gtq(i,k,iti)) * delp(i,k) * gravi
!        pr_ice = pr_ice  + ( cp*gtt(i,k) + el*gtq(i,k,1) - emelt*gtq(i,k,iti) )   &
!                    * delp(i,k)*gravi
!      enddo
       !pr_ice = max( min(pr_tot, pr_ice / (emelt)),zero)
!      pr_ice = pr_ice / emelt
!      pr_liq = pr_tot - pr_ice
!    END DO
!    print *,'precip1',pr_tot*86400.,gprcp(1,1)*86400.,gsnwp(1,1)*86400.
!    print *,'precip2',pr_tot*86400.,pr_liq*86400.,pr_ice*86400.
! endif

   DO K = 1, KMAX
     DO I = ISTS, IENS
       GPRCPF( I,K ) = 0.5*( GPRCP( I,K )+GPRCP( I,K+1 ) )
       GSNWPF( I,K ) = 0.5*( GSNWP( I,K )+GSNWP( I,K+1 ) )
     END DO
   END DO

!
!   do i=ISTS,IENS
!      GPRCC( I,1 ) = GPRCP( I,1 )
!      GSNWC( I   ) = GSNWP( I,1 )
!   enddo

!  adjust sfc precip consistently with vertically integrated
!     temperature and moisture tendencies

   do k=1,kmax+1
     do i=ISTS,IENS
       GTPRP(I,k) = GPRCP(I,k) + GSNWP(I,k)
     enddo
   enddo
!
!DD provide GFS with a separate downdraft mass flux
     if(do_aw) then
     DO K = 1, KMAX+1
        DO I = ISTS, IENS           
           fsigma(i,k)     = one - sigma(i,k)
           GMFX0( I,K ) = GMFX0( I,K ) * fsigma(i,k)
           GMFLX( I,K ) = GMFLX( I,K ) * fsigma(i,k)
        END DO
     END DO
     endif
     DO K = 1, KMAX+1
        DO I = ISTS, IENS           
           GMFX1( I,K ) = GMFX0( I,K ) - GMFLX( I,K )
        END DO
     END DO
     
   if (allocated(gprcc)) deallocate(gprcc)
     
!
      END SUBROUTINE CS_CUMLUS
!> @}
!***********************************************************************
!>\ingroup cs_scheme
!! This subroutine calculates cloud base properties.
      SUBROUTINE CUMBAS                            & ! cloud base
               ( IJSDIM, KMAX  ,                   & !DD dimensions
                 KB    , GCYM  , KBMX  ,           & ! output
                 ntr   , ntrq  ,                   &
                 GCHB  , GCWB  , GCUB  , GCVB  ,   & ! output
                 GCIB  , gctrb,                    & ! output
                 GDH   , GDW   , GDHS  , GDQS  ,   & ! input
                 GDQI  , GDU   , GDV   , GDZM  ,   & ! input
                 GDPM  , FDQS  , GAM   ,           & ! input
                 lprnt,  ipr,                      &
                 ISTS  , IENS , gctbl, gcqbl ,gdq, &
                 gcwbl, gcqlbl, gcqibl, gctrbl   )   ! input  !DDsigmadiag add updraft profiles below cloud base
!
!
      IMPLICIT NONE
!     integer, parameter  :: crtrh=0.80
      integer, parameter  :: crtrh=0.70
      INTEGER, INTENT(IN) :: IJSDIM, KMAX , ntr, ntrq  ! DD, for GFS, pass in
      integer  ipr
      logical  lprnt
!
!   [OUTPUT]
      INTEGER    KB    (IJSDIM)         ! cloud base
      REAL(kind_phys)   GCYM  (IJSDIM, KMAX)   ! norm. mass flux (half lev)
      INTEGER    KBMX
      REAL(kind_phys)   GCHB  (IJSDIM)         ! cloud base MSE
      REAL(kind_phys)   GCWB  (IJSDIM)         ! cloud base total water
      REAL(kind_phys)   GCUB  (IJSDIM)         ! cloud base U
      REAL(kind_phys)   GCVB  (IJSDIM)         ! cloud base V
      REAL(kind_phys)   GCIB  (IJSDIM)         ! cloud base ice
      REAL(kind_phys)   GCtrB (IJSDIM,ntrq:ntr)     ! cloud base tracer

!DDsigma added to arglist for AW, subcloud updraft profiles: temperature, water vapor
!                               total water, cloud water, and cloud ice respectively
      REAL(kind_phys), dimension(ijsdim,kmax)     :: gctbl, gcqbl, gcwbl, gcqlbl, gcqibl   !>DDsigmadiag
      REAL(kind_phys), dimension(ijsdim,kmax,ntrq:ntr) :: gctrbl   !DDsigmadiag
!
!   [INPUT]
      REAL(kind_phys)   GDH   (IJSDIM, KMAX)        ! moist static energy
      REAL(kind_phys)   GDW   (IJSDIM, KMAX)        ! total water
      REAL(kind_phys)   GDq   (IJSDIM, KMAX, ntr)   ! water vapor  and tracer
      REAL(kind_phys)   GDHS  (IJSDIM, KMAX)        ! saturate MSE
      REAL(kind_phys)   GDQS  (IJSDIM, KMAX)        ! saturate humidity
      REAL(kind_phys)   GDQI  (IJSDIM, KMAX)        ! cloud ice
      REAL(kind_phys)   GDU   (IJSDIM, KMAX)        ! u-velocity
      REAL(kind_phys)   GDV   (IJSDIM, KMAX)        ! v-velocity
      REAL(kind_phys)   GDZM  (IJSDIM, KMAX+1)      ! Altitude (half lev)
      REAL(kind_phys)   GDPM  (IJSDIM, KMAX+1)      ! pressure (half lev)
      REAL(kind_phys)   FDQS  (IJSDIM, KMAX)
      REAL(kind_phys)   GAM   (IJSDIM, KMAX)
      INTEGER    ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(kind_phys)   CBASE (IJSDIM)              ! one over cloud base height
!     REAL(kind_phys)   CBASEP(IJSDIM)              ! cloud base pressure
      REAL(kind_phys)   DELZ, GAMX, wrk
!     REAL(kind_phys)   DELZ, QSL, GAMX, wrk
!     REAL(kind_phys), dimension(ijsdim,kmax) :: gchbl   !DDsigmadiag
      real(kind_phys), dimension(ijsdim)      :: gcqb, tx1, spbl, qsl
      INTEGER    I, K, kp1, n
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
!     KBMAX  = KMAXM1              ! cloud base max
      KBMAX  = KMAX/2              ! cloud base max
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
!         KB(I) = KBMAX
          KB(I) = -1
          tx1(i) = one / gdpm(i,1)
        ENDDO
        DO K=KLCLB+1,KBMAX-1
          DO I=ISTS,IENS
            GAMX    = FDQS(I,K) / (one+GAM(I,K)) * oneocp
            QSL(i)  = GDQS(I,K) + GAMX * (GDH(I,KLCLB)-GDHS(I,K))
            spbl(i) = one - gdpm(i,k) * tx1(i)
            IF (GDW(I,KLCLB) >= QSL(i) .and. kb(i) < 0              &
                                       .and. spbl(i) >= spblmin) THEN
!             .and. spbl(i) >= spblcrit .and. spbl(i) < spblcrit*10.0) THEN
              KB(I) = K + KBOFS
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
      do i=ists,iens
        tx1(i) = zero
        qsl(i) = zero
      enddo
      do k=1,kbmax
        do i=ists,iens
          if (k < kb(i)) then
            tx1(i) = tx1(i) + gdw(i,k)  * (GDZM(i,k+1)-GDZM(i,k))
            qsl(i) = qsl(i) + gdqs(i,k) * (GDZM(i,k+1)-GDZM(i,k))
          endif
        enddo
      enddo
      do i=ists,iens
        if (qsl(i) > zero) tx1(i) = tx1(i) / qsl(i)
        if (tx1(i) < crtrh) kb(i) = -1
      enddo
          
!
      KBMX = 1
      DO I=ISTS,IENS
        KBMX = MAX(KBMX, KB(I))
        if (kb(i) > 0) then
          CBASE (I) = one / (GDZM(I,KB(I)) - GDZM(I,1))
!         CBASEP(I) = GDPM(I,KB(I))
        endif
      ENDDO
!
      DO K=2,KBMX
        DO I=ISTS,IENS
          IF (K <= KB(I)) THEN
!           GCYM(I,K) = sqrt((GDZM(I,K)-GDZM(I,1))*CBASE(i))
            GCYM(I,K) = (GDZM(I,K)-GDZM(I,1))*CBASE(i)
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
      do n = ntrq,ntr
        DO I=ISTS,IENS
          GCtrB(I,n) = zero
        enddo
      enddo
      do k=1,kmax
        do i=ists,iens
!         GChbl(i,k)  = zero
          gcqbl(i,k)  = zero
          gcqlbl(i,k) = zero
          gcqibl(i,k) = zero
          gctbl(i,k)  = zero
          gcwbl(i,k)  = zero
        enddo
      enddo
!     do n=ntrq,ntr
!       do k=1,kmax
!         do i=ists,iens
!           gtrqbl(i,k,n)  = zero
!         enddo
!       enddo
!     enddo
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
            GCqB(I) = GCqB(I) + DELZ * GDQ (I,K,1)
!           do n = ntrq,ntr
!             GCtrB(I,n) = GCtrB(I,n) + DELZ * GDQ (I,K,n)
!           enddo
! get subcloud profiles to pass out and do AW eddy flux tendencies
!   removing the normalized mass flux weighting
            wrk       = one /  gcym(i,kp1)
!           gchbl(i,kp1)  = gchb(i) * wrk
            gcqbl(i,kp1)  = gcqb(i) * wrk
            gcqibl(i,kp1) = gcib(i) * wrk
            gcwbl(i,kp1)  = gcwb(i) * wrk
            gcqlbl(i,kp1) = gcwbl(i,kp1) - (gcqibl(i,kp1)+gcqbl(i,kp1))
!           gctbl(i,kp1)  = (gchbl(i,kp1) - grav*gdzm(i,kp1) - el*gcqbl(i,kp1)) * oneocp
            gctbl(i,kp1)  = (gchb(i)*wrk - grav*gdzm(i,kp1) - el*gcqbl(i,kp1)) * oneocp
! tracers
            do n=ntrq,ntr
              GCtrB(I,n)      = GCtrB(I,n) + DELZ * GDQ (I,K,n)
              GCtrBl(I,kp1,n) = gctrb(i,n) * wrk
            enddo
          ENDIF
        ENDDO
      ENDDO
!
      END SUBROUTINE CUMBAS
!***********************************************************************
!>\ingroup cs_scheme
!! This subroutine calculates in-cloud properties.
      SUBROUTINE CUMUP                              & !! in-cloud properties
               ( IJSDIM, KMAX  , NTR   , ntrq  ,    & !DD dimensions
                 ACWF  ,                            & ! output
                 GCLZ  , GCIZ  , GPRCIZ, GSNWIZ,    & ! output
                 GCYT  , GCHT  , GCQT  ,            & ! output
                 GCLT  , GCIT  , GTPRT ,            & ! output
                 GCUT  , GCVT  , gctrt ,            & ! output
                 KT    , KTMX  ,                    & ! output
                 GCYM  ,                            & ! modified
                 wcv   ,                            & ! !DDsigma new output
                 GCHB  , GCWB  , GCUB  , GCVB  ,    & ! input   !DDsigmadiag
                 GCIB  , gctrb ,                    & ! input
                 GDU   , GDV   , GDH   , GDW   ,    & ! input
                 GDHS  , GDQS  , GDT   , GDTM  ,    & ! input
                 GDQ   , GDQI  , GDZ   , GDZM  ,    & ! input
                 GDPM  , FDQS  , GAM   , GDZTR ,    & ! input
                 CPRES , WCB   ,                    & ! input
!                CPRES , WCB   , ERMR  ,            & ! input
                 KB    , CTP   , ISTS  , IENS,      & ! input
                 gctm  , gcqm  , gcwm  , gchm, gcwt,&
                 gclm,   gcim  , gctrm , lprnt, ipr )
!
!DD AW the above line of arguments were previously local, and often scalars.
!  Dimensions were added to them to save profiles for each grid point.
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, KMAX, NTR, ipr , ntrq    ! DD, for GFS, pass in
      logical :: lprnt
!
!   [OUTPUT]
      REAL(kind_phys)   ACWF  (IJSDIM)             !< cloud work function
      REAL(kind_phys)   GCLZ  (IJSDIM, KMAX)       !< cloud liquid water*eta
      REAL(kind_phys)   GCIZ  (IJSDIM, KMAX)       !< cloud ice*eta
      REAL(kind_phys)   GPRCIZ(IJSDIM, KMAX+1)     !< rain generation*eta
      REAL(kind_phys)   GSNWIZ(IJSDIM, KMAX+1)     !< snow generation*eta
      REAL(kind_phys)   GCYT  (IJSDIM)             !< norm. mass flux @top
      REAL(kind_phys)   GCHT  (IJSDIM)             !< cloud top MSE*eta
      REAL(kind_phys)   GCQT  (IJSDIM)             !< cloud top moisture*eta
      REAL(kind_phys)   GCLT  (IJSDIM)             !< cloud top liquid water*eta
      REAL(kind_phys)   GCIT  (IJSDIM)             !< cloud top ice*eta
      REAL(kind_phys)   GCtrT (IJSDIM, ntrq:ntr)   !< cloud top tracer*eta
      REAL(kind_phys)   GTPRT (IJSDIM)             !< cloud top (rain+snow)*eta
      REAL(kind_phys)   GCUT  (IJSDIM)             !< cloud top u*eta
      REAL(kind_phys)   GCVT  (IJSDIM)             !< cloud top v*eta
      REAL(kind_phys)   GCwT  (IJSDIM)             !< cloud top v*eta
      INTEGER    KT    (IJSDIM)                    !< cloud top
      INTEGER    KTMX                              !< max of cloud top
      REAL(kind_phys)   WCV   (IJSDIM, KMAX+1)     !< updraft velocity (half lev) !DD sigma make output
!
!   [MODIFIED]
      REAL(kind_phys)   GCYM  (IJSDIM, KMAX)       !< norm. mass flux
!
!   [INPUT]
      REAL(kind_phys)   GCHB  (IJSDIM)             !< cloud base Moist Static Energy
      REAL(kind_phys)   GCWB  (IJSDIM)             !< cloud base total water
      REAL(kind_phys)   GCUB  (IJSDIM)             !< cloud base U
      REAL(kind_phys)   GCVB  (IJSDIM)             !< cloud base V
      REAL(kind_phys)   GCIB  (IJSDIM)             !< cloud base ice
      REAL(kind_phys)   GCtrB  (IJSDIM,ntrq:ntr)   !< cloud base tracers
      REAL(kind_phys)   GDU   (IJSDIM, KMAX)       !< U
      REAL(kind_phys)   GDV   (IJSDIM, KMAX)       !< V
      REAL(kind_phys)   GDH   (IJSDIM, KMAX)       !< moist static energy
      REAL(kind_phys)   GDW   (IJSDIM, KMAX)       !< total water
      REAL(kind_phys)   GDHS  (IJSDIM, KMAX)       !< saturation MSE
      REAL(kind_phys)   GDQS  (IJSDIM, KMAX)       !< saturation q
      REAL(kind_phys)   GDT   (IJSDIM, KMAX)       !< T
      REAL(kind_phys)   GDTM  (IJSDIM, KMAX+1)     !< T (half lev)
      REAL(kind_phys)   GDQ   (IJSDIM, KMAX, NTR)  !< q  !!DDsigmadiag
      REAL(kind_phys)   GDQI  (IJSDIM, KMAX)       !< cloud ice
      REAL(kind_phys)   GDZ   (IJSDIM, KMAX)       !< z
      REAL(kind_phys)   GDZM  (IJSDIM, KMAX+1)     !< z (half lev)
      REAL(kind_phys)   GDPM  (IJSDIM, KMAX+1)     !< p (half lev)
      REAL(kind_phys)   FDQS  (IJSDIM, KMAX)
      REAL(kind_phys)   GAM   (IJSDIM, KMAX)
      REAL(kind_phys)   GDZTR (IJSDIM)             !< tropopause height
      REAL(kind_phys)   CPRES                      !< pres. fac. for cum. fric.
      REAL(kind_phys)   WCB(ijsdim)                !< cloud base updraft velocity**2
!     REAL(kind_phys)   ERMR                       !< entrainment rate (ASMODE)
      INTEGER    KB    (IJSDIM)
      INTEGER    CTP, ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(kind_phys)     ACWFK (IJSDIM,KMAX)      !< cloud work function
      REAL(kind_phys)     ACWFN (IJSDIM,KMAX)      !< negative part of cloud work function
      REAL(kind_phys)     myGCHt                   !< cloud top h *eta (half lev)
      REAL(kind_phys)     GCHMZ (IJSDIM, KMAX)     !< cloud h *eta (half lev)
      REAL(kind_phys)     GCWMZ (IJSDIM, KMAX)     !< cloud Qt*eta (half lev)
      REAL(kind_phys)     GCUMZ (IJSDIM, KMAX)     !< cloud U *eta (half lev)
      REAL(kind_phys)     GCVMZ (IJSDIM, KMAX)     !< cloud V *eta (half lev)
      REAL(kind_phys)     GCqMZ (IJSDIM      )     !< cloud qv*eta (half lev)
      REAL(kind_phys)     GCIMZ (IJSDIM, KMAX)     !< cloud Qi*eta (half lev)
      REAL(kind_phys)     GCtrMZ(IJSDIM, KMAX,ntrq:ntr)!< cloud tracer*eta (half lev)
      REAL(kind_phys)     GTPRMZ(IJSDIM, KMAX)     !< rain+snow *eta (half lev)
!
      REAL(kind_phys)     BUOY  (IJSDIM, KMAX)     !< buoyancy
      REAL(kind_phys)     BUOYM (IJSDIM, KMAX)     !< buoyancy (half lev)
      REAL(kind_phys)     WCM   (IJSDIM      )     !< updraft velocity**2 (half lev)
!     REAL(kind_phys)     WCM   (IJSDIM, KMAX)     !< updraft velocity**2 (half lev)
!DD sigma make output     REAL(kind_phys)     WCV   ( IJSDIM, KMAX+1 )   !! updraft velocity (half lev)
      REAL(kind_phys)     GCY   (IJSDIM, KMAX)     !< norm. mass flux
!     REAL(kind_phys)     ELAR  (IJSDIM, KMAX)     !< entrainment rate
      REAL(kind_phys)     ELAR                     !< entrainment rate at mid layer
!
      REAL(kind_phys)     GCHM  (IJSDIM, KMAX+1)   !< cloud MSE (half lev)
      REAL(kind_phys)     GCWM  (IJSDIM, KMAX+1)   !< cloud Qt  (half lev)  !DDsigmadiag
      REAL(kind_phys)     GCTM  (IJSDIM, KMAX+1)   !< cloud T (half lev)   !DDsigmadiag make output
      REAL(kind_phys)     GCQM  (IJSDIM, KMAX+1)   !< cloud q (half lev)   !DDsigmadiag make output
      REAL(kind_phys)     GCLM  (IJSDIM, KMAX+1)   !< cloud liquid ( half lev)
      REAL(kind_phys)     GCIM  (IJSDIM, KMAX+1)   !< cloud ice (half lev)
      REAL(kind_phys)     GCUM  (IJSDIM, KMAX)     !< cloud U (half lev)
      REAL(kind_phys)     GCVM  (IJSDIM, KMAX)     !< cloud V (half lev)
      REAL(kind_phys)     GCtrM (IJSDIM, KMAX,ntrq:ntr) !< cloud tracer (half lev)
!
      REAL(kind_phys), dimension(IJSDIM) :: WCM_, ELARM1, GDZMKB
      REAL(kind_phys)     GDQSM, GDHSM, GDQM, GDSM, GDCM, FDQSM, GCCM,         &
                   DELZ, ELADZ, DCTM , CPGMI, DELC, FICE, ELARM2,GCCMZ, &
                   PRECR, GTPRIZ, DELZL, GCCT, DCT, WCVX, PRCZH, wrk
      INTEGER      K, I, kk, km1, kp1, n

!     CHARACTER    CTNUM*2
!
!DD#ifdef OPT_CUMBGT
!DD   REAL(kind_phys)     HBGT  (IJSDIM)           ! heat budget
!DD   REAL(kind_phys)     WBGT  (IJSDIM)           ! water budget
!DD   REAL(kind_phys)     PBGT  (IJSDIM)           ! precipitation budget
!DD   REAL(kind_phys)     MBGT  (IJSDIM)           ! mass budget
!DD   REAL(kind_phys)     GTPRX (IJSDIM)           ! (rain+snow)*eta at top
!DD   REAL(kind_phys)     GSNWT (IJSDIM)           ! cloud top snow*eta
!DD   REAL(kind_phys)     HBMX, WBMX, PBMX, MBMX
!DD   SAVE       HBMX, WBMX, PBMX, MBMX
!DD#endif
!
!   [INTERNAL PARAM]

      REAL(kind_phys), parameter  ::  ZTREF  = 1._kind_phys,   ztrefi = one/ztref, &
                               ELAMIN = zero,    ELAMAX = 4.e-3   ! min and max entrainment rate
      REAL(kind_phys) ::  PB      = 1.0_kind_phys
!m    REAL(kind_phys) ::  TAUZ    = 5.0e3_kind_phys
      REAL(kind_phys) ::  TAUZ    = 1.0e4_kind_phys
!m    REAL(kind_phys) ::  ELMD    = 2.4e-3     ! for Neggers and Siebesma (2002)
!m    REAL(kind_phys) ::  ELAMAX  = 5.e-3      ! max. of entrainment rate
!     REAL(kind_phys) ::  WCCRT   = zero
!m    REAL(kind_phys) ::  WCCRT   = 0.01
      REAL(kind_phys) ::  WCCRT   = 1.0e-6_kind_phys, wvcrt=1.0e-3_kind_phys
      REAL(kind_phys) ::  TSICE   = 273.15_kind_phys  ! compatible with macrop_driver
      REAL(kind_phys) ::  TWICE   = 233.15_kind_phys  ! compatible with macrop_driver
      REAL(kind_phys) ::  c1t

!     REAL(kind_phys) ::  wfn_neg = 0.1
      REAL(kind_phys) ::  wfn_neg = 0.15
!     REAL(kind_phys) ::  wfn_neg = 0.25
!     REAL(kind_phys) ::  wfn_neg = 0.30
!     REAL(kind_phys) ::  wfn_neg = 0.35
      
      REAL(kind_phys) ::  esat, tem
!     REAL(kind_phys) ::  esat, tem, rhs_h, rhs_q
!
!   [INTERNAL FUNC]
      REAL(kind_phys)     FPREC   ! precipitation ratio in condensate
      REAL(kind_phys)     FRICE   ! ice ratio in cloud water
      REAL(kind_phys)     Z       ! altitude
      REAL(kind_phys)     ZH      ! scale height
      REAL(kind_phys)     T       ! temperature
!
      FPREC(Z,ZH) = MIN(MAX(one-EXP(-(Z-PRECZ0)/ZH), zero), one)
      FRICE(T)    = MIN(MAX((TSICE-T)/(TSICE-TWICE), zero), one)
!
! Note: iteration is not made to diagnose cloud ice for simplicity
!
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
      do k=1,kmax+1
        do i=ists,iens
          GPRCIZ(I,k) = zero
          GSNWIZ(I,k) = zero
        enddo
      enddo
      do k=1,kmax
        do i=ists,iens
          WCV   (I,k) = unset_kind_phys
          GCLM  (I,k) = unset_kind_phys
          GCIM  (I,k) = unset_kind_phys
        enddo
      enddo
      do k=1,kmax
        do i=ists,iens
          ACWFK (I,k) = unset_kind_phys
          ACWFN (I,k) = unset_kind_phys
          GCLZ  (I,k) = zero
          GCIZ  (I,k) = zero
!
          GCHMZ (I,k) = zero
          GCWMZ (I,k) = zero
          GCIMZ (I,k) = zero
          GCUMZ (I,k) = zero
          GCVMZ (I,k) = zero
          GTPRMZ(I,k) = zero
!
          BUOY  (I,k) = unset_kind_phys
          BUOYM (I,k) = unset_kind_phys
          GCY   (I,k) = unset_kind_phys
!
          GCHM  (I,k) = unset_kind_phys
          GCWM  (I,k) = unset_kind_phys
          GCTM  (I,k) = unset_kind_phys
          GCQM  (I,k) = unset_kind_phys
          GCUM  (I,k) = unset_kind_phys
          GCVM  (I,k) = unset_kind_phys
        enddo
      enddo
      do i=ists,iens
        GCqMZ(I) = zero
        WCM(I)   = unset_kind_phys
        WCM_(I)  = zero
      enddo
!  tracers
      do n=ntrq,ntr
        do i=ists,iens
          GCtrT(I,n) = zero
        enddo
        do k=1,kmax
          do i=ists,iens
            GCTRM(I,k,n) = unset_kind_phys
          enddo
        enddo
      enddo

!     DO I=ISTS,IENS
!       if (kb(i) > 0) then
!         GDZMKB(I) = GDZM(I,KB(I))     ! cloud base height
!       endif
!     ENDDO
!
!     < cloud base properties >
!
      DO I=ISTS,IENS
        K = KB(I)
        if (k > 0) then
          GDZMKB(I) = GDZM(I,K)     ! cloud base height
          GCHM(I,K) = GCHB(I)
          GCWM(I,K) = GCWB(I)
          WCM_(I)   = WCB(i)
          GCUM(I,K) = GCUB(I)
          GCVM(I,K) = GCVB(I)
          do n = ntrq,ntr
            GCtrM(I,K,n) = GCtrB(I,n)
          enddo
!
          esat  = min(gdpm(i,k), fpvs(gdtm(i,k)))
          GDQSM = min(EPSV*esat/max(gdpm(i,k)+epsvm1*esat, 1.0e-10), 0.1)
          gdsm  = CP*GDTM(I,K) + GRAV*GDZMKB(I)        ! dse
          GDHSM = gdsm + EL*GDQSM                      ! saturated mse
!         FDQSM = FDQSAT(GDTM(I,K), GDQSM)
          tem   = one / GDTM(I,K)
          FDQSM = GDQSM * tem * (fact1 + fact2*tem)    ! calculate d(qs)/dT
!
          tem       = one / (CP+EL*FDQSM)
          DCTM      = (GCHB(I) - GDHSM) * tem
          GCQM(I,K) = min(GDQSM + FDQSM*DCTM, GCWM(I,K))
          GCCM      = MAX(GCWM(I,K)-GCQM(I,K), zero)
!         GCTM(I,K) = GDT(I,K) + DCTM                  ! old
!         GCTM(I,K) = (GCHB(I) - gdsm - el*gcqm(i,k)) * oneocp + dctm  ! new
          GCTM(I,K) = (GCHB(I) - grav*gdzm(i,k) - el*gcqm(i,k)) * oneocp + dctm  ! new
!
          GCIM(I,K) = FRICE(GCTM(I,K)) * GCCM          ! cloud base ice
          GCLM(I,K) = MAX(GCCM-GCIM(I,K), zero)        ! cloud base liquid
          GCHM(I,K) = GCHM(I,K) + EMELT * (GCIM(I,K)-GCIB(I))
          DCTM      = (GCHM(I,K) - GDHSM) * tem
!         GCTM(I,K) = dctm + GDT(I,K) + gocp*gdzm(i,k) !DD old AW convert to DSE
          GCTM(I,K) = dctm + (GCHB(I) - el*gcqm(i,k)) * oneocp ! new, make dse
!
          GDQM  = half * (GDQ(I,K,1)     + GDQ(I,K-1,1))
          GDCM  = half * (GDQ(I,K,ITL)   + GDQI(I,K)                &
                       +  GDQ(I,K-1,ITL) + GDQI(I,K-1))
                       
!
          BUOYM(I,K) = (DCTM*tem + EPSVT*(GCQM(I,K)-GDQM) - GCCM + GDCM )*GRAV
!
          ACWFK(I,K) = zero
          ACWFN(I,K) = zero
!
!DD#ifdef OPT_ASMODE
!DD     ELARM1(I) = ERMR
!DD#elif defined OPT_NS02
!DD     ELARM1(I) = ELMD / SQRT(WCM(I,K))
!DD#else
!         ELARM1(I) = CLMD*PA*BUOYM(I,K)/WCM(I,K)
!         ELARM1(I) = min(max(CLMD*PA*BUOYM(I,K)/WCM_(I), ELAMIN), ELAMAX)
!DD#endif
!         ELARM1(I) = MIN(MAX(ELARM1(I), ELAMIN), ELAMAX)
!
          GCHMZ(I,K) = GCHM(I,K)
          GCWMZ(I,K) = GCWM(I,K)
          GCUMZ(I,K) = GCUM(I,K)
          GCVMZ(I,K) = GCVM(I,K)
          GCqMZ(I)   = GCqM(I,K)
          GCIMZ(I,K) = GCIM(I,K)
          do n = ntrq,ntr
            GCtrMZ(I,K,n) = GCtrM(I,K,n)
          enddo
        endif
      ENDDO
!
!     < in-cloud properties >
!
      DO K=3,KMAX
        km1 = k - 1
        DO I=ISTS,IENS
          IF (kb(i) > 0 .and. K > KB(I) .AND. WCM_(I) > WCCRT) THEN
            WCV(I,KM1) = SQRT(MAX(WCM_(I), zero))
            DELZ       = GDZM(I,K) - GDZM(I,KM1)
            ELARM1(I)  = min(max(CLMDPA*BUOYM(I,KM1)/WCM_(I), ELAMIN), ELAMAX)
            GCYM(I,K)  = GCYM(I,KM1) * EXP(ELARM1(I)*DELZ)
            ELADZ      = GCYM(I,K) - GCYM(I,KM1)
!
            GCHMZ(I,K) = GCHMZ(I,KM1) + GDH(I,KM1)*ELADZ
            GCWMZ(I,K) = GCWMZ(I,KM1) + GDW(I,KM1)*ELADZ
!
            esat  = min(gdpm(i,k), fpvs(gdtm(i,k)))
            GDQSM = min(EPSV*esat/max(gdpm(i,k)+epsvm1*esat, 1.0e-10), 0.1)
            GDHSM = CP*GDTM(I,K ) + GRAV*GDZM(I,K) + EL*GDQSM
!           FDQSM = FDQSAT(GDTM(I,K), GDQSM)
            tem   = one / GDTM(I,K)
            FDQSM = GDQSM * tem * (fact1 + fact2*tem) ! calculate d(qs)/dT
            CPGMI = one / (CP + EL*FDQSM)

!
            wrk   = one / GCYM(I,K)
            DCTM          = (GCHMZ(I,K)*wrk - GDHSM) * CPGMI
            GCQMZ(i)      = min((GDQSM+FDQSM*DCTM)*GCYM(I,K), GCWMZ(I,K))
            if(PRECZH > zero) then
              PRCZH = PRECZH * MIN(GDZTR(I)*ZTREFI, one)
              PRECR = FPREC(GDZM(I,K)-GDZMKB(I), PRCZH )
              GTPRMZ(I,K)   = PRECR * (GCWMZ(I,K)-GCQMZ(i))
            else
              DELC=GDZ(I,K)-GDZ(I,KM1)
              if(gdtm(i,k)>TSICE) then
                c1t=c0t*delc
              else
                c1t=c0t*exp(d0t*(gdtm(i,k)-TSICE))*delc
              end if
              c1t=min(c1t, one)
              GTPRMZ(I,K)   = c1t * (GCWMZ(I,K)-GCQMZ(i))
            end if
            GTPRMZ(I,K)   = MAX(GTPRMZ(I,K), GTPRMZ(I,KM1))
            GCCMZ         = GCWMZ(I,K) - GCQMZ(i) - GTPRMZ(I,K )
            DELC          = MIN(GCCMZ, zero)
            GCCMZ         = GCCMZ    - DELC
            GCQMZ(i)      = GCQMZ(i) + DELC
!
            FICE          = FRICE(GDTM(I,K)+DCTM )
            GCIMZ(I,K)    = FICE * GCCMZ
            GSNWIZ(I,KM1) = FICE * (GTPRMZ(I,K)-GTPRMZ(I,KM1))
            GCHMZ(I,K)    = GCHMZ(I,K) + EMELT * (GCIMZ(I,K) + GSNWIZ(I,KM1) &
                                       - GCIMZ(I,KM1) - GDQI(I,KM1)*ELADZ)
            DCTM          = (GCHMZ(I,K)*wrk - GDHSM) * CPGMI
!
            GDQM          = half * (GDQ(I,K,1)     + GDQ(I,KM1,1))
            GDCM          = half * (GDQ(I,K,ITL)   + GDQI(I,K)          &
                                 +  GDQ(I,KM1,ITL) + GDQI(I,KM1))
            GCQM(I,K)     = wrk * GCQMZ(i)
            GCCM          = wrk * GCCMZ
!
            BUOYM(I,K)    = (DCTM*tem + EPSVT*(GCQM(I,K)-GDQM)-GCCM+GDCM) * GRAV
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
              WCM(I) = (WCM_(I) + CLMP*DELZ*BUOY(I,KM1)) / (one + DELZ/TAUZ)
            ELSE
              WCM(I) = (WCM_(I) + PA*(DELZ+DELZ)*BUOY(I,KM1) ) &
                     / (one + DELZ/TAUZ + (DELZ+DELZ)*ELAMIN )
            ENDIF
!DD#endif
!
!DD#ifdef OPT_ASMODE
!DD         ELARM2 = ERMR
!DD#elif OPT_NS02
!DD         ELARM2 = ELMD/SQRT(MAX(WCM(I,K), EPSln))
!DD#else
!           ELARM2        = CLMD*PA*BUOYM(I,K) / MAX(WCM(I), EPSln)
!DD#endif
            if (WCM(I) > zero) then
              ELARM2 = min(max(CLMDPA*BUOYM(I,K)/WCM(I),ELAMIN), ELAMAX)
            else
              ELARM2 = zero
            endif
            ELAR       = half * (ELARM1(I) + ELARM2)
            GCYM(I,K)  = GCYM(I,KM1) * EXP(ELAR*DELZ)
            ELADZ      = GCYM(I,K) - GCYM(I,KM1)
!
            GCHMZ(I,K) = GCHMZ(I,KM1) + GDH(I,KM1)*ELADZ
            GCWMZ(I,K) = GCWMZ(I,KM1) + GDW(I,KM1)*ELADZ
            GCUMZ(I,K) = GCUMZ(I,KM1) + GDU(I,KM1)*ELADZ
            GCVMZ(I,K) = GCVMZ(I,KM1) + GDV(I,KM1)*ELADZ
            do n = ntrq,ntr
              GCtrMZ(I,K,n) = GCtrMZ(I,KM1,n) + GDq(I,KM1,n)*ELADZ
            enddo
!
            wrk           = one / GCYM(I,K)
            DCTM          = (GCHMZ(I,K)*wrk - GDHSM) * CPGMI
            GCQMZ(i)      = min((GDQSM+FDQSM*DCTM)*GCYM(I,K), GCWMZ(I,K))
            if(PRECZH > zero) then
              GTPRMZ(I,K)   = PRECR * (GCWMZ(I,K)-GCQMZ(i))
            else
              GTPRMZ(I,K)   = c1t * (GCWMZ(I,K)-GCQMZ(i))
            end if
            GTPRMZ(I,K)   = MAX(GTPRMZ(I,K), GTPRMZ(I,KM1))
            GCCMZ         = GCWMZ(I,K) - GCQMZ(i) - GTPRMZ(I,K)
            DELC          = MIN(GCCMZ, zero)
            GCCMZ         = GCCMZ    - DELC
            GCQMZ(i)      = GCQMZ(i) + DELC
            GCCM          = wrk * GCCMZ
            GCQM(I,K)     = wrk * GCQMZ(i)
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
            GCWM(I,K)     = GCWMZ(I,K) * wrk
            GCUM(I,K)     = GCUMZ(I,K) * wrk
            GCVM(I,K)     = GCVMZ(I,K) * wrk
            do n = ntrq,ntr
             GCtrM(I,K,n) = GCtrMZ(I,K,n) * wrk
            enddo
            DELZL         = GDZ(I,KM1)-GDZM(I,KM1)
            GCY (I,KM1)   = GCYM(I,KM1) * EXP(ELAR*DELZL)
            GCLZ(I,KM1)   = half * (GCLM(I,K) + GCLM(I,KM1)) * GCY(I,KM1)
            GCIZ(I,KM1)   = half * (GCIM(I,K) + GCIM(I,KM1)) * GCY(I,KM1)

!
            BUOYM(I,K)    = (DCTM*tem + EPSVT*(GCQM(I,K)-GDQM)-GCCM+GDCM) * GRAV
            BUOY(I,KM1)   = half * (BUOYM(I,K)+BUOYM(I,KM1))
!
            IF (BUOY(I,KM1) > zero) THEN
              WCM(I)   = (WCM_(I) + CLMP*DELZ*BUOY(I,KM1)) / (one + DELZ/TAUZ)
            ELSE
              WCM(I)   = (WCM_(I) + PA*(DELZ+DELZ)*BUOY(I,KM1) ) &
                       / (one + DELZ/TAUZ + (DELZ+DELZ)*ELAMIN )
            ENDIF

!
!           IF (BUOY(I,KM1) > zero) THEN
!             ACWF(I)     = ACWF(I) + BUOY(I,KM1)*GCY(I,KM1)*DELZ
!           ENDIF
!           ACWF(I)       = ACWF(I) + BUOY(I,KM1)*GCY(I,KM1)*DELZ
!!!         wrk           = BUOY(I,KM1)*GCY(I,KM1)*DELZ
!!!         ACWFK(I,K)    = ACWFK(I,KM1) + wrk
!!!         ACWFN(I,K)    = ACWFN(I,KM1) - min(wrk,0.0)
!           ACWFN(I,K)    = ACWFN(I,KM1) + abs(min(wrk,0.0))
!

            wrk        = BUOY(I,KM1)*GCY(I,KM1)*DELZ
            ACWFK(I,K) = ACWFK(I,KM1) + wrk
            ACWFN(I,K) = ACWFN(I,KM1) - min(wrk,0.0)

            WCM_(I)    = WCM(I)

!      if (lprnt .and. i == ipr) write(0,*) ' in cumup k=',k,' km1=',km1,' WCM_=',WCM_(I),' gcy=',gcy(i,km1),' buoym=',buoym(i,km1)

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
          if (kb(i) > 0 .and.  k > kb(i) .and. ACWFK(I,K) > 1.0e-10) then
            wrk = ACWFN(I,K) / ACWFK(I,K)
            IF (KT(I)  == -1 .and. wrk < wfn_neg .AND. WCV(I,K) > WVCRT) THEN
              KT(I)   = K
              ACWF(I) = ACWFK(I,K)
            ENDIF
          endif
        ENDDO
      ENDDO
!     if (lprnt .and. kt(ipr) > 0) write(0,*) ' in cumup kt=',kt(ipr),' gcy=',gcy(ipr,kt(ipr))
!
      KTMX = 2
      DO I=ISTS,IENS
        kt(i) = min(kt(i), kmax-1)
        KTMX  = max(ktmx,  KT(I))
      ENDDO
!
      DO I=ISTS,IENS
        kk = max(1, kt(i)+1)
        do k=kk,kmax
          GCYM  (I,K) = zero
          GCLZ  (I,K) = zero 
          GCIZ  (I,K) = zero
          GPRCIZ(I,K) = zero
          GSNWIZ(I,K) = zero
        enddo
      ENDDO
!     if (lprnt .and. kt(ipr) > 0) write(0,*) ' in cumup2 kt=',kt(ipr),' gcy=',gcy(ipr,kt(ipr))
!
!     < cloud top properties >
!
      DO I=ISTS,IENS
        IF (kb(i) > 0 .and. KT(I) > kb(i)) THEN
          K   = KT(I)
          kp1 = k + 1
          GCYT(I) = GCY(I,K)
          ELADZ   = GCYT(I) - GCYM(I,K)
!
          GCHT(I) = GCHMZ(I,K) + GDH(I,K)*ELADZ
          GCWT(i) = GCWMZ(I,K) + GDW(I,K)*ELADZ
          GCUT(I) = GCUMZ(I,K) + GDU(I,K)*ELADZ
          GCVT(I) = GCVMZ(I,K) + GDV(I,K)*ELADZ
          do n = ntrq,NTR
            GCtrT(I,n) = GCtrMZ(I,K,n) + GDq(I,K,n)*ELADZ
          enddo
!
          wrk         = one / gcyt(i)
          DCT         = (GCHT(I)*wrk - GDHS(I,K)) / (CP*(one + GAM(I,K)))
          GCQT(I)     = min((GDQS(I,K) + FDQS(I,K)*DCT) * GCYT(I), GCWT(i))
          if(PRECZH > zero) then
            PRCZH       = PRECZH * MIN(GDZTR(I)*ZTREFI, one)
            GTPRT(I)    = FPREC(GDZ(I,K)-GDZMKB(I), PRCZH) * (GCWT(i)-GCQT(I))
          else
            DELC=GDZ(I,K)-GDZ(I,K-1)
            if(gdtm(i,k)>TSICE) then
              c1t=c0t*delc
            else
              c1t=c0t*exp(d0t*(gdtm(i,k)-TSICE))*delc
            end if
            c1t=min(c1t, one)
            GTPRT(I)    = c1t * (GCWT(i)-GCQT(I))
          end if
          GTPRT(I)    = MAX(GTPRT(I), GTPRMZ(I,K))
          GCCT        = GCWT(i) - GCQT(I) - GTPRT(I)
          DELC        = MIN(GCCT, zero)
          GCCT        = GCCT    - DELC
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
          do n = ntrq,NTR
!           GCtrT(I,n) = GCtrT(I,n)*(one-CPRES) + GCY(I,K)*GDq(I,K,n)*CPRES
            GCtrT(I,n) = GCtrT(I,n) + GCY(I,K)*GDq(I,K,n)
          enddo
          GCLZ(I,K)   = GCLT(I)
          GCIZ(I,K)   = GCIT(I)

!DD AW get  the cloud top values denormalized and put into profile
          mygcht      = gcht(I) - el*(gcwt(i) - gcqt(i))

          gctm(i,kp1) = wrk * (mygcht - el*gcqt(i)) * oneocp
!Moorthi  gcqm(i,kp1) = gcqt(i)
          gcqm(i,kp1) = gcqt(i)*wrk     ! check this - oct17 2016
          gcim(i,kp1) = gcit(i)*wrk
          gclm(i,kp1) = gclt(i)*wrk
          do n = ntrq,NTR
            gctrm(i,kp1,n) = gctrt(i,n)*wrk
          enddo
!
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
!>\ingroup cs_scheme
!! This subroutine computes cloud base mass flux.
      SUBROUTINE CUMBMX                    & !! cloud base mass flux
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
      REAL(kind_phys)     CBMFX (IJSDIM)          !< cloud base mass flux
!
!   [INPUT]
      REAL(kind_phys)     ACWF  (IJSDIM)          !< cloud work function
      REAL(kind_phys)     GCYT  (IJSDIM)          !< norm mass flux @top
      REAL(kind_phys)     GDZM  (IJSDIM, KMAX+1)  !< height
      REAL(kind_phys)     GDW   (IJSDIM, KMAX)    !< total water
      REAL(kind_phys)     GDQS  (IJSDIM, KMAX)    !< saturate humidity
      REAL(kind_phys)     DELP  (IJSDIM, KMAX)    !< delt pressure
      INTEGER      KT    (IJSDIM)                 !< cloud top
      INTEGER      KTMX                           !< max. of cloud top
      INTEGER      KB    (IJSDIM)                 !< cloud base
      REAL(kind_phys)     DELT                    !< time step
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(kind_phys), dimension(ijsdim) :: QX, QSX, RHM
      INTEGER      I, K
      REAL(kind_phys)     ALP, FMAX1, wrk
!
!   [INTERNAL PARAM]
      REAL(kind_phys) :: FMAX   = 1.5e-2_kind_phys         ! maximum flux
!     REAL(kind_phys) :: RHMCRT = zero                     ! critical val. of cloud mean RH
      REAL(kind_phys) :: RHMCRT = 0.25_kind_phys           ! critical val. of cloud mean RH
!     REAL(kind_phys) :: RHMCRT = 0.50_kind_phys           ! critical val. of cloud mean RH
      REAL(kind_phys) :: ALP1   = zero
      REAL(kind_phys) :: TAUD   = 1.e3_kind_phys
!     REAL(kind_phys) :: TAUD   = 6.e2_kind_phys
      REAL(kind_phys) :: ZFMAX  = 3.5e3_kind_phys
      REAL(kind_phys) :: ZDFMAX = 5.e2_kind_phys
!     REAL(kind_phys) :: FMAXP  = 2._kind_phys
!
      do i=ists,iens
        qx(i)  = zero
        qsx(i) = zero
      enddo
!
      DO K=1,KTMX
        DO I=ISTS,IENS
          IF (kb(i) > 0 .and. K >= KB(I) .AND. K <= KT(I)) THEN
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
        k = kb(i)
        IF (k > 0 .and. KT(I) > K .AND. RHM(I) >= RHMCRT) THEN
          cbmfx(i) = max(cbmfx(i), zero)
          ALP      = ALP0 + ALP1 * (GDZM(I,KT(I))-GDZM(I,K))
          FMAX1    = (one - TANH((GDZM(I,1)-ZFMAX)/ZDFMAX)) * half
!         FMAX1    = FMAX * FMAX1**FMAXP
          FMAX1    = FMAX * FMAX1 * FMAX1
!         CBMFX(I) = CBMFX(I) + MAX(ACWF(I), zero)/(ALP+ALP)*DELT
!         CBMFX(I) = CBMFX(I) / (one + DELT/(TAUD+TAUD))
          CBMFX(I) = (CBMFX(I) + ACWF(I)*(delt/(ALP+ALP))) * wrk
          CBMFX(I) = MIN(max(CBMFX(I), zero), FMAX1/GCYT(I))
        ELSE
          CBMFX(I) = zero
        ENDIF
      ENDDO
!
      END SUBROUTINE CUMBMX
!***********************************************************************
!>\ingroup cs_scheme
!! This subroutine computes cloud mass flux & precip.
      SUBROUTINE CUMFLX                                   & !! cloud mass flux
                      ( IM    , IJSDIM, KMAX  ,           & !DD dimensions
                        GMFLX , GPRCI , GSNWI , CMDET,    & ! output
                        QLIQ  , QICE  , GTPRC0,           & ! output
                        CBMFX , GCYM  , GPRCIZ, GSNWIZ,   & ! input
                        GTPRT , GCLZ  , GCIZ  , GCYT,     & ! input
                        KB    , KT    , KTMX  ,           & ! input
                        ISTS  , IENS                    )   ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, KMAX, IM            !! DD, for GFS, pass in
!
!   [OUTPUT]
      REAL(kind_phys)     GMFLX (IJSDIM, KMAX+1)   !< mass flux
      REAL(kind_phys)     CMDET (IJSDIM, KMAX)     !< detrainment mass flux
      REAL(kind_phys)     GPRCI (IJSDIM, KMAX)     !< rainfall generation
      REAL(kind_phys)     GSNWI (IJSDIM, KMAX)     !< snowfall generation
      REAL(kind_phys)     QLIQ  (IJSDIM, KMAX)     !< cloud liquid
      REAL(kind_phys)     QICE  (IJSDIM, KMAX)     !< cloud ice
      REAL(kind_phys)     GTPRC0(IJSDIM)           !< precip. before evap.
!
!   [INPUT]
      REAL(kind_phys)     CBMFX (IJSDIM)           !< cloud base mass flux
      REAL(kind_phys)     GCYM  (IJSDIM, KMAX)     !< normalized mass flux
      REAL(kind_phys)     GCYT  (IJSDIM)           !< detraining mass flux
      REAL(kind_phys)     GPRCIZ(IJSDIM, KMAX+1)   !< precipitation/M
      REAL(kind_phys)     GSNWIZ(IJSDIM, KMAX+1)   !< snowfall/M
      REAL(kind_phys)     GTPRT (IJSDIM)           !< rain+snow @top
      REAL(kind_phys)     GCLZ  (IJSDIM, KMAX)     !< cloud liquid/M
      REAL(kind_phys)     GCIZ  (IJSDIM, KMAX)     !< cloud ice/M
      INTEGER      KB    (IJSDIM)           !< cloud base
      INTEGER      KT    (IJSDIM)           !< cloud top
      INTEGER      KTMX                     !< max of cloud top
      INTEGER      ISTS, IENS, I, K
!
      DO K=1,KTMX
        DO I=ISTS,IENS
          if (kb(i) > 0) then
            GMFLX(I,K) = GMFLX(I,K) + CBMFX(I) * GCYM(I,K)
            GPRCI(I,K) = GPRCI(I,K) + CBMFX(I) * GPRCIZ(I,K)
            GSNWI(I,K) = GSNWI(I,K) + CBMFX(I) * GSNWIZ(I,K)
            QLIQ(I,K)  = QLIQ (I,K) + CBMFX(I) * GCLZ(I,K)
            QICE(I,K)  = QICE (I,K) + CBMFX(I) * GCIZ(I,K)
          endif
        ENDDO
      ENDDO
!
      DO I= ISTS,IENS
        k = kt(i)
        if (kb(i) > 0 .and. k > kb(i)) then
          GTPRC0(I)  = GTPRC0(I)  + CBMFX(I) * GTPRT(I)
          CMDET(I,K) = CMDET(I,K) + CBMFX(I) * GCYT(I)
        endif
      ENDDO
!
      END SUBROUTINE CUMFLX
!***********************************************************************
!>\ingroup cs_scheme
!! This subroutine calculates cloud detrainment heating.
      SUBROUTINE CUMDET                                    & !! detrainment
               ( im    , IJSDIM, KMAX  , NTR   , ntrq  ,   & !DD dimensions
                 GTT   , GTQ   ,         GTU   , GTV   ,   & ! modified
                 GDH   , GDQ   ,         GDU   , GDV   ,   & ! input
!                GTT   , GTQ   , GTCFRC, GTU   , GTV   ,   & ! modified
!                GDH   , GDQ   , GDCFRC, GDU   , GDV   ,   & ! input
                 CBMFX , GCYT  , DELPI , GCHT  , GCQT  ,   & ! input
                 GCLT  , GCIT  , GCUT  , GCVT  , GDQI  ,   & ! input
                 gctrt,                                    &
                 KT    , ISTS  , IENS  , nctp  )             ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: im, IJSDIM, KMAX, NTR, nctp, ntrq !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(kind_phys)     GTT   (IJSDIM, KMAX)   !< temperature tendency
      REAL(kind_phys)     GTQ   (IJSDIM, KMAX, NTR)   !< moisture tendency
!     REAL(kind_phys)     GTCFRC(IJSDIM, KMAX)   !< cloud fraction tendency
      REAL(kind_phys)     GTU   (IJSDIM, KMAX)   !< u tendency
      REAL(kind_phys)     GTV   (IJSDIM, KMAX)   !< v tendency
!
!   [INPUT]
      REAL(kind_phys)     GDH   (IJSDIM, KMAX)      !< moist static energy
      REAL(kind_phys)     GDQ   (IJSDIM, KMAX, NTR) !< humidity qv
!     REAL(kind_phys)     GDCFRC(IJSDIM, KMAX)      !< cloud fraction
      REAL(kind_phys)     GDU   (IJSDIM, KMAX)
      REAL(kind_phys)     GDV   (IJSDIM, KMAX)
      REAL(kind_phys)     DELPI (IJSDIM, KMAX)
      REAL(kind_phys)     CBMFX (IM,     NCTP)      !< cloud base mass flux
      REAL(kind_phys)     GCYT  (IJSDIM, NCTP)      !< detraining mass flux
      REAL(kind_phys)     GCHT  (IJSDIM, NCTP)      !< detraining MSE
      REAL(kind_phys)     GCQT  (IJSDIM, NCTP)      !< detraining qv
      REAL(kind_phys)     GCLT  (IJSDIM, NCTP)      !< detraining ql
      REAL(kind_phys)     GCIT  (IJSDIM, NCTP)      !< detraining qi
      REAL(kind_phys)     GCtrT (IJSDIM, ntrq:ntr, NCTP)!< detraining tracer
      REAL(kind_phys)     GCUT  (IJSDIM, NCTP)      !< detraining u
      REAL(kind_phys)     GCVT  (IJSDIM, NCTP)      !< detraining v
      REAL(kind_phys)     GDQI  (IJSDIM, KMAX)      !< cloud ice
      INTEGER      KT    (IJSDIM, NCTP)      !< cloud top
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(kind_phys)     GTHCI, GTQVCI, GTXCI
      integer      I, K, CTP, kk,n 
!

      DO CTP=1,NCTP
        DO I=ISTS,IENS
          K = KT(I,CTP)
          IF (K > 0) THEN
            GTXCI  = DELPI(I,K)*CBMFX(I,CTP)

            GTHCI  = GTXCI * (GCHT(I,CTP) - GCYT(I,CTP)*GDH(I,K))
            GTQVCI = GTXCI * (GCQT(I,CTP) - GCYT(I,CTP)*GDQ(I,K,1))
!
            GTQ(I,K,1)   = GTQ(I,K,1)   + GTQVCI
            GTT(I,K)     = GTT(I,K)     + (GTHCI - EL*GTQVCI) * oneocp
! ql tendency by detrainment is treated by stratiform scheme
            GTQ(I,K,ITL) = GTQ(I,K,ITL) + GTXCI * (GCLT(I,CTP) - GCYT(I,CTP)*GDQ(I,K,ITL))
! qi tendency by detrainment is treated by stratiform scheme
            GTQ(I,K,ITI) = GTQ(I,K,ITI) + GTXCI * (GCIT(I,CTP) - GCYT(I,CTP)*GDQI(I,K))
            do n = ntrq,NTR
              GTQ(I,K,n) = GTQ(I,K,n)   + GTXCI * (GCtrT(I,n,CTP) - GCYT(I,CTP)*GDQ(I,K,n))
            enddo

!           GTCFRC(I,K)  = GTCFRC(I,K) + GTXCI * (GCYT(I,CTP) - GCYT(I,CTP)*GDCFRC(I,K))
            GTU(I,K)     = GTU(I,K)    + GTXCI * (GCUT(I,CTP) - GCYT(I,CTP)*GDU(I,K))
            GTV(I,K)     = GTV(I,K)    + GTXCI * (GCVT(I,CTP) - GCYT(I,CTP)*GDV(I,K))
          ENDIF
        ENDDO
      ENDDO
!
      END SUBROUTINE CUMDET
!***********************************************************************
!>\ingroup cs_scheme
      SUBROUTINE CUMSBH                             & !! adiabat. descent
               ( IM    , IJSDIM, KMAX  , NTR, ntrq, & !DD dimensions
                 GTT   , GTQ   ,                    & ! modified
                 GTU   , GTV   ,                    & ! modified
                 GDH   , GDQ   , GDQI  ,            & ! input
                 GDU   , GDV   ,                    & ! input
                 DELPI , GMFLX , GMFX0 ,            & ! input
                 KTMX  , CPRES , KB, ISTS  , IENS )   ! input
!
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, IM, KMAX, NTR, ntrq      !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(kind_phys)     GTT   (IJSDIM, KMAX)      !< Temperature tendency
      REAL(kind_phys)     GTQ   (IJSDIM, KMAX, NTR) !< Moisture etc tendency
      REAL(kind_phys)     GTU   (IJSDIM, KMAX)      !< u tendency
      REAL(kind_phys)     GTV   (IJSDIM, KMAX)      !< v tendency
!
!   [INPUT]
      REAL(kind_phys)     GDH   (IJSDIM, KMAX)
      REAL(kind_phys)     GDQ   (IJSDIM, KMAX, NTR) !< humidity etc
      REAL(kind_phys)     GDQI  (IJSDIM, KMAX)
      REAL(kind_phys)     GDU   (IJSDIM, KMAX)
      REAL(kind_phys)     GDV   (IJSDIM, KMAX)
      REAL(kind_phys)     DELPI (IJSDIM, KMAX)
      REAL(kind_phys)     GMFLX (IJSDIM, KMAX+1)    !< mass flux (updraft+downdraft)
      REAL(kind_phys)     GMFX0 (IJSDIM, KMAX+1)    !< mass flux (updraft only)
      INTEGER      KB(IJSDIM)                !< cloud base index - negative means no convection
      INTEGER      KTMX
      REAL(kind_phys)     CPRES                     !< pressure factor for cumulus friction
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(kind_phys)     SBH0, SBQ0, SBL0, SBI0, SBC0,  SBS0,        &
                   SBH1, SBQ1, SBL1, SBI1, SBC1,  SBS1,   FX1, &
                   SBU0, SBV0, SBU1, SBV1, GTHCI, GTQVCI,      &
                   GTQLCI, GTQICI, GTM2CI, GTM3CI, wrk, wrk1
      REAL(kind_phys)     FX(ISTS:IENS)

      REAL(kind_phys), dimension(IJSDIM, KMAX)  :: GTLSBH, GTISBH
      integer   :: I, K, KM, KP, n
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
          if (kb(i) > 0) then
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
            IF (GMFLX(I,K) > GMFLX(I,KP)) THEN
               FX1 = half
            ELSE
               FX1 = zero
            ENDIF
!
            wrk    = DELPI(I,K)
            wrk1   = one - FX(I)
            GTHCI  = wrk * (wrk1*SBH0 + FX1 *SBH1)
            GTQVCI = wrk * (wrk1*SBQ0 + FX1 *SBQ1)
            GTQLCI = wrk * (wrk1*SBL0 + FX1 *SBL1)
            GTQICI = wrk * (wrk1*SBI0 + FX1 *SBI1)
!
            GTT (I,K    ) = GTT(I,K)     + (GTHCI-EL*GTQVCI)*oneocp
            GTQ (I,K,1  ) = GTQ(I,K,1)   +  GTQVCI
            GTQ (I,K,ITL) = GTQ(I,K,ITL) +  GTQLCI
            GTQ (I,K,ITI) = GTQ(I,K,ITI) +  GTQICI
            GTU (I,K)     = GTU(I,K)     +  wrk * (wrk1*SBU0 + FX1*SBU1)
            GTV (I,K)     = GTV(I,K)     +  wrk * (wrk1*SBV0 + FX1*SBV1)
            DO n = ntrq, ntr
              GTQ (I,K,n) = GTQ(I,K,n) + wrk                                   &
                          * ( wrk1 * (GMFLX(I,KP) * (GDQ(I,KP,n)-GDQ(I,K ,n))) &
                             + FX1 * (GMFLX(I,K ) * (GDQ(I,K ,n)-GDQ(I,KM,n))) )
            ENDDO

            GTLSBH(I,K)   = GTQLCI
            GTISBH(I,K)   = GTQICI
!
!           SBC0 = GMFLX(I,K+1) * (GDQ(I,KP,IMU2)-GDQ(I,K,IMU2))
!           SBS0 = GMFLX(I,K+1) * (GDQ(I,KP,IMU3)-GDQ(I,K,IMU3))
!           SBC1 = GMFLX(I,K)   * (GDQ(I,K,IMU2)-GDQ(I,KM,IMU2))
!           SBS1 = GMFLX(I,K)   * (GDQ(I,K,IMU3)-GDQ(I,KM,IMU3))
!           GTM2CI = DELPI(I,K) * (( one-FX(I))*SBC0 + FX1 *SBC1)
!           GTM3CI = DELPI(I,K) * ((one-FX(I))*SBS0 + FX1 *SBS1)
!           GTQ(I,K,IMU2) = GTQ(I,K,IMU2) + GTM2CI
!           GTQ(I,K,IMU3) = GTQ(I,K,IMU3) + GTM3CI
!
            FX(I) = FX1
          endif
        enddo
      enddo
!
      END SUBROUTINE CUMSBH
!***********************************************************************
!
!***********************************************************************
!>\ingroup cs_scheme
!! This subroutine calculate cloud subsidence heating.
      SUBROUTINE CUMSBW                         & !! adiabat. descent
               ( IM    , IJSDIM, KMAX  ,        & !DD dimensions
                 GTU   , GTV   ,                & ! modified
                 GDU   , GDV   ,                & ! input
                 DELPI , GMFLX , GMFX0 ,        & ! input
                 KTMX  , CPRES , KB, ISTS  , IENS )   ! input
!
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSDIM, IM, KMAX!! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(kind_phys)     GTU   (IJSDIM, KMAX)      !< u tendency
      REAL(kind_phys)     GTV   (IJSDIM, KMAX)      !< v tendency
!
!   [INPUT]
      REAL(kind_phys)     GDU   (IJSDIM, KMAX)
      REAL(kind_phys)     GDV   (IJSDIM, KMAX)
      REAL(kind_phys)     DELPI (IJSDIM, KMAX)
      REAL(kind_phys)     GMFLX (IJSDIM, KMAX+1)    !< mass flux (updraft+downdraft)
      REAL(kind_phys)     GMFX0 (IJSDIM, KMAX+1)    !< mass flux (updraft only)
      INTEGER      KB(IJSDIM)                !< cloud base index - negative means no convection
      INTEGER      KTMX, ISTS, IENS
      REAL(kind_phys)     CPRES                     !< pressure factor for cumulus friction
!
!   [INTERNAL WORK]
      REAL(kind_phys)     FX1, SBU0, SBV0, SBU1, SBV1, wrk, wrk1
      REAL(kind_phys)     FX(ISTS:IENS)

      integer   :: I, K, KM, KP
!
!
      FX     = zero
!
      DO K=KTMX,1,-1
        KM = MAX(K-1, 1)
        KP = MIN(K+1, KMAX)
        DO I=ISTS,IENS
          if (kb(i) > 0) then
            SBU0 = GMFLX(I,KP) * (GDU(I,KP)-GDU(I,K))           &
                 - GMFX0(I,KP) * (GDU(I,KP)-GDU(I,K))*CPRES
            SBV0 = GMFLX(I,KP) * (GDV(I,KP)-GDV(I,K))           &
                 - GMFX0(I,KP) * (GDV(I,KP)-GDV(I,K))*CPRES
!
            SBU1 = GMFLX(I,K) * (GDU(I,K)-GDU(I,KM))            &
                 - GMFX0(I,K) * (GDU(I,K)-GDU(I,KM))*CPRES
            SBV1 = GMFLX(I,K) * (GDV(I,K)-GDV(I,KM))            &
                 - GMFX0(I,K) * (GDV(I,K)-GDV(I,KM))*CPRES
!
            IF (GMFLX(I,K) > GMFLX(I,KP)) THEN
               FX1 = half
            ELSE
               FX1 = zero
            ENDIF
!
            wrk      = DELPI(I,K)
            wrk1     = one - FX(I)
!
            GTU(I,K) = GTU(I,K) +  wrk * (wrk1*SBU0 + FX1*SBU1)
            GTV(I,K) = GTV(I,K) +  wrk * (wrk1*SBV0 + FX1*SBV1)
!
            FX(I) = FX1
          endif
        enddo
      enddo
!
      END SUBROUTINE CUMSBW
!***********************************************************************
!>\ingroup cs_scheme
!! This subroution calculates freeze, melt and evaporation in cumulus downdraft.
      SUBROUTINE CUMDWN                            & ! Freeze & Melt & Evaporation
               ( IM    , IJSDIM, KMAX  , NTR,ntrq,nctp, & !DD dimensions
                 GTT   , GTQ   , GTU   , GTV   ,        & ! modified
                         GMFLX ,                        & ! modified
                 GPRCP , GSNWP , GTEVP , GMDD  ,        & ! output
                 GPRCI , GSNWI ,                        & ! input
                 GDH   , GDW   , GDQ   , GDQI  ,        & ! input
                 GDQS  , GDS   , GDHS  , GDT   ,        & ! input
                 GDU   , GDV   , GDZ   ,                & ! input
                 GDZM  ,         FDQS  , DELP  ,        & ! input
                 DELPI ,                                &
                 sigmad, do_aw , do_awdd, flx_form,     & !DDsigma input
                 gtmelt, gtevap, gtsubl,                & !DDsigma input
                 dtdwn , dqvdwn, dqldwn, dqidwn,        & !DDsigma input
                 dtrdwn,                                &
                 KB    , KTMX  , ISTS  , IENS    )        ! input
!
! DD AW : modify to get eddy fluxes and microphysical tendencies for AW
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX, NTR , ntrq, nctp   !! DD, for GFS, pass in
      logical, intent(in) :: do_aw, do_awdd, flx_form
!
!   [MODIFY]
      REAL(kind_phys)     GTT   (IJSDIM, KMAX)       !< Temperature tendency
      REAL(kind_phys)     GTQ   (IJSDIM, KMAX, NTR)  !< Moisture etc tendency
      REAL(kind_phys)     GTU   (IJSDIM, KMAX)       !< u tendency
      REAL(kind_phys)     GTV   (IJSDIM, KMAX)       !< v tendency
      REAL(kind_phys)     GMFLX (IJSDIM, KMAX+1)     !< mass flux
!
!   [OUTPUT]
      REAL(kind_phys)     GPRCP (IJSDIM, KMAX+1)     !< rainfall flux
      REAL(kind_phys)     GSNWP (IJSDIM, KMAX+1)     !< snowfall flux
      REAL(kind_phys)     GTEVP (IJSDIM, KMAX)       !< evaporation+sublimation
      REAL(kind_phys)     GMDD  (IJSDIM, KMAX+1)     !< downdraft mass flux

!AW microphysical tendencies
      REAL(kind_phys)     gtmelt (IJSDIM, KMAX)      !< t tendency ice-liq
      REAL(kind_phys)     gtevap (IJSDIM, KMAX)      !< t tendency liq-vapor
      REAL(kind_phys)     gtsubl (IJSDIM, KMAX)      !< t tendency ice-vapor
!AW eddy flux tendencies
      REAL(kind_phys)     dtdwn  (IJSDIM, KMAX)      !< t tendency downdraft detrainment
      REAL(kind_phys)     dqvdwn (IJSDIM, KMAX)      !< qv tendency downdraft detrainment
      REAL(kind_phys)     dqldwn (IJSDIM, KMAX)      !< ql tendency downdraft detrainment
      REAL(kind_phys)     dqidwn (IJSDIM, KMAX)      !< qi tendency downdraft detrainment
      REAL(kind_phys)     dtrdwn (IJSDIM, KMAX, ntrq:ntr) !< tracer tendency downdraft detrainment

!   [INPUT]
      REAL(kind_phys)     GPRCI (IJSDIM, KMAX)       !< rainfall generation
      REAL(kind_phys)     GSNWI (IJSDIM, KMAX)       !< snowfall generation
      REAL(kind_phys)     GDH   (IJSDIM, KMAX)       !< moist static energy
      REAL(kind_phys)     GDW   (IJSDIM, KMAX)       !< total water
      REAL(kind_phys)     GDQ   (IJSDIM, KMAX, NTR)  !< humidity etc
      REAL(kind_phys)     GDQI  (IJSDIM, KMAX)       !< cloud ice
      REAL(kind_phys)     GDQS  (IJSDIM, KMAX)       !< saturate humidity
      REAL(kind_phys)     GDS   (IJSDIM, KMAX)       !< dry static energy
      REAL(kind_phys)     GDHS  (IJSDIM, KMAX)       !< saturate moist static energy
      REAL(kind_phys)     GDT   (IJSDIM, KMAX)       !< air temperature T
      REAL(kind_phys)     GDU   (IJSDIM, KMAX)       !< u-velocity
      REAL(kind_phys)     GDV   (IJSDIM, KMAX)       !< v-velocity
      REAL(kind_phys)     GDZ   (IJSDIM, KMAX)       !< altitude
      REAL(kind_phys)     GDZM  (IJSDIM, KMAX+1)     !< altitude (half lev)
      REAL(kind_phys)     FDQS  (IJSDIM, KMAX)
      REAL(kind_phys)     DELP  (IJSDIM, KMAX)
      REAL(kind_phys)     DELPI (IJSDIM, KMAX)
      INTEGER      KB    (IJSDIM)
      INTEGER      KTMX, ISTS, IENS
      REAL(kind_phys)     sigmad (IM,KMAX+1)         !< DDsigma cloud downdraft area fraction

!
!   [INTERNAL WORK]
! Note: Some variables have 3-dimensions for the purpose of budget check.
      REAL(kind_phys)     EVAPD (IJSDIM, KMAX)        !< evap. in downdraft
      REAL(kind_phys)     SUBLD (IJSDIM, KMAX)        !< subl. in downdraft
      REAL(kind_phys)     EVAPE (IJSDIM, KMAX)        !< evap. in environment
      REAL(kind_phys)     SUBLE (IJSDIM, KMAX)        !< subl. in environment
      REAL(kind_phys)     EVAPX (IJSDIM, KMAX)        !< evap. env. to DD
      REAL(kind_phys)     SUBLX (IJSDIM, KMAX)        !< subl. env. to DD
      REAL(kind_phys)     GMDDE (IJSDIM, KMAX)        !< downdraft entrainment
      REAL(kind_phys)     SNMLT (IJSDIM, KMAX)        !< melt - freeze
      REAL(kind_phys)     GCHDD (IJSDIM, KMAX)        !< MSE detrainment
      REAL(kind_phys)     GCWDD (IJSDIM, KMAX)        !< water detrainment
      REAL(kind_phys)     GTTEV (IJSDIM, KMAX)        !< T tendency by evaporation
      REAL(kind_phys)     GTQEV (IJSDIM, KMAX)        !< q tendency by evaporation
      REAL(kind_phys)     GCHD  (ISTS:IENS)           !< downdraft MSE
      REAL(kind_phys)     GCWD  (ISTS:IENS)           !< downdraft q
! profiles of downdraft variables for AW flux tendencies
      REAL(kind_phys)     GCdseD(ISTS:IENS, KMAX)     !< downdraft dse
      REAL(kind_phys)     GCqvD (ISTS:IENS, KMAX)     !< downdraft qv
      REAL(kind_phys)     GCqlD (ISTS:IENS, KMAX)     !< downdraft ql
      REAL(kind_phys)     GCqiD (ISTS:IENS, KMAX)     !< downdraft qi
      REAL(kind_phys)     GCtrD (ISTS:IENS, ntrq:ntr) !< downdraft tracer

      REAL(kind_phys)     GCUD  (ISTS:IENS)           !< downdraft u
      REAL(kind_phys)     GCVD  (ISTS:IENS)           !< downdraft v
      REAL(kind_phys)     FSNOW (ISTS:IENS)
      REAL(kind_phys)     GMDDD (ISTS:IENS)
      INTEGER      I, K
      REAL(kind_phys)     GDTW
      REAL(kind_phys)     GCHX, GCTX, GCQSX, GTPRP, EVSU, GTEVE, LVIC
      REAL(kind_phys)     DQW, DTW, GDQW, DZ, GCSD, FDET, GDHI
      REAL(kind_phys)     GMDDX, GMDDMX
      REAL(kind_phys)     GCHDX, GCWDX
      REAL(kind_phys)     GCUDD, GCVDD
      REAL(kind_phys)     GTHCI, GTQVCI, GTQLCI, GTQICI
!M    REAL(kind_phys)     GTHCI, GTQVCI, GTQLCI, GTQICI, GTUCI, GTVCI
      real(kind_phys)     wrk, fmelt, fevp,  gctrdd(ntrq:ntr)
!DD#ifdef OPT_CUMBGT
      REAL(kind_phys)     WBGT  ( ISTS:IENS )         !! water budget
      REAL(kind_phys)     HBGT  ( ISTS:IENS )         !! energy budget
      REAL(kind_phys)     DDWBGT( ISTS:IENS )         !! downdraft water budget
      REAL(kind_phys)     DDHBGT( ISTS:IENS )         !! downdraft energy budget
      REAL(kind_phys)     WMX, HMX, DDWMX, DDHMX, tx1, wrk1, wrk2, wrk3, wrk4, wrkn
      REAL(kind_phys) dp_above, dp_below
      real(kind_phys) fsigma
      integer ij, n, kp1
!DD#endif
!
!   [INTERNAL PARM]
      REAL(kind_phys), parameter :: TWSNOW = 273.15_kind_phys   !< wet-bulb temp. rain/snow
      REAL(kind_phys), parameter :: FTMLT  = 4._kind_phys       !< temp. factor for melt
      REAL(kind_phys), parameter :: GMFLXC = 5.e-2_kind_phys    !< critical mass flux
      REAL(kind_phys), parameter :: VTERMS = 2._kind_phys       !< terminal velocity of snowflake
!     REAL(kind_phys), parameter :: MELTAU = 10._kind_phys      !< melting timescale
      REAL(kind_phys), parameter :: MELTAU = 20._kind_phys      !< melting timescale       ! Moorthi june 30, 2017
!
!     REAL(kind_phys), parameter :: EVAPR  = 0.4_kind_phys      !< evaporation factor      ! Moorthi June 28, 2017
      REAL(kind_phys), parameter :: EVAPR  = 0.3_kind_phys      !< evaporation factor
!     REAL(kind_phys), parameter :: EVAPR  = 0._kind_phys       !< evaporation factor
      REAL(kind_phys), parameter :: REVPDD = 1._kind_phys       !< max rate of DD to evapolation
      REAL(kind_phys), parameter :: RDDR   = 5.e-4_kind_phys    !< DD rate (T0 R0 W0)^-1
!     REAL(kind_phys), parameter :: RDDR   = 0._kind_phys       !< DD rate (T0 R0 W0)^-1
      REAL(kind_phys), parameter :: RDDMX  = 0.5_kind_phys      !< norm. flux of downdraft
      REAL(kind_phys), parameter :: VTERM  = 5._kind_phys       !< term. vel. of precip.
!     REAL(kind_phys), parameter :: VTERM  = 4._kind_phys       !< term. vel. of precip.   ! Moorthi June 28, 2017
      REAL(kind_phys), parameter :: EVATAU = 2._kind_phys       !< evaporation/sublimation timescale
      REAL(kind_phys), parameter :: ZDMIN  = 5.e2_kind_phys     !< min altitude of downdraft detrainment
      real(kind_phys), parameter :: evapovtrm=EVAPR/VTERM

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
!         GCqlD (I,k) = zero 
!         GCqiD (I,k) = zero 
          gtevap(I,k) = zero 
          gtmelt(I,k) = zero 
          gtsubl(I,k) = zero 
        enddo
      enddo

!  These are zeroed by the calling routine, cs_cumlus
!          do k=1,kmax
!            do i=ists,iens
!              dtdwn (I,k) = zero
!              dqvdwn(I,k) = zero
!              dqldwn(I,k) = zero
!              dqidwn(I,k) = zero
!            enddo
!          enddo
!          do n=ntrq,ntr
!            do k=1,kmax
!              do i=ists,iens
!                dtrdwn(i,k,n) = zero
!              enddo
!            enddo
!          enddo
!
      do i=ists,iens
        GCHD(I) = zero
        GCWD(I) = zero
        GCUD(I) = zero
        GCVD(I) = zero
      enddo
      do n=ntrq,ntr
        do i=ists,iens
          GCtrD (I,n) = zero 
        enddo
      enddo
!
      DO K=KTMX,1,-1   ! loop A
        kp1 = min(k+1,kmax)
!
!     < precipitation melt & freeze >
!
        DO I=ISTS,IENS
          if (kb(i) > 0) then
            GTPRP = GPRCP(I,KP1) + GSNWP(I,KP1)
            IF (GTPRP > zero) THEN
              FSNOW(I) = GSNWP(I,KP1) / GTPRP
            ELSE
              FSNOW(I) = zero
            ENDIF
            LVIC = ELocp + EMELTocp*FSNOW(I)
            GDTW = GDT(I,K) - LVIC*(GDQS(I,K) - GDQ(I,K,1)) &
                            / (one + LVIC*FDQS(I,K))

            DZ = GDZM(I,KP1) - GDZM(I,K)
            FMELT      = (one + FTMLT*(GDTW - TWSNOW))     &
                      * (one - TANH(GMFLX(I,KP1)/GMFLXC)) &
                       * (one - TANH(VTERMS*MELTAU/DZ))
            IF (GDTW  < TWSNOW) THEN
              SNMLT(I,K) = GPRCP(I,KP1)*min(max(FMELT, one), zero)
            ELSE
              SNMLT(I,K) = GSNWP(I,KP1)*max(min(FMELT, one), zero)
            ENDIF
            GSNWP(I,K) = GSNWP(I,KP1)+GSNWI(I,K) - SNMLT(I,K)
            GPRCP(I,K) = GPRCP(I,KP1)+GPRCI(I,K) + SNMLT(I,K)
            GTTEV(I,K) = -EMELToCP * SNMLT(I,K) * DELPI(I,K)
!DD heating rate due to precip melting for AW
            gtmelt(i,k) = gtmelt(i,k) + GTTEV(I,K)
          endif
        ENDDO
!
!     < downdraft >
!
        DO I=ISTS,IENS   ! loop B
          if (kb(i) > 0) then
            wrk  = delpi(i,k)
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
                  GCHD(I) = GCHDX
                  GCWD(I) = GCWDX
                  GCUD(I) = GCUD(I) + GDU(I,K)*GMDDE(I,K)
                  GCVD(I) = GCVD(I) + GDV(I,K)*GMDDE(I,K)
                  do n = ntrq,ntr
                    GCtrD(I,n) = GCtrD(I,n) + GDq(I,K,n)*GMDDE(I,K)
                  enddo
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
              do n = ntrq,ntr
                GCtrDD(n)    = FDET*GCtrD(I,n)
              enddo
!
              GTHCI  =  wrk * (GCHDD(I,K) - GMDDD(I)*GDH(I,K))
              GTQVCI =  wrk * (GCWDD(I,K) - GMDDD(I)*GDQ(I,K,1))
!
              GTT (I,K)     = GTT(I,K)     + (GTHCI - EL*GTQVCI)*oneoCP
              GTQ (I,K,1)   = GTQ(I,K,1)   + GTQVCI
              GTQ (I,K,ITL) = GTQ(I,K,ITL) - wrk * GMDDD(I)*GDQ(I,K,ITL)
              GTQ (I,K,ITI) = GTQ(I,K,ITI) - wrk * GMDDD(I)*GDQI(I,K)

              do n = ntrq,ntr
                GTQ (I,K,n)  = GTQ(I,K,n)   + wrk * (GCtrDD(n) - GMDDD(I)*GDQ(I,K,n))
                GCtrD(I,n)   = GCtrD(I,n)   - GCtrDD(n)
              enddo

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
          endif
        ENDDO   ! loop B
!
      ENDDO   ! loop A
!
      DO K=1,KTMX
        kp1 = min(k+1,kmax)
        DO I=ISTS,IENS
          if (kb(i) > 0) then
            wrk = DELPI(I,k)
            tx1 = DELPI(I,kp1)
              
            GTTEV(I,K) = GTTEV(I,K) - wrk                              &
                       * (ELocp*EVAPE(I,K)+(ELocp+EMELTocp)*SUBLE(I,K))
            GTT(I,K)   = GTT(I,K) + GTTEV(I,K)
!
            GTQEV(I,K) = GTQEV(I,K) + (EVAPE(I,K)+SUBLE(I,K)) * wrk
            GTQ(I,K,1) = GTQ(I,K,1) + GTQEV(I,K)
!
            GMFLX(I,K) = GMFLX(I,K) - GMDD(I,K)
          endif
        ENDDO   ! end of i loop
      ENDDO     ! end of k loop

! AW tendencies due to vertical divergence of eddy fluxes
      DO K=2,KTMX
         kp1 = min(k+1,kmax)
        DO I=ISTS,IENS
          if (kb(i) > 0) then
            if (k > 1 .and. flx_form) then
              fsigma        = one - sigmad(i,kp1)
              dp_below      = wrk * (one - sigmad(i,k))
              dp_above      = tx1 * (one - sigmad(i,kp1))

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
              do n = ntrq,ntr
                wrkn            = gmdd(i,kp1) * gdq(i,k,n)
                dtrdwn(i,k,n)   = dtrdwn(i,k,n)   + dp_below * wrkn
                dtrdwn(i,kp1,n) = dtrdwn(i,kp1,n) - dp_above * wrkn
              enddo
            endif

          endif
        ENDDO   ! end of i loop
      ENDDO     ! end of k loop
!
      END SUBROUTINE CUMDWN
!***********************************************************************
!>\ingroup cs_scheme
!> This subroutine computes cumulus cloudiness.
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
      REAL(kind_phys)     CUMFRC(IJSDIM)          !< cumulus cloud fraction
!
!   [MODIFY]
      REAL(kind_phys)     CUMCLW(IJSDIM, KMAX)    !< cloud water in cumulus
      REAL(kind_phys)     QLIQ  (IJSDIM, KMAX)    !< cloud liquid
      REAL(kind_phys)     QICE  (IJSDIM, KMAX)    !< cloud ice
      REAL(kind_phys)     FLIQC (IJSDIM, KMAX)    !< liquid ratio in cumulus
!
!   [INPUT]
      REAL(kind_phys)     GMFLX (IJSDIM, KMAX+1)  ! cumulus mass flux
      INTEGER      KTMX
      INTEGER      ISTS, IENS
!
!   [WORK]
      INTEGER      I, K
      REAL(kind_phys)     CUMF, QC, wrk
      LOGICAL, SAVE :: OFIRST = .TRUE.
!
!   [INTERNAL PARAM]
      REAL(kind_phys) :: FACLW  = 0.1_kind_phys     !> Mc->CLW
      REAL(kind_phys) :: CMFMIN = 2.e-3_kind_phys   !> Mc->cloudiness
      REAL(kind_phys) :: CMFMAX = 3.e-1_kind_phys   !> Mc->cloudiness
      REAL(kind_phys) :: CLMIN  = 1.e-3_kind_phys   !> cloudiness Min.
      REAL(kind_phys) :: CLMAX  = 0.1_kind_phys     !> cloudiness Max.
      REAL(kind_phys), SAVE :: FACLF
!     
      IF ( OFIRST ) THEN
         FACLF = (CLMAX-CLMIN)/LOG(CMFMAX/CMFMIN)
         OFIRST = .FALSE.
      END IF
                       
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
!>\ingroup cs_scheme
!! This subroutine calculates
      SUBROUTINE CUMUPR                                    & !! Tracer Updraft
               ( im    , IJSDIM, KMAX  , NTR   ,           & !DD dimensions
                 GTR   , GPRCC ,                           & ! modified
                 GDR   , CBMFX ,                           & ! input
                 GCYM  , GCYT  , GCQT  , GCLT  , GCIT  ,   & ! input
                 GTPRT , GTEVP , GTPRC0,                   & ! input
                 KB    , KBMX  , KT    , KTMX  , KTMXT ,   & ! input
                 DELPI , OTSPT , ISTS  , IENS,             & ! input
                 fscav, fswtr, nctp)
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: im, IJSDIM, KMAX, NTR, nctp             !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(kind_phys)     GTR   (IJSDIM, KMAX, NTR)
      REAL(kind_phys)     GPRCC (IJSDIM, NTR)
!
!   [INPUT]
      REAL(kind_phys)     GDR   (IJSDIM, KMAX, NTR)
      REAL(kind_phys)     CBMFX (IM, NCTP)
      REAL(kind_phys)     GCYM  (IJSDIM, KMAX, nctp)
      REAL(kind_phys)     GCYT  (IJSDIM, NCTP)
      REAL(kind_phys)     GCQT  (IJSDIM, NCTP)
      REAL(kind_phys)     GCLT  (IJSDIM, NCTP)
      REAL(kind_phys)     GCIT  (IJSDIM, NCTP)
      REAL(kind_phys)     GTPRT (IJSDIM, NCTP)
      REAL(kind_phys)     GTEVP (IJSDIM, KMAX)
      REAL(kind_phys)     GTPRC0(IJSDIM)   !! precip. before evap.
      real(kind_phys)     fscav(ntr), fswtr(ntr)
      INTEGER      KB    (IJSDIM )
      INTEGER      KBMX
      INTEGER      KT    (IJSDIM, NCTP)
      INTEGER      KTMX  (NCTP)
      INTEGER      KTMXT
      REAL(kind_phys)     DELPI (IJSDIM, KMAX)
      LOGICAL      OTSPT (NTR)              !! transport with this routine?
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER      I, K, LT, TP, CTP
      REAL(kind_phys)                            :: GCRTD, SCAV, GCWT, GPRCR, evpf, cbmfxl
      REAL(kind_phys), dimension(ists:iens)      :: GCRB, GCRT,  DR,   gtprc0i
!     REAL(kind_phys), dimension(ists:iens,kmax) :: DGCB, DZ,    RDZM,  EVPF
!     REAL(kind_phys), dimension(ists:iens,nctp) :: DZT,  RGCWT, MASK1, MASK2
      REAL(kind_phys), dimension(ists:iens,nctp) ::       RGCWT, MASK1
!
      do i=ists,iens
        if (gtprc0(i) > zero) then
          gtprc0i(i) = one / gtprc0(i)
        else
          gtprc0i(i) = zero
        endif
      enddo
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
          IF (kb(i) > 0 .and. K > KB(I)) THEN
            MASK1(I,CTP) = one
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
          DO CTP=1,NCTP
            DO I=ISTS,IENS
              GCRB(i)  = zero
              DR(i)    = zero
            enddo
            DO K=1,KBMX
              DO I=ISTS,IENS
                IF (kb(i) > 0 .and. K < KB(I)) THEN
                  GCRB(I) = GCRB(I) + (GCYM(I,K+1,ctp)-GCYM(I,K,ctp))* GDR(I,K,LT)
                ENDIF
              ENDDO
            ENDDO
!
            DO K=2,KTMX(CTP)
              DO I=ISTS,IENS
                IF (kb(i) > 0 .and. K >= KB(I) .AND.  K < KT(I,CTP)) THEN
                  DR(I) = DR(I) + (GCYM(I,K+1,ctp)-GCYM(I,K,ctp)) * GDR(I,K,LT)
                ENDIF
              ENDDO
            ENDDO
!
            DO I=ISTS,IENS
              K = KT(I,CTP)
              if (kb(i) > 0 .and. k > kb(i)) then
                DR(I)   = DR(I) + (GCYT(I,CTP) - GCYM(I,K,ctp)) * GDR (I,K,LT) &
                                             * MASK1(I,CTP)
                GCRT(I) = (GCRB(I) + DR(I))  * MASK1(I,CTP)
!
                SCAV    = FSCAV(LT)*GTPRT(I,CTP) + FSWTR(LT)*GTPRT(I,CTP)*RGCWT(I,CTP)
                SCAV    = MIN(SCAV, one)
                GCRTD   = GCRT(I) * (one - SCAV)
                cbmfxl  = max(zero, CBMFX(I,CTP))
                GPRCR   = SCAV * GCRT(I) * CBMFXl

                GTR(I,K,LT) = GTR(I,K,LT) + DELPI(I,K) * CBMFXl          &
                                * (GCRTD - GCYT(I,CTP) * GDR(I,K,LT))
                GPRCC(I,LT) = GPRCC(I,LT) + GPRCR

!               GPRCR   = SCAV * GCRT(I) * CBMFX(I,CTP)
!               GTR(I,K,LT) = GTR(I,K,LT) + DELPI(I,K) * CBMFX(I,CTP) &
!                             * (GCRTD - GCYT(I,CTP) * GDR(I,K,LT)) * MASK2(I,CTP)
!               GPRCC(I,LT) = GPRCC(I,LT) + GPRCR * MASK2(I,CTP)
              endif
            ENDDO
          ENDDO
!
          DO K=KTMXT,1,-1
            DO I=ISTS,IENS
              evpf = GTEVP(i,k) * gtprc0i(i)
              GTR(I,K,LT) = GTR(I,K,LT) + DELPI(I,K) * GPRCC(I,LT) * EVPF
              GPRCC(I,LT) = GPRCC(I,LT) * (one - EVPF)
!             GTR(I,K,LT) = GTR(I,K,LT) + DELPI(I,K) * GPRCC(I,LT) * EVPF(I,K)
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
!>\ingroup cs_scheme
      SUBROUTINE CUMDNR                                 & !! Tracer Downdraft
                      ( IM    , IJSDIM, KMAX  , NTR   , & !DD dimensions
                        GTR   ,                         & ! modified
                        GDR   , GMDD  , DELPI ,         & ! input
                        KTMX  , OTSPT , ISTS  , IENS )    ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX, NTR             !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(kind_phys)     GTR   (IJSDIM, KMAX, NTR)   ! Temperature tendency
!
!   [INPUT]
      REAL(kind_phys)     GDR   (IJSDIM, KMAX, NTR)
      REAL(kind_phys)     GMDD  (IJSDIM, KMAX)        ! downdraft mass flux
      REAL(kind_phys)     DELPI (IJSDIM, KMAX  )
      LOGICAL      OTSPT (NTR)
      INTEGER      KTMX, ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(kind_phys)     GCRD  (ISTS:IENS)           ! downdraft q
      REAL(kind_phys)     GMDDE, GMDDD, GCRDD
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
                GTR(I,K,LT) = GTR(I,K,LT) + DELPI(I,K) &
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
!>\ingroup cs_scheme
      SUBROUTINE CUMSBR                                      & !! Tracer Subsidence
                      ( IM    , IJSDIM, KMAX, NTR, NCTP,     & !DD dimensions
                        GTR   ,                              & ! modified
                        GDR   , DELP  ,                      & ! input
                        GMFLX , KTMX  , OTSPT ,              & ! input
                        sigmai        , sigma ,              & !DDsigma input
                        ISTS, IENS )                           ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX, NTR, nctp       !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(kind_phys)     GTR   (IJSDIM, KMAX, NTR)   !! tracer tendency
!
!   [INPUT]
      REAL(kind_phys)     GDR   (IJSDIM, KMAX, NTR)   !! tracer
      REAL(kind_phys)     DELP  (IJSDIM, KMAX)
      REAL(kind_phys)     GMFLX (IJSDIM, KMAX+1)      !! mass flux
      INTEGER      KTMX
      LOGICAL      OTSPT (NTR)                 !! tracer transport on/off
      INTEGER      ISTS, IENS
      REAL(kind_phys)     sigmai (IM,KMAX+1,NCTP), sigma(IM,KMAX+1)   !!DDsigma cloud updraft fraction
!
!   [INTERNAL WORK]
      INTEGER      I, K, KM, KP, LT
      REAL(kind_phys)     SBR0, SBR1, FX1
      REAL(kind_phys)     FX(ISTS:IENS)
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
              SBR0 = GMFLX(I,K+1) * (GDR(I,KP,LT) - GDR(I,K,LT))
              SBR1 = GMFLX(I,K)   * (GDR(I,K,LT)  - GDR(I,KM,LT))
              IF (GMFLX(I,K) > GMFLX(I,K+1)) THEN
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
!>\ingroup cs_scheme
!! This subroutine calculates tracer mass fixer without detrainment.
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
      REAL(kind_phys)     GTR   (IJSDIM, KMAX, NTR)   ! tracer tendency
!
!   [INPUT]
      REAL(kind_phys)     GDR   (IJSDIM, KMAX, NTR)   ! tracer
      REAL(kind_phys)     DELP  (IJSDIM, KMAX)
      REAL(kind_phys)     DELTA                       ! time step
      INTEGER      KTMX
      INTEGER      IMFXR (NTR)
        ! 0: mass fixer is not applied
        !    tracers which may become negative values
        !    e.g. subgrid-PDFs
        ! 1: mass fixer is applied, total mass may change through cumulus scheme
        !    e.g. moisture, liquid cloud, ice cloud, aerosols
        ! 2: mass fixer is applied, total mass never change through cumulus scheme
        !    e.g. CO2
        !DD add new CASE
        ! 3: just fill holes, no attempt to conserve
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(kind_phys)     GDR1
      REAL(kind_phys)     GDR2  (ISTS:IENS, KMAX)
      REAL(kind_phys), dimension(ISTS:IENS) :: TOT0, TOT1, TRAT
      REAL(kind_phys)     FWAT
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
          CASE (3)
          CASE DEFAULT
            EXIT
        END SELECT
!
        DO I=ISTS,IENS
          TOT0(I) = zero
          TOT1(I) = zero
          TRAT(I) = one
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
        if(imfxr(LT) .ne. 3) then
        DO I=ISTS,IENS
          IF (TOT1(I) > zero ) THEN
            TRAT(I) = MAX(TOT0(I), zero) / TOT1(I)
          ENDIF
        ENDDO
        endif
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
!>\ingroup cs_scheme
      SUBROUTINE CUMFXR1                                   & ! Tracer mass fixer
               ( IM    , IJSDIM, KMAX  ,nctp,              & !DD dimensions
                 GTR   ,                                   & ! modified
                 GDR   , DELP  , DELTA , KTMX  , IMFXR ,   & ! input
                 ISTS  , IENS                            )   ! input
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IM, IJSDIM, KMAX, nctp           !! DD, for GFS, pass in
!
!   [MODIFY]
      REAL(kind_phys)     GTR   (IJSDIM, KMAX)      ! tracer tendency
!
!   [INPUT]
      REAL(kind_phys)     GDR   (IJSDIM, KMAX)      ! tracer
      REAL(kind_phys)     DELP  (IJSDIM, KMAX)
      REAL(kind_phys)     DELTA                     ! time step
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
      REAL(kind_phys)     GDR1
      REAL(kind_phys)     GDR2  (ISTS:IENS, KMAX)
      REAL(kind_phys), dimension(ISTS:IENS) :: TOT0, TOT1, TRAT
      REAL(kind_phys)     FWAT
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
!>\ingroup cs_scheme
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
      REAL(kind_phys)     GTT   (IJSDIM, KMAX)      ! heating rate
      REAL(kind_phys)     GTQ   (IJSDIM, KMAX, NTR) ! change in q
      REAL(kind_phys)     GTU   (IJSDIM, KMAX)      ! tendency of u
      REAL(kind_phys)     GTV   (IJSDIM, KMAX)      ! tendency of v
      REAL(kind_phys)     GPRCC (IJSDIM, NTR )      ! rainfall
      REAL(kind_phys)     GSNWC (IJSDIM)            ! snowfall
      REAL(kind_phys)     CUMCLW(IJSDIM, KMAX)      ! cloud water in cumulus
      REAL(kind_phys)     CUMFRC(IJSDIM)            ! cumulus cloud fraction
      REAL(kind_phys)     GTCFRC(IJSDIM, KMAX)      ! change in cloud fraction
      REAL(kind_phys)     FLIQC (IJSDIM, KMAX)      ! liquid ratio in cumulus
      REAL(kind_phys)     GTPRP (IJSDIM, KMAX)      ! rain+snow flux
!
      INTEGER    ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER    I, K
!
!   [INTERNAL PARM]
      REAL(kind_phys) :: GTTMAX  = 1.e-2_kind_phys
      REAL(kind_phys) :: GTQVMAX = 1.e-4_kind_phys
      REAL(kind_phys) :: GTQLMAX = 1.e-5_kind_phys
      REAL(kind_phys) :: GTUMAX  = 1.e-2_kind_phys
      REAL(kind_phys) :: GTVMAX  = 1.e-2_kind_phys
      REAL(kind_phys) :: GTCFMAX = 1.e-3_kind_phys
      REAL(kind_phys) :: PRCCMAX = 1.e-2_kind_phys
      REAL(kind_phys) :: SNWCMAX = 1.e-2_kind_phys
      REAL(kind_phys) :: CLWMAX  = 1.e-3_kind_phys
      REAL(kind_phys) :: TPRPMAX = 1.e-2_kind_phys
      REAL(kind_phys) :: GTQIMAX = 1.e-5_kind_phys
      !REAL(kind_phys) :: GTM2MAX = 1._kind_phys
      !REAL(kind_phys) :: GTM3MAX = 1._kind_phys
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

end module cs_conv

   module canopy_levs_mod
   contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   subroutine canopy_levs_init(im, ix, km,         &
              ntrac1,  ntqv, ntke,                 &
              zi, zl, zm,                          & ! in: 3D meters
              prsl, prsi,                          & ! in: 3D (Pa)
              dv, du, tdt, rtg,                    & ! in: 3D
              U1, V1, T1, Q1,                      & ! in: 3D " 4D q1(ix,km,ntrac1) kg kg-1
              dens, dkt, dku,                      & ! in: 3D
              dtend,                               & ! in: 4D
              zmom_can3, zmid_can3,                & !out: 3D
              sigmom_can, sigmid_can,              & !out
              ZH_CAN, ZF_CAN,                      & !out
              PRSL_CAN, PRSI_CAN,                  & !out
              dv_can, du_can, tdt_can, rtg_can,    & !out: 3D
              T1_CAN, QV_CAN, DENS_CAN,            & !out
              WS_CAN, DKT_CAN, DKU_CAN,            & !out
              Q1_CAN, Q1_2M,                       & !out
              DTEND_CAN )

   use machine , only : kind_phys
! Allocated in mfpbltq_mod:  q1(ix,km,ntrac1)  t1(ix,km) u1(ix,km), v1(ix,km)
   use mfpbltq_mod
   use canopy_mask_mod

   IMPLICIT NONE

!...Arguments:
! ntrac1 = ntrac - 1
   integer, intent(in)  :: im, ix, km, ntrac1, ntqv, ntke

   real(kind=kind_phys), intent(in) ::  zi(:,:),   zl(:,:),   zm(:,:), &
                                      prsi(:,:), prsl(:,:)
   real(kind=kind_phys), intent(in) :: dv(:,:),    du(:,:),         &
                                      tdt(:,:),   rtg(:,:,:)
   real(kind=kind_phys), intent(in) ::  u1(:,:),   v1(:,:), t1(:,:)
   real(kind=kind_phys), intent(in) :: dens(:,:), dkt(:,:),  dku(:,:)
   real(kind=kind_phys), intent(in) :: dtend(:,:,:)

! ** Q1 is concentration field (including gas and aerosol variables) mass mixing ratio kg kg-1
!   NB. mfpbltq_mod: q1(ix,km,ntrac1)
   real(kind=kind_phys), intent(in) :: Q1(:,:,:)  ! consider only gas-phase species (NO aerosol species)

   real(kind=kind_phys), intent(out) ::  &
! tendencies
             DTEND_CAN (:, :, :)     , & ! dim(km , ndtend)
             dv_can    (:, :)        , & ! dim(km)
             du_can    (:, :)        , & ! dim(km)
             tdt_can   (:, :)        , & ! dim(km)
! tendencies all gas-phase species & TKE
             RTG_CAN   (:, :, :)     , & ! dim(km )

! met3d arrays
             ZH_CAN    (:, :)        , & ! dim(nkt)
             ZF_CAN    (:, :)        , & ! dim(nkt)
             T1_CAN    (:, :)        , & ! dim(nkt)
             QV_CAN    (:, :)        , & ! dim(nkt)
             WS_CAN    (:, :)        , & ! dim(nkt)
             PRSL_CAN  (:, :)        , & ! dim(nkt)
             PRSI_CAN  (:, :)        , & ! dim(nkt+1)
             DENS_CAN  (:, :)        , & ! dim(nkt)
             DKT_CAN   (:, :)        , & ! dim(nkt)
             DKU_CAN   (:, :)        , & ! dim(nkt)
! all gas-phase species array
             Q1_CAN    (:, :, :)     , & ! dim(nkt)
             Q1_2M     (:, :)        , & ! dim(nkt)
! canopy layers height arrays
             zmom_can3 (:, :)        , & ! dim(nkt+1)    ! Paul's sigmcan(:,nkt)
             zmid_can3 (:, :)        , & ! dim(nkt)      ! Paul's sigtcan(:,nkt)
             sigmom_can(:, :)        , & ! dim(nkt) ~ prsi(:,km+1)
             sigmid_can(:, :)            ! dim(nkt) ~ prsl(:,km)

!...local variables

   character(256) :: errmsg
   integer        :: errflg

   integer        :: k, kc

! Initialize with values before in-canopy diffusion

! Layers height
   zmom_can3(:,:) = 0.
   zmid_can3(:,:) = 0.
   sigmom_can(:,:) = 0.
   sigmid_can(:,:) = 0.

! met3d arrays
   ZH_CAN   (:,:) = 0.
   ZF_CAN   (:,:) = 0.

! Zero in-canopy tendencies
   dtend_can(:, :, : ) = 0.0

! Tracers
   Q1_2M (:,           :) = Q1(:,1,    :)       ! kg kg-1

! Subset (km combined layers minus top nkc layers)
   do k = 1, km-nkc
     ! km     is top combined subset
     ! nkc+1  is bot combined
     kc= nkc+k     ! 4th from top (nkt) to nkc+1 combined canopy plus resolved model layer

! Tendencies
     DU_CAN  (:,kc)    = DU  (:,k)       ! m s-2
     DV_CAN  (:,kc)    = DV  (:,k)       ! m s-2
     TDT_CAN (:,kc)    = TDT (:,k)       ! K s-1

     RTG_CAN (:,kc, 1:ntrac1) = RTG (:,k, 1:ntrac1) ! kg kg-1 s-1
     RTG_CAN (:,kc,   ntke  ) = RTG (:,k,   ntke  ) ! J   s-1 s-1

   end do

! All combined canopy plus resolved layers
   do k = 1, km
     ! nkc+km is top (nkt) combined
     ! nkc+1  is bot       combined
     kc= nkc+k     ! top (nkt) to nkc+1 combined canopy plus resolved model layer

! Height
     zh_can  (:,kc) = zl  (:,k)
     zf_can  (:,kc) = zm  (:,k)

! Pressure & temperature
     prsl_can(:,kc) = prsl(:,k)  ! km  combined canopy plus resolved layers
     prsi_can(:,kc) = prsi(:,k)  ! km  combined canopy plus resolved layers
     T1_CAN  (:,kc) = T1  (:,k)
     DENS_CAN(:,kc) = DENS(:,k)

! Diffusivities
     DKT_CAN (:,kc) = DKT    (:,k)       ! m2 s-1
     DKU_CAN (:,kc) = DKU    (:,k)       ! m2 s-1

! Wind
     WS_CAN  (:,kc) = sqrt(u1(:,k)**2+v1(:,k)**2)       ! m s-1

! Mass tracers
     Q1_CAN  (:,kc, 1:ntrac1) = Q1 (:,k, 1:ntrac1) ! all tracers ntrac1

! TKE tracer
     Q1_CAN  (:,kc,   ntke  ) = Q1 (:,k,   ntke  ) ! ntke=198  TKE tracer

! Humidity
     QV_CAN(:,kc) = Q1(:,k, ntqv) ! ntqv=1

   end do
   prsi_can(:,nkt+1 ) = prsi(:,km+1)  ! km  combined canopy plus resolved layers

! Canopy layers
   do kc = 1, nkc ! 3-nkc canopy layers

! Tendencies
     DU_CAN  (:,kc)    = DU  (:,1)       ! m s-2
     DV_CAN  (:,kc)    = DV  (:,1)       ! m s-2
     TDT_CAN (:,kc)    = TDT (:,1)       ! K s-1

     RTG_CAN (:,kc, 1:ntrac1) = RTG (:,1, 1:ntrac1) ! kg kg-1 s-1
     RTG_CAN (:,kc,   ntke  ) = RTG (:,1,   ntke  ) ! J   s-1 s-1

! Height
     zh_can  (:,kc) = zl  (:,1)
     zf_can  (:,kc) = zm  (:,1)

! Pressure & temperature
     prsl_can(:,kc) = prsl(:,1)  ! km  combined canopy plus resolved layers
     prsi_can(:,kc) = prsi(:,1)  ! km  combined canopy plus resolved layers
     T1_CAN  (:,kc) = T1  (:,1)
     DENS_CAN(:,kc) = DENS(:,1)

! Diffusivities
     DKT_CAN (:,kc) = DKT    (:,1)       ! m2 s-1
     DKU_CAN (:,kc) = DKU    (:,1)       ! m2 s-1

! Wind
     WS_CAN  (:,kc) = sqrt(u1(:,1)**2+v1(:,1)**2)       ! m s-1

! Mass tracers
     Q1_CAN  (:,kc, 1:ntrac1) = Q1 (:,1, 1:ntrac1) ! all tracers ntrac1

! TKE tracer
     Q1_CAN  (:,kc,   ntke  ) = Q1 (:,1,   ntke  ) ! ntke=198  TKE tracer

! Water vapor
     QV_CAN  (:,kc) = Q1(:,1, ntqv) ! ntqv=1

   end do


   end subroutine canopy_levs_init

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   subroutine canopy_levs_run(im, ix, km,           &
              ntrac1, ntqv, ntke,                   & ! in
              RDGAS, PI,                            & ! in ?? units ??
              zi, zl, zm,                           & ! in: 1D    zm(i,k) = zi(i,k+1)
              prsl, prsi, psfc,                     & ! in: 3D 3D 2D  (Pa)
              cfch,                                 & ! in: 2D
              garea, u10m, v10m, fm,fh,             & ! in: 2D
              rbsoil,                               & ! in: 2D
              T2M, Q2M,                             & ! in: 2D
              stress, spd1,                         & !zol, & ! in: 2D
              dv, du, tdt, rtg,                     & ! in: 3D
              U1, V1, T1, Q1,                       & ! in: 3D " 4D
              DENS, dkt, dku,                       & ! in  3D
              FRT_MASK,                             & ! in  2D
              kmod, kcan3,                          & ! out
              zmom_can3, zmid_can3,                 & ! out    zmom_can3 (:, nkt)   zmid_can3(im, nkt)
              sigmom_can, sigmid_can,               & ! out 3D sigmom_can(:, nkt)  sigmid_can(im, nkt)
              ZH_CAN, ZF_CAN,                       & ! out 3D
              PRSL_CAN, PRSI_CAN,                   & ! out 3D prsi_can (:, nkt+1)
              dv_can, du_can, tdt_can, rtg_can,    & ! out: 3D
              T1_CAN,  QV_CAN, DENS_CAN,            & ! out 3D
              WS_CAN, DKT_CAN,  DKU_CAN,            & ! out 3D
              Q1_CAN, Q1_2M)            !out

   use machine , only : kind_phys
! Allocated in mfpbltq_mod:  q1(ix,km,ntrac1)  t1(ix,km) u1(ix,km), v1(ix,km)
   use mfpbltq_mod
!   use physcons, grav => con_g, cp => con_cp, &
!                   rd => con_rd
   use canopy_mask_mod

   IMPLICIT NONE

! Includes:

!...Arguments:

   integer, intent(in)  :: im, ix, km, ntrac1, ntqv, ntke
   real(kind=kind_phys), intent(in) :: RDGAS, PI
! NB. zi   (im, km+1), zl   (im, km),  zm(im,km)
!     prsi (im, km+1), prsl (im, km)
   real(kind=kind_phys), intent(in) :: dv(:,:),    du(:,:),         &
                                      tdt(:,:),   rtg(:,:,:)
   real(kind=kind_phys), intent(in) :: zi(:,:),   zl(:,:), zm(:,:), &
                                     prsi(:,:), prsl(:,:)
   real(kind=kind_phys), intent(in) ::  psfc(:)   ! Pa
   real(kind=kind_phys), intent(in) :: cfch(:), garea(:), u10m(:), v10m(:),  &
                                       spd1(:),stress(:),                    &
                                        t2m(:),   q2m(:),   fm(:),   fh(:),  &
                                     rbsoil(:)

! Allocated in mfpbltq_mod:  q1(ix,km,ntrac1)  t1(ix,km) u1(ix,km), v1(ix,km)
!     ** Q1 is concentration field (including gas and aerosol variables) kg kg-1
   real(kind=kind_phys), intent(in) :: u1(:,:), v1(:,:), t1(:,:), q1(:,:,:)

   real(kind=kind_phys), intent(in) :: dens(:,:), dkt(:,:), dku(:,:)

   real(kind=kind_phys), intent(in) :: FRT_mask(:)

   integer, intent(out) ::          &
               kmod   (:, :)      , &
               kcan3  (:, :)

   real(kind=kind_phys), intent(out) ::  &
! tendencies
             dv_can    (:, :)        , & ! dim(km)
             du_can    (:, :)        , & ! dim(km)
             tdt_can   (:, :)        , & ! dim(km)
! tendencies all gas-phase species & TKE
             RTG_CAN   (:, :, :)     , & ! dim(km )
! met3d arrays
             ZH_CAN    (:, :)        , & ! dim(nkt)
             ZF_CAN    (:, :)        , & ! dim(nkt)
             T1_CAN    (:, :)        , & ! dim(nkt)
             QV_CAN    (:, :)        , & ! dim(nkt)
             WS_CAN    (:, :)        , & ! dim(nkt)
             PRSL_CAN  (:, :)        , & ! dim(nkt)
             PRSI_CAN  (:, :)        , & ! dim(nkt+1)
             DENS_CAN  (:, :)        , & ! dim(nkt)
             DKT_CAN   (:, :)        , & ! dim(nkt)
             DKU_CAN   (:, :)        , & ! dim(nkt)
! all gas-phase species array
             Q1_CAN    (:, :, :)     , & ! dim(nkt)
             Q1_2M     (:, :)        , & ! dim(nkt)
! canopy layers height arrays
             zmom_can3 (:, :)        , & ! dim(nkt+1)    ! Paul's sigmcan(:,nkt)
             zmid_can3 (:, :)        , & ! dim(nkt)      ! Paul's sigtcan(:,nkt)
             sigmom_can(:, :)        , & ! dim(nkt) ~ prsi(:,km+1)
             sigmid_can(:, :)            ! dim(nkt) ~ prsl(:,km)

!...Local arrays:

   character(256) :: errmsg
   integer        :: errflg

   integer(kind=4) :: kcan_top
   real   (kind=kind_phys) :: hcan

   logical      :: sfcflg(im)

   integer(kind=4) :: &
               ka        (im)             , &
               kl        (im)

   real(kind=kind_phys) ::               &
                                zmid3    (km)  , &
                                zmom3    (km)  , & ! Paul's zfull
                                sigmom3  (km+1)  , &
                                z2       (km+1), & ! Paul's    z2(:,chm_nk+1)
                                sigmid2  (km+1), & ! Paul's sigt2(:,chm_nk+1)
             zcan3     (nkc),                    &
             ta_can3   (nkt),   ta3      (km)  , &
             qv_can3   (nkt),   qv3      (km)  , &
             ws_can3   (nkt),   ws3      (km)  , &
             dkt_can3  (nkt),   dkt3     (km)  , &
             dku_can3  (nkt),   dku3     (km)  , &
             prsl_can3 (nkt),   prsl3    (km)  , &
             prsi_can3 (nkt+1), prsi3    (km+1), &
             dens_can3 (nkt),   dens3    (km)  , &
                                mol3     (km)  , &
             klower_can(nkc)

   real(kind=kind_phys) ::                  &
             dxdy      (im),  ustar   (im), &
             ws10m     (im),                &
             zol       (im),  ilmo    (im), &
                    safe_inv_mo_length(im)

!...local variables

   INTEGER          :: i,L

   logical(kind=4) :: flag_error
   integer(kind=4) :: k, kk, kc, k2, II, npass

   real(kind=kind_phys) :: tmp
   real(kind=kind_phys) :: hol, a1, b1, c1, rat

! del:  Minimum allowable distance between a resolved model layer and a canopy layer
!       (fraction of canopy layer height)
   real(kind=kind_phys),    parameter :: del = 0.2
   real(kind=kind_phys),    parameter :: min_kt = 0.1
   real(kind=kind_phys),    parameter :: zfmin=1.e-8
   real(kind=kind_phys),    parameter :: rimin=-100.
   real(kind=kind_phys),    parameter :: karman=0.4            ! von karman constant
   real(kind=kind_phys),    parameter :: THRESHOLD = 1.e06 ! MOL threshold, similar to mach_plumerise

   real(kind=kind_phys)    :: zm2, zr, td, hd, ddel
   real(kind=kind_phys)    :: uh, uspr, wndr, sigw, tl, ktr, kur

!  Assign the fractional heights of the canopy layers (fraction of canopy height)
   real(kind=kind_phys), dimension(3), parameter   :: can_frac = (/1.0, 0.5, 0.2/)

   logical(kind=4)                         :: local_dbg

   local_dbg = (.false.)

   kmod (:,:) = -999
   kcan3(:,:) = -999

! Initializations

! NB. mfpbltq_mod: q1(ix,km,ntrac1) kg kg-1
   Q1_2M (:,   :) = Q1(:,1,:)       ! kg kg-1

! Subset (km combined layers minus top nkc layers)
   do k = 1, km-nkc
     ! km     is top combined subset
     ! nkc+1  is bot combined
     kc= nkc+k     ! 4th from top (nkt) to nkc+1 combined canopy plus resolved model layer

! PBL Tendencies are declared in CCPP_typdefs as dim(im,km) instead (im, nkt)
     DU_CAN  (:,kc)    = DU  (:,k)       ! m s-2
     DV_CAN  (:,kc)    = DV  (:,k)       ! m s-2
     TDT_CAN (:,kc)    = TDT (:,k)       ! K s-1

     RTG_CAN (:,kc, 1:ntrac1) = RTG (:,k, 1:ntrac1) ! kg kg-1 s-1
     RTG_CAN (:,kc,   ntke  ) = RTG (:,k,   ntke  ) ! J   s-1 s-1

   end do

! All combined canopy plus resolved layers
   do k = km, 1, -1 ! top to 1hy model layer
     ! nkc+km is top (nkt) combined
     ! nkc+1  is bot       combined
     kc= nkc+k     ! top (nkt) to nkc+1 combined canopy plus resolved model layer

! Pressure & Temperature
     prsl_can(:,kc) = prsl(:,k)  ! km  combined canopy plus resolved layers
     prsi_can(:,kc) = prsi(:,k)  ! km  combined canopy plus resolved layers
     T1_CAN  (:,kc) = T1  (:,k)
     DENS_CAN(:,kc) = DENS(:,k)

! Diffusivities
     DKT_CAN (:,kc) = DKT    (:,k)       ! m2 s-1
     DKU_CAN (:,kc) = DKU    (:,k)       ! m2 s-1

! Wind
     WS_CAN  (:,kc) = sqrt(u1(:,k)**2+v1(:,k)**2)       ! m s-1

! Mass tracers
     Q1_CAN  (:,kc, 1:ntrac1) = Q1 (:,k, 1:ntrac1) ! all tracers ntrac1

! TKE tracer
     Q1_CAN  (:,kc,   ntke  ) = Q1 (:,k,   ntke  ) ! ntke=198  TKE tracer

! Water vapor
     QV_CAN(:,kc) = Q1(:,k, ntqv) ! ntqv=1

   end do
   prsi_can(:,nkt+1 ) = prsi(:,km+1)  ! km  combined canopy plus resolved layers

! Canopy layers
   do kc = 1, nkc ! 3-nkc canopy layers

     DU_CAN  (:,kc)    = DU  (:,1)       ! m s-2
     DV_CAN  (:,kc)    = DV  (:,1)       ! m s-2
     TDT_CAN (:,kc)    = TDT (:,1)       ! K s-1

     RTG_CAN (:,kc, 1:ntrac1) = RTG (:,1, 1:ntrac1) ! kg kg-1 s-1
     RTG_CAN (:,kc,   ntke  ) = RTG (:,1,   ntke  ) ! J   s-1 s-1

     prsl_can(:,kc) = prsl(:,1)  ! km  combined canopy plus resolved layers
     prsi_can(:,kc) = prsi(:,1)  ! km  combined canopy plus resolved layers
     T1_CAN  (:,kc) = T1  (:,1)
     DENS_CAN(:,kc) = DENS(:,1)

     DKT_CAN (:,kc) = DKT    (:,1)       ! m2 s-1
     DKU_CAN (:,kc) = DKU    (:,1)       ! m2 s-1
     WS_CAN  (:,kc) = sqrt(u1(:,1)**2+v1(:,1)**2)       ! m s-1

! Mass tracers
     Q1_CAN  (:,kc, 1:ntrac1) = Q1 (:,1, 1:ntrac1) ! all tracers ntrac1

! TKE tracer
     Q1_CAN  (:,kc,   ntke  ) = Q1 (:,1,   ntke  ) ! ntke=198  TKE tracer

! Water vapor
     QV_CAN  (:,kc) = Q1(:,1, ntqv) ! ntqv=1

   end do


   DO i = 1, im

      sfcflg(i)= .true.
      if(rbsoil(i) > 0.) sfcflg(i) = .false.

      dxdy(i) = garea( i ) !  dx*dy ~1.6E+8 m2

      ustar(i) = sqrt(stress(i))
!      ws10m(i)  = sqrt(u10m(i)**2+v10m(i)**2)

!> ## Compute Monin-Obukhov similarity parameters
!!  - Calculate the Monin-Obukhov nondimensional stability paramter, commonly
!!    referred to as \f$\zeta\f$ using the following equation from Businger et al.(1971) \cite businger_et_al_1971
!!    (eqn 28):
!!    \f[
!!    \zeta = Ri_{sfc}\frac{F_m^2}{F_h} = \frac{z}{L}
!!    \f]
!!    where \f$F_m\f$ and \f$F_h\f$ are surface Monin-Obukhov stability functions calculated in sfc_diff.f and
!!    \f$L\f$ is the Obukhov length.
      zol(i) = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
      if(sfcflg(i)) then
         zol(i) = min(zol(i),-zfmin)
      else
         zol(i) = max(zol(i),zfmin)
      endif
! Inverse of Monin-Obukhov length
      ilmo(i) = 1./zol(i)

!!!! Non-Canopy columns
   IF (FRT_mask(i) <= 0.) THEN

!!!!! Start all columns!!!!! canopy & non-canopy (canopy columns are overwritten below if FRT_MASK > 0.)
      do k = 1, km     ! from bottom to top
         II = km + 1 - k  ! from top to bottom of resolved model layers
!!! Paul's zmom is our zmom
! zmom(1)     = ZFULL(km) is top    model  layer height
! zmom(km) = ZFULL(1)     is bottom model  layer height
         ! NB. zm(:,k) = zi(:,k+1)
         ! zmom3(II) = zi(i,k) ! ZFULL(i,k) Mar24, 2025 replace zi with zm
         zmom3(II) = zm(i,k) ! ZFULL(i,k)
!! Heights of the original model layers for the canopy columns are extracted to the zmom array.

         ! Create temperature & humidity array on reversed layer order for interpolation
         ta3  (II) = T1  (i,k)               ! K
         qv3  (II) = Q1  (i,k,1)      ! 1=water vapor kg kg-1
         prsl3(II) = PRSL(i,k)  ! Pa  mean layer pressure
         dens3(II) = DENS(i,k)  ! kg m-3
         ws3  (II) = sqrt(u1(i,k)**2+v1(i,k)**2)  ! rename wspd3 ???
         dkt3 (II) = DKT (i,k)  ! m2 s-2
         dku3 (II) = DKU (i,k)  ! m2 s-2
      end do ! k = 1, km     ! from bottom to top

      do k = 1, km+1     ! from bottom to top
         II = (km + 1) + 1 - k  ! from top to bottom of resolved model layers

         prsi3  (II) = PRSI(i,k)  ! Pa  air pressure at model layer interfaces
! ! [pgr] surface air pressure meta var
         sigmom3(II) = PRSI(i, k) / psfc(i) ! PRES_FULL
      end do ! k = 1, km+1

!  First, carry over original model values for the matching layers
      do k = 1, km ! from bottom to top of resolved model layers
        ! kmod(1)     is 1  top model  layer
        ! kmod(km) is 65 top canopy layer (modified after mono adj.)
         kk = kmod(i,k)

! to do
!        zmom_can3 (i,kk) = zmom3  (k) ! full layer height [m]
         sigmom_can(i,kk) = sigmom3(k) !

         ta_can3  (kk) = ta3  (k)   !          TA  (i, k) ! temperature [K]
         qv_can3  (kk) = qv3  (k)   ! Met_Data%QV  (i, k) ! spec. humidity
         prsl_can3(kk) = prsl3(k)   ! Met_Data%PRES(i, k) ! Pa
         prsi_can3(kk) = prsi3(k)   !
         dens_can3(kk) = dens3(k)   ! Met_Data%DENS(i, k) ! kg m-3
         ws_can3  (kk) = ws3  (k)   !                     ! m s-1
         dkt_can3 (kk) = dkt3 (k)   !          DKT (i, k) ! m2 s-2 atmos. thermal diffus.
         dku_can3 (kk) = dku3 (k)   !          DKU (i, k) ! m2 s-2 atmos. momentum diffus.
      end do

      do kc = 1, nkc       ! from top to bottom of canopy layers
         ! kk = 65 = kcan3(1) = km + 1
         ! kk = 66 = kcan3(2) = km + 2
         ! kk = 67 = kcan3(3) = km + 3
         kk = kcan3(i,kc)

!        zmom_can3 (i,kk) = zmom3  (km) ! full layer height [m]
         sigmom_can(i,kk) = sigmom3(km) !

         ta_can3  (kk) = ta3  (km)   !          TA  (i, k) ! temperature [K]
         qv_can3  (kk) = qv3  (km)   ! Met_Data%QV  (i, k) ! spec. humidity
         prsl_can3(kk) = prsl3(km)   ! Met_Data%PRES(i, k) ! Pa
         prsi_can3(kk) = prsi3(km)
         dens_can3(kk) = dens3(km)   ! Met_Data%DENS(i, k) ! kg m-3
         ws_can3  (kk) = ws3  (km)   !                     ! m s-1
         dkt_can3 (kk) = dkt3 (km)   !          DKT (i, k) ! m2 s-2 atmos. thermal diffus.
         dku_can3 (kk) = dku3 (km)   !          DKU (i, k) ! m2 s-2 atmos. momentum diffus.
      end do
! Lower interface at surface
      prsi_can3 (   nkt+1) = prsi3(km+1)
      sigmom_can(i, nkt+1) = 1.0

!!!!! End all columns!!!!!

   ! Continuous forest canopy
   ELSE IF (FRT_mask(i) > 0.) THEN

!      print*, 'CANOPY_LEVS: ZOL ILMO= ', i, zol(i), ilmo(i)

      hcan = cfch( i )
!!! Extract the canopy height (FCH)

! Generate initial canopy levels, as altitude above sea level
!
      do kc = 1, nkc
         zcan3(kc) = hcan * can_frac(kc)  ! Paul's hc is our hcan
                                          ! Paul's zcan is our zcan3
!!! Set the initial values of the heights of the inserted canopy layers to hc, 0.5 hc, and 0.2 hc
!!!
!!! NB. zcan3(1) is hc, top of canopy
!!!     zcan3(2) is 0.5 * hc
!!!     zcan3(3) is 0.2 * hc (bottom canopy level)

!        print*,'canopy_levs: ZCAN = ', i, kc, zcan3(kc)
      end do

! 1 = bottom (1st) model layer
! km= top model layer
      do k = 1, km     ! from bottom to top
         II = km + 1 - k  ! from top to bottom of resolved model layers
! zl is height of layer center
! zmid3(1)  = zl(km) is top    model  layer height
! zmid3(km) = zl(1)  is bottom model  layer height
        ! Paul's zt is our zmid
        zmid3(II) = ZL(i,k) ! mid layer height [m]
!!! Heights of the original model layers for the canopy columns are extracted to the zmid array.

!       write(errmsg,*) 'canopy_levs: ZMID = ', i, II, zmid3(II)

        ! Paul's sigt2 is our sigmid2
        sigmid2(II) = prsl(i,k)/ psfc(i)

! 65 1.0   (surface) !! Set to 1 here !!
! 64  0.997329666888429  (1hy model layer)
! 63  0.994572224115356
! 62  0.987953350646348
! 61  0.980671961372880
!  ...
! 4  2.651067835355327E-003
! 3  1.751135612532654E-003
! 2  9.570774376723687E-004
! 1  3.757488135785848E-004 (top model )
!        print*,'canopy_levs: sigmid2= ', i, II, sigmid2(II)
      end do
      sigmid2(km+1) = 1.0

      do k = 1, km     ! from bottom to top
         II = km + 1 - k  ! from top to bottom of resolved model layers
!!! Paul's zmom is our zmom
! zmom(1)     = ZFULL(km) is top    model  layer height
! zmom(km) = ZFULL(1)     is bottom model  layer height
         ! NB. zm(:,k) = zi(:,k+1)
         ! zmom(II) = zi(i,k) ! ZFULL(i,k) Mar24, 2025 replace zi with zm
         zmom3(II) = zm(i,k) ! ZFULL(i,k)
!! Heights of the original model layers for the canopy columns are extracted to the zmom array.

         ! Create temperature & humidity array on reversed layer order for interpolation
         ta3  (II) = T1  (i,k)               ! K
         qv3  (II) = Q1  (i,k,1)      ! 1=water vapor kg kg-1
         prsl3(II) = PRSL(i,k)  ! Pa  mean layer pressure
         dens3(II) = DENS(i,k)  ! kg m-3
         ws3  (II) = sqrt(u1(i,k)**2+v1(i,k)**2)
         dkt3 (II) = DKT (i,k)  ! m2 s-2
         dku3 (II) = DKU (i,k)  ! m2 s-2

! From satmedmfvdifq.F:
! MOL  = zol(i)/zl(i,k)  !Monin-Obukhov Length in layer
! ZL is mid layer height [m]
         mol3(II) = zol(i)/ZL(i,k)  !Monin-Obukhov Length in layer
      end do

      do k = 1, km+1     ! from bottom to top
         II = (km + 1) + 1 - k  ! from top to bottom of resolved model layers

         prsi3  (II) = PRSI(i,k)  ! Pa  air pressure at model layer interfaces
! Paul's SIGM does not include surface layer lower interface (1.0) !!!
         sigmom3(II) = PRSI(i, k)/ psfc(i) ! PRES_FULL(i, k) / psfc(i)

! prsi (km+1) =>  prsi3( 1)   Top    model layer upper interface
!
! 65  1.00000000000000  93074.3428508980 mb (km+1) surface bottom model layer interface
! 64  0.994671010591796 92578.3506636700 mb
! 63  0.988632406984791 92016.3116012109
! ...
! 3  1.367636992545316E-003 137.789993286133
! 2  6.376847405122714E-004 64.2470016479492
! 1  1.985103504149681E-004 20.0000000000000 mb (top model layer)
!
!         print*,'canopy_levs: sigmom3= ', i, II, sigmom3(II),prsi3 (II)

      end do

!!! Find the resolved model level which lies above the top of the forest canopy,
!!! in each canopy column.  Usually the canopy is within the km or km-1
!!! level of the original model structure.
!
! The model level above the tallest canopy in grid
      kcan_top = 2         ! initialize to 2nd top model layer
      do L = km, 3, -1  ! from bottom to top model layer (going up)
         ! Mid-layer height m, zmid
         if (zmid3(L) > hcan) then ! Paul's zt is our zmid
            kcan_top = L - 1  ! level above the tallest canopy
            exit
         end if
      end do
! kcan_top = 62 or 63
!      print*,'canopy_levs: kcan_top = ', i, kcan_top

! MV2D_ILMO: Aggregated Inverse of Monin-Obukhov length
!  Setup of Monin-Obhukov Length similar to plumerise for upper limit:
! from satmedmfvdifq.F: ! MOL  = zol(i)/zl(i,k)  !Monin-Obukhov Length in layer
      safe_inv_mo_length(i) = ilmo(i)
      if (abs(ilmo(i)) > THRESHOLD) then
         safe_inv_mo_length(i) = sign(THRESHOLD, ilmo(i))
      end if
!
! Adjust the canopy levels: we don't want canopy levels to get closer than del (0.2m)
! to the model levels to prevent possible differencing errors in the diffusion.
! If zcan3 > zmid3 but is too close to zmid3, move zcan3 up by ddel.  If zcan3 < zmid3
! but is too close to zmid3, move zcan3 down by ddel.  The net result will be that the
! canopy levels are never closer than del from the original model levels.
      do k = kcan_top, km! from model layer above the canopy to bottom of model layer
         do kc = 1, nkc   ! from top to bottom of canopy
            if (abs(zmid3(k) - zcan3(kc)) < del) then
               ddel = max(0.0, del - abs(zcan3(kc) - zmid3(k)))
               zcan3(kc) = zcan3(kc) + sign(ddel, zcan3(kc) - zmid3(k))

!!! The reason why this section is necessary:  while it would be preferable for
!!! the canopy levels to stick with values of hc, 0.5 hc and 0.2 hc, somewhere in a
!!! large domain, there may be an overlap where one of these canopy levels is very
!!! close to or on top of an existing model level.  Which we dont want!
!!! What is done here, if the canopy levels come within "ddel" of an original model level,
!!! is to shift the canopy level in question a bit, to avoid overlaps.
            end if
         end do
      end do

!!! Starts the creation of the local array with the heights of the thermodynamic
!!! levels (layer midpoints) for the combined canopy + no canopy layers.
!!! Note that zmid_can at this point does not have these layers sorted in the
!!! correct order - the canopy layers have been tacked onto the bottom of the
!!! zmid_can array, but the values of zmid_can are not monotonically increasing with
!!! decreasing height index.
!
!  Set the initial values of the combined height array:
!
! Note that here, zmid_can is created, but the heights within each column have
! yet to be sorted to rearrange the layers in the correct order.
      do k = 1, km ! from top to bottom model layers
         zmid_can3(i,k)    = zmid3(k)
         ! Paul's zthrmcan is our zmid_can
!
      end do

! Add zcan3 additional thermo levels into zmid_can array for later sorting
      do kc = 1, nkc ! from top to bottom canopy layers
         zmid_can3 (i,km+kc) = zcan3(kc)
      end do

!
!  Determine locations of canopy and resolved model levels within
!  the combined array for the canopy columns:
!
!!! This section sorts the zmid_can array to make sure that the new layers are
!!! all ordered so they monotonically increase with decreasing height.

! Top canopy layer height (km+1) is higher than the bottom model layer height (km)
      if (zmid_can3(i, km) < zmid_can3 (i,km+1)) then
!
!  Non-trivial case:  the ancilliary and original array levels intermingle.
!  Sort the combined height array to get the right order of the the heights:
!
!  zmid_can is the height locations of the combined array, which needs to be  sorted:
!  since there are only NKC levels in the canopy, and both zcan3 and z
!  decrease monotonically, only nkc+1 passes are needed to sort the combined array:
         do npass = 1, nkc+1
            flag_error = .false.
            do k = nkt, 2, -1
! Top    canopy layer height (nkt-2) is larger than the bottom model  layer height (nkt-3 = km)
! Middle canopy layer height (nkt-1) is larger than the top    canopy layer height (nkt-2)
! Bottom canopy layer height (nkt)   is larger than the middle canopy layer height (nkt-1)
               if (zmid_can3(i, k) > zmid_can3(i, k-1)) then
!  The combined array heights are out of order, sort them:
                  tmp = zmid_can3(i, k-1)
                  zmid_can3(i, k-1) = zmid_can3(i,k)
                  zmid_can3(i, k)   = tmp
                  flag_error = .true.
               end if
            end do
         end do
         if (flag_error) then
!           write(errmsg,*) 'NKC+1 passes insufficient to sort canopy array '
!           write(errmsg,*) 'in can_levs_defn.F90.  Scream and die.'
! ABORT!
            return
         end if
      end if

!
!  Heights in zmid_can should now be monotonically decreasing.

! Print
!      do k = nkt, 1, -1 ! sfc to top model layer
! 67   3.71699981689453
! 66   9.29249954223633
! 65   18.5849990844727
! 64   22.5893670351600
!         print*,'canopy_levs: zmid_can = ', i, k, zmid_can3(i, k)
!      end do

!  Next, identify the locations of the vertical levels in the combined
!  array relative to the resolved model array and canopy array
!
!!! Now that the heights in zmid_can are in the right order, we can use them to
!!! identify the values of kcan and kmod:  the vertical locations of the canopy and
!!! original model layers in the augmented canopy layer code.
      do kc = 1, nkc         ! from top to bottom canopy layers
         do kk = nkt, 1, -1  ! from bottom to top of combined canopy and resolved model levels
            if (zmid_can3 (i, kk) == zcan3(kc)) then
               kcan3(i,kc) = kk
               exit
            endif
         end do
      end do

! k=1        is top    model layer
! k=km       is bottom model layer
      do k = 1, km ! from top to bottom model layers
         do kk = k, nkt! from bottom to top of combined resolved plus canopy layers

! zmid_can3(1)    = zmid3(1)     is top    model  layer height
! ...
! zmid_can3(km)= zmid3(km) is bottom model  layer height
            if (zmid_can3(i, kk) == zmid3(k)) then

! kmod(1)       is 1      , top model  layer
! kmod(km-1) is km-1, 2nd model  layer
! kmod(km)   is          top canopy layer (modified after monotonic adj.)
               kmod(i,k) = kk
               exit
            endif
         end do
      end do

      if (local_dbg) then
      do kc = 1, nkc
         if (kcan3(i,kc) < 1) then
!           write(errmsg,*) 'get_can_levs: kcan undefined: ', kc, kcan3(i,kc)
            !ABORT
            return
         end if
      end do
      do k = 1,km
         if (kmod(i,k) < 1) then
!           write(errmsg,*) 'get_can_levs: kmod undefined: ',k, kmod(i,k)
            !ABORT
            return
         end if
      end do
      end if


!  Create the corresponding momentum height array
!
!  The original methodology adopted made use of the at2m array and the  thermodynamic heights determined above.
!  However, this methodology resulted in momentum levels which did not match the original model levels
!  above the region modified for canopy layers.  Here, the thermodynamic layers will be used to
!  (1) Determine whether the original model and canopy thermodynamic layers coincide, and if so,
!  (2) Use the existing model layer values for the momentum layers, while if not,
!  (3) Assign the new momentum layers as being 1/2 way between the canopy layers
! Note that these changes only exist inside the chemistry part of GEM-MACH and do not affect the model physics
!!!
!!! Create the momentum height (layer interface) array.  The original momentum layers are used above the canopy height.
!!! Below the canopy height, the "momentum"layers are assumed to be ½ way between the thermodynamiclayers.

! Default case:  all added canopy thermodynamic layers are below the lowest resolved model thermodynamic layer
! kcan_top is either 2nd or 3rd (63 or 62) resolved model layer
      do k = 1, kcan_top - 1 ! from top model layer to model layer above the canopy
!          zmom(1) is top    model  layer height
! zmom(kcan_top-1) is model layer above the canopy < 234.061m
         zmom_can3(i,k) = zmom3(k) ! full layer height [m]
      end do

      ka(i) =km
      inner0: do k = kcan_top, km-1  ! from resolved model layer above the canopy to top model layer
!  Starting from the top, scan down through the original and combined mid layer heights, to see when
!  they first deviate from each other
         !Paul's zthrmcan is our zmid_can
         !Paul's zt is our zmid
         if (zmid_can3(i,k) == zmid3(k) .and. zmid_can3(i,k+1) == zmid3(k+1)) then
            ! Paul's zmom is our zmom
            ! Paul's zmomcan(nkt+1) is our zmom_can
            zmom_can3(i,k) = zmom3(k)  ! full layer height [m]
         else
            ka(i) = k
            exit inner0
         end if
      end do  inner0

! ka is 63 or 64
!      print*,'canopy_levs: ka = ', i, ka(i)

! ka  is the lower-most layer for which the combined layer zmom_can = zmom resolved model layer
      ! Paul's zmom is our zmom
      ! Paul's zmomcan is our zmom_can
      zmom_can3(i,ka(i)) = zmom3(ka(i))
      do k = ka(i)+1, nkt! from ka to bottom combined canopy and resolved layers
         zmom_can3(i,k) = (zmid_can3(i,k-1) + zmid_can3(i,k)) * 0.5

      end do
! Oct31:      zmom_can3(i, nkt+ 1) = 0.


! Print
!      do k = nkt, 1, -1 ! sfc to top model layer
! 67   6.62900018692017
! 66   14.2050004005432
! 65   21.0653651654053
! 64   46.3814595935061  1hy
! 63   99.2328891021972  2hy
!         print*,'canopy_levs: zmom_can = ', i, k, zmom_can3(i, k), zmom3 (k)
!      end do

!########################################################################

! create original model arrays of z and sigma-t which include the surface, to
!  allow interpolation:
      !Paul's sigtcan is our sigmid_can
      sigmid_can(:,:) = 0.0
      do k = 1, km! from top to bottom of resolved model layers

! zmid3(1)  is top    model  layer height
! zmid3(km) is bottom model  layer height
         z2(k)    = zmid3(k)

!  Fill in the thermodynamic sigma levels (Pre-existing levels first):
! kmod(1)     is 1  top  (last) model layer
! kmod(km) is 64 bottom (1st)  model layer
         kk = kmod(i,k)
         sigmid_can(i, kk) = sigmid2(k)

!     sigmid_can             zmid_can3
! 1  3.875425449149410E-004 54904.9550581240 m
! 2  9.844331193192971E-004 47732.0690652646 m
! ...
! 62  0.985167158577051       125.175103771062 m
! 63  0.991717417180879        70.5363577077242
! 64  0.997329666888429        22.4844313034714
!
!         print*,'canopy_levs: sigmid_can = ', i, kk, sigmid_can(i, kk), &
!                                                      zmid_can3(i, kk)

      end do
      klower_can(:) = -999
      z2(km+1) = 0.0

!
!  fill in the remaining sigma levels by interpolating in z:
      do kc = 1, nkc   ! from top to bottom canopy layers
         do k2 = kcan_top, km+1  ! from resolved model layer above the canopy to top model layer
            if (zcan3(kc) > z2(k2) .and. zcan3(kc) <= z2(k2-1)) then

! k2 is either 64 or 65
! 64  0.997509580701422 0.991549245511511       5.960335189910571E-003
!    23.4420505707344 73.6016275069086        -50.1595769361742
!    23.6420505707344                         -49.9595769361742
!
! 64  0.997359509095134 0.991637283835972       5.722225259162883E-003
!    23.5479167685719  73.9801156184914       -50.4321988499195
!    23.7479167685719                         -50.2321988499195
!
! 65   1.00000000000000        0.997352976969389       2.647023030610929E-003
!      0.000000000000000E+000 22.3611756580077       -22.3611756580077
!      2.73199996948242                              -19.6291756885253
!
! 65   1.00000000000000         0.997352976969389       2.647023030610929E-003
!      0.000000000000000E+000  22.3611756580077       -22.3611756580077
!     13.6599998474121                                 -8.70117581059563
!
!               print*, 'canopy_levs: sigmid_can (1) = ', i, k2,         &
!                sigmid2(k2), sigmid2(k2-1), sigmid2(k2) - sigmid2(k2-1),&
!                     z2(k2),      z2(k2-1),      z2(k2) -      z2(k2-1),&
!                  zcan3(kc),                  zcan3(kc) -      z2(k2-1)

! Interpolate in sigma
               sigmid_can(i, kcan3(i,kc)) = sigmid2(k2-1)  +          &
                                 (sigmid2(k2) - sigmid2(k2-1)) /      &
                                   (   z2(k2) -    z2(k2-1)) *        &
                                   (zcan3(kc) -    z2(k2-1))

! Store grid locations for use in later interpolations
               klower_can(kc) = k2
            end if

         end do ! do k2=kcan_top, km+1

! Print
!
! kcan3    sigmid_can             zmid_can3
! 65 0.999628269764443        3.13000011444092
! 66 0.999814134882221        1.56500005722046
! 67 0.999925653952889       0.626000022888184
!
! 65 0.997117582813635        24.1049995422363
! 66 0.998648976933277        11.8999996185303
! 67 0.999459590773311        4.75999984741211
!
!          print*,'canopy_levs: sigmid_can (2) = ', i, kc, kcan3(i,kc),  &
!                                            sigmid_can(i, kcan3(i,kc)), &
!                                             zmid_can3(i, kcan3(i,kc))
!
!
         if (klower_can(kc) < 1) then
!           write(errmsg,*) 'get_can_levs:  klower_can is unassigned at i, kc: ', i, kc
!           write(errmsg,*) 'get_can_levs:  zcan3(kc): ',zcan3(kc)
            do kk = kcan_top, km+1
!              write(errmsg,*) 'get_can_levs: kk z2(kk) which should bracket the above zcan3: ',kk, z2(kk)
            end do
            do kk = 1, km+1
!              write(errmsg,*) 'get_can_levs:  kk z2(kk) full set of z2 values: ', kk, z2(kk)
            end do
            do kk = 1,nkc
!              write(errmsg,*) 'get_can_levs:  kc zcan3(kc) hcan fr(kc) for full set of zcan3 values: ',kk, zcan3(kk), hcan, can_frac(kk)
            end do
            return
         end if
      end do

! NB.
! klower_can(1) is 64 or 65
! klower_can(2) is 65 except for individual grid points near West coast
! klower_can(3) is 65 uniformly
!
!
   if (local_dbg) then
!  Check on klower_can for NaN or out of bounds:
   do kc = 1,nkc
      if ((klower_can(kc) /= klower_can(kc)) .or. &
          (klower_can(kc) <= 0)              .or. &
          (klower_can(kc) > km+ 1) ) then
!        write(errmsg,*) 'get_can_levs: klower_can after creation NaN or <=0 or >km+1 : ', &
!                      kc, klower_can(kk)
         return
      end if
   end do
   end if
!
!  Create sigma coordinate  momentum levels:
!
!  As above, the existing momentum levels and the canopy values are used to create SIGM levels
!
!  (1) Determine whether the original model and canopy thermodynamic layers coincide, and if so,
!  (2) Use the existing model layer values for the momentum layers, while if not,
!  (3) Assign the new momentum layers as being 1/2 way between the canopy layers
! Note that these changes only exist inside the chemistry part of GEM-MACH and do not affect the
! model physics


!      Default case:  all added canopy half layers are
!      below the lowest resolved model half layer
      ka(i) = km
      inner2:   do k = 1, km-1
         if (sigmid_can(i, k) == sigmid2(k) .and. sigmid_can(i, k+1) == sigmid2(k+1) ) then
            sigmom_can(i, k) = sigmom3(k)
         else
            ka(i) = k
            exit inner2
         end if
      end do inner2
! ka  is the last layer for which sigmom_can= sigmom3(k)
      sigmom_can(i, ka(i)) = sigmom3(ka(i))   ! Jul23: sigmid2(ka(i))
      do k = ka(i)+1,nkt
         sigmom_can(i, k) = (sigmid_can(i, k-1) + sigmid_can(i, k)) * 0.5
      end do
! Jul24, 2025
      sigmom_can(i, nkt+1) = 1.0

! Print
!      do k = 1,nkt+1 ! from top to bottom

! 1  1.985103504149681E-004  prsi3(1) = 20.0000000000000 mb
! 2  6.376847405122714E-004
! ...
!
! 62  0.981799237332539
! 63  0.988632335800729
! 64  0.994671160237943
! 65  0.997541255229605
! 66  0.998374138576117
! 67  0.999241264668854
! 68  1.0   set to 1.0 above
!
!        print*,'canopy_levs: sigmom_can =',i, k, sigmom_can(i, k)
!      end do ! nkt+1


!
!  Next, do a sort of all of the variables in the original METV3D array into canopy.  Note that
!  the declaration of the met arrays for the new canopy subdomain has occurred earlie in the code.
!  Three-D variables are a bit more complicated, in that one must make decisions regarding
!  the values of the met variables in the canopy region.
!  The code which follows is based on chm_load_metvar.ftn90
!
!  First, carry over original model values for the matching layers
   do k = 1, km ! from bottom to top of resolved model layers
      ! kmod(1)     is 1  top model  layer
      ! kmod(km) is 65 top canopy layer (modified after mono adj.)
      kk = kmod(i,k)
      ta_can3  (kk) = ta3  (k)   !          TA  (i, k) ! temperature [K]
      qv_can3  (kk) = qv3  (k)   ! Met_Data%QV  (i, k) ! spec. humidity
      prsl_can3(kk) = prsl3(k)   ! Met_Data%PRES(i, k) ! Pa
      prsi_can3(kk) = prsi3(k)   !
      dens_can3(kk) = dens3(k)   ! Met_Data%DENS(i, k) ! kg m-3
      ws_can3  (kk) = ws3  (k)   !                     ! m s-1
      dkt_can3 (kk) = dkt3 (k)   !          DKT (i, k) ! m2 s-2 atmos. thermal diffus.
      dku_can3 (kk) = dku3 (k)   !          DKU (i, k) ! m2 s-2 atmos. momentum diffus.

! Print
! (km+1) (68=nkc+km +1) prsi3( 1)   Top    model layer upper interface  prsi_can3(1)
! i = 1
! 1   20.0000000000000
! 2   64.2470016479492
! 3   137.789993286133
!...

! 62          62    96311.7483321220        96981.9123946220
! 63          63    96981.9123946220        97574.2952071220
!             64    --> in kcan3 loop: 64 97551.5096832975
! 64          65    97574.2952071220        98097.0373946220
!
!      print*,'canopy_levs: prsi_can3 kmod=', i, k, kk, prsi_can3(kk), prsi3(k+1)

   end do ! km

!----------------------------------------------------------------------------
!  Canopy region:  next, go through each variable to work out canopy values.
!
!  (1) Do those variables for which special canopy formulae will NOT be used:
      do kc = 1, nkc       ! from top to bottom of canopy layers

!  Each of the following 2 variables have a screen height (2m) value in the 2D met arrays
!          Temperature:         TA, T2M
!          Specific humidity:   Q,  Q2M

!   kcan3(1) = 65
!   kcan3(2) = 66
!   kcan3(3) = 67
         kk = kcan3(i,kc)
         if (klower_can(kc) <= km) then
!  Level is above first resolved model level

            k2 = klower_can(kc)
            zm2 = (zcan3(kc) - z2(k2-1)) / (z2(k2) - z2(k2-1))

            td = ( ta3(k2)  - ta3(k2-1)) * zm2
            hd = ( qv3(k2)  - qv3(k2-1)) * zm2
            ta_can3(kk) = ta3(k2-1) + td
            qv_can3(kk) = qv3(k2-1) + hd

         else
!  Level is below first resolved model level

            if (zcan3(kc) - z2(km+1) >= 2.0) then
            !  Level is below first resolved model level but above screen height

               zm2 = (zcan3(kc) - z2(km+1) - 2.0) / (z2(km) - z2(km+1) - 2.0)

               td = (ta3(km)  - T2M( i ) )  * zm2
               hd = (qv3(km)  - Q2M( i ) )  * zm2
               ta_can3(kk)  = T2M( i ) + td
               qv_can3(kk)  = Q2M( i ) + hd

            else
            ! Level in canopy is below screen height; assume constant values below screen height

               ta_can3(kk)  = T2M( i ) ! 2-m  temperature [K]
               qv_can3(kk)  = Q2M( i ) ! 2-m  spec. humidity
            end if

         end if

!  Evaluate the air density in canopy columns using values determined above
!
! NB. PRSL  is air pressure on ZL, formerly ZH (mid-layers)
!     PRSI  is air pressure on ZI, formerly ZF (interfaces)
!     psfc   is surface air pressure psfc

! get pressure from sigma levels in Pa
         prsl_can3(kk) = sigmid_can(i, kk) * psfc(i)  ! ~zl mid-layers centers
         prsi_can3(kk) = sigmom_can(i, kk) * psfc(i)  ! ~zm/zi layers interfaces

! Print
! 1 64 97551.5096832975
!   65 --> in kmod loop : 65    97574.2952071220
! 2 66 97892.5615950123
! 3 67 97999.3464530241
!
!         print*,'canopy_levs: prsi_can3 kcan3=', i, kc, kk, prsi_can3(kk)


! aqm_methods: dens: buffer(k) = stateIn % prl(c,r,l) / ( rdgas * stateIn % temp(c,r,l) )
         ! dens_can3(1)       is top model  layer
         ! ...
         ! dens_can3(km)   is 1hy model  layer
         ! dens_can3(km+1) is top canopy layer
         ! dens_can3(nkt)   is 1st canopy layer
         dens_can3(kk) = prsl_can3(kk) / ( RDGAS * ta_can3(kk))  ! kg m-3


!  The following variables are assumed to have uniform values throughout the
!  lowest resolved model layer:
!
!    Cloud liquid water mass mixing ratio (QCPLUS)
!    Total cloud fraction (FTOT)
!    Stratospheric cloud fraction (FXP)
!    Convective cloud fraction (FDC)
!    Total liquid water flux (RNFLX)
!    Total solid water flux (SNOFLX)
!    Precipitation evaporation (FEVP)
!    Cloud to rain collection tendency (PPRO)
!  Search over the original model layers (k).  Note that the outer loop above this
!  one is over the canopy layers kc:  we are looking for the values to assign the
!  canopy layers in the combined canopy+resolved scale space.  For these variables,
!  the resolved scale values will be used, hence the aim is to determine the
!  resolved scale layer in which the canopy layer resides, and assign the
!  corresponding values to the locations of the canopy layers in the combined
!  canopy + resolved scale space (kk).


      end do ! kc = 1,nkc
! Surface layer lower interface
      prsi_can3(nkt+1) = prsi3(km+1)

   if (local_dbg) then
! Several checks for suspicious values:
      do kk = 1,nkt
         if ( ta_can3(kk) < 150.0) then
            write(errmsg,*) 'get_can_levs:  suspicious temperature detected in get_can_levs after creation (kk value): ',&
                        i, kk, ta_can3(kk)
            do kc = 1, nkc
               write(errmsg,*) 'get_can_levs: value of zcan(kc) z2(km+1) and difference  at this value of ic for kk: ',&
                            kc,' are: ',zcan3(kc),z2(km+1), zcan3(kc)-z2(km+1)
            end do

            do k = 1, nkt
               write(errmsg,*) 'get_can_levs: value of zmid_can for = ', i,' at k = ',k,' is: ',zmid_can3(i,k)
            end do

            do kc = 1,nkc
               write(errmsg,*) 'get_can_levs:  values of kcan zcan and original zcan for = ', i,' at kc = ',kc,' are: ',&
                           kcan3(i,kc), zcan3(kc), hcan * can_frac(kc)
            end do

            do k = 1,km
               write(errmsg,*) 'get_can_levs:  value of kmod and z for = ', i,' at k = ',k,' are: ',kmod(i,k), zmid3(k)
            end do

            do kc = 1,nkc
               write(errmsg,*) 'get_can_levs: value of klower_can at this grid point for kc: ',kc,' is: ',klower_can(kc)
            end do

            return
         end if
      end do
   end if

!  (2) For the last few variables, the value at the lowest resolved model layer and typical profiles for that variable
!  within the canopy will be used to create the canopy values:
      do kc = 1, nkc
         kk = kcan3(i,kc)
!  Ratio of lowest model level to canopy height:
!
         zr = (zmid3(km) - z2(km+1)) / hcan
!
!  Horizontal wind and KT profiles are from Raupach, Quarterly Journal
!  of the Royal Meteorological Society, vol 115, pp 609-632, 1989, examples
!  from page 626, equations (48) through (51).
!
!  Wind speed (equation 51), assumed to scale similarly in each horizontal dimension:
!
!   U(z) = ustar/karman * ln((z - d) / z0), where
!   k = 0.4
!   d = 0.75 hc
!   z0 = 0.07530 hc
! The next few lines calculate the average value of u(z), v(z), Raupach's eqn 51,
! at the first resolved level model height
! Paul's UE is our ustar,  surface friction velocity
         uh = ustar(i) * 3.0
         if (zr >= 1.0) then
         ! Paul's zt is our zmid (i.e. zmid(km) is zt(i,chm_nk))
         ! Paul's hc is our hcan
            uspr = ustar(i) / karman * &
                   alog((zmid3(km) - z2(km+1) - 0.75 * hcan) / &
                   (0.07530 * hcan))
         else
            uspr = uh * exp(- 2.0 * ( 1.0 - zr))
         end if
!  wndr is the ratio of the wind to Raupach's average us(), eqn 51.
!  This is used to scale the wind speed with height values from eqn 51 to the current grid square
         ! Paul's WS(nk) is our spd1, wind speed at lowest model level m s-1
         wndr = spd1(i) / uspr
!  Using Raupach's formulae for wind speed, multiplied by the above ratio, for the canopy layers:
!
         zr = (zcan3(kc) - z2(km+1)) / hcan
         if (zr >= 1.0) then
            uspr = alog((zcan3(kc) - z2(km+1) - 0.75 * hcan) / &
                   (0.07530 * hcan)) * ustar(i)
         else
            uspr = uh * exp(- 2.0 * (1.0 - (zcan3(kc) - z2(km+1)) / hcan))
         end if

         ws_can3(kk) = wndr * uspr
!
!  Coefficients of diffusivity:
!  Find value of K at first model level from raupach's sigw and TL formulae (eqns 48, 49)
         zr = (zmid3(km) - z2(km+1)) / hcan
!  Gradient in stability under the canopy is reduced for higher stability conditions
!  in accord with Shaw, den Hartog and Neumann, BLM 45, 391-409, 1988, Fig 16.
         ! Paul's zl is our hol (as in satmedmfvdifq.F)
         hol = hcan * safe_inv_mo_length(i)
! Unstable:
         if(hol < -0.1) then
             a1 = 0.75
             b1 = 0.5
             c1 = 1.25
         end if
! Neutral:
         if(hol >= -0.1 .and. hol < 0.1) then
             a1 = 0.625
             b1 = 0.375
             c1 = 1.0
         end if
! Stable:
         if(hol >= 0.1 .and. hol < 0.9) then
             rat = 4.375 - 3.75 * hol
             a1 = 0.125 * rat + 0.125
             b1 = 0.125 * rat - 0.125
             c1 = 0.25 * rat
         end if
! Very stable (from extrapolation of Shaw et al's values at 0.1 and 0.5:
         ! Paul's MV3D_KT(nk) is our dkt3(km) m2 s-1 atmospheric heat diffusivity (thermal vertical diffusion coefficient)
         !  1st (bottom) model layer
         if(hol >= 0.9 .or. dkt3(km) <= min_kt) then
             a1 = 0.25
             b1 = 0.0
             c1 = 0.25
         end if
!  Raupach's originals:
!         if (zr >= 1.0) then
!            sigw = ustar(i) * 1.25
!         else
!            sigw = ustar(i) * ( 0.75 + 0.5 * cos(pi * (1.0 - (zmid3(km) - z2(km+1))/hcan) ) )
!         end if
!  Replace Raupach's originals with fit to Patton et al and Shaw et al 1988
         if(zr < 0.175) then
               sigw = ustar(i) * 0.25
         else
           if(zr < 1.25) then
               sigw = ustar(i) * ( a1 + b1 * cos(pi / 1.06818 * &
                      (1.25 - (zmid3(km) - z2(km+1)) / hcan)))
           else
               sigw = ustar(i) * c1
           end if
         end if

         tl = hcan / ustar(i)  * &
              (0.256 * ((zmid3(km) - z2(km+1) - 0.75 * hcan) / hcan) + &
               0.492 * exp (-(0.256 * ((zmid3(km) - z2(km+1)) / hcan) / 0.492)))
! ktr is the ratio of the resolved model diffusivity at the lowest resolved
! model level to that derived by Raupach's formula
!
         ktr =  dkt3(km) / (sigw * sigw * tl)
         kur =  dku3(km) / (sigw * sigw * tl)
!         print*, 'CANOPY_LEVS: KTR= ', i, ktr, dkt3(km), kk, kc
!
!  Use Raupach's formulae for diffusivity, multiplied by the above ratio, for the canopy layers:
!
         zr = (zcan3(kc) - z2(km+1)) / hcan
!  Gradient in stability under the canopy is reduced for higher stability conditions
!  in accord with Shaw, den Hartog and Neumann, BLM 45, 391-409, 1988, Fig 16.
!  Raupach's original:
!         if (zr >= 1.0) then
!            sigw = ustar(i) * 1.25
!         else
!            sigw = ustar(i) * ( 0.75 + 0.5 * cos(pi * (1.0 - (zcan3(kc) - z2(km+1))/hcan) ) )
!         end if
         if(zr < 0.175) then
               sigw = ustar(i) * 0.25
         else
           if(zr < 1.25) then
               sigw = ustar(i) * ( a1 + b1 * cos(pi / 1.06818 * &
                      (1.25 - (zcan3(kc) - z2(km+1))/hcan)))
           else
               sigw = ustar(i) * c1
           end if
         end if
!
         tl = hcan / ustar(i) *  &
              (0.256 * ( (zcan3(kc) - z2(km+1) - 0.75 * hcan) / hcan) + &
              (0.492 * exp (-(0.256 * (zcan3(kc) - z2(km+1)) / hcan) / 0.492) ) )

         dkt_can3(kk)  = (sigw * sigw * tl) * ktr
         dku_can3(kk)  = (sigw * sigw * tl) * kur

!  DKT_CAN=0.178022242775362      54.2361811640303        1.11225899578581               64           1
!  DKT_CAN=7.201550034628344E-002 47.9798060091286       0.161019598920152               66           2
!  DKT_CAN=3.982132984178101E-002 46.0438951730293       4.724674166464671E-002          67           3
!        print*, 'CANOPY_LEVS: DKT_CAN= ', i, sigw, tl, dkt_can3(kk), kk, kc
      end do ! kc = 1,nkc
!
   if (local_dbg) then
      do kc = 1, nkc
         flag_error = .false.
         if (kcan3(i, kc) == 0) then
            write(6,*) 'kcan zero inside canopy_levs at i kc = ', &
                        i, kc
            flag_error = .true.
            return
         end if
      end do
   end if
!
      do k = 1, nkt! from top to bottom of combined layers
         II = nkt + 1 - k  ! from bottom to top of combined layer

         ! Flip back meteo arrays on combined layers in same layer order as original model layer
         ! nkt  is top model layer         <= 1
         ! ...
         ! (4) is 1st (bottom) model layer   <= km
         ! (3) is 3rd (top) canopy layer     <= nkt-2
         ! (2) is 2nd canopy layer           <= nkt-1
         ! (1) is 1st (bottom) canopy layer  <=nkt
         ZH_CAN  (i,II) = zmid_can3(i, k)
         ZF_CAN  (i,II) = zmom_can3(i, k)
         PRSL_CAN(i,II) = prsl_can3(k)

         T1_CAN  (i,II) = ta_can3  (k)
         QV_CAN  (i,II) = qv_can3  (k)
         DENS_CAN(i,II) = dens_can3(k)
         WS_CAN  (i,II) = ws_can3  (k)
         DKT_CAN (i,II) = dkt_can3 (k)
         DKU_CAN (i,II) = dku_can3 (k)

! Pressure at layers centers
! 1   37.9003337896498 96.3881049029277
! 2   96.3881049029277 176.687747254452
! 3   176.687747254452 267.236282600406
! ...
! 63   99570.0993392892 100118.892141721
! 64   100118.892141721 100129.946869981
! 65   100129.946869981 100257.714673645
! 66   100257.714673645 100341.141349630
! 67   100341.141349630
!         print*,'canopy_levs: prsl_can3 =',i,k, &
!                                  prsl_can3(k),  prsl_can3(k+1)
      end do ! k = 1, nkt

! Pressure at layers interfaces
      do k = 1, nkt+1  ! from top to bottom of combined layers
         II = (nkt+1) + 1 - k  ! from bottom to top of combined layer

! Pressure at layers interfaces:
! 1    20.0000000000000
! 2    64.2470016479492
! 3   137.789993286133
! 4   221.957992553711
! ...
! 65   97574.2952071220
! 66   97892.5615950123
! 67   97999.3464530241
! 68   98097.0373946220
!
!         print*,'canopy_levs: prsi_can3 =',i,k, &
!                                  prsi_can3(k)

! (km+1) (68=nkc+km +1) prsi3( 1)   Top    model layer upper interface  prsi_can3(1)
! (km)   (67=nkc+km   ) prsi3( 2)
! ...
! (2)  (5 =nkc    +2)   prsi3(km)   Bottom model  layer upper interface  prsi_can3(km)
!      (4 =nkc    +1)               Top    canopy layer upper interface  prsi_can3(km+1)
!      (3)                          Mid    canopy layer upper interface
!      (2)                          Bottom canopy layer upper interface prsi_can3(nkt)
! (1)  (1)              prsi3(km+1) Bottom model  layer LOWER interface prsi_can3(nkt+1)
!
         PRSI_CAN(i,II) = prsi_can3(k)

      end do ! k = 1, nkt+1


   END IF ! Continuous forest canopy: FRT_MASK == 1.

! ... have not finished Paul's code ...


   END DO  !I-index

   end subroutine canopy_levs_run

   end module canopy_levs_mod

   module canopy_transfer_mod
   contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   subroutine canopy_transfer_init( im, km, nkc, nkt,   & !in
              massair_can, massair,                       & !out
              mmr_o3_can,                           & !inout
              nfrct, ifrct,                               & !out
              frctr2c, frctc2r,                           & !out
              errmsg, errflg )

!  Input/Output variables, original Horizontal coordinate
!
!  Local variables:
!  massair_can(:,nkt)  :  mass of air in canopy layers (kg)
!  massair    (:, km)  :  mass of air in model layers (kg)
!                                        (gathered canopy + resolved scale columns)
!   nfrct  (nkt,   :)   :  Number of original model levels contributing to canopy level k
!   ifrct  (nkt, 2,:)   :  Index of the original model level contributing to canopy level k
!   frctr2c(nkt, 2,:)   :  Fractional contribution of the original model level to canopy level k
!   frctc2r(nkt, 2,:)   :  Fractional contribution of the canopy level to the original model level
!
!=============================================================================

   use machine , only : kind_phys

   IMPLICIT NONE

!...Arguments:

   integer, intent(in)  :: im, km, nkc, nkt

   integer, intent(out) ::                &
                nfrct  (km+nkc,    im)  , &
                ifrct  (km+nkc, 2, im)

   real(kind=kind_phys), intent(out) ::  &
            massair_can(im, km+nkc)    , &
            massair    (im, km)        , &
            mmr_o3_can (im, km+nkc)    , &
             frctr2c   (km+nkc, 2, im) , &
             frctc2r   (km+nkc, 2, im)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!...local variables

! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   massair_can(:,:) = 0.
   massair    (:,:) = 0.
   mmr_o3_can (:,:) = 0.

   nfrct  (:,:)   = 0
   ifrct  (:,:,:) = 0
   frctr2c(:,:,:) = 1.
   frctc2r(:,:,:) = 1.

   return
   end subroutine canopy_transfer_init

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   subroutine canopy_transfer_run( im, km, nkc, nkt, &
              ntrac, ntoz,                           &
              GAREA,                                 &
              zi, zl, zm,                            &
              Q1, DENS,                              & !in: kg kg-1
              FLAG,                                  & !in
              FRT_MASK,                              & !in
              kmod, kcan3,                           & !in
              zmom_can, zmid_can,                    & !in
              PRES_CAN, DENS_CAN,                    & !in
              Q1_MOD, Q1_CAN, Q1_2M,                 & !inout kg kg-1
              massair_can, massair,                  & !inout
              mmr_o3_can,                            & !inout
              nfrct, ifrct,                          & !inout
              frctr2c, frctc2r,                      & !inout
              errmsg, errflg )

! Arguments:
! Input variables
!-----------------------------------------------------------------------------
! Array dimensions:
!  nkc                                :  number of canopy levels
!  nkt= km + nkc                :  number of levels in gathered canopy + resolved scale columns
!  met???_CAN(:.:, nkt):            :  met 3d variables, gathered canopy + resolved scale columns
!  kmod(km)                        :  Vertical index location of original ungathered model layer in combined
!                                        canopy + resolved scale column
!  flag                               : 0 -> resolved_to_canopy
!                                       1 -> canopy_to_resolved
!
!  Input/Output variables, original Horizontal coordinate
!  Q1_CAN(:,nkt, NSPCSD)     :  Chemical tracers concentrations kg kg-1 combined canopy and resolved model layers
!  Q1_MOD(:,km,  NSPCSD)     :  Chemical tracers concentrations kg kg-1 on model levels (copy of CONC)
!  Q1    (:,km,  NSPCSD)     :  Chemical tracers concentrations kg kg-1 on model levels
!  Q1_2M (:,     NSPCSD)
!
!  Local variables:
!  massair_can(:, nkt)    :  mass of air in canopy layers (kg)
!  massair    (:, km)    :  mass of air in model layers (kg)
!                                        (gathered canopy + resolved scale columns)
!   nfrct  (nkt,  :)        :  Number of original model levels contributing to canopy level k
!   ifrct  (nkt,2,:)        :  Index of the original model level contributing to canopy level k
!   frctr2c(nkt,2,:)        :  Fractional contribution of the original model level to canopy level k
!   frctc2r(nkt,2,:)        :  Fractional contribution of the canopy level to the original model level
!
!=============================================================================

   use machine , only : kind_phys

   IMPLICIT NONE

!...Arguments:

   integer, intent(in)  :: im, km, nkc, nkt, ntrac, ntoz
   integer, intent(in)  :: flag
   real(kind=kind_phys), intent(in) :: zi(im, km+1),  zl(im, km),   zm(im, km)
   real(kind=kind_phys), intent(in) :: GAREA(im)

! ** Q1 is concentration field (including gas and aerosol variables) mass mixing ratio kg kg-1
   real(kind=kind_phys), intent(in) ::   Q1(im, km, ntrac)

   real(kind=kind_phys), intent(in) :: DENS(im, km)

   integer, intent(in) :: kmod   (im, km), kcan3  (im, nkc)

   real(kind=kind_phys), intent(inout) :: zmom_can  (im, nkt+1) , &
                                          zmid_can  (im, nkt)

   real(kind=kind_phys), intent(in) ::  FRT_MASK  (im)        , &
! met3d arrays
                                        PRES_CAN  (im, nkt)   , &
                                        DENS_CAN  (im, nkt)

! all gas-phase species array
   real(kind=kind_phys), intent(inout) ::  Q1_MOD (im,  km, ntrac), &
                                           Q1_CAN (im, nkt, ntrac)
   real(kind=kind_phys), intent(inout) ::  Q1_2M  (im,      ntrac)

   integer, intent(inout) ::  nfrct  (km+nkc,    im) , &
                              ifrct  (km+nkc, 2, im)

   real(kind=kind_phys), intent(inout) ::  massair_can(im, km+nkc), &
                                           massair    (im, km)    , &
                                           mmr_o3_can (im, km+nkc), &
                                           frctr2c (km+nkc, 2, im), &
                                           frctc2r (km+nkc, 2, im)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!...Local arrays:

   real(kind=kind_phys) ::      zmid     (km)  , &
                                zmom     (km+1), & ! Same as zfull !
                                z2       (km+1), &
                                sigmid2  (km+1), &
                                zcan3    (nkc)  ,&
             pres_can3 (nkt),   pres3    (km)  , &
             dens_can3 (nkt),   dens3    (km)  , &
             klower_can(nkc)            , &
             dxdy      (im)

   real(kind=kind_phys) ::         &
              mass_canopy  (nkt), &
              mmr_canopy   (nkt), &
              vmr_canopy   (nkt), &
              mmr_resolved (km + 1), &
              vmr_resolved (km + 1), &
              mass_resolved(km), &
              conc3        (km), &
              conc_can3    (nkc)

!...local variables

   INTEGER          :: i, S, IS

   INTEGER          :: LEV, L

   INTEGER          :: KOUNT

! Diagnostic height is the assumed height above ground of the sampling for observations
   real(kind=kind_phys),    parameter              :: diag_hgt = 2.0
   real(kind=kind_phys),    parameter              :: epsilon = 1.e-10

!--------------
!hrinit.F: ...set scale factor for [ppm] -> [kg/kg]
!
! CGRID to CHEM  Species conversion factor
!         FORWARD_CONV( N ) = 1.0E-3 * MWAIR / SPECIES_MOLWT( N )  ! ug kg-1 to ppm
! CHEM  to CGRID Species conversion factor
!         REVERSE_CONV( N ) = 1.0E+3 / MWAIR * SPECIES_MOLWT( N )  ! ppm    to ug kg-1
!--------------
! Conversion factor from units in [kg kg-1] to [ug kg-1]
   REAL(kind=kind_phys), PARAMETER :: FORWARD_CONV = 1.E-9 ! ug kg-1 -> kg kg-1
   REAL(kind=kind_phys), PARAMETER :: REVERSE_CONV = 1.E+9 ! kg kg-1 -> ug kg-1

   real(kind=kind_phys) :: mmr_diag

   logical(kind=4) :: chm_error_l = .false.

   integer(kind=4) :: k, kk, kc, k2, II, npass

   logical(kind=4)                         :: local_dbg
   local_dbg = .true.

! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   conc_can3(:)=0.
   conc3    (:)=0.
   mass_canopy(:) = 0.
   mmr_canopy (:) = 0.
   vmr_canopy (:) = 0.
   mmr_resolved(:) = 0.
   vmr_resolved(:) = 0.
   mass_resolved(:) = 0.

!      RELWTEM( ICG ) = CONVMW / NR_MOLWT( SP_INDX )

   DO i = 1, im !i-index

!!! Non-Canopy columns
   IF (FRT_mask(i) <= 0.) THEN

!!!!! Start all columns!!!!! canopy & non-canopy (canopy columns are overwritten below)
      do k = 1, km       ! from bottom to top
         II = km + 1 - k  ! from top to bottom of resolved model layers km+1 ???
!!! Paul's zmom is our zmom
! zmom(1)     = ZFULL(km) is top    model  layer height
! zmom(km) = ZFULL(1)     is bottom model  layer height
         zmom (II) = zm(i,k)    ! ZFULL(i,k)
         dens3(II) = DENS(i,k)  ! kg/m**3
!! Heights of the original model layers for the canopy columns are extracted to the zmom array.

      end do

!  Calculate mass of air in model levels
      !Paul's zmom is our zmom
      zmom(km + 1) = 0.0
      do k = km, 1, -1
         ! Paul's massairmod is our massair
         massair(i, k) = dens3(k) * GAREA (i) * &
                             (zmom(k) - zmom(k + 1))
      end do

!  First, carry over original model values for the matching layers
      do k = 1, km ! from bottom to top of resolved model layers
         massair_can(i, k) = massair(i, k) ! full layer height [m]
      end do

      do kc = 1, nkc       ! from top to bottom of canopy layers
         massair_can(i, km+kc) = massair(i, km)
      end do ! kc = 1, nkc
!!!!! Non-Canopy columns !!!!!

!!!! Continuous forest canopy
   ELSE IF (FRT_mask(i) > 0.) THEN

! Put vars on combined layers in layer order as in Paul's code (GEM-MACH)
! 1      <=  nkt is top model layer
! ...
! km<= (4) is 1st (bottom) model layer
! nkt-2 <= (3) is 3rd (top) canopy layer
! nkt-1 <= (2) is 2nd canopy layer
! nkt   <= (1) is 1st (bottom) canopy layer

      do k = 1, nkt
         II = nkt + 1 - k
         pres_can3(II) = PRES_CAN(i,k)
         dens_can3(II) = DENS_CAN(i,k)  ! kg/m**3
      end do

! Calculate mass of air on combined levels
      !Paul's zmomcan is our zmom_can(im, nkt+ 1)
      ! Layers in reverse order!
      ! zmom_can(:,:,1)    is top resolved layer
      ! zmom_can(:,:,km)   is 1hy resolved layer
      ! zmom_can(:,:,nkt)  is 1st canopy layer
      zmom_can(i, nkt+ 1) = 0.0
      do k = nkt, 1, -1
         ! Paul's massaircan is our massair_can
         massair_can(i, k) = dens_can3(k) * GAREA (i) * &
                              (zmom_can(i, k) - zmom_can(i, k + 1))
      end do

      do k = 1, km       ! from bottom to top
         II = km + 1 - k  ! from top to bottom of resolved model layers km+1 ???
!!! Paul's zmom is our zmom
! zmom(1)     = ZFULL(km) is top    model  layer height
! zmom(km) = ZFULL(1)     is bottom model  layer height
         zmom (II) = zm(i,k)    ! ZFULL(i,k)
         dens3(II) = DENS(i,k)  ! kg/m**3
!! Heights of the original model layers for the canopy columns are extracted to the zmom array.
      end do

!  Calculate mass of air in model levels
      !Paul's zmom is our zmom
      zmom(km + 1) = 0.0
      do k = km, 1, -1
         ! Paul's massairmod is our massair
         massair(i, k) = dens3(k) * GAREA (i) * &
                             (zmom(k) - zmom(k + 1))
      end do

!  Next, we need a set of arrays which track mass transfer from resolved to model layers;
!  how much of the original (aka "resolved") model layer mass goes into each canopy layer,
!  given the above level structure.  The three arrays are:
!    nfrct(k,  :) : the number of resolved model levels contributing to canopy level k
!    ifrct(k,n,:) : the index of the resolved model level contributing to canopy level k (n is at most 2)
!  frctr2c(k,n,:) : the fractional contribution of the resolved model level to canopy level k
!  frctc2r(k,n,:) : the fractional contribution of the canopy model level to the resolved model level
!
!  Check for coincident layers first:
!
      inner: do k = 1, km
!  If the following IF statement is true, then the canopy and resolved
!  model layer upper and lower boundaries coincide, and the entire resolved model
!  model layer contributes to the combined model layer (trivial case).
         if (zmom_can(i, k) == zmom(k) .and. zmom_can(i, k+1) == zmom(k+1)) then
            nfrct(k,    i) = 1
            ifrct(k, 1, i) = k
            frctr2c(k, 1, i) = 1.0
            frctc2r(k, 1, i) = 1.0
         else
            exit inner
         end if
      end do inner
!
!  "k" is the first layer where boundary levels do not match on output from the above loops.
! Determine fractions of original model layer structure contributing to canopy model layers.
      k2 = k
      do k = k2, nkt
         do kk = k2, km
! (1) Upper boundaries of combined and resolved model layers coincide,
! lower boundary of combined layer is within resolved layer, so canopy
! layer resides entirely within resolved layer, and shares an upper boundary
! with the resolved layer:
            if ((zmom_can(i, k) == zmom(kk) .and. zmom_can(i, k+1) > zmom(kk+1)) .or. &
! (2) Lower boundaries coincide, upper boundary of combined layer is within resolved layer,
! so canopy layer resides entirely within the resolved layer, and shares a lower boundary
! with the canopy layer.
                (zmom_can(i, k+1) == zmom(kk+1) .and. zmom_can(i, k) > zmom(kk)) .or. &
! (3) Both canopy layer boundaries exist inside a resolved layer, with no shared boundaries:
!  nfrct(km + 1) = 1
!  ifrct(km + 1) = 64
                (zmom_can(i, k) < zmom(kk) .and. zmom_can(i, k+1) >= zmom(kk+1))) then
               nfrct(k,    i) = 1
               ifrct(k, 1, i) = kk
!              frctr2c(k, 1, i) = (zmom_can(i, k) - zmom_can(i, k+1)) / max(zmom(kk) - zmom(kk+1), epsilon)
               frctr2c(k, 1, i) = (zmom_can(i, k) - zmom_can(i, k+1)) / (zmom(kk) - zmom(kk+1))
               frctc2r(k, 1, i) = 1.0  ! canopy layer resides within resolved model layer
            end if
!  Resolved layer boundary splits a combined canopy layer:
!  This case arises if, due to the use of the momentum levels in the canopy column
!  sometimes being half-way between the thermodynamic levels, a resolved model
!  momentum layer falls within the canopy layer.  Since the resolved model layers are
!  defacto thicker than the canopy layers, this means that there can at most be two
!  resolved model layers contributing to the canopy layer (only case where nfrct = 2).
            if (zmom_can(i, k+1) < zmom(kk) .and. zmom_can(i, k) > zmom(kk)) then
               nfrct(k,    i) = 2
               ifrct(k, 1, i) = kk
               ifrct(k, 2, i) = kk-1
!  Fraction of resolved model layer contributing to canopy layer:
!              frctr2c(k, 1, i) = (zmom(kk) - zmom_can(i, k+1)) / max(zmom(kk) - zmom(kk+1), epsilon)
!              frctr2c(k, 2, i) = (zmom_can(i, k) - zmom(kk)) / max(zmom(kk-1) - zmom(kk), epsilon)
               frctr2c(k, 1, i) = (zmom(kk) - zmom_can(i, k+1)) / (zmom(kk) - zmom(kk+1))
               frctr2c(k, 2, i) = (zmom_can(i, k) - zmom(kk)) / (zmom(kk-1) - zmom(kk) )
!  Fraction of canopy layer contributing to resolved model layer:
!              frctc2r(k, 1, i) = (zmom(kk) - zmom_can(i, k+1)) / max(zmom_can(i, k) - zmom_can(i, k+1), epsilon)
!              frctc2r(k, 2, i) = (zmom_can(i, k) - zmom(kk)) / max(zmom_can(i, k) - zmom_can(i, k+1), epsilon)
               frctc2r(k, 1, i) = (zmom(kk) - zmom_can(i, k+1)) / (zmom_can(i, k) - zmom_can(i, k+1))
               frctc2r(k, 2, i) = (zmom_can(i, k) - zmom(kk)) / (zmom_can(i, k) - zmom_can(i, k+1) )
            end if
         end do
      end do

!
!  massair_can thus contains the mass of air in the canopy layers in kg, while massair contains the
!  mass of air in the original model layers, at the canopy columns (i)
!
   END IF ! Continuous forest canopy: FRT_MASK == 1.


   END DO  !i = 1, im  !I-index


!  return tracers to resolved scale model layers:

   if (flag == 1) then  ! "canopy_to_resolved"

!  At this point, the model mass is distributed over the combined layers,
!  and the tracer concentration arrays are both in the combined layer system.
!
   DO i = 1, im  !I-index

   KOUNT = 0

   !  loop over canopy columns
   IF (FRT_mask(i) > 0.) THEN

! Q1_MOD/Q1_CAN:
!   Assigned/Initilized in canopy_levs FIRSTIME

!...fetch all species in units kg kg-1 mass mixing ratio
      do S = 1, ntrac-1          ! ntrac1= 197 (ntrac=ntke=198)

! Flip resolved layer arrays into a new array for use here
         do k = 1, km        ! from bottom to top
            II = km + 1 - k  ! from top to bottom of resolved model layers
            ! conc3(1)     is top model layer
            ! conc3(km) is 1st (bottom) model layer
            ! Paul's chem_tr is our conc3 = vmr_resolved
! NB. mfpbltq_mod: q1(ix,km,ntrac1) kg kg-1
            ! Paul's chem_tr is our vmr_resolved =conc3
            vmr_resolved(II) = Q1_MOD(i, k, S)          ! kg kg-1
         end do

! Flip combined layer arrays into a new array for use here
         do k = 1, nkt          ! from top to bottom
            II = nkt + 1 - k  ! from bottom to top of resolved model layers
            ! Paul's trppm is our vmr_canopy (conc_can)
            ! (km) is top model layer
            ! (1)     is 1hy model layer
            vmr_canopy(II) = Q1_CAN(i, k, S) !kg kg-1
         end do

! (ii): Canopy shaded layers
         do kc = 1, nkc
            k = kcan3(i, kc) ! kcan3(1,2,3) = 65,66,67

            ! Paul's tracers_can is our conc_can3 array
            conc_can3 (kc) = vmr_canopy(k)    ! kg kg-1
         end do

!  (1) We start off by converting these mass mixing ratio [kg kg-1] to mass in [ug]:

         do k = 1, km
            ! kmod(1)     is 1  top model  layer
            ! kmod(km) is 65 top canopy layer (modified after mono adj.)
            kk = kmod(i, k)

! ...fetch gas mass mix. ratios [kg kg-1] and convert to [ug kg-1]
            ! Paul's conc is our mmr_canopy
            !mmr_canopy(kk) = REVERSE_CONV * conc3(k)        ! ug kg-1
            mmr_canopy(kk) = REVERSE_CONV * vmr_resolved(k)  ! ug kg-1
         end do

         do k = 1, nkc
            ! kcan3(k=1,2,3) = 65,66,67
            kc = kcan3(i, k)

! ...fetch gas mass mix. ratios [kg kg-1] and convert to [ug kg-1]
            mmr_canopy(kc) = REVERSE_CONV * conc_can3(k)  ! ug kg-1
         end do

! (2) Array "mass_canopy" now holds the mass of the tracer in each of the combined levels.
! This mass must be added back to the resolved levels:
         ! Paul's masscan is our mass_canopy
         ! Paul's mass_resolved is our mass_resolved
         mass_resolved(:) = 0.
         do k = 1, nkt

! Output diag
            if(S == 11) mmr_o3_can(i,k) = mmr_canopy(k) ! nto3=11 "canopy_to_resolved"

            mass_canopy(k) = mmr_canopy(k) * massair_can(i, k)  ! ug
            do kk = 1, nfrct(k, i)
               kc = ifrct(k, kk, i)
               mass_resolved(kc) = mass_resolved(kc) + mass_canopy(k) * frctc2r(k,kk,i)  ! ug
            end do
         end do

!
!  Check:  total mass in the column should be the same
         if (local_dbg) then
             call canopy_mass_check(mass_canopy, mass_resolved, i, flag, nkc, nkt, errmsg, errflg)
             if (errflg /= 0) return
         end if
!
!  (3) The masses in [ug] need to be converted back to [kg kg-1]
         do k = 1, km
!
            ! Paul's massairmod is our massair
            ! Paul's mass_resolved is our mass_resolved
!           mmr_resolved(k) = mass_resolved(k) / max(massair(i, k), epsilon)  ! ug kg-1
            mmr_resolved(k) = mass_resolved(k) / (massair(i, k))  ! ug kg-1

! (3a) Convert back m.m.r. [ug kg-1] to [kg kg-1]
            ! NB. This is Q1_MOD to be used in gas-phase hrdriver call on canopy columns
            ! Paul's chem_tr is our conc3 = vmr_resolved
            vmr_resolved(k)            = FORWARD_CONV *  mmr_resolved(k)    ! kg kg-1

         end do

         do k = 1, km ! from bottom to top
            II = km + 1 - k  ! from top to bottom of resolved model layers
            ! zmid(1)     = ZM(km) is top    model  layer height
            ! zmid(km) = ZM(1)     is bottom model  layer height
            ! Paul's zt (or ZPLUS) is our zmid
            zmid(II) = ZL(i,k) ! mid layer height [m]
!!! Heights of the original model layers for the canopy columns are extracted to the zmid array.
         end do

!
!  (4) Evaluate the diagnostic level concentration
!  Find the bounding layers above and below the diagnostic height:
!  kk'th layer is the layer above the inlet height
            kk = nkt
            do k = nkt, nkt-8, -1
               ! Paul's zt_can (MV3D_ZPLUS) is our zmid
               if (diag_hgt <= zmid_can(i, k-1) .and. &
                  diag_hgt > zmid_can(i, k)) then
                  kk = k - 1
               end if
            end do
!  If the diagnostic height is less than the lowest level, then use that level
!  for the concentration.
            if (kk == nkt) then
               mmr_diag =  mmr_canopy(nkt) ! ug kg-1
               vmr_resolved      (km + 1)      = FORWARD_CONV * mmr_canopy(nkt) ! kg kg-1

            else
! Diagnostic height 2m is always above the lowest model hybrid level ~42m
! The lines below never executed
               mmr_diag =  &
                        mmr_canopy(kk) +                        &
                       (mmr_canopy(kk) -  mmr_canopy(kk + 1)) / &
!                  max(zmid_can(i, kk) - zmid_can(i, kk + 1), epsilon) * &
                      (zmid_can(i, kk) - zmid_can(i, kk + 1)) * &
                          (diag_hgt    - zmid_can(i, kk + 1))        ! ug kg-1
               vmr_resolved      (km + 1)      = FORWARD_CONV * mmr_diag       ! kg kg-1

            end if

! Flip back resolved layers arrays for gas-phase integration (hrdriver)
         do k = 1, km           ! from top to bottom
            II = km + 1 - k  ! from bottom to top of resolved model layers
            ! Paul's trppm is our conc_can (vmr_canopy)
            ! (km) is top model layer
            ! (1)     is 1hy model layer
            Q1_MOD(i, II, S) = vmr_resolved(k)             ! kg kg-1
         end do

! 2M Diagnostics
         Q1_2M (i,     S) = vmr_resolved(km+1)             ! kg kg-1

      end do ! number of species loop s = 1, NUMB_MECH_SPC

! Print up to KOUNT number of canopy columns
      KOUNT = KOUNT + 1
!
   END IF ! loop over canopy columns FRT_MASK == 1.


   END DO  !I = 1, im !I-index

!  Done transferring from combined canopy + resolved scale back to resolved scale.  :)
!
! ========================================================================
   else ! if (flag == 0) then (canopy_transfer == "resolved_to_canopy") then
!
! In: Q1_MOD
!
! Out: Q1_CAN (vmr_canopy)
! NB. ! Paul's trppm (mach_gas_canopy) is our vmr_canopy
! ========================================================================
!

   DO i = 1, im  !I-index

   KOUNT = 0

   IF (FRT_mask(i) > 0.) THEN

!...fetch all species and convert to kg kg-1 mass mixing ratio
      DO S = 1, NTRAC-1  ! ntrac1= 197 (ntrac=ntke=198)
!     DO ISP = 1, 1     ! ntqv=1 ntoz=7 nto3=11

         ! S = CGRID_INDEX( ISP )

! Flip resolved layer arrays into a new array for use here
! (i): Model resolved layers
      do k = 1, km
         II = km + 1 - k  ! from top to bottom of resolved model layers km+1 ???
         ! Paul's chem_tr is our conc3 = vmr_resolved (q1_mod)
         ! conc3(1)     is top model layer
         ! conc3(km) is 1st (bottom) model layer
         ! conc3(II) = Q1 (i, k, S)           ! kg kg-1
         conc3(II) = Q1_MOD(i, k, S)        ! kg kg-1
         vmr_resolved(II) = Q1_MOD(i, k, S) ! kg kg-1
      end do

!  (1) We start off by converting these mass mixing ratio [kg kg-1]to mass in [ug]:
      do k = 1, km
! ...fetch gas mass mixing ratios [kg kg-1] and convert to [ug kg-1]
         ! Paul's conc is our mmr_resolved
         mmr_resolved(k) = REVERSE_CONV * conc3(k)      ! ug kg-1
      end do

!  (1) Convert the original model domain values in the current column to mass from mass mixing ratio:
!  mass_resolved = Mass mixing ratio * (density) / (volume of original model layer)  (ug)
      do k = 1, km
         mass_resolved(k) = mmr_resolved(k) * massair(i, k) ! ug
      end do

!  (2) Use the array fractions defined earlier to divide the resolved layer masses into the canopy layers,
!  and convert back to mixing ratios.  Note that the frctr2c fractions are vertical extent of the
!  contribution of the resolved layer into the canopy layer, hence the mass/volume can be divided up
!  this way:
!  mmr_canopy = sum of masses contributed / (density * volume of canopy model layeri)
            ! Paul's mmr_canopy is our mmr_canopy in ug kg-1
            ! Paul's masscan is our mass_canopy
            mmr_canopy(:) = 0.
            mass_canopy(:) = 0.
            do k = 1, nkt
               do kk = 1, nfrct(k, i)
                  kc = ifrct(k, kk, i)
                  mass_canopy(k) = mass_canopy(k) + mass_resolved(kc) * frctr2c(k, kk, i) ! ug
               end do
            end do

!
!  Check:  total mass in the column should be the same
            if (local_dbg) then
               call canopy_mass_check(mass_canopy, mass_resolved, i, flag, nkc, nkt, errmsg, errflg)
               if (errflg /= 0) return
            end if
!
            do k = 1, nkt
               ! Paul's massaircan is our massair_can
!              mmr_canopy(k) = mass_canopy(k) / max(massair_can(i, k), epsilon)  ! ug kg-1
               mmr_canopy(k) = mass_canopy(k) /    (massair_can(i, k))  ! ug kg-1

! Output diags
!               ! if(S == 11) mmr_o3_can(i,k) = mmr_canopy(k) ! nto3=11 "resolved_to_canopy"
!               if(S == 11) mmr_o3_can(i,k) = frctr2c(k, 1, i) ! "resolved_to_canopy"
               if(S == 11) mmr_o3_can(i,k) = frctr2c(k, 2, i)
            end do


!
!  (3) Replace the original model layer values with the corresponding canopy layer values, when
!  a canopy exists:
            do kk = 1, km
               k = kmod(i, kk)
               ! Paul's chem_tr is our conc3 = vmr_resolved (q1_mod) <================
!              conc3(kk)         = FORWARD_CONV * mmr_canopy(k)  ! kg kg-1
               vmr_resolved (kk) = FORWARD_CONV * mmr_canopy(k)  ! kg kg-1
            end do

! (i): Model resolved layers: for hrdriver (trppm from mach_gas_canopy)
            do kk = 1, km
               ! kmod(1)     is 1  top model  layer
               ! kmod(km) is 65 top canopy layer (modified after mono adj.)
               k = kmod(i, kk)

               ! Paul's trppm is our vmr_canopy (conc_can)
!              vmr_canopy(k) = conc3(kk)                     ! kg kg-1
               vmr_canopy(k) = vmr_resolved(kk)              ! kg kg-1
            end do
!
!  (4) Fill the canopy layers with the new mass mixing ratios
            do kc = 1, nkc
               k  = kcan3(i, kc)
               ! Paul's tracers_can is our conc_can3              <====================
               conc_can3(kc)  = FORWARD_CONV * mmr_canopy(k) ! kg kg-1
            end do

! (ii): Canopy shaded layers (for hrdriver) (trppm from mach_gas_canopy)
            do kc = 1, nkc
               ! Paul's trppm is our vmr_canopy (conc_can)
               ! kcan3(1) = 65
               ! kcan3(2) = 66
               ! kcan3(3) = 67
               k = kcan3(i, kc)
               vmr_canopy(k) = conc_can3(kc)                      ! kg kg-1
            end do

! Prepare array for gas-phase chemical integration. (Paul's mach_gas_canopy)
!
! Flip back augmented canopy+resolved arrays for gas-phase integration (hrdriver)
         do k = 1, nkt          ! from top to bottom
            II = nkt + 1 - k  ! from bottom to top of resolved model layers
            ! (nkt) is top model layer
            ! (4)     is 1hy model layer
            ! (1-3)   are canopy layers
            ! Paul's trppm is our vmr_canopy (conc_can)
            Q1_CAN(i, II, S) = vmr_canopy(k)                    ! kg kg-1
         end do

      end do !species index loop S (formerly isp)

! Print up to KOUNT number of canopy columns
      KOUNT = KOUNT + 1

!  loop over canopy columns
   END IF ! loop over canopy columns FRT_MASK == 1.
!
   END DO  !i = 1, im  !I-index
!
   end if  ! 1="canopy_to_resolved" 0= "resolved_to_canopy"

   return

   contains

   subroutine canopy_mass_check(mass_canopy, mass_model, i, flag, nkc, nkt, errmsg, errflg)
      implicit none
      integer(kind=4),   intent(in) :: flag, i, nkc, nkt
      real(kind=kind_phys),      intent(in) :: mass_canopy(nkt), mass_model(km)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      character(len=18) :: mode_transfer
      real(kind=kind_phys) :: masstotcan, masstotres, massrat
      real(kind=kind_phys) :: sum2can(nkt), sum2res(nkt)

      masstotcan = 0.
      masstotres = 0.
      do k = 1, nkt
         masstotcan = masstotcan + mass_canopy(k)
      end do
      do k = 1, km
         masstotres = masstotres + mass_model(k)
      end do

      if (flag == 1) then
         mode_transfer = "canopy_to_resolved"
      else
         mode_transfer = "resolved_to_canopy"
      end if

!     if (masstotres > epsilon) then
      if (masstotres > 0.0 ) then
         massrat = masstotcan / masstotres
         if (massrat > 1.001 .or. massrat < 0.999) then
            write(errmsg,fmt='(*(a,f10.4,a,f10.4))') 'Conversion of mass in ccpp_canopy_transfer not conserved ' // &
                              'during ' // mode_transfer // ' evaluation. masstotcan = ', masstotcan, &
                              ' and masstotres = ', masstotres
            errflg = 1
            return
         end if
      end if
!
!  Check on the values of the fractions:  they should sum to unity across the number
!  of original model levels!
      sum2can = 0.
      sum2res = 0.
      do k = nkt, 1, -1
         do kk = 1, nfrct(k, i)
            kc = ifrct(k, kk, i)
            sum2can(kc) = sum2can(kc) + frctr2c(k, kk, i)
            sum2res(k) = sum2res(k) + frctc2r(k, kk, i)
         end do
      end do

      do k = km , 1, -1
         if (sum2can(k) < 0.999 .or. sum2can(k) > 1.001) then
            write(errmsg,fmt='(*(a,i0,a,i0,a,f10.4))') 'layer mismatch in canopy level setup in resolved to canopy indexing: ' // &
               'column ', i, ' layer ', k, ' sum=', sum2can(k)
            errflg = 1
            return
         end if
      end do
      do k = nkt, 1, -1
         if (sum2res(k) < 0.999 .or. sum2res(k) > 1.001) then
            write(errmsg,fmt='(*(a,i0,a,i0,a,f10.4))') 'layer mismatch in canopy level setup in canopy to resolved indexing: ' // &
               'column ', i, ' layer ', k, ' sum=', sum2res(k)
            errflg = 1
            return
         end if
      end do

!
   return
   end subroutine canopy_mass_check

   end subroutine canopy_transfer_run

   end module canopy_transfer_mod

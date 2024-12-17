!>\file dcyc2t3.f
!! This file contains the CCPP-compliant dcyc2t3 codes that fits
!! radiative fluxes and heating rates from a coarse radiation
!! calculation time interval into model's more frequent time steps.

!! This module contains the CCPP-compliant dcyc2t3 codes that fits
!! radiative fluxes and heating rates from a coarse radiation
!! calculation time interval into model's more frequent time steps.
      module dcyc2t3

      implicit none

      private

      public :: dcyc2t3_run

      contains

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!    dcyc2t3 fits radiative fluxes and heating rates from a coarse      !
!    radiation calc time interval into model's more frequent time steps.!
!    solar heating rates and fluxes are scaled by the ratio of cosine   !
!    of zenith angle at the current time to the mean value used in      !
!    radiation calc.  surface downward lw flux is scaled by the ratio   !
!    of current surface air temperature (temp**4) to the corresponding  !
!    temperature saved during lw radiation calculation. upward lw flux  !
!    at the surface is computed by current ground surface temperature.  !
!    surface emissivity effect will be taken in other part of the model.!
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call dcyc2t3                                                       !
!      inputs:                                                          !
!          ( solhr,slag,sdec,cdec,sinlat,coslat,                        !
!            xlon,coszen,tsfc_lnd,tsfc_ice,tsfc_wat,                    !
!            tf,tsflw,sfcemis_lnd,sfcemis_ice,sfcemis_wat,              !
!            sfcdsw,sfcnsw,sfcdlw,sfculw,swh,swhc,hlw,hlwc,             !
!            sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   !
!            sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   !
!            im, levs, deltim, fhswr,                                   !
!            dry, icy, wet                                              !
!      input/output:                                                    !
!            dtdt,dtdtnp,                                               !
!      outputs:                                                         !
!            adjsfcdsw,adjsfcnsw,adjsfcdlw,                             !
!            adjsfculw_lnd,adjsfculw_ice,adjsfculw_wat,xmu,xcosz,       !
!            adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   !
!            adjdnnbmd,adjdnndfd,adjdnvbmd,adjdnvdfd)                   !
!                                                                       !
!                                                                       !
!                                                                       
!  inputs:                                                              
!     solhr        - real, forecast time in 24-hour form (hr)           
!     slag         - real, equation of time in radians                  
!     sdec, cdec   - real, sin and cos of the solar declination angle   
!     sinlat(im), coslat(im):                                           !
!                  - real, sin and cos of latitude                      !
!     xlon   (im)  - real, longitude in radians                         !
!     coszen (im)  - real, avg of cosz over daytime sw call interval    !
!     tsfc_lnd  (im) - real, bottom surface temperature over land (k)   !
!     tsfc_ice  (im) - real, bottom surface temperature over ice (k)    !
!     tsfc_wat  (im) - real, bottom surface temperature over ocean (k)  !
!     tf     (im)  - real, surface air (layer 1) temperature (k)        !
!     sfcemis_lnd(im) - real, surface emissivity (fraction) o. land (k) !
!     sfcemis_ice(im) - real, surface emissivity (fraction) o. ice (k)  !
!     sfcemis_wat(im) - real, surface emissivity (fraction) o. ocean (k)!
!     tsflw  (im)  - real, sfc air (layer 1) temp in k saved in lw call !
!     sfcdsw (im)  - real, total sky sfc downward sw flux ( w/m**2 )    !
!     sfcnsw (im)  - real, total sky sfc net sw into ground (w/m**2)    !
!     sfcdlw (im)  - real, total sky sfc downward lw flux ( w/m**2 )    !
!     sfculw (im)  - real, total sky sfc upward lw flux ( w/m**2 )    !
!     swh(im,levs) - real, total sky sw heating rates ( k/s )           !
!     swhc(im,levs) - real, clear sky sw heating rates ( k/s )          !
!     hlw(im,levs) - real, total sky lw heating rates ( k/s )           !
!     hlwc(im,levs) - real, clear sky lw heating rates ( k/s )          !
!     sfcnirbmu(im)- real, tot sky sfc nir-beam sw upward flux (w/m2)   !
!     sfcnirdfu(im)- real, tot sky sfc nir-diff sw upward flux (w/m2)   !
!     sfcvisbmu(im)- real, tot sky sfc uv+vis-beam sw upward flux (w/m2)!
!     sfcvisdfu(im)- real, tot sky sfc uv+vis-diff sw upward flux (w/m2)!
!     sfcnirbmd(im)- real, tot sky sfc nir-beam sw downward flux (w/m2) !
!     sfcnirdfd(im)- real, tot sky sfc nir-diff sw downward flux (w/m2) !
!     sfcvisbmd(im)- real, tot sky sfc uv+vis-beam sw dnward flux (w/m2)!
!     sfcvisdfd(im)- real, tot sky sfc uv+vis-diff sw dnward flux (w/m2)!
!     im           - integer, horizontal dimension                      !
!     levs         - integer, vertical layer dimension                  !
!     deltim       - real, physics time step in seconds                 !
!     fhswr        - real, Short wave radiation time step in seconds    !
!     fhlwr        - real, Long wave radiation time step in seconds     !
!     lslwr        - logical, flag for sw radiation calls               !
!     fluxr (im, nfxr) - real, time-accum 2d rad diag fields            !
!     topflw (im)  - topfsw_type, lw radiation fluxes at toa            !
!     sfcflw (im)  - sfcflw_type, lw radiation fluxes at sfc            !
!     lsswr        - logical, flag for lw radiation calls               !
!     coszdg       - real, Cosine(SZA), daytime                         !
!     topfsw (im)  - topfsw_type, sw radiation fluxes at TOA            !
!     sfcfsw (im)  - sfcfsw_type, sw radiation fluxes at sfc            !
!     scmpsw (im)  - cmpfsw_type, special components of sw down fluxes  !
!     raddt        - real, radiation timestep                           !
!     cldsa (im,5) - real, frac of clouds in low, mid, high, total, BL  !
!     mtopa (im,3) - integer, Vertical indices for low, middle and high !
!                    cloud bases                                        !
!     mbota (im,3) - integer, Vertical indices for low, middle and high !
!                    cloud tops                                         !
!     cldtausw (im, im+LTP) - real, approx .55mu band layer cloud       !
!                             optical depth                             !
!     cldtaulw (im, im+LTP) - real, approx 10 mu band layer cloud       ! 
!                             optical depth                             !
!     tgrs (im,levs) - real, model layer mean temperature               !
!     aerodp (im, nspc1) - real, Vertical integrated optical depth for  ! 
!                          various aerosol species                      !
!     nfxr         - integer, num variables stored in fluxr array       !
!     lm           - integer, number vertical layers for rad calc       !
!     ltp          - integer, extra top layers                          !
!     nday         - integer, daytime points dimension                  !
!     kb           - integer, vertical index diff. b/w layer and lower  !
!                    bound                                              !
!     kd           - integer, vertical index diff. b/w in/out and local ! 
!     kt           - integer, vertical index diff. b/w in/out and upper !
!                    bound
!     lssav        - logical, flag for storing radiation diagnostics    !
!     nspc1        - logical, num. species for optical dept plus total  !
!     dry          - logical, true over land                            !
!     icy          - logical, true over ice                             !
!     wet          - logical, true over water                           !
!                                                                       !
!  input/output:                                                        !
!     dtdt(im,levs)- real, model time step adjusted total radiation     !
!                          heating rates ( k/s )                        !
!     dtdtnp(im,levs)- real, heating rate adjustment for SPPT           !
!                                                                       !
!  outputs:                                                             !
!     adjsfcdsw(im)- real, time step adjusted sfc dn sw flux (w/m**2)   !
!     adjsfcnsw(im)- real, time step adj sfc net sw into ground (w/m**2)!
!     adjsfcdlw(im)- real, time step adjusted sfc dn lw flux (w/m**2)   !
!     adjsfculw_lnd(im)- real, sfc upw. lw flux at current time (w/m**2)!
!     adjsfculw_ice(im)- real, sfc upw. lw flux at current time (w/m**2)!
!     adjsfculw_wat(im)- real, sfc upw. lw flux at current time (w/m**2)!
!     adjnirbmu(im)- real, t adj sfc nir-beam sw upward flux (w/m2)     !
!     adjnirdfu(im)- real, t adj sfc nir-diff sw upward flux (w/m2)     !
!     adjvisbmu(im)- real, t adj sfc uv+vis-beam sw upward flux (w/m2)  !
!     adjvisdfu(im)- real, t adj sfc uv+vis-diff sw upward flux (w/m2)  !
!     adjnirbmd(im)- real, t adj sfc nir-beam sw downward flux (w/m2)   !
!     adjnirdfd(im)- real, t adj sfc nir-diff sw downward flux (w/m2)   !
!     adjvisbmd(im)- real, t adj sfc uv+vis-beam sw dnward flux (w/m2)  !
!     adjvisdfd(im)- real, t adj sfc uv+vis-diff sw dnward flux (w/m2)  !
!     xmu   (im)   - real, time step zenith angle adjust factor for sw  !
!     xcosz (im)   - real, cosine of zenith angle at current time step  !
!                                                                       !
!  ====================    end of description    =====================  !

!>\defgroup dcyc2t3_mod RRTMG dcyc2t3 Module
!! This module contains the CCPP-compliant dcyc2t3 codes that fits
!! radiative fluxes and heating rates from a coarse radiation
!! calculation time interval into model's more frequent time steps.
!!
!! Solar heating rates and fluxes are scaled by the ratio of cosine
!! of zenith angle at the current time to the mean value used in
!! radiation calculation. Surface downward LW flux is scaled by the
!! ratio of current surface air temperature to the corresponding
!! temperature saved during LW radiation calculation. Upward LW flux
!! at the surface is computed by current ground surface temperature.
!! Surface emissivity effect will be taken in other part of the model.
!!
!! program history:
!!-          198?  nmc mrf    - created, similar as treatment in gfdl
!!                             radiation treatment
!!-          1994  y. hou     - modified solar zenith angle calculation
!!-     nov  2004  x. wu      - add sfc sw downward flux to the variable
!!                             list for sea-ice model
!!-     mar  2008  y. hou     - add cosine of zenith angle as output for
!!                             sunshine duration time calc.
!!-     sep  2008  y. hou     - separate net sw and downward lw in slrad,
!!                changed the sign of sfc net sw to consistent with
!!                 other parts of the mdl (positive value defines from
!!                 atmos to the ground). rename output fluxes as adjusted
!!                 fluxes. other minor changes such as renaming some of
!!                 passing argument names to be consistent with calling
!!                 program.
!!-     apr  2009  y. hou     - integrated with the new parallel model
!!                 along with other modifications
!!-     mar  2011  y. hou     - minor modification including rearrange
!!                 loop orders and loop structures to improve efficiency
!!-     mar  2014  x. wu      - add sfc nir/vis bm/df to the variable
!!                             list for the coupled model input
!!-     jul  2014  s moorthi  - merge gfs and nems versions
!!-     jun  2014  y. hou     - revised to include both up and down sw
!!                 spectral component fluxes
!!-     Oct  2014  y. hous s. moorthi - add emissivity contribution to
!!                             upward longwave flux
!!-     Mar  2019  s. moorthi - modify xmu calculation in a time centered
!!                             way and add more accuracy when physics
!!                             time step is close to radiation time step
!!-     Dec  2024  l. reames - move fluxr calculations from *_post
!!                             routines to dcyc3t3_run  
!> \section arg_table_dcyc2t3_run Argument Table
!! \htmlinclude dcyc2t3_run.html
!!
!!\section dcyc2t3_general RRTMG dcyc2t3 General Algorithm
!> @{
      subroutine dcyc2t3_run                                            &
!  ---  inputs:
     &     ( solhr,slag,sdec,cdec,sinlat,coslat,                        &
     &       con_g, con_cp, con_pi, con_sbc,                            &
     &       xlon,coszen,tsfc_lnd,tsfc_ice,tsfc_wat,tf,tsflw,tsfc,      &
     &       sfcemis_lnd, sfcemis_ice, sfcemis_wat,                     &
     &       sfcdsw,sfcnsw,sfcdlw,swh,swhc,hlw,hlwc,                    &
     &       sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   &
     &       sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   &
     &       im, levs, deltim, fhswr,                                   &
     &       dry, icy, wet, damp_LW_fluxadj, lfnc_k, lfnc_p0,           &
     &       use_LW_jacobian, sfculw, use_med_flux, sfculw_med,         &
     &       fluxlwUP_jac, t_lay, p_lay, p_lev, flux2D_lwUP,            &
     &       flux2D_lwDOWN,pert_radtend,do_sppt,ca_global,tsfc_radtime, &
     &       lslwr,fluxr,fhlwr,topflw,sfcflw,lsswr,coszdg,topfsw,sfcfsw,&
     &       scmpsw,raddt,cldsa,mtopa,mbota,cldtausw,cldtaulw,tgrs,     &
     &       aerodp, nfxr, lm, ltp, nday,  kb, kd , kt, lssav, nspc1,   &
!    &       dry, icy, wet, lprnt, ipr,                                 &
!  ---  input/output:
     &       dtdt,dtdtnp,htrlw,                                         &
!  ---  outputs:
     &       adjsfcdsw,adjsfcnsw,adjsfcdlw,                             &
     &       adjsfculw_lnd,adjsfculw_ice,adjsfculw_wat,xmu,xcosz,       &
     &       adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   &
     &       adjnirbmd,adjnirdfd,adjvisbmd,adjvisdfd,                   &
     &       errmsg,errflg                                              &
     &     )
!
      use machine,         only : kind_phys
      use module_radsw_parameters,            only: topfsw_type,        &
     &                                              sfcfsw_type,        &
     &                                               cmpfsw_type
      use module_radlw_parameters,            only: topflw_type,        &
     &                                              sfcflw_type
      implicit none
!
!  ---  constant parameters:
      real(kind=kind_phys), parameter :: f_eps  = 0.0001_kind_phys,     &
     &                                   zero   = 0.0d0, one = 1.0d0,   &
     &                                   hour12 = 12.0_kind_phys,       &
     &                                   f3600  = one/3600.0_kind_phys, &
     &                                   f7200  = one/7200.0_kind_phys, &
     &                                   czlimt = 0.0001_kind_phys        ! ~ cos(89.99427)

!  ---  inputs:
      integer, intent(in) :: im, kb, kd, kt, levs, nday, nspc1

!     integer, intent(in) :: ipr
!     logical lprnt
      logical, dimension(:), intent(in) :: dry, icy, wet
      logical, intent(in) :: use_LW_jacobian, damp_LW_fluxadj,          &
     &     pert_radtend, use_med_flux, lssav
      logical, intent(in) :: do_sppt,ca_global
      real(kind=kind_phys),   intent(in) :: solhr, slag, cdec, sdec,    &
     &     deltim, fhswr, lfnc_k, lfnc_p0

      real(kind=kind_phys), dimension(:), intent(in) ::                 &
     &      sinlat, coslat, xlon, coszen, tf, tsflw, sfcdlw,            &
     &      sfcdsw, sfcnsw, sfculw, tsfc
      real(kind=kind_phys), dimension(:), intent(in), optional ::       &
     &      sfculw_med, tsfc_radtime
      real(kind=kind_phys), dimension(:), intent(in) ::                 &
     &                         tsfc_lnd, tsfc_ice, tsfc_wat,            &
     &                         sfcemis_lnd, sfcemis_ice, sfcemis_wat

      real(kind=kind_phys), dimension(:), intent(in) ::                 &
     &      sfcnirbmu, sfcnirdfu, sfcvisbmu, sfcvisdfu,                 &
     &      sfcnirbmd, sfcnirdfd, sfcvisbmd, sfcvisdfd

      real(kind=kind_phys), dimension(:,:), intent(in) :: swh, hlw,     &
     &                                     swhc, hlwc, p_lay, t_lay

      real(kind=kind_phys), dimension(:,:), intent(in) :: p_lev
      real(kind=kind_phys), dimension(:,:), intent(in), optional ::     &
     &     flux2D_lwUP, flux2D_lwDOWN, fluxlwUP_jac

      real(kind_phys),           intent(in   ) :: con_g, con_cp,        &
     &     con_pi, con_sbc

      real(kind_phys)  :: pid12
      
      logical,              intent(in) :: lsswr, lslwr
      real(kind=kind_phys), intent(in) :: raddt, fhlwr
      integer,              intent(in) :: lm, ltp, nfxr         
      real(kind=kind_phys), dimension(im),        intent(in) :: coszdg
      
      real(kind=kind_phys), dimension(im,5),      intent(in) :: cldsa
      integer,              dimension(im,3),      intent(in) :: mbota
      integer,              dimension(im,3),      intent(in) :: mtopa
   
      real(kind=kind_phys), dimension(im,lm+LTP),intent(in) :: cldtausw
      real(kind=kind_phys), dimension(im,lm+LTP),intent(in) :: cldtaulw

      real(kind=kind_phys), dimension(im,levs),  intent(in) :: tgrs
      
      type(sfcflw_type), dimension(im), intent(in) :: sfcflw
      type(sfcfsw_type), dimension(im), intent(in) :: sfcfsw
      type(cmpfsw_type), dimension(im), intent(in) :: scmpsw
      type(topflw_type), dimension(im), intent(in) :: topflw
      type(topfsw_type), dimension(im), intent(in) :: topfsw
    
      real(kind=kind_phys), dimension(im,NSPC1),  intent(in) :: aerodp

!  ---  input/output:
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dtdt
      real(kind=kind_phys), dimension(:,:), intent(inout), optional ::  &
     &      dtdtnp, htrlw
      real(kind=kind_phys), dimension(im,nfxr), intent(inout) :: fluxr
!  ---  outputs:
      real(kind=kind_phys), dimension(:), intent(out) ::                &
     &      adjsfcdsw, adjsfcnsw, adjsfcdlw, xmu, xcosz,                &
     &      adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                 &
     &      adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd

      real(kind=kind_phys), dimension(:), intent(out) ::                &
     &      adjsfculw_lnd, adjsfculw_ice, adjsfculw_wat

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      integer :: i, j, k, nstp, nstl, it, istsun(im),iSFC,iTOA,itop,ibtc
      real(kind=kind_phys) :: cns,  coszn, tem1, tem2, anginc,          &
     &                        rstl, solang, dT, tem0d
      real(kind=kind_phys), dimension(im,levs+1) :: flxlwup_adj,        &
     &     flxlwdn_adj
      real(kind=kind_phys) :: fluxlwnet_adj,fluxlwnet,dT_sfc,           &
     &fluxlwDOWN_jac,lfnc,c1
      ! Length scale for flux-adjustment scaling
      real(kind=kind_phys), parameter ::                                &
     &     L = 1.
      ! Scaling factor for downwelling LW Jacobian profile.
      real(kind=kind_phys), parameter ::                                &
     &     gamma = 0.2
!
!===> ...  begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!     Vertical ordering?
      if (p_lev(1,1) .lt.  p_lev(1, levs)) then 
         iSFC = levs + 1
         iTOA = 1
      else
         iSFC = 1
         iTOA = levs + 1
      endif

      tem1 = fhswr / deltim
      nstp = max(6, nint(tem1))
      nstl = max(1, nint(nstp/tem1))
      pid12  = con_pi / hour12
      
!  - For time averaged output quantities (including total-sky and
!    clear-sky SW and LW fluxes at TOA and surface; conventional
!    3-domain cloud amount, cloud top and base pressure, and cloud top
!    temperature; aerosols AOD, etc.), store computed results in
!    corresponding slots of array fluxr with appropriate time weights.

!  --- ...  collect the fluxr data for wrtsfc

      if (lssav) then
        !if (lsswr) then
          do i=1,im
!            fluxr(i,34) = fluxr(i,34) + fhswr*aerodp(i,1)  ! total aod at 550nm
!            fluxr(i,35) = fluxr(i,35) + fhswr*aerodp(i,2)  ! DU aod at 550nm
!            fluxr(i,36) = fluxr(i,36) + fhswr*aerodp(i,3)  ! BC aod at 550nm
!            fluxr(i,37) = fluxr(i,37) + fhswr*aerodp(i,4)  ! OC aod at 550nm
!            fluxr(i,38) = fluxr(i,38) + fhswr*aerodp(i,5)  ! SU aod at 550nm
!            fluxr(i,39) = fluxr(i,39) + fhswr*aerodp(i,6)  ! SS aod at 550nm
             fluxr(i,34) = aerodp(i,1)  ! total aod at 550nm
             fluxr(i,35) = aerodp(i,2)  ! DU aod at 550nm
             fluxr(i,36) = aerodp(i,3)  ! BC aod at 550nm
             fluxr(i,37) = aerodp(i,4)  ! OC aod at 550nm
             fluxr(i,38) = aerodp(i,5)  ! SU aod at 550nm
             fluxr(i,39) = aerodp(i,6)  ! SS aod at 550nm
          enddo
        !endif

!  ---  save lw toa and sfc fluxes
        !if (lslwr) then
          do i=1,im
!  ---  lw total-sky fluxes
            fluxr(i,1 ) = fluxr(i,1 ) + fhlwr * topflw(i)%upfxc   ! total sky top lw up
            fluxr(i,19) = fluxr(i,19) + fhlwr * sfcflw(i)%dnfxc   ! total sky sfc lw dn
            fluxr(i,20) = fluxr(i,20) + fhlwr * sfcflw(i)%upfxc   ! total sky sfc lw up
!  ---  lw clear-sky fluxes
            fluxr(i,28) = fluxr(i,28) + fhlwr * topflw(i)%upfx0   ! clear sky top lw up
            fluxr(i,30) = fluxr(i,30) + fhlwr * sfcflw(i)%dnfx0   ! clear sky sfc lw dn
            fluxr(i,33) = fluxr(i,33) + fhlwr * sfcflw(i)%upfx0   ! clear sky sfc lw up
          enddo
        !endif

!  ---  save sw toa and sfc fluxes with proper diurnal sw wgt. coszen=mean cosz over daylight
!       part of sw calling interval, while coszdg= mean cosz over entire interval
        !if (lsswr) then
          do i = 1, IM
            if (coszen(i) > 0.) then
!  ---                                  sw total-sky fluxes
!                                       -------------------
              tem0d = fhswr * coszdg(i) / coszen(i)
              fluxr(i,2 ) = fluxr(i,2)  +  topfsw(i)%upfxc * tem0d  ! total sky top sw up
              fluxr(i,3 ) = fluxr(i,3)  +  sfcfsw(i)%upfxc * tem0d  ! total sky sfc sw up
              fluxr(i,4 ) = fluxr(i,4)  +  sfcfsw(i)%dnfxc * tem0d  ! total sky sfc sw dn
!  ---                                  sw uv-b fluxes
!                                       --------------
              fluxr(i,21) = fluxr(i,21) + scmpsw(i)%uvbfc * tem0d          ! total sky uv-b sw dn
              fluxr(i,22) = fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d          ! clear sky uv-b sw dn
!  ---                                  sw toa incoming fluxes
!                                       ----------------------
              fluxr(i,23) = fluxr(i,23) + topfsw(i)%dnfxc * tem0d     ! top sw dn
!  ---                                  sw sfc flux components
!                                       ----------------------
              fluxr(i,24) = fluxr(i,24) + scmpsw(i)%visbm * tem0d          ! uv/vis beam sw dn
              fluxr(i,25) = fluxr(i,25) + scmpsw(i)%visdf * tem0d          ! uv/vis diff sw dn
              fluxr(i,26) = fluxr(i,26) + scmpsw(i)%nirbm * tem0d          ! nir beam sw dn
              fluxr(i,27) = fluxr(i,27) + scmpsw(i)%nirdf * tem0d          ! nir diff sw dn
!  ---                                  sw clear-sky fluxes
!                                       -------------------
              fluxr(i,29) = fluxr(i,29) + topfsw(i)%upfx0 * tem0d  ! clear sky top sw up
              fluxr(i,31) = fluxr(i,31) + sfcfsw(i)%upfx0 * tem0d  ! clear sky sfc sw up
              fluxr(i,32) = fluxr(i,32) + sfcfsw(i)%dnfx0 * tem0d  ! clear sky sfc sw dn
            endif
          enddo
        !endif

!  ---  save total and boundary layer clouds

        !if (lsswr .or. lslwr) then
          do i=1,im
            fluxr(i,17) = fluxr(i,17) + raddt * cldsa(i,4)
            fluxr(i,18) = fluxr(i,18) + raddt * cldsa(i,5)
          enddo

!  ---  save cld frac,toplyr,botlyr and top temp, note that the order
!       of h,m,l cloud is reversed for the fluxr output.
!  ---  save interface pressure (pa) of top/bot

          do j = 1, 3
            do i = 1, IM
              tem0d = raddt * cldsa(i,j)
              itop  = mtopa(i,j) - kd
              ibtc  = mbota(i,j) - kd
              fluxr(i, 8-j) = fluxr(i, 8-j) + tem0d
              fluxr(i,11-j) = fluxr(i,11-j) + tem0d * p_lev(i,itop+kt)
              fluxr(i,14-j) = fluxr(i,14-j) + tem0d * p_lev(i,ibtc+kb)
              fluxr(i,17-j) = fluxr(i,17-j) + tem0d * tgrs(i,itop)
            enddo
          enddo

!       Anning adds optical depth and emissivity output
          !if (lsswr .and. (nday > 0)) then
          if (nday > 0) then
            do j = 1, 3
              do i = 1, IM
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j) - kd
                ibtc  = mbota(i,j) - kd
                tem1 = 0.
                do k=ibtc,itop
                  tem1 = tem1 + cldtausw(i,k)      ! approx .55 um channel
                enddo
                fluxr(i,43-j) = fluxr(i,43-j) + tem0d * tem1
              enddo
            enddo
          endif

          !if (lslwr) then
            do j = 1, 3
              do i = 1, IM
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j) - kd
                ibtc  = mbota(i,j) - kd
                tem2 = 0.
                do k=ibtc,itop
                  tem2 = tem2 + cldtaulw(i,k)      ! approx 10. um channel
                enddo
                fluxr(i,46-j) = fluxr(i,46-j) + tem0d * (1.0-exp(-tem2))
              enddo
            enddo
          !endif

        !endif

      endif                                ! end_if_lssav
      
      
!
!  --- ...  sw time-step adjustment for current cosine of zenith angle
!           ----------------------------------------------------------
      if (nstl == 1) then
        cns = pid12 * (solhr + deltim*f7200 - hour12) + slag
        do i = 1, IM
          xcosz(i) = sdec*sinlat(i) + cdec*coslat(i)*cos(cns+xlon(i))
        enddo
      elseif (nstl == nstp) then
        do i = 1, IM
          xcosz(i) = coszen(i)
        enddo
      else
        rstl = one / float(nstl)
        solang = pid12 * (solhr - hour12)         
        anginc = pid12 * deltim * f3600 * rstl
        do i = 1, im
          xcosz(i)  = zero
          istsun(i) = zero
        enddo
        do it=1,nstl
          cns = solang + (float(it)-0.5_kind_phys)*anginc + slag
          do i = 1, IM
            coszn    = sdec*sinlat(i) + cdec*coslat(i)*cos(cns+xlon(i))
            xcosz(i) = xcosz(i) + max(zero, coszn)
            if (coszn > czlimt) istsun(i) = istsun(i) + 1
          enddo
        enddo
        do i = 1, IM
          if (istsun(i) > 0) xcosz(i) = xcosz(i) / istsun(i)  ! mean cosine of solar zenith angle at current time
        enddo
      endif
!

      do i = 1, im
!> - LW time-step adjustment:
         tem1 = tf(i) / tsflw(i)
         tem2 = tem1 * tem1
         adjsfcdlw(i) = sfcdlw(i) * tem2 * tem2
!!  - adjust \a sfc downward LW flux to account for t changes in the lowest model layer.
!! compute 4th power of the ratio of \c tf in the lowest model layer over the mean value \c tsflw.
         if (dry(i)) then
            tem2 = tsfc_lnd(i) * tsfc_lnd(i)
            adjsfculw_lnd(i) =  sfcemis_lnd(i) * con_sbc * tem2 * tem2
     &                        + (one - sfcemis_lnd(i)) * adjsfcdlw(i)
         endif
         if (icy(i)) then
            tem2 = tsfc_ice(i) * tsfc_ice(i)
            adjsfculw_ice(i) =  sfcemis_ice(i) * con_sbc * tem2 * tem2
     &                        + (one - sfcemis_ice(i)) * adjsfcdlw(i)
         endif
         if (wet(i)) then
            tem2 = tsfc_wat(i) * tsfc_wat(i)
            adjsfculw_wat(i) =  sfcemis_wat(i) * con_sbc *
     &                        tem2 * tem2
     &                        + (one - sfcemis_wat(i)) * adjsfcdlw(i)
!>  - replace upward longwave flux provided by the mediator (zero over lakes)
            if (use_med_flux) then
               if (sfculw_med(i) > f_eps) then
                  adjsfculw_wat(i) = sfculw_med(i)
               end if
            end if
         endif

!     if (lprnt .and. i == ipr) write(0,*)' in dcyc3: dry==',dry(i)
!    &,' wet=',wet(i),' icy=',icy(i),' tsfc3=',tsfc3(i,:)
!    &,' sfcemis=',sfcemis(i,:)
!

!>  - normalize by average value over radiation period for daytime.
        if ( xcosz(i) > f_eps .and. coszen(i) > f_eps ) then
          xmu(i) = xcosz(i) / coszen(i)
        else
          xmu(i) = zero
        endif

!>  - adjust \a sfc net and downward SW fluxes for zenith angle changes.
!      note: sfc emiss effect will not be appied here

        adjsfcnsw(i) = sfcnsw(i)    * xmu(i)
        adjsfcdsw(i) = sfcdsw(i)    * xmu(i)

        adjnirbmu(i) = sfcnirbmu(i) * xmu(i)
        adjnirdfu(i) = sfcnirdfu(i) * xmu(i)
        adjvisbmu(i) = sfcvisbmu(i) * xmu(i)
        adjvisdfu(i) = sfcvisdfu(i) * xmu(i)

        adjnirbmd(i) = sfcnirbmd(i) * xmu(i)
        adjnirdfd(i) = sfcnirdfd(i) * xmu(i)
        adjvisbmd(i) = sfcvisbmd(i) * xmu(i)
        adjvisdfd(i) = sfcvisdfd(i) * xmu(i)
      enddo

      ! Adjust the LW and SW heating-rates.
      ! For LW, optionally scale using the Jacobian of the upward LW flux. *RRTMGP ONLY*
      ! For SW, adjust heating rates with zenith angle change.
      if (use_LW_jacobian) then
         ! Compute adjusted net LW flux foillowing Hogan and Bozzo 2015 (10.1002/2015MS000455)
         ! Here we assume that the profile of the downwelling LW Jaconiam has the same shape
         ! as the upwelling, but scaled and offset.
         ! The scaling factor is 0.2
         ! The profile of the downwelling Jacobian (J) is offset so that
         !     J_dn_sfc / J_up_sfc = scaling_factor
         !     J_dn_toa / J_up_sfc = 0
         !
         ! Optionally, the flux adjustment can be damped with height using a logistic function
         ! fx ~ L / (1 + exp(-k*dp)), where dp = p - p0
         ! L  = 1, fix scale between 0-1.      - Fixed
         ! k  = 1 / pressure decay length (Pa) - Controlled by namelist
         ! p0 = Transition pressure (Pa)       - Controlled by namelsit
         do i = 1, im
            c1 = fluxlwUP_jac(i,iTOA) / fluxlwUP_jac(i,iSFC)
            !dT_sfc = t_lev2(i,iSFC) - t_lev(i,iSFC)
            dT_sfc = tsfc(i) - tsfc_radtime(i)
            do k = 1, levs
               ! LW net flux
               fluxlwnet = (flux2D_lwUP(i,  k+1) - flux2D_lwUP(i,  k) - &
     &                      flux2D_lwDOWN(i,k+1) + flux2D_lwDOWN(i,k))
               ! Downward LW Jacobian (Eq. 9)
               fluxlwDOWN_jac = gamma *                                 &
     &              (fluxlwUP_jac(i,k)/fluxlwUP_jac(i,iSFC) - c1) /     &
     &              (1 - c1)
               ! Adjusted LW net flux(Eq. 10)
               fluxlwnet_adj = fluxlwnet + dT_sfc*                      &
     &              (fluxlwUP_jac(i,k)/fluxlwUP_jac(i,iSFC) -           &
     &              fluxlwDOWN_jac)
               ! Adjusted LW heating rate
               htrlw(i,k) = fluxlwnet_adj * con_g /                     &
     &              (con_cp * (p_lev(i,k+1) - p_lev(i,k)))

               ! Add radiative heating rates to physics heating rate. Optionally, scaled w/ height
               ! using a logistic function
               if (damp_LW_fluxadj) then
                  lfnc = L / (1+exp(-(p_lev(i,k) - lfnc_p0)/lfnc_k))
               else
                  lfnc = 1.
               endif
               dtdt(i,k) = dtdt(i,k) + swh(i,k)*xmu(i) +                &
     &              htrlw(i,k)*lfnc + (1.-lfnc)*hlw(i,k)
            enddo
         enddo
      else
         do k = 1, levs
            do i = 1, im
               dtdt(i,k)  = dtdt(i,k)  + swh(i,k)*xmu(i)  + hlw(i,k)
            enddo
         enddo
      endif

      if (do_sppt .or. ca_global) then
         if (pert_radtend) then
! clear sky
           do k = 1, levs
             do i = 1, im
               dtdtnp(i,k) = dtdtnp(i,k) + swhc(i,k)*xmu(i) + hlwc(i,k)  
             enddo
           enddo
         else
! all sky
           do k = 1, levs
             do i = 1, im
               dtdtnp(i,k) = dtdtnp(i,k) + swh(i,k)*xmu(i) + hlw(i,k)
             enddo
           enddo
         endif
      endif
!
      return
!...................................
      end subroutine dcyc2t3_run
!> @}
!-----------------------------------
      end module dcyc2t3

!>\file dcyc2.f
!! This file contains the CCPP-compliant dcyc2t3 codes that fits
!! radiative fluxes and heating rates from a coarse radiation
!! calculation time interval into model's more frequent time steps.

!! This module contains the CCPP-compliant dcyc2t3 codes that fits
!! radiative fluxes and heating rates from a coarse radiation
!! calculation time interval into model's more frequent time steps.
      module dcyc2t3

      implicit none

      private

      public :: dcyc2t3_init, dcyc2t3_run, dcyc2t3_finalize

      contains

      subroutine dcyc2t3_init()
      end subroutine dcyc2t3_init

      subroutine dcyc2t3_finalize()
      end subroutine dcyc2t3_finalize

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
!            adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,                   !
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
!> \section arg_table_dcyc2t3_run Argument Table
!! \htmlinclude dcyc2t3_run.html
!!
!!\section dcyc2t3_general RRTMG dcyc2t3 General Algorithm
!> @{
      subroutine dcyc2t3_run                                            &
!  ---  inputs:
     &     ( solhr,slag,sdec,cdec,sinlat,coslat,                        &
     &       xlon,coszen,tsfc_lnd,tsfc_ice,tsfc_wat,tf,tsflw,           &
     &       sfcemis_lnd, sfcemis_ice, sfcemis_wat,                     &
     &       sfcdsw,sfcnsw,sfcdlw,swh,swhc,hlw,hlwc,                    &
     &       sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   &
     &       sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   &
     &       im, levs, deltim, fhswr,                                   &
     &       dry, icy, wet,                                             &
     &       use_LW_jacobian, sfculw, sfculw_jac,                       &
     &       pert_radtend, do_sppt,ca_global,                           &
!    &       dry, icy, wet, lprnt, ipr,                                 &
!  ---  input/output:
     &       dtdt,dtdtnp,                                               &
!  ---  outputs:
     &       adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,                   &
     &       adjsfculw_lnd,adjsfculw_ice,adjsfculw_wat,xmu,xcosz,       &
     &       adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   &
     &       adjnirbmd,adjnirdfd,adjvisbmd,adjvisdfd,                   &
     &       errmsg,errflg                                              &
     &     )
!
      use machine,         only : kind_phys
      use physcons,        only : con_pi, con_sbc

      implicit none
!
!  ---  constant parameters:
      real(kind=kind_phys), parameter :: f_eps  = 0.0001_kind_phys,     &
     &                                   zero   = 0.0d0, one = 1.0d0,   &
     &                                   hour12 = 12.0_kind_phys,       &
     &                                   f3600  = one/3600.0_kind_phys, &
     &                                   f7200  = one/7200.0_kind_phys, &
     &                                   czlimt = 0.0001_kind_phys,     &    ! ~ cos(89.99427)
     &                                   pid12  = con_pi / hour12

!  ---  inputs:
      integer, intent(in) :: im, levs

!     integer, intent(in) :: ipr
!     logical lprnt
      logical, dimension(im), intent(in) :: dry, icy, wet
      logical, intent(in) :: use_LW_jacobian, pert_radtend
      logical, intent(in) :: do_sppt,ca_global
      real(kind=kind_phys),   intent(in) :: solhr, slag, cdec, sdec,    &
     &                                      deltim, fhswr

      real(kind=kind_phys), dimension(im), intent(in) ::                &
     &      sinlat, coslat, xlon, coszen, tf, tsflw, sfcdlw,            &
     &      sfcdsw, sfcnsw, sfculw, sfculw_jac

      real(kind=kind_phys), dimension(im), intent(in) ::                &
     &                         tsfc_lnd, tsfc_ice, tsfc_wat,            &
     &                         sfcemis_lnd, sfcemis_ice, sfcemis_wat

      real(kind=kind_phys), dimension(im), intent(in) ::                &
     &      sfcnirbmu, sfcnirdfu, sfcvisbmu, sfcvisdfu,                 &
     &      sfcnirbmd, sfcnirdfd, sfcvisbmd, sfcvisdfd

      real(kind=kind_phys), dimension(im,levs), intent(in) :: swh,  hlw &
     &,                                                       swhc, hlwc

!  ---  input/output:
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: dtdt 
      real(kind=kind_phys), dimension(:,:),     intent(inout) :: dtdtnp

!  ---  outputs:
      real(kind=kind_phys), dimension(im), intent(out) ::               &
     &      adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw, xmu, xcosz,     &
     &      adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                 &
     &      adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd

      real(kind=kind_phys), dimension(im), intent(out) ::               &
     &      adjsfculw_lnd, adjsfculw_ice, adjsfculw_wat

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      integer :: i, k, nstp, nstl, it, istsun(im)
      real(kind=kind_phys) :: cns,  coszn, tem1, tem2, anginc,          &
     &                        rstl, solang, dT
!
!===> ...  begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      tem1 = fhswr / deltim
      nstp = max(6, nint(tem1))
      nstl = max(1, nint(nstp/tem1))
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
         tem1 = tf(i) / tsflw(i)
         tem2 = tem1 * tem1
         adjsfcdlw(i) = sfcdlw(i) * tem2 * tem2
!> - LW time-step adjustment:
         if (use_LW_Jacobian) then
            ! F_adj = F_o + (dF/dT) * dT	
            dT           = tf(i) - tsflw(i)
            adjsfculw(i) = sfculw(i) + sfculw_jac(i) * dT
         else
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
             adjsfculw_wat(i) =  sfcemis_wat(i) * con_sbc * tem2 * tem2
     &                        + (one - sfcemis_wat(i)) * adjsfcdlw(i)
          endif
        endif  
!     if (lprnt .and. i == ipr) write(0,*)' in dcyc3: dry==',dry(i)
!    &,' wet=',wet(i),' icy=',icy(i),' tsfc3=',tsfc3(i,:)
!    &,' sfcemis=',sfcemis(i,:),' adjsfculw=',adjsfculw(i,:)
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

!>  - adjust SW heating rates with zenith angle change and
!! add with LW heating to temperature tendency.

      do k = 1, levs
        do i = 1, im
          dtdt(i,k)  = dtdt(i,k)  + swh(i,k)*xmu(i)  + hlw(i,k)
        enddo
      enddo
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

      module dcyc2t3

      contains

!! \section arg_table_dcyc2t3_init Argument Table
!!
      subroutine dcyc2t3_init()
      end subroutine dcyc2t3_init


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
!            xlon,coszen,tsea,tf,tsflw,sfcemis,                         !
!            sfcdsw,sfcnsw,sfcdlw,swh,swhc,hlw,hlwc,                    !
!            sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   !
!            sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   !
!            ix, im, levs,                                              !
!      input/output:                                                    !
!            dtdt,dtdtc,                                                !
!      outputs:                                                         !
!            adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,xmu,xcosz,         !
!            adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   !
!            adjdnnbmd,adjdnndfd,adjdnvbmd,adjdnvdfd)                   !
!                                                                       !
!                                                                       !
!  program history:                                                     !
!          198?  nmc mrf    - created, similar as treatment in gfdl     !
!                             radiation treatment                       !
!          1994  y. hou     - modified solar zenith angle calculation   !
!     nov  2004  x. wu      - add sfc sw downward flux to the variable  !
!                             list for sea-ice model                    !
!     mar  2008  y. hou     - add cosine of zenith angle as output for  !
!                             sunshine duration time calc.              !
!     sep  2008  y. hou     - separate net sw and downward lw in slrad, !
!                 changed the sign of sfc net sw to consistent with     !
!                 other parts of the mdl (positive value defines from   !
!                 atmos to the ground). rename output fluxes as adjusted!
!                 fluxes. other minor changes such as renaming some of  !
!                 passing argument names to be consistent with calling  !
!                 program.                                              !
!     apr  2009  y. hou     - integrated with the new parallel model    !
!                 along with other modifications                        !
!     mar  2011  y. hou     - minor modification including rearrange    !
!                 loop orders and loop structures to improve efficiency !
!     mar  2014  x. wu      - add sfc nir/vis bm/df to the variable     !
!                             list for the coupled model input          !
!     jul  2014  s moorthi  - merge gfs and nems versions               !
!     jun  2014  y. hou     - revised to include both up and down sw    !
!                 spectral component fluxes                             !
!     Oct  2014  y. hous s. moorthi - add emissivity contribution to    !
!                             upward longwave flux                      !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     solhr        - real, forecast time in 24-hour form (hr)           !
!     slag         - real, equation of time in radians                  !
!     sdec, cdec   - real, sin and cos of the solar declination angle   !
!     sinlat(im), coslat(im):                                           !
!                  - real, sin and cos of latitude                      !
!     xlon   (im)  - real, longitude in radians                         !
!     coszen (im)  - real, avg of cosz over daytime sw call interval    !
!     tsea   (im)  - real, ground surface temperature (k)               !
!     tf     (im)  - real, surface air (layer 1) temperature (k)        !
!     sfcemis(im)  - real, surface emissivity (fraction)                !
!     tsflw  (im)  - real, sfc air (layer 1) temp in k saved in lw call !
!     sfcdsw (im)  - real, total sky sfc downward sw flux ( w/m**2 )    !
!     sfcnsw (im)  - real, total sky sfc net sw into ground (w/m**2)    !
!     sfcdlw (im)  - real, total sky sfc downward lw flux ( w/m**2 )    !
!     swh(ix,levs) - real, total sky sw heating rates ( k/s )           !
!     swhc(ix,levs) - real, clear sky sw heating rates ( k/s )          !
!     hlw(ix,levs) - real, total sky lw heating rates ( k/s )           !
!     hlwc(ix,levs) - real, clear sky lw heating rates ( k/s )          !
!     sfcnirbmu(im)- real, tot sky sfc nir-beam sw upward flux (w/m2)   !
!     sfcnirdfu(im)- real, tot sky sfc nir-diff sw upward flux (w/m2)   !
!     sfcvisbmu(im)- real, tot sky sfc uv+vis-beam sw upward flux (w/m2)!
!     sfcvisdfu(im)- real, tot sky sfc uv+vis-diff sw upward flux (w/m2)!
!     sfcnirbmd(im)- real, tot sky sfc nir-beam sw downward flux (w/m2) !
!     sfcnirdfd(im)- real, tot sky sfc nir-diff sw downward flux (w/m2) !
!     sfcvisbmd(im)- real, tot sky sfc uv+vis-beam sw dnward flux (w/m2)!
!     sfcvisdfd(im)- real, tot sky sfc uv+vis-diff sw dnward flux (w/m2)!
!     ix, im       - integer, horiz. dimention and num of used points   !
!     levs         - integer, vertical layer dimension                  !
!                                                                       !
!  input/output:                                                        !
!     dtdt(im,levs)- real, model time step adjusted total radiation     !
!                          heating rates ( k/s )                        !
!     dtdtc(im,levs)- real, model time step adjusted clear sky radiation!
!                          heating rates ( k/s )                        !
!                                                                       !
!  outputs:                                                             !
!     adjsfcdsw(im)- real, time step adjusted sfc dn sw flux (w/m**2)   !
!     adjsfcnsw(im)- real, time step adj sfc net sw into ground (w/m**2)!
!     adjsfcdlw(im)- real, time step adjusted sfc dn lw flux (w/m**2)   !
!     adjsfculw(im)- real, sfc upward lw flux at current time (w/m**2)  !
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

!-----------------------------------
!! \section arg_table_dcyc2t3_run Argument Table
!! | local var name | longname                                                                                 | description                                                                         | units   | rank | type    | kind      | intent | optional |
!! |----------------|------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|---------|------|---------|-----------|--------|----------|
!! | solhr          | forecast_hour                                                                            | forecast time in 24-hour form                                                       | hr      | 0    | real    | kind_phys | in     | F        |
!! | slag           | equation_of_time                                                                         | equation of time                                                                    | radians | 0    | real    | kind_phys | in     | F        |
!! | sdec           | sine_of_solar_declination_angle                                                          | sine of solar declination angle                                                     | none    | 0    | real    | kind_phys | in     | F        |
!! | cdec           | cosine_of_solar_declination_angle                                                        | cosine of solar declination angle                                                   | none    | 0    | real    | kind_phys | in     | F        |
!! | sinlat         | sine_of_latitude                                                                         | sine of latitude                                                                    | none    | 1    | real    | kind_phys | in     | F        |
!! | coslat         | cosine_of_latitude                                                                       | cosine of latitude                                                                  | none    | 1    | real    | kind_phys | in     | F        |
!! | xlon           | longitude_of_grid_point                                                                  | longitude of grid point                                                             | radians | 1    | real    | kind_phys | in     | F        |
!! | coszen         | cosine_of_zenith_angle                                                                   | average of cosine of zenith angle over daytime shortwave call time interval         | none    | 1    | real    | kind_phys | in     | F        |
!! | tsea           | ground_surface_temperature                                                               | ground surface temperature                                                          | K       | 1    | real    | kind_phys | in     | F        |
!! | tf             | surface_midlayer_air_temperature                                                         | surface (first layer) air temperature                                               | K       | 1    | real    | kind_phys | in     | F        |
!! | tsflw          | surface_midlayer_air_temperature_in_longwave_radiation                                   | surface (first layer) air temperature saved in longwave radiation call              | K       | 1    | real    | kind_phys | in     | F        |
!! | sfcemis        | surface_longwave_emissivity                                                              | surface emissivity                                                                  | frac    | 1    | real    | kind_phys | in     | F        |
!! | sfcdsw         | downwelling_shortwave_flux_at_surface                                                    | total sky surface downward shortwave flux                                           | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcnsw         | net_shortwave_flux_into_ground_at_surface                                                | total sky surface net shortwave flux into ground                                    | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcdlw         | downwelling_longwave_flux_at_surface                                                     | total sky surface downward longwave flux                                            | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | swh            | tendency_of_air_temperature_due_to_shortwave_heating                                     | total sky shortwave heating rate                                                    | K s-1   | 2    | real    | kind_phys | in     | F        |
!! | swhc           | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky                  | clear sky shortwave heating rate                                                    | K s-1   | 2    | real    | kind_phys | in     | F        |
!! | hlw            | tendency_of_air_temperature_due_to_longwave_heating                                      | total sky longwave heating rate                                                     | K s-1   | 2    | real    | kind_phys | in     | F        |
!! | hlwc           | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky                   | clear sky longwave heating rates                                                    | K s-1   | 2    | real    | kind_phys | in     | F        |
!! | sfcnirbmu      | direct_upwelling_near_infrared_shortwave_flux_at_surface                                 | total sky surface near-infrared beam shortwave upward flux                          | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcnirdfu      | diffuse_upwelling_near_infrared_shortwave_flux_at_surface                                | total sky surface near-infrared diffuse shortwave upward flux                       | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcvisbmu      | direct_upwelling_ultraviolet_and_visible_shortwave_flux_at_surface                       | total sky surface ultraviolet plus visible beam shortwave upward flux               | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcvisdfu      | diffuse_upwelling_ultraviolet_and_visible_shortwave_flux_at_surface                      | total sky surface ultraviolet plus visible diffuse shortwave upward flux            | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcnirbmd      | direct_downwelling_near_infrared_shortwave_flux_at_surface                               | total sky surface near-infrared beam shortwave downward flux                        | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcnirdfd      | diffuse_downwelling_near_infrared_shortwave_flux_at_surface                              | total sky surface near-infrared diffuse shortwave downward flux                     | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcvisbmd      | direct_downwelling_ultraviolet_and_visible_shortwave_flux_at_surface                     | total sky surface ultraviolet plus visible beam shortwave downward flux             | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | sfcvisdfd      | diffuse_downwelling_ultraviolet_and_visible_shortwave_flux_at_surface                    | total sky surface ultraviolet plus visible diffuse shortwave downward flux          | W m-2   | 1    | real    | kind_phys | in     | F        |
!! | ix             | horizontal_dimension                                                                     | horizontal dimension                                                                | index   | 0    | integer | default   | in     | F        |
!! | im             | horizontal_loop_extent                                                                   | horizontal loop extent                                                              | index   | 0    | integer | default   | in     | F        |
!! | levs           | vertical_dimension                                                                       | number of vertical layers                                                           | index   | 0    | integer | default   | in     | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics                                         | model time-step-adjusted total radiation heating rate                               | K s-1   | 2    | real    | kind_phys | inout  | F        |
!! | dtdtc          | tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky                  | model time-step-adjusted clear sky radiative (shortwave+longwave) heating rate      | K s-1   | 2    | real    | kind_phys | inout  | F        |
!! | adjsfcdsw      | time_step_adjusted_downwelling_shortwave_flux_at_surface                                 | time-step-adjusted surface downward shortwave flux                                  | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjsfcnsw      | time_step_adjusted_net_shortwave_flux_into_ground_at_surface                             | time-step-adjusted surface net shortwave flux into ground                           | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjsfcdlw      | time_step_adjusted_downwelling_longwave_flux_at_surface                                  | time-step-adjusted surface downward longwave flux                                   | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjsfculw      | upwelling_longwave_flux_at_surface                                                       | surface upward longwave flux at current time                                        | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | xmu            | time_step_zenith_angle_adjust_factor_for_sw                                              | time-step zenith angle adjustment factor for shortwave                              | none    | 1    | real    | kind_phys | out    | F        |
!! | xcosz          | cosine_of_zenith_angle_at_current_time                                                   | cosine of zenith angle at current time step                                         | none    | 1    | real    | kind_phys | out    | F        |
!! | adjnirbmu      | time_step_adjusted_direct_upwelling_near_infrared_shortwave_flux_at_surface              | time-step-adjusted surface near-infrared beam shortwave upward flux                 | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjnirdfu      | time_step_adjusted_diffuse_upwelling_near_infrared_shortwave_flux_at_surface             | time-step-adjusted surface near-infrared diffuse shortwave upward flux              | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjvisbmu      | time_step_adjusted_direct_upwelling_ultraviolet_and_visible_shortwave_flux_at_surface    | time-step-adjusted surface ultraviolet plus visible beam shortwave upward flux      | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjvisdfu      | time_step_adjusted_diffuse_upwelling_ultraviolet_and_visible_shortwave_flux_at_surface   | time-step-adjusted surface ultraviolet plus visible diffuse shortwave upward flux   | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjnirbmd      | time_step_adjusted_direct_downwelling_near_infrared_shortwave_flux_at_surface            | time-step-adjusted surface near-infrared beam shortwave downward flux               | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjnirdfd      | time_step_adjusted_diffuse_downwelling_near_infrared_shortwave_flux_at_surface           | time-step-adjusted surface near-infrared diffuse shortwave downward flux            | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjvisbmd      | time_step_adjusted_direct_downwelling_ultraviolet_and_visible_shortwave_flux_at_surface  | time-step-adjusted surface ultraviolet plus visible beam shortwave downward flux    | W m-2   | 1    | real    | kind_phys | out    | F        |
!! | adjvisdfd      | time_step_adjusted_diffuse_downwelling_ultraviolet_and_visible_shortwave_flux_at_surface | time-step-adjusted surface ultraviolet plus visible diffuse shortwave downward flux | W m-2   | 1    | real    | kind_phys | out    | F        |
!!
      subroutine dcyc2t3_run                                            &
!...................................
!  ---  inputs:
     &     ( solhr,slag,sdec,cdec,sinlat,coslat,                        &
     &       xlon,coszen,tsea,tf,tsflw,sfcemis,                         &
     &       sfcdsw,sfcnsw,sfcdlw,swh,swhc,hlw,hlwc,                    &
     &       sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   &
     &       sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   &
     &       ix, im, levs,                                              &
!  ---  input/output:
     &       dtdt,dtdtc,                                                &
!  ---  outputs:
     &       adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,xmu,xcosz,         &
     &       adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   &
     &       adjnirbmd,adjnirdfd,adjvisbmd,adjvisdfd                    &
     &     )
!
      use machine,         only : kind_phys
      use physcons,        only : con_pi, con_sbc

      implicit none
!
!  ---  constant parameters:
      real(kind=kind_phys), parameter :: f_eps  = 0.0001, hour12 = 12.0

!  ---  inputs:
      integer, intent(in) :: ix, im, levs

      real(kind=kind_phys), intent(in) :: solhr, slag, cdec, sdec

      real(kind=kind_phys), dimension(im), intent(in) ::                &
     &      sinlat, coslat, xlon, coszen, tsea, tf, tsflw, sfcdlw,      &
     &      sfcdsw, sfcnsw, sfcemis
      real(kind=kind_phys), dimension(im), intent(in) ::                &
     &      sfcnirbmu, sfcnirdfu, sfcvisbmu, sfcvisdfu,                 &
     &      sfcnirbmd, sfcnirdfd, sfcvisbmd, sfcvisdfd

      real(kind=kind_phys), dimension(ix,levs), intent(in) :: swh,  hlw
     &,                                                       swhc, hlwc&

!  ---  input/output:
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: dtdt   &
     &,                                                          dtdtc

!  ---  outputs:
      real(kind=kind_phys), dimension(im), intent(out) ::               &
     &      adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw, xmu, xcosz,     &
     &      adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                 &
     &      adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd

!  ---  locals:
      integer :: i, k
      real(kind=kind_phys) :: cns, ss, cc, ch, tem1, tem2
!
!===> ...  begin here
!
      cns = con_pi * (solhr - hour12) / hour12 + slag
!
      do i = 1, im

!  --- ...  lw time-step adjustment
!           -----------------------
!  --- ...  adjust sfc downward lw flux to account for t changes in layer 1
!           compute 4th power of the ratio of layer 1 tf over the mean value tsflw

        tem1 = tf(i) / tsflw(i)
        tem2  = tem1 * tem1
        adjsfcdlw(i) = sfcdlw(i) * tem2 * tem2

!  --- ...  compute sfc upward lw flux from current sfc temp,
!      note: sfc emiss effect is not appied here, and will be dealt in other place

        tem2 = tsea(i) * tsea(i)
        adjsfculw(i) =  sfcemis(i) * con_sbc * tem2 * tem2
     &               + (1.0 - sfcemis(i)) * adjsfcdlw(i)
!
!  --- ...  sw time-step adjustment
!           -----------------------

        ss     = sinlat(i) * sdec
        cc     = coslat(i) * cdec
        ch     = cc * cos( xlon(i)+cns )
        xcosz(i) = ch + ss        ! cosine of solar zenith angle at current time

!  --- ...  normalize by average value over radiation period for daytime.

        if ( xcosz(i) > f_eps .and. coszen(i) > f_eps ) then
          xmu(i) = xcosz(i) / coszen(i)
        else
          xmu(i) = 0.0
        endif

!  --- ...  adjust sfc net and downward sw fluxes for zenith angle changes
!      note: sfc emiss effect will not be appied here

        adjsfcnsw(i) = sfcnsw(i)   * xmu(i)
        adjsfcdsw(i) = sfcdsw(i)   * xmu(i)

        adjnirbmu(i)  = sfcnirbmu(i) * xmu(i)
        adjnirdfu(i)  = sfcnirdfu(i) * xmu(i)
        adjvisbmu(i)  = sfcvisbmu(i) * xmu(i)
        adjvisdfu(i)  = sfcvisdfu(i) * xmu(i)

        adjnirbmd(i)  = sfcnirbmd(i) * xmu(i)
        adjnirdfd(i)  = sfcnirdfd(i) * xmu(i)
        adjvisbmd(i)  = sfcvisbmd(i) * xmu(i)
        adjvisdfd(i)  = sfcvisdfd(i) * xmu(i)
      enddo

!  --- ...  adjust sw heating rates with zenith angle change and
!           add with lw heating to temperature tendency

      do k = 1, levs
        do i = 1, im
          dtdt(i,k)  = dtdt(i,k)  + swh(i,k)*xmu(i)  + hlw(i,k)
          dtdtc(i,k) = dtdtc(i,k) + swhc(i,k)*xmu(i) + hlwc(i,k)
        enddo
      enddo
!
      return
!...................................
      end subroutine dcyc2t3_run
!-----------------------------------


!! \section arg_table_dcyc2t3_finalize Argument Table
!!
      subroutine dcyc2t3_finalize()
      end subroutine dcyc2t3_finalize

      end module dcyc2t3

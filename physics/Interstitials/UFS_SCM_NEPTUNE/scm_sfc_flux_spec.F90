!> \file scm_sfc_flux_spec.F90
!!  Contains code to calculate parameters needed by the rest of the GFS physics suite given specified surface fluxes.

module scm_sfc_flux_spec

  implicit none

  private

!----------------
! Public entities
!----------------
  public  scm_sfc_flux_spec_init, scm_sfc_flux_spec_run

  CONTAINS
!*******************************************************************************************

!!
!! \section arg_table_scm_sfc_flux_spec_init Argument Table
!! \htmlinclude scm_sfc_flux_spec_init.html
!!
  subroutine scm_sfc_flux_spec_init(lheatstrg, errmsg, errflg)
    
    logical, intent(in) :: lheatstrg
    
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    
    if (lheatstrg) then
      errmsg = 'Using specified surface fluxes is not compatible with canopy heat storage (lheatstrg) being true. Stopping.'
      errflg = 1
      return
    end if
  end subroutine scm_sfc_flux_spec_init

!> \brief This routine calculates surface-related parameters given specified sensible and latent heat fluxes and a roughness length. Most of the calculation
!! is "backing out" parameters that are calculated in sfc_dff.f from the known surface heat fluxes and roughness length.
!!
!! \section arg_table_scm_sfc_flux_spec_run Argument Table
!! \htmlinclude scm_sfc_flux_spec_run.html
!!
!! \section general_sfc_flux_spec General Algorithm
!!  -# Compute friction velocity from the wind speed at the lowest model layer, the height about the ground, and the roughness length.
!!  -# Compute the surface stress from the friction velocity.
!!  -# Calculate the surface drag coefficient for momentum given the surface stress and wind on the lowest model layer.
!!  -# Calculate the Monin-Obukhov similarity funciton for momentum from the surface drag coefficient.
!!  -# Calculate the Obukhov length from the friction velocity, surface virtual potential temperature, and surface vertical virtual potential temperature flux.
!!  -# Calculate the bulk Richardson number at the lowest model layer.
!!  -# Calculate the Monin-Obukhov similarity function for heat and moisture from the bulk Richardson number and diagnosed similarity function for momentum.
!!  -# Calculate the surface drag coefficient for heat and moisture.
!!  -# Calculate the u and v wind at 10m.
  subroutine scm_sfc_flux_spec_run (im, u1, v1, z1, t1, q1, p1, roughness_length, spec_sh_flux, spec_lh_flux, &
    exner_inverse, T_surf, cp, grav, hvap, rd, fvirt, vonKarman, tgice, islmsk, dry, frland, cice, icy, tisfc,&
    oceanfrac, min_seaice, cplflx, cplice, flag_cice, wet, min_lakeice, tsfcl, tsfc_wat, slmsk, lakefrac, lkm,&
    lakedepth, use_lake_model, sh_flux, lh_flux, sh_flux_chs, u_star, sfc_stress, cm, ch, &
    fm, fh, rb, u10m, v10m, wind1, qss, t2m, q2m, errmsg, errflg)

    use machine,             only: kind_phys
    
    integer, intent(in)    :: im, lkm
    integer, intent(inout) :: islmsk(:)
    logical, intent(in)    :: cplflx, cplice
    logical, intent(inout) :: dry(:), icy(:), flag_cice(:), wet(:), use_lake_model(:)
    real(kind=kind_phys), intent(in)    :: cp, grav, hvap, rd, fvirt, vonKarman, min_seaice, tgice, min_lakeice
    real(kind=kind_phys), intent(in)    :: u1(:), v1(:), z1(:), t1(:), q1(:), p1(:), roughness_length(:), &
      spec_sh_flux(:), spec_lh_flux(:), exner_inverse(:), T_surf(:), oceanfrac(:), lakefrac(:), lakedepth(:)
    real(kind=kind_phys), intent(inout) :: cice(:), tisfc(:), tsfcl(:), tsfc_wat(:), slmsk(:)
    real(kind=kind_phys), intent(out)   :: sh_flux(:), lh_flux(:), u_star(:), sfc_stress(:), &
      cm(:), ch(:), fm(:), fh(:), rb(:), u10m(:), v10m(:), wind1(:), qss(:), t2m(:), q2m(:), &
      sh_flux_chs(:), frland(:)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i

    real(kind=kind_phys) :: rho, q1_non_neg, w_thv1, rho_cp_inverse, rho_hvap_inverse, Obukhov_length, thv1, tvs, &
      dtv, adtv, wind10m, u_fraction, roughness_length_m
    
    real(kind=kind_phys), parameter :: timin = 173.0_kind_phys  ! minimum temperature allowed for snow/ice

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
  
!     !--- set control properties (including namelist read)
  !calculate u_star from wind profiles (need roughness length, and wind and height at lowest model level)
  do i=1, im
    sh_flux(i) = spec_sh_flux(i)
    lh_flux(i) = spec_lh_flux(i)
    sh_flux_chs(i) = sh_flux(i)
    
    roughness_length_m = 0.01*roughness_length(i)

    wind1(i) = sqrt(u1(i)*u1(i) + v1(i)*v1(i))
    u_star(i) = vonKarman*wind1(i)/(log(z1(i)/roughness_length_m))

    !calculate variables related to u_star
    sfc_stress(i) = u_star(i)*u_star(i)

    if(wind1(i) > 0.0) then
      cm(i) = sfc_stress(i)/(wind1(i)*wind1(i))
    else
      cm(i) = 0.0
    end if
    fm(i) = sqrt((vonKarman*vonKarman)/cm(i))

    !calculate the Obukhov length from the specified surface fluxes
    q1_non_neg       = max( q1(i), 1.0e-8 )
    rho      = p1(i) / (rd*t1(i)*(1.0 + fvirt*q1_non_neg))
    rho_cp_inverse = 1.0/(rho*cp)
    rho_hvap_inverse = 1.0/(rho*hvap)

    thv1    = t1(i) * exner_inverse(i) * (1.0 + fvirt*q1_non_neg)
    ! sh_flux = rho_cp_inverse*sh_flux_wm2(i)
    ! lh_flux = rho_hvap_inverse*lh_flux_wm2(i)
    w_thv1 = sh_flux(i)*exner_inverse(i) + fvirt*t1(i)*lh_flux(i)

    Obukhov_length = -u_star(i)*u_star(i)*u_star(i)*thv1/(vonKarman*grav*w_thv1)

    !calculate the bulk Richardson number and the M-O function for scalars
    tvs     = T_surf(i)*(1.0 + fvirt*q1_non_neg)

    dtv     = thv1 - tvs
    adtv    = max(abs(dtv),0.001)
    dtv     = sign(1._kind_phys,dtv) * adtv
    rb(i)   = max(-5000.0, (grav+grav) * dtv * z1(i) / ((thv1 + tvs) * wind1(i) * wind1(i)))

    fh(i) = rb(i)*fm(i)*fm(i)*Obukhov_length/z1(i)
    ch(i) = vonKarman*vonKarman/(fm(i)*fh(i))

    !calculate sfc_diag variables (bypassing t2m and q2m for now)
    !should calculate qss, but it is not needed in moninedmf (q2m depends on qss - will implement later)

    wind10m = u_star(i)/vonKarman*log(10.0/roughness_length_m)
    u_fraction = sqrt(1.0 - v1(i)**2/wind1(i)**2)
    u10m(i) = u_fraction*wind10m
    v10m(i) = sqrt(wind10m**2 - u10m(i)**2)

    qss(i) = 0.0
    t2m(i) = 0.0
    q2m(i) = 0.0
  end do
  
  !GJF: The following code is from GFS_surface_composites.F90; only statements that are used in physics schemes outside of surface schemes are kept
  !GJF: Adding this code means that this scheme should be called before dcyc2t3
  do i = 1, im
    if (islmsk(i) == 1) then
      dry(i)    = .true.
      frland(i) = 1.0_kind_phys
      cice(i)   = 0.0_kind_phys
      icy(i)    = .false.
    else
      frland(i) = 0.0_kind_phys
      if (oceanfrac(i) > 0.0_kind_phys) then
        if (cice(i) >= min_seaice) then
          icy(i)   = .true.
          ! This cplice namelist option was added to deal with the
          ! situation of the FV3ATM-HYCOM coupling without an active sea
          ! ice (e.g., CICE6) component. By default, the cplice is true
          ! when cplflx is .true. (e.g., for the S2S application).
          ! Whereas, for the HAFS FV3ATM-HYCOM coupling, cplice is set as
          ! .false.. In the future HAFS FV3ATM-MOM6 coupling, the cplflx
          ! could be .true., while cplice being .false..
          if (cplice .and. cplflx)  then
            flag_cice(i)   = .true.
          else
            flag_cice(i)   = .false.
          endif
          islmsk(i) = 2
        else
          cice(i)        = 0.0_kind_phys
          flag_cice(i)   = .false.
          islmsk(i)      = 0
          icy(i)         = .false.
        endif
        if (cice(i) < 1.0_kind_phys) then
          wet(i) = .true. ! some open ocean
        endif
      else
        if (cice(i) >= min_lakeice) then
          icy(i) = .true.
          islmsk(i) = 2
        else
          cice(i)   = 0.0_kind_phys
          islmsk(i) = 0
          icy(i)    = .false.
        endif
        flag_cice(i)   = .false.
        if (cice(i) < 1.0_kind_phys) then
          wet(i) = .true. ! some open lake
        endif
      endif
    endif
    if (nint(slmsk(i)) /= 1) slmsk(i)  = islmsk(i)
  enddo
  
  do i = 1, im
    if (wet(i)) then
      tsfc_wat(i) = T_surf(i)
    end if
    if (dry(i)) then
      tsfcl(i)  = T_surf(i)
    end if
    if (icy(i)) then
      tisfc(i) = T_surf(i)
      tisfc(i) = max(timin, min(tisfc(i), tgice))
    end if
  end do

! to prepare to separate lake from ocean under water category
  do i = 1, im
    if ((wet(i) .or. icy(i)) .and. lakefrac(i) > 0.0_kind_phys) then
      if (lkm == 1 .and. lakefrac(i) >= 0.15 .and. lakedepth(i) > 1.0_kind_phys) then
        use_lake_model(i) = .true.
      else
        use_lake_model(i) = .false.
      endif
    else
      use_lake_model(i) = .false.
    endif
  enddo
!
  
  end subroutine scm_sfc_flux_spec_run

end module scm_sfc_flux_spec

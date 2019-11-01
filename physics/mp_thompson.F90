!>\file mp_thompson.F90
!! This file contains aerosol-aware Thompson MP scheme.


!>\defgroup aathompson Aerosol-Aware Thompson MP Module
!! This module contains the aerosol-aware Thompson microphysics scheme.
module mp_thompson

      use machine, only : kind_phys

      use module_mp_thompson, only : thompson_init, mp_gt_driver, thompson_finalize

      implicit none

      public :: mp_thompson_init, mp_thompson_run, mp_thompson_finalize

      private

      logical :: is_initialized = .False.

   contains

!> This subroutine is a wrapper around the actual mp_gt_driver().
#if 0
!! \section arg_table_mp_thompson_init Argument Table
!! \htmlinclude mp_thompson_init.html
!!
#endif
      subroutine mp_thompson_init(ncol, nlev, is_aerosol_aware, &
                                       nwfa2d, nifa2d, nwfa, nifa,   &
                                       mpicomm, mpirank, mpiroot,    &
                                       imp_physics,                  &
                                       imp_physics_thompson,         &
                                       threads, errmsg, errflg)

         implicit none

         ! Interface variables
         integer,                   intent(in)    :: ncol
         integer,                   intent(in)    :: nlev

         logical,                   intent(in)    :: is_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nwfa2d(1:ncol)
         real(kind_phys), optional, intent(inout) :: nifa2d(1:ncol)
         real(kind_phys), optional, intent(inout) :: nwfa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nifa(1:ncol,1:nlev)
         integer,                   intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         integer,                   intent(in)    :: threads
         integer,                   intent(in)    :: imp_physics
         integer,                   intent(in)    :: imp_physics_thompson
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables: dimensions used in thompson_init
         integer               :: ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return

         ! DH* temporary
         if (mpirank==mpiroot) then
            write(0,*) ' ----------------------------------------------------------------------------------------------------------------'
            write(0,*) ' --- WARNING --- the CCPP Thompson MP scheme is currently under development, use at your own risk --- WARNING ---'
            write(0,*) ' ----------------------------------------------------------------------------------------------------------------'
         end if
         ! *DH temporary

         if (imp_physics/=imp_physics_thompson) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Thompson MP"
            errflg = 1
            return
         end if

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

         if (is_aerosol_aware .and. present(nwfa2d) &
                              .and. present(nifa2d) &
                              .and. present(nwfa)   &
                              .and. present(nifa)   ) then
            ! Call init
            call thompson_init(nwfa2d=nwfa2d, nifa2d=nifa2d, nwfa=nwfa, nifa=nifa,   &
                               ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                               its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, &
                               mpicomm=mpicomm, mpirank=mpirank, mpiroot=mpiroot,    &
                               threads=threads, errmsg=errmsg, errflg=errflg)
            if (errflg /= 0) return
         else if (is_aerosol_aware) then
            write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_init:',                    &
                                       ' aerosol-aware microphysics require all of the following', &
                                       ' optional arguments: nifa2d, nwfa2d, nwfa, nifa'
            errflg = 1
            return
         else
            call thompson_init(ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                               its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, &
                               mpicomm=mpicomm, mpirank=mpirank, mpiroot=mpiroot,    &
                               threads=threads, errmsg=errmsg, errflg=errflg)
            if (errflg /= 0) return
         end if

         is_initialized = .true.

      end subroutine mp_thompson_init


#if 0
!> \section arg_table_mp_thompson_run Argument Table
!! \htmlinclude mp_thompson_run.html
!!
#endif
!>\ingroup aathompson
!>\section gen_thompson_hrrr Thompson MP General Algorithm
!>@{
      subroutine mp_thompson_run(ncol, nlev, con_g, con_rd,         &
                              spechum, qc, qr, qi, qs, qg, ni, nr,       &
                              is_aerosol_aware, nc, nwfa, nifa,          &
                              nwfa2d, nifa2d,                            &
                              tgrs, prsl, phii, omega, dtp,              &
                              prcp, rain, graupel, ice, snow, sr,        &
                              refl_10cm, do_radar_ref,                   &
                              re_cloud, re_ice, re_snow,                 &
                              mpicomm, mpirank, mpiroot,                 &
                              errmsg, errflg)

         implicit none

         ! Interface variables

         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         real(kind_phys),           intent(in   ) :: con_g
         real(kind_phys),           intent(in   ) :: con_rd
         ! Hydrometeors
         real(kind_phys),           intent(inout) :: spechum(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qc(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qr(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qi(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qs(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qg(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: ni(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nr(1:ncol,1:nlev)
         ! Aerosols
         logical,                   intent(in)    :: is_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nc(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nwfa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nifa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(in   ) :: nwfa2d(1:ncol)
         real(kind_phys), optional, intent(in   ) :: nifa2d(1:ncol)
         ! State variables and timestep information
         real(kind_phys),           intent(inout) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: phii(1:ncol,1:nlev+1)
         real(kind_phys),           intent(in   ) :: omega(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: dtp
         ! Precip/rain/snow/graupel fall amounts and fraction of frozen precip
         real(kind_phys),           intent(  out) :: prcp(1:ncol)
         real(kind_phys),           intent(  out) :: rain(1:ncol)
         real(kind_phys),           intent(  out) :: graupel(1:ncol)
         real(kind_phys),           intent(  out) :: ice(1:ncol)
         real(kind_phys),           intent(  out) :: snow(1:ncol)
         real(kind_phys),           intent(  out) :: sr(1:ncol)
         ! Radar reflectivity
         real(kind_phys),           intent(  out) :: refl_10cm(1:ncol,1:nlev)
         logical,         optional, intent(in   ) :: do_radar_ref
         ! Cloud effective radii
         real(kind_phys), optional, intent(  out) :: re_cloud(1:ncol,1:nlev)
         real(kind_phys), optional, intent(  out) :: re_ice(1:ncol,1:nlev)
         real(kind_phys), optional, intent(  out) :: re_snow(1:ncol,1:nlev)
         ! MPI information
         integer,                   intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables

         ! Air density
         real(kind_phys) :: rho(1:ncol,1:nlev)              !< kg m-3
         ! Hydrometeors
         real(kind_phys) :: qv_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qc_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qr_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qi_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qs_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qg_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         ! Vertical velocity and level width
         real(kind_phys) :: w(1:ncol,1:nlev)                !< m s-1
         real(kind_phys) :: dz(1:ncol,1:nlev)               !< m
         ! Rain/snow/graupel fall amounts
         real(kind_phys) :: rain_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: graupel_mp(1:ncol)              ! mm, dummy, not used
         real(kind_phys) :: ice_mp(1:ncol)                  ! mm, dummy, not used
         real(kind_phys) :: snow_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: delta_rain_mp(1:ncol)           ! mm
         real(kind_phys) :: delta_graupel_mp(1:ncol)        ! mm
         real(kind_phys) :: delta_ice_mp(1:ncol)            ! mm
         real(kind_phys) :: delta_snow_mp(1:ncol)           ! mm
         ! Radar reflectivity
         logical         :: diagflag                        ! must be true if do_radar_ref is true, not used otherwise
         integer         :: do_radar_ref_mp                 ! integer instead of logical do_radar_ref
         ! Effective cloud radii
         logical         :: do_effective_radii
         integer         :: has_reqc
         integer         :: has_reqi
         integer         :: has_reqs
         ! Dimensions used in mp_gt_driver
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Check initialization state
         if (.not.is_initialized) then
            write(errmsg, fmt='((a))') 'mp_thompson_run called before mp_thompson_init'
            errflg = 1
            return
         end if

         !> - Convert specific humidity/moist mixing ratios to dry mixing ratios
         qv_mp = spechum/(1.0_kind_phys-spechum)
         qc_mp = qc/(1.0_kind_phys-spechum)
         qr_mp = qr/(1.0_kind_phys-spechum)
         qi_mp = qi/(1.0_kind_phys-spechum)
         qs_mp = qs/(1.0_kind_phys-spechum)
         qg_mp = qg/(1.0_kind_phys-spechum)

         if (is_aerosol_aware .and. .not. (present(nc)     .and. &
                                           present(nwfa)   .and. &
                                           present(nifa)   .and. &
                                           present(nwfa2d) .and. &
                                           present(nifa2d)       )) then
            write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_run:',  &
                                       ' aerosol-aware microphysics require all of the', &
                                       ' following optional arguments:', &
                                       ' nc, nwfa, nifa, nwfa2d, nifa2d'
            errflg = 1
            return
         end if

         !> - Density of air in kg m-3
         rho = prsl/(con_rd*tgrs)

         !> - Convert omega in Pa s-1 to vertical velocity w in m s-1
         w = -omega/(rho*con_g)

         !> - Layer width in m from geopotential in m2 s-2
         dz = (phii(:,2:nlev+1) - phii(:,1:nlev)) / con_g

         ! Accumulated values inside Thompson scheme, not used;
         ! only use delta and add to inout variables (different units)
         rain_mp          = 0
         graupel_mp       = 0
         ice_mp           = 0
         snow_mp          = 0
         delta_rain_mp    = 0
         delta_graupel_mp = 0
         delta_ice_mp     = 0
         delta_snow_mp    = 0

         ! Flags for calculating radar reflectivity; diagflag is redundant
         if (do_radar_ref) then
             diagflag = .true.
             do_radar_ref_mp = 1
         else
             diagflag = .false.
             do_radar_ref_mp = 0
         end if

         if (present(re_cloud) .and. present(re_ice) .and. present(re_snow)) then
             do_effective_radii = .true.
             has_reqc = 1
             has_reqi = 1
             has_reqs = 1
             ! Initialize to zero, intent(out) variables
             re_cloud = 0
             re_ice   = 0
             re_snow  = 0
         else if (.not.present(re_cloud) .and. .not.present(re_ice) .and. .not.present(re_snow)) then
             do_effective_radii = .false.
             has_reqc = 0
             has_reqi = 0
             has_reqs = 0
         else
             write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_run:',  &
                                        ' all or none of the following optional', &
                                        ' arguments are required: re_cloud, re_ice, re_snow'
             errflg = 1
             return
         end if

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


         !> - Call mp_gt_driver() with or without aerosols
         if (is_aerosol_aware) then
            call mp_gt_driver(qv=qv_mp, qc=qc_mp, qr=qr_mp, qi=qi_mp, qs=qs_mp, qg=qg_mp,    &
                              ni=ni, nr=nr, nc=nc,                                           &
                              nwfa=nwfa, nifa=nifa, nwfa2d=nwfa2d, nifa2d=nifa2d,            &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtp,                        &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              icenc=ice_mp, icencv=delta_ice_mp,                             &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm,                                           &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow,             &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,          &
                              errmsg=errmsg, errflg=errflg)

         else
            call mp_gt_driver(qv=qv_mp, qc=qc_mp, qr=qr_mp, qi=qi_mp, qs=qs_mp, qg=qg_mp,    &
                              ni=ni, nr=nr, nc=nc,                                           &
                              tt=tgrs, p=prsl, w=w, dz=dz, dt_in=dtp,                        &
                              rainnc=rain_mp, rainncv=delta_rain_mp,                         &
                              snownc=snow_mp, snowncv=delta_snow_mp,                         &
                              icenc=ice_mp, icencv=delta_ice_mp,                             &
                              graupelnc=graupel_mp, graupelncv=delta_graupel_mp, sr=sr,      &
                              refl_10cm=refl_10cm,                                           &
                              diagflag=diagflag, do_radar_ref=do_radar_ref_mp,               &
                              re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow,             &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,          &
                              errmsg=errmsg, errflg=errflg)
         end if
         if (errflg/=0) return

         !> - Convert dry mixing ratios to specific humidity/moist mixing ratios
         spechum = qv_mp/(1.0_kind_phys+qv_mp)
         qc      = qc_mp/(1.0_kind_phys+qv_mp)
         qr      = qr_mp/(1.0_kind_phys+qv_mp)
         qi      = qi_mp/(1.0_kind_phys+qv_mp)
         qs      = qs_mp/(1.0_kind_phys+qv_mp)
         qg      = qg_mp/(1.0_kind_phys+qv_mp)


         !> - Convert rainfall deltas from mm to m (on physics timestep); add to inout variables
         ! "rain" in Thompson MP refers to precipitation (total of liquid rainfall+snow+graupel+ice)
         prcp    = max(0.0, delta_rain_mp/1000.0_kind_phys)
         graupel = max(0.0, delta_graupel_mp/1000.0_kind_phys)
         ice     = max(0.0, delta_ice_mp/1000.0_kind_phys)
         snow    = max(0.0, delta_snow_mp/1000.0_kind_phys)
         rain    = max(0.0, (delta_rain_mp - (delta_graupel_mp + delta_ice_mp + delta_snow_mp))/1000.0_kind_phys)

      end subroutine mp_thompson_run
!>@}

#if 0
!! \section arg_table_mp_thompson_finalize Argument Table
!! \htmlinclude mp_thompson_finalize.html
!!
#endif
      subroutine mp_thompson_finalize(errmsg, errflg)

         implicit none

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (.not.is_initialized) return

         call thompson_finalize()

         is_initialized = .false.

      end subroutine mp_thompson_finalize

end module mp_thompson

!>\file mp_nsslg.F90
!! This file contains NSSL 2-moment MP scheme.


!>\defgroup aansslg NSSL MP Module
!! This module contains the NSSL microphysics scheme.
module mp_nsslg

    use machine, only : kind_phys, kind_real
    use module_mp_nssl_2mom, only : nssl_2mom_init, nssl_2mom_driver

    implicit none

    public :: mp_nsslg_init, mp_nsslg_run, mp_nsslg_finalize

    private
    logical :: is_initialized = .False.
    real :: nssl_qccn

    contains

!> This subroutine is a wrapper around the nssl_2mom_init().
!! \section arg_table_mp_nsslg_init Argument Table
!! \htmlinclude mp_nsslg_init.html
!!
    subroutine mp_nsslg_init(ncol, nlev, errflg, errmsg,threads, &
                                       mpicomm, mpirank, mpiroot,    &
                                       imp_physics,                  &
                                       imp_physics_nssl2m, imp_physics_nssl2mccn,         &
                                       nssl_cccn, nssl_alphah, nssl_alphahl, nssl_hail_on)

        implicit none
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         integer,                   intent(in)    :: ncol
         integer,                   intent(in)    :: nlev

         integer,                   intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         integer,                   intent(in)    :: threads
         integer,                   intent(in)    :: imp_physics
         integer,                   intent(in)    :: imp_physics_nssl2m, imp_physics_nssl2mccn
         real(kind_phys),           intent(in)    :: nssl_cccn, nssl_alphah, nssl_alphahl
         logical,                   intent(in)    :: nssl_hail_on

         ! Local variables: dimensions used in nssl_init
         integer               :: ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte
         real :: nssl_params(20)
         integer :: ihailv
         


        errflg = 0
        errmsg = ''


        if (is_initialized) return

         if (mpirank==mpiroot) then
            write(0,*) ' ----------------------------------------------------------------------------------------------------------------'
            write(0,*) ' --- WARNING! --- the CCPP NSSL MP scheme is currently under development, use at your own risk --- WARNING ---'
            write(0,*) ' ----------------------------------------------------------------------------------------------------------------'
            write(6,*) ' ----------------------------------------------------------------------------------------------------------------'
            write(6,*) ' --- WARNING! --- the CCPP NSSL MP scheme is currently under development, use at your own risk --- WARNING ---'
            write(6,*) ' ----------------------------------------------------------------------------------------------------------------'
         end if
        
!        IF ( kind_phys /= kind_real ) THEN
!          errflg = 1
!          write(errmsg,'(a)') 'NSSL MP does not yet work for double precision. Compile for single precision'
!          return
!        ENDIF

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

         is_initialized = .true.

         nssl_params(:) = 0.0
         nssl_params(1)  = nssl_cccn
         nssl_params(2)  = nssl_alphah
         nssl_params(3)  = nssl_alphahl
         nssl_params(4)  = 4.e5 ! nssl_cnoh
         nssl_params(5)  = 4.e4 ! nssl_cnohl
         nssl_params(6)  = 4.e5 ! nssl_cnor
         nssl_params(7)  = 4.e6 ! nssl_cnos
         nssl_params(8)  = 500. ! nssl_rho_qh
         nssl_params(9)  = 800. ! nssl_rho_qhl
         nssl_params(10) = 100. ! nssl_rho_qs
         nssl_params(11) = 0 ! nssl_ipelec_tmp
         nssl_params(12) = 11 ! nssl_isaund
         nssl_params(13) = 0 ! 1= turn on cccna; 0 = turn off
         
         nssl_qccn = nssl_cccn/1.225
         if (mpirank==mpiroot) then
          write(*,*) 'nssl_init: nssl_qccn = ',nssl_qccn
         endif
         
         IF ( nssl_hail_on ) THEN
           ihailv = 1
         ELSE
           ihailv = -1
         ENDIF

         IF ( imp_physics == imp_physics_nssl2m ) THEN
!           write(0,*) 'call nssl_2mom_init'
         CALL nssl_2mom_init(ims,ime, jms,jme, kms,kme,nssl_params,ipctmp=5,mixphase=0,ihvol=ihailv)
!           write(0,*) 'done nssl_2mom_init'
         ELSEIF ( imp_physics == imp_physics_nssl2mccn ) THEN
!           write(0,*) 'call nssl_2mom_init ccn'
         CALL nssl_2mom_init(ims,ime, jms,jme, kms,kme,nssl_params,ipctmp=5,mixphase=0,ihvol=ihailv)
!           write(0,*) 'done nssl_2mom_init ccn'
         ELSE
!           write(0,*) 'call nssl_2mom_init ccn: imp_physics, imp_physics_nssl2mccn = ',imp_physics, imp_physics_nssl2mccn
         CALL nssl_2mom_init(ims,ime, jms,jme, kms,kme,nssl_params,ipctmp=5,mixphase=0,ihvol=ihailv)
!           write(0,*) 'done nssl_2mom_init ccn'
         ENDIF

    end subroutine mp_nsslg_init

!>\ingroup aansslg
!>\section gen_nsslg NSSL MP General Algorithm
!>@{
!> \section arg_table_mp_nsslg_run Argument Table
!! \htmlinclude mp_nsslg_run.html
!!
    subroutine mp_nsslg_run(ncol, nlev, con_g, con_rd, &
!                             spechum, cccn, qc, qr, qi, qs, qh, qhl,         &
                             spechum, cccn, cccna, qc, qr, qi, qs, qh, qhl,         &
                             ccw, crw, cci, csw, chw, chl, vh, vhl,          &
                              tgrs, prslk, prsl, phii, omega, dtp,           &
                              prcp, rain, graupel, ice, snow, sr,            &
                             refl_10cm, do_radar_ref, first_time_step,       &
                             re_cloud, re_ice, re_snow, re_rain,             &
                             imp_physics,                                    &
                             imp_physics_nssl2m, imp_physics_nssl2mccn,      &
                             nssl_hail_on, nssl_invertccn, ntccn, ntccna,    &
                             errflg, errmsg)
        implicit none
        integer, intent(in) :: ncol, nlev
         real(kind_phys),           intent(in   ) :: con_g
         real(kind_phys),           intent(in   ) :: con_rd
         ! Hydrometeors
         real(kind_phys),           intent(inout) :: spechum(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: cccn(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: cccna(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qc(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qr(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qi(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qs(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qh(1:ncol,1:nlev)  ! graupel
         real(kind_phys),           intent(inout) :: qhl(1:ncol,1:nlev) ! hail
         real(kind_phys),           intent(inout) :: ccw(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: crw(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: cci(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: csw(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: chw(1:ncol,1:nlev) ! graupel number
         real(kind_phys),           intent(inout) :: chl(1:ncol,1:nlev) ! hail number
         real(kind_phys),           intent(inout) :: vh(1:ncol,1:nlev)  ! graupel volume
         real(kind_phys),           intent(inout) :: vhl(1:ncol,1:nlev) ! hail volume
         ! State variables and timestep information
         real(kind_phys),           intent(inout) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prslk(1:ncol,1:nlev)
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
         logical,                   intent(in   ) :: do_radar_ref, first_time_step
         ! Cloud effective radii
         real(kind_phys), optional, intent(  out) :: re_cloud(1:ncol,1:nlev)
         real(kind_phys), optional, intent(  out) :: re_ice(1:ncol,1:nlev)
         real(kind_phys), optional, intent(  out) :: re_snow(1:ncol,1:nlev)
         real(kind_phys), optional, intent(  out) :: re_rain(1:ncol,1:nlev)
         integer,                   intent(in)    :: imp_physics
         integer,                   intent(in)    :: imp_physics_nssl2m, imp_physics_nssl2mccn
         logical,                   intent(in)    :: nssl_hail_on, nssl_invertccn
         integer,                   intent(in)    :: ntccn, ntccna
        
        integer,          intent(out)   :: errflg
        character(len=*), intent(out)   :: errmsg


         ! Local variables

         ! Air density
         real(kind_phys) :: rho(1:ncol,1:nlev)              !< kg m-3
         ! Hydrometeors
         real(kind_phys) :: qv_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qc_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qr_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qi_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qs_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qh_mp(1:ncol,1:nlev)            !< kg kg-1 (graupel dry mixing ratio)
         real(kind_phys) :: qhl_mp(1:ncol,1:nlev)           !< kg kg-1 (hail dry mixing ratio)
         real(kind_phys) :: cn_mp(1:ncol,1:nlev) 
         real(kind_phys) :: cna_mp(1:ncol,1:nlev) 
         ! create temporaries for hail in case it does not exist
         real(kind_phys) :: chl_mp(1:ncol,1:nlev)           !< kg-1 (number mixing ratio)
         real(kind_phys) :: vhl_mp(1:ncol,1:nlev)           !< m3 kg-1 (volume mixing ratio)
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

         real(kind_phys) :: xrain_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: xgraupel_mp(1:ncol)              ! mm, dummy, not used
         real(kind_phys) :: xice_mp(1:ncol)                  ! mm, dummy, not used
         real(kind_phys) :: xsnow_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: xdelta_rain_mp(1:ncol)           ! mm
         real(kind_phys) :: xdelta_graupel_mp(1:ncol)        ! mm
         real(kind_phys) :: xdelta_ice_mp(1:ncol)            ! mm
         real(kind_phys) :: xdelta_snow_mp(1:ncol)           ! mm

         ! Radar reflectivity
         logical         :: diagflag                        ! must be true if do_radar_ref is true, not used otherwise
         integer         :: do_radar_ref_mp                 ! integer instead of logical do_radar_ref
         ! Effective cloud radii
         logical         :: do_effective_radii
         real(kind_phys) :: re_cloud_mp(1:ncol,1:nlev)      ! m
         real(kind_phys) :: re_ice_mp(1:ncol,1:nlev)        ! m
         real(kind_phys) :: re_snow_mp(1:ncol,1:nlev)       ! m
         integer         :: has_reqc
         integer         :: has_reqi
         integer         :: has_reqs
         ! Dimensions used in driver
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte, i,j,k
         integer :: itimestep = 0 ! timestep counter
         integer :: ntmul, n
         real, parameter    :: dtpmax = 300. ! 600. ! 120.
         real(kind_phys)    :: dtptmp
         integer, parameter :: ndebug = 0
         logical, parameter :: convertdry = .true.
         logical :: invertccn
         


        errflg = 0
        errmsg = ''

        IF ( ndebug > 1 ) write(0,*) 'In physics nsslg_run'


         ! Check initialization state
         if (.not.is_initialized) then
            write(errmsg, fmt='((a))') 'mp_nssl_run called before mp_nssl_init'
            errflg = 1
            return
         end if
         
         invertccn = nssl_invertccn

         !> - Convert specific humidity/moist mixing ratios to dry mixing ratios
         qv_mp = spechum/(1.0_kind_phys-spechum)
         IF ( convertdry ) THEN
         qc_mp = qc/(1.0_kind_phys-spechum)
         qr_mp = qr/(1.0_kind_phys-spechum)
         qi_mp = qi/(1.0_kind_phys-spechum)
         qs_mp = qs/(1.0_kind_phys-spechum)
         qh_mp = qh/(1.0_kind_phys-spechum)
         IF ( nssl_hail_on ) THEN
           qhl_mp = qhl/(1.0_kind_phys-spechum)
         ENDIF
         ELSE
!         qv_mp = spechum ! /(1.0_kind_phys-spechum)
         qc_mp = qc ! /(1.0_kind_phys-spechum)
         qr_mp = qr ! /(1.0_kind_phys-spechum)
         qi_mp = qi ! /(1.0_kind_phys-spechum)
         qs_mp = qs ! /(1.0_kind_phys-spechum)
         qh_mp = qh ! /(1.0_kind_phys-spechum)
         IF ( nssl_hail_on ) THEN
           qhl_mp = qhl ! /(1.0_kind_phys-spechum)
         ENDIF
         
         ENDIF
         
         IF ( nssl_hail_on ) THEN
           chl_mp = chl
           vhl_mp = vhl
         ELSE
           qhl_mp = 0
           chl_mp = 0
           vhl_mp = 0
         ENDIF

         
         !> - Density of air in kg m-3
         rho = prsl/(con_rd*tgrs)

         !> - Convert omega in Pa s-1 to vertical velocity w in m s-1
         w = -omega/(rho*con_g)

         !> - Layer width in m from geopotential in m2 s-2
         dz = (phii(:,2:nlev+1) - phii(:,1:nlev)) / con_g

         ! Accumulated values inside scheme, not used;
         ! only use delta and add to inout variables (different units)
         rain_mp          = 0
         graupel_mp       = 0
         ice_mp           = 0
         snow_mp          = 0
         delta_rain_mp    = 0
         delta_graupel_mp = 0
         delta_ice_mp     = 0
         delta_snow_mp    = 0
         xrain_mp          = 0
         xgraupel_mp       = 0
         xice_mp           = 0
         xsnow_mp          = 0
         xdelta_rain_mp    = 0
         xdelta_graupel_mp = 0
         xdelta_ice_mp     = 0
         xdelta_snow_mp    = 0

         IF ( ndebug >= 1 ) THEN
         write(*,*) 'Max q before micro'
         write(*,*) 'qc = ',1000.*maxval(qc_mp)
         write(*,*) 'qr = ',1000.*maxval(qr_mp)
         write(*,*) 'qi = ',1000.*maxval(qi_mp)
         write(*,*) 'qs = ',1000.*maxval(qs_mp)
         write(*,*) 'qh = ',1000.*maxval(qh_mp)
         IF ( nssl_hail_on ) write(*,*) 'qhl = ',1000.*maxval(qhl_mp)
         write(*,*) 'ccw = ',1.e-6*maxval(ccw*rho)
         ENDIF

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
         else if (.not.present(re_cloud) .and. .not.present(re_ice) .and. .not.present(re_snow)) then
             do_effective_radii = .false.
             has_reqc = 0
             has_reqi = 0
             has_reqs = 0
         else
             write(errmsg,fmt='(*(a))') 'Logic error in mp_nssl_run:',  &
                                        ' all or none of the following optional', &
                                        ' arguments are required: re_cloud, re_ice, re_snow'
             errflg = 1
             return
         end if
         ! Initialize to zero, intent(out) variables
         re_cloud_mp = 0
         re_ice_mp   = 0
         re_snow_mp  = 0

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


       IF ( ndebug > 1 )  write(0,*) 'call nssl_2mom_driver'

        IF ( dtp > 1.5*dtpmax ) THEN
           ntmul = Nint( dtp/dtpmax )
           dtptmp = dtp/ntmul
        ELSE
           dtptmp = dtp
           ntmul = 1
        ENDIF
        
        IF ( first_time_step ) THEN
          itimestep = 2
          IF ( imp_physics == imp_physics_nssl2mccn ) THEN
            IF ( invertccn ) THEN
              cccn = 0
              !cccn = nssl_qccn
            ELSE
              cccn = nssl_qccn
            ENDIF
          ENDIF
        ELSE
          itimestep = 2
        ENDIF
 
 
       IF ( imp_physics == imp_physics_nssl2mccn ) THEN
         IF ( invertccn ) THEN
!            cn_mp = Max(0.0, nssl_qccn - Max(0.0,cccn))
              DO k = 1,nlev
               DO i = 1,ncol
                 cn_mp(i,k) = Max(0.0, nssl_qccn - Max(0.0, cccn(i,k)) )
!                 cn_mp(i,k) = Min(nssl_qccn, nssl_qccn - cccn(i,k) ) 
               ENDDO
              ENDDO
            !  DO k = 1,nlev
            !   DO i = 1,ncol
            !     cccn(i,k) = Max(0.0, nssl_qccn - cn_mp(i,k) )
            !     cn_mp(i,k) = cccn(i,k)
            !   ENDDO
            !  ENDDO
         ELSE
            cn_mp = cccn
         ENDIF
          IF ( ntccna > 0 ) THEN
!         cna_mp = cccna
          ELSE 
            cna_mp = 0
          ENDIF
        ENDIF
       
       
        DO n = 1,ntmul
        
        itimestep = itimestep + 1



         IF ( imp_physics == imp_physics_nssl2mccn ) THEN


         CALL nssl_2mom_driver(                          &
                    ITIMESTEP=itimestep,                &
                  !   TH=th,                              &
                     tt=tgrs,                          &
                     QV=qv_mp,                         &
                     QC=qc_mp,                         &
                     QR=qr_mp,                         &
                     QI=qi_mp,                         &
                     QS=qs_mp,                         &
                     QH=qh_mp,                         &
                     QHL=qhl_mp,                        &
                     CCW=ccw,                    &
                     CRW=crw,                       &
                     CCI=cci,                       &
                     CSW=csw,                       &
                     CHW=chw,                       &
                     CHL=chl_mp,                       &
                     VHW=vh,                     &
                     VHL=vhl_mp,                     &
                     cn=cn_mp,                        &
!                     cna=cna_mp, f_cna=( ntccna > 0 ),  & ! for future use
                      cna=cna_mp, f_cna=.false. ,           &
                    PII=prslk,                         &
                     P=prsl,                                &
                     W=w,                                &
                     DZ=dz,                              &
                     DTP=dtptmp,                         &
                     DN=rho,                             &
                     rainnc=xrain_mp, rainncv=xdelta_rain_mp,                         &
                     snownc=xsnow_mp, snowncv=xdelta_snow_mp,                         &
!                     icenc=ice_mp, icencv=delta_ice_mp,                             &
                     GRPLNC=xgraupel_mp, GRPLNCV=xdelta_graupel_mp, sr=sr,      &
                     dbz      = refl_10cm,               &
!                     nssl_progn=.false.,                       &
                     diagflag = diagflag,                &
                     re_cloud=re_cloud_mp,                  &
                     re_ice=re_ice_mp,                      &
                     re_snow=re_snow_mp,                    &
                     has_reqc=has_reqc,                  & ! ala G. Thompson
                     has_reqi=has_reqi,                  & ! ala G. Thompson
                     has_reqs=has_reqs,                  & ! ala G. Thompson
                  IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde, &
                  IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme, &
                  ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte  &
                                                                    )


           ELSE

         CALL nssl_2mom_driver(                          &
                    ITIMESTEP=itimestep,                &
                  !   TH=th,                              &
                     tt=tgrs,                          &
                     QV=qv_mp,                         &
                     QC=qc_mp,                         &
                     QR=qr_mp,                         &
                     QI=qi_mp,                         &
                     QS=qs_mp,                         &
                     QH=qh_mp,                         &
                     QHL=qhl_mp,                        &
!                     CCW=qnc_mp,                       &
                     CCW=ccw,                    &
                     CRW=crw,                       &
                     CCI=cci,                       &
                     CSW=csw,                       &
                     CHW=chw,                       &
                     CHL=chl_mp,                       &
                     VHW=vh,                     &
                     VHL=vhl_mp,                     &
                !     cn=cccn,                        &
                     PII=prslk,                         &
                     P=prsl,                                &
                     W=w,                                &
                     DZ=dz,                              &
                     DTP=dtptmp,                         &
                     DN=rho,                             &
                     rainnc=xrain_mp, rainncv=xdelta_rain_mp,                         &
                     snownc=xsnow_mp, snowncv=xdelta_snow_mp,                         &
!                     icenc=ice_mp, icencv=delta_ice_mp,                             &
                     GRPLNC=xgraupel_mp, GRPLNCV=xdelta_graupel_mp, sr=sr,      &
                     dbz      = refl_10cm,               &
!                     nssl_progn=.false.,                       &
                     diagflag = diagflag,                &
                     re_cloud=re_cloud_mp,                  &
                     re_ice=re_ice_mp,                      &
                     re_snow=re_snow_mp,                    &
                     has_reqc=has_reqc,                  & ! ala G. Thompson
                     has_reqi=has_reqi,                  & ! ala G. Thompson
                     has_reqs=has_reqs,                  & ! ala G. Thompson
                  IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde, &
                  IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme, &
                  ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte  &
                                                                    )
           
           ENDIF
           
           
           DO i = 1,ncol
             delta_rain_mp(i) = delta_rain_mp(i) + xdelta_rain_mp(i)
             delta_graupel_mp(i) = delta_graupel_mp(i) + xdelta_graupel_mp(i)
             delta_ice_mp(i) = delta_ice_mp(i) + xdelta_ice_mp(i)
             delta_snow_mp(i) = delta_snow_mp(i) + xdelta_snow_mp(i)
           ENDDO

          ENDDO


         IF ( imp_physics == imp_physics_nssl2mccn ) THEN
           IF ( invertccn ) THEN
              !cccn = Max(0.0, nssl_qccn - cn_mp )
              DO k = 1,nlev
               DO i = 1,ncol
!                 cccn(i,k) = Max(0.0, nssl_qccn - cn_mp(i,k) )
                 cccn(i,k) = nssl_qccn - cn_mp(i,k) 
               ENDDO
              ENDDO
           ELSE
              cccn = cn_mp
           ENDIF
!           cccna = cna_mp
          ENDIF
          
! test code
!          IF ( ntccna > 1 .and. do_effective_radii ) THEN
!            cccna = re_ice_mp*1.0E6_kind_phys
!          ENDIF

        IF ( ndebug > 1 ) write(0,*) 'done nssl_2mom_driver'

         if (errflg/=0) return

         IF ( ndebug >= 1 ) THEN
         write(*,*) 'Max q after micro'
         write(*,*) 'qc = ',1000.*maxval(qc_mp)
         write(*,*) 'qr = ',1000.*maxval(qr_mp)
         write(*,*) 'qi = ',1000.*maxval(qi_mp)
         write(*,*) 'qs = ',1000.*maxval(qs_mp)
         write(*,*) 'qh = ',1000.*maxval(qh_mp)
         IF ( nssl_hail_on ) THEN
           write(*,*) 'qhl = ',1000.*maxval(qhl_mp)
         ENDIF
         write(*,*) 'ccw = ',1.e-6*maxval(ccw*rho)
           IF ( 1000.*maxval(qc_mp) > 0.5 .or. 1000.*maxval(qi_mp) > 0.09 .or. 1000.*maxval(qs_mp) > 0.1 ) THEN
             IF ( imp_physics == imp_physics_nssl2mccn ) THEN
             write(*,*) 'qc, ccn, ccw, tt, qi+qs by height'
             DO k = 1,nlev
               write(*,*) qc_mp(1,k)*1000., cccn(1,k)*rho(1,k)*1.e-6, ccw(1,k)*rho(1,k)*1.e-6, tgrs(1,k), (qs_mp(1,k)+qi_mp(1,k))*1000. ! cccn(1,k)*1.e-6
             ENDDO
             ELSE
             write(*,*) 'qc, ccn, ccw, tt, qi+qs by height'
             DO k = 1,nlev
               write(*,*) qc_mp(1,k)*1000., cccn(1,k)*rho(1,k)*1.e-6, 0.0, tgrs(1,k), (qs_mp(1,k)+qi_mp(1,k))*1000. ! cccn(1,k)*1.e-6
             ENDDO
             ENDIF
           ENDIF
         ENDIF

         IF ( nssl_hail_on ) THEN
           chl = chl_mp
           vhl = vhl_mp
         ENDIF

         !> - Convert dry mixing ratios to specific humidity/moist mixing ratios
         spechum = qv_mp/(1.0_kind_phys+qv_mp)
         IF ( convertdry ) THEN
         qc      = qc_mp/(1.0_kind_phys+qv_mp)
         qr      = qr_mp/(1.0_kind_phys+qv_mp)
         qi      = qi_mp/(1.0_kind_phys+qv_mp)
         qs      = qs_mp/(1.0_kind_phys+qv_mp)
         qh      = qh_mp/(1.0_kind_phys+qv_mp)
         IF ( nssl_hail_on ) THEN
          qhl     = qhl_mp/(1.0_kind_phys+qv_mp)
         ENDIF
         ELSE
!         spechum = qv_mp ! /(1.0_kind_phys+qv_mp)
         qc      = qc_mp ! /(1.0_kind_phys+qv_mp)
         qr      = qr_mp ! /(1.0_kind_phys+qv_mp)
         qi      = qi_mp ! /(1.0_kind_phys+qv_mp)
         qs      = qs_mp ! /(1.0_kind_phys+qv_mp)
         qh      = qh_mp ! /(1.0_kind_phys+qv_mp)
         IF ( nssl_hail_on ) THEN
          qhl     = qhl_mp ! /(1.0_kind_phys+qv_mp)
         ENDIF
         
         ENDIF

!        write(0,*) 'mp_nsslg: done q'

         !> - Convert rainfall deltas from mm to m (on physics timestep); add to inout variables
         ! "rain" in NSSL MP refers to precipitation (total of liquid rainfall+snow+graupel+ice)
         
         prcp    = max(0.0, delta_rain_mp/1000.0_kind_phys)
         graupel = max(0.0, delta_graupel_mp/1000.0_kind_phys)
         ice     = max(0.0, delta_ice_mp/1000.0_kind_phys)
         snow    = max(0.0, delta_snow_mp/1000.0_kind_phys)
         rain    = max(0.0, delta_rain_mp - (delta_graupel_mp + delta_ice_mp + delta_snow_mp)/1000.0_kind_phys)

!        write(0,*) 'mp_nsslg: done precip'

         if (do_effective_radii) then
            ! Convert m to micron
            re_cloud = re_cloud_mp*1.0E6_kind_phys
            re_ice   = re_ice_mp*1.0E6_kind_phys
            re_snow  = re_snow_mp*1.0E6_kind_phys
!            re_rain  = 1.0E3_kind_phys
         end if

        IF ( ndebug > 1 ) write(0,*) 'mp_nsslg: end'

    end subroutine mp_nsslg_run
!>@}

#if 0
!! \section arg_table_mp_nsslg_finalize Argument Table
!! \htmlinclude mp_nsslg_finalize.html
!!
#endif
    subroutine mp_nsslg_finalize(errflg, errmsg)
        implicit none
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

        errflg = 0
        errmsg = ''


    end subroutine mp_nsslg_finalize

end module mp_nsslg

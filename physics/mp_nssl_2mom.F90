!>\file mp_nssl_2mom.F90
!! This file contains the NSSL 2-moment microphysics scheme.
!!  Added by Chunxi Zhang and Tim Supinie of CAPS in May 2020


!>\defgroup nssl2mom NSSL 2-moment microphysics scheme
!! This module contains the NSSL 2-moment microphysics scheme.
module mp_nssl_2mom

      use machine, only : kind_phys

      use module_mp_nssl_2mom, only : nssl_2mom_init, nssl_2mom_driver

      implicit none

      public :: mp_nssl_2mom_init, mp_nssl_2mom_run, mp_nssl_2mom_finalize

      private

      logical :: is_initialized = .False.

! #####################################################################
contains

!> This subroutine is a wrapper around the actual nssl_2mom_init.
!! \section arg_table_mp_nssl_2mom_init Argument Table
!! \htmlinclude mp_nssl_2mom_init.html
!!
  subroutine mp_nssl_2mom_init(mpirank, mpiroot, con_pi, con_g, con_rd, con_cp, con_rv, con_t0c,  &
                               con_cliq, con_csol, con_eps, imp_physics, imp_physics_nssl, errmsg, errflg)

  implicit none

   integer,                   intent(in)    :: mpirank
   integer,                   intent(in)    :: mpiroot
   real(kind=kind_phys),      intent(in)    :: con_pi, con_g, con_rd, con_cp, con_rv, con_t0c, &
                                               con_cliq, con_csol, con_eps
   integer,                   intent(in)    :: imp_physics
   integer,                   intent(in)    :: imp_physics_nssl
   character(len=*),          intent(  out) :: errmsg
   integer,                   intent(  out) :: errflg

   integer :: ipctmp,mixphase,ihvol
     double precision :: arg
     real(kind=kind_phys)    :: temq
     integer :: igam
     integer :: i,il,j,l
     integer :: ltmp
     integer :: isub
     real(kind=kind_phys)    :: bxh,bxhl

     real(kind=kind_phys)    :: alp,ratio,x,y,y7
     real(kind=kind_phys), dimension(20) :: nssl_params

 ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

       if (is_initialized) return

       if (mpirank==mpiroot) then
            write(0,*) '----------------------------------------------------------------------------------------------------------------'
            write(0,*) ' --- WARNING --- the CCPP NSSL MP scheme is currently under development, use at your own risk --- WARNING ---'
            write(0,*) '----------------------------------------------------------------------------------------------------------------'
       end if

       if (imp_physics/=imp_physics_nssl) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from NSSL"
            errflg = 1
            return
       end if

! set some values here
       ipctmp   = 5
       mixphase = 0
       ihvol    = 1

       nssl_params = 0.
       nssl_params(1)  = 0.7e9 
       nssl_params(2)  = 0.
       nssl_params(3)  = 2.
       nssl_params(4)  = 4.e5
       nssl_params(5)  = 4.e4
       nssl_params(6)  = 8.e5
       nssl_params(7)  = 3.e6
       nssl_params(8)  = 500.
       nssl_params(9)  = 900.
       nssl_params(10) = 100.
       nssl_params(11) = 0.
       nssl_params(12) = 0.

       CALL nssl_2mom_init(nssl_params,ipctmp,mixphase,ihvol, con_pi, con_g, con_rd, con_cp, con_rv, con_t0c, &
                           con_cliq, con_csol, con_eps, errmsg, errflg)
       if (errflg /= 0) return

       is_initialized = .true.

     end subroutine mp_nssl_2mom_init

!! \section arg_table_mp_nssl_2mom_finalize Argument Table
!! \htmlinclude mp_nssl_2mom_finalize.html
!!
     subroutine mp_nssl_2mom_finalize(errmsg, errflg)

         implicit none

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

      end subroutine mp_nssl_2mom_finalize

!> \section arg_table_mp_nssl_2mom_run Argument Table
!! \htmlinclude mp_nssl_2mom_run.html
!!
!>\ingroup nssl2mom
!>\section nssl2mom_run NSSL 2-moment microphysics
!>@{
      subroutine mp_nssl_2mom_run(ncol, nlev, kdt, con_g, con_rd,        &
                              spechum, qc, qr, qi, qs, qg, qh,           &
                              ni, nr, nc, ns, ng, nh, nqvolg, nqvolh,    &
                              tgrs, prsl, prslk, phii, omega, dtp,       &
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
         integer,                   intent(in   ) :: kdt
         real(kind_phys),           intent(in   ) :: con_g
         real(kind_phys),           intent(in   ) :: con_rd
         ! Hydrometeors
         real(kind_phys),           intent(inout) :: spechum(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qc(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qr(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qi(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qs(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qg(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qh(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: ni(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nr(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nc(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: ns(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: ng(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nh(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nqvolg(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nqvolh(1:ncol,1:nlev)
         ! State variables and timestep information
         real(kind_phys),           intent(inout) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prslk(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: phii(1:ncol,1:nlev+1)
         real(kind_phys),           intent(in   ) :: omega(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: dtp
         ! Precip/rain/snow/graupel fall amounts and fraction of frozen precip
         real(kind_phys),           intent(inout) :: prcp(1:ncol)
         real(kind_phys),           intent(inout) :: rain(1:ncol)
         real(kind_phys),           intent(inout) :: graupel(1:ncol)
         real(kind_phys),           intent(inout) :: ice(1:ncol)
         real(kind_phys),           intent(inout) :: snow(1:ncol)
         real(kind_phys),           intent(inout) :: sr(1:ncol)
         ! Radar reflectivity
         real(kind_phys),           intent(inout) :: refl_10cm(1:ncol,1:nlev)
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
         real(kind_phys) :: rho(1:ncol,1:nlev)              ! kg m-3
         ! Hydrometeors
         real(kind_phys) :: qv_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qc_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qr_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qi_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qs_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qg_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qh_mp(1:ncol,1:nlev)            ! kg kg-1 (dry mixing ratio)
         ! Rain/snow/graupel fall amounts
         real(kind_phys) :: delta_rain_mp(1:ncol)           ! mm
         real(kind_phys) :: delta_graupel_mp(1:ncol)        ! mm
         real(kind_phys) :: delta_ice_mp(1:ncol)            ! mm
         real(kind_phys) :: delta_snow_mp(1:ncol)           ! mm
         real(kind_phys) :: delta_hail_mp(1:ncol)           ! mm
         ! Vertical velocity and level width
         real(kind_phys) :: w(1:ncol,1:nlev)                ! m s-1
         real(kind_phys) :: dz(1:ncol,1:nlev)               ! m
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

         ! Dimensions used in mp_gt_driver
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte

         real(kind_phys) :: hail(1:ncol)
          ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Check initialization state
         if (.not.is_initialized) then
            write(errmsg, fmt='((a))') 'mp_nssl_2mom_run called before mp_nssl_2mom_init'
            errflg = 1
            return
         end if

         ! Density of air in kg m-3
         rho = prsl/(con_rd*tgrs)

         ! Convert omega in Pa s-1 to vertical velocity w in m s-1
         w = -omega/(rho*con_g)

         ! Layer width in m from geopotential in m2 s-2
         dz = (phii(:,2:nlev+1) - phii(:,1:nlev)) / con_g

         ! Flags for calculating radar reflectivity; diagflag is redundant
         if (do_radar_ref) then
             diagflag = .true.
         else
             diagflag = .false.
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
             write(errmsg,fmt='(*(a))') 'Logic error in mp_nssl_2mom_run:',  &
                                        ' all or none of the following optional', &
                                        ' arguments are required: re_cloud, re_ice, re_snow'
             errflg = 1
             return
         end if

         ! Convert specific humidity/moist mixing ratios to dry mixing ratios
         qv_mp = spechum/(1.0_kind_phys-spechum)
         qc_mp = qc/(1.0_kind_phys-spechum)
         qr_mp = qr/(1.0_kind_phys-spechum)
         qi_mp = qi/(1.0_kind_phys-spechum)
         qs_mp = qs/(1.0_kind_phys-spechum)
         qg_mp = qg/(1.0_kind_phys-spechum)
         qh_mp = qh/(1.0_kind_phys-spechum)

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
 
         call  nssl_2mom_driver(       &
                     ITIMESTEP=kdt + 1,    & ! TAS: XXX this is dumb
                     TH=tgrs,          &
                     QV=qv_mp,         &
                     QC=qc_mp,         &
                     QR=qr_mp,         &
                     QI=qi_mp,         &
                     QS=qs_mp,         &
                     QH=qg_mp,         &
                     QHL=qh_mp,        &
                     CCW=nc,           &
                     CRW=nr,           &
                     CCI=ni,           &
                     CSW=ns,           &
                     CHW=ng,           &
                     CHL=nh,           &
                     VHW=nqvolg,       &
                     VHL=nqvolh,       &
                     PII=prslk,                       &
                     P=prsl,                          &
                     W=w,                             &
                     DZ=dz,                           &
                     DTP=dtp,                         &
                     DN=rho,                          &
                     RAINNCV  = delta_rain_mp,        &
                     ICENCV   = delta_ice_mp,         &
                     SNOWNCV  = delta_snow_mp,        &
                     GRPLNCV  = delta_graupel_mp,     &
                     HAILNCV  = delta_hail_mp,        &
                     SR       = sr,                   &
                     dbz      = refl_10cm,            &
                     nssl_progn= .false.,             &
                     diagflag  = diagflag,            &
                     re_cloud = re_cloud_mp,          &
                     re_ice   = re_ice_mp,            &
                     re_snow  = re_snow_mp,           &
                     has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                     errmsg=errmsg, errflg=errflg,                                  &
                     ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,          &
                     ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,          &
                     its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte  )
         
         if (errflg/=0) return

         ! convert dry mixing ratios to specific humidity/moist mixing ratios
         spechum = qv_mp/(1.0_kind_phys+qv_mp)
         qc      = qc_mp/(1.0_kind_phys+qv_mp)
         qr      = qr_mp/(1.0_kind_phys+qv_mp)
         qi      = qi_mp/(1.0_kind_phys+qv_mp)
         qs      = qs_mp/(1.0_kind_phys+qv_mp)
         qg      = qg_mp/(1.0_kind_phys+qv_mp)         
         qh      = qh_mp/(1.0_kind_phys+qv_mp)         

         ! Convert mm output from microphysics to m expected by CCPP
         prcp    = max(0.0, delta_rain_mp/1000.0_kind_phys)
         graupel = max(0.0, (delta_hail_mp + delta_graupel_mp)/1000.0_kind_phys)
         ice     = max(0.0, delta_ice_mp/1000.0_kind_phys)
         snow    = max(0.0, delta_snow_mp/1000.0_kind_phys)
         rain    = max(0.0, (delta_rain_mp - (delta_graupel_mp + delta_ice_mp + delta_snow_mp + delta_hail_mp))/1000.0_kind_phys)

         if (do_effective_radii) then
            re_cloud = re_cloud_mp
            re_ice   = re_ice_mp
            re_snow  = re_snow_mp
         end if

        end subroutine mp_nssl_2mom_run
!>@}

end module mp_nssl_2mom

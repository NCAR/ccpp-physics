!> \file gfdl_cloud_microphys.F90
!! This file contains the CCPP entry point for the column GFDL cloud microphysics ( Chen and Lin (2013)
!! \cite chen_and_lin_2013 ).

!> This module contains the CCPP entry point for the column GFDL cloud microphysics ( Chen and Lin (2013)
!! \cite chen_and_lin_2013 ).
module gfdl_cloud_microphys

   use gfdl_cloud_microphys_mod, only: gfdl_cloud_microphys_mod_init,   &
                                       gfdl_cloud_microphys_mod_driver, &
                                       gfdl_cloud_microphys_mod_end,    &
                                       cloud_diagnosis

   implicit none

   private

   public gfdl_cloud_microphys_run, gfdl_cloud_microphys_init, gfdl_cloud_microphys_finalize

   logical :: is_initialized = .false.

contains

! -----------------------------------------------------------------------
! CCPP entry points for gfdl cloud microphysics
! -----------------------------------------------------------------------

!>\brief The subroutine initializes the GFDL
!! cloud microphysics.
!!
!> \section arg_table_gfdl_cloud_microphys_init Argument Table
!! \htmlinclude gfdl_cloud_microphys_init.html
!!
   subroutine gfdl_cloud_microphys_init (me, master, nlunit, input_nml_file, logunit, fn_nml, &
                                         imp_physics, imp_physics_gfdl, do_shoc, errmsg, errflg)

       implicit none

       integer, intent (in) :: me
       integer, intent (in) :: master
       integer, intent (in) :: nlunit
       integer, intent (in) :: logunit
       character(len=*), intent (in) :: fn_nml
       character(len=*), intent (in) :: input_nml_file(:)
       integer,          intent( in) :: imp_physics
       integer,          intent( in) :: imp_physics_gfdl
       logical,          intent( in) :: do_shoc
       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (is_initialized) return

       if (imp_physics/=imp_physics_gfdl) then
          write(errmsg,'(*(a))') 'Namelist option for microphysics does not match choice in suite definition file'
          errflg = 1
          return
       end if

       if (do_shoc) then
           write(errmsg,'(*(a))') 'SHOC is not currently compatible with GFDL MP'
           errflg = 1
           return
       endif

       call gfdl_cloud_microphys_mod_init(me, master, nlunit, input_nml_file, logunit, fn_nml, errmsg, errflg)

       is_initialized = .true.

   end subroutine gfdl_cloud_microphys_init

! =======================================================================
!>\brief The subroutine 'gfdl_cloud_microphys_finalize' terminates the GFDL
!! cloud microphysics.
!!
!! \section arg_table_gfdl_cloud_microphys_finalize  Argument Table
!! \htmlinclude gfdl_cloud_microphys_finalize.html
!!
   subroutine gfdl_cloud_microphys_finalize(errmsg, errflg)

       implicit none

       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg

   ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (.not.is_initialized) return

       call gfdl_cloud_microphys_mod_end()

       is_initialized = .false.

   end subroutine gfdl_cloud_microphys_finalize

!>\defgroup gfdlmp  GFDL Cloud Microphysics Module
!! This is cloud microphysics package for GFDL global cloud resolving model.
!! The algorithms are originally derived from Lin et al. (1983) \cite lin_et_al_1983.
!! Most of the key elements have been simplified/improved. This code at this stage
!! bears little to no similarity to the original Lin MP.
!! Therefore, it is best to be called GFDL microphysics (GFDL MP) .
!!
!>\brief The module contains the GFDL cloud
!! microphysics (Chen and Lin (2013) \cite chen_and_lin_2013 ).
!> The module is paired with \ref fast_sat_adj, which performs the "fast"
!! processes.
!!
!>\brief The subroutine executes the full GFDL cloud microphysics.
!! \section arg_table_gfdl_cloud_microphys_run Argument Table
!! \htmlinclude gfdl_cloud_microphys_run.html
!!
   subroutine gfdl_cloud_microphys_run(                                            &
      levs, im, rainmin, con_g, con_fvirt, con_rd, con_eps, frland, garea, islmsk, &
      gq0, gq0_ntcw, gq0_ntrw, gq0_ntiw, gq0_ntsw, gq0_ntgl, gq0_ntclamt,          &
      gt0, gu0, gv0, vvl, prsl, phii, del,                                         &
      rain0, ice0, snow0, graupel0, prcp0, sr,                                     &
      dtp, hydrostatic, phys_hydrostatic, lradar, refl_10cm,                       &
      reset, effr_in, rew, rei, rer, res, reg,                                     &
      cplchm, pfi_lsan, pfl_lsan, ten_t, ten_u, ten_v, ten_qv, ten_ql, ten_qr,     &
      ten_qi, ten_qs, ten_qg, ten_cldfrc, ten_q, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! DH* TODO: CLEANUP, all of these should be coming in through the argument list
      ! parameters
      real(kind=kind_phys), parameter :: one = 1.0d0
      real(kind=kind_phys), parameter :: con_p001= 0.001d0
      real(kind=kind_phys), parameter :: con_day = 86400.d0
      !real(kind=kind_phys), parameter :: rainmin = 1.0d-13
      ! *DH

      ! interface variables
      integer,              intent(in   ) :: levs, im
      real(kind=kind_phys), intent(in   ) :: con_g, con_fvirt, con_rd, con_eps, rainmin
      real(kind=kind_phys), intent(in   ), dimension(:)     :: frland, garea
      integer,              intent(in   ), dimension(:)     :: islmsk
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: gq0, gq0_ntcw, gq0_ntrw, gq0_ntiw, &
                                                               gq0_ntsw, gq0_ntgl, gq0_ntclamt
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: gt0, gu0, gv0
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: vvl, prsl, del
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: phii

      ! rain/snow/ice/graupel/precip amounts, fraction of frozen precip
      real(kind_phys),      intent(out  ), dimension(:) :: rain0
      real(kind_phys),      intent(out  ), dimension(:) :: snow0
      real(kind_phys),      intent(out  ), dimension(:) :: ice0
      real(kind_phys),      intent(out  ), dimension(:) :: graupel0
      real(kind_phys),      intent(out  ), dimension(:) :: prcp0
      real(kind_phys),      intent(out  ), dimension(:) :: sr
      
      real(kind_phys),      intent(out  ), dimension(:,:)    :: ten_t, ten_u, ten_v, ten_qv, ten_ql, ten_qr, ten_qi, ten_qs, ten_qg, ten_cldfrc
      real(kind_phys),      intent(out  ), dimension(:,:,:)  :: ten_q 

      real(kind_phys),      intent(in) :: dtp ! physics time step
      logical, intent (in) :: hydrostatic, phys_hydrostatic

      logical, intent (in) :: lradar
      real(kind=kind_phys), intent(inout), dimension(:,:) :: refl_10cm
      logical, intent (in) :: reset, effr_in
      real(kind=kind_phys), intent(inout), dimension(:,:), optional :: rew, rei, rer, res, reg
      logical, intent (in) :: cplchm
      ! ice and liquid water 3d precipitation fluxes - only allocated if cplchm is .true.
      real(kind=kind_phys), intent(inout), dimension(:,:), optional :: pfi_lsan, pfl_lsan

      character(len=*), intent(out) :: errmsg
      integer, intent(out)          :: errflg

      ! local variables
      integer :: iis, iie, jjs, jje, kks, kke, kbot, ktop
      integer :: i, k, kk
      real(kind=kind_phys), dimension(1:im,1:levs) :: delp, dz, uin, vin, pt, qv1, ql1, qr1, qg1, qa1, qn1, qi1,    &
                                                      qs1, pt_dt, qa_dt, u_dt, v_dt, w, qv_dt, ql_dt, qr_dt, qi_dt, &
                                                      qs_dt, qg_dt, p123, refl, new_qv, new_ql, new_qi, new_qr,     &
                                                      new_qs, new_qg, new_t 
      real(kind=kind_phys), dimension(1:im,1,1:levs) :: pfils, pflls
      real(kind=kind_phys), dimension(:,:), allocatable :: den
      real(kind=kind_phys) :: onebg
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      ten_t = 0.0
      ten_u = 0.0
      ten_v = 0.0
      ten_q = 0.0  !set tendency of entire tracer array to zero to make sure that those tracers not affected by this scheme do not change when tendencies are applied
      ten_qv = 0.0
      ten_ql = 0.0
      ten_qr = 0.0
      ten_qi = 0.0
      ten_qs = 0.0
      ten_qg = 0.0
      ten_cldfrc = 0.0
      
      iis = 1
      iie = im
      jjs = 1
      jje = 1
      kks = 1
      kke = levs
      ! flipping of vertical direction
      ktop = 1
      kbot = levs

      onebg = one/con_g

      do k = 1, levs
         kk = levs-k+1
         do i = 1, im
            qv_dt(i,k) = 0.0
            ql_dt(i,k) = 0.0
            qr_dt(i,k) = 0.0
            qi_dt(i,k) = 0.0
            qs_dt(i,k) = 0.0
            qg_dt(i,k) = 0.0
            qa_dt(i,k) = 0.0
            pt_dt(i,k) = 0.0
            u_dt(i,k)  = 0.0
            v_dt(i,k)  = 0.0
            qn1(i,k)   = 0.0
            pfils(i,1,k) = 0.0
            pflls(i,1,k) = 0.0
            ! flip vertical (k) coordinate
            qv1(i,k)  = gq0(i,kk)
            ql1(i,k)  = gq0_ntcw(i,kk)
            qr1(i,k)  = gq0_ntrw(i,kk)
            qi1(i,k)  = gq0_ntiw(i,kk)
            qs1(i,k)  = gq0_ntsw(i,kk)
            qg1(i,k)  = gq0_ntgl(i,kk)
            qa1(i,k)  = gq0_ntclamt(i,kk)
            pt(i,k)   = gt0(i,kk)
            w(i,k)    = -vvl(i,kk) * (one+con_fvirt * gq0(i,kk))   &
                          *  gt0(i,kk) / prsl(i,kk) * (con_rd*onebg)
            uin(i,k)  = gu0(i,kk)
            vin(i,k)  = gv0(i,kk)
            delp(i,k) = del(i,kk)
            dz(i,k)   = (phii(i,kk)-phii(i,kk+1))*onebg
            p123(i,k) = prsl(i,kk)
            refl(i,k) = refl_10cm(i,kk)
         enddo
      enddo

      ! reset precipitation amounts to zero
      rain0     = 0
      ice0      = 0
      snow0     = 0
      graupel0  = 0

      call gfdl_cloud_microphys_mod_driver(iis, iie, jjs, jje, kks, kke, ktop, kbot, &
                 qv1, ql1, qr1, qi1, qs1, qg1, qa1, qn1, qv_dt, ql_dt, qr_dt, qi_dt, &
                 qs_dt, qg_dt, qa_dt, pt_dt, pt, w,  uin, vin, u_dt, v_dt, dz, delp, &
                 garea, dtp, frland, rain0, snow0, ice0, graupel0, hydrostatic,      &
                 phys_hydrostatic, p123, lradar, refl, reset, pfils, pflls)
      tem   = dtp*con_p001/con_day

      ! fix negative values
      do i = 1, im
        !rain0(i)    = max(con_d00, rain0(i))
        !snow0(i)    = max(con_d00, snow0(i))
        !ice0(i)     = max(con_d00, ice0(i))
        !graupel0(i) = max(con_d00, graupel0(i))
        if(rain0(i)*tem < rainmin) then
          rain0(i) = 0.0
        endif
        if(ice0(i)*tem < rainmin) then
          ice0(i) = 0.0
        endif
        if(snow0(i)*tem < rainmin) then
          snow0(i) = 0.0
        endif
        if(graupel0(i)*tem < rainmin) then
          graupel0(i) = 0.0
        endif
      enddo

      ! calculate fraction of frozen precipitation using unscaled
      ! values of rain0, ice0, snow0, graupel0 (for bit-for-bit)
      do i=1,im
        prcp0(i) = (rain0(i)+snow0(i)+ice0(i)+graupel0(i)) * tem
        if ( prcp0(i) > rainmin ) then
          sr(i) = (snow0(i) + ice0(i)  + graupel0(i)) &
                      / (rain0(i) + snow0(i) + ice0(i) + graupel0(i))
        else
          sr(i) = 0.0
        endif
      enddo

      ! convert rain0, ice0, snow0, graupel0 from mm per day to m per physics timestep
      rain0    = rain0*tem
      ice0     = ice0*tem
      snow0    = snow0*tem
      graupel0 = graupel0*tem

      ! flip vertical coordinate back
      do k=1,levs
        kk = levs-k+1
        do i=1,im
            ten_qv(i,k)      = qv_dt(i,kk)
            ten_ql(i,k)      = ql_dt(i,kk)
            ten_qr(i,k)      = qr_dt(i,kk)
            ten_qi(i,k)      = qi_dt(i,kk)
            ten_qs(i,k)      = qs_dt(i,kk)
            ten_qg(i,k)      = qg_dt(i,kk)
            ten_cldfrc(i,k)  = qa_dt(i,kk)
            
            ten_t(i,k)       = pt_dt(i,kk)
            ten_u(i,k)       = u_dt(i,kk)
            ten_v(i,k)       = v_dt(i,kk)
            
            refl_10cm(i,k)   = refl(i,kk)
            
            !new values needed for cloud_diagnosis below
            new_qv(i,k)      = qv1(i,kk) + qv_dt(i,kk) * dtp
            new_ql(i,k)      = ql1(i,kk) + ql_dt(i,kk) * dtp
            new_qr(i,k)      = qr1(i,kk) + qr_dt(i,kk) * dtp
            new_qi(i,k)      = qi1(i,kk) + qi_dt(i,kk) * dtp
            new_qs(i,k)      = qs1(i,kk) + qs_dt(i,kk) * dtp
            new_qg(i,k)      = qg1(i,kk) + qg_dt(i,kk) * dtp
            new_t(i,k)       = gt0(i,k)  + pt_dt(i,kk) * dtp
        enddo
      enddo

      ! output ice and liquid water 3d precipitation fluxes if requested
      if (cplchm) then
        do k=1,levs
          kk = levs-k+1
          do i=1,im
            pfi_lsan(i,k) = pfils(i,1,kk)
            pfl_lsan(i,k) = pflls(i,1,kk)
          enddo
        enddo
      endif
            
      if(effr_in) then
         allocate(den(1:im,1:levs))
         do k=1,levs
            do i=1,im
               den(i,k)=con_eps*prsl(i,k)/(con_rd*new_t(i,k)*(new_qv(i,k)+con_eps))
            enddo
         enddo
         call cloud_diagnosis (1, im, 1, levs, den(1:im,1:levs), &
            del(1:im,1:levs),      islmsk(1:im),                 &
            new_ql(1:im,1:levs), new_qi(1:im,1:levs),            &
            new_qr(1:im,1:levs),                                 &
            new_qs(1:im,1:levs) + new_qg(1:im,1:levs),           &
            new_qg(1:im,1:levs)*0.0, new_t(1:im,1:levs),         &
            rew(1:im,1:levs), rei(1:im,1:levs), rer(1:im,1:levs),&
            res(1:im,1:levs), reg(1:im,1:levs))
         deallocate(den)
      endif

   end subroutine gfdl_cloud_microphys_run

end module gfdl_cloud_microphys

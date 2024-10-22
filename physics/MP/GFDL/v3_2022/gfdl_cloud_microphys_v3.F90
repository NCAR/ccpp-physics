!> \file gfdl_cloud_microphys_v3.F90
!! This file contains the CCPP entry point for the column GFDL cloud microphysics version 3 ( Chen and Lin (2013)
!! \cite chen_and_lin_2013 ).
module gfdl_cloud_microphys_v3

   use module_gfdl_cloud_microphys_v3, only: module_gfdl_cloud_microphys_v3_init,         &
                                             module_gfdl_cloud_microphys_v3_driver,       &
                                             module_gfdl_cloud_microphys_v3_end,          &
                                             rad_ref, cld_eff_rad    

   implicit none

   private

   public gfdl_cloud_microphys_v3_run, gfdl_cloud_microphys_v3_init, gfdl_cloud_microphys_v3_finalize

   logical :: is_initialized = .false.

contains

! -----------------------------------------------------------------------
! CCPP entry points for gfdl cloud microphysics
! -----------------------------------------------------------------------

!>\brief The subroutine initializes the GFDL
!! cloud microphysics.
!!
!> \section arg_table_gfdl_cloud_microphys_v3_init Argument Table
!! \htmlinclude gfdl_cloud_microphys_v3_init.html
!!

   subroutine gfdl_cloud_microphys_v3_init (me, master, nlunit, input_nml_file, logunit, &
                            fn_nml, imp_physics, imp_physics_gfdl, do_shoc,  &
                            hydrostatic, errmsg, errflg)

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
       logical,          intent( in) :: hydrostatic
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
           write(errmsg,'(*(a))') 'SHOC is not currently compatible with GFDL MP v3'
           errflg = 1
           return
       endif

       call module_gfdl_cloud_microphys_v3_init(me, master, nlunit, input_nml_file, logunit, fn_nml, hydrostatic, errmsg, errflg)

       is_initialized = .true.

   end subroutine gfdl_cloud_microphys_v3_init


! =======================================================================
!>\brief The subroutine 'gfdl_cloud_microphys_v3_finalize' terminates the GFDL
!! cloud microphysics.
!!
!! \section arg_table_gfdl_cloud_microphys_v3_finalize  Argument Table
!! \htmlinclude gfdl_cloud_microphys_v3_finalize.html
!!
   subroutine gfdl_cloud_microphys_v3_finalize(errmsg, errflg)

       implicit none

       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg

   ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (.not.is_initialized) return

       call module_gfdl_cloud_microphys_v3_end()

       is_initialized = .false.

   end subroutine gfdl_cloud_microphys_v3_finalize

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
!! \section arg_table_gfdl_cloud_microphys_v3_run Argument Table
!! \htmlinclude gfdl_cloud_microphys_v3_run.html
!!
   subroutine gfdl_cloud_microphys_v3_run(                                         &
      imp_physics, imp_physics_gfdl, fast_mp_consv,                                &
      levs, im, rainmin, con_g, con_fvirt, con_rd, con_eps, garea, slmsk, snowd,   &
      gq0, gq0_ntcw, gq0_ntrw, gq0_ntiw, gq0_ntsw, gq0_ntgl, gq0_ntclamt, aerfld,  &
      gt0, gu0, gv0, vvl, prsl, phii, del,                                         &
      rain0, ice0, snow0, graupel0, prcp0, sr, oro,                                &
      dtp, hydrostatic, lradar, refl_10cm,                                         &
      reset, effr_in, rew, rei, rer, res, reg,                                     &
      cplchm, pfi_lsan, pfl_lsan, errmsg, errflg)

      use machine, only: kind_phys, kind_dyn 

      implicit none

      ! DH* TODO: CLEANUP, all of these should be coming in through the argument list
      ! parameters
      real(kind=kind_phys), parameter :: one = 1.0d0
      real(kind=kind_phys), parameter :: con_p001= 0.001d0
      real(kind=kind_phys), parameter :: con_day = 86400.d0
      !real(kind=kind_phys), parameter :: rainmin = 1.0d-13
      ! *DH

      ! interface variables
      integer,              intent(in   ) :: imp_physics
      integer,              intent(in   ) :: imp_physics_gfdl
      integer,              intent(in   ) :: levs, im
      real(kind=kind_phys), intent(in   ) :: con_g, con_fvirt, con_rd, con_eps, rainmin
      real(kind=kind_phys), intent(in   ), dimension(:)     :: garea, slmsk, snowd, oro 
      real(kind=kind_phys), intent(inout), dimension(:,:)   :: gq0, gq0_ntcw, gq0_ntrw, gq0_ntiw, &
                                                               gq0_ntsw, gq0_ntgl, gq0_ntclamt
      real(kind_phys),      intent(in   ), dimension(:,:,:) :: aerfld
      real(kind=kind_phys), intent(inout), dimension(:,:)   :: gt0, gu0, gv0
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: vvl, prsl, del
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: phii

      ! rain/snow/ice/graupel/precip amounts, fraction of frozen precip
      !real(kind_phys),                     dimension(:) :: water0
      real(kind_phys),      intent(out  ), dimension(:), optional :: rain0
      real(kind_phys),      intent(out  ), dimension(:), optional :: snow0
      real(kind_phys),      intent(out  ), dimension(:), optional :: ice0
      real(kind_phys),      intent(out  ), dimension(:), optional :: graupel0
      real(kind_phys),      intent(out  ), dimension(:) :: prcp0
      real(kind_phys),      intent(out  ), dimension(:) :: sr

      real(kind_phys),      intent(in) :: dtp ! physics time step
      logical, intent (in) :: hydrostatic, fast_mp_consv

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
      real(kind=kind_phys), dimension(1:im,1:levs) :: delp, dz, uin, vin, pt, qv1, ql1, qi1, qr1, qs1, qg1,    &  
                                                      qa1, qnl, qni, pt_dt, qa_dt, u_dt, v_dt, w, qv_dt, ql_dt,&
                                                      qr_dt, qi_dt, qs_dt, qg_dt, p123, refl
      real(kind=kind_phys), dimension(1:im,1:levs) :: q_con, cappa !for inline MP option  
      real(kind=kind_phys), dimension(1:im,1,1:levs) :: pfils, pflls
      real(kind=kind_phys), dimension(1:im,1,1:levs) :: adj_vmr, te
      real(kind=kind_phys), dimension(1:im,1:levs) :: prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
      real(kind=kind_phys), dimension(1:im) :: dte, hs, gsize  
      !real(kind=kind_phys), dimension(:,:), allocatable :: den
      real(kind=kind_phys), dimension(1:im) :: water0
      real(kind=kind_phys) :: onebg
      real(kind=kind_phys) :: tem
      logical last_step, do_inline_mp

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

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
            qnl(i,k)   = aerfld(i,kk,11) ! sulfate 
            pfils(i,1,k) = 0.0
            pflls(i,1,k) = 0.0
            prefluxw(i,k) =0.0 
            prefluxi(i,k) =0.0 
            prefluxr(i,k) =0.0 
            prefluxs(i,k) =0.0 
            prefluxg(i,k) =0.0 

            ! flip vertical (k) coordinate top =1 
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
            qni(i,k)  = 10.
            q_con(i,k) = 0.0 
            cappa(i,k) = 0.0 
         enddo
      enddo

      ! reset precipitation amounts to zero
      water0     = 0
      rain0     = 0
      ice0      = 0
      snow0     = 0
      graupel0  = 0
 
      if(imp_physics == imp_physics_gfdl) then 

        last_step = .false.
        do_inline_mp = .false. 
        hs = oro(:) * con_g 
        gsize = sqrt(garea(:)) 

        call module_gfdl_cloud_microphys_v3_driver( qv1, ql1, qr1, qi1, qs1, qg1, qa1, qnl, qni, pt, w,&
                  uin, vin, dz, delp, gsize, dtp, hs, water0, rain0,                        &
                  ice0, snow0, graupel0, hydrostatic, iis, iie, kks, kke, q_con, cappa,     &
                  fast_mp_consv, adj_vmr, te, dte, prefluxw, prefluxr, prefluxi, prefluxs,  &
                  prefluxg, last_step, do_inline_mp ) 

      else 
        write(errmsg,'(*(a))') 'Invalid imp_physics option for GFDL MP v3'
        errflg = 1
        return
      endif 
      tem   = dtp*con_p001/con_day

      ! fix negative values
      do i = 1, im
        !rain0(i)    = max(con_d00, rain0(i))
        !snow0(i)    = max(con_d00, snow0(i))
        !ice0(i)     = max(con_d00, ice0(i))
        !graupel0(i) = max(con_d00, graupel0(i))
        if(water0(i)*tem < rainmin) then
          water0(i) = 0.0
        endif
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
      water0    = water0*tem
      rain0    = rain0*tem
      ice0     = ice0*tem
      snow0    = snow0*tem
      graupel0 = graupel0*tem

      ! flip vertical coordinate back
      do k=1,levs
        kk = levs-k+1
        do i=1,im

          if (imp_physics == imp_physics_gfdl) then 
            gq0(i,k)         = qv1(i,kk)
            gq0_ntcw(i,k)    = ql1(i,kk) 
            gq0_ntrw(i,k)    = qr1(i,kk)
            gq0_ntiw(i,k)    = qi1(i,kk)
            gq0_ntsw(i,k)    = qs1(i,kk)
            gq0_ntgl(i,k)    = qg1(i,kk)
            gq0_ntclamt(i,k) = qa1(i,kk)
            gt0(i,k)         = pt(i,kk)  
            gu0(i,k)         = uin(i,kk) 
            gv0(i,k)         = vin(i,kk) 
            refl_10cm(i,k)   = refl(i,kk)
           else 
            write(errmsg,'(*(a))') 'Invalid imp_physics option for GFDL MP v3' 
            errflg = 1
            return
           endif 
        enddo
      enddo

      ! output ice and liquid water 3d precipitation fluxes if requested
      if (cplchm) then
        do k=1,levs
          kk = levs-k+1
          do i=1,im
            if (imp_physics==imp_physics_gfdl) then 
              pfi_lsan(i,k) = prefluxi (i,kk) + prefluxs (i,kk) + prefluxg (i,kk)
              pfl_lsan(i,k) = prefluxr (i,kk) 
            else 
              write(errmsg,'(*(a))') 'Invalid imp_physics option for GFDL MP v3'
              errflg = 1
              return
            endif 
          enddo
        enddo
      endif

      if(effr_in) then
         call cld_eff_rad (1, im, 1, levs, slmsk(1:im),          &
            prsl(1:im,1:levs),  del(1:im,1:levs),                &
            gt0(1:im,1:levs), gq0(1:im,1:levs),                  &
            gq0_ntcw(1:im,1:levs), gq0_ntiw(1:im,1:levs),        &
            gq0_ntrw(1:im,1:levs), gq0_ntsw(1:im,1:levs),        &
            gq0_ntgl(1:im,1:levs), gq0_ntclamt(1:im,1:levs),     &
            rew(1:im,1:levs), rei(1:im,1:levs), rer(1:im,1:levs),&
            res(1:im,1:levs), reg(1:im,1:levs),snowd(1:im))
      endif

      if(lradar) then 
         call rad_ref (1, im, 1, 1, qv1(1:im,1:levs), qr1(1:im,1:levs), &
	    qs1(1:im,1:levs),qg1(1:im,1:levs),pt(1:im,1:levs),          & 
            delp(1:im,1:levs), dz(1:im,1:levs), refl(1:im,1:levs), levs, hydrostatic,  & 
            do_inline_mp, 1)

         do k=1,levs
           kk = levs-k+1
           do i=1,im
              refl_10cm(i,k)   = max(-35.,refl(i,kk))
           enddo 
         enddo 
      endif 

   end subroutine gfdl_cloud_microphys_v3_run

end module gfdl_cloud_microphys_v3

!> \file GFS_PBL_generic.F90
!!  Contains code related to PBL schemes to be used within the GFS physics suite.

      module GFS_PBL_generic_common

      implicit none

      private

      public :: set_aerosol_tracer_index

      contains

      subroutine set_aerosol_tracer_index(imp_physics, imp_physics_wsm6,          &
                                          imp_physics_thompson, ltaerosol,        &
                                          imp_physics_mg, ntgl, imp_physics_gfdl, &
                                          imp_physics_zhao_carr, kk,              &
                                          errmsg, errflg)
      implicit none
      !
      integer, intent(in )          :: imp_physics, imp_physics_wsm6,          &
                                       imp_physics_thompson,                   &
                                       imp_physics_mg, ntgl, imp_physics_gfdl, &
                                       imp_physics_zhao_carr
      logical, intent(in )          :: ltaerosol
      integer, intent(out)          :: kk
      character(len=*), intent(out) :: errmsg
      integer, intent(out)          :: errflg

      errflg = 0

! Set Interstitial%kk = last index in diffused tracer array before chemistry-aerosol tracers
      if (imp_physics == imp_physics_wsm6) then
! WSM6
        kk = 4
      elseif (imp_physics == imp_physics_thompson) then
! Thompson
        if(ltaerosol) then
          kk = 10
        else
          kk = 7
        endif
! MG
      elseif (imp_physics == imp_physics_mg) then
        if (ntgl > 0) then
          kk = 12
        else
          kk = 10
        endif
      elseif (imp_physics == imp_physics_gfdl) then
! GFDL MP
        kk = 7
      elseif (imp_physics == imp_physics_zhao_carr) then
! Zhao/Carr/Sundqvist
        kk = 3
      else
        write(errmsg,'(*(a))') 'Logic error: unknown microphysics option in set_aerosol_tracer_index'
        kk = -999
        errflg = 1
        return
      endif

      end subroutine set_aerosol_tracer_index

      end module GFS_PBL_generic_common


      module GFS_PBL_generic_pre

      contains

      subroutine GFS_PBL_generic_pre_init ()
      end subroutine GFS_PBL_generic_pre_init

      subroutine GFS_PBL_generic_pre_finalize()
      end subroutine GFS_PBL_generic_pre_finalize

!> \brief This scheme sets up the vertically diffused tracer array for any PBL scheme based on the microphysics scheme chosen
!! \section arg_table_GFS_PBL_generic_pre_run Argument Table
!! \htmlinclude GFS_PBL_generic_pre_run.html
!!
      subroutine GFS_PBL_generic_pre_run (im, levs, nvdiff, ntrac, rtg_ozone_index,      &
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc,                 &
        ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef, trans_aero, ntchs, ntchm,          &
        imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,           &
        imp_physics_zhao_carr, imp_physics_mg, imp_physics_fer_hires, cplchm, ltaerosol, &
        hybedmf, do_shoc, satmedmf, qgrs, vdftra, save_u, save_v, save_t, save_q,        &
        flag_for_pbl_generic_tend, ldiag3d, qdiag3d, lssav, ugrs, vgrs, tgrs, errmsg, errflg)
        
      use machine,                only : kind_phys
      use GFS_PBL_generic_common, only : set_aerosol_tracer_index

      implicit none

      integer, parameter  :: kp = kind_phys
      integer, intent(out) :: rtg_ozone_index
      integer, intent(in) :: im, levs, nvdiff, ntrac
      integer, intent(in) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc
      integer, intent(in) :: ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef,ntchs, ntchm
      logical, intent(in) :: trans_aero, ldiag3d, qdiag3d, lssav
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      integer, intent(in) :: imp_physics_zhao_carr, imp_physics_mg, imp_physics_fer_hires
      logical, intent(in) :: cplchm, ltaerosol, hybedmf, do_shoc, satmedmf, flag_for_pbl_generic_tend

      real(kind=kind_phys), dimension(:,:,:), intent(in) :: qgrs
      real(kind=kind_phys), dimension(:,:), intent(in) :: ugrs, vgrs, tgrs
      real(kind=kind_phys), dimension(:,:, :), intent(inout) :: vdftra
      real(kind=kind_phys), dimension(:,:), intent(out) :: save_u, save_v, save_t
      real(kind=kind_phys), dimension(:,:, :), intent(out) :: save_q

      ! CCPP error handling variables
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real (kind=kind_phys), parameter :: zero = 0.0_kp, one=1.0_kp

      ! Local variables
      integer :: i, k, kk, k1, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      rtg_ozone_index=-1
!DH: dvdftra is only used if nvdiff != ntrac or (nvdiff == ntrac .and. )
      if (nvdiff == ntrac .and. (hybedmf .or. do_shoc .or. satmedmf)) then
        vdftra = qgrs
        rtg_ozone_index = ntoz
      else
        if (imp_physics == imp_physics_wsm6) then
  ! WSM6
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntiw)
              vdftra(i,k,4) = qgrs(i,k,ntoz)
            enddo
          enddo
          rtg_ozone_index = 4

  ! Ferrier-Aligo
        elseif (imp_physics == imp_physics_fer_hires) then
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntiw)
              vdftra(i,k,4) = qgrs(i,k,ntrw)
              vdftra(i,k,5) = qgrs(i,k,nqrimef)
              vdftra(i,k,6) = qgrs(i,k,ntoz)
            enddo
          enddo
          rtg_ozone_index = 6
        
        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          if(ltaerosol) then
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = qgrs(i,k,ntinc)
                vdftra(i,k,9)  = qgrs(i,k,ntrnc)
                vdftra(i,k,10) = qgrs(i,k,ntoz)
                vdftra(i,k,11) = qgrs(i,k,ntwa)
                vdftra(i,k,12) = qgrs(i,k,ntia)
              enddo
            enddo
            rtg_ozone_index = 10
          else
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = qgrs(i,k,ntqv)
                vdftra(i,k,2) = qgrs(i,k,ntcw)
                vdftra(i,k,3) = qgrs(i,k,ntiw)
                vdftra(i,k,4) = qgrs(i,k,ntrw)
                vdftra(i,k,5) = qgrs(i,k,ntsw)
                vdftra(i,k,6) = qgrs(i,k,ntgl)
                vdftra(i,k,7) = qgrs(i,k,ntinc)
                vdftra(i,k,8) = qgrs(i,k,ntrnc)
                vdftra(i,k,9) = qgrs(i,k,ntoz)
              enddo
            enddo
            rtg_ozone_index = 9
          endif
  ! MG
        elseif (imp_physics == imp_physics_mg) then        ! MG3/2
          if (ntgl > 0) then                               ! MG3
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = qgrs(i,k,ntinc)
                vdftra(i,k,9)  = qgrs(i,k,ntrnc)
                vdftra(i,k,10) = qgrs(i,k,ntsnc)
                vdftra(i,k,11) = qgrs(i,k,ntgnc)
                vdftra(i,k,12) = qgrs(i,k,ntoz)
              enddo
            enddo
            rtg_ozone_index = 12
          else                                             ! MG2
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntlnc)
                vdftra(i,k,7)  = qgrs(i,k,ntinc)
                vdftra(i,k,8)  = qgrs(i,k,ntrnc)
                vdftra(i,k,9)  = qgrs(i,k,ntsnc)
                vdftra(i,k,10) = qgrs(i,k,ntoz)
              enddo
            enddo
            rtg_ozone_index = 10
          endif
        elseif (imp_physics == imp_physics_gfdl) then
  ! GFDL MP
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntiw)
              vdftra(i,k,4) = qgrs(i,k,ntrw)
              vdftra(i,k,5) = qgrs(i,k,ntsw)
              vdftra(i,k,6) = qgrs(i,k,ntgl)
              vdftra(i,k,7) = qgrs(i,k,ntoz)
            enddo
          enddo
          rtg_ozone_index = 7
        elseif (imp_physics == imp_physics_zhao_carr) then
! Zhao/Carr/Sundqvist
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntoz)
            enddo
          enddo
          rtg_ozone_index = 3
        endif
!
        if (trans_aero) then
          call set_aerosol_tracer_index(imp_physics, imp_physics_wsm6,          &
                                        imp_physics_thompson, ltaerosol,        &
                                        imp_physics_mg, ntgl, imp_physics_gfdl, &
                                        imp_physics_zhao_carr, kk,              &
                                        errmsg, errflg)
          if (.not.errflg==1) return
          !
          k1 = kk
          do n=ntchs,ntchm+ntchs-1
            k1 = k1 + 1
            do k=1,levs
              do i=1,im
                vdftra(i,k,k1) = qgrs(i,k,n)
              enddo
            enddo
          enddo
        endif
!
        if (ntke>0) then
          do k=1,levs
            do i=1,im
              vdftra(i,k,ntkev) = qgrs(i,k,ntke)
            enddo
          enddo
        endif
!
      endif

      if(ldiag3d .and. lssav .and. flag_for_pbl_generic_tend) then
        do k=1,levs
          do i=1,im
            save_t(i,k) = tgrs(i,k)
            save_u(i,k) = ugrs(i,k)
            save_v(i,k) = vgrs(i,k)
          enddo
        enddo
        if(qdiag3d) then
          do k=1,levs
            do i=1,im
              save_q(i,k,ntqv) = qgrs(i,k,ntqv)
              save_q(i,k,ntoz) = qgrs(i,k,ntoz)
            enddo
          enddo
          if(ntke>0) then
            do k=1,levs
              do i=1,im
                save_q(i,k,ntke) = qgrs(i,k,ntke)
              enddo
            enddo
          endif
        endif
      endif

    end subroutine GFS_PBL_generic_pre_run

    end module GFS_PBL_generic_pre


    module GFS_PBL_generic_post

    contains

    subroutine GFS_PBL_generic_post_init ()
    end subroutine GFS_PBL_generic_post_init

    subroutine GFS_PBL_generic_post_finalize ()
    end subroutine GFS_PBL_generic_post_finalize

!> \section arg_table_GFS_PBL_generic_post_run Argument Table
!! \htmlinclude GFS_PBL_generic_post_run.html
!!
      subroutine GFS_PBL_generic_post_run (im, levs, nvdiff, ntrac,                                                            &
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc, ntwa, ntia, ntgl, ntoz, ntke, ntkev,nqrimef,          &
        trans_aero, ntchs, ntchm,                                                                                              &
        imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6, imp_physics_zhao_carr, imp_physics_mg,          &
        imp_physics_fer_hires,                                                                                                 &
        ltaerosol, cplflx, cplchm, lssav, flag_for_pbl_generic_tend, ldiag3d, lsidea, hybedmf, do_shoc, satmedmf,              &
        shinhong, do_ysu, dvdftra, dusfc1, dvsfc1, dtsfc1, dqsfc1, dtf, dudt, dvdt, dtdt, htrsw, htrlw, xmu,                   &
        dqdt, dusfc_cpl, dvsfc_cpl, dtsfc_cpl, dtend, dtidx, index_of_temperature, index_of_x_wind, index_of_y_wind,           &
        index_of_process_pbl, dqsfc_cpl, dusfci_cpl, dvsfci_cpl, dtsfci_cpl, dqsfci_cpl, dusfc_diag, dvsfc_diag, dtsfc_diag,   &
        dqsfc_diag, dusfci_diag, dvsfci_diag, dtsfci_diag, dqsfci_diag,                                                        &
        rd, cp, fvirt, hvap, t1, q1, prsl, hflx, ushfsfci, oceanfrac, kdt, dusfc_cice, dvsfc_cice,                             &
        dtsfc_cice, dqsfc_cice, wet, dry, icy, wind, stress_wat, hflx_wat, evap_wat, ugrs1, vgrs1, dkt_cpl, dkt, hffac, hefac, &
        ugrs, vgrs, tgrs, qgrs, save_u, save_v, save_t, save_q, errmsg, errflg)

      use machine,                only : kind_phys
      use GFS_PBL_generic_common, only : set_aerosol_tracer_index

      implicit none

      integer, parameter  :: kp = kind_phys
      integer, intent(in) :: im, levs, nvdiff, ntrac, ntchs, ntchm
      integer, intent(in) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc, ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef
      logical, intent(in) :: trans_aero
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      integer, intent(in) :: imp_physics_zhao_carr, imp_physics_mg, imp_physics_fer_hires
      logical, intent(in) :: ltaerosol, cplflx, cplchm, lssav, ldiag3d, lsidea
      logical, intent(in) :: hybedmf, do_shoc, satmedmf, shinhong, do_ysu
      integer, intent(in) :: kdt

      logical, intent(in) :: flag_for_pbl_generic_tend      
      real(kind=kind_phys), dimension(:,:), intent(in) :: save_u, save_v, save_t
      real(kind=kind_phys), dimension(:,:, :), intent(in) :: save_q

      real(kind=kind_phys), intent(in) :: dtf
      real(kind=kind_phys), intent(in) :: rd, cp, fvirt, hvap
      real(kind=kind_phys), dimension(:), intent(in) :: t1, q1, hflx, oceanfrac
      real(kind=kind_phys), dimension(:,:), intent(in) :: prsl
      real(kind=kind_phys), dimension(:), intent(in) :: dusfc_cice, dvsfc_cice, dtsfc_cice, dqsfc_cice, &
          wind, stress_wat, hflx_wat, evap_wat, ugrs1, vgrs1

      real(kind=kind_phys), dimension(:,:, :), intent(in) :: qgrs
      real(kind=kind_phys), dimension(:,:), intent(in) :: ugrs, vgrs, tgrs

      real(kind=kind_phys), dimension(:,:, :), intent(in) :: dvdftra
      real(kind=kind_phys), dimension(:), intent(in) :: dusfc1, dvsfc1, dtsfc1, dqsfc1, xmu
      real(kind=kind_phys), dimension(:,:), intent(in) :: dudt, dvdt, dtdt, htrsw, htrlw

      real(kind=kind_phys), dimension(:,:, :), intent(inout) :: dqdt

      ! The following arrays may not be allocated, depending on certain flags (cplflx, ...).
      ! Since Intel 15 crashes when passing unallocated arrays to arrays defined with explicit shape,
      ! use assumed-shape arrays. Note that Intel 18 and GNU 6.2.0-8.1.0 tolerate explicit-shape arrays
      ! as long as these do not get used when not allocated.
      real(kind=kind_phys), dimension(:),   intent(inout) :: dusfc_cpl, dvsfc_cpl, dtsfc_cpl, dqsfc_cpl, dusfci_cpl, dvsfci_cpl, &
        dtsfci_cpl, dqsfci_cpl, dusfc_diag, dvsfc_diag, dtsfc_diag, dqsfc_diag, dusfci_diag, dvsfci_diag, dtsfci_diag, dqsfci_diag
      real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:)
      integer, intent(in) :: index_of_temperature, index_of_x_wind, index_of_y_wind, index_of_process_pbl

      logical, dimension(:),intent(in) :: wet, dry, icy
      real(kind=kind_phys), dimension(:), intent(out) ::  ushfsfci

      real(kind=kind_phys), dimension(:,:), intent(inout) :: dkt_cpl
      real(kind=kind_phys), dimension(:,:), intent(in)    :: dkt

      ! From canopy heat storage - reduction factors in latent/sensible heat flux due to surface roughness
      real(kind=kind_phys), dimension(:), intent(in) :: hffac, hefac

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      real(kind=kind_phys), parameter :: zero  = 0.0_kp, one = 1.0_kp
      real(kind=kind_phys), parameter :: huge  = 9.9692099683868690E36 ! NetCDF float FillValue, same as in GFS_typedefs.F90
      real(kind=kind_phys), parameter :: qmin  = 1.0e-8_kp
      integer :: i, k, kk, k1, n
      real(kind=kind_phys) :: tem, rho
      integer :: idtend
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!GJF: dvdftra is only used if nvdiff != ntrac or (nvdiff == ntrac .and. )
      if (nvdiff == ntrac .and. (hybedmf .or. do_shoc .or. satmedmf)) then
        dqdt = dvdftra
      elseif (nvdiff /= ntrac .and. .not. shinhong .and. .not. do_ysu) then
!
        if (ntke>0) then
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntke)  = dvdftra(i,k,ntkev)
            enddo
          enddo
        endif
!
        if (trans_aero) then
          ! Set kk if chemistry-aerosol tracers are diffused
          call set_aerosol_tracer_index(imp_physics, imp_physics_wsm6,          &
                                        imp_physics_thompson, ltaerosol,        &
                                        imp_physics_mg, ntgl, imp_physics_gfdl, &
                                        imp_physics_zhao_carr, kk,              &
                                        errmsg, errflg)
          if (.not.errflg==1) return
          !
          k1 = kk
          do n=ntchs,ntchm+ntchs-1
            k1 = k1 + 1
            do k=1,levs
              do i=1,im
                dqdt(i,k,n) = dvdftra(i,k,k1)
              enddo
            enddo
          enddo
        endif
!
        if (imp_physics == imp_physics_wsm6) then
  ! WSM6
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntqv)  = dvdftra(i,k,1)
              dqdt(i,k,ntcw)  = dvdftra(i,k,2)
              dqdt(i,k,ntiw)  = dvdftra(i,k,3)
              dqdt(i,k,ntoz)  = dvdftra(i,k,4)
            enddo
          enddo

        elseif (imp_physics == imp_physics_fer_hires) then
  ! Ferrier-Aligo 
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntqv)    = dvdftra(i,k,1)
              dqdt(i,k,ntcw)    = dvdftra(i,k,2)
              dqdt(i,k,ntiw)    = dvdftra(i,k,3)
              dqdt(i,k,ntrw)    = dvdftra(i,k,4)
              dqdt(i,k,nqrimef) = dvdftra(i,k,5)
              dqdt(i,k,ntoz)    = dvdftra(i,k,6)
            enddo
          enddo

        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          if(ltaerosol) then
            do k=1,levs
              do i=1,im
                dqdt(i,k,ntqv)  = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntrw)  = dvdftra(i,k,4)
                dqdt(i,k,ntsw)  = dvdftra(i,k,5)
                dqdt(i,k,ntgl)  = dvdftra(i,k,6)
                dqdt(i,k,ntlnc) = dvdftra(i,k,7)
                dqdt(i,k,ntinc) = dvdftra(i,k,8)
                dqdt(i,k,ntrnc) = dvdftra(i,k,9)
                dqdt(i,k,ntoz)  = dvdftra(i,k,10)
                dqdt(i,k,ntwa)  = dvdftra(i,k,11)
                dqdt(i,k,ntia)  = dvdftra(i,k,12)
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                dqdt(i,k,ntqv)  = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntrw)  = dvdftra(i,k,4)
                dqdt(i,k,ntsw)  = dvdftra(i,k,5)
                dqdt(i,k,ntgl)  = dvdftra(i,k,6)
                dqdt(i,k,ntinc) = dvdftra(i,k,7)
                dqdt(i,k,ntrnc) = dvdftra(i,k,8)
                dqdt(i,k,ntoz)  = dvdftra(i,k,9)
              enddo
            enddo
          endif
        elseif (imp_physics == imp_physics_mg) then          ! MG3/2
          if (ntgl > 0) then                                 ! MG
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)     = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntrw)  = dvdftra(i,k,4)
                dqdt(i,k,ntsw)  = dvdftra(i,k,5)
                dqdt(i,k,ntgl)  = dvdftra(i,k,6)
                dqdt(i,k,ntlnc) = dvdftra(i,k,7)
                dqdt(i,k,ntinc) = dvdftra(i,k,8)
                dqdt(i,k,ntrnc) = dvdftra(i,k,9)
                dqdt(i,k,ntsnc) = dvdftra(i,k,10)
                dqdt(i,k,ntgnc) = dvdftra(i,k,11)
                dqdt(i,k,ntoz)  = dvdftra(i,k,12)
              enddo
            enddo
          else                                               ! MG2
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)     = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntrw)  = dvdftra(i,k,4)
                dqdt(i,k,ntsw)  = dvdftra(i,k,5)
                dqdt(i,k,ntlnc) = dvdftra(i,k,6)
                dqdt(i,k,ntinc) = dvdftra(i,k,7)
                dqdt(i,k,ntrnc) = dvdftra(i,k,8)
                dqdt(i,k,ntsnc) = dvdftra(i,k,9)
                dqdt(i,k,ntoz)  = dvdftra(i,k,10)
              enddo
            enddo
          endif
        elseif (imp_physics == imp_physics_gfdl) then        ! GFDL MP
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntqv) = dvdftra(i,k,1)
              dqdt(i,k,ntcw) = dvdftra(i,k,2)
              dqdt(i,k,ntiw) = dvdftra(i,k,3)
              dqdt(i,k,ntrw) = dvdftra(i,k,4)
              dqdt(i,k,ntsw) = dvdftra(i,k,5)
              dqdt(i,k,ntgl) = dvdftra(i,k,6)
              dqdt(i,k,ntoz) = dvdftra(i,k,7)
            enddo
          enddo
        elseif (imp_physics == imp_physics_zhao_carr) then
          do k=1,levs
            do i=1,im
              dqdt(i,k,1)    = dvdftra(i,k,1)
              dqdt(i,k,ntcw) = dvdftra(i,k,2)
              dqdt(i,k,ntoz) = dvdftra(i,k,3)
            enddo
          enddo
        endif

      endif ! nvdiff == ntrac

      if (cplchm) then
        do i = 1, im
          tem  = prsl(i,1) / (rd*t1(i)*(one+fvirt*max(q1(i), qmin)))
          ushfsfci(i) = -cp * tem * hflx(i) ! upward sensible heat flux
        enddo
        dkt_cpl(1:im,1:levs) = dkt(1:im,1:levs)
      endif


!  --- ...  coupling insertion

      if (cplflx) then
        do i=1,im
          if (oceanfrac(i) > zero) then ! Ocean only, NO LAKES
            if ( .not. wet(i)) then ! no open water
              if ( kdt > 1 ) then !use results from CICE
                dusfci_cpl(i) = dusfc_cice(i)
                dvsfci_cpl(i) = dvsfc_cice(i)
                dtsfci_cpl(i) = dtsfc_cice(i)
                dqsfci_cpl(i) = dqsfc_cice(i)
              else !use PBL fluxes when CICE fluxes is unavailable
                dusfci_cpl(i) = dusfc1(i)
                dvsfci_cpl(i) = dvsfc1(i)
                dtsfci_cpl(i) = dtsfc1(i)
                dqsfci_cpl(i) = dqsfc1(i)
              end if
            elseif (icy(i) .or. dry(i)) then ! use stress_ocean from sfc_diff for opw component at mixed point
              rho = prsl(i,1) / (rd*t1(i)*(one+fvirt*max(q1(i), qmin)))
              if (wind(i) > zero) then
                tem = - rho * stress_wat(i) / wind(i)
                dusfci_cpl(i) = tem * ugrs1(i)   ! U-momentum flux
                dvsfci_cpl(i) = tem * vgrs1(i)   ! V-momentum flux
              else
                dusfci_cpl(i) = zero
                dvsfci_cpl(i) = zero
              endif
              dtsfci_cpl(i) = cp   * rho * hflx_wat(i) ! sensible heat flux over open ocean
              dqsfci_cpl(i) = hvap * rho * evap_wat(i) ! latent heat flux over open ocean
            else                                       ! use results from PBL scheme for 100% open ocean
              dusfci_cpl(i) = dusfc1(i)
              dvsfci_cpl(i) = dvsfc1(i)
              dtsfci_cpl(i) = dtsfc1(i)*hffac(i)
              dqsfci_cpl(i) = dqsfc1(i)*hefac(i)
            endif
!
            dusfc_cpl (i) = dusfc_cpl(i) + dusfci_cpl(i) * dtf
            dvsfc_cpl (i) = dvsfc_cpl(i) + dvsfci_cpl(i) * dtf
            dtsfc_cpl (i) = dtsfc_cpl(i) + dtsfci_cpl(i) * dtf
            dqsfc_cpl (i) = dqsfc_cpl(i) + dqsfci_cpl(i) * dtf
!
          else
            dusfc_cpl(i) = huge
            dvsfc_cpl(i) = huge
            dtsfc_cpl(i) = huge
            dqsfc_cpl(i) = huge
!!
          endif ! Ocean only, NO LAKES
        enddo
      endif

!-------------------------------------------------------lssav if loop ----------
      if (lssav) then
        do i=1,im
          dusfc_diag (i) = dusfc_diag(i) + dusfc1(i) * dtf
          dvsfc_diag (i) = dvsfc_diag(i) + dvsfc1(i) * dtf
          dusfci_diag(i) = dusfc1(i)
          dvsfci_diag(i) = dvsfc1(i)
          dtsfci_diag(i) = dtsfc1(i)*hffac(i)
          dqsfci_diag(i) = dqsfc1(i)*hefac(i)
          dtsfc_diag (i) = dtsfc_diag(i) + dtsfci_diag(i) * dtf
          dqsfc_diag (i) = dqsfc_diag(i) + dqsfci_diag(i) * dtf
        enddo

        if (ldiag3d .and. flag_for_pbl_generic_tend) then
          if (lsidea) then
            idtend = dtidx(index_of_temperature, index_of_process_pbl)
            if(idtend>=1) then
              dtend(1:im,1:levs,idtend) = dtend(1:im,1:levs,idtend) + dtdt(1:im,1:levs)*dtf
            endif
          else
            idtend = dtidx(index_of_temperature, index_of_process_pbl)
            if(idtend>=1) then
              dtend(1:im,1:levs,idtend) = dtend(1:im,1:levs,idtend) + (tgrs(1:im,1:levs) - save_t(1:im,1:levs))
            endif
          endif
          idtend = dtidx(index_of_x_wind, index_of_process_pbl)
          if(idtend>=1) then
            dtend(1:im,1:levs,idtend) = dtend(1:im,1:levs,idtend) + (ugrs(1:im,1:levs) - save_u(1:im,1:levs))
          endif
          idtend = dtidx(index_of_y_wind, index_of_process_pbl)
          if(idtend>=1) then
            dtend(1:im,1:levs,idtend) = dtend(1:im,1:levs,idtend) + (vgrs(1:im,1:levs) - save_v(1:im,1:levs))
          endif
          idtend = dtidx(100+ntqv, index_of_process_pbl)
          if(idtend>=1) then
            dtend(1:im,1:levs,idtend) = dtend(1:im,1:levs,idtend) + qgrs(1:im,1:levs,ntqv) - save_q(1:im,1:levs,ntqv)
          endif
          idtend = dtidx(100+ntoz, index_of_process_pbl)
          if(idtend>=1) then
            dtend(1:im,1:levs,idtend) = dtend(1:im,1:levs,idtend) + qgrs(1:im,1:levs,ntoz) - save_q(1:im,1:levs,ntoz)
          endif
          idtend = dtidx(100+ntke, index_of_process_pbl)
          if(idtend>=1) then
            dtend(1:im,1:levs,idtend) = dtend(1:im,1:levs,idtend) + (qgrs(1:im,1:levs,ntke) - save_q(1:im,1:levs,ntke))
          endif
        endif

      endif   ! end if_lssav

      end subroutine GFS_PBL_generic_post_run

      end module GFS_PBL_generic_post

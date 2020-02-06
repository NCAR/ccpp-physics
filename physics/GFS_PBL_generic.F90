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
      subroutine GFS_PBL_generic_pre_run (im, levs, nvdiff, ntrac,                       &
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc,                 &
        ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef, trans_aero, ntchs, ntchm,                   &
        imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,           &
        imp_physics_zhao_carr, imp_physics_mg, imp_physics_fer_hires, cplchm, ltaerosol, hybedmf, do_shoc,      &
        satmedmf, qgrs, vdftra, errmsg, errflg)

      use machine,                only : kind_phys
      use GFS_PBL_generic_common, only : set_aerosol_tracer_index

      implicit none

      integer, intent(in) :: im, levs, nvdiff, ntrac
      integer, intent(in) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc
      integer, intent(in) :: ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef,ntchs, ntchm
      logical, intent(in) :: trans_aero
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      integer, intent(in) :: imp_physics_zhao_carr, imp_physics_mg, imp_physics_fer_hires
      logical, intent(in) :: cplchm, ltaerosol, hybedmf, do_shoc, satmedmf

      real(kind=kind_phys), dimension(im, levs, ntrac), intent(in) :: qgrs
      real(kind=kind_phys), dimension(im, levs, nvdiff), intent(inout) :: vdftra

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      !local variables
      integer :: i, k, kk, k1, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!DH: dvdftra is only used if nvdiff != ntrac or (nvdiff == ntrac .and. )
      if (nvdiff == ntrac .and. (hybedmf .or. do_shoc .or. satmedmf)) then
        vdftra = qgrs
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
                vdftra(i,k,6)  = qgrs(i,k,ntlnc)
                vdftra(i,k,7)  = qgrs(i,k,ntinc)
                vdftra(i,k,8)  = qgrs(i,k,ntoz)
                vdftra(i,k,9)  = qgrs(i,k,ntwa)
                vdftra(i,k,10) = qgrs(i,k,ntia)
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = qgrs(i,k,ntqv)
                vdftra(i,k,2) = qgrs(i,k,ntcw)
                vdftra(i,k,3) = qgrs(i,k,ntiw)
                vdftra(i,k,4) = qgrs(i,k,ntrw)
                vdftra(i,k,5) = qgrs(i,k,ntsw)
                vdftra(i,k,6) = qgrs(i,k,ntinc)
                vdftra(i,k,7) = qgrs(i,k,ntoz)
              enddo
            enddo
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
        elseif (imp_physics == imp_physics_zhao_carr) then
! Zhao/Carr/Sundqvist
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntoz)
            enddo
          enddo
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
        ltaerosol, cplflx, cplchm, lssav, flag_for_pbl_generic_tend, ldiag3d, qdiag3d, lsidea, hybedmf, do_shoc, satmedmf, shinhong, do_ysu,              &
        dvdftra, dusfc1, dvsfc1, dtsfc1, dqsfc1, dtf, dudt, dvdt, dtdt, htrsw, htrlw, xmu,                                     &
        dqdt, dusfc_cpl, dvsfc_cpl, dtsfc_cpl,                                                                                 &
        dqsfc_cpl, dusfci_cpl, dvsfci_cpl, dtsfci_cpl, dqsfci_cpl, dusfc_diag, dvsfc_diag, dtsfc_diag, dqsfc_diag,             &
        dusfci_diag, dvsfci_diag, dtsfci_diag, dqsfci_diag, dt3dt, du3dt_PBL, du3dt_OGWD, dv3dt_PBL, dv3dt_OGWD, dq3dt,        &
        dq3dt_ozone, rd, cp,fvirt, hvap, t1, q1, prsl, hflx, ushfsfci, oceanfrac, fice, dusfc_cice, dvsfc_cice, dtsfc_cice,    &
        dqsfc_cice, wet, dry, icy, wind, stress_ocn, hflx_ocn, evap_ocn, ugrs1, vgrs1, dkt_cpl, dkt, errmsg, errflg)

      use machine,                only : kind_phys
      use GFS_PBL_generic_common, only : set_aerosol_tracer_index

      implicit none

      integer, intent(in) :: im, levs, nvdiff, ntrac, ntchs, ntchm
      integer, intent(in) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc, ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef
      logical, intent(in) :: trans_aero
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      integer, intent(in) :: imp_physics_zhao_carr, imp_physics_mg, imp_physics_fer_hires
      logical, intent(in) :: ltaerosol, cplflx, cplchm, lssav, ldiag3d, qdiag3d, lsidea
      logical, intent(in) :: hybedmf, do_shoc, satmedmf, shinhong, do_ysu
      logical, intent(in) :: flag_for_pbl_generic_tend

      real(kind=kind_phys), intent(in) :: dtf
      real(kind=kind_phys), intent(in) :: rd, cp, fvirt, hvap
      real(kind=kind_phys), dimension(:), intent(in) :: t1, q1, hflx, oceanfrac, fice
      real(kind=kind_phys), dimension(:,:), intent(in) :: prsl
      real(kind=kind_phys), dimension(:), intent(in) :: dusfc_cice, dvsfc_cice, dtsfc_cice, dqsfc_cice, &
          wind, stress_ocn, hflx_ocn, evap_ocn, ugrs1, vgrs1
      real(kind=kind_phys), dimension(im, levs, nvdiff), intent(in) :: dvdftra
      real(kind=kind_phys), dimension(im), intent(in) :: dusfc1, dvsfc1, dtsfc1, dqsfc1, xmu
      real(kind=kind_phys), dimension(im, levs), intent(in) :: dudt, dvdt, dtdt, htrsw, htrlw

      real(kind=kind_phys), dimension(im, levs, ntrac), intent(inout) :: dqdt

      ! The following arrays may not be allocated, depending on certain flags (cplflx, ...).
      ! Since Intel 15 crashes when passing unallocated arrays to arrays defined with explicit shape,
      ! use assumed-shape arrays. Note that Intel 18 and GNU 6.2.0-8.1.0 tolerate explicit-shape arrays
      ! as long as these do not get used when not allocated.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dt3dt, du3dt_PBL, du3dt_OGWD, dv3dt_PBL, dv3dt_OGWD, dq3dt, dq3dt_ozone
      real(kind=kind_phys), dimension(:), intent(inout) :: dusfc_cpl, dvsfc_cpl, dtsfc_cpl, dqsfc_cpl, dusfci_cpl, dvsfci_cpl, &
        dtsfci_cpl, dqsfci_cpl, dusfc_diag, dvsfc_diag, dtsfc_diag, dqsfc_diag, dusfci_diag, dvsfci_diag, dtsfci_diag, dqsfci_diag

      logical, dimension(:),intent(in) :: wet, dry, icy
      real(kind=kind_phys), dimension(:), intent(out) ::  ushfsfci

      real(kind=kind_phys), dimension(:,:), intent(inout) :: dkt_cpl
      real(kind=kind_phys), dimension(:,:), intent(in) :: dkt

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k, kk, k1, n
      real(kind=kind_phys) :: tem, tem1, rho

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
              dqdt(i,k,ntqv)  = dvdftra(i,k,1)
              dqdt(i,k,ntcw)  = dvdftra(i,k,2)
              dqdt(i,k,ntiw)  = dvdftra(i,k,3)
              dqdt(i,k,ntrw)  = dvdftra(i,k,4)
              dqdt(i,k,nqrimef) = dvdftra(i,k,5)
              dqdt(i,k,ntoz)  = dvdftra(i,k,6)
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
                dqdt(i,k,ntlnc) = dvdftra(i,k,6)
                dqdt(i,k,ntinc) = dvdftra(i,k,7)
                dqdt(i,k,ntoz)  = dvdftra(i,k,8)
                dqdt(i,k,ntwa)  = dvdftra(i,k,9)
                dqdt(i,k,ntia)  = dvdftra(i,k,10)
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
                dqdt(i,k,ntinc) = dvdftra(i,k,6)
                dqdt(i,k,ntoz)  = dvdftra(i,k,7)
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
          tem1 = max(q1(i), 1.e-8)
          tem  = prsl(i,1) / (rd*t1(i)*(1.0+fvirt*tem1))
          ushfsfci(i) = -cp * tem * hflx(i) ! upward sensible heat flux
        enddo
        ! dkt_cpl has dimensions (1:im,1:levs), but dkt has (1:im,1:levs-1)
        dkt_cpl(1:im,1:levs-1) = dkt(1:im,1:levs-1)
      endif


!  --- ...  coupling insertion

      if (cplflx) then
        do i=1,im
          if (oceanfrac(i) > 0.0) then ! Ocean only, NO LAKES
!            if (fice(i) == ceanfrac(i)) then ! use results from CICE
!              dusfci_cpl(i) = dusfc_cice(i)
!              dvsfci_cpl(i) = dvsfc_cice(i)
!              dtsfci_cpl(i) = dtsfc_cice(i)
!              dqsfci_cpl(i) = dqsfc_cice(i)
!            elseif (dry(i) .or. icy(i)) then   ! use stress_ocean from sfc_diff for opw component at mixed point
            if (wet(i)) then                   ! use stress_ocean from sfc_diff for opw component at mixed point
              if (icy(i) .or. dry(i)) then
                tem1 = max(q1(i), 1.e-8)
                rho = prsl(i,1) / (rd*t1(i)*(1.0+fvirt*tem1))
                if (wind(i) > 0.0) then
                  tem = - rho * stress_ocn(i) / wind(i)
                  dusfci_cpl(i) = tem * ugrs1(i)   ! U-momentum flux
                  dvsfci_cpl(i) = tem * vgrs1(i)   ! V-momentum flux
                else
                  dusfci_cpl(i) = 0.0
                  dvsfci_cpl(i) = 0.0
                endif
                dtsfci_cpl(i) = cp   * rho * hflx_ocn(i) ! sensible heat flux over open ocean
                dqsfci_cpl(i) = hvap * rho * evap_ocn(i) ! latent heat flux over open ocean
              else                                                    ! use results from PBL scheme for 100% open ocean
                dusfci_cpl(i) = dusfc1(i)
                dvsfci_cpl(i) = dvsfc1(i)
                dtsfci_cpl(i) = dtsfc1(i)
                dqsfci_cpl(i) = dqsfc1(i)
              endif
            endif
!
            dusfc_cpl (i) = dusfc_cpl(i) + dusfci_cpl(i) * dtf
            dvsfc_cpl (i) = dvsfc_cpl(i) + dvsfci_cpl(i) * dtf
            dtsfc_cpl (i) = dtsfc_cpl(i) + dtsfci_cpl(i) * dtf
            dqsfc_cpl (i) = dqsfc_cpl(i) + dqsfci_cpl(i) * dtf
!!
          endif ! Ocean only, NO LAKES
        enddo
      endif

!-------------------------------------------------------lssav if loop ----------
      if (lssav) then
        do i=1,im
          dusfc_diag (i) = dusfc_diag(i) + dusfc1(i)*dtf
          dvsfc_diag (i) = dvsfc_diag(i) + dvsfc1(i)*dtf
          dtsfc_diag (i) = dtsfc_diag(i) + dtsfc1(i)*dtf
          dqsfc_diag (i) = dqsfc_diag(i) + dqsfc1(i)*dtf
          dusfci_diag(i) = dusfc1(i)
          dvsfci_diag(i) = dvsfc1(i)
          dtsfci_diag(i) = dtsfc1(i)
          dqsfci_diag(i) = dqsfc1(i)
        enddo
  !       if (lprnt) then
  !         write(0,*)' dusfc=',dusfc(ipr),' dusfc1=',dusfc1(ipr),' dtf=',
  !    &     dtf,' kdt=',kdt,' lat=',lat
  !       endif

        if (ldiag3d .and. flag_for_pbl_generic_tend) then
          if (lsidea) then
            dt3dt(1:im,:) = dt3dt(1:im,:) + dtdt(1:im,:)*dtf
          else
            do k=1,levs
              do i=1,im
                tem  = dtdt(i,k) - (htrlw(i,k)+htrsw(i,k)*xmu(i))
                dt3dt(i,k) = dt3dt(i,k) + tem*dtf
              enddo
            enddo
          endif
          do k=1,levs
            do i=1,im
              du3dt_PBL(i,k) = du3dt_PBL(i,k) + dudt(i,k) * dtf
              du3dt_OGWD(i,k) = du3dt_OGWD(i,k) - dudt(i,k) * dtf
              dv3dt_PBL(i,k) = dv3dt_PBL(i,k) + dvdt(i,k) * dtf
              dv3dt_OGWD(i,k) = dv3dt_OGWD(i,k) - dvdt(i,k) * dtf
            enddo
          enddo
          if (qdiag3d) then
             do k=1,levs
                do i=1,im
                   tem  = dqdt(i,k,ntqv) * dtf
                   dq3dt(i,k) = dq3dt(i,k) + tem
                enddo
             enddo
             if (ntoz > 0) then
                do k=1,levs
                   do i=1,im
                      dq3dt_ozone(i,k) = dq3dt_ozone(i,k) + dqdt(i,k,ntoz) * dtf
                   enddo
                enddo
             endif
          endif
        endif

      endif   ! end if_lssav

      end subroutine GFS_PBL_generic_post_run

      end module GFS_PBL_generic_post

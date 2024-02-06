!> \file GFS_PBL_generic_pre.F90
!!  Contains code related to PBL schemes to be called prior to PBL schemes within GFS-based physics suites.

      module GFS_PBL_generic_pre

      contains

!> \brief This scheme sets up the vertically diffused tracer array for any PBL scheme based on the microphysics scheme chosen
!! \section arg_table_GFS_PBL_generic_pre_run Argument Table
!! \htmlinclude GFS_PBL_generic_pre_run.html
!!
      subroutine GFS_PBL_generic_pre_run (im, levs, nvdiff, ntrac, rtg_ozone_index,      &
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntsnc, ntgnc,                 &
        ntwa, ntia, ntgl, ntoz, ntke, ntkev, nqrimef, trans_aero, ntchs, ntchm,          &
        ntccn, nthl, nthnc, ntgv, nthv, ntrz, ntgz, nthz,                                &
        imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,           &
        imp_physics_zhao_carr, imp_physics_mg, imp_physics_fer_hires, imp_physics_nssl,  &
        ltaerosol, mraerosol, nssl_ccn_on, nssl_hail_on, nssl_3moment,                   &
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
      integer, intent(in) :: ntccn, nthl, nthnc, ntgv, nthv, ntrz, ntgz, nthz
      logical, intent(in) :: trans_aero, ldiag3d, qdiag3d, lssav
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      integer, intent(in) :: imp_physics_zhao_carr, imp_physics_mg, imp_physics_fer_hires
      logical, intent(in) :: ltaerosol, hybedmf, do_shoc, satmedmf, flag_for_pbl_generic_tend, mraerosol
      integer, intent(in) :: imp_physics_nssl
      logical, intent(in) :: nssl_hail_on, nssl_ccn_on, nssl_3moment

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
          elseif(mraerosol) then
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
        elseif (imp_physics == imp_physics_nssl ) then
  ! nssl
            IF ( nssl_hail_on ) THEN
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,nthl)
                vdftra(i,k,8)  = qgrs(i,k,ntlnc)
                vdftra(i,k,9)  = qgrs(i,k,ntinc)
                vdftra(i,k,10) = qgrs(i,k,ntrnc)
                vdftra(i,k,11) = qgrs(i,k,ntsnc)
                vdftra(i,k,12) = qgrs(i,k,ntgnc)
                vdftra(i,k,13) = qgrs(i,k,nthnc)
                vdftra(i,k,14) = qgrs(i,k,ntgv)
                vdftra(i,k,15) = qgrs(i,k,nthv)
                vdftra(i,k,16) = qgrs(i,k,ntoz)
                n = 16
                IF ( nssl_ccn_on ) THEN
                 vdftra(i,k,n+1)  = qgrs(i,k,ntccn)
                 n = n+1
                ENDIF
                IF ( nssl_3moment ) THEN
                 vdftra(i,k,n+1)  = qgrs(i,k,ntrz)
                 vdftra(i,k,n+2)  = qgrs(i,k,ntgz)
                 vdftra(i,k,n+3)  = qgrs(i,k,nthz)
                 n = n+3
                ENDIF
              enddo
            enddo
            
            ELSE
            ! no hail
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
                vdftra(i,k,12) = qgrs(i,k,ntgv)
                vdftra(i,k,13) = qgrs(i,k,ntoz)
                 n = 13
                IF ( nssl_ccn_on ) THEN
                 vdftra(i,k,n+1) = qgrs(i,k,ntccn)
                 n = n+1
                ENDIF
                IF ( nssl_3moment ) THEN
                 vdftra(i,k,n+1) = qgrs(i,k,ntrz)
                 vdftra(i,k,n+2) = qgrs(i,k,ntgz)
                 n = n+2
                ENDIF
              enddo
            enddo
            
            ENDIF


        endif
!
        if (trans_aero) then
          call set_aerosol_tracer_index(imp_physics, imp_physics_wsm6,          &
                                        imp_physics_thompson, ltaerosol,mraerosol, &
                                        imp_physics_mg, ntgl, imp_physics_gfdl, &
                                        imp_physics_zhao_carr, imp_physics_nssl,&
                                        nssl_hail_on, nssl_ccn_on, kk,          &
                                        errmsg, errflg)
          if (errflg /= 0) return
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

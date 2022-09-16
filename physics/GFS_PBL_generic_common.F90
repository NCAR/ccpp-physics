!> \file GFS_PBL_generic_common.F90
!!  Contains code used in both pre/post PBL-related interstitial schemes to be used within the GFS physics suite.

      module GFS_PBL_generic_common

      implicit none

      private

      public :: set_aerosol_tracer_index

      contains

      subroutine set_aerosol_tracer_index(imp_physics, imp_physics_wsm6,          &
                                          imp_physics_thompson, ltaerosol,mraerosol,   &
                                          imp_physics_mg, ntgl, imp_physics_gfdl, &
                                          imp_physics_zhao_carr, imp_physics_nssl,&
                                          nssl_hail_on, nssl_ccn_on, kk,          &
                                          errmsg, errflg)
      implicit none
      !
      integer, intent(in )          :: imp_physics, imp_physics_wsm6,          &
                                       imp_physics_thompson,                   &
                                       imp_physics_mg, ntgl, imp_physics_gfdl, &
                                       imp_physics_zhao_carr,imp_physics_nssl
      logical, intent(in )          :: ltaerosol, mraerosol, nssl_hail_on, nssl_ccn_on
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
          kk = 12
        else if(mraerosol) then
          kk = 10
        else
          kk = 9
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
      elseif (imp_physics == imp_physics_nssl) then
        IF ( nssl_hail_on ) THEN
          kk = 16
        ELSE
          kk = 13
        ENDIF
        IF ( nssl_ccn_on ) kk = kk + 1
      else
        write(errmsg,'(*(a))') 'Logic error: unknown microphysics option in set_aerosol_tracer_index'
        kk = -999
        errflg = 1
        return
      endif

      end subroutine set_aerosol_tracer_index

      end module GFS_PBL_generic_common

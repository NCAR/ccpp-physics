!>\file rrtmg_lw_pre.f90
!! This file contains a call to module_radiation_surface::setemis() to
!! setup surface emissivity for LW radiation.
      module rrtmg_lw_pre
      contains

!>\defgroup rrtmg_lw_pre GFS RRTMG scheme pre
!! @{
      subroutine rrtmg_lw_pre_init ()
      end subroutine rrtmg_lw_pre_init 

!> \section arg_table_rrtmg_lw_pre_run Argument Table
!! \htmlinclude rrtmg_lw_pre_run.html
!!
      subroutine rrtmg_lw_pre_run (im, lslwr, kdt, lsm, lsm_noahmp, lsm_ruc, vtype, &
        xlat, xlon, slmsk, snowd, sncovr, sncovr_ice, zorl, hprime, tsfg, tsfa,     &
        semis_lnd, semis_ice, semisbase, semis, errmsg, errflg)
    
      use machine,                   only: kind_phys
      use module_radiation_surface,  only: setemis

      implicit none
      
      integer,                              intent(in)  :: im
      logical,                              intent(in)  :: lslwr
      integer, intent(in) :: kdt, lsm, lsm_noahmp, lsm_ruc

      real(kind=kind_phys), dimension(im),  intent(in)  :: xlat, xlon, vtype, slmsk,&
        snowd, sncovr, sncovr_ice, zorl, hprime, tsfg, tsfa
      real(kind=kind_phys), dimension(:),   intent(in)  :: semis_lnd 
      real(kind=kind_phys), dimension(:),   intent(in)  :: semis_ice 
      real(kind=kind_phys), dimension(im),  intent(out) :: semisbase
      real(kind=kind_phys), dimension(im),  intent(out) :: semis
      character(len=*),                     intent(out) :: errmsg
      integer,                              intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lslwr) then
!>  - Call module_radiation_surface::setemis(),to setup surface
!! emissivity for LW radiation.
        call setemis (kdt, lsm, lsm_noahmp, lsm_ruc, vtype, xlon, xlat, slmsk, &
                      snowd, sncovr, sncovr_ice, zorl, tsfg, tsfa,             &
                      hprime, semis_lnd, semis_ice, im,                        & !  ---  inputs
                      semisbase, semis)                                          !  ---  outputs

      endif

      end subroutine rrtmg_lw_pre_run

       subroutine rrtmg_lw_pre_finalize ()
       end subroutine rrtmg_lw_pre_finalize
!! @}
       end module rrtmg_lw_pre

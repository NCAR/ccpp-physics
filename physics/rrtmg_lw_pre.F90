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
      subroutine rrtmg_lw_pre_run (im, lslwr, xlat, xlon, landfrac, snowd, sncovr,&
                                   fice, zorl, hprime, tsfg, tsfa, semis, emiss, errmsg, errflg)
    
      use machine,                   only: kind_phys
      use module_radiation_surface,  only: setemis

      implicit none
      
      integer,                              intent(in)  :: im
      logical,                              intent(in)  :: lslwr
      real(kind=kind_phys), dimension(:),   intent(in)  :: xlat, xlon, landfrac,  &
        snowd, sncovr, fice, zorl, hprime, tsfg, tsfa, emiss

      real(kind=kind_phys), dimension(:),   intent(out) :: semis
      character(len=*),                     intent(out) :: errmsg
      integer,                              intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lslwr) then
!>  - Call module_radiation_surface::setemis(),to setup surface
!! emissivity for LW radiation.
        call setemis (xlon, xlat, landfrac, snowd, sncovr, fice, zorl, tsfg, tsfa, &
                      hprime, emiss, im,                                           & !  ---  inputs
                      semis)                              !  ---  outputs
      endif

      end subroutine rrtmg_lw_pre_run

       subroutine rrtmg_lw_pre_finalize ()
       end subroutine rrtmg_lw_pre_finalize
!! @}
       end module rrtmg_lw_pre

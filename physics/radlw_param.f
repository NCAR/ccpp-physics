!>  \file radlw_param.f
!!  This file contains LW band parameters setup.

!!!!!  ==============================================================  !!!!!
!!!!!             lw-rrtm3 radiation package description               !!!!!
!!!!!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtm-lw radiation    !
!   code from aer inc.                                                     !
!                                                                          !
!    the rrtm3 package includes these parts:                               !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    the 'radlw_rrtm3_param.f' contains:                                   !
!                                                                          !
!       'module_radlw_parameters'  -- band parameters set up               !
!                                                                          !
!    the 'radlw_rrtm3_datatb.f' contains:                                  !
!                                                                          !
!       'module_radlw_avplank'     -- plank flux data                      !
!       'module_radlw_ref'         -- reference temperature and pressure   !
!       'module_radlw_cldprlw'     -- cloud property coefficients          !
!       'module_radlw_kgbnn'       -- absorption coeffients for 16         !
!                                     bands, where nn = 01-16              !
!                                                                          !
!    the 'radlw_rrtm3_main.f' contains:                                    !
!                                                                          !
!       'rrtmg_lw'        -- main lw radiation transfer                    !
!                                                                          !
!    in the main module 'rrtmg_lw' there are only two                      !
!    externally callable subroutines:                                      !
!                                                                          !
!       'lwrad'     -- main rrtm3 lw radiation routine                     !
!       'rlwinit'   -- to initialize rrtm3 lw radiation                    !
!                                                                          !
!    all the lw radiation subprograms become contained subprograms         !
!    in module 'module_radlw_rrtm' and many of them are not directly       !
!    accessable from places outside the module.                            !
!                                                                          !
!    compilation sequence is:                                              !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    and all should be put in front of routines that use lw modules        !
!                                                                          !
!    ncep modifications history log:                                       !
!                                                                          !
!       see list in program "radlw_rrtm3_main.f"                           !
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!

!>\ingroup module_radlw_main
!> This module contains LW band parameters set up.
      module module_radlw_parameters     !
!........................................!

!! \section arg_table_module_radlw_parameters
!! \htmlinclude module_radlw_parameters.html
!!

      use physparam,               only : kind_phys

      implicit none
!
      public
!
!> derived type for LW fluxes at top of atmosphere
!! \section arg_table_topflw_type
!! \htmlinclude topflw_type.html
!!
      type topflw_type                      !< define type construct for radiation fluxes at toa
        real (kind=kind_phys) :: upfxc      !< total sky upward flux at toa
        real (kind=kind_phys) :: upfx0      !< clear sky upward flux at toa
      end type topflw_type
!
!> derived type for LW fluxes at surface
!! \section arg_table_sfcflw_type
!! \htmlinclude sfcflw_type.html
!!
      type sfcflw_type                      !< define type construct for radiation fluxes at surface
        real (kind=kind_phys) :: upfxc      !< total sky upward flux at sfc
        real (kind=kind_phys) :: upfx0      !< clear sky upward flux at sfc
        real (kind=kind_phys) :: dnfxc      !< total sky downward flux at sfc
        real (kind=kind_phys) :: dnfx0      !< clear sky downward flux at sfc
      end type sfcflw_type
!
      type proflw_type                      !< define type construct for optional radiation flux profiles
        real (kind=kind_phys) :: upfxc      !< level up flux for total sky
        real (kind=kind_phys) :: dnfxc      !< level down flux for total sky
        real (kind=kind_phys) :: upfx0      !< level up for clear sky
        real (kind=kind_phys) :: dnfx0      !< level down flux for clear sky
      end type proflw_type
!
! Parameter constants for LW band structures
      integer, parameter :: NBANDS = 16     !< number of total spectral bands
      integer, parameter :: NGPTLW = 140    !< number of total g-points   
      integer, parameter :: NTBL   = 10000  !< lookup table dimension 
      integer, parameter :: MAXGAS = 7      !< maximum number of absorbing gases
      integer, parameter :: MAXXSEC= 4      !< number of halocarbon gases
      integer, parameter :: NRATES = 6      !< number of ref rates of binary species   
      integer, parameter :: NPLNK  = 181    !< dimension for plank function table

      integer, parameter :: NBDLW  = NBANDS

!> Number of g-point in each band
      integer  :: NG01, NG02, NG03, NG04, NG05, NG06, NG07, NG08,       &
     &            NG09, NG10, NG11, NG12, NG13, NG14, NG15, NG16
      parameter (NG01=10, NG02=12, NG03=16, NG04=14, NG05=16, NG06=08,  &
     &           NG07=12, NG08=08, NG09=12, NG10=06, NG11=08, NG12=08,  &
     &           NG13=04, NG14=02, NG15=02, NG16=02)

!> Begining index of each band
      integer  :: NS01, NS02, NS03, NS04, NS05, NS06, NS07, NS08,       &
     &            NS09, NS10, NS11, NS12, NS13, NS14, NS15, NS16
      parameter (NS01=00, NS02=10, NS03=22, NS04=38, NS05=52, NS06=68,  &
     &           NS07=76, NS08=88, NS09=96, NS10=108, NS11=114,         &
     &           NS12=122, NS13=130, NS14=134, NS15=136, NS16=138)

!> band indices for each g-point
      integer, dimension(NGPTLW) :: NGB
      data NGB(:) / 10*1, 12*2, 16*3, 14*4, 16*5,  8*6, 12*7,  8*8,     & ! band  1- 8
     &              12*9, 6*10, 8*11, 8*12, 4*13, 2*14, 2*15, 2*16 /      ! band  9-16

!> Band spectrum structures (wavenumber is 1/cm)
      real (kind=kind_phys) :: wvnlw1(NBANDS), wvnlw2(NBANDS)
      data wvnlw1  /                                                    &
     &         10.,  350.,  500.,  630.,  700.,  820.,  980., 1080.,    &
     &       1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600. /
      data wvnlw2  /                                                    &
     &        350.,  500.,  630.,  700.,  820.,  980., 1080., 1180.,    &
     &       1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250. /

      real (kind=kind_phys) :: delwave(nbands)
      data delwave / 340., 150., 130.,  70., 120., 160., 100., 100.,    &
     &               210.,  90., 320., 280., 170., 130., 220., 650. /

!........................................!
      end module module_radlw_parameters !
!========================================!

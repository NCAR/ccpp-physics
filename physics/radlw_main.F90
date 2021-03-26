!>  \file radlw_main.f
!!  This file contains NCEP's modifications of the rrtmg-lw radiation
!!  code from AER.

!!!!!  ==============================================================  !!!!!
!!!!!               lw-rrtm3 radiation package description             !!!!!
!!!!!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtmg-lw radiation   !
!   code from aer inc.                                                     !
!                                                                          !
!    the lw-rrtm3 package includes these parts:                            !
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
!                                                                          !
!       'lwrad'     -- main lw radiation routine                           !
!          inputs:                                                         !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                         !
!            clouds,icseed,aerosols,sfemis,sfgtmp,                         !
!            dzlyr,delpin,de_lgth,alpha,                                   !
!            npts, nlay, nlp1, lprnt,                                      !
!          outputs:                                                        !
!            hlwc,topflx,sfcflx,cldtau,                                    !
!!         optional outputs:                                               !
!            HLW0,HLWB,FLXPRF)                                             !
!                                                                          !
!       'rlwinit'   -- initialization routine                              !
!          inputs:                                                         !
!           ( me )                                                         !
!          outputs:                                                        !
!           (none)                                                         !
!                                                                          !
!    all the lw radiation subprograms become contained subprograms         !
!    in module 'rrtmg_lw' and many of them are not directly                !
!    accessable from places outside the module.                            !
!                                                                          !
!    derived data type constructs used:                                    !
!                                                                          !
!     1. radiation flux at toa: (from module 'module_radlw_parameters')    !
!          topflw_type   -  derived data type for toa rad fluxes           !
!            upfxc              total sky upward flux at toa               !
!            upfx0              clear sky upward flux at toa               !
!                                                                          !
!     2. radiation flux at sfc: (from module 'module_radlw_parameters')    !
!          sfcflw_type   -  derived data type for sfc rad fluxes           !
!            upfxc              total sky upward flux at sfc               !
!            upfx0              clear sky upward flux at sfc               !
!            dnfxc              total sky downward flux at sfc             !
!            dnfx0              clear sky downward flux at sfc             !
!                                                                          !
!     3. radiation flux profiles(from module 'module_radlw_parameters')    !
!          proflw_type    -  derived data type for rad vertical prof       !
!            upfxc              level upward flux for total sky            !
!            dnfxc              level downward flux for total sky          !
!            upfx0              level upward flux for clear sky            !
!            dnfx0              level downward flux for clear sky          !
!                                                                          !
!    external modules referenced:                                          !
!                                                                          !
!       'module physparam'                                                 !
!       'module physcons'                                                  !
!       'mersenne_twister'                                                 !
!                                                                          !
!    compilation sequence is:                                              !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    and all should be put in front of routines that use lw modules        !
!                                                                          !
!==========================================================================!
!                                                                          !
!    the original aer program declarations:                                !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
! Copyright (c) 2002-2020, Atmospheric & Environmental Research, Inc. (AER)    !
! All rights reserved.                                                         !
!                                                                              !
! Redistribution and use in source and binary forms, with or without           !
! modification, are permitted provided that the following conditions are met:  !
!  * Redistributions of source code must retain the above copyright            !
!    notice, this list of conditions and the following disclaimer.             !
!  * Redistributions in binary form must reproduce the above copyright         !
!    notice, this list of conditions and the following disclaimer in the       !
!    documentation and/or other materials provided with the distribution.      !
!  * Neither the name of Atmospheric & Environmental Research, Inc., nor       !
!    the names of its contributors may be used to endorse or promote products  !
!    derived from this software without specific prior written permission.     !
!                                                                              !
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  !
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    !
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   !
! ARE DISCLAIMED. IN NO EVENT SHALL ATMOSPHERIC & ENVIRONMENTAL RESEARCH, INC.,!
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR       !
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         !
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     !
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      !
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      !
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF       !
! THE POSSIBILITY OF SUCH DAMAGE.                                              !
!                        (http://www.rtweb.aer.com/)                           !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
! ************************************************************************ !
!                                                                          !
!                              rrtmg_lw                                    !
!                                                                          !
!                                                                          !
!                   a rapid radiative transfer model                       !
!                       for the longwave region                            !
!             for application to general circulation models                !
!                                                                          !
!                                                                          !
!            atmospheric and environmental research, inc.                  !
!                        131 hartwell avenue                               !
!                        lexington, ma 02421                               !
!                                                                          !
!                           eli j. mlawer                                  !
!                        jennifer s. delamere                              !
!                         michael j. iacono                                !
!                         shepard a. clough                                !
!                                                                          !
!                                                                          !
!                       email:  miacono@aer.com                            !
!                       email:  emlawer@aer.com                            !
!                       email:  jdelamer@aer.com                           !
!                                                                          !
!        the authors wish to acknowledge the contributions of the          !
!        following people:  steven j. taubman, karen cady-pereira,         !
!        patrick d. brown, ronald e. farren, luke chen, robert bergstrom.  !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    references:                                                           !
!    (rrtmg_lw/rrtm_lw):                                                   !
!      iacono, m.j., j.s. delamere, e.j. mlawer, m.w. shepard,             !
!      s.a. clough, and w.d collins, radiative forcing by long-lived       !
!      greenhouse gases: calculations with the aer radiative transfer      !
!      models, j, geophys. res., 113, d13103, doi:10.1029/2008jd009944,    !
!      2008.                                                               !
!                                                                          !
!      clough, s.a., m.w. shephard, e.j. mlawer, j.s. delamere,            !
!      m.j. iacono, k. cady-pereira, s. boukabara, and p.d. brown:         !
!      atmospheric radiative transfer modeling: a summary of the aer       !
!      codes, j. quant. spectrosc. radiat. transfer, 91, 233-244, 2005.    !
!                                                                          !
!      mlawer, e.j., s.j. taubman, p.d. brown, m.j. iacono, and s.a.       !
!      clough:  radiative transfer for inhomogeneous atmospheres: rrtm,    !
!      a validated correlated-k model for the longwave.  j. geophys. res., !
!      102, 16663-16682, 1997.                                             !
!                                                                          !
!    (mcica):                                                              !
!      pincus, r., h. w. barker, and j.-j. morcrette: a fast, flexible,    !
!      approximation technique for computing radiative transfer in         !
!      inhomogeneous cloud fields, j. geophys. res., 108(d13), 4376,       !
!      doi:10.1029/2002JD003322, 2003.                                     !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    aer's revision history:                                               !
!     this version of rrtmg_lw has been modified from rrtm_lw to use a     !
!     reduced set of g-points for application to gcms.                     !
!                                                                          !
! --  original version (derived from rrtm_lw), reduction of g-points,      !
!     other revisions for use with gcms.                                   !
!        1999: m. j. iacono, aer, inc.                                     !
! --  adapted for use with ncar/cam3.                                      !
!        may 2004: m. j. iacono, aer, inc.                                 !
! --  revised to add mcica capability.                                     !
!        nov 2005: m. j. iacono, aer, inc.                                 !
! --  conversion to f90 formatting for consistency with rrtmg_sw.          !
!        feb 2007: m. j. iacono, aer, inc.                                 !
! --  modifications to formatting to use assumed-shape arrays.             !
!        aug 2007: m. j. iacono, aer, inc.                                 !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    ncep modifications history log:                                       !
!                                                                          !
!       nov 1999,  ken campana       -- received the original code from    !
!                    aer (1998 ncar ccm version), updated to link up with  !
!                    ncep mrf model                                        !
!       jun 2000,  ken campana       -- added option to switch random and  !
!                    maximum/random cloud overlap                          !
!           2001,  shrinivas moorthi -- further updates for mrf model      !
!       may 2001,  yu-tai hou        -- updated on trace gases and cloud   !
!                    property based on rrtm_v3.0 codes.                    !
!       dec 2001,  yu-tai hou        -- rewritten code into fortran 90 std !
!                    set ncep radiation structure standard that contains   !
!                    three plug-in compatable fortran program files:       !
!                    'radlw_param.f', 'radlw_datatb.f', 'radlw_main.f'     !
!                    fixed bugs in subprograms taugb14, taugb2, etc. added !
!                    out-of-bounds protections. (a detailed note of        !
!                    up_to_date modifications/corrections by ncep was sent !
!                    to aer in 2002)                                       !
!       jun 2004,  yu-tai hou        -- added mike iacono's apr 2004       !
!                    modification of variable diffusivity angles.          !
!       apr 2005,  yu-tai hou        -- minor modifications on module      !
!                    structures include rain/snow effect (this version of  !
!                    code was given back to aer in jun 2006)               !
!       mar 2007,  yu-tai hou        -- added aerosol effect for ncep      !
!                    models using the generallized aerosol optical property!
!                    scheme for gfs model.                                 !
!       apr 2007,  yu-tai hou        -- added spectral band heating as an  !
!                    optional output to support the 500 km gfs model's     !
!                    upper stratospheric radiation calculations. and       !
!                    restructure optional outputs for easy access by       !
!                    different models.                                     !
!       oct 2008,  yu-tai hou        -- modified to include new features   !
!                    from aer's newer release v4.4-v4.7, including the     !
!                    mcica sub-grid cloud option. add rain/snow optical    !
!                    properties support to cloudy sky calculations.        !
!                    correct errors in mcica cloud optical properties for  !
!                    ebert & curry scheme (ilwcice=1) that needs band      !
!                    index conversion. simplified and unified sw and lw    !
!                    sub-column cloud subroutines into one module by using !
!                    optional parameters.                                  !
!       mar 2009,  yu-tai hou        -- replaced the original random number!
!                    generator coming from the original code with ncep w3  !
!                    library to simplify the program and moved sub-column  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       oct 2009,  yu-tai hou        -- modified subrtines "cldprop" and   !
!                    "rlwinit" according updats from aer's rrtmg_lw v4.8.  !
!       nov 2009,  yu-tai hou        -- modified subrtine "taumol" according
!                    updats from aer's rrtmg_lw version 4.82. notice the   !
!                    cloud ice/liquid are assumed as in-cloud quantities,  !
!                    not as grid averaged quantities.                      !
!       jun 2010,  yu-tai hou        -- optimized code to improve efficiency
!       apr 2012,  b. ferrier and y. hou -- added conversion factor to fu's!
!                    cloud-snow optical property scheme.                   !
!       nov 2012,  yu-tai hou        -- modified control parameters thru   !
!                     module 'physparam'.                                  !  
!       FEB 2017    A.Cheng   - add odpth output, effective radius input   !
!       jun 2018,  h-m lin/y-t hou   -- added new option of cloud overlap  !
!                     method 'de-correlation-length' for mcica application !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    additional aer revision history:                                      !
!       jul 2020,  m.j. iacono   -- added new mcica cloud overlap options  !
!                     exponential and exponential-random. each method can  !
!                     use either a constant or a latitude-varying and      !
!                     day-of-year varying decorrelation length selected    !
!                     with parameter "idcor".                              !
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!

!> This module contains the CCPP-compliant NCEP's modifications of the 
!! rrtmg-lw radiation code from aer inc.
      module rrtmg_lw 
!
      use physparam,        only : ilwrate, ilwrgas, ilwcliq, ilwcice,  &
     &                             isubclw, icldflg, iovr,  ivflip
      use physcons,         only : con_g, con_cp, con_avgd, con_amd,    &
     &                             con_amw, con_amo3
      use mersenne_twister, only : random_setseed, random_number,       &
     &                             random_stat
!mz
      use machine,          only : kind_phys,                           &
     &                             im => kind_io4, rb => kind_phys

      use module_radlw_parameters
!
      use module_radlw_avplank, only : totplnk
      use module_radlw_ref,     only : preflog, tref, chi_mls
!
      implicit none
!
      private
!
!  ...  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGLW='NCEP LW v5.1  Nov 2012 -RRTMG-LW v4.82  '
!    &   VTAGLW='NCEP LW v5.0  Aug 2012 -RRTMG-LW v4.82  '
!    &   VTAGLW='RRTMG-LW v4.82  Nov 2009  '
!    &   VTAGLW='RRTMG-LW v4.8   Oct 2009  '
!    &   VTAGLW='RRTMG-LW v4.71  Mar 2009  '
!    &   VTAGLW='RRTMG-LW v4.4   Oct 2008  '
!    &   VTAGLW='RRTM-LW v2.3g   Mar 2007  '
!    &   VTAGLW='RRTM-LW v2.3g   Apr 2004  '

!  ---  constant values
      real (kind=kind_phys), parameter :: eps     = 1.0e-6
      real (kind=kind_phys), parameter :: oneminus= 1.0-eps
      real (kind=kind_phys), parameter :: cldmin  = tiny(cldmin)
      real (kind=kind_phys), parameter :: bpade   = 1.0/0.278  ! pade approx constant
      real (kind=kind_phys), parameter :: stpfac  = 296.0/1013.0
      real (kind=kind_phys), parameter :: wtdiff  = 0.5        ! weight for radiance to flux conversion
      real (kind=kind_phys), parameter :: tblint  = ntbl       ! lookup table conversion factor
      real (kind=kind_phys), parameter :: f_zero  = 0.0
      real (kind=kind_phys), parameter :: f_one   = 1.0

!  ...  atomic weights for conversion from mass to volume mixing ratios
      real (kind=kind_phys), parameter :: amdw    = con_amd/con_amw
      real (kind=kind_phys), parameter :: amdo3   = con_amd/con_amo3

!  ...  band indices
      integer, dimension(nbands) :: nspa, nspb

      data nspa / 1, 1, 9, 9, 9, 1, 9, 1, 9, 1, 1, 9, 9, 1, 9, 9 /
      data nspb / 1, 1, 5, 5, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0 /

!  ...  band wavenumber intervals
!     real (kind=kind_phys) :: wavenum1(nbands), wavenum2(nbands)
!     data wavenum1/                                                    &
!    &         10.,  350.,  500.,  630.,  700.,  820.,  980., 1080.,    &
!err &       1180., 1390., 1480., 1800., 2080., 2250., 2390., 2600. /
!    &       1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600. /
!     data wavenum2/                                                    &
!    &        350.,  500.,  630.,  700.,  820.,  980., 1080., 1180.,    &
!err &       1390., 1480., 1800., 2080., 2250., 2390., 2600., 3250. /
!    &       1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250. /
!     real (kind=kind_phys) :: delwave(nbands)
!     data delwave / 340., 150., 130.,  70., 120., 160., 100., 100.,    &
!    &               210.,  90., 320., 280., 170., 130., 220., 650. /

!  ---  reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
!       and 1.80) as a function of total column water vapor.  the function
!       has been defined to minimize flux and cooling rate errors in these bands
!       over a wide range of precipitable water values.
      real (kind=kind_phys), dimension(nbands) :: a0, a1, a2

      data a0 / 1.66,  1.55,  1.58,  1.66,  1.54, 1.454,  1.89,  1.33,  &
     &         1.668,  1.66,  1.66,  1.66,  1.66,  1.66,  1.66,  1.66 /
      data a1 / 0.00,  0.25,  0.22,  0.00,  0.13, 0.446, -0.10,  0.40,  &
     &        -0.006,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /
      data a2 / 0.00, -12.0, -11.7,  0.00, -0.72,-0.243,  0.19,-0.062,  &
     &         0.414,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /

!! ---  logical flags for optional output fields

      logical :: lhlwb  = .false.
      logical :: lhlw0  = .false.
      logical :: lflxprf= .false.

!  ---  those data will be set up only once by "rlwinit"

!  ...  fluxfac, heatfac are factors for fluxes (in w/m**2) and heating
!       rates (in k/day, or k/sec set by subroutine 'rlwinit')
!       semiss0 are default surface emissivity for each bands

      real (kind=kind_phys) :: fluxfac, heatfac, semiss0(nbands)
      data semiss0(:) / nbands*1.0 /

      real (kind=kind_phys) :: tau_tbl(0:ntbl)  !< clr-sky opt dep (for cldy transfer)
      real (kind=kind_phys) :: exp_tbl(0:ntbl)  !< transmittance lookup table
      real (kind=kind_phys) :: tfn_tbl(0:ntbl)  !< tau transition function; i.e. the
                                                !< transition of planck func from mean lyr
                                                !< temp to lyr boundary temp as a func of
                                                !< opt dep. "linear in tau" method is used.

!  ---  the following variables are used for sub-column cloud scheme

      integer, parameter :: ipsdlw0 = ngptlw     ! initial permutation seed

!  ---  public accessable subprograms

      public rrtmg_lw_init, rrtmg_lw_run, rrtmg_lw_finalize, rlwinit


! ================
      contains
! ================

         subroutine rrtmg_lw_init ()
         end subroutine rrtmg_lw_init

!> \defgroup module_radlw_main GFS RRTMG Longwave Module 
!! \brief This module includes NCEP's modifications of the RRTMG-LW radiation
!! code from AER.
!!
!! The RRTMG-LW package includes three files:
!! - radlw_param.f, which contains:
!!  - module_radlw_parameters: band parameters set up
!! - radlw_datatb.f, which contains modules:
!!  - module_radlw_avplank: plank flux data
!!  - module_radlw_ref: reference temperature and pressure
!!  - module_radlw_cldprlw: cloud property coefficients
!!  - module_radlw_kgbnn: absorption coeffients for 16 bands, where nn = 01-16
!! - radlw_main.f, which contains:
!!  - rrtmg_lw_run(): the main LW radiation routine
!!  - rlwinit(): the initialization routine
!!
!!\version NCEP LW v5.1  Nov 2012 -RRTMG-LW v4.82
!!
!!\copyright  2002-2007, Atmospheric & Environmental Research, Inc. (AER).
!!  This software may be used, copied, or redistributed as long as it is
!!  not sold and this copyright notice is reproduced on each copy made.
!!  This model is provided as is without any express or implied warranties.
!!  (http://www.rtweb.aer.com/)
!! \section arg_table_rrtmg_lw_run Argument Table
!! \htmlinclude rrtmg_lw_run.html
!!
!> \section gen_lwrad RRTMG Longwave Radiation Scheme General Algorithm
!> @{
      subroutine rrtmg_lw_run                                           &
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr_co2, gasvmr_n2o,      &   !  ---  inputs
     &       gasvmr_ch4, gasvmr_o2, gasvmr_co, gasvmr_cfc11,            &
     &       gasvmr_cfc12, gasvmr_cfc22, gasvmr_ccl4,                   &
     &       icseed,aeraod,aerssa,sfemis,sfgtmp,                        &
     &       dzlyr,delpin,de_lgth,alpha,                                &
     &       npts, nlay, nlp1, lprnt, cld_cf, lslwr,                    &
     &       hlwc,topflx,sfcflx,cldtau,                                 &   !  ---  outputs
     &       HLW0,HLWB,FLXPRF,                                          &   !  ---  optional
     &       cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice,                &
     &       cld_rwp,cld_ref_rain, cld_swp, cld_ref_snow,               &
     &       cld_od, errmsg, errflg                                     &
     &     )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!     plyr (npts,nlay) : layer mean pressures (mb)                      !
!     plvl (npts,nlp1) : interface pressures (mb)                       !
!     tlyr (npts,nlay) : layer mean temperature (k)                     !
!     tlvl (npts,nlp1) : interface temperatures (k)                     !
!     qlyr (npts,nlay) : layer specific humidity (gm/gm)   *see inside  !
!     olyr (npts,nlay) : layer ozone concentration (gm/gm) *see inside  !
!     gasvmr(npts,nlay,:): atmospheric gases amount:                    !
!                       (check module_radiation_gases for definition)   !
!       gasvmr(:,:,1)  -   co2 volume mixing ratio                      !
!       gasvmr(:,:,2)  -   n2o volume mixing ratio                      !
!       gasvmr(:,:,3)  -   ch4 volume mixing ratio                      !
!       gasvmr(:,:,4)  -   o2  volume mixing ratio                      !
!       gasvmr(:,:,5)  -   co  volume mixing ratio                      !
!       gasvmr(:,:,6)  -   cfc11 volume mixing ratio                    !
!       gasvmr(:,:,7)  -   cfc12 volume mixing ratio                    !
!       gasvmr(:,:,8)  -   cfc22 volume mixing ratio                    !
!       gasvmr(:,:,9)  -   ccl4  volume mixing ratio                    !
!     clouds(npts,nlay,:): layer cloud profiles:                        !
!                       (check module_radiation_clouds for definition)  !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer in-cloud liq water path   (g/m**2)     !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)     !
!       clouds(:,:,4)  -   layer in-cloud ice water path   (g/m**2)     !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)     !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)     !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)     !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)     !
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)     !
!     icseed(npts)   : auxiliary special cloud related array            !
!                      when module variable isubclw=2, it provides      !
!                      permutation seed for each column profile that    !
!                      are used for generating random numbers.          !
!                      when isubclw /=2, it will not be used.           !
!     aerosols(npts,nlay,nbands,:) : aerosol optical properties         !
!                       (check module_radiation_aerosols for definition)!
!        (:,:,:,1)     - optical depth                                  !
!        (:,:,:,2)     - single scattering albedo                       !
!        (:,:,:,3)     - asymmetry parameter                            !
!     sfemis (npts)  : surface emissivity                               !
!     sfgtmp (npts)  : surface ground temperature (k)                   !
!     dzlyr(npts,nlay) : layer thickness (km)                           !
!     delpin(npts,nlay): layer pressure thickness (mb)                  !
!     de_lgth(npts)    : cloud decorrelation length (km)                !
!     alpha(npts,nlay) : EXP/ER cloud overlap decorrelation parameter   !
!     npts           : total number of horizontal points                !
!     nlay, nlp1     : total number of vertical layers, levels          !
!     lprnt          : cntl flag for diagnostic print out               !
!                                                                       !
!  output variables:                                                    !
!     hlwc  (npts,nlay): total sky heating rate (k/day or k/sec)        !
!     topflx(npts)     : radiation fluxes at top, component:            !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux at top (w/m2)          !
!        upfx0           - clear sky upward flux at top (w/m2)          !
!     sfcflx(npts)     : radiation fluxes at sfc, component:            !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux at sfc (w/m2)          !
!        upfx0           - clear sky upward flux at sfc (w/m2)          !
!        dnfxc           - total sky downward flux at sfc (w/m2)        !
!        dnfx0           - clear sky downward flux at sfc (w/m2)        !
!     cldtau(npts,nlay): approx 10mu band layer cloud optical depth     !
!                                                                       !
!! optional output variables:                                           !
!     hlwb(npts,nlay,nbands): spectral band total sky heating rates     !
!     hlw0  (npts,nlay): clear sky heating rate (k/day or k/sec)        !
!     flxprf(npts,nlp1): level radiative fluxes (w/m2), components:     !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux                        !
!        dnfxc           - total sky dnward flux                        !
!        upfx0           - clear sky upward flux                        !
!        dnfx0           - clear sky dnward flux                        !
!                                                                       !
!  external module variables:  (in physparam)                           !
!   ilwrgas - control flag for rare gases (ch4,n2o,o2,cfcs, etc.)       !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   ilwcliq - control flag for liq-cloud optical properties             !
!           =1: input cld liqp & reliq, hu & stamnes (1993)             !
!           =2: not used                                                !
!   ilwcice - control flag for ice-cloud optical properties             !
!           =1: input cld icep & reice, ebert & curry (1997)            !
!           =2: input cld icep & reice, streamer (1996)                 !
!           =3: input cld icep & reice, fu (1998)                       !
!   isubclw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   iovr  - cloud overlapping control flag                              !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud (used for isubclw>0 only)         !
!           =3: decorrelation-length overlap (for isubclw>0 only)       !
!           =4: exponential cloud overlap (AER)                         !
!           =5: exponential-random cloud overlap (AER)                  !
!   ivflip  - control flag for vertical index direction                 !
!           =0: vertical index from toa to surface                      !
!           =1: vertical index from surface to toa                      !
!                                                                       !
!  module parameters, control variables:                                !
!     nbands           - number of longwave spectral bands              !
!     maxgas           - maximum number of absorbing gaseous            !
!     maxxsec          - maximum number of cross-sections               !
!     ngptlw           - total number of g-point subintervals           !
!     ng##             - number of g-points in band (##=1-16)           !
!     ngb(ngptlw)      - band indices for each g-point                  !
!     bpade            - pade approximation constant (1/0.278)          !
!     nspa,nspb(nbands)- number of lower/upper ref atm's per band       !
!     delwave(nbands)  - longwave band width (wavenumbers)              !
!     ipsdlw0          - permutation seed for mcica sub-col clds        !
!                                                                       !
!  major local variables:                                               !
!     pavel  (nlay)         - layer pressures (mb)                      !
!     delp   (nlay)         - layer pressure thickness (mb)             !
!     tavel  (nlay)         - layer temperatures (k)                    !
!     tz     (0:nlay)       - level (interface) temperatures (k)        !
!     semiss (nbands)       - surface emissivity for each band          !
!     wx     (nlay,maxxsec) - cross-section molecules concentration     !
!     coldry (nlay)         - dry air column amount                     !
!                                   (1.e-20*molecules/cm**2)            !
!     cldfrc (0:nlp1)       - layer cloud fraction                      !
!     taucld (nbands,nlay)  - layer cloud optical depth for each band   !
!     cldfmc (ngptlw,nlay)  - layer cloud fraction for each g-point     !
!     tauaer (nbands,nlay)  - aerosol optical depths                    !
!     fracs  (ngptlw,nlay)  - planck fractions                          !
!     tautot (ngptlw,nlay)  - total optical depths (gaseous+aerosols)   !
!     colamt (nlay,maxgas)  - column amounts of absorbing gases         !
!                             1-maxgas are for watervapor, carbon       !
!                             dioxide, ozone, nitrous oxide, methane,   !
!                             oxigen, carbon monoxide, respectively     !
!                             (molecules/cm**2)                         !
!     pwvcm                 - column precipitable water vapor (cm)      !
!     secdiff(nbands)       - variable diffusivity angle defined as     !
!                             an exponential function of the column     !
!                             water amount in bands 2-3 and 5-9.        !
!                             this reduces the bias of several w/m2 in  !
!                             downward surface flux in high water       !
!                             profiles caused by using the constant     !
!                             diffusivity angle of 1.66.         (mji)  !
!     facij  (nlay)         - indicator of interpolation factors        !
!                             =0/1: indicate lower/higher temp & height !
!     selffac(nlay)         - scale factor for self-continuum, equals   !
!                          (w.v. density)/(atm density at 296K,1013 mb) !
!     selffrac(nlay)        - factor for temp interpolation of ref      !
!                             self-continuum data                       !
!     indself(nlay)         - index of the lower two appropriate ref    !
!                             temp for the self-continuum interpolation !
!     forfac (nlay)         - scale factor for w.v. foreign-continuum   !
!     forfrac(nlay)         - factor for temp interpolation of ref      !
!                             w.v. foreign-continuum data               !
!     indfor (nlay)         - index of the lower two appropriate ref    !
!                             temp for the foreign-continuum interp     !
!     laytrop               - tropopause layer index at which switch is !
!                             made from one conbination kew species to  !
!                             another.                                  !
!     jp(nlay),jt(nlay),jt1(nlay)                                       !
!                           - lookup table indexes                      !
!     totuflux(0:nlay)      - total-sky upward longwave flux (w/m2)     !
!     totdflux(0:nlay)      - total-sky downward longwave flux (w/m2)   !
!     htr(nlay)             - total-sky heating rate (k/day or k/sec)   !
!     totuclfl(0:nlay)      - clear-sky upward longwave flux (w/m2)     !
!     totdclfl(0:nlay)      - clear-sky downward longwave flux (w/m2)   !
!     htrcl(nlay)           - clear-sky heating rate (k/day or k/sec)   !
!     fnet    (0:nlay)      - net longwave flux (w/m2)                  !
!     fnetc   (0:nlay)      - clear-sky net longwave flux (w/m2)        !
!                                                                       !
!                                                                       !
!  ======================    end of definitions    ===================  !

!  ---  inputs:
      integer, intent(in) :: npts, nlay, nlp1
      integer, intent(in) :: icseed(npts)

      logical,  intent(in) :: lprnt

      real (kind=kind_phys), dimension(npts,nlp1), intent(in) :: plvl,  &
     &       tlvl
      real (kind=kind_phys), dimension(npts,nlay), intent(in) :: plyr,  &
     &       tlyr, qlyr, olyr, dzlyr, delpin

      real (kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_co2,&
     &     gasvmr_n2o, gasvmr_ch4, gasvmr_o2, gasvmr_co, gasvmr_cfc11,  &
     &     gasvmr_cfc12, gasvmr_cfc22, gasvmr_ccl4

      real (kind=kind_phys), dimension(npts,nlay),intent(in):: cld_cf
      real (kind=kind_phys), dimension(npts,nlay),intent(in),optional:: &
     &       cld_lwp, cld_ref_liq,  cld_iwp, cld_ref_ice,               &
     &       cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow,              &
     &       cld_od

      real (kind=kind_phys), dimension(npts), intent(in) :: sfemis,     &
     &       sfgtmp, de_lgth
      real (kind=kind_phys), dimension(npts,nlay), intent(in) :: alpha

      real (kind=kind_phys), dimension(npts,nlay,nbands),intent(in)::   &
     &       aeraod, aerssa

!mz* HWRF -- OUTPUT from mcica_subcol_lw
      real(kind=kind_phys),dimension(ngptlw,npts,nlay)  :: cldfmcl     ! Cloud fraction
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=kind_phys),dimension(ngptlw,npts,nlay)  :: ciwpmcl     ! In-cloud ice water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=kind_phys),dimension(ngptlw,npts,nlay)  :: clwpmcl     ! In-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=kind_phys),dimension(ngptlw,npts,nlay)  :: cswpmcl     ! In-cloud snow   water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=kind_phys),dimension(npts,nlay) :: relqmcl   ! Cloud water drop  effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=kind_phys),dimension(npts,nlay) :: reicmcl   ! Cloud ice  effective size (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=kind_phys),dimension(npts,nlay) :: resnmcl   ! Snow effective size (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=kind_phys),dimension(ngptlw,npts,nlay) :: taucmcl     ! In-cloud optical depth
                                                      !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=kind_phys),dimension(npts,nlay,nbands) :: tauaer      ! Aerosol optical  depth
!                                                      !    Dimensions: (ncol,nlay,nbndlw)
!mz* output from cldprmc
      integer :: ncbands        ! number of cloud  spectral bands
      real(kind=kind_phys),dimension(ngptlw,nlay) :: taucmc     ! cloud optical depth [mcica]        
                                                      !    Dimensions: (ngptlw,nlayers)    
!mz

!  ---  outputs:
      real (kind=kind_phys), dimension(npts,nlay), intent(inout) :: hlwc
      real (kind=kind_phys), dimension(npts,nlay), intent(inout) ::     &
     &       cldtau

      type (topflw_type),    dimension(npts), intent(inout) :: topflx
      type (sfcflw_type),    dimension(npts), intent(inout) :: sfcflx

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!! ---  optional outputs:
      real (kind=kind_phys), dimension(npts,nlay,nbands),optional,      &
     &       intent(inout) :: hlwb
      real (kind=kind_phys), dimension(npts,nlay),       optional,      &
     &       intent(inout) :: hlw0
      type (proflw_type),    dimension(npts,nlp1),       optional,      &
     &       intent(inout) :: flxprf
      logical, intent(in) :: lslwr

!  ---  locals:
! mz* - Add height of each layer for exponential-random cloud overlap
! This will be derived below from the dzlyr in each layer
      real (kind=kind_phys), dimension( npts,nlay )  ::   hgt
      real (kind=kind_phys) :: dzsum

      real (kind=kind_phys), dimension(0:nlp1) :: cldfrc

      real (kind=kind_phys), dimension(0:nlay) :: totuflux, totdflux,   &
     &       totuclfl, totdclfl, tz

      real (kind=kind_phys), dimension(nlay)   :: htr, htrcl

      real (kind=kind_phys), dimension(nlay)   :: pavel, tavel, delp,   &
     &       clwp, ciwp, relw, reiw, cda1, cda2, cda3, cda4,            &
     &       coldry, colbrd, h2ovmr, o3vmr, fac00, fac01, fac10, fac11, &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       scaleminorn2, temcol, dz

!mz*
      real(kind=rb),dimension(0:nlay,nbands) :: planklay,planklev 
      real(kind=rb),dimension(0:nlay) :: pz 

!      real(kind=rb) :: plankbnd(nbndlw)
      real (kind=kind_phys), dimension(nbands,0:nlay) :: pklev, pklay

      real (kind=kind_phys), dimension(nlay,nbands) :: htrb
      real (kind=kind_phys), dimension(nbands,nlay) :: taucld, tauaer
      real (kind=kind_phys), dimension(nbands,npts,nlay) :: taucld3
      real (kind=kind_phys), dimension(ngptlw,nlay) :: fracs, tautot
      real (kind=kind_phys), dimension(nlay,ngptlw) :: fracs_r
!mz rtrnmc_mcica
      real (kind=kind_phys), dimension(nlay,ngptlw) :: taut 
!mz* Atmosphere/clouds - cldprop
      real(kind=kind_phys), dimension(ngptlw,nlay) :: cldfmc,    &
     &                                                cldfmc_save  ! cloud fraction [mcica]
                                                                   !    Dimensions: (ngptlw,nlay)
      real(kind=kind_phys), dimension(ngptlw,nlay) :: ciwpmc       ! in-cloud ice water path [mcica]
                                                                   !    Dimensions: (ngptlw,nlay)
      real(kind=kind_phys), dimension(ngptlw,nlay) :: clwpmc       ! in-cloud liquid water path [mcica]
                                                                   !    Dimensions: (ngptlw,nlay)
      real(kind=kind_phys), dimension(ngptlw,nlay) :: cswpmc       ! in-cloud snow path [mcica]
                                                                   !    Dimensions: (ngptlw,nlay)
      real(kind=kind_phys), dimension(nlay) :: relqmc              ! liquid particle effective radius (microns)
                                                                   !    Dimensions: (nlay)
      real(kind=kind_phys), dimension(nlay) :: reicmc              ! ice particle effective size (microns)
                                                                   !    Dimensions: (nlay)
      real(kind=kind_phys), dimension(nlay) :: resnmc              ! snow effective size (microns)
                                                                   !    Dimensions: (nlay)


      real (kind=kind_phys), dimension(nbands) :: semiss, secdiff

!  ---  column amount of absorbing gases:
!       (:,m) m = 1-h2o, 2-co2, 3-o3, 4-n2o, 5-ch4, 6-o2, 7-co
      real (kind=kind_phys) :: colamt(nlay,maxgas)

!  ---  column cfc cross-section amounts:
!       (:,m) m = 1-ccl4, 2-cfc11, 3-cfc12, 4-cfc22
      real (kind=kind_phys) :: wx(nlay,maxxsec)

!  ---  reference ratios of binary species parameter in lower atmosphere:
!       (:,m,:) m = 1-h2o/co2, 2-h2o/o3, 3-h2o/n2o, 4-h2o/ch4, 5-n2o/co2, 6-o3/co2
      real (kind=kind_phys) :: rfrate(nlay,nrates,2)

      real (kind=kind_phys) :: tem0, tem1, tem2, pwvcm, summol, stemp,  &
     &                         delgth
      real (kind=kind_phys), dimension(nlay) :: alph

      integer, dimension(npts) :: ipseed
      integer, dimension(nlay) :: jp, jt, jt1, indself, indfor, indminor
      integer                  :: laytrop, iplon, i, j, k, k1
      ! mz* added local arrays for RRTMG
      integer                  :: irng, permuteseed,ig
      integer                  :: inflglw, iceflglw, liqflglw
      logical :: lcf1
      integer :: istart              ! beginning band of calculation
      integer :: iend                ! ending band of calculation
      integer :: iout                ! output option flag (inactive)


!
!===> ... begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!mz*
! For passing in cloud physical properties; cloud optics parameterized
! in RRTMG:
      inflglw = 2
      iceflglw = 3
      liqflglw = 1
      istart = 1
      iend = 16
      iout = 0

!
      if (.not. lslwr) return

!  --- ...  initialization

      lhlwb  = present ( hlwb )
      lhlw0  = present ( hlw0 )
      lflxprf= present ( flxprf )

      colamt(:,:) = f_zero
      cldtau(:,:) = f_zero

!! --- check for optional input arguments, depending on cloud method
      if (ilwcliq > 0) then    ! use prognostic cloud method
        if ( .not.present(cld_lwp) .or. .not.present(cld_ref_liq) .or.  &
     &       .not.present(cld_iwp) .or. .not.present(cld_ref_ice) .or.  &
     &       .not.present(cld_rwp) .or. .not.present(cld_ref_rain) .or. &
     &       .not.present(cld_swp) .or. .not.present(cld_ref_snow)) then
          write(errmsg,'(*(a))')                                        &
     &               'Logic error: ilwcliq>0 requires the following',   &
     &               ' optional arguments to be present:',              &
     &               ' cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice,',    &
     &               ' cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow'
          errflg = 1
          return
        end if
      else                     ! use diagnostic cloud method
        if ( .not.present(cld_od) ) then
          write(errmsg,'(*(a))')                                        &
     &               'Logic error: ilwcliq<=0 requires the following',  &
     &               ' optional argument to be present: cld_od'
          errflg = 1
          return
        end if
      endif                    ! end if_ilwcliq

!> -# Change random number seed value for each radiation invocation
!!    (isubclw =1 or 2).

      if     ( isubclw == 1 ) then     ! advance prescribed permutation seed
        do i = 1, npts
          ipseed(i) = ipsdlw0 + i
        enddo
      elseif ( isubclw == 2 ) then     ! use input array of permutaion seeds
        do i = 1, npts
          ipseed(i) = icseed(i)
        enddo
      endif

!     if ( lprnt ) then
!       print *,'  In rrtmg_lw, isubclw, ipsdlw0,ipseed =',             &
!    &          isubclw, ipsdlw0, ipseed
!     endif

!  --- ...  loop over horizontal npts profiles

      lab_do_iplon : do iplon = 1, npts

!> -# Read surface emissivity.
        if (sfemis(iplon) > eps .and. sfemis(iplon) <= 1.0) then  ! input surface emissivity
          do j = 1, nbands
            semiss(j) = sfemis(iplon)
          enddo
        else                                                      ! use default values
          do j = 1, nbands
            semiss(j) = semiss0(j)
          enddo
        endif

        stemp = sfgtmp(iplon)          ! surface ground temp
        if (iovr == 3) delgth= de_lgth(iplon)    ! clouds decorr-length

! mz*: HWRF
        if (iovr == 4 ) then

!Add layer height needed for exponential (icld=4) and
! exponential-random (icld=5) overlap options  

         !iplon = 1
         irng = 0
         permuteseed = 150

!mz* Derive height 
         dzsum =0.0
         do k = 1,nlay
         hgt(iplon,k)= dzsum+0.5*dzlyr(iplon,k)*1000.   !km->m
         dzsum =  dzsum+ dzlyr(iplon,k)*1000.   
         enddo

! Zero out cloud optical properties here; not used when passing physical properties
! to radiation and taucld is calculated in radiation 
            do k = 1, nlay 
               do j = 1, nbands
                  taucld3(j,iplon,k) = 0.0
               enddo
            enddo

          call mcica_subcol_lw(1, iplon, nlay, iovr, permuteseed,       &
     &                 irng, plyr, hgt,                                 &
     &                 cld_cf, cld_iwp, cld_lwp,cld_swp,                &
     &                 cld_ref_ice, cld_ref_liq,                        &
     &                 cld_ref_snow, taucld3,                           &
     &                 cldfmcl,                                         &  !--output
     &                 ciwpmcl, clwpmcl, cswpmcl, reicmcl, relqmcl,     &
     &                 resnmcl, taucmcl)     

       endif
!mz* end

!> -# Prepare atmospheric profile for use in rrtm.
!           the vertical index of internal array is from surface to top

!  --- ...  molecular amounts are input or converted to volume mixing ratio
!           and later then converted to molecular amount (molec/cm2) by the
!           dry air column coldry (in molec/cm2) which is calculated from the
!           layer pressure thickness (in mb), based on the hydrostatic equation
!  --- ...  and includes a correction to account for h2o in the layer.

        if (ivflip == 0) then       ! input from toa to sfc

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd
          tz(0) = tlvl(iplon,nlp1)

          do k = 1, nlay
            k1 = nlp1 - k
            pavel(k)= plyr(iplon,k1)
            delp(k) = delpin(iplon,k1)
            tavel(k)= tlyr(iplon,k1)
            tz(k)   = tlvl(iplon,k1)
            dz(k)   = dzlyr(iplon,k1)
            if (iovr == 4 .or. iovr == 5) alph(k) = alpha(iplon,k) ! alpha decorrelation

!> -# Set absorber amount for h2o, co2, and o3.

!test use
!           h2ovmr(k)= max(f_zero,qlyr(iplon,k1)*amdw)                  ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qlyr(iplon,k1))                       ! input vol mixing ratio
!           o3vmr (k)= max(f_zero,olyr(iplon,k1))                       ! input vol mixing ratio
!ncep model use
            h2ovmr(k)= max(f_zero,qlyr(iplon,k1)                        &
     &                           *amdw/(f_one-qlyr(iplon,k1)))          ! input specific humidity
            o3vmr (k)= max(f_zero,olyr(iplon,k1)*amdo3)                 ! input mass mixing ratio

!  --- ...  tem0 is the molecular weight of moist air
            tem0 = (f_one - h2ovmr(k))*con_amd + h2ovmr(k)*con_amw
            coldry(k) = tem2*delp(k) / (tem1*tem0*(f_one+h2ovmr(k)))
            temcol(k) = 1.0e-12 * coldry(k)

            colamt(k,1) = max(f_zero,    coldry(k)*h2ovmr(k))          ! h2o
            colamt(k,2) = max(temcol(k), coldry(k)*gasvmr_co2(iplon,k1)) ! co2
            colamt(k,3) = max(temcol(k), coldry(k)*o3vmr(k))           ! o3
          enddo

!> -# Set up column amount for rare gases n2o,ch4,o2,co,ccl4,cf11,cf12,
!!    cf22, convert from volume mixing ratio to molec/cm2 based on
!!    coldry (scaled to 1.0e-20).

          if (ilwrgas > 0) then
            do k = 1, nlay
              k1 = nlp1 - k
              colamt(k,4)=max(temcol(k), coldry(k)*gasvmr_n2o(iplon,k1))  ! n2o
              colamt(k,5)=max(temcol(k), coldry(k)*gasvmr_ch4(iplon,k1))  ! ch4
              colamt(k,6)=max(f_zero,    coldry(k)*gasvmr_o2(iplon,k1))   ! o2
              colamt(k,7)=max(f_zero,    coldry(k)*gasvmr_co(iplon,k1))   ! co

              wx(k,1) = max( f_zero, coldry(k)*gasvmr_ccl4(iplon,k1) )    ! ccl4
              wx(k,2) = max( f_zero, coldry(k)*gasvmr_cfc11(iplon,k1) )   ! cf11
              wx(k,3) = max( f_zero, coldry(k)*gasvmr_cfc12(iplon,k1) )   ! cf12
              wx(k,4) = max( f_zero, coldry(k)*gasvmr_cfc22(iplon,k1) )   ! cf22
            enddo
          else
            do k = 1, nlay
              colamt(k,4) = f_zero     ! n2o
              colamt(k,5) = f_zero     ! ch4
              colamt(k,6) = f_zero     ! o2
              colamt(k,7) = f_zero     ! co

              wx(k,1) = f_zero
              wx(k,2) = f_zero
              wx(k,3) = f_zero
              wx(k,4) = f_zero
            enddo
          endif

!> -# Set aerosol optical properties.

          do k = 1, nlay
            k1 = nlp1 - k
            do j = 1, nbands
              tauaer(j,k) = aeraod(iplon,k1,j)                          &
     &                    * (f_one - aerssa(iplon,k1,j))
            enddo
          enddo

!> -# Read cloud optical properties.
          if (ilwcliq > 0) then    ! use prognostic cloud method
            do k = 1, nlay
              k1 = nlp1 - k
              cldfrc(k)= cld_cf(iplon,k1)
              clwp(k)  = cld_lwp(iplon,k1)
              relw(k)  = cld_ref_liq(iplon,k1)
              ciwp(k)  = cld_iwp(iplon,k1)
              reiw(k)  = cld_ref_ice(iplon,k1)
              cda1(k)  = cld_rwp(iplon,k1)
              cda2(k)  = cld_ref_rain(iplon,k1)
              cda3(k)  = cld_swp(iplon,k1)
              cda4(k)  = cld_ref_snow(iplon,k1)
            enddo
            ! HWRF RRMTG
            if (iovr == 4) then   !mz  HWRF 
               do k = 1, nlay
                  k1 = nlp1 - k
               do ig = 1, ngptlw
                   cldfmc(ig,k) = cldfmcl(ig,iplon,k1)
                   taucmc(ig,k) = taucmcl(ig,iplon,k1)
                   ciwpmc(ig,k) = ciwpmcl(ig,iplon,k1)
                   clwpmc(ig,k) = clwpmcl(ig,iplon,k1)
              !mz     cswpmc(ig,k) = cswpmcl(ig,iplon,k1)
                   cswpmc(ig,k) = 0.0
               enddo
                   reicmc(k) = reicmcl(iplon,k1)
                   relqmc(k) = relqmcl(iplon,k1)
                   resnmc(k) = resnmcl(iplon,k1)
               enddo
            endif
          else                       ! use diagnostic cloud method
            do k = 1, nlay
              k1 = nlp1 - k
              cldfrc(k)= cld_cf(iplon,k1)
              cda1(k)  = cld_od(iplon,k1)
            enddo
          endif                      ! end if_ilwcliq

          cldfrc(0)    = f_one       ! padding value only
          cldfrc(nlp1) = f_zero      ! padding value only

!> -# Compute precipitable water vapor for diffusivity angle adjustments.

          tem1 = f_zero
          tem2 = f_zero
          do k = 1, nlay
            tem1 = tem1 + coldry(k) + colamt(k,1)
            tem2 = tem2 + colamt(k,1)
          enddo

          tem0 = 10.0 * tem2 / (amdw * tem1 * con_g)
          pwvcm = tem0 * plvl(iplon,nlp1)

        else                        ! input from sfc to toa

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd
          tz(0) = tlvl(iplon,1)

          do k = 1, nlay
            pavel(k)= plyr(iplon,k)
            delp(k) = delpin(iplon,k)
            tavel(k)= tlyr(iplon,k)
            tz(k)   = tlvl(iplon,k+1)
            dz(k)   = dzlyr(iplon,k)
            if (iovr == 4 .or. iovr == 5) alph(k) = alpha(iplon,k) ! alpha decorrelation

!  --- ...  set absorber amount
!test use
!           h2ovmr(k)= max(f_zero,qlyr(iplon,k)*amdw)                   ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qlyr(iplon,k))                        ! input vol mixing ratio
!           o3vmr (k)= max(f_zero,olyr(iplon,k))                        ! input vol mixing ratio
!ncep model use
            h2ovmr(k)= max(f_zero,qlyr(iplon,k)                         &
     &                           *amdw/(f_one-qlyr(iplon,k)))           ! input specific humidity
            o3vmr (k)= max(f_zero,olyr(iplon,k)*amdo3)                  ! input mass mixing ratio

!  --- ...  tem0 is the molecular weight of moist air
            tem0 = (f_one - h2ovmr(k))*con_amd + h2ovmr(k)*con_amw
            coldry(k) = tem2*delp(k) / (tem1*tem0*(f_one+h2ovmr(k)))
            temcol(k) = 1.0e-12 * coldry(k)

            colamt(k,1) = max(f_zero,    coldry(k)*h2ovmr(k))          ! h2o
            colamt(k,2) = max(temcol(k), coldry(k)*gasvmr_co2(iplon,k))! co2
            colamt(k,3) = max(temcol(k), coldry(k)*o3vmr(k))           ! o3
          enddo

!  --- ...  set up col amount for rare gases, convert from volume mixing ratio
!           to molec/cm2 based on coldry (scaled to 1.0e-20)

          if (ilwrgas > 0) then
            do k = 1, nlay
              colamt(k,4)=max(temcol(k), coldry(k)*gasvmr_n2o(iplon,k))  ! n2o
              colamt(k,5)=max(temcol(k), coldry(k)*gasvmr_ch4(iplon,k))  ! ch4
              colamt(k,6)=max(f_zero,    coldry(k)*gasvmr_o2(iplon,k))   ! o2
              colamt(k,7)=max(f_zero,    coldry(k)*gasvmr_co(iplon,k))   ! co

              wx(k,1) = max( f_zero, coldry(k)*gasvmr_ccl4(iplon,k) )    ! ccl4
              wx(k,2) = max( f_zero, coldry(k)*gasvmr_cfc11(iplon,k) )   ! cf11
              wx(k,3) = max( f_zero, coldry(k)*gasvmr_cfc12(iplon,k) )   ! cf12
              wx(k,4) = max( f_zero, coldry(k)*gasvmr_cfc22(iplon,k) )   ! cf22
            enddo
          else
            do k = 1, nlay
              colamt(k,4) = f_zero     ! n2o
              colamt(k,5) = f_zero     ! ch4
              colamt(k,6) = f_zero     ! o2
              colamt(k,7) = f_zero     ! co

              wx(k,1) = f_zero
              wx(k,2) = f_zero
              wx(k,3) = f_zero
              wx(k,4) = f_zero
            enddo
          endif

!  --- ...  set aerosol optical properties

          do j = 1, nbands
            do k = 1, nlay
              tauaer(j,k) = aeraod(iplon,k,j)                           &
     &                    * (f_one - aerssa(iplon,k,j))
            enddo
          enddo

          if (ilwcliq > 0) then    ! use prognostic cloud method
            do k = 1, nlay
              cldfrc(k)= cld_cf(iplon,k)
              clwp(k)  = cld_lwp(iplon,k)
              relw(k)  = cld_ref_liq(iplon,k)
              ciwp(k)  = cld_iwp(iplon,k)
              reiw(k)  = cld_ref_ice(iplon,k)
              cda1(k)  = cld_rwp(iplon,k)
              cda2(k)  = cld_ref_rain(iplon,k)
              cda3(k)  = cld_swp(iplon,k)
              cda4(k)  = cld_ref_snow(iplon,k)
            enddo
            if (iovr == 4) then
!mz* Move incoming GCM cloud arrays to RRTMG cloud arrays.
!For GCM input, incoming reicmcl is defined based on selected 
!ice parameterization (inflglw)
            do k = 1, nlay
            do ig = 1, ngptlw
               cldfmc(ig,k) = cldfmcl(ig,iplon,k)
               taucmc(ig,k) = taucmcl(ig,iplon,k)
               ciwpmc(ig,k) = ciwpmcl(ig,iplon,k)
               clwpmc(ig,k) = clwpmcl(ig,iplon,k)
              !mz cswpmc(ig,k) = cswpmcl(ig,iplon,k)
               cswpmc(ig,k) = 0.0
            enddo
               reicmc(k) = reicmcl(iplon,k)
               relqmc(k) = relqmcl(iplon,k)
               resnmc(k) = resnmcl(iplon,k)
            enddo
            endif
          else                       ! use diagnostic cloud method
            do k = 1, nlay
              cldfrc(k)= cld_cf(iplon,k)
              cda1(k)  = cld_od(iplon,k)
            enddo
          endif                      ! end if_ilwcliq

          cldfrc(0)    = f_one       ! padding value only
          cldfrc(nlp1) = f_zero      ! padding value only

!  --- ...  compute precipitable water vapor for diffusivity angle adjustments

          tem1 = f_zero
          tem2 = f_zero
          do k = 1, nlay
            tem1 = tem1 + coldry(k) + colamt(k,1)
            tem2 = tem2 + colamt(k,1)
          enddo

          tem0 = 10.0 * tem2 / (amdw * tem1 * con_g)
          pwvcm = tem0 * plvl(iplon,1)

        endif                       ! if_ivflip

!> -# Compute column amount for broadening gases.

        do k = 1, nlay
          summol = f_zero
          do i = 2, maxgas
            summol = summol + colamt(k,i)
          enddo
          colbrd(k) = coldry(k) - summol
        enddo

!> -# Compute diffusivity angle adjustments.

        tem1 = 1.80
        tem2 = 1.50
        do j = 1, nbands
          if (j==1 .or. j==4 .or. j==10) then
            secdiff(j) = 1.66
          else
            secdiff(j) = min( tem1, max( tem2,                          &
     &                   a0(j)+a1(j)*exp(a2(j)*pwvcm) ))
          endif
        enddo

!     if (lprnt) then
!      print *,'  coldry',coldry
!      print *,' wx(*,1) ',(wx(k,1),k=1,NLAY)
!      print *,' wx(*,2) ',(wx(k,2),k=1,NLAY)
!      print *,' wx(*,3) ',(wx(k,3),k=1,NLAY)
!      print *,' wx(*,4) ',(wx(k,4),k=1,NLAY)
!      print *,' iplon ',iplon
!      print *,'  pavel ',pavel
!      print *,'  delp ',delp
!      print *,'  tavel ',tavel
!      print *,'  tz ',tz
!      print *,' h2ovmr ',h2ovmr
!      print *,' o3vmr ',o3vmr
!     endif

!> -# For cloudy atmosphere, call cldprop() to set cloud optical
!!    properties.

        lcf1 = .false.
        lab_do_k0 : do k = 1, nlay
          if ( cldfrc(k) > eps ) then
            lcf1 = .true.
            exit lab_do_k0
          endif
        enddo  lab_do_k0

        if ( lcf1 ) then

          !mz* for HWRF, save cldfmc with mcica
          if (iovr == 4) then
               do k = 1, nlay
               do ig = 1, ngptlw
                  cldfmc_save(ig,k)=cldfmc (ig,k)
               enddo
               enddo
          endif

          call cldprop                                                  &
!  ---  inputs:
     &     ( cldfrc,clwp,relw,ciwp,reiw,cda1,cda2,cda3,cda4,            &
     &       nlay, nlp1, ipseed(iplon), dz, delgth, iovr, alph,         &
!  ---  outputs:
     &       cldfmc, taucld                                             &
     &     )

          if (iovr == 4) then
          !mz for HWRF, still using mcica cldfmc
               do k = 1, nlay
               do ig = 1, ngptlw
                  cldfmc(ig,k)=cldfmc_save(ig,k)
               enddo
               enddo
          endif

!  --- ...  save computed layer cloud optical depth for output
!           rrtm band-7 is apprx 10mu channel (or use spectral mean of bands 6-8)

          if (ivflip == 0) then       ! input from toa to sfc
            do k = 1, nlay
              k1 = nlp1 - k
              cldtau(iplon,k1) = taucld( 7,k)
            enddo
          else                        ! input from sfc to toa
            do k = 1, nlay
              cldtau(iplon,k) = taucld( 7,k)
            enddo
          endif                       ! end if_ivflip_block

        else
          cldfmc = f_zero
          taucld = f_zero
        endif

!mz* HWRF: calculate taucmc with mcica
        if (iovr == 4) then 
          call cldprmc(nlay, inflglw, iceflglw, liqflglw,               &
     &                 cldfmc, ciwpmc,                                  &
     &                 clwpmc, cswpmc, reicmc, relqmc, resnmc,          &
     &                 ncbands, taucmc, errmsg, errflg)
          ! return immediately if cldprmc throws an error
          if (errflg/=0) return
        endif

!     if (lprnt) then
!      print *,' after cldprop'
!      print *,' clwp',clwp
!      print *,' ciwp',ciwp
!      print *,' relw',relw
!      print *,' reiw',reiw
!      print *,' taucl',cda1
!      print *,' cldfrac',cldfrc
!     endif

!> -# Calling setcoef() to compute various coefficients needed in
!!    radiative transfer calculations.
        call setcoef                                                    &
!  ---  inputs:
     &     ( pavel,tavel,tz,stemp,h2ovmr,colamt,coldry,colbrd,          &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       laytrop,pklay,pklev,jp,jt,jt1,                             &
     &       rfrate,fac00,fac01,fac10,fac11,                            &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor                 &
     &     )

!     if (lprnt) then
!      print *,'laytrop',laytrop
!      print *,'colh2o',(colamt(k,1),k=1,NLAY)
!      print *,'colco2',(colamt(k,2),k=1,NLAY)
!      print *,'colo3', (colamt(k,3),k=1,NLAY)
!      print *,'coln2o',(colamt(k,4),k=1,NLAY)
!      print *,'colch4',(colamt(k,5),k=1,NLAY)
!      print *,'fac00',fac00
!      print *,'fac01',fac01
!      print *,'fac10',fac10
!      print *,'fac11',fac11
!      print *,'jp',jp
!      print *,'jt',jt
!      print *,'jt1',jt1
!      print *,'selffac',selffac
!      print *,'selffrac',selffrac
!      print *,'indself',indself
!      print *,'forfac',forfac
!      print *,'forfrac',forfrac
!      print *,'indfor',indfor
!     endif

!> -# Call taumol() to calculte the gaseous optical depths and Plank
!! fractions for each longwave spectral band.

        call taumol                                                     &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay,                                                      &
!  ---  outputs:
     &       fracs, tautot                                              &
     &     )

!     if (lprnt) then
!     print *,' after taumol'
!     do k = 1, nlay
!       write(6,121) k
!121    format(' k =',i3,5x,'FRACS')
!       write(6,122) (fracs(j,k),j=1,ngptlw)
!122    format(10e14.7)
!       write(6,123) k
!123    format(' k =',i3,5x,'TAUTOT')
!       write(6,122) (tautot(j,k),j=1,ngptlw)
!     enddo
!     endif

!> -# Call the radiative transfer routine based on cloud scheme
!!    selection. Compute the upward/downward radiative fluxes, and
!!    heating rates for both clear or cloudy atmosphere.
!!\n  - call rtrn(): clouds are assumed as randomly overlaping in a
!!                   vertical column
!!\n  - call rtrnmr(): clouds are assumed as in maximum-randomly
!!                     overlaping in a vertical column;
!!\n  - call rtrnmc(): clouds are treated with the mcica stochastic
!!                     approach.

        if (isubclw <= 0) then

          if (iovr <= 0) then

            call rtrn                                                   &
!  ---  inputs:
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdiff,nlay,nlp1,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
     &     )

          else

            call rtrnmr                                                 &
!  ---  inputs:
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdiff,nlay,nlp1,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
     &     )

          endif   ! end if_iovr_block

        else

          call rtrnmc                                                   &
!  ---  inputs:
     &     ( semiss,delp,cldfmc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdiff,nlay,nlp1,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
     &     )

        endif   ! end if_isubclw_block

!> -# Save outputs.

        topflx(iplon)%upfxc = totuflux(nlay)
        topflx(iplon)%upfx0 = totuclfl(nlay)

        sfcflx(iplon)%upfxc = totuflux(0)
        sfcflx(iplon)%upfx0 = totuclfl(0)
        sfcflx(iplon)%dnfxc = totdflux(0)
        sfcflx(iplon)%dnfx0 = totdclfl(0)

        if (ivflip == 0) then       ! output from toa to sfc

!! --- ...  optional fluxes
          if ( lflxprf ) then
            do k = 0, nlay
              k1 = nlp1 - k
              flxprf(iplon,k1)%upfxc = totuflux(k)
              flxprf(iplon,k1)%dnfxc = totdflux(k)
              flxprf(iplon,k1)%upfx0 = totuclfl(k)
              flxprf(iplon,k1)%dnfx0 = totdclfl(k)
            enddo
          endif

          do k = 1, nlay
            k1 = nlp1 - k
            hlwc(iplon,k1) = htr(k)
          enddo

!! --- ...  optional clear sky heating rate
          if ( lhlw0 ) then
            do k = 1, nlay
              k1 = nlp1 - k
              hlw0(iplon,k1) = htrcl(k)
            enddo
          endif

!! --- ...  optional spectral band heating rate
          if ( lhlwb ) then
            do j = 1, nbands
            do k = 1, nlay
              k1 = nlp1 - k
              hlwb(iplon,k1,j) = htrb(k,j)
            enddo
            enddo
          endif

        else                        ! output from sfc to toa

!! --- ...  optional fluxes
          if ( lflxprf ) then
            do k = 0, nlay
              flxprf(iplon,k+1)%upfxc = totuflux(k)
              flxprf(iplon,k+1)%dnfxc = totdflux(k)
              flxprf(iplon,k+1)%upfx0 = totuclfl(k)
              flxprf(iplon,k+1)%dnfx0 = totdclfl(k)
            enddo
          endif

          do k = 1, nlay
            hlwc(iplon,k) = htr(k)
          enddo

!! --- ...  optional clear sky heating rate
          if ( lhlw0 ) then
            do k = 1, nlay
              hlw0(iplon,k) = htrcl(k)
            enddo
          endif

!! --- ...  optional spectral band heating rate
          if ( lhlwb ) then
            do j = 1, nbands
            do k = 1, nlay
              hlwb(iplon,k,j) = htrb(k,j)
            enddo
            enddo
          endif

        endif                       ! if_ivflip

      enddo  lab_do_iplon

!...................................
      end subroutine rrtmg_lw_run
!-----------------------------------
!> @}
      subroutine rrtmg_lw_finalize ()
      end subroutine rrtmg_lw_finalize 



!> \ingroup module_radlw_main
!> \brief This subroutine performs calculations necessary for the initialization
!! of the longwave model, which includes non-varying model variables, conversion
!! factors, and look-up tables  
!!
!! Lookup tables are computed for use in the lw
!! radiative transfer, and input absorption coefficient data for each
!! spectral band are reduced from 256 g-point intervals to 140.
!!\param me        print control for parallel process
!!\section rlwinit_gen rlwinit General Algorithm
!! @{
      subroutine rlwinit                                                &
     &     ( me ) !  ---  inputs
!  ---  outputs: (none)

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  initialize non-varying module variables, conversion factors,!
! and look-up tables.                                                   !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!    me       - print control for parallel process                      !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  external module variables:  (in physparam)                            !
!   ilwrate - heating rate unit selections                              !
!           =1: output in k/day                                         !
!           =2: output in k/second                                      !
!   ilwrgas - control flag for rare gases (ch4,n2o,o2,cfcs, etc.)       !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   ilwcliq - liquid cloud optical properties contrl flag               !
!           =0: input cloud opt depth from diagnostic scheme            !
!           >0: input cwp,rew, and other cloud content parameters       !
!   isubclw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   icldflg - cloud scheme control flag                                 !
!           =0: diagnostic scheme gives cloud tau, omiga, and g.        !
!           =1: prognostic scheme gives cloud liq/ice path, etc.        !
!   iovr  - clouds vertical overlapping control flag                    !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud (isubcol>0 only)                  !
!           =3: decorrelation-length overlap (for isubclw>0 only)       !
!           =4: exponential overlap cloud
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:       michael j. iacono; july, 1998                !
!  first revision for ncar ccm:               september, 1998           !
!  second revision for rrtm_v3.0:             september, 2002           !
!                                                                       !
!  this subroutine performs calculations necessary for the initialization
!  of the longwave model.  lookup tables are computed for use in the lw !
!  radiative transfer, and input absorption coefficient data for each   !
!  spectral band are reduced from 256 g-point intervals to 140.         !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
! definitions:                                                          !
!   arrays for 10000-point look-up tables:                              !
!   tau_tbl - clear-sky optical depth (used in cloudy radiative transfer!
!   exp_tbl - exponential lookup table for tansmittance                 !
!   tfn_tbl - tau transition function; i.e. the transition of the Planck!
!             function from that for the mean layer temperature to that !
!             for the layer boundary temperature as a function of optical
!             depth. the "linear in tau" method is used to make the table
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: me

!  ---  outputs: none

!  ---  locals:
      real (kind=kind_phys), parameter :: expeps = 1.e-20

      real (kind=kind_phys) :: tfn, pival, explimit

      integer               :: i

!
!===> ... begin here
!
      if ( iovr<0 .or. iovr>4 ) then
        print *,'  *** Error in specification of cloud overlap flag',   &
     &          ' IOVR=',iovr,' in RLWINIT !!'
        stop
      elseif ( (iovr==2 .or. iovr==3) .and. isubclw==0 ) then
        if (me == 0) then
          print *,'  *** IOVR=',iovr,' is not available for',           &
     &          ' ISUBCLW=0 setting!!'
          print *,'      The program uses maximum/random overlap',      &
     &          ' instead.'
        endif

        iovr = 1
      endif

      if (me == 0) then
        print *,' - Using AER Longwave Radiation, Version: ', VTAGLW

        if (ilwrgas > 0) then
          print *,'   --- Include rare gases N2O, CH4, O2, CFCs ',      &
     &            'absorptions in LW'
        else
          print *,'   --- Rare gases effect is NOT included in LW'
        endif

        if ( isubclw == 0 ) then
          print *,'   --- Using standard grid average clouds, no ',     &
     &            'sub-column clouds approximation applied'
        elseif ( isubclw == 1 ) then
          print *,'   --- Using MCICA sub-colum clouds approximation ', &
     &            'with a prescribed sequence of permutaion seeds'
        elseif ( isubclw == 2 ) then
          print *,'   --- Using MCICA sub-colum clouds approximation ', &
     &            'with provided input array of permutation seeds'
        else
          print *,'  *** Error in specification of sub-column cloud ',  &
     &            ' control flag isubclw =',isubclw,' !!'
          stop
        endif
      endif

!> -# Check cloud flags for consistency.

      if ((icldflg == 0 .and. ilwcliq /= 0) .or.                        &
     &    (icldflg == 1 .and. ilwcliq == 0)) then
        print *,'  *** Model cloud scheme inconsistent with LW',        &
     &          ' radiation cloud radiative property setup !!'
        stop
      endif

!> -# Setup default surface emissivity for each band.

      semiss0(:) = f_one

!> -# Setup constant factors for flux and heating rate
!! the 1.0e-2 is to convert pressure from mb to \f$N/m^2\f$.

      pival = 2.0 * asin(f_one)
      fluxfac = pival * 2.0d4
!     fluxfac = 62831.85307179586                   ! = 2 * pi * 1.0e4

      if (ilwrate == 1) then
!       heatfac = 8.4391
!       heatfac = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
        heatfac = con_g * 864.0 / con_cp            !   (in k/day)
      else
        heatfac = con_g * 1.0e-2 / con_cp           !   (in k/second)
      endif

!> -# Compute lookup tables for transmittance, tau transition
!! function, and clear sky tau (for the cloudy sky radiative
!! transfer).  tau is computed as a function of the tau
!! transition function, transmittance is calculated as a
!! function of tau, and the tau transition function is
!! calculated using the linear in tau formulation at values of
!! tau above 0.01.  tf is approximated as tau/6 for tau < 0.01.
!! all tables are computed at intervals of 0.001.  the inverse
!! of the constant used in the pade approximation to the tau
!! transition function is set to b.

      tau_tbl(0) = f_zero
      exp_tbl(0) = f_one
      tfn_tbl(0) = f_zero

      tau_tbl(ntbl) = 1.e10
      exp_tbl(ntbl) = expeps
      tfn_tbl(ntbl) = f_one

      explimit = aint( -log(tiny(exp_tbl(0))) )

      do i = 1, ntbl-1
!org    tfn = float(i) / float(ntbl)
!org    tau_tbl(i) = bpade * tfn / (f_one - tfn)
        tfn = real(i, kind_phys) / real(ntbl-i, kind_phys)
        tau_tbl(i) = bpade * tfn
        if (tau_tbl(i) >= explimit) then
          exp_tbl(i) = expeps
        else
          exp_tbl(i) = exp( -tau_tbl(i) )
        endif

        if (tau_tbl(i) < 0.06) then
          tfn_tbl(i) = tau_tbl(i) / 6.0
        else
          tfn_tbl(i) = f_one - 2.0*( (f_one / tau_tbl(i))               &
     &               - ( exp_tbl(i) / (f_one - exp_tbl(i)) ) )
        endif
      enddo

!...................................
      end subroutine rlwinit
!! @}
!-----------------------------------


!>\ingroup module_radlw_main
!> \brief This subroutine computes the cloud optical depth(s) for each cloudy
!! layer and g-point interval.
!!\param cfrac           layer cloud fraction
!!\n     ---  for  ilwcliq > 0 (prognostic cloud scheme)  - - -
!!\param cliqp           layer in-cloud liq water path (\f$g/m^2\f$)
!!\param reliq           mean eff radius for liq cloud (micron)
!!\param cicep           layer in-cloud ice water path (\f$g/m^2\f$)
!!\param reice           mean eff radius for ice cloud (micron)
!!\param cdat1           layer rain drop water path (\f$g/m^2\f$)
!!\param cdat2           effective radius for rain drop (micron)
!!\param cdat3           layer snow flake water path(\f$g/m^2\f$)
!!\param cdat4           mean effective radius for snow flake(micron)
!!\n     ---  for ilwcliq = 0  (diagnostic cloud scheme)  - - -
!!\param cliqp           not used
!!\param cicep           not used
!!\param reliq           not used
!!\param reice           not used
!!\param cdat1           layer cloud optical depth
!!\param cdat2           layer cloud single scattering albedo
!!\param cdat3           layer cloud asymmetry factor
!!\param cdat4           optional use
!!\param nlay            number of layer number
!!\param nlp1            number of veritcal levels
!!\param ipseed          permutation seed for generating random numbers (isubclw>0)
!!\param dz              layer thickness (km) 
!!\param de_lgth         layer cloud decorrelation length (km)  
!!\param alpha           EXP/ER cloud overlap decorrelation parameter
!!\param cldfmc          cloud fraction for each sub-column
!!\param taucld          cloud optical depth for bands (non-mcica)
!!\section gen_cldprop cldprop General Algorithm
!> @{
      subroutine cldprop                                                &
     &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     & !  ---  inputs
     &       nlay, nlp1, ipseed, dz, de_lgth, iovr, alpha,              &
     &       cldfmc, taucld                                             & !  ---  outputs
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the cloud optical depth(s) for each cloudy layer    !
! and g-point interval.                                                 !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       -size- !
!    cfrac - real, layer cloud fraction                          0:nlp1 !
!        .....  for ilwcliq > 0  (prognostic cloud sckeme)  - - -       !
!    cliqp - real, layer in-cloud liq water path (g/m**2)          nlay !
!    reliq - real, mean eff radius for liq cloud (micron)          nlay !
!    cicep - real, layer in-cloud ice water path (g/m**2)          nlay !
!    reice - real, mean eff radius for ice cloud (micron)          nlay !
!    cdat1 - real, layer rain drop water path  (g/m**2)            nlay !
!    cdat2 - real, effective radius for rain drop (microm)         nlay !
!    cdat3 - real, layer snow flake water path (g/m**2)            nlay !
!    cdat4 - real, effective radius for snow flakes (micron)       nlay !
!        .....  for ilwcliq = 0  (diagnostic cloud sckeme)  - - -       !
!    cdat1 - real, input cloud optical depth                       nlay !
!    cdat2 - real, layer cloud single scattering albedo            nlay !
!    cdat3 - real, layer cloud asymmetry factor                    nlay !
!    cdat4 - real, optional use                                    nlay !
!    cliqp - not used                                              nlay !
!    reliq - not used                                              nlay !
!    cicep - not used                                              nlay !
!    reice - not used                                              nlay !
!                                                                       !
!    dz     - real, layer thickness (km)                           nlay !
!    de_lgth- real, layer cloud decorrelation length (km)             1 !
!    alpha  - real, EXP/ER decorrelation parameter                 nlay !
!    nlay  - integer, number of vertical layers                      1  !
!    nlp1  - integer, number of vertical levels                      1  !
!    ipseed- permutation seed for generating random numbers (isubclw>0) !
!                                                                       !
!  outputs:                                                             !
!    cldfmc - real, cloud fraction for each sub-column       ngptlw*nlay!
!    taucld - real, cld opt depth for bands (non-mcica)      nbands*nlay!
!                                                                       !
!  explanation of the method for each value of ilwcliq, and ilwcice.    !
!    set up in module "module_radlw_cntr_para"                          !
!                                                                       !
!     ilwcliq=0  : input cloud optical property (tau, ssa, asy).        !
!                  (used for diagnostic cloud method)                   !
!     ilwcliq>0  : input cloud liq/ice path and effective radius, also  !
!                  require the user of 'ilwcice' to specify the method  !
!                  used to compute aborption due to water/ice parts.    !
!  ...................................................................  !
!                                                                       !
!     ilwcliq=1:   the water droplet effective radius (microns) is input!
!                  and the opt depths due to water clouds are computed  !
!                  as in hu and stamnes, j., clim., 6, 728-742, (1993). !
!                  the values for absorption coefficients appropriate for
!                  the spectral bands in rrtm have been obtained for a  !
!                  range of effective radii by an averaging procedure   !
!                  based on the work of j. pinto (private communication).
!                  linear interpolation is used to get the absorption   !
!                  coefficients for the input effective radius.         !
!                                                                       !
!     ilwcice=1:   the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in ebert and curry, jgr, 97,  !
!                  3831-3836 (1992).  the spectral regions in this work !
!                  have been matched with the spectral bands in rrtm to !
!                  as great an extent as possible:                      !
!                     e&c 1      ib = 5      rrtm bands 9-16            !
!                     e&c 2      ib = 4      rrtm bands 6-8             !
!                     e&c 3      ib = 3      rrtm bands 3-5             !
!                     e&c 4      ib = 2      rrtm band 2                !
!                     e&c 5      ib = 1      rrtm band 1                !
!     ilwcice=2:   the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in rt code, streamer v3.0     !
!                  (ref: key j., streamer user's guide, cooperative     !
!                  institute for meteorological satellite studies, 2001,!
!                  96 pp.) valid range of values for re are between 5.0 !
!                  and 131.0 micron.                                    !
!     ilwcice=3:   the ice generalized effective size (dge) is input and!
!                  the optical properties, are calculated as in q. fu,  !
!                  j. climate, (1998). q. fu provided high resolution   !
!                  tales which were appropriately averaged for the bands!
!                  in rrtm_lw. linear interpolation is used to get the  !
!                  coeff from the stored tables. valid range of values  !
!                  for deg are between 5.0 and 140.0 micron.            !
!                                                                       !
!  other cloud control module variables:                                !
!     isubclw =0: standard cloud scheme, no sub-col cloud approximation !
!             >0: mcica sub-col cloud scheme using ipseed as permutation!
!                 seed for generating rundom numbers                    !
!                                                                       !
!  ======================  end of description block  =================  !
!
      use module_radlw_cldprlw

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ipseed, iovr

      real (kind=kind_phys), dimension(0:nlp1), intent(in) :: cfrac
      real (kind=kind_phys), dimension(nlay),   intent(in) :: cliqp,    &
     &       reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4, dz
      real (kind=kind_phys),                    intent(in) :: de_lgth
      real (kind=kind_phys), dimension(nlay),   intent(in) :: alpha

!  ---  outputs:
      real (kind=kind_phys), dimension(ngptlw,nlay),intent(out):: cldfmc
      real (kind=kind_phys), dimension(nbands,nlay),intent(out):: taucld

!  ---  locals:
      real (kind=kind_phys), dimension(nbands) :: tauliq, tauice
      real (kind=kind_phys), dimension(nlay)   :: cldf

      real (kind=kind_phys) :: dgeice, factor, fint, tauran, tausnw,    &
     &       cldliq, refliq, cldice, refice

      logical :: lcloudy(ngptlw,nlay)
      integer :: ia, ib, ig, k, index

!
!===> ...  begin here
!
      do k = 1, nlay
        do ib = 1, nbands
          taucld(ib,k) = f_zero
        enddo
      enddo

      do k = 1, nlay
        do ig = 1, ngptlw
          cldfmc(ig,k) = f_zero
        enddo
      enddo

!> -# Compute cloud radiative properties for a cloudy column:
!!\n  - Compute cloud radiative properties for rain and snow (tauran,tausnw)
!!\n  - Calculation of absorption coefficients due to water clouds(tauliq)
!!\n  - Calculation of absorption coefficients due to ice clouds (tauice).
!!\n  - For prognostic cloud scheme: sum up the cloud optical property:
!!\n    \f$ taucld=tauice+tauliq+tauran+tausnw \f$

!  --- ...  compute cloud radiative properties for a cloudy column

      lab_if_ilwcliq : if (ilwcliq > 0) then

        lab_do_k : do k = 1, nlay
          lab_if_cld : if (cfrac(k) > cldmin) then

            tauran = absrain * cdat1(k)                      ! ncar formula
!!          tausnw = abssnow1 * cdat3(k)                     ! ncar formula
!  ---  if use fu's formula it needs to be normalized by snow density
!       !not use snow density = 0.1 g/cm**3 = 0.1 g/(mu * m**2)
!       use ice density = 0.9167 g/cm**3 = 0.9167 g/(mu * m**2)
!       factor 1.5396=8/(3*sqrt(3)) converts reff to generalized ice particle size
!       use newer factor value 1.0315
!       1/(0.9167*1.0315) = 1.05756
            if (cdat3(k)>f_zero .and. cdat4(k)>10.0_kind_phys) then
              tausnw = abssnow0*1.05756*cdat3(k)/cdat4(k)      ! fu's formula
            else
              tausnw = f_zero
            endif

            cldliq = cliqp(k)
            cldice = cicep(k)
!           refliq = max(2.5e0, min(60.0e0, reliq(k) ))
!           refice = max(5.0e0, reice(k) )
            refliq = reliq(k)
            refice = reice(k)

!  --- ...  calculation of absorption coefficients due to water clouds.

            if ( cldliq <= f_zero ) then
              do ib = 1, nbands
                tauliq(ib) = f_zero
              enddo
            else
              if ( ilwcliq == 1 ) then

                factor = refliq - 1.5
                index  = max( 1, min( 57, int( factor ) ))
                fint   = factor - float(index)

                do ib = 1, nbands
                  tauliq(ib) = max(f_zero, cldliq*(absliq1(index,ib)    &
     &              + fint*(absliq1(index+1,ib)-absliq1(index,ib)) ))
                enddo
              endif   ! end if_ilwcliq_block
            endif   ! end if_cldliq_block

!  --- ...  calculation of absorption coefficients due to ice clouds.

            if ( cldice <= f_zero ) then
              do ib = 1, nbands
                tauice(ib) = f_zero
              enddo
            else

!  --- ...  ebert and curry approach for all particle sizes though somewhat
!           unjustified for large ice particles

              if ( ilwcice == 1 ) then
                refice = min(130.0, max(13.0, real(refice) ))

                do ib = 1, nbands
                  ia = ipat(ib)             ! eb_&_c band index for ice cloud coeff
                  tauice(ib) = max(f_zero, cldice*(absice1(1,ia)        &
     &                         + absice1(2,ia)/refice) )
                enddo

!  --- ...  streamer approach for ice effective radius between 5.0 and 131.0 microns
!           and ebert and curry approach for ice eff radius greater than 131.0 microns.
!           no smoothing between the transition of the two methods.

              elseif ( ilwcice == 2 ) then

                factor = (refice - 2.0) / 3.0
                index  = max( 1, min( 42, int( factor ) ))
                fint   = factor - float(index)

                do ib = 1, nbands
                  tauice(ib) = max(f_zero, cldice*(absice2(index,ib)    &
     &              + fint*(absice2(index+1,ib) - absice2(index,ib)) ))
                enddo

!  --- ...  fu's approach for ice effective radius between 4.8 and 135 microns
!           (generalized effective size from 5 to 140 microns)

              elseif ( ilwcice == 3 ) then

!               dgeice = max(5.0, 1.5396*refice)              ! v4.4 value
                dgeice = max(5.0, 1.0315*refice)              ! v4.71 value
                factor = (dgeice - 2.0) / 3.0
                index  = max( 1, min( 45, int( factor ) ))
                fint   = factor - float(index)

                do ib = 1, nbands
                  tauice(ib) = max(f_zero, cldice*(absice3(index,ib)    &
     &              + fint*(absice3(index+1,ib) - absice3(index,ib)) ))
                enddo

              endif   ! end if_ilwcice_block
            endif   ! end if_cldice_block

            do ib = 1, nbands
              taucld(ib,k) = tauice(ib) + tauliq(ib) + tauran + tausnw
            enddo

          endif  lab_if_cld
        enddo  lab_do_k

      else  lab_if_ilwcliq

        do k = 1, nlay
          if (cfrac(k) > cldmin) then
            do ib = 1, nbands
              taucld(ib,k) = cdat1(k)
            enddo
          endif
        enddo

      endif  lab_if_ilwcliq

!> -# if physparam::isubclw > 0, call mcica_subcol() to distribute
!!    cloud properties to each g-point.

      if ( isubclw > 0 ) then      ! mcica sub-col clouds approx
        do k = 1, nlay
          if ( cfrac(k) < cldmin ) then
            cldf(k) = f_zero
          else
            cldf(k) = cfrac(k)
          endif
        enddo

!  --- ...  call sub-column cloud generator

!mz*
      if (iovr .ne. 4) then
        call mcica_subcol                                               &
!  ---  inputs:
     &     ( cldf, nlay, ipseed, dz, de_lgth, alpha,                    &
!  ---  output:
     &       lcloudy                                                    &
     &     )

        do k = 1, nlay
          do ig = 1, ngptlw
            if ( lcloudy(ig,k) ) then
              cldfmc(ig,k) = f_one
            else
              cldfmc(ig,k) = f_zero
            endif
          enddo
        enddo
      endif  !iovr

      endif   ! end if_isubclw_block

      return
! ..................................
      end subroutine cldprop
! ----------------------------------
!> @}

!>\ingroup module_radlw_main
!>\brief This suroutine computes sub-colum cloud profile flag array.
!!\param cldf        layer cloud fraction
!!\param nlay        number of model vertical layers
!!\param ipseed      permute seed for random num generator
!!\param dz          layer thickness
!!\param de_lgth     layer cloud decorrelation length (km)
!!\param alpha       EXP/ER cloud overlap decorrelation parameter
!!\param lcloudy     sub-colum cloud profile flag array
!!\section mcica_subcol_gen mcica_subcol General Algorithm
!! @{
      subroutine mcica_subcol                                           &
     &    ( cldf, nlay, ipseed, dz, de_lgth, alpha,                     & !  ---  inputs
     &      lcloudy                                                     & !  ---  outputs
     &    )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                size !
!   cldf    - real, layer cloud fraction                           nlay !
!   nlay    - integer, number of model vertical layers               1  !
!   ipseed  - integer, permute seed for random num generator         1  !
!    ** note : if the cloud generator is called multiple times, need    !
!              to permute the seed between each call; if between calls  !
!              for lw and sw, use values differ by the number of g-pts. !
!   dz      - real, layer thickness (km)                           nlay !
!   de_lgth - real, layer cloud decorrelation length (km)            1  !
!    alpha  - real, EXP/ER decorrelation parameter                 nlay !
!                                                                       !
!  output variables:                                                    !
!   lcloudy - logical, sub-colum cloud profile flag array    ngptlw*nlay!
!                                                                       !
!  other control flags from module variables:                           !
!     iovr    : control flag for cloud overlapping method               !
!                 =0:random; =1:maximum/random: =2:maximum; =3:decorr   !
!                                                                       !
!  =====================    end of definitions    ====================  !

      implicit none

!  ---  inputs:
      integer, intent(in) :: nlay, ipseed

      real (kind=kind_phys), dimension(nlay), intent(in) :: cldf, dz
      real (kind=kind_phys),                  intent(in) :: de_lgth
      real (kind=kind_phys), dimension(nlay), intent(in) :: alpha

!  ---  outputs:
      logical, dimension(ngptlw,nlay), intent(out) :: lcloudy

!  ---  locals:
      real (kind=kind_phys) :: cdfunc(ngptlw,nlay), rand1d(ngptlw),     &
     &       rand2d(nlay*ngptlw), tem1, fac_lcf(nlay),                  &
     &       cdfun2(ngptlw,nlay)

      type (random_stat) :: stat          ! for thread safe random generator

      integer :: k, n, k1
!
!===> ...  begin here
!
!> -# Call random_setseed() to advance randum number generator by ipseed values.

      call random_setseed                                               &
!  ---  inputs:
     &    ( ipseed,                                                     &
!  ---  outputs:
     &      stat                                                        &
     &    )

!> -# Sub-column set up according to overlapping assumption:
!!  - For random overlap, pick a random value at every level 
!!  - For max-random overlap, pick a random value at every level
!!  - For maximum overlap, pick same random numebr at every level

      select case ( iovr )

        case( 0 )        ! random overlap, pick a random value at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptlw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(n,k) = rand2d(k1)
            enddo
          enddo

        case( 1 )        ! max-ran overlap

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptlw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(n,k) = rand2d(k1)
            enddo
          enddo

!  ---  first pick a random number for bottom (or top) layer.
!       then walk up the column: (aer's code)
!       if layer below is cloudy, use the same rand num in the layer below
!       if layer below is clear,  use a new random number

!  ---  from bottom up
          do k = 2, nlay
            k1 = k - 1
            tem1 = f_one - cldf(k1)

            do n = 1, ngptlw
              if ( cdfunc(n,k1) > tem1 ) then
                cdfunc(n,k) = cdfunc(n,k1)
              else
                cdfunc(n,k) = cdfunc(n,k) * tem1
              endif
            enddo
          enddo

!  ---  or walk down the column: (if use original author's method)
!       if layer above is cloudy, use the same rand num in the layer above
!       if layer above is clear,  use a new random number

!  ---  from top down
!         do k = nlay-1, 1, -1
!           k1 = k + 1
!           tem1 = f_one - cldf(k1)

!           do n = 1, ngptlw
!             if ( cdfunc(n,k1) > tem1 ) then
!               cdfunc(n,k) = cdfunc(n,k1)
!             else
!               cdfunc(n,k) = cdfunc(n,k) * tem1
!             endif
!           enddo
!         enddo

        case( 2 )        !<  - For maximum overlap, pick same random numebr at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand1d, stat )

          do n = 1, ngptlw
            tem1 = rand1d(n)

            do k = 1, nlay
              cdfunc(n,k) = tem1
            enddo
          enddo

        case( 3 )        ! decorrelation length overlap

!  ---  compute overlapping factors based on layer midpoint distances
!       and decorrelation depths

          do k = nlay, 2, -1
            fac_lcf(k) = exp( -0.5 * (dz(k)+dz(k-1)) / de_lgth )
          enddo

!  ---  setup 2 sets of random numbers

          call random_number ( rand2d, stat )

          k1 = 0
          do k = 1, nlay
            do n = 1, ngptlw
              k1 = k1 + 1
              cdfunc(n,k) = rand2d(k1)
            enddo
          enddo

          call random_number ( rand2d, stat )

          k1 = 0
          do k = 1, nlay
            do n = 1, ngptlw
              k1 = k1 + 1
              cdfun2(n,k) = rand2d(k1)
            enddo
          enddo

!  ---  then working from the top down:
!       if a random number (from an independent set -cdfun2) is smaller then the
!       scale factor: use the upper layer's number,  otherwise use a new random
!       number (keep the original assigned one).

          do k = nlay-1, 1, -1
            k1 = k + 1

            do n = 1, ngptlw
              if ( cdfun2(n,k) <= fac_lcf(k1) ) then
                cdfunc(n,k) = cdfunc(n,k1)
              endif
            enddo
          enddo

        case( 4:5 )        ! exponential and exponential-random cloud overlap

!  ---  Use previously derived decorrelation parameter, alpha, to specify
!       the exponenential transition of cloud correlation in the vertical column.
!
!       For exponential cloud overlap, the correlation is applied across layers
!       without regard to the configuration of clear and cloudy layers.

!       For exponential-random cloud overlap, a new exponential transition is 
!       performed within each group of adjacent cloudy layers and blocks of 
!       cloudy layers with clear layers between them are correlated randomly. 
!
!       NOTE: The code below is identical for case (4) and (5) because the 
!       distinction in the vertical correlation between EXP and ER is already 
!       built into the specification of alpha (in subroutine get_alpha_exp). 

!  ---  setup 2 sets of random numbers

          call random_number ( rand2d, stat )

          k1 = 0
          do k = 1, nlay
            do n = 1, ngptlw
              k1 = k1 + 1
              cdfunc(n,k) = rand2d(k1)
            enddo
          enddo

          call random_number ( rand2d, stat )

          k1 = 0
          do k = 1, nlay
            do n = 1, ngptlw
              k1 = k1 + 1
              cdfun2(n,k) = rand2d(k1)
            enddo
          enddo

!  ---  then working upward from the surface:
!       if a random number (from an independent set: cdfun2) is smaller than 
!       alpha, then use the previous layer's number, otherwise use a new random
!       number (keep the originally assigned one in cdfunc for that layer).

          do k = 2, nlay
            k1 = k - 1
            do n = 1, ngptlw
              if ( cdfun2(n,k) < alpha(k) ) then
                   cdfunc(n,k) = cdfunc(n,k1)
              endif
            enddo
          enddo

      end select

!> -# Generate subcolumns for homogeneous clouds.

      do k = 1, nlay
        tem1 = f_one - cldf(k)

        do n = 1, ngptlw
          lcloudy(n,k) = cdfunc(n,k) >= tem1
        enddo
      enddo

      return
! ..................................
      end subroutine mcica_subcol
!! @}
! ----------------------------------

!>\ingroup module_radlw_main
!> This subroutine computes various coefficients needed in radiative
!! transfer calculations.
!!\param pavel           layer pressure (mb)
!!\param tavel           layer temperature (K)
!!\param tz              level(interface) temperatures (K)
!!\param stemp           surface ground temperature (K)
!!\param h2ovmr          layer w.v. volumn mixing ratio (kg/kg)
!!\param colamt           column amounts of absorbing gases.
!! 2nd indices range: 1-maxgas, for watervapor,carbon dioxide, ozone,
!! nitrous oxide, methane,oxigen, carbon monoxide,etc. \f$(mol/cm^2)\f$
!!\param coldry          dry air column amount
!!\param colbrd          column amount of broadening gases
!!\param nlay            total number of vertical layers
!!\param nlp1            total number of vertical levels
!!\param laytrop         tropopause layer index (unitless)
!!\param pklay           integrated planck func at lay temp
!!\param pklev           integrated planck func at lev temp
!!\param jp              indices of lower reference pressure
!!\param jt, jt1         indices of lower reference temperatures
!!\param rfrate          ref ratios of binary species param
!!\n                     (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,
!!                                4-h2o/ch4,5-n2o/co2,6-o3/co2
!!\n                     (:,:,n)n=1,2: the rates of ref press at
!!                                the 2 sides of the layer
!!\param fac00,fac01,fac10,fac11           factors multiply the reference ks, i,j=0/1 for
!!                       lower/higher of the 2 appropriate temperatures
!!                       and altitudes.
!!\param selffac         scale factor for w. v. self-continuum equals
!!                       (w. v. density)/(atmospheric density at 296k and 1013 mb)
!!\param selffrac        factor for temperature interpolation of
!!                       reference w. v. self-continuum data
!!\param indself         index of lower ref temp for selffac
!!\param forfac          scale factor for w. v. foreign-continuum
!!\param forfrac         factor for temperature interpolation of
!!                       reference w.v. foreign-continuum data
!!\param indfor          index of lower ref temp for forfac
!!\param minorfrac       factor for minor gases
!!\param scaleminor,scaleminorn2         scale factors for minor gases
!!\param indminor        index of lower ref temp for minor gases
!>\section setcoef_gen setcoef General Algorithm
!> @{
      subroutine setcoef                                                &
     &     ( pavel,tavel,tz,stemp,h2ovmr,colamt,coldry,colbrd,          & !  ---  inputs:
     &       nlay, nlp1,                                                &
     &       laytrop,pklay,pklev,jp,jt,jt1,                             & !  ---  outputs:
     &       rfrate,fac00,fac01,fac10,fac11,                            &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor                 &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute various coefficients needed in radiative transfer   !
!    calculations.                                                      !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       -size- !
!   pavel     - real, layer pressures (mb)                         nlay !
!   tavel     - real, layer temperatures (k)                       nlay !
!   tz        - real, level (interface) temperatures (k)         0:nlay !
!   stemp     - real, surface ground temperature (k)                1   !
!   h2ovmr    - real, layer w.v. volum mixing ratio (kg/kg)        nlay !
!   colamt    - real, column amounts of absorbing gases      nlay*maxgas!
!                 2nd indices range: 1-maxgas, for watervapor,          !
!                 carbon dioxide, ozone, nitrous oxide, methane,        !
!                 oxigen, carbon monoxide,etc. (molecules/cm**2)        !
!   coldry    - real, dry air column amount                        nlay !
!   colbrd    - real, column amount of broadening gases            nlay !
!   nlay/nlp1 - integer, total number of vertical layers, levels    1   !
!                                                                       !
!  outputs:                                                             !
!   laytrop   - integer, tropopause layer index (unitless)          1   !
!   pklay     - real, integrated planck func at lay temp   nbands*0:nlay!
!   pklev     - real, integrated planck func at lev temp   nbands*0:nlay!
!   jp        - real, indices of lower reference pressure          nlay !
!   jt, jt1   - real, indices of lower reference temperatures      nlay !
!   rfrate    - real, ref ratios of binary species param   nlay*nrates*2!
!     (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,4-h2o/ch4,5-n2o/co2,6-o3/co2!
!     (:,:,n)n=1,2: the rates of ref press at the 2 sides of the layer  !
!   facij     - real, factors multiply the reference ks,           nlay !
!                 i,j=0/1 for lower/higher of the 2 appropriate         !
!                 temperatures and altitudes.                           !
!   selffac   - real, scale factor for w. v. self-continuum        nlay !
!                 equals (w. v. density)/(atmospheric density           !
!                 at 296k and 1013 mb)                                  !
!   selffrac  - real, factor for temperature interpolation of      nlay !
!                 reference w. v. self-continuum data                   !
!   indself   - integer, index of lower ref temp for selffac       nlay !
!   forfac    - real, scale factor for w. v. foreign-continuum     nlay !
!   forfrac   - real, factor for temperature interpolation of      nlay !
!                 reference w.v. foreign-continuum data                 !
!   indfor    - integer, index of lower ref temp for forfac        nlay !
!   minorfrac - real, factor for minor gases                       nlay !
!   scaleminor,scaleminorn2                                             !
!             - real, scale factors for minor gases                nlay !
!   indminor  - integer, index of lower ref temp for minor gases   nlay !
!                                                                       !
!  ======================    end of definitions    ===================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nlay,maxgas),intent(in):: colamt
      real (kind=kind_phys), dimension(0:nlay),     intent(in):: tz

      real (kind=kind_phys), dimension(nlay), intent(in) :: pavel,      &
     &       tavel, h2ovmr, coldry, colbrd

      real (kind=kind_phys), intent(in) :: stemp

!  ---  outputs:
      integer, dimension(nlay), intent(out) :: jp, jt, jt1, indself,    &
     &       indfor, indminor

      integer, intent(out) :: laytrop

      real (kind=kind_phys), dimension(nlay,nrates,2), intent(out) ::   &
     &       rfrate
      real (kind=kind_phys), dimension(nbands,0:nlay), intent(out) ::   &
     &       pklev, pklay

      real (kind=kind_phys), dimension(nlay),          intent(out) ::   &
     &       fac00, fac01, fac10, fac11, selffac, selffrac, forfac,     &
     &       forfrac, minorfrac, scaleminor, scaleminorn2

!  ---  locals:
      real (kind=kind_phys) :: tlvlfr, tlyrfr, plog, fp, ft, ft1,       &
     &       tem1, tem2

      integer :: i, k, jp1, indlev, indlay
!
!===> ... begin here
!
!> -# Calculate information needed by the radiative transfer routine
!! that is specific to this atmosphere, especially some of the
!! coefficients and indices needed to compute the optical depths
!! by interpolating data from stored reference atmospheres.

      indlay = min(180, max(1, int(stemp-159.0) ))
      indlev = min(180, max(1, int(tz(0)-159.0) ))
      tlyrfr = stemp - int(stemp)
      tlvlfr = tz(0) - int(tz(0))
      do i = 1, nbands
        tem1 = totplnk(indlay+1,i) - totplnk(indlay,i)
        tem2 = totplnk(indlev+1,i) - totplnk(indlev,i)
        pklay(i,0) = delwave(i) * (totplnk(indlay,i) + tlyrfr*tem1)
        pklev(i,0) = delwave(i) * (totplnk(indlev,i) + tlvlfr*tem2)
      enddo

!  --- ...  begin layer loop
!> -# Calculate the integrated Planck functions for each band at the
!! surface, level, and layer temperatures.

      laytrop = 0

      do k = 1, nlay

        indlay = min(180, max(1, int(tavel(k)-159.0) ))
        tlyrfr = tavel(k) - int(tavel(k))

        indlev = min(180, max(1, int(tz(k)-159.0) ))
        tlvlfr = tz(k) - int(tz(k))

!  --- ...  begin spectral band loop

        do i = 1, nbands
          pklay(i,k) = delwave(i) * (totplnk(indlay,i) + tlyrfr         &
     &               * (totplnk(indlay+1,i) - totplnk(indlay,i)) )
          pklev(i,k) = delwave(i) * (totplnk(indlev,i) + tlvlfr         &
     &               * (totplnk(indlev+1,i) - totplnk(indlev,i)) )
        enddo

!> -# Find the two reference pressures on either side of the
!! layer pressure. store them in jp and jp1. store in fp the
!! fraction of the difference (in ln(pressure)) between these
!! two values that the layer pressure lies.

        plog = log(pavel(k))
        jp(k)= max(1, min(58, int(36.0 - 5.0*(plog+0.04)) ))
        jp1  = jp(k) + 1
!  --- ...  limit pressure extrapolation at the top
        fp   = max(f_zero, min(f_one, 5.0*(preflog(jp(k))-plog) ))
!org    fp   = 5.0 * (preflog(jp(k)) - plog)

!> -# Determine, for each reference pressure (jp and jp1), which
!! reference temperature (these are different for each
!! reference pressure) is nearest the layer temperature but does
!! not exceed it. store these indices in jt and jt1, resp.
!! store in ft (resp. ft1) the fraction of the way between jt
!! (jt1) and the next highest reference temperature that the
!! layer temperature falls.

        tem1 = (tavel(k)-tref(jp(k))) / 15.0
        tem2 = (tavel(k)-tref(jp1  )) / 15.0
        jt (k) = max(1, min(4, int(3.0 + tem1) ))
        jt1(k) = max(1, min(4, int(3.0 + tem2) ))
!  --- ...  restrict extrapolation ranges by limiting abs(det t) < 37.5 deg
        ft  = max(-0.5, min(1.5, tem1 - float(jt (k) - 3) ))
        ft1 = max(-0.5, min(1.5, tem2 - float(jt1(k) - 3) ))
!org    ft  = tem1 - float(jt (k) - 3)
!org    ft1 = tem2 - float(jt1(k) - 3)

!> -# We have now isolated the layer ln pressure and temperature,
!! between two reference pressures and two reference temperatures
!!(for each reference pressure).  we multiply the pressure
!! fraction fp with the appropriate temperature fractions to get
!! the factors that will be needed for the interpolation that yields
!! the optical depths (performed in routines taugbn for band n).

        tem1 = f_one - fp
        fac10(k) = tem1 * ft
        fac00(k) = tem1 * (f_one - ft)
        fac11(k) = fp * ft1
        fac01(k) = fp * (f_one - ft1)

        forfac(k) = pavel(k)*stpfac / (tavel(k)*(1.0 + h2ovmr(k)))
        selffac(k) = h2ovmr(k) * forfac(k)

!> -# Set up factors needed to separately include the minor gases
!! in the calculation of absorption coefficient.

        scaleminor(k) = pavel(k) / tavel(k)
        scaleminorn2(k) = (pavel(k) / tavel(k))                         &
     &                  * (colbrd(k)/(coldry(k) + colamt(k,1)))
        tem1 = (tavel(k) - 180.8) / 7.2
        indminor(k) = min(18, max(1, int(tem1)))
        minorfrac(k) = tem1 - float(indminor(k))

!> -# If the pressure is less than ~100mb, perform a different
!! set of species interpolations.

        if (plog > 4.56) then

          laytrop =  laytrop + 1

          tem1 = (332.0 - tavel(k)) / 36.0
          indfor(k) = min(2, max(1, int(tem1)))
          forfrac(k) = tem1 - float(indfor(k))

!> -# Set up factors needed to separately include the water vapor
!! self-continuum in the calculation of absorption coefficient.

          tem1 = (tavel(k) - 188.0) / 7.2
          indself(k) = min(9, max(1, int(tem1)-7))
          selffrac(k) = tem1 - float(indself(k) + 7)

!> -# Setup reference ratio to be used in calculation of binary
!! species parameter in lower atmosphere.

          rfrate(k,1,1) = chi_mls(1,jp(k)) / chi_mls(2,jp(k))
          rfrate(k,1,2) = chi_mls(1,jp(k)+1) / chi_mls(2,jp(k)+1)

          rfrate(k,2,1) = chi_mls(1,jp(k)) / chi_mls(3,jp(k))
          rfrate(k,2,2) = chi_mls(1,jp(k)+1) / chi_mls(3,jp(k)+1)

          rfrate(k,3,1) = chi_mls(1,jp(k)) / chi_mls(4,jp(k))
          rfrate(k,3,2) = chi_mls(1,jp(k)+1) / chi_mls(4,jp(k)+1)

          rfrate(k,4,1) = chi_mls(1,jp(k)) / chi_mls(6,jp(k))
          rfrate(k,4,2) = chi_mls(1,jp(k)+1) / chi_mls(6,jp(k)+1)

          rfrate(k,5,1) = chi_mls(4,jp(k)) / chi_mls(2,jp(k))
          rfrate(k,5,2) = chi_mls(4,jp(k)+1) / chi_mls(2,jp(k)+1)

        else

          tem1 = (tavel(k) - 188.0) / 36.0
          indfor(k) = 3
          forfrac(k) = tem1 - f_one

          indself(k) = 0
          selffrac(k) = f_zero

!> -# Setup reference ratio to be used in calculation of binary
!! species parameter in upper atmosphere.

          rfrate(k,1,1) = chi_mls(1,jp(k)) / chi_mls(2,jp(k))
          rfrate(k,1,2) = chi_mls(1,jp(k)+1) / chi_mls(2,jp(k)+1)

          rfrate(k,6,1) = chi_mls(3,jp(k)) / chi_mls(2,jp(k))
          rfrate(k,6,2) = chi_mls(3,jp(k)+1) / chi_mls(2,jp(k)+1)

        endif

!> -# Rescale \a selffac and \a forfac for use in taumol.

        selffac(k) = colamt(k,1) * selffac(k)
        forfac(k)  = colamt(k,1) * forfac(k)

      enddo   ! end do_k layer loop

      return
! ..................................
      end subroutine setcoef
!> @}
! ----------------------------------

!>\ingroup module_radlw_main
!> This subroutine computes the upward/downward radiative fluxes, and
!! heating rates for both clear or cloudy atmosphere. Clouds assumed as
!! randomly overlaping in a vertical column.
!!\brief Original Code Description: this program calculates the upward
!! fluxes, downward fluxes, and heating rates for an arbitrary clear or
!! cloudy atmosphere. The input to this program is the atmospheric
!! profile, all Planck function information, and the cloud fraction by
!! layer. A variable diffusivity angle (secdif) is used for the angle
!! integration. Bands 2-3 and 5-9 use a value for secdif that varies
!! from 1.50 to 1.80 as a function of the column water vapor, and other
!! bands use a value of 1.66. The gaussian weight appropriate to this
!! angle (wtdiff =0.5) is applied here. Note that use of the emissivity
!! angle for the flux integration can cause errors of 1 to 4 \f$W/m^2\f$
!! within cloudy layers. Clouds are treated with a random cloud overlap
!! method.
!!\param semiss      lw surface emissivity
!!\param delp        layer pressure thickness (mb)
!!\param cldfrc      layer cloud fraction
!!\param taucld      layer cloud opt depth
!!\param tautot      total optical depth (gas+aerosols)
!!\param pklay       integrated planck function at lay temp
!!\param pklev       integrated planck func at lev temp
!!\param fracs       planck fractions
!!\param secdif      secant of diffusivity angle
!!\param nlay        number of vertical layers
!!\param nlp1        number of vertical levels (interfaces)
!!\param totuflux    total sky upward flux \f$(w/m^2)\f$
!!\param totdflux    total sky downward flux \f$(w/m^2)\f$
!!\param htr         total sky heating rate (k/sec or k/day)
!!\param totuclfl    clear sky upward flux \f$(w/m^2)\f$
!!\param totdclfl    clear sky downward flux \f$(w/m^2)\f$
!!\param htrcl       clear sky heating rate (k/sec or k/day)
!!\param htrb        spectral band lw heating rate (k/day)
!>\section gen_rtrn rtrn General Algorithm
!! @{
! ----------------------------------
      subroutine rtrn                                                   &
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              & !  ---  inputs
     &       fracs,secdif, nlay,nlp1,                                   &
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       & !  ---  outputs
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are assumed as     !
! randomly overlaping in a vertical colum.                              !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfrc  - real, layer cloud fraction                         0:nlp1 !
!   taucld  - real, layer cloud opt depth                    nbands,nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw,nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw,nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    efclrfr- real, effective clear sky fraction (1-efcldfr)       nlay !
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck fn         1  !
!    totfac - real, gas+cld pade factor, used for planck fn          1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas only            nlay!
!    totsrcu- real, upwd source radiance due to gas+cld             nlay!
!    gassrcd- real, dnwd source radiance due to gas only             1  !
!    totsrcd- real, dnwd source radiance due to gas+cld              1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all Planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 W/m2 within cloudy layers.                !
!  clouds are treated with a random cloud overlap method.               !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(0:nlp1), intent(in) :: cldfrc
      real (kind=kind_phys), dimension(nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(nbands,nlay),intent(in):: taucld
      real (kind=kind_phys), dimension(ngptlw,nlay),intent(in):: fracs, &
     &       tautot

      real (kind=kind_phys), dimension(nbands,0:nlay), intent(in) ::    &
     &       pklev, pklay

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(nlay,nbands),intent(out) :: htrb

      real (kind=kind_phys), dimension(0:nlay), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      real (kind=kind_phys), dimension(0:nlay,nbands) :: clrurad,       &
     &       clrdrad, toturad, totdrad

      real (kind=kind_phys), dimension(nlay)   :: gassrcu, totsrcu,     &
     &       trngas, efclrfr, rfdelp
      real (kind=kind_phys), dimension(0:nlay) :: fnet, fnetc

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, radtotd, radclrd, rad0,          &
     &       clfr, trng, gasu

      integer :: ittot, itgas, ib, ig, k
!
!===> ...  begin here
!
      do ib = 1, NBANDS
        do k = 0, NLAY
          toturad(k,ib) = f_zero
          totdrad(k,ib) = f_zero
          clrurad(k,ib) = f_zero
          clrdrad(k,ib) = f_zero
        enddo
      enddo

      do k = 0, nlay
        totuflux(k) = f_zero
        totdflux(k) = f_zero
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

!  --- ...  loop over all g-points

      do ig = 1, ngptlw
        ib = ngb(ig)

        radtotd = f_zero
        radclrd = f_zero

!> -# Downward radiative transfer loop.

        do k = nlay, 1, -1

!!\n  - clear sky, gases contribution

          odepth = max( f_zero, secdif(ib)*tautot(ig,k) )
          if (odepth <= 0.06) then
            atrgas = odepth - 0.5*odepth*odepth
            trng   = f_one - atrgas
            gasfac = rec_6 * odepth
          else
            tblind = odepth / (bpade + odepth)
            itgas = tblint*tblind + 0.5
            trng  = exp_tbl(itgas)
            atrgas = f_one - trng
            gasfac = tfn_tbl(itgas)
            odepth = tau_tbl(itgas)
          endif

          plfrac = fracs(ig,k)
          blay = pklay(ib,k)

          dplnku = pklev(ib,k  ) - blay
          dplnkd = pklev(ib,k-1) - blay
          bbdgas = plfrac * (blay + dplnkd*gasfac)
          bbugas = plfrac * (blay + dplnku*gasfac)
          gassrcd= bbdgas * atrgas
          gassrcu(k)= bbugas * atrgas
          trngas(k) = trng

!!\n  - total sky, gases+clouds contribution

          clfr = cldfrc(k)
          if (clfr >= eps) then
!!\n  - cloudy layer

            odcld = secdif(ib) * taucld(ib,k)
            efclrfr(k) = f_one-(f_one - exp(-odcld))*clfr
            odtot = odepth + odcld
            if (odtot < 0.06) then
              totfac = rec_6 * odtot
              atrtot = odtot - 0.5*odtot*odtot
            else
              tblind = odtot / (bpade + odtot)
              ittot  = tblint*tblind + 0.5
              totfac = tfn_tbl(ittot)
              atrtot = f_one - exp_tbl(ittot)
            endif

            bbdtot = plfrac * (blay + dplnkd*totfac)
            bbutot = plfrac * (blay + dplnku*totfac)
            totsrcd= bbdtot * atrtot
            totsrcu(k)= bbutot * atrtot

!  --- ...  total sky radiance
            radtotd = radtotd*trng*efclrfr(k) + gassrcd                 &
     &              + clfr*(totsrcd - gassrcd)
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotd = radtotd*trng + gassrcd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

!> -# Compute spectral emissivity & reflectance, include the
!!    contribution of spectrally varying longwave emissivity and
!!     reflection from the surface to the upward radiative transfer.

!     note: spectral and Lambertian reflection are identical for the
!           diffusivity angle flux integration used here.

        reflct = f_one - semiss(ib)
        rad0 = semiss(ib) * fracs(ig,1) * pklay(ib,0)

!> -# Compute total sky radiance.
        radtotu = rad0 + reflct*radtotd
        toturad(0,ib) = toturad(0,ib) + radtotu

!> -# Compute clear sky radiance
        radclru = rad0 + reflct*radclrd
        clrurad(0,ib) = clrurad(0,ib) + radclru

!> -# Upward radiative transfer loop.

        do k = 1, nlay
          clfr = cldfrc(k)
          trng = trngas(k)
          gasu = gassrcu(k)

          if (clfr >= eps) then
!  --- ...  cloudy layer

!  --- ... total sky radiance
            radtotu = radtotu*trng*efclrfr(k) + gasu                    &
     &            + clfr*(totsrcu(k) - gasu)
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          else
!  --- ...  clear layer

!  --- ... total sky radiance
            radtotu = radtotu*trng + gasu
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

      enddo   ! end do_ig_loop

!> -# Process longwave output from band for total and clear streams.
!!    Calculate upward, downward, and net flux.

      flxfac = wtdiff * fluxfac

      do k = 0, nlay
        do ib = 1, nbands
          totuflux(k) = totuflux(k) + toturad(k,ib)
          totdflux(k) = totdflux(k) + totdrad(k,ib)
          totuclfl(k) = totuclfl(k) + clrurad(k,ib)
          totdclfl(k) = totdclfl(k) + clrdrad(k,ib)
        enddo

        totuflux(k) = totuflux(k) * flxfac
        totdflux(k) = totdflux(k) * flxfac
        totuclfl(k) = totuclfl(k) * flxfac
        totdclfl(k) = totdclfl(k) * flxfac
      enddo

!  --- ...  calculate net fluxes and heating rates
      fnet(0) = totuflux(0) - totdflux(0)

      do k = 1, nlay
        rfdelp(k) = heatfac / delp(k)
        fnet(k) = totuflux(k) - totdflux(k)
        htr (k) = (fnet(k-1) - fnet(k)) * rfdelp(k)
      enddo

!! --- ...  optional clear sky heating rates
      if ( lhlw0 ) then
        fnetc(0) = totuclfl(0) - totdclfl(0)

        do k = 1, nlay
          fnetc(k) = totuclfl(k) - totdclfl(k)
          htrcl(k) = (fnetc(k-1) - fnetc(k)) * rfdelp(k)
        enddo
      endif

!! --- ...  optional spectral band heating rates
      if ( lhlwb ) then
        do ib = 1, nbands
          fnet(0) = (toturad(0,ib) - totdrad(0,ib)) * flxfac

          do k = 1, nlay
            fnet(k) = (toturad(k,ib) - totdrad(k,ib)) * flxfac
            htrb(k,ib) = (fnet(k-1) - fnet(k)) * rfdelp(k)
          enddo
        enddo
      endif

! ..................................
      end subroutine rtrn
!! @}
! ----------------------------------


!>\ingroup module_radlw_main
!> This subroutine computes the upward/downward radiative fluxes, and
!! heating rates for both clear or cloudy atmosphere. Clouds are
!! assumed as in maximum-randomly overlaping in a vertical column.
!!\param semiss        lw surface emissivity
!!\param delp          layer pressure thickness (mb)
!!\param cldfrc        layer cloud fraction
!!\param taucld        layer cloud opt depth
!!\param tautot        total optical depth (gas+aerosols)
!!\param pklay         integrated planck func at lay temp
!!\param pklev         integrated planck func at lev temp
!!\param fracs         planck fractions
!!\param secdif        secant of diffusivity angle
!!\param nlay          number of vertical layers
!!\param nlp1          number of vertical levels (interfaces)
!!\param totuflux      total sky upward flux (\f$w/m^2\f$)
!!\param totdflux      total sky downward flux (\f$w/m^2\f$)
!!\param htr           total sky heating rate (k/sec or k/day)
!!\param totuclfl      clear sky upward flux (\f$w/m^2\f$)
!!\param totdclfl      clear sky downward flux (\f$w/m^2\f$)
!!\param htrcl         clear sky heating rate (k/sec or k/day)
!!\param htrb          spectral band lw heating rate (k/day)
!!\section gen_rtrnmr rtrnmr General Algorithm
!> @{
! ----------------------------------
      subroutine rtrnmr                                                 &
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &!  ---  inputs
     &       fracs,secdif, nlay,nlp1,                                   &
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &  !  ---  outputs:
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are assumed as in  !
! maximum-randomly overlaping in a vertical colum.                      !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfrc  - real, layer cloud fraction                         0:nlp1 !
!   taucld  - real, layer cloud opt depth                    nbands,nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw,nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw,nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck fn         1  !
!    totfac - real, gas+cld pade factor, used for planck fn          1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas only            nlay!
!    totsrcu- real, upwd source radiance due to gas + cld           nlay!
!    gassrcd- real, dnwd source radiance due to gas only             1  !
!    totsrcd- real, dnwd source radiance due to gas + cld            1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all Planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 W/m2 within cloudy layers.                !
!  clouds are treated with a maximum-random cloud overlap method.       !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(0:nlp1), intent(in) :: cldfrc
      real (kind=kind_phys), dimension(nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(nbands,nlay),intent(in):: taucld
      real (kind=kind_phys), dimension(ngptlw,nlay),intent(in):: fracs, &
     &       tautot

      real (kind=kind_phys), dimension(nbands,0:nlay), intent(in) ::    &
     &       pklev, pklay

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(nlay,nbands),intent(out) :: htrb

      real (kind=kind_phys), dimension(0:nlay), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      real (kind=kind_phys), dimension(0:nlay,nbands) :: clrurad,       &
     &       clrdrad, toturad, totdrad

      real (kind=kind_phys), dimension(nlay)   :: gassrcu, totsrcu,     &
     &       trngas, trntot, rfdelp
      real (kind=kind_phys), dimension(0:nlay) :: fnet, fnetc

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, radtotd, radclrd, rad0, rad,     &
     &       totradd, clrradd, totradu, clrradu, fmax, fmin, rat1, rat2,&
     &       radmod, clfr, trng, trnt, gasu, totu

      integer :: ittot, itgas, ib, ig, k

!  dimensions for cloud overlap adjustment
      real (kind=kind_phys), dimension(nlp1) :: faccld1u, faccld2u,     &
     &        facclr1u, facclr2u, faccmb1u, faccmb2u
      real (kind=kind_phys), dimension(0:nlay) :: faccld1d, faccld2d,   &
     &        facclr1d, facclr2d, faccmb1d, faccmb2d

      logical :: lstcldu(nlay), lstcldd(nlay)
!
!===> ...  begin here
!
      do k = 1, nlp1
        faccld1u(k) = f_zero
        faccld2u(k) = f_zero
        facclr1u(k) = f_zero
        facclr2u(k) = f_zero
        faccmb1u(k) = f_zero
        faccmb2u(k) = f_zero
      enddo

      lstcldu(1) = cldfrc(1) > eps
      rat1 = f_zero
      rat2 = f_zero

      do k = 1, nlay-1

        lstcldu(k+1) = cldfrc(k+1)>eps .and. cldfrc(k)<=eps

        if (cldfrc(k) > eps) then

!> -# Setup maximum/random cloud overlap.

          if (cldfrc(k+1) >= cldfrc(k)) then
            if (lstcldu(k)) then
              if (cldfrc(k) < f_one) then
                facclr2u(k+1) = (cldfrc(k+1) - cldfrc(k))               &
     &                        / (f_one - cldfrc(k))
              endif
              facclr2u(k) = f_zero
              faccld2u(k) = f_zero
            else
              fmax = max(cldfrc(k), cldfrc(k-1))
              if (cldfrc(k+1) > fmax) then
                facclr1u(k+1) = rat2
                facclr2u(k+1) = (cldfrc(k+1) - fmax)/(f_one - fmax)
              elseif (cldfrc(k+1) < fmax) then
                facclr1u(k+1) = (cldfrc(k+1) - cldfrc(k))               &
     &                        / (cldfrc(k-1) - cldfrc(k))
              else
                facclr1u(k+1) = rat2
              endif
            endif

            if (facclr1u(k+1)>f_zero .or. facclr2u(k+1)>f_zero) then
              rat1 = f_one
              rat2 = f_zero
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          else
            if (lstcldu(k)) then
              faccld2u(k+1) = (cldfrc(k) - cldfrc(k+1)) / cldfrc(k)
              facclr2u(k) = f_zero
              faccld2u(k) = f_zero
            else
              fmin = min(cldfrc(k), cldfrc(k-1))
              if (cldfrc(k+1) <= fmin) then
                faccld1u(k+1) = rat1
                faccld2u(k+1) = (fmin - cldfrc(k+1)) / fmin
              else
                faccld1u(k+1) = (cldfrc(k) - cldfrc(k+1))               &
     &                        / (cldfrc(k) - fmin)
              endif
            endif

            if (faccld1u(k+1)>f_zero .or. faccld2u(k+1)>f_zero) then
              rat1 = f_zero
              rat2 = f_one
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          endif

          faccmb1u(k+1) = facclr1u(k+1) * faccld2u(k) * cldfrc(k-1)
          faccmb2u(k+1) = faccld1u(k+1) * facclr2u(k)                   &
     &                  * (f_one - cldfrc(k-1))
        endif

      enddo

      do k = 0, nlay
        faccld1d(k) = f_zero
        faccld2d(k) = f_zero
        facclr1d(k) = f_zero
        facclr2d(k) = f_zero
        faccmb1d(k) = f_zero
        faccmb2d(k) = f_zero
      enddo

      lstcldd(nlay) = cldfrc(nlay) > eps
      rat1 = f_zero
      rat2 = f_zero

      do k = nlay, 2, -1

        lstcldd(k-1) = cldfrc(k-1) > eps .and. cldfrc(k)<=eps

        if (cldfrc(k) > eps) then

          if (cldfrc(k-1) >= cldfrc(k)) then
            if (lstcldd(k)) then
              if (cldfrc(k) < f_one) then
                facclr2d(k-1) = (cldfrc(k-1) - cldfrc(k))               &
     &                        / (f_one - cldfrc(k))
              endif

              facclr2d(k) = f_zero
              faccld2d(k) = f_zero
            else
              fmax = max(cldfrc(k), cldfrc(k+1))

              if (cldfrc(k-1) > fmax) then
                facclr1d(k-1) = rat2
                facclr2d(k-1) = (cldfrc(k-1) - fmax) / (f_one - fmax)
              elseif (cldfrc(k-1) < fmax) then
                facclr1d(k-1) = (cldfrc(k-1) - cldfrc(k))               &
     &                        / (cldfrc(k+1) - cldfrc(k))
              else
                facclr1d(k-1) = rat2
              endif
            endif

            if (facclr1d(k-1)>f_zero .or. facclr2d(k-1)>f_zero) then
              rat1 = f_one
              rat2 = f_zero
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          else
            if (lstcldd(k)) then
              faccld2d(k-1) = (cldfrc(k) - cldfrc(k-1)) / cldfrc(k)
              facclr2d(k) = f_zero
              faccld2d(k) = f_zero
            else
              fmin = min(cldfrc(k), cldfrc(k+1))

              if (cldfrc(k-1) <= fmin) then
                faccld1d(k-1) = rat1
                faccld2d(k-1) = (fmin - cldfrc(k-1)) / fmin
              else
                faccld1d(k-1) = (cldfrc(k) - cldfrc(k-1))               &
     &                        / (cldfrc(k) - fmin)
              endif
            endif

            if (faccld1d(k-1)>f_zero .or. faccld2d(k-1)>f_zero) then
              rat1 = f_zero
              rat2 = f_one
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          endif

          faccmb1d(k-1) = facclr1d(k-1) * faccld2d(k) * cldfrc(k+1)
          faccmb2d(k-1) = faccld1d(k-1) * facclr2d(k)                   &
     &                  * (f_one - cldfrc(k+1))
        endif

      enddo

!> -# Initialize for radiative transfer

      do ib = 1, NBANDS
        do k = 0, NLAY
          toturad(k,ib) = f_zero
          totdrad(k,ib) = f_zero
          clrurad(k,ib) = f_zero
          clrdrad(k,ib) = f_zero
        enddo
      enddo

      do k = 0, nlay
        totuflux(k) = f_zero
        totdflux(k) = f_zero
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

!  --- ...  loop over all g-points

      do ig = 1, ngptlw
        ib = ngb(ig)

        radtotd = f_zero
        radclrd = f_zero

!> -# Downward radiative transfer loop:

        do k = nlay, 1, -1

!  --- ...  clear sky, gases contribution

          odepth = max( f_zero, secdif(ib)*tautot(ig,k) )
          if (odepth <= 0.06) then
            atrgas = odepth - 0.5*odepth*odepth
            trng   = f_one - atrgas
            gasfac = rec_6 * odepth
          else
            tblind = odepth / (bpade + odepth)
            itgas = tblint*tblind + 0.5
            trng  = exp_tbl(itgas)
            atrgas = f_one - trng
            gasfac = tfn_tbl(itgas)
            odepth = tau_tbl(itgas)
          endif

          plfrac = fracs(ig,k)
          blay = pklay(ib,k)

          dplnku = pklev(ib,k  ) - blay
          dplnkd = pklev(ib,k-1) - blay
          bbdgas = plfrac * (blay + dplnkd*gasfac)
          bbugas = plfrac * (blay + dplnku*gasfac)
          gassrcd   = bbdgas * atrgas
          gassrcu(k)= bbugas * atrgas
          trngas(k) = trng

!  --- ...  total sky, gases+clouds contribution

          clfr = cldfrc(k)
          if (lstcldd(k)) then
            totradd = clfr * radtotd
            clrradd = radtotd - totradd
            rad = f_zero
          endif

          if (clfr >= eps) then
!>  - cloudy layer

            odcld = secdif(ib) * taucld(ib,k)
            odtot = odepth + odcld
            if (odtot < 0.06) then
              totfac = rec_6 * odtot
              atrtot = odtot - 0.5*odtot*odtot
              trnt   = f_one - atrtot
            else
              tblind = odtot / (bpade + odtot)
              ittot  = tblint*tblind + 0.5
              totfac = tfn_tbl(ittot)
              trnt   = exp_tbl(ittot)
              atrtot = f_one - trnt
            endif

            bbdtot = plfrac * (blay + dplnkd*totfac)
            bbutot = plfrac * (blay + dplnku*totfac)
            totsrcd   = bbdtot * atrtot
            totsrcu(k)= bbutot * atrtot
            trntot(k) = trnt

            totradd = totradd*trnt + clfr*totsrcd
            clrradd = clrradd*trng + (f_one - clfr)*gassrcd

!>  - total sky radiance
            radtotd = totradd + clrradd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!>  - clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

            radmod = rad*(facclr1d(k-1)*trng + faccld1d(k-1)*trnt)      &
     &             - faccmb1d(k-1)*gassrcd + faccmb2d(k-1)*totsrcd

            rad = -radmod + facclr2d(k-1)*(clrradd + radmod)            &
     &                    - faccld2d(k-1)*(totradd - radmod)
            totradd = totradd + rad
            clrradd = clrradd - rad

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotd = radtotd*trng + gassrcd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

!> -# Compute spectral emissivity & reflectance, include the
!!    contribution of spectrally varying longwave emissivity and
!!    reflection from the surface to the upward radiative transfer.

!     note: spectral and Lambertian reflection are identical for the
!           diffusivity angle flux integration used here.

        reflct = f_one - semiss(ib)
        rad0 = semiss(ib) * fracs(ig,1) * pklay(ib,0)

!> -# Compute total sky radiance.
        radtotu = rad0 + reflct*radtotd
        toturad(0,ib) = toturad(0,ib) + radtotu

!> -# Compute clear sky radiance.
        radclru = rad0 + reflct*radclrd
        clrurad(0,ib) = clrurad(0,ib) + radclru

!> -# Upward radiative transfer loop:

        do k = 1, nlay

          clfr = cldfrc(k)
          trng = trngas(k)
          gasu = gassrcu(k)

          if (lstcldu(k)) then
            totradu = clfr * radtotu
            clrradu = radtotu - totradu
            rad = f_zero
          endif

          if (clfr >= eps) then
!>  - cloudy layer radiance

            trnt = trntot(k)
            totu = totsrcu(k)
            totradu = totradu*trnt + clfr*totu
            clrradu = clrradu*trng + (f_one - clfr)*gasu

!>  - total sky radiance
            radtotu = totradu + clrradu
            toturad(k,ib) = toturad(k,ib) + radtotu

!>  - clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

            radmod = rad*(facclr1u(k+1)*trng + faccld1u(k+1)*trnt)      &
     &             - faccmb1u(k+1)*gasu + faccmb2u(k+1)*totu
            rad = -radmod + facclr2u(k+1)*(clrradu + radmod)            &
     &                    - faccld2u(k+1)*(totradu - radmod)
            totradu = totradu + rad
            clrradu = clrradu - rad

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotu = radtotu*trng + gasu
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ...  clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

      enddo   ! end do_ig_loop

!> -# Process longwave output from band for total and clear streams.
!! calculate upward, downward, and net flux.

      flxfac = wtdiff * fluxfac

      do k = 0, nlay
        do ib = 1, nbands
          totuflux(k) = totuflux(k) + toturad(k,ib)
          totdflux(k) = totdflux(k) + totdrad(k,ib)
          totuclfl(k) = totuclfl(k) + clrurad(k,ib)
          totdclfl(k) = totdclfl(k) + clrdrad(k,ib)
        enddo

        totuflux(k) = totuflux(k) * flxfac
        totdflux(k) = totdflux(k) * flxfac
        totuclfl(k) = totuclfl(k) * flxfac
        totdclfl(k) = totdclfl(k) * flxfac
      enddo

!  --- ...  calculate net fluxes and heating rates
      fnet(0) = totuflux(0) - totdflux(0)

      do k = 1, nlay
        rfdelp(k) = heatfac / delp(k)
        fnet(k) = totuflux(k) - totdflux(k)
        htr (k) = (fnet(k-1) - fnet(k)) * rfdelp(k)
      enddo

!! --- ...  optional clear sky heating rates
      if ( lhlw0 ) then
        fnetc(0) = totuclfl(0) - totdclfl(0)

        do k = 1, nlay
          fnetc(k) = totuclfl(k) - totdclfl(k)
          htrcl(k) = (fnetc(k-1) - fnetc(k)) * rfdelp(k)
        enddo
      endif

!! --- ...  optional spectral band heating rates
      if ( lhlwb ) then
        do ib = 1, nbands
          fnet(0) = (toturad(0,ib) - totdrad(0,ib)) * flxfac

          do k = 1, nlay
            fnet(k) = (toturad(k,ib) - totdrad(k,ib)) * flxfac
            htrb(k,ib) = (fnet(k-1) - fnet(k)) * rfdelp(k)
          enddo
        enddo
      endif

! .................................
      end subroutine rtrnmr
! ---------------------------------
!> @}

!>\ingroup module_radlw_main
!> \brief This subroutine computes the upward/downward radiative fluxes, and
!! heating rates for both clear or cloudy atmosphere.Clouds are treated
!! with the mcica stochastic approach.
!!
!!\param semiss       lw surface emissivity
!!\param delp         layer pressure thickness (mb)
!!\param cldfmc       layer cloud fraction (sub-column)
!!\param taucld       layer cloud opt depth
!!\param tautot       total optical depth (gas+aerosols)
!!\param pklay        integrated planck func at lay temp
!!\param pklev        integrated planck func at lev temp
!!\param fracs        planck fractions
!!\param secdif       secant of diffusivity angle
!!\param nlay         number of vertical layers
!!\param nlp1         number of vertical levels (interfaces)
!!\param totuflux     total sky upward flux \f$(w/m^2)\f$
!!\param totdflux     total sky downward flux \f$(w/m^2)\f$
!!\param htr          total sky heating rate (k/sec or k/day)
!!\param totuclfl     clear sky upward flux \f$(w/m^2)\f$
!!\param totdclfl     clear sky downward flux \f$(w/m^2)\f$
!!\param htrcl        clear sky heating rate (k/sec or k/day)
!!\param htrb         spectral band lw heating rate (k/day)
!!\section gen_rtrnmc rtrnmc General Algorithm
!> @{
! ---------------------------------
      subroutine rtrnmc                                                 &
     &     ( semiss,delp,cldfmc,taucld,tautot,pklay,pklev,              & !  ---  inputs:
     &       fracs,secdif, nlay,nlp1,                                   &
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       & !  ---  outputs:
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are treated with   !
! the mcica stochastic approach.                                        !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfmc  - real, layer cloud fraction (sub-column)        ngptlw*nlay!
!   taucld  - real, layer cloud opt depth                    nbands*nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw*nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw*nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    efclrfr- real, effective clear sky fraction (1-efcldfr)        nlay!
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck function   1  !
!    totfac - real, gas and cloud pade factor, used for planck fn    1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas                 nlay!
!    totsrcu- real, upwd source radiance due to gas+cld             nlay!
!    gassrcd- real, dnwd source radiance due to gas                  1  !
!    totsrcd- real, dnwd source radiance due to gas+cld              1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all Planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 W/m2 within cloudy layers.                !
!  clouds are treated with the mcica stochastic approach and            !
!  maximum-random cloud overlap.                                        !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(nbands,nlay),intent(in):: taucld
      real (kind=kind_phys), dimension(ngptlw,nlay),intent(in):: fracs, &
     &       tautot, cldfmc

      real (kind=kind_phys), dimension(nbands,0:nlay), intent(in) ::    &
     &       pklev, pklay

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(nlay,nbands),intent(out) :: htrb

      real (kind=kind_phys), dimension(0:nlay), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      real (kind=kind_phys), dimension(0:nlay,nbands) :: clrurad,       &
     &       clrdrad, toturad, totdrad

      real (kind=kind_phys), dimension(nlay)   :: gassrcu, totsrcu,     &
     &       trngas, efclrfr, rfdelp
      real (kind=kind_phys), dimension(0:nlay) :: fnet, fnetc

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, radtotd, radclrd, rad0,          &
     &       clfm, trng, gasu

      integer :: ittot, itgas, ib, ig, k
!
!===> ...  begin here
!
      do ib = 1, NBANDS
        do k = 0, NLAY
          toturad(k,ib) = f_zero
          totdrad(k,ib) = f_zero
          clrurad(k,ib) = f_zero
          clrdrad(k,ib) = f_zero
        enddo
      enddo

      do k = 0, nlay
        totuflux(k) = f_zero
        totdflux(k) = f_zero
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

!  --- ...  loop over all g-points

      do ig = 1, ngptlw
        ib = ngb(ig)

        radtotd = f_zero
        radclrd = f_zero

!> -# Downward radiative transfer loop.
!!\n  - Clear sky, gases contribution
!!\n  - Total sky, gases+clouds contribution
!!\n  - Cloudy layer
!!\n  - Total sky radiance
!!\n  - Clear sky radiance

        do k = nlay, 1, -1

!  --- ...  clear sky, gases contribution

          odepth = max( f_zero, secdif(ib)*tautot(ig,k) )
          if (odepth <= 0.06) then
            atrgas = odepth - 0.5*odepth*odepth
            trng   = f_one - atrgas
            gasfac = rec_6 * odepth
          else
            tblind = odepth / (bpade + odepth)
            itgas = tblint*tblind + 0.5
            trng  = exp_tbl(itgas)
            atrgas = f_one - trng
            gasfac = tfn_tbl(itgas)
            odepth = tau_tbl(itgas)
          endif

          plfrac = fracs(ig,k)
          blay = pklay(ib,k)

          dplnku = pklev(ib,k  ) - blay
          dplnkd = pklev(ib,k-1) - blay
          bbdgas = plfrac * (blay + dplnkd*gasfac)
          bbugas = plfrac * (blay + dplnku*gasfac)
          gassrcd= bbdgas * atrgas
          gassrcu(k)= bbugas * atrgas
          trngas(k) = trng

!  --- ...  total sky, gases+clouds contribution

          clfm = cldfmc(ig,k)
          if (clfm >= eps) then
!  --- ...  cloudy layer

            odcld = secdif(ib) * taucld(ib,k)
            efclrfr(k) = f_one - (f_one - exp(-odcld))*clfm
            odtot = odepth + odcld
            if (odtot < 0.06) then
              totfac = rec_6 * odtot
              atrtot = odtot - 0.5*odtot*odtot
            else
              tblind = odtot / (bpade + odtot)
              ittot  = tblint*tblind + 0.5
              totfac = tfn_tbl(ittot)
              atrtot = f_one - exp_tbl(ittot)
            endif

            bbdtot = plfrac * (blay + dplnkd*totfac)
            bbutot = plfrac * (blay + dplnku*totfac)
            totsrcd= bbdtot * atrtot
            totsrcu(k)= bbutot * atrtot

!  --- ...  total sky radiance
            radtotd = radtotd*trng*efclrfr(k) + gassrcd                 &
     &              + clfm*(totsrcd - gassrcd)
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotd = radtotd*trng + gassrcd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          endif   ! end if_clfm_block

        enddo   ! end do_k_loop

!> -# Compute spectral emissivity & reflectance, include the
!!    contribution of spectrally varying longwave emissivity and
!!    reflection from the surface to the upward radiative transfer.

!     note: spectral and Lambertian reflection are identical for the
!           diffusivity angle flux integration used here.

        reflct = f_one - semiss(ib)
        rad0 = semiss(ib) * fracs(ig,1) * pklay(ib,0)

!> -# Compute total sky radiance.
        radtotu = rad0 + reflct*radtotd
        toturad(0,ib) = toturad(0,ib) + radtotu

!> -# Compute clear sky radiance.
        radclru = rad0 + reflct*radclrd
        clrurad(0,ib) = clrurad(0,ib) + radclru

!> -# Upward radiative transfer loop.
!!\n  - Compute total sky radiance
!!\n  - Compute clear sky radiance

!          toturad holds summed radiance for total sky stream
!          clrurad holds summed radiance for clear sky stream

        do k = 1, nlay
          clfm = cldfmc(ig,k)
          trng = trngas(k)
          gasu = gassrcu(k)

          if (clfm > eps) then
!  --- ...  cloudy layer

!  --- ... total sky radiance
            radtotu = radtotu*trng*efclrfr(k) + gasu                    &
     &              + clfm*(totsrcu(k) - gasu)
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          else
!  --- ...  clear layer

!  --- ... total sky radiance
            radtotu = radtotu*trng + gasu
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          endif   ! end if_clfm_block

        enddo   ! end do_k_loop

      enddo   ! end do_ig_loop

!> -# Process longwave output from band for total and clear streams.
!!    Calculate upward, downward, and net flux.

      flxfac = wtdiff * fluxfac

      do k = 0, nlay
        do ib = 1, nbands
          totuflux(k) = totuflux(k) + toturad(k,ib)
          totdflux(k) = totdflux(k) + totdrad(k,ib)
          totuclfl(k) = totuclfl(k) + clrurad(k,ib)
          totdclfl(k) = totdclfl(k) + clrdrad(k,ib)
        enddo

        totuflux(k) = totuflux(k) * flxfac
        totdflux(k) = totdflux(k) * flxfac
        totuclfl(k) = totuclfl(k) * flxfac
        totdclfl(k) = totdclfl(k) * flxfac
      enddo

!> -# Calculate net fluxes and heating rates.
      fnet(0) = totuflux(0) - totdflux(0)

      do k = 1, nlay
        rfdelp(k) = heatfac / delp(k)
        fnet(k) = totuflux(k) - totdflux(k)
        htr (k) = (fnet(k-1) - fnet(k)) * rfdelp(k)
      enddo

!> -# Optional clear sky heating rates.
      if ( lhlw0 ) then
        fnetc(0) = totuclfl(0) - totdclfl(0)

        do k = 1, nlay
          fnetc(k) = totuclfl(k) - totdclfl(k)
          htrcl(k) = (fnetc(k-1) - fnetc(k)) * rfdelp(k)
        enddo
      endif

!> -# Optional spectral band heating rates.
      if ( lhlwb ) then
        do ib = 1, nbands
          fnet(0) = (toturad(0,ib) - totdrad(0,ib)) * flxfac

          do k = 1, nlay
            fnet(k) = (toturad(k,ib) - totdrad(k,ib)) * flxfac
            htrb(k,ib) = (fnet(k-1) - fnet(k)) * rfdelp(k)
          enddo
        enddo
      endif

! ..................................
      end subroutine rtrnmc
! ----------------------------------
!> @}

!>\ingroup module_radlw_main
!>\brief This subroutine contains optical depths developed for the rapid
!! radiative transfer model.
!!
!! It contains the subroutines \a taugbn (where n goes from
!! 1 to 16). \a taugbn calculates the optical depths and planck fractions
!! per g-value and layer for band n.
!!\param laytrop          tropopause layer index (unitless) layer at
!!                        which switch is made for key species
!!\param pavel            layer pressures (mb)
!!\param coldry           column amount for dry air \f$(mol/cm^2)\f$
!!\param colamt           column amounts of h2o, co2, o3, n2o, ch4,o2,
!!                        co \f$(mol/cm^2)\f$
!!\param colbrd           column amount of broadening gases
!!\param wx               cross-section amounts \f$(mol/cm^2)\f$
!!\param tauaer           aerosol optical depth
!!\param rfrate           reference ratios of binary species parameter
!!\n                      (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,4-h2o/ch4,
!!                                 5-n2o/co2,6-o3/co2
!!\n                      (:,:,n)n=1,2: the rates of ref press at the 2
!!                                 sides of the layer
!!\param fac00,fac01,fac10,fac11            factors multiply the reference ks, i,j of 0/1
!!                        for lower/higher of the 2 appropriate
!!                        temperatures and altitudes
!!\param jp               index of lower reference pressure
!!\param jt, jt1          indices of lower reference temperatures for
!!                        pressure levels jp and jp+1, respectively
!!\param selffac          scale factor for water vapor self-continuum
!!                        equals (water vapor density)/(atmospheric
!!                        density at 296k and 1013 mb)
!!\param selffrac         factor for temperature interpolation of
!!                        reference water vapor self-continuum data
!!\param indself          index of lower reference temperature for the
!!                        self-continuum interpolation
!!\param forfac           scale factor for w. v. foreign-continuum
!!\param forfrac          factor for temperature interpolation of
!!                        reference w.v. foreign-continuum data
!!\param indfor           index of lower reference temperature for the
!!                        foreign-continuum interpolation
!!\param minorfrac        factor for minor gases
!!\param scaleminor,scaleminorn2       scale factors for minor gases
!!\param indminor         index of lower reference temperature for
!!                        minor gases
!!\param nlay             total number of layers
!!\param fracs            planck fractions
!!\param tautot           total optical depth (gas+aerosols)
!>\section taumol_gen taumol General Algorithm
!! @{
!! subprograms called:  taugb## (## = 01 -16) 
      subroutine taumol                                                 &
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              & !  ---  inputs
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay,                                                      &
     &       fracs, tautot                                              & !  ---  outputs
     &     )

!  ************    original subprogram description    ***************   !
!                                                                       !
!                  optical depths developed for the                     !
!                                                                       !
!                rapid radiative transfer model (rrtm)                  !
!                                                                       !
!            atmospheric and environmental research, inc.               !
!                        131 hartwell avenue                            !
!                        lexington, ma 02421                            !
!                                                                       !
!                           eli j. mlawer                               !
!                         jennifer delamere                             !
!                         steven j. taubman                             !
!                         shepard a. clough                             !
!                                                                       !
!                       email:  mlawer@aer.com                          !
!                       email:  jdelamer@aer.com                        !
!                                                                       !
!        the authors wish to acknowledge the contributions of the       !
!        following people:  karen cady-pereira, patrick d. brown,       !
!        michael j. iacono, ronald e. farren, luke chen,                !
!        robert bergstrom.                                              !
!                                                                       !
!  revision for g-point reduction: michael j. iacono; aer, inc.         !
!                                                                       !
!     taumol                                                            !
!                                                                       !
!     this file contains the subroutines taugbn (where n goes from      !
!     1 to 16).  taugbn calculates the optical depths and planck        !
!     fractions per g-value and layer for band n.                       !
!                                                                       !
!  *******************************************************************  !
!  ==================   program usage description   ==================  !
!                                                                       !
!    call  taumol                                                       !
!       inputs:                                                         !
!          ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              !
!            rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  !
!            selffac,selffrac,indself,forfac,forfrac,indfor,            !
!            minorfrac,scaleminor,scaleminorn2,indminor,                !
!            nlay,                                                      !
!       outputs:                                                        !
!            fracs, tautot )                                            !
!                                                                       !
!  subprograms called:  taugb## (## = 01 -16)                           !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!     laytrop   - integer, tropopause layer index (unitless)        1   !
!                   layer at which switch is made for key species       !
!     pavel     - real, layer pressures (mb)                       nlay !
!     coldry    - real, column amount for dry air (mol/cm2)        nlay !
!     colamt    - real, column amounts of h2o, co2, o3, n2o, ch4,       !
!                   o2, co (mol/cm**2)                       nlay*maxgas!
!     colbrd    - real, column amount of broadening gases          nlay !
!     wx        - real, cross-section amounts(mol/cm2)      nlay*maxxsec!
!     tauaer    - real, aerosol optical depth               nbands*nlay !
!     rfrate    - real, reference ratios of binary species parameter    !
!     (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,4-h2o/ch4,5-n2o/co2,6-o3/co2!
!     (:,:,n)n=1,2: the rates of ref press at the 2 sides of the layer  !
!                                                          nlay*nrates*2!
!     facij     - real, factors multiply the reference ks, i,j of 0/1   !
!                   for lower/higher of the 2 appropriate temperatures  !
!                   and altitudes                                  nlay !
!     jp        - real, index of lower reference pressure          nlay !
!     jt, jt1   - real, indices of lower reference temperatures    nlay !
!                   for pressure levels jp and jp+1, respectively       !
!     selffac   - real, scale factor for water vapor self-continuum     !
!                   equals (water vapor density)/(atmospheric density   !
!                   at 296k and 1013 mb)                           nlay !
!     selffrac  - real, factor for temperature interpolation of         !
!                   reference water vapor self-continuum data      nlay !
!     indself   - integer, index of lower reference temperature for     !
!                   the self-continuum interpolation               nlay !
!     forfac    - real, scale factor for w. v. foreign-continuum   nlay !
!     forfrac   - real, factor for temperature interpolation of         !
!                   reference w.v. foreign-continuum data          nlay !
!     indfor    - integer, index of lower reference temperature for     !
!                   the foreign-continuum interpolation            nlay !
!     minorfrac - real, factor for minor gases                     nlay !
!     scaleminor,scaleminorn2                                           !
!               - real, scale factors for minor gases              nlay !
!     indminor  - integer, index of lower reference temperature for     !
!                   minor gases                                    nlay !
!     nlay      - integer, total number of layers                   1   !
!                                                                       !
!  outputs:                                                             !
!     fracs     - real, planck fractions                     ngptlw,nlay!
!     tautot    - real, total optical depth (gas+aerosols)   ngptlw,nlay!
!                                                                       !
!  internal variables:                                                  !
!     ng##      - integer, number of g-values in band ## (##=01-16) 1   !
!     nspa      - integer, for lower atmosphere, the number of ref      !
!                   atmos, each has different relative amounts of the   !
!                   key species for the band                      nbands!
!     nspb      - integer, same but for upper atmosphere          nbands!
!     absa      - real, k-values for lower ref atmospheres (no w.v.     !
!                   self-continuum) (cm**2/molecule)  nspa(##)*5*13*ng##!
!     absb      - real, k-values for high ref atmospheres (all sources) !
!                   (cm**2/molecule)               nspb(##)*5*13:59*ng##!
!     ka_m'mgas'- real, k-values for low ref atmospheres minor species  !
!                   (cm**2/molecule)                          mmn##*ng##!
!     kb_m'mgas'- real, k-values for high ref atmospheres minor species !
!                   (cm**2/molecule)                          mmn##*ng##!
!     selfref   - real, k-values for w.v. self-continuum for ref atmos  !
!                   used below laytrop (cm**2/mol)               10*ng##!
!     forref    - real, k-values for w.v. foreign-continuum for ref atmos
!                   used below/above laytrop (cm**2/mol)          4*ng##!
!                                                                       !
!  ******************************************************************   !

!  ---  inputs:
      integer, intent(in) :: nlay, laytrop

      integer, dimension(nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2

      real (kind=kind_phys), dimension(nlay,maxgas), intent(in):: colamt
      real (kind=kind_phys), dimension(nlay,maxxsec),intent(in):: wx

      real (kind=kind_phys), dimension(nbands,nlay), intent(in):: tauaer

      real (kind=kind_phys), dimension(nlay,nrates,2), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ngptlw,nlay), intent(out) ::     &
     &       fracs, tautot

!  ---  locals
      real (kind=kind_phys), dimension(ngptlw,nlay) :: taug

      integer :: ib, ig, k
!
!===> ...  begin here
!
      call taugb01
      call taugb02
      call taugb03
      call taugb04
      call taugb05
      call taugb06
      call taugb07
      call taugb08
      call taugb09
      call taugb10
      call taugb11
      call taugb12
      call taugb13
      call taugb14
      call taugb15
      call taugb16

!  ---  combine gaseous and aerosol optical depths

      do ig = 1, ngptlw
        ib = ngb(ig)

        do k = 1, nlay
          tautot(ig,k) = taug(ig,k) + tauaer(ib,k)
        enddo
      enddo

! =================
      contains
! =================

!>\ingroup module_radlw_main
!> band 1:  10-350 cm-1 (low key - h2o; low minor - n2);
!!  (high key - h2o; high minor - n2)
! ----------------------------------
      subroutine taugb01
! ..................................

!  ------------------------------------------------------------------  !
!  written by eli j. mlawer, atmospheric & environmental research.     !
!  revised by michael j. iacono, atmospheric & environmental research. !
!                                                                      !
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)             !
!                          (high key - h2o; high minor - n2)           !
!                                                                      !
!  compute the optical depth by interpolating in ln(pressure) and      !
!  temperature.  below laytrop, the water vapor self-continuum and     !
!  foreign continuum is interpolated (in temperature) separately.      !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb01

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &           indm, indmp, ig

      real (kind=kind_phys) :: pp, corradj, scalen2, tauself, taufor,   &
     &       taun2
!
!===> ...  begin here
!
!  ---  minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(1) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(1) + 1
        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1

        pp = pavel(k)
        scalen2 = colbrd(k) * scaleminorn2(k)
        if (pp < 250.0) then
          corradj = f_one - 0.15 * (250.0-pp) / 154.4
        else
          corradj = f_one
        endif

        do ig = 1, ng01
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) -  forref(ig,indf)))
          taun2   = scalen2 * (ka_mn2(ig,indm) + minorfrac(k)           &
     &            * (ka_mn2(ig,indmp) - ka_mn2(ig,indm)))

          taug(ig,k) = corradj * (colamt(k,1)                           &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor + taun2)

          fracs(ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(1) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(1) + 1
        indf = indfor(k)
        indm = indminor(k)

        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indfp = indf + 1
        indmp = indm + 1

        scalen2 = colbrd(k) * scaleminorn2(k)
        corradj = f_one - 0.15 * (pavel(k) / 95.6)

        do ig = 1, ng01
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))
          taun2  = scalen2 * (kb_mn2(ig,indm) + minorfrac(k)            &
     &           * (kb_mn2(ig,indmp) - kb_mn2(ig,indm)))

          taug(ig,k) = corradj * (colamt(k,1)                           &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))   &
     &           + taufor + taun2)

          fracs(ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb01
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
! ----------------------------------
      subroutine taugb02
! ..................................

!  ------------------------------------------------------------------  !
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)            !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb02

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &           ig

      real (kind=kind_phys) :: corradj, tauself, taufor
!
!===> ...  begin here
!
!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(2) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(2) + 1
        inds = indself(k)
        indf = indfor(k)

        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1

        corradj = f_one - 0.05 * (pavel(k) - 100.0) / 900.0

        do ig = 1, ng02
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns02+ig,k) = corradj * (colamt(k,1)                      &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor)

          fracs(ns02+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(2) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(2) + 1
        indf = indfor(k)

        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indfp = indf + 1

        do ig = 1, ng02
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns02+ig,k) = colamt(k,1)                                 &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))   &
     &           + taufor

          fracs(ns02+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb02
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o);
!!                        (high key - h2o,co2; high minor - n2o)
! ----------------------------------
      subroutine taugb03
! ..................................

!  ------------------------------------------------------------------  !
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)       !
!                           (high key - h2o,co2; high minor - n2o)     !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb03

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmn2o, jmn2op,   &
     &       id001, id011, id101, id111, id201, id211, jpl, jplp,       &
     &       ig, js, js1

      real (kind=kind_phys) ::  absn2o, ratn2o, adjfac, adjcoln2o,      &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mn2o,  specparm_mn2o,  specmult_mn2o,  fmn2o,      &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b,   &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      tau_major, tau_major1, tauself, taufor, n2om1, n2om2,       &
     &      p, p4, fk0, fk1, fk2
!
!===> ...  begin here
!
!  --- ...  minor gas mapping levels:
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

      refrat_planck_a = chi_mls(1,9)/chi_mls(2,9)    ! P = 212.725 mb
      refrat_planck_b = chi_mls(1,13)/chi_mls(2,13)  ! P = 95.58   mb
      refrat_m_a      = chi_mls(1,3)/chi_mls(2,3)    ! P = 706.270 mb
      refrat_m_b      = chi_mls(1,13)/chi_mls(2,13)  ! P = 95.58   mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(3) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(3) + js1

        speccomb_mn2o = colamt(k,1) + refrat_m_a*colamt(k,2)
        specparm_mn2o = colamt(k,1) / speccomb_mn2o
        specmult_mn2o = 8.0 * min(specparm_mn2o, oneminus)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = mod(specmult_mn2o, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jmn2op= jmn2o+ 1
        jplp  = jpl  + 1

!  --- ...  in atmospheres where the amount of n2O is too great to be considered
!           a minor species, adjust the column amount of n2O by an empirical factor
!           to obtain the proper contribution.

        p = coldry(k) * chi_mls(4,jp(k)+1)
        ratn2o = colamt(k,4) / p
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o - 0.5)**0.65
          adjcoln2o = adjfac * p
        else
          adjcoln2o = colamt(k,4)
        endif

        if (specparm < 0.125) then
          p = fs - f_one
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        else if (specparm > 0.875) then
          p = -fs
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk0 = f_one - fs
          fk1 = fs
          fk2 = f_zero
          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk0*fac00(k)
        fac100 = fk1*fac00(k)
        fac200 = fk2*fac00(k)
        fac010 = fk0*fac10(k)
        fac110 = fk1*fac10(k)
        fac210 = fk2*fac10(k)

        if (specparm1 < 0.125) then
          p = fs1 - f_one
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p = -fs1
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk0 = f_one - fs1
          fk1 = fs1
          fk2 = f_zero
          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk0*fac01(k)
        fac101 = fk1*fac01(k)
        fac201 = fk2*fac01(k)
        fac011 = fk0*fac11(k)
        fac111 = fk1*fac11(k)
        fac211 = fk2*fac11(k)

        do ig = 1, ng03
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          n2om1   = ka_mn2o(ig,jmn2o,indm) + fmn2o                      &
     &            * (ka_mn2o(ig,jmn2op,indm) - ka_mn2o(ig,jmn2o,indm))
          n2om2   = ka_mn2o(ig,jmn2o,indmp) + fmn2o                     &
     &            * (ka_mn2o(ig,jmn2op,indmp) - ka_mn2o(ig,jmn2o,indmp))
          absn2o  = n2om1 + minorfrac(k) * (n2om2 - n2om1)

          tau_major = speccomb                                          &
     &              * (fac000*absa(ig,id000) + fac010*absa(ig,id010)    &
     &              +  fac100*absa(ig,id100) + fac110*absa(ig,id110)    &
     &              +  fac200*absa(ig,id200) + fac210*absa(ig,id210))

          tau_major1 = speccomb1                                        &
     &              * (fac001*absa(ig,id001) + fac011*absa(ig,id011)    &
     &              +  fac101*absa(ig,id101) + fac111*absa(ig,id111)    &
     &              +  fac201*absa(ig,id201) + fac211*absa(ig,id211))

          taug(ns03+ig,k) = tau_major + tau_major1                      &
     &                    + tauself + taufor + adjcoln2o*absn2o

          fracs(ns03+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo     ! end do_k_loop
      enddo   ! end do_ig_loop

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 4.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-13)*5 + (jt(k)-1)) * nspb(3) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 4.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(3) + js1

        speccomb_mn2o = colamt(k,1) + refrat_m_b*colamt(k,2)
        specparm_mn2o = colamt(k,1) / speccomb_mn2o
        specmult_mn2o = 4.0 * min(specparm_mn2o, oneminus)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = mod(specmult_mn2o, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_b*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 4.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        indf = indfor(k)
        indm = indminor(k)
        indfp = indf + 1
        indmp = indm + 1
        jmn2op= jmn2o+ 1
        jplp  = jpl  + 1

        id000 = ind0
        id010 = ind0 + 5
        id100 = ind0 + 1
        id110 = ind0 + 6
        id001 = ind1
        id011 = ind1 + 5
        id101 = ind1 + 1
        id111 = ind1 + 6

!  --- ...  in atmospheres where the amount of n2o is too great to be considered
!           a minor species, adjust the column amount of N2O by an empirical factor
!           to obtain the proper contribution.

        p = coldry(k) * chi_mls(4,jp(k)+1)
        ratn2o = colamt(k,4) / p
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o - 0.5)**0.65
          adjcoln2o = adjfac * p
        else
          adjcoln2o = colamt(k,4)
        endif

        fk0 = f_one - fs
        fk1 = fs
        fac000 = fk0*fac00(k)
        fac010 = fk0*fac10(k)
        fac100 = fk1*fac00(k)
        fac110 = fk1*fac10(k)

        fk0 = f_one - fs1
        fk1 = fs1
        fac001 = fk0*fac01(k)
        fac011 = fk0*fac11(k)
        fac101 = fk1*fac01(k)
        fac111 = fk1*fac11(k)

        do ig = 1, ng03
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))
          n2om1  = kb_mn2o(ig,jmn2o,indm) + fmn2o                       &
     &           * (kb_mn2o(ig,jmn2op,indm) - kb_mn2o(ig,jmn2o,indm))
          n2om2  = kb_mn2o(ig,jmn2o,indmp) + fmn2o                      &
     &           * (kb_mn2o(ig,jmn2op,indmp) - kb_mn2o(ig,jmn2o,indmp))
          absn2o = n2om1 + minorfrac(k) * (n2om2 - n2om1)

          tau_major = speccomb                                          &
     &              * (fac000*absb(ig,id000) + fac010*absb(ig,id010)    &
     &              +  fac100*absb(ig,id100) + fac110*absb(ig,id110))

          tau_major1 = speccomb1                                        &
     &              * (fac001*absb(ig,id001) + fac011*absb(ig,id011)    &
     &              +  fac101*absb(ig,id101) + fac111*absb(ig,id111))

          taug(ns03+ig,k) = tau_major + tau_major1                      &
     &                    + taufor + adjcoln2o*absn2o

          fracs(ns03+ig,k) = fracrefb(ig,jpl) + fpl                     &
     &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
        enddo
      enddo

! ..................................
      end subroutine taugb03
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
! ----------------------------------
      subroutine taugb04
! ..................................

!  ------------------------------------------------------------------  !
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)     !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb04

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, jpl, jplp,    &
     &       id000, id010, id100, id110, id200, id210, ig, js, js1,     &
     &       id001, id011, id101, id111, id201, id211

      real (kind=kind_phys) :: tauself, taufor, p, p4, fk0, fk1, fk2,   &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      refrat_planck_a, refrat_planck_b, tau_major, tau_major1
!
!===> ...  begin here
!
      refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)     ! P = 142.5940 mb
      refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)     ! P = 95.58350 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(4) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = ( jp(k)*5 + (jt1(k)-1)) * nspa(4) + js1

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, 1.0)

        inds = indself(k)
        indf = indfor(k)
        indsp = inds + 1
        indfp = indf + 1
        jplp  = jpl  + 1

        if (specparm < 0.125) then
          p = fs - f_one
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p = -fs
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk0 = f_one - fs
          fk1 = fs
          fk2 = f_zero
          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk0*fac00(k)
        fac100 = fk1*fac00(k)
        fac200 = fk2*fac00(k)
        fac010 = fk0*fac10(k)
        fac110 = fk1*fac10(k)
        fac210 = fk2*fac10(k)

        if (specparm1 < 0.125) then
          p = fs1 - f_one
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p = -fs1
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk0 = f_one - fs1
          fk1 = fs1
          fk2 = f_zero
          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk0*fac01(k)
        fac101 = fk1*fac01(k)
        fac201 = fk2*fac01(k)
        fac011 = fk0*fac11(k)
        fac111 = fk1*fac11(k)
        fac211 = fk2*fac11(k)

        do ig = 1, ng04
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          tau_major = speccomb                                          &
     &              * (fac000*absa(ig,id000) + fac010*absa(ig,id010)    &
     &              +  fac100*absa(ig,id100) + fac110*absa(ig,id110)    &
     &              +  fac200*absa(ig,id200) + fac210*absa(ig,id210))

          tau_major1 = speccomb1                                        &
     &              * (fac001*absa(ig,id001) + fac011*absa(ig,id011)    &
     &              +  fac101*absa(ig,id101) + fac111*absa(ig,id111)    &
     &              +  fac201*absa(ig,id201) + fac211*absa(ig,id211))

          taug(ns04+ig,k) = tau_major + tau_major1 + tauself + taufor

          fracs(ns04+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo     ! end do_k_loop
      enddo   ! end do_ig_loop

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        speccomb = colamt(k,3) + rfrate(k,6,1)*colamt(k,2)
        specparm = colamt(k,3) / speccomb
        specmult = 4.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-13)*5 + (jt(k)-1)) * nspb(4) + js

        speccomb1 = colamt(k,3) + rfrate(k,6,2)*colamt(k,2)
        specparm1 = colamt(k,3) / speccomb1
        specmult1 = 4.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(4) + js1

        speccomb_planck = colamt(k,3) + refrat_planck_b*colamt(k,2)
        specparm_planck = colamt(k,3) / speccomb_planck
        specmult_planck = 4.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)
        jplp = jpl + 1

        id000 = ind0
        id010 = ind0 + 5
        id100 = ind0 + 1
        id110 = ind0 + 6
        id001 = ind1
        id011 = ind1 + 5
        id101 = ind1 + 1
        id111 = ind1 + 6

        fk0 = f_one - fs
        fk1 = fs
        fac000 = fk0*fac00(k)
        fac010 = fk0*fac10(k)
        fac100 = fk1*fac00(k)
        fac110 = fk1*fac10(k)

        fk0 = f_one - fs1
        fk1 = fs1
        fac001 = fk0*fac01(k)
        fac011 = fk0*fac11(k)
        fac101 = fk1*fac01(k)
        fac111 = fk1*fac11(k)

        do ig = 1, ng04
          tau_major =  speccomb                                         &
     &              * (fac000*absb(ig,id000) + fac010*absb(ig,id010)    &
     &              +  fac100*absb(ig,id100) + fac110*absb(ig,id110))
          tau_major1 = speccomb1                                        &
     &              * (fac001*absb(ig,id001) + fac011*absb(ig,id011)    &
     &              +  fac101*absb(ig,id101) + fac111*absb(ig,id111))

          taug(ns04+ig,k) =  tau_major + tau_major1

          fracs(ns04+ig,k) = fracrefb(ig,jpl) + fpl                     &
     &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
        enddo

!  --- ...  empirical modification to code to improve stratospheric cooling rates
!           for co2. revised to apply weighting for g-point reduction in this band.

        taug(ns04+ 8,k) = taug(ns04+ 8,k) * 0.92
        taug(ns04+ 9,k) = taug(ns04+ 9,k) * 0.88
        taug(ns04+10,k) = taug(ns04+10,k) * 1.07
        taug(ns04+11,k) = taug(ns04+11,k) * 1.1
        taug(ns04+12,k) = taug(ns04+12,k) * 0.99
        taug(ns04+13,k) = taug(ns04+13,k) * 0.88
        taug(ns04+14,k) = taug(ns04+14,k) * 0.943
      enddo

! ..................................
      end subroutine taugb04
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!!                       (high key - o3,co2)
! ----------------------------------
      subroutine taugb05
! ..................................

!  ------------------------------------------------------------------  !
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)  !
!                           (high key - o3,co2)                        !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb05

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmo3, jmo3p,     &
     &       id001, id011, id101, id111, id201, id211, jpl, jplp,       &
     &       ig, js, js1

      real (kind=kind_phys)  :: tauself, taufor, o3m1, o3m2, abso3,     &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mo3,   specparm_mo3,   specmult_mo3,   fmo3,       &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_planck_b, refrat_m_a,               &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower/upper atmosphere.

      refrat_planck_a = chi_mls(1,5)/chi_mls(2,5)      ! P = 473.420 mb
      refrat_planck_b = chi_mls(3,43)/chi_mls(2,43)    ! P = 0.2369  mb
      refrat_m_a = chi_mls(1,7)/chi_mls(2,7)           ! P = 317.348 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(5) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(5) + js1

        speccomb_mo3 = colamt(k,1) + refrat_m_a*colamt(k,2)
        specparm_mo3 = colamt(k,1) / speccomb_mo3
        specmult_mo3 = 8.0 * min(specparm_mo3, oneminus)
        jmo3 = 1 + int(specmult_mo3)
        fmo3 = mod(specmult_mo3, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmo3p = jmo3 + 1

        if (specparm < 0.125) then
          p0   = fs - f_one
          p40  = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p0   = -fs
          p40  = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        if (specparm1 < 0.125) then
          p1   = fs1 - f_one
          p41  = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p1   = -fs1
          p41  = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng05
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          o3m1    = ka_mo3(ig,jmo3,indm) + fmo3                         &
     &            * (ka_mo3(ig,jmo3p,indm) -  ka_mo3(ig,jmo3,indm))
          o3m2    = ka_mo3(ig,jmo3,indmp) + fmo3                        &
     &            * (ka_mo3(ig,jmo3p,indmp) - ka_mo3(ig,jmo3,indmp))
          abso3   = o3m1 + minorfrac(k)*(o3m2 - o3m1)

          taug(ns05+ig,k) = speccomb                                    &
     &            * (fac000*absa(ig,id000) + fac010*absa(ig,id010)      &
     &            +  fac100*absa(ig,id100) + fac110*absa(ig,id110)      &
     &            +  fac200*absa(ig,id200) + fac210*absa(ig,id210))     &
     &            +     speccomb1                                       &
     &            * (fac001*absa(ig,id001) + fac011*absa(ig,id011)      &
     &            +  fac101*absa(ig,id101) + fac111*absa(ig,id111)      &
     &            +  fac201*absa(ig,id201) + fac211*absa(ig,id211))     &
     &            + tauself + taufor+abso3*colamt(k,3)+wx(k,1)*ccl4(ig)

          fracs(ns05+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        speccomb = colamt(k,3) + rfrate(k,6,1)*colamt(k,2)
        specparm = colamt(k,3) / speccomb
        specmult = 4.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-13)*5 + (jt(k)-1)) * nspb(5) + js

        speccomb1 = colamt(k,3) + rfrate(k,6,2)*colamt(k,2)
        specparm1 = colamt(k,3) / speccomb1
        specmult1 = 4.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(5) + js1

        speccomb_planck = colamt(k,3) + refrat_planck_b*colamt(k,2)
        specparm_planck = colamt(k,3) / speccomb_planck
        specmult_planck = 4.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)
        jplp= jpl + 1

        id000 = ind0
        id010 = ind0 + 5
        id100 = ind0 + 1
        id110 = ind0 + 6
        id001 = ind1
        id011 = ind1 + 5
        id101 = ind1 + 1
        id111 = ind1 + 6

        fk00 = f_one - fs
        fk10 = fs

        fk01 = f_one - fs1
        fk11 = fs1

        fac000 = fk00 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac100 = fk10 * fac00(k)
        fac110 = fk10 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac101 = fk11 * fac01(k)
        fac111 = fk11 * fac11(k)

        do ig = 1, ng05
          taug(ns05+ig,k) = speccomb                                    &
     &                * (fac000*absb(ig,id000) + fac010*absb(ig,id010)  &
     &                +  fac100*absb(ig,id100) + fac110*absb(ig,id110)) &
     &                +     speccomb1                                   &
     &                * (fac001*absb(ig,id001) + fac011*absb(ig,id011)  &
     &                +  fac101*absb(ig,id101) + fac111*absb(ig,id111)) &
     &                + wx(k,1) * ccl4(ig)

          fracs(ns05+ig,k) = fracrefb(ig,jpl) + fpl                     &
     &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
        enddo
      enddo

! ..................................
      end subroutine taugb05
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!!                       (high key - none; high minor - cfc11, cfc12)
! ----------------------------------
      subroutine taugb06
! ..................................

!  ------------------------------------------------------------------  !
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)           !
!                           (high key - none; high minor - cfc11, cfc12)
!  ------------------------------------------------------------------  !

      use module_radlw_kgb06

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: ratco2, adjfac, adjcolco2, tauself,      &
     &      taufor, absco2, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(6) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(6) + 1

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        ind0p = ind0 + 1
        ind1p = ind1 + 1

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.77
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        do ig = 1, ng06
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          absco2  = ka_mco2(ig,indm) + minorfrac(k)                     &
     &            * (ka_mco2(ig,indmp) - ka_mco2(ig,indm))

          taug(ns06+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            +  tauself + taufor + adjcolco2*absco2                &
     &            +  wx(k,2)*cfc11adj(ig) + wx(k,3)*cfc12(ig)

          fracs(ns06+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop
!           nothing important goes on above laytrop in this band.

      do k = laytrop+1, nlay
        do ig = 1, ng06
          taug(ns06+ig,k) = wx(k,2)*cfc11adj(ig) + wx(k,3)*cfc12(ig)

          fracs(ns06+ig,k) = fracrefa(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb06
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!!                        (high key - o3; high minor - co2)
! ----------------------------------
      subroutine taugb07
! ..................................

!  ------------------------------------------------------------------  !
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)       !
!                            (high key - o3; high minor - co2)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb07

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, indm, indmp,     &
     &       id001, id011, id101, id111, id201, id211, jmco2, jmco2p,   &
     &       jpl, jplp, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, co2m1, co2m2, absco2,   &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mco2,  specparm_mco2,  specmult_mco2,  fmco2,      &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_m_a, ratco2, adjfac, adjcolco2,     &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - co2, p = 706.2620 mbar, t= 278.94 k
!     upper - co2, p = 12.9350 mbar, t = 234.01 k

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower atmosphere.

      refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)     ! P = 706.2620 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(3,3)          ! P = 706.2720 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,2,1)*colamt(k,3)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(7) + js

        speccomb1 = colamt(k,1) + rfrate(k,2,2)*colamt(k,3)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(7) + js1

        speccomb_mco2 = colamt(k,1) + refrat_m_a*colamt(k,3)
        specparm_mco2 = colamt(k,1) / speccomb_mco2
        specmult_mco2 = 8.0 * min(specparm_mco2, oneminus)
        jmco2 = 1 + int(specmult_mco2)
        fmco2 = mod(specmult_mco2, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,3)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmco2p= jmco2+ 1
        ind0p = ind0 + 1
        ind1p = ind1 + 1

!  --- ...  in atmospheres where the amount of CO2 is too great to be considered
!           a minor species, adjust the column amount of CO2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 3.0 + (ratco2-3.0)**0.79
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        if (specparm < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        if (specparm1 < 0.125) then
          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng07
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          co2m1   = ka_mco2(ig,jmco2,indm) + fmco2                      &
     &            * (ka_mco2(ig,jmco2p,indm) - ka_mco2(ig,jmco2,indm))
          co2m2   = ka_mco2(ig,jmco2,indmp) + fmco2                     &
     &            * (ka_mco2(ig,jmco2p,indmp) - ka_mco2(ig,jmco2,indmp))
          absco2  = co2m1 + minorfrac(k) * (co2m2 - co2m1)

          taug(ns07+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor + adjcolco2*absco2

          fracs(ns07+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

      do k = laytrop+1, nlay
        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.79
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(7) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(7) + 1

        indm = indminor(k)
        indmp = indm + 1
        ind0p = ind0 + 1
        ind1p = ind1 + 1

        do ig = 1, ng07
          absco2 = kb_mco2(ig,indm) + minorfrac(k)                      &
     &           * (kb_mco2(ig,indmp) - kb_mco2(ig,indm))

          taug(ns07+ig,k) = colamt(k,3)                                 &
     &            * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)   &
     &            +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))  &
     &            + adjcolco2 * absco2

          fracs(ns07+ig,k) = fracrefb(ig)
        enddo

!  --- ...  empirical modification to code to improve stratospheric cooling rates
!           for o3.  revised to apply weighting for g-point reduction in this band.

        taug(ns07+ 6,k) = taug(ns07+ 6,k) * 0.92
        taug(ns07+ 7,k) = taug(ns07+ 7,k) * 0.88
        taug(ns07+ 8,k) = taug(ns07+ 8,k) * 1.07
        taug(ns07+ 9,k) = taug(ns07+ 9,k) * 1.1
        taug(ns07+10,k) = taug(ns07+10,k) * 0.99
        taug(ns07+11,k) = taug(ns07+11,k) * 0.855
      enddo

! ..................................
      end subroutine taugb07
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!!                         (high key - o3; high minor - co2, n2o)
! ----------------------------------
      subroutine taugb08
! ..................................

!  ------------------------------------------------------------------  !
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)  !
!                             (high key - o3; high minor - co2, n2o)   !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb08

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: tauself, taufor, absco2, abso3, absn2o,  &
     &      ratco2, adjfac, adjcolco2, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(8) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(8) + 1

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.65
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        do ig = 1, ng08
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          absco2  = (ka_mco2(ig,indm) + minorfrac(k)                    &
     &            * (ka_mco2(ig,indmp) - ka_mco2(ig,indm)))
          abso3   = (ka_mo3(ig,indm) + minorfrac(k)                     &
     &            * (ka_mo3(ig,indmp) - ka_mo3(ig,indm)))
          absn2o  = (ka_mn2o(ig,indm) + minorfrac(k)                    &
     &            * (ka_mn2o(ig,indmp) - ka_mn2o(ig,indm)))

          taug(ns08+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself+taufor + adjcolco2*absco2                   &
     &            + colamt(k,3)*abso3 + colamt(k,4)*absn2o              &
     &            + wx(k,3)*cfc12(ig) + wx(k,4)*cfc22adj(ig)

          fracs(ns08+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(8) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(8) + 1

        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indmp = indm + 1

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.65
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        do ig = 1, ng08
          absco2 = (kb_mco2(ig,indm) + minorfrac(k)                     &
     &           * (kb_mco2(ig,indmp) - kb_mco2(ig,indm)))
          absn2o = (kb_mn2o(ig,indm) + minorfrac(k)                     &
     &           * (kb_mn2o(ig,indmp) - kb_mn2o(ig,indm)))

          taug(ns08+ig,k) = colamt(k,3)                                 &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))   &
     &           + adjcolco2*absco2 + colamt(k,4)*absn2o                &
     &           + wx(k,3)*cfc12(ig) + wx(k,4)*cfc22adj(ig)

          fracs(ns08+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb08
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!!                         (high key - ch4; high minor - n2o)
! ----------------------------------
      subroutine taugb09
! ..................................

!  ------------------------------------------------------------------  !
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)     !
!                             (high key - ch4; high minor - n2o)       !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb09

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, indm, indmp,     &
     &       id001, id011, id101, id111, id201, id211, jmn2o, jmn2op,   &
     &       jpl, jplp, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, n2om1, n2om2, absn2o,   &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mn2o,  specparm_mn2o,  specmult_mn2o,  fmn2o,     &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, ratn2o, adjfac, adjcoln2o,    &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower/upper atmosphere.

      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)       ! P = 212 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)            ! P = 706.272 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,4,1)*colamt(k,5)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(9) + js

        speccomb1 = colamt(k,1) + rfrate(k,4,2)*colamt(k,5)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(9) + js1

        speccomb_mn2o = colamt(k,1) + refrat_m_a*colamt(k,5)
        specparm_mn2o = colamt(k,1) / speccomb_mn2o
        specmult_mn2o = 8.0 * min(specparm_mn2o, oneminus)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = mod(specmult_mn2o, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,5)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmn2op= jmn2o+ 1

!  --- ...  in atmospheres where the amount of n2o is too great to be considered
!           a minor species, adjust the column amount of n2o by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(4,jp(k)+1)
        ratn2o = colamt(k,4) / temp
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o-0.5)**0.65
          adjcoln2o = adjfac * temp
        else
          adjcoln2o = colamt(k,4)
        endif

        if (specparm < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        if (specparm1 < 0.125) then
          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng09
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          n2om1   = ka_mn2o(ig,jmn2o,indm) + fmn2o                      &
     &            * (ka_mn2o(ig,jmn2op,indm) - ka_mn2o(ig,jmn2o,indm))
          n2om2   = ka_mn2o(ig,jmn2o,indmp) + fmn2o                     &
     &            * (ka_mn2o(ig,jmn2op,indmp) - ka_mn2o(ig,jmn2o,indmp))
          absn2o  = n2om1 + minorfrac(k) * (n2om2 - n2om1)

          taug(ns09+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor + adjcoln2o*absn2o

          fracs(ns09+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(9) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(9) + 1

        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indmp = indm + 1

!  --- ...  in atmospheres where the amount of n2o is too great to be considered
!           a minor species, adjust the column amount of n2o by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(4,jp(k)+1)
        ratn2o = colamt(k,4) / temp
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o - 0.5)**0.65
          adjcoln2o = adjfac * temp
        else
          adjcoln2o = colamt(k,4)
        endif

        do ig = 1, ng09
          absn2o = kb_mn2o(ig,indm) + minorfrac(k)                      &
     &           * (kb_mn2o(ig,indmp) - kb_mn2o(ig,indm))

          taug(ns09+ig,k) = colamt(k,5)                                 &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))   &
     &           + adjcoln2o*absn2o

          fracs(ns09+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb09
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
! ----------------------------------
      subroutine taugb10
! ..................................

!  ------------------------------------------------------------------  !
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb10

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       ig

      real (kind=kind_phys) :: tauself, taufor
!
!===> ...  begin here
!
!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(10) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(10) + 1

        inds = indself(k)
        indf = indfor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1

        do ig = 1, ng10
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns10+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor

          fracs(ns10+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(10) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(10) + 1

        indf = indfor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indfp = indf + 1

        do ig = 1, ng10
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns10+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)   &
     &            +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))  &
     &            + taufor

          fracs(ns10+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb10
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!!                          (high key - h2o; high minor - o2)
! ----------------------------------
      subroutine taugb11
! ..................................

!  ------------------------------------------------------------------  !
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)             !
!                              (high key - h2o; high minor - o2)       !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb11

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: scaleo2, tauself, taufor, tauo2
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(11) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(11) + 1

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1

        scaleo2 = colamt(k,6) * scaleminor(k)

        do ig = 1, ng11
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          tauo2   = scaleo2 * (ka_mo2(ig,indm) + minorfrac(k)           &
     &            * (ka_mo2(ig,indmp) - ka_mo2(ig,indm)))

          taug(ns11+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor + tauo2

          fracs(ns11+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(11) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(11) + 1

        indf = indfor(k)
        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indfp = indf + 1
        indmp = indm + 1

        scaleo2 = colamt(k,6) * scaleminor(k)

        do ig = 1, ng11
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))
          tauo2  = scaleo2 * (kb_mo2(ig,indm) + minorfrac(k)            &
     &           * (kb_mo2(ig,indmp) - kb_mo2(ig,indm)))

          taug(ns11+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)   &
     &            +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))  &
     &            + taufor + tauo2

          fracs(ns11+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb11
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
! ----------------------------------
      subroutine taugb12
! ..................................

!  ------------------------------------------------------------------  !
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb12

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, jpl, jplp,    &
     &       id000, id010, id100, id110, id200, id210, ig, js, js1,     &
     &       id001, id011, id101, id111, id201, id211

      real (kind=kind_phys) :: tauself, taufor, refrat_planck_a,        &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21
!
!===> ...  begin here
!
!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower/upper atmosphere.

      refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)      ! P =   174.164 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(12) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(12) + js1

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        if (specparm_planck >= oneminus) specparm_planck=oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indsp = inds + 1
        indfp = indf + 1
        jplp  = jpl  + 1

        if (specparm < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        if (specparm1 < 0.125) then
          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng12
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns12+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor

          fracs(ns12+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     *(fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        do ig = 1, ng12
          taug(ns12+ig,k) = f_zero
          fracs(ns12+ig,k) = f_zero
        enddo
      enddo

! ..................................
      end subroutine taugb12
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 13:  2080-2250 cm-1 (low key-h2o,n2o; high minor-o3 minor)
! ----------------------------------
      subroutine taugb13
! ..................................

!  ------------------------------------------------------------------  !
!     band 13:  2080-2250 cm-1 (low key-h2o,n2o; high minor-o3 minor)  !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb13

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmco2, jpl,      &
     &       id001, id011, id101, id111, id201, id211, jmco2p, jplp,    &
     &       jmco, jmcop, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, co2m1, co2m2, absco2,   &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mco2,  specparm_mco2,  specmult_mco2,  fmco2,     &
     &       speccomb_mco,   specparm_mco,   specmult_mco,   fmco,      &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, refrat_m_a3, ratco2,          &
     &       adjfac, adjcolco2, com1, com2, absco, abso3,               &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping levels :
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - co, p = 706 mb, t = 278.94 k
!     upper - o3, p = 95.5835 mb, t = 215.7 k

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower/upper atmosphere.

      refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)        ! P = 473.420 mb (Level 5)
      refrat_m_a = chi_mls(1,1)/chi_mls(4,1)             ! P = 1053. (Level 1)
      refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)            ! P = 706. (Level 3)

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,3,1)*colamt(k,4)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(13) + js

        speccomb1 = colamt(k,1) + rfrate(k,3,2)*colamt(k,4)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(13) + js1

        speccomb_mco2 = colamt(k,1) + refrat_m_a*colamt(k,4)
        specparm_mco2 = colamt(k,1) / speccomb_mco2
        specmult_mco2 = 8.0 * min(specparm_mco2, oneminus)
        jmco2 = 1 + int(specmult_mco2)
        fmco2 = mod(specmult_mco2, f_one)

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        speccomb_mco = colamt(k,1) + refrat_m_a3*colamt(k,4)
        specparm_mco = colamt(k,1) / speccomb_mco
        specmult_mco = 8.0 * min(specparm_mco, oneminus)
        jmco = 1 + int(specmult_mco)
        fmco = mod(specmult_mco, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,4)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmco2p= jmco2+ 1
        jmcop = jmco + 1

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * 3.55e-4
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.68
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        if (specparm < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        if (specparm1 < 0.125) then
          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng13
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          co2m1   = ka_mco2(ig,jmco2,indm) + fmco2                      &
     &            * (ka_mco2(ig,jmco2p,indm) - ka_mco2(ig,jmco2,indm))
          co2m2   = ka_mco2(ig,jmco2,indmp) + fmco2                     &
     &            * (ka_mco2(ig,jmco2p,indmp) - ka_mco2(ig,jmco2,indmp))
          absco2  = co2m1 + minorfrac(k) * (co2m2 - co2m1)
          com1    = ka_mco(ig,jmco,indm) + fmco                         &
     &            * (ka_mco(ig,jmcop,indm) - ka_mco(ig,jmco,indm))
          com2    = ka_mco(ig,jmco,indmp) + fmco                        &
     &            * (ka_mco(ig,jmcop,indmp) - ka_mco(ig,jmco,indmp))
          absco   = com1 + minorfrac(k) * (com2 - com1)

          taug(ns13+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor + adjcolco2*absco2             &
     &                + colamt(k,7)*absco

          fracs(ns13+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        indm = indminor(k)
        indmp = indm + 1

        do ig = 1, ng13
          abso3 = kb_mo3(ig,indm) + minorfrac(k)                        &
     &          * (kb_mo3(ig,indmp) - kb_mo3(ig,indm))

          taug(ns13+ig,k) = colamt(k,3)*abso3

          fracs(ns13+ig,k) =  fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb13
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 14:  2250-2380 cm-1 (low - co2; high - co2)
! ----------------------------------
      subroutine taugb14
! ..................................

!  ------------------------------------------------------------------  !
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)                 !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb14

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       ig

      real (kind=kind_phys) :: tauself, taufor
!
!===> ...  begin here
!
!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(14) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(14) + 1

        inds = indself(k)
        indf = indfor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1

        do ig = 1, ng14
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns14+ig,k) = colamt(k,2)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor

          fracs(ns14+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(14) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(14) + 1

        ind0p = ind0 + 1
        ind1p = ind1 + 1

        do ig = 1, ng14
          taug(ns14+ig,k) = colamt(k,2)                                 &
     &             * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)  &
     &             +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))

          fracs(ns14+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb14
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!!                          (high - nothing)
! ----------------------------------
      subroutine taugb15
! ..................................

!  ------------------------------------------------------------------  !
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)         !
!                              (high - nothing)                        !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb15

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jpl, jplp,       &
     &       id001, id011, id101, id111, id201, id211, jmn2, jmn2p,     &
     &       ig, js, js1

      real (kind=kind_phys) :: scalen2, tauself, taufor,                &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mn2,   specparm_mn2,   specmult_mn2,   fmn2,      &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, n2m1, n2m2, taun2,            &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - nitrogen continuum, P = 1053., T = 294.

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower atmosphere.

      refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)      ! P = 1053. mb (Level 1)
      refrat_m_a = chi_mls(4,1)/chi_mls(2,1)           ! P = 1053. mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,4) + rfrate(k,5,1)*colamt(k,2)
        specparm = colamt(k,4) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(15) + js

        speccomb1 = colamt(k,4) + rfrate(k,5,2)*colamt(k,2)
        specparm1 = colamt(k,4) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(15) + js1

        speccomb_mn2 = colamt(k,4) + refrat_m_a*colamt(k,2)
        specparm_mn2 = colamt(k,4) / speccomb_mn2
        specmult_mn2 = 8.0 * min(specparm_mn2, oneminus)
        jmn2 = 1 + int(specmult_mn2)
        fmn2 = mod(specmult_mn2, f_one)

        speccomb_planck = colamt(k,4) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,4) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        scalen2 = colbrd(k) * scaleminor(k)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmn2p = jmn2 + 1

        if (specparm < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        if (specparm1 < 0.125) then
          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng15
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          n2m1    = ka_mn2(ig,jmn2,indm) + fmn2                         &
     &            * (ka_mn2(ig,jmn2p,indm) - ka_mn2(ig,jmn2,indm))
          n2m2    = ka_mn2(ig,jmn2,indmp) + fmn2                        &
     &            * (ka_mn2(ig,jmn2p,indmp) - ka_mn2(ig,jmn2,indmp))
          taun2   = scalen2 * (n2m1 + minorfrac(k) * (n2m2 - n2m1))

          taug(ns15+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor + taun2

          fracs(ns15+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        do ig = 1, ng15
          taug(ns15+ig,k) = f_zero

          fracs(ns15+ig,k) = f_zero
        enddo
      enddo

! ..................................
      end subroutine taugb15
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
! ----------------------------------
      subroutine taugb16
! ..................................

!  ------------------------------------------------------------------  !
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)      !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb16

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, jpl, jplp,       &
     &       id001, id011, id101, id111, id201, id211, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, refrat_planck_a,        &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21
!
!===> ...  begin here
!
!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower atmosphere.

      refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)        ! P = 387. mb (Level 6)

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,4,1)*colamt(k,5)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(16) + js

        speccomb1 = colamt(k,1) + rfrate(k,4,2)*colamt(k,5)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(16) + js1

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,5)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indsp = inds + 1
        indfp = indf + 1
        jplp  = jpl  + 1

        if (specparm < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        if (specparm1 < 0.125) then
          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng16
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns16+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor

          fracs(ns16+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(16) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(16) + 1

        ind0p = ind0 + 1
        ind1p = ind1 + 1

        do ig = 1, ng16
          taug(ns16+ig,k) = colamt(k,5)                                 &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))

          fracs(ns16+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb16
! ----------------------------------

! ..................................
      end subroutine taumol
!! @}
!-----------------------------------

!mz* exponential cloud overlapping subroutines
!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------
! mz* - Add height needed for exponential and exponential-random cloud overlap methods (icld=4 and 5, respectively)
      subroutine mcica_subcol_lw(iplon, ncol, nlay, icld, permuteseed,  &
     &                 irng, play, hgt,                                 &
     &                 cldfrac, ciwp, clwp, cswp, rei, rel, res, tauc,  &
     &                 cldfmcl,                                         &
     &                 ciwpmcl, clwpmcl, cswpmcl, reicmcl, relqmcl,     &
     &                 resnmcl, taucmcl)

      use machine, only : im => kind_io4, rb => kind_phys
! ----- Input -----
! Control
      integer(kind=im), intent(in) :: iplon           ! column/longitude index
      integer(kind=im), intent(in) :: ncol            ! number of  columns
      integer(kind=im), intent(in) :: nlay            ! number of model layers
      integer(kind=im), intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer(kind=im), intent(in) :: permuteseed     ! if the cloud generator is called multiple times, 
                                                      ! permute the seed between each call.
                                                      ! between calls for LW and SW, recommended
                                                      ! permuteseed differes by 'ngpt'
      integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                      !  0 = kissvec
                                                      !  1 = Mersenne
                                                      !  Twister

! Atmosphere
      real(kind=rb), intent(in) :: play(:,:)          ! layer pressures (mb) 
                                                      !    Dimensions: (ncol,nlay)

! mji - Add height
      real(kind=rb), intent(in) :: hgt(:,:)           ! layer height (m)
                                                      !    Dimensions: (ncol,nlay)

! Atmosphere/clouds - cldprop
      real(kind=rb), intent(in) :: cldfrac(:,:)       ! layer cloud fraction
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: tauc(:,:,:)        ! in-cloud optical depth
                                                      !    Dimensions: (nbndlw,ncol,nlay)
!      real(kind=rb), intent(in) :: ssac(:,:,:)       ! in-cloud single scattering albedo
                                                      !    Dimensions: (nbndlw,ncol,nlay)
!      real(kind=rb), intent(in) :: asmc(:,:,:)       ! in-cloud asymmetry parameter
                                                      !    Dimensions: (nbndlw,ncol,nlay)
      real(kind=rb), intent(in) :: ciwp(:,:)          ! in-cloud ice water path
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: clwp(:,:)          ! in-cloud liquid water path
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: cswp(:,:)          ! in-cloud snow path
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: rei(:,:)           ! cloud ice particle size
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: rel(:,:)           ! cloud liquid particle size
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: res(:,:)           ! snow particle size
                                                      !    Dimensions: (ncol,nlay)

! ----- Output -----                                                                                                          
! Atmosphere/clouds - cldprmc [mcica]                                                                                                
      real(kind=rb), intent(out) :: cldfmcl(:,:,:)    ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(out) :: ciwpmcl(:,:,:)    ! in-cloud ice water path [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(out) :: clwpmcl(:,:,:)    ! in-cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(out) :: cswpmcl(:,:,:)    ! in-cloud snow path [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(out) :: relqmcl(:,:)      ! liquid particle size (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(out) :: reicmcl(:,:)      ! ice partcle size (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(out) :: resnmcl(:,:)      ! snow partcle size (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(out) :: taucmcl(:,:,:)    ! in-cloud optical depth [mcica]
!mz*
                                                      !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=rb), intent(out) :: ssacmcl(:,:,:)   ! in-cloud single scattering albedo [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=rb), intent(out) :: asmcmcl(:,:,:)   ! in-cloud asymmetry parameter [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
! ----- Local -----

! Stochastic cloud generator variables [mcica]
      integer(kind=im), parameter :: nsubclw = ngptlw ! number of sub-columns (g-point intervals)
      integer(kind=im) :: ilev                        ! loop index

      real(kind=rb) :: pmid(ncol, nlay)               ! layer pressures (Pa)
!      real(kind=rb) :: pdel(ncol, nlay)              ! layer pressure thickness (Pa)
!      real(kind=rb) :: qi(ncol, nlay)                ! ice water (specific humidity)
!      real(kind=rb) :: ql(ncol, nlay)                ! liq water (specific humidity)

! Return if clear sky
      if (icld.eq.0) return

! NOTE: For GCM mode, permuteseed must be offset between LW and SW by at least the number of subcolumns


! Pass particle sizes to new arrays, no subcolumns for these properties yet
! Convert pressures from mb to Pa

      reicmcl(:ncol,:nlay) = rei(:ncol,:nlay)
      relqmcl(:ncol,:nlay) = rel(:ncol,:nlay)
      resnmcl(:ncol,:nlay) = res(:ncol,:nlay)
      pmid(:ncol,:nlay) = play(:ncol,:nlay)*1.e2_rb

!  Generate the stochastic subcolumns of cloud optical properties for
!  the longwave
      call generate_stochastic_clouds (ncol, nlay, nsubclw, icld, irng, &
     &                      pmid, hgt, cldfrac, clwp, ciwp, cswp, tauc, &
     &                         cldfmcl, clwpmcl, ciwpmcl, cswpmcl,      &
     &                         taucmcl, permuteseed)

      end subroutine mcica_subcol_lw
!-------------------------------------------------------------------------------------------------
      subroutine generate_stochastic_clouds(ncol, nlay, nsubcol, icld,  &
     &                    irng, pmid, hgt, cld, clwp, ciwp, cswp, tauc, &
     &                             cld_stoch, clwp_stoch, ciwp_stoch,   &
     &                              cswp_stoch, tauc_stoch, changeSeed)  
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
! Contact: Cecile Hannay (hannay@ucar.edu)
!
! Original code: Based on Raisanen et al., QJRMS, 2004.
!
! Modifications:
!   1) Generalized for use with RRTMG and added Mersenne Twister as the default
!   random number generator, which can be changed to the optional kissvec random number generator
!   with flag 'irng'. Some extra functionality has been commented or removed.
!   Michael J. Iacono, AER, Inc., February 2007
!   2) Activated exponential and exponential/random cloud overlap method
!   Michael J. Iacono, AER, November 2017
!
! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one
! and uniform cloud liquid and cloud ice concentration.
! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer
! and obeys an overlap assumption in the vertical.
!
! Overlap assumption:
!  The cloud are consistent with 5 overlap assumptions: random, maximum, maximum-random, exponential and exponential random.
!  The default option is maximum-random (option 2)
!  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap, 5=exp/random 
!  This is set with the variable "overlap"
!  The exponential overlap uses also a length scale, Zo. (real,  parameter  :: Zo = 2500. )
!
! Seed:
!  If the stochastic cloud generator is called several times during the same timestep,
!  one should change the seed between the call to insure that the
!  subcolumns are different.                                        
!  This is done by changing the argument 'changeSeed'                                                                              
!  For example, if one wants to create a set of columns for the
!  shortwave and another set for the longwave ,
!  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call

! PDF assumption:
!  We can use arbitrary complicated PDFS.
!  In the present version, we produce homogeneuous clouds (the simplest case).
!  Future developments include using the PDF scheme of Ben Johnson.
!
! History file:
!  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
!  nsubcol = number of subcolumns
!  overlap = overlap type (1-3)
!  Zo = length scale                                               
!  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)  
!  CLDLIQ_S = mean of the subcolumn cloud water
!  CLDICE_S = mean of the subcolumn cloud ice
!
! Note:
!   Here: we force that the cloud condensate to be consistent with the cloud fraction
!   i.e we only have cloud condensate when the cell is cloudy.
!   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations
!   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction
!   without cloud condensate or the opposite).
!-----------------------------------------------------------------

      use mcica_random_numbers
! The Mersenne Twister random number engine
      use MersenneTwister, only: randomNumberSequence,                  &
     &                    new_RandomNumberSequence, getRandomReal
      use machine ,only : im => kind_io4, rb => kind_phys

      type(randomNumberSequence) :: randomNumbers

! -- Arguments

      integer(kind=im), intent(in) :: ncol            ! number of columns
      integer(kind=im), intent(in) :: nlay            ! number of layers
      integer(kind=im), intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                      !  0 = kissvec
                                                      !  1 = Mersenne Twister
      integer(kind=im), intent(in) :: nsubcol         ! number of sub-columns (g-point intervals)
      integer(kind=im), optional, intent(in) :: changeSeed     ! allows permuting seed 

! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state
      real(kind=rb), intent(in) :: pmid(:,:)          ! layer pressure (Pa)
                                                      !    Dimensions: (ncol,nlay)

      real(kind=rb), intent(in) :: hgt(:,:)           ! layer height (m)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: cld(:,:)           ! cloud fraction
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: clwp(:,:)          ! in-cloud liquid water path
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: ciwp(:,:)          ! in-cloud ice water path
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: cswp(:,:)          ! in-cloud snow path
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: tauc(:,:,:)        ! in-cloud optical depth
                                                      !    Dimensions:(nbndlw,ncol,nlay)
!      real(kind=rb), intent(in) :: ssac(:,:,:)       ! in-cloud single scattering albedo
                                                      !    Dimensions: (nbndlw,ncol,nlay)
                                                      !   inactive - for future expansion
!      real(kind=rb), intent(in) :: asmc(:,:,:)       ! in-cloud asymmetry parameter
                                                      !    Dimensions: (nbndlw,ncol,nlay)
                                                      !   inactive - for future expansion

      real(kind=rb), intent(out) :: cld_stoch(:,:,:)  ! subcolumn cloud fraction
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(out) :: clwp_stoch(:,:,:) ! subcolumn in-cloud liquid water path
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(out) :: ciwp_stoch(:,:,:) ! subcolumn in-cloud ice water path
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(out) :: cswp_stoch(:,:,:) ! subcolumn in-cloud snow path
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(out) :: tauc_stoch(:,:,:) ! subcolumn in-cloud optical depth
                                                      !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=rb), intent(out) :: ssac_stoch(:,:,:)! subcolumn in-cloud single scattering albedo
                                                      !    Dimensions: (ngptlw,ncol,nlay)
                                                      !   inactive - for future expansion
!      real(kind=rb), intent(out) :: asmc_stoch(:,:,:)! subcolumn in-cloud asymmetry parameter
                                                      !    Dimensions: (ngptlw,ncol,nlay)
                                                      !   inactive - for future expansion

! -- Local variables
      real(kind=rb) :: cldf(ncol,nlay)                ! cloud fraction

! Mean over the subcolumns (cloud fraction, cloud water , cloud ice) - inactive
!      real(kind=rb) :: mean_cld_stoch(ncol, nlay)    ! cloud fraction
!      real(kind=rb) :: mean_clwp_stoch(ncol, nlay)   ! cloud water
!      real(kind=rb) :: mean_ciwp_stoch(ncol, nlay)   ! cloud ice
!      real(kind=rb) :: mean_tauc_stoch(ncol, nlay)   ! cloud optical depth
!      real(kind=rb) :: mean_ssac_stoch(ncol, nlay)   ! cloud single scattering albedo 
!      real(kind=rb) :: mean_asmc_stoch(ncol, nlay)   ! cloud asymmetry parameter 

! Set overlap
      integer(kind=im) :: overlap                     ! 1 = random overlap, 2 = maximum-random,
                                                      ! 3 = maximum overlap, 4 = exponential, 
                                                      ! 5 = exponential-random
      real(kind=rb), parameter  :: Zo = 2500._rb      ! length scale (m)
      real(kind=rb), dimension(ncol,nlay) :: alpha    ! overlap parameter

! Constants (min value for cloud fraction and cloud water and ice)
      real(kind=rb), parameter :: cldmin = 1.0e-20_rb ! min cloud fraction
!      real(kind=rb), parameter :: qmin   = 1.0e-10_rb   ! min cloud water and cloud ice (not used)

! Variables related to random number and seed
      real(kind=rb), dimension(nsubcol, ncol, nlay) :: CDF, CDF2      !random numbers
      integer(kind=im), dimension(ncol) :: seed1, seed2, seed3, seed4 !seed to create random number (kissvec)
      real(kind=rb), dimension(ncol) :: rand_num      ! random number (kissvec)
      integer(kind=im) :: iseed                       ! seed to create random number (Mersenne Teister)
      real(kind=rb) :: rand_num_mt                    ! random number (Mersenne Twister)

! Flag to identify cloud fraction in subcolumns
      logical,  dimension(nsubcol, ncol, nlay) :: iscloudy   ! flag that says whether a gridbox is cloudy

! Indices
      integer(kind=im) :: ilev, isubcol, i, n         ! indices

!-------------------------------------------------------------------

! Check that irng is in bounds; if not, set to default
      if (irng .ne. 0) irng = 1

! Pass input cloud overlap setting to local variable
      overlap = icld

! Ensure that cloud fractions are in bounds
      do ilev = 1, nlay
         do i = 1, ncol
            cldf(i,ilev) = cld(i,ilev)
            if (cldf(i,ilev) < cldmin) then
               cldf(i,ilev) = 0._rb
            endif
         enddo
      enddo

! ----- Create seed  --------

! Advance randum number generator by changeseed values
      if (irng.eq.0) then 
! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.
! Must use pmid from bottom four layers.
         do i=1,ncol
            if (pmid(i,1).lt.pmid(i,2)) then
               stop 'MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID &
     &               FROM BOTTOM FOUR LAYERS.'
            endif
            seed1(i) = (pmid(i,1) - int(pmid(i,1)))  * 1000000000_im
            seed2(i) = (pmid(i,2) - int(pmid(i,2)))  * 1000000000_im
            seed3(i) = (pmid(i,3) - int(pmid(i,3)))  * 1000000000_im
            seed4(i) = (pmid(i,4) - int(pmid(i,4)))  * 1000000000_im
          enddo
         do i=1,changeSeed
            call kissvec(seed1, seed2, seed3, seed4, rand_num)
         enddo
      elseif (irng.eq.1) then
         randomNumbers = new_RandomNumberSequence(seed = changeSeed)
      endif

! ------ Apply overlap assumption --------

! generate the random numbers

      select case (overlap) 

      case(1)
! Random overlap
! i) pick a random value at every level

         if (irng.eq.0) then
            do isubcol = 1,nsubcol
               do ilev = 1,nlay
                  call kissvec(seed1, seed2, seed3, seed4, rand_num)  ! we get different random number for each level
                  CDF(isubcol,:,ilev) = rand_num
               enddo    
            enddo
         elseif (irng.eq.1) then
            do isubcol = 1, nsubcol
               do i = 1, ncol
                  do ilev = 1, nlay
                     rand_num_mt = getRandomReal(randomNumbers)
                     CDF(isubcol,i,ilev) = rand_num_mt
                  enddo
               enddo
             enddo
         endif

      case(2)
! Maximum-Random overlap 
! i) pick a random number for top layer.
! ii) walk down the column:
!    - if the layer above is cloudy, we use the same random number than in the layer above
!    - if the layer above is clear, we use a new random number

         if (irng.eq.0) then
            do isubcol = 1,nsubcol
               do ilev = 1,nlay
                  call kissvec(seed1, seed2, seed3, seed4, rand_num)
                  CDF(isubcol,:,ilev) = rand_num
               enddo
            enddo
         elseif (irng.eq.1) then
            do isubcol = 1, nsubcol
               do i = 1, ncol
                  do ilev = 1, nlay
                     rand_num_mt = getRandomReal(randomNumbers)
                     CDF(isubcol,i,ilev) = rand_num_mt
                  enddo
               enddo
             enddo
         endif

         do ilev = 2,nlay
            do i = 1, ncol
               do isubcol = 1, nsubcol
                  if (CDF(isubcol, i, ilev-1) > 1._rb - cldf(i,ilev-1) )&
     &             then
                     CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev-1)
                  else
                     CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev) * (1._rb &
     &               - cldf(i,ilev-1))
                  endif
               enddo
            enddo
         enddo

      case(3)
! Maximum overlap
! i) pick the same random numebr at every level

         if (irng.eq.0) then
            do isubcol = 1,nsubcol
               call kissvec(seed1, seed2, seed3, seed4, rand_num)
               do ilev = 1,nlay
                  CDF(isubcol,:,ilev) = rand_num
               enddo
            enddo
         elseif (irng.eq.1) then
            do isubcol = 1, nsubcol
               do i = 1, ncol
                  rand_num_mt = getRandomReal(randomNumbers)
                  do ilev = 1, nlay
                     CDF(isubcol,i,ilev) = rand_num_mt
                  enddo
               enddo
             enddo
         endif

! mji - Activate exponential cloud overlap option
         case(4)
            ! Exponential overlap: weighting between maximum and random overlap increases with the distance.
            ! The random numbers for exponential overlap verify:
            ! j=1   RAN(j)=RND1
            ! j>1   if RND1 < alpha(j,j-1) => RAN(j) = RAN(j-1)
            !                                 RAN(j) = RND2
            ! alpha is obtained from the equation
            ! alpha = exp(-(Z(j)-Z(j-1))/Zo) where Zo is a characteristic length scale

            ! compute alpha
            do i = 1, ncol
               alpha(i, 1) = 0._rb
               do ilev = 2,nlay
                  alpha(i, ilev) = exp( -( hgt (i, ilev) -              &
     &                  hgt (i, ilev-1)) / Zo)
               enddo
            enddo

            ! generate 2 streams of random numbers
            if (irng.eq.0) then
               do isubcol = 1,nsubcol
                  do ilev = 1,nlay
                     call kissvec(seed1, seed2, seed3, seed4, rand_num)
                     CDF(isubcol, :, ilev) = rand_num
                     call kissvec(seed1, seed2, seed3, seed4, rand_num)
                     CDF2(isubcol, :, ilev) = rand_num
                  enddo
               enddo
            elseif (irng.eq.1) then
               do isubcol = 1, nsubcol
                  do i = 1, ncol
                     do ilev = 1, nlay
                        rand_num_mt = getRandomReal(randomNumbers)
                        CDF(isubcol,i,ilev) = rand_num_mt
                        rand_num_mt = getRandomReal(randomNumbers)
                        CDF2(isubcol,i,ilev) = rand_num_mt
                     enddo
                  enddo
               enddo
            endif

            ! generate random numbers
            do ilev = 2,nlay
               where (CDF2(:, :, ilev) < spread(alpha (:,ilev),         &
     &               dim=1,nCopies=nsubcol) )
                  CDF(:,:,ilev) = CDF(:,:,ilev-1)
               end where
            end do

! Activate exponential-random cloud overlap option
         case(5)
            ! Exponential-random overlap:
!mz*            call wrf_error_fatal("Cloud Overlap case 5: ER has not yet  &
!                            been implemented. Stopping...")

      end select

! -- generate subcolumns for homogeneous clouds ----- 
      do ilev = 1,nlay 
         iscloudy(:,:,ilev) = (CDF(:,:,ilev) >= 1._rb -                 &
     &        spread(cldf(:,ilev), dim=1, nCopies=nsubcol) ) 
      enddo

! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
! where there is a cloud, define the subcolumn cloud properties,
! otherwise set these to zero

      do ilev = 1,nlay
         do i = 1, ncol
            do isubcol = 1, nsubcol
               if (iscloudy(isubcol,i,ilev) ) then
                  cld_stoch(isubcol,i,ilev) = 1._rb
                  clwp_stoch(isubcol,i,ilev) = clwp(i,ilev)
                  ciwp_stoch(isubcol,i,ilev) = ciwp(i,ilev)
!mz  
!                  cswp_stoch(isubcol,i,ilev) = cswp(i,ilev)
                   cswp_stoch(isubcol,i,ilev) = 0._rb
                  n = ngb(isubcol)
                  tauc_stoch(isubcol,i,ilev) = tauc(n,i,ilev)
!                  ssac_stoch(isubcol,i,ilev) = ssac(n,i,ilev)
!                  asmc_stoch(isubcol,i,ilev) = asmc(n,i,ilev)
               else
                  cld_stoch(isubcol,i,ilev) = 0._rb
                  clwp_stoch(isubcol,i,ilev) = 0._rb
                  ciwp_stoch(isubcol,i,ilev) = 0._rb
                  cswp_stoch(isubcol,i,ilev) = 0._rb
                  tauc_stoch(isubcol,i,ilev) = 0._rb
!                  ssac_stoch(isubcol,i,ilev) = 1._rb
!                  asmc_stoch(isubcol,i,ilev) = 1._rb
               endif
            enddo
         enddo
      enddo

! -- compute the means of the subcolumns ---
!      mean_cld_stoch(:,:) = 0._rb
!      mean_clwp_stoch(:,:) = 0._rb
!      mean_ciwp_stoch(:,:) = 0._rb
!      mean_tauc_stoch(:,:) = 0._rb
!      mean_ssac_stoch(:,:) = 0._rb
!      mean_asmc_stoch(:,:) = 0._rb
!      do i = 1, nsubcol
!         mean_cld_stoch(:,:) =  cld_stoch(i,:,:) + mean_cld_stoch(:,:) 
!         mean_clwp_stoch(:,:) =  clwp_stoch( i,:,:) + mean_clwp_stoch(:,:)
!         mean_ciwp_stoch(:,:) =  ciwp_stoch( i,:,:) + mean_ciwp_stoch(:,:)
!         mean_tauc_stoch(:,:) =  tauc_stoch( i,:,:) + mean_tauc_stoch(:,:)
!         mean_ssac_stoch(:,:) =  ssac_stoch( i,:,:) + mean_ssac_stoch(:,:)
!         mean_asmc_stoch(:,:) =  asmc_stoch( i,:,:) + mean_asmc_stoch(:,:)
!      end do 
!      mean_cld_stoch(:,:) = mean_cld_stoch(:,:) / nsubcol
!      mean_clwp_stoch(:,:) = mean_clwp_stoch(:,:) / nsubcol
!      mean_ciwp_stoch(:,:) = mean_ciwp_stoch(:,:) / nsubcol
!      mean_tauc_stoch(:,:) = mean_tauc_stoch(:,:) / nsubcol
!      mean_ssac_stoch(:,:) = mean_ssac_stoch(:,:) / nsubcol
!      mean_asmc_stoch(:,:) = mean_asmc_stoch(:,:) / nsubcol

      end subroutine generate_stochastic_clouds

!------------------------------------------------------------------
! Private subroutines
!------------------------------------------------------------------  

!----------------------------------------------------------------- 
      subroutine kissvec(seed1,seed2,seed3,seed4,ran_arr) 
!----------------------------------------------------------------

! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04 for NCAR CAM
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004

! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59  
!  Overall period>2^123;                                                                                                             
      real(kind=rb), dimension(:), intent(inout)  :: ran_arr
      integer(kind=im), dimension(:), intent(inout) :: seed1,seed2,seed3&
     &                                                 ,seed4
      integer(kind=im) :: i,sz,kiss
      integer(kind=im) :: m, k, n

! inline function  
      m(k, n) = ieor (k, ishft (k, n) )

      sz = size(ran_arr)
      do i = 1, sz 
         seed1(i) = 69069_im * seed1(i) + 1327217885_im
         seed2(i) = m (m (m (seed2(i), 13_im), - 17_im), 5_im)
         seed3(i) = 18000_im * iand (seed3(i), 65535_im) +              &
     &              ishft (seed3(i), - 16_im)
         seed4(i) = 30903_im * iand (seed4(i), 65535_im) +              &
     &              ishft (seed4(i), - 16_im)
         kiss = seed1(i) + seed2(i) + ishft (seed3(i), 16_im) + seed4(i)
         ran_arr(i) = kiss*2.328306e-10_rb + 0.5_rb
      end do 

      end subroutine kissvec
!
      subroutine rtrnmc_mcica(nlayers, istart, iend, iout, pz, semiss,  &
     &       ncbands,  cldfmc, taucmc, planklay, planklev,              &!plankbnd,    &
     &       pwvcm, fracs, taut,                                        &
     &                   totuflux, totdflux,  htr,                      &
     &                   totuclfl, totdclfl,  htrc )
!---------------------------------------------------------------
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90:  Michael J. Iacono; June, 2006
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity 
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
!  use of the emissivity angle for the flux integration can cause errors of 
!  1 to 4 W/m2 within cloudy layers.  
!  Clouds are treated with the McICA stochastic approach and maximum-random               
!  cloud overlap.                                                                         
!***************************************************************************              
                                                                                          
! ------- Declarations -------                                                            
                                                                                          
! ----- Input -----                                                                       
      integer(kind=im), intent(in) :: nlayers         ! total number of layers            
      integer(kind=im), intent(in) :: istart          ! beginning band of calculation     
      integer(kind=im), intent(in) :: iend            ! ending band of calculation        
      integer(kind=im), intent(in) :: iout            ! output option flag                
                                                                                          
! Atmosphere                                                                              
      real(kind=rb), intent(in) :: pz(0:)             ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlayers)        
      real(kind=rb), intent(in) :: pwvcm              ! precipitable water vapor (cm)     
      real(kind=rb), intent(in) :: semiss(:)          ! lw surface emissivity             
                                                      !    Dimensions: (nbndlw)           
!mz
      real(kind=rb), intent(in) :: planklay(0:,:)      !                                   
                                                      !    Dimensions: (nlayers,nbndlw)   
      real(kind=rb), intent(in) :: planklev(0:,:)     !                                   
                                                      !    Dimensions: (0:nlayers,nbndlw) 
!      real(kind=rb), intent(in) :: plankbnd(:)        !                                   
                                                      !    Dimensions: (nbndlw)           
      real(kind=rb), intent(in) :: fracs(:,:)         !                                   
                                                      !    Dimensions: (nlayers,ngptw)    
      real(kind=rb), intent(in) :: taut(:,:)          ! gaseous + aerosol optical depths  
                                                      !    Dimensions: (nlayers,ngptlw)   
                                                                                          
! Clouds                                                                                  
      integer(kind=im), intent(in) :: ncbands         ! number of cloud spectral bands    
      real(kind=rb), intent(in) :: cldfmc(:,:)        ! layer cloud fraction [mcica]      
                                                      !    Dimensions: (ngptlw,nlayers)   
      real(kind=rb), intent(in) :: taucmc(:,:)        ! layer cloud optical depth [mcica] 
                                                      !    Dimensions: (ngptlw,nlayers)   
                                                                                          
! ----- Output -----                                                                      
      real(kind=rb), intent(out) :: totuflux(0:)      ! upward longwave flux (w/m2)       
                                                      !    Dimensions: (0:nlayers)        
      real(kind=rb), intent(out) :: totdflux(0:)      ! downward longwave flux (w/m2)     
                                                      !    Dimensions: (0:nlayers)        
!mz* real(kind=rb), intent(out) :: fnet(0:)          ! net longwave flux (w/m2)          
                                                      !    Dimensions: (0:nlayers)        
         real(kind=rb), intent(out) :: htr(:)
!mz      real(kind=rb), intent(out) :: htr(0:)           ! longwave heating rate (k/day)     
                                                      !    Dimensions: (0:nlayers)        
      real(kind=rb), intent(out) :: totuclfl(0:)      ! clear sky upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)        
      real(kind=rb), intent(out) :: totdclfl(0:)      ! clear sky downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)        
!mz*real(kind=rb), intent(out) :: fnetc(0:)         ! clear sky net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)        
       real(kind=rb), intent(out) :: htrc(:) 
!      real(kind=rb), intent(out) :: htrc(0:)          ! clear sky longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)        
                                                                                          
! ----- Local -----                                                                       
! Declarations for radiative transfer                                                     
      real (kind=kind_phys), dimension(0:nlayers) :: fnet, fnetc
      real(kind=rb) :: abscld(nlayers,ngptlw)                                             
      real(kind=rb) :: atot(nlayers)                                                      
      real(kind=rb) :: atrans(nlayers)                                                    
      real(kind=rb) :: bbugas(nlayers)                                                    
      real(kind=rb) :: bbutot(nlayers)                                                    
      real(kind=rb) :: clrurad(0:nlayers)                                                 
      real(kind=rb) :: clrdrad(0:nlayers)                                                 
      real(kind=rb) :: efclfrac(nlayers,ngptlw)                                           
      real(kind=rb) :: uflux(0:nlayers)                                                   
      real(kind=rb) :: dflux(0:nlayers)                                                   
      real(kind=rb) :: urad(0:nlayers)                                                    
      real(kind=rb) :: drad(0:nlayers)                                                    
      real(kind=rb) :: uclfl(0:nlayers)                                                   
      real(kind=rb) :: dclfl(0:nlayers)                                                   
      real(kind=rb) :: odcld(nlayers,ngptlw)                                              
                                                                                          
                                                                                          
      real(kind=rb) :: secdiff(nbands)                 ! secant of diffusivity angle      
      real(kind=rb) :: transcld, radld, radclrd, plfrac, blay, dplankup,&
     &                 dplankdn         
      real(kind=rb) :: odepth, odtot, odepth_rec, odtot_rec, gassrc                       
      real(kind=rb) :: tblind, tfactot, bbd, bbdtot, tfacgas, transc,   &
     &                 tausfac             
      real(kind=rb) :: rad0, reflect, radlu, radclru                                      
                                                                                          
      integer(kind=im) :: icldlyr(nlayers)                  ! flag for cloud in layer     
      integer(kind=im) :: ibnd, ib, iband, lay, lev, l, ig  ! loop indices                
      integer(kind=im) :: igc                               ! g-point interval counter    
      integer(kind=im) :: iclddn                            ! flag for cloud in down path 
      integer(kind=im) :: ittot, itgas, itr                 ! lookup table indices        
!mz*
      real (kind=kind_phys), parameter :: rec_6 = 0.166667
      ! The cumulative sum of new g-points for each band
      integer(kind=im) :: ngs(nbands)
      ngs(:) = (/10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,   &
     &          140/)
                                                                                          
! ------- Definitions -------                                                             
! input                                                                                   
!    nlayers                      ! number of model layers                                
!    ngptlw                       ! total number of g-point subintervals                  
!    nbndlw                       ! number of longwave spectral bands                     
!    ncbands                      ! number of spectral bands for clouds                   
!    secdiff                      ! diffusivity angle                                     
!    wtdiff                       ! weight for radiance to flux conversion                
!    pavel                        ! layer pressures (mb)                                  
!    pz                           ! level (interface) pressures (mb)                      
!    tavel                        ! layer temperatures (k)                                
!    tz                           ! level (interface) temperatures(mb)                    
!    tbound                       ! surface temperature (k)                               
!    cldfrac                      ! layer cloud fraction                                  
!    taucloud                     ! layer cloud optical depth                             
!    itr                          ! integer look-up table index                           
!    icldlyr                      ! flag for cloudy layers                                
!    iclddn                       ! flag for cloud in column at any layer                 
!    semiss                       ! surface emissivities for each band                    
!    reflect                      ! surface reflectance                                   
!    bpade                        ! 1/(pade constant)                                     
!    tau_tbl                      ! clear sky optical depth look-up table                 
!    exp_tbl                      ! exponential look-up table for transmittance           
!    tfn_tbl                      ! tau transition function look-up table                 
                                                                                          
! local                                                                                   
!    atrans                       ! gaseous absorptivity                                  
!    abscld                       ! cloud absorptivity                                    
!    atot                         ! combined gaseous and cloud absorptivity               
!    odclr                        ! clear sky (gaseous) optical depth                     
!    odcld                        ! cloud optical depth                                   
!    odtot                        ! optical depth of gas and cloud                        
!    tfacgas                      ! gas-only pade factor, used for planck fn              
!    tfactot                      ! gas and cloud pade factor, used for planck fn         
!    bbdgas                       ! gas-only planck function for downward rt              
!    bbugas                       ! gas-only planck function for upward rt                
!    bbdtot                       ! gas and cloud planck function for downward rt         
!    bbutot                       ! gas and cloud planck function for upward calc.        
!    gassrc                       ! source radiance due to gas only                       
!    efclfrac                     ! effective cloud fraction                              
!    radlu                        ! spectrally summed upward radiance                     
!    radclru                      ! spectrally summed clear sky upward radiance           
!    urad                         ! upward radiance by layer                              
!    clrurad                      ! clear sky upward radiance by layer                    
!    radld                        ! spectrally summed downward radiance                   
!    radclrd                      ! spectrally summed clear sky downward radiance         
!    drad                         ! downward radiance by layer                            
!    clrdrad                      ! clear sky downward radiance by layer                  
                                                                                          
                                                                                          
! output                                                                                  
!    totuflux                     ! upward longwave flux (w/m2)                           
!    totdflux                     ! downward longwave flux (w/m2)                         
!    fnet                         ! net longwave flux (w/m2)                              
!    htr                          ! longwave heating rate (k/day)                         
!    totuclfl                     ! clear sky upward longwave flux (w/m2)                 
!    totdclfl                     ! clear sky downward longwave flux (w/m2)               
!    fnetc                        ! clear sky net longwave flux (w/m2)                    
!    htrc                         ! clear sky longwave heating rate (k/day)               
                                                                                          
                                                                                          
!jm not thread safe      hvrrtc = '$Revision: 1.3 $'                                      
                                                                                          
      do ibnd = 1,nbands!mz*nbndlw                                                                  
         if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then                               
           secdiff(ibnd) = 1.66_rb                                                        
         else                                                                             
           secdiff(ibnd) = a0(ibnd) + a1(ibnd)*exp(a2(ibnd)*pwvcm)                        
           if (secdiff(ibnd) .gt. 1.80_rb) secdiff(ibnd) = 1.80_rb                        
           if (secdiff(ibnd) .lt. 1.50_rb) secdiff(ibnd) = 1.50_rb                        
         endif                                                                            
      enddo                                                                               
                                                                                          
      urad(0) = 0.0_rb                                                                    
      drad(0) = 0.0_rb                                                                    
      totuflux(0) = 0.0_rb                                                                
      totdflux(0) = 0.0_rb                                                                
      clrurad(0) = 0.0_rb                                                                 
      clrdrad(0) = 0.0_rb                                                                 
      totuclfl(0) = 0.0_rb                                                                
      totdclfl(0) = 0.0_rb                                                                
                                                                                          
      do lay = 1, nlayers                                                                 
         urad(lay) = 0.0_rb                                                               
         drad(lay) = 0.0_rb                                                               
         totuflux(lay) = 0.0_rb                                                           
         totdflux(lay) = 0.0_rb                                                           
         clrurad(lay) = 0.0_rb                                                            
         clrdrad(lay) = 0.0_rb                                                            
         totuclfl(lay) = 0.0_rb                                                           
         totdclfl(lay) = 0.0_rb                                                           
         icldlyr(lay) = 0                                                                 
                                                                 
! Change to band loop?                                                                    
         do ig = 1, ngptlw                                                                
            if (cldfmc(ig,lay) .eq. 1._rb) then                                           
               ib = ngb(ig)                                                               
               odcld(lay,ig) = secdiff(ib) * taucmc(ig,lay)                               
               transcld = exp(-odcld(lay,ig))                                             
               abscld(lay,ig) = 1._rb - transcld                                          
               efclfrac(lay,ig) = abscld(lay,ig) * cldfmc(ig,lay)                         
               icldlyr(lay) = 1                                                           
            else                                                                          
               odcld(lay,ig) = 0.0_rb                                                     
               abscld(lay,ig) = 0.0_rb                                                    
               efclfrac(lay,ig) = 0.0_rb                                                  
            endif                                                                         
         enddo                                                                            
                                                                                          
      enddo                                                                               
                                                                                          
      igc = 1                                                                             
! Loop over frequency bands.                                                              
      do iband = istart, iend                                                             
                                                                                          
! Reinitialize g-point counter for each band if output for each band is requested.        
         if (iout.gt.0.and.iband.ge.2) igc = ngs(iband-1)+1                               
                                                                                          
! Loop over g-channels.                                                                   
 1000    continue                                                                         
                                                                                          
! Radiative transfer starts here.                                                         
         radld = 0._rb                                                                    
         radclrd = 0._rb                                                                  
         iclddn = 0                                                                       
                                                                                          
! Downward radiative transfer loop.                                                       
                                                                                          
         do lev = nlayers, 1, -1                                                          
               plfrac = fracs(lev,igc)                                                    
               blay = planklay(lev,iband)                                                 
               dplankup = planklev(lev,iband) - blay                                      
               dplankdn = planklev(lev-1,iband) - blay                                    
               odepth = secdiff(iband) * taut(lev,igc)                                    
               if (odepth .lt. 0.0_rb) odepth = 0.0_rb              
!  Cloudy layer                                                                  
               if (icldlyr(lev).eq.1) then  
                  iclddn = 1              
                  odtot = odepth + odcld(lev,igc) 
                  if (odtot .lt. 0.06_rb) then            
                     atrans(lev) = odepth - 0.5_rb*odepth*odepth                          
                     odepth_rec = rec_6*odepth
                 gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)               
                                                                                          
                     atot(lev) =  odtot - 0.5_rb*odtot*odtot                              
                     odtot_rec = rec_6*odtot                                              
                     bbdtot =  plfrac * (blay+dplankdn*odtot_rec)                         
                     bbd = plfrac*(blay+dplankdn*odepth_rec)                              
                     radld = radld - radld * (atrans(lev) +             &
     &                    efclfrac(lev,igc) * (1. - atrans(lev))) +     &
     &                    gassrc + cldfmc(igc,lev) *                    &
     &                    (bbdtot * atot(lev) - gassrc)                                    
                     drad(lev-1) = drad(lev-1) + radld                                    
                                                                                          
                     bbugas(lev) =  plfrac * (blay+dplankup*odepth_rec)                   
                     bbutot(lev) =  plfrac * (blay+dplankup*odtot_rec)                    
                                                                                          
                  elseif (odepth .le. 0.06_rb) then                                       
                     atrans(lev) = odepth - 0.5_rb*odepth*odepth                          
                     odepth_rec = rec_6*odepth                                            
                 gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)               
                                                                                          
                     odtot = odepth + odcld(lev,igc)                                      
                     tblind = odtot/(bpade+odtot)                                         
                     ittot = tblint*tblind + 0.5_rb                                       
                     tfactot = tfn_tbl(ittot)                                             
                     bbdtot = plfrac * (blay + tfactot*dplankdn)                          
                     bbd = plfrac*(blay+dplankdn*odepth_rec)                              
                     atot(lev) = 1. - exp_tbl(ittot)                                      
                                                                                          
                     radld = radld - radld * (atrans(lev) +             &
     &                   efclfrac(lev,igc) * (1._rb - atrans(lev))) +   &
     &                   gassrc + cldfmc(igc,lev) *                     &
     &                   (bbdtot * atot(lev) - gassrc)             
                     drad(lev-1) = drad(lev-1) + radld                                    
                                                                                          
                     bbugas(lev) = plfrac * (blay + dplankup*odepth_rec)                  
                     bbutot(lev) = plfrac * (blay + tfactot * dplankup)                   
                                                                                          
                  else                                                                    
                                                                                          
                     tblind = odepth/(bpade+odepth)                                       
                     itgas = tblint*tblind+0.5_rb                                         
                     odepth = tau_tbl(itgas)                                              
                     atrans(lev) = 1._rb - exp_tbl(itgas)                                 
                     tfacgas = tfn_tbl(itgas)                                             
              gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)            
                                                                                          
                     odtot = odepth + odcld(lev,igc)                                      
                     tblind = odtot/(bpade+odtot)                                         
                     ittot = tblint*tblind + 0.5_rb                                       
                     tfactot = tfn_tbl(ittot)                                             
                     bbdtot = plfrac * (blay + tfactot*dplankdn)                          
                     bbd = plfrac*(blay+tfacgas*dplankdn)                                 
                     atot(lev) = 1._rb - exp_tbl(ittot)                                   
                                                                                          
                  radld = radld - radld * (atrans(lev) +                &
     &               efclfrac(lev,igc) * (1._rb - atrans(lev))) +       &
     &               gassrc + cldfmc(igc,lev) *                         &
     &               (bbdtot * atot(lev) - gassrc)                                         
                  drad(lev-1) = drad(lev-1) + radld                                       
                  bbugas(lev) = plfrac * (blay + tfacgas * dplankup)                      
                  bbutot(lev) = plfrac * (blay + tfactot * dplankup)                      
                  endif                                                                   
!  Clear layer                                                                            
               else                                                                       
                  if (odepth .le. 0.06_rb) then                                           
                     atrans(lev) = odepth-0.5_rb*odepth*odepth                            
                     odepth = rec_6*odepth                                                
                     bbd = plfrac*(blay+dplankdn*odepth)                                  
                     bbugas(lev) = plfrac*(blay+dplankup*odepth)                          
                  else                                                                    
                     tblind = odepth/(bpade+odepth)                                       
                     itr = tblint*tblind+0.5_rb                                           
                     transc = exp_tbl(itr)                                                
                     atrans(lev) = 1._rb-transc                                           
                     tausfac = tfn_tbl(itr)                                               
                     bbd = plfrac*(blay+tausfac*dplankdn)                                 
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)                   
                  endif                                                                   
                  radld = radld + (bbd-radld)*atrans(lev)                                 
                  drad(lev-1) = drad(lev-1) + radld                                       
               endif                                                                      
!  Set clear sky stream to total sky stream as long as layers                             
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),                     
!  and clear sky stream must be computed separately from that point.                      
                  if (iclddn.eq.1) then                                                   
                     radclrd = radclrd + (bbd-radclrd) * atrans(lev)                      
                     clrdrad(lev-1) = clrdrad(lev-1) + radclrd                            
                  else                                                                    
                     radclrd = radld                                                      
                     clrdrad(lev-1) = drad(lev-1)                                         
                  endif                                                                   
            enddo                                                                         
                                                                                          
! Spectral emissivity & reflectance                                                       
!  Include the contribution of spectrally varying longwave emissivity                     
!  and reflection from the surface to the upward radiative transfer.                      
!  Note: Spectral and Lambertian reflection are identical for the                         
!  diffusivity angle flux integration used here.                                          
                                                                                          
!mz*
!         rad0 = fracs(1,igc) * plankbnd(iband)                                            
          rad0 = semiss(iband) * fracs(1,igc) * planklay(0,iband)
!mz
!  Add in specular reflection of surface downward radiance.                               
         reflect = 1._rb - semiss(iband)                                                  
         radlu = rad0 + reflect * radld                                                   
         radclru = rad0 + reflect * radclrd                                               
                                                                                          
                                                                                          
! Upward radiative transfer loop.                                                         
         urad(0) = urad(0) + radlu                                                        
         clrurad(0) = clrurad(0) + radclru                                                
                                                                                          
         do lev = 1, nlayers                                                              
!  Cloudy layer                                                                           
            if (icldlyr(lev) .eq. 1) then 
               gassrc = bbugas(lev) * atrans(lev)
               radlu = radlu - radlu * (atrans(lev) +                   &
     &             efclfrac(lev,igc) * (1._rb - atrans(lev))) +         &
     &              gassrc + cldfmc(igc,lev) *                          &
     &              (bbutot(lev) * atot(lev) - gassrc)                                     
               urad(lev) = urad(lev) + radlu                                              
!  Clear layer                                                                            
            else                                                                          
               radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)                            
               urad(lev) = urad(lev) + radlu                                              
            endif                                                                         
!  Set clear sky stream to total sky stream as long as all layers                         
!  are clear (iclddn=0).  Streams must be calculated separately at                        
!  all layers when a cloud is present (ICLDDN=1), because surface                         
!  reflectance is different for each stream.                         
               if (iclddn.eq.1) then                                                      
                  radclru = radclru + (bbugas(lev)-radclru)*atrans(lev)                   
                  clrurad(lev) = clrurad(lev) + radclru                                   
               else                                                                       
                  radclru = radlu                                                         
                  clrurad(lev) = urad(lev)                                                
               endif                                                                      
         enddo                                                                            
                                                                                          
! Increment g-point counter                                                               
         igc = igc + 1                                                                    
! Return to continue radiative transfer for all g-channels in present band                
         if (igc .le. ngs(iband)) go to 1000                                              
                                                                                          
! Process longwave output from band for total and clear streams.                          
! Calculate upward, downward, and net flux.                                               
         do lev = nlayers, 0, -1                                                          
            uflux(lev) = urad(lev)*wtdiff                                                 
            dflux(lev) = drad(lev)*wtdiff                                                 
            urad(lev) = 0.0_rb                                                            
            drad(lev) = 0.0_rb                                                            
            totuflux(lev) = totuflux(lev) + uflux(lev) * delwave(iband)                   
            totdflux(lev) = totdflux(lev) + dflux(lev) * delwave(iband)                   
            uclfl(lev) = clrurad(lev)*wtdiff                                              
            dclfl(lev) = clrdrad(lev)*wtdiff                                              
            clrurad(lev) = 0.0_rb                                                         
            clrdrad(lev) = 0.0_rb                                                         
            totuclfl(lev) = totuclfl(lev) + uclfl(lev) * delwave(iband)                   
            totdclfl(lev) = totdclfl(lev) + dclfl(lev) * delwave(iband)                   
         enddo                                                                            
                                                                                          
! End spectral band loop                                                                  
      enddo                                                                               
                                                                                          
! Calculate fluxes at surface                                                             
      totuflux(0) = totuflux(0) * fluxfac                                                 
      totdflux(0) = totdflux(0) * fluxfac                                                 
      fnet(0) = totuflux(0) - totdflux(0)                                                 
      totuclfl(0) = totuclfl(0) * fluxfac                                                 
      totdclfl(0) = totdclfl(0) * fluxfac                                                 
      fnetc(0) = totuclfl(0) - totdclfl(0)                                                
                                                                                          
! Calculate fluxes at model levels                                                        
      do lev = 1, nlayers                                                                 
         totuflux(lev) = totuflux(lev) * fluxfac                                          
         totdflux(lev) = totdflux(lev) * fluxfac                                          
         fnet(lev) = totuflux(lev) - totdflux(lev)                                        
         totuclfl(lev) = totuclfl(lev) * fluxfac                                          
         totdclfl(lev) = totdclfl(lev) * fluxfac                                          
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)                                       
         l = lev - 1                                                                      
                                                                          
! Calculate heating rates at model layers                                                 
         htr(l)=heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev))                               
         htrc(l)=heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev))                            
      enddo                                                                               
                                                                                          
! Set heating rate to zero in top layer                                                   
      htr(nlayers) = 0.0_rb                                                               
      htrc(nlayers) = 0.0_rb                                                              
                                                                                          
      end subroutine rtrnmc_mcica             

! ------------------------------------------------------------------------------
      subroutine cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc,     &
     &  ciwpmc, clwpmc, cswpmc, reicmc, relqmc, resnmc, ncbands, taucmc, errmsg, errflg)
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers         ! total number of layers
      integer(kind=im), intent(in) :: inflag          ! see definitions
      integer(kind=im), intent(in) :: iceflag         ! see definitions
      integer(kind=im), intent(in) :: liqflag         ! see definitions

      real(kind=rb), intent(in) :: cldfmc(:,:)        ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: ciwpmc(:,:)        ! cloud ice water path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: clwpmc(:,:)        ! cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: cswpmc(:,:)        ! cloud snow path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: relqmc(:)          ! liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: reicmc(:)          ! ice particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: resnmc(:)          ! snow particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of reicmc depends on setting of iceflag:
                                                      ! iceflag = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec must be >= 10.0 microns
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]

! ------- Output -------

      integer(kind=im), intent(out)   :: ncbands      ! number of cloud spectral bands
      real(kind=rb),    intent(inout) :: taucmc(:,:)  ! cloud optical depth [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      character(len=*), intent(inout) :: errmsg
      integer,          intent(inout) :: errflg

! ------- Local -------

      integer(kind=im) :: lay                         ! Layer index
      integer(kind=im) :: ib                          ! spectral band index
      integer(kind=im) :: ig                          ! g-point interval index
      integer(kind=im) :: index                                                                                     
      integer(kind=im) :: icb(nbands)                                                                               
      real(kind=rb) , dimension(2) :: absice0
      real(kind=rb) , dimension(2,5) :: absice1
      real(kind=rb) , dimension(43,16) :: absice2
      real(kind=rb) , dimension(46,16) :: absice3
      real(kind=rb) :: absliq0
      real(kind=rb) , dimension(58,16) :: absliq1
                                                                                                                    
      real(kind=rb) :: abscoice(ngptlw)               ! ice absorption coefficients                                 
      real(kind=rb) :: abscoliq(ngptlw)               ! liquid absorption coefficients                              
      real(kind=rb) :: abscosno(ngptlw)               ! snow absorption coefficients                                
      real(kind=rb) :: cwp                            ! cloud water path                                            
      real(kind=rb) :: radice                         ! cloud ice effective size (microns)                          
      real(kind=rb) :: factor                         !                                                             
      real(kind=rb) :: fint                           !                                                             
      real(kind=rb) :: radliq                         ! cloud liquid droplet radius (microns)                       
      real(kind=rb) :: radsno                         ! cloud snow effective size (microns)                         
      real(kind=rb), parameter :: eps = 1.e-6_rb      ! epsilon                                                     
      real(kind=rb), parameter :: cldmin = 1.e-20_rb  ! minimum value for cloud quantities                          
                                                                                                                    
! ------- Definitions -------                                                                                       
                                                                                                                    
!     Explanation of the method for each value of INFLAG.  Values of                                                
!     0 or 1 for INFLAG do not distingish being liquid and ice clouds.                                              
!     INFLAG = 2 does distinguish between liquid and ice clouds, and                                                
!     requires further user input to specify the method to be used to                                               
!     compute the aborption due to each.                                                                            
!     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)                                             
!                  optical depth are input.                                                                         
!     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud                                              
!                  water path (g/m2) are input.  The (gray) cloud optical                                           
!                  depth is computed as in CCM2.                                                                    
!     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud                                                 
!                  water path (g/m2), and cloud ice fraction are input.                                             
!       ICEFLAG = 0:  The ice effective radius (microns) is input and the                                           
!                     optical depths due to ice clouds are computed as in CCM3.                                     
!       ICEFLAG = 1:  The ice effective radius (microns) is input and the                                           
!                     optical depths due to ice clouds are computed as in                                           
!                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The                                              
!                     spectral regions in this work have been matched with                                          
!                     the spectral bands in RRTM to as great an extent                                              
!                     as possible:                                                                                  
!                     E&C 1      IB = 5      RRTM bands 9-16                                                        
!                     E&C 2      IB = 4      RRTM bands 6-8                                                         
!                     E&C 3      IB = 3      RRTM bands 3-5                                                         
!                     E&C 4      IB = 2      RRTM band 2                                                            
!                     E&C 5      IB = 1      RRTM band 1                                                            
!       ICEFLAG = 2:  The ice effective radius (microns) is input and the
!                     optical properties due to ice clouds are computed from
!                     the optical properties stored in the RT code,
!                     STREAMER v3.0 (Reference: Key. J., Streamer 
!                     User's Guide, Cooperative Institute for
!                     Meteorological Satellite Studies, 2001, 96 pp.).
!                     Valid range of values for re are between 5.0 and
!                     131.0 micron.
!       ICEFLAG = 3: The ice generalized effective size (dge) is input
!                    and the optical properties, are calculated as in
!                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
!                    tables which were appropriately averaged for the
!                    bands in RRTM_LW.  Linear interpolation is used to
!                    get the coefficients from the stored tables.
!                    Valid range of values for dge are between 5.0 and
!                    140.0 micron.
!       LIQFLAG = 0:  The optical depths due to water clouds are computed as
!                     in CCM3.
!       LIQFLAG = 1:  The water droplet effective radius (microns) is input 
!                     and the optical depths due to water clouds are computed 
!                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!                     The values for absorption coefficients appropriate for
!                     the spectral bands in RRTM have been obtained for a 
!                     range of effective radii by an averaging procedure 
!                     based on the work of J. Pinto (private communication).
!                     Linear interpolation is used to get the absorption 
!                     coefficients for the input effective radius.

      data icb /1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5/
! Everything below is for INFLAG = 2.

! ABSICEn(J,IB) are the parameters needed to compute the liquid water 
! absorption coefficient in spectral region IB for ICEFLAG=n.  The units
! of ABSICEn(1,IB) are m2/g and ABSICEn(2,IB) has units (microns (m2/g)).
! For ICEFLAG = 0.

      absice0(:)= (/0.005_rb,  1.0_rb/)

! For ICEFLAG = 1.
      absice1(1,:) = (/0.0036_rb, 0.0068_rb, 0.0003_rb, 0.0016_rb,      &
     &                0.0020_rb/)
      absice1(2,:) = (/1.136_rb , 0.600_rb , 1.338_rb , 1.166_rb ,      &
     &                1.118_rb /)

! For ICEFLAG = 2.  In each band, the absorption
! coefficients are listed for a range of effective radii from 5.0
! to 131.0 microns in increments of 3.0 microns.
! Spherical Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice2(:,1) = (/ &
! band 1
       7.798999e-02_rb,6.340479e-02_rb,5.417973e-02_rb,4.766245e-02_rb,4.272663e-02_rb, &
       3.880939e-02_rb,3.559544e-02_rb,3.289241e-02_rb,3.057511e-02_rb,2.855800e-02_rb, &
       2.678022e-02_rb,2.519712e-02_rb,2.377505e-02_rb,2.248806e-02_rb,2.131578e-02_rb, &
       2.024194e-02_rb,1.925337e-02_rb,1.833926e-02_rb,1.749067e-02_rb,1.670007e-02_rb, &
       1.596113e-02_rb,1.526845e-02_rb,1.461739e-02_rb,1.400394e-02_rb,1.342462e-02_rb, &
       1.287639e-02_rb,1.235656e-02_rb,1.186279e-02_rb,1.139297e-02_rb,1.094524e-02_rb, &
       1.051794e-02_rb,1.010956e-02_rb,9.718755e-03_rb,9.344316e-03_rb,8.985139e-03_rb, &
       8.640223e-03_rb,8.308656e-03_rb,7.989606e-03_rb,7.682312e-03_rb,7.386076e-03_rb, &
       7.100255e-03_rb,6.824258e-03_rb,6.557540e-03_rb/)
      absice2(:,2) = (/ &
! band 2
       2.784879e-02_rb,2.709863e-02_rb,2.619165e-02_rb,2.529230e-02_rb,2.443225e-02_rb, &
       2.361575e-02_rb,2.284021e-02_rb,2.210150e-02_rb,2.139548e-02_rb,2.071840e-02_rb, &
       2.006702e-02_rb,1.943856e-02_rb,1.883064e-02_rb,1.824120e-02_rb,1.766849e-02_rb, &
       1.711099e-02_rb,1.656737e-02_rb,1.603647e-02_rb,1.551727e-02_rb,1.500886e-02_rb, &
       1.451045e-02_rb,1.402132e-02_rb,1.354084e-02_rb,1.306842e-02_rb,1.260355e-02_rb, &
       1.214575e-02_rb,1.169460e-02_rb,1.124971e-02_rb,1.081072e-02_rb,1.037731e-02_rb, &
       9.949167e-03_rb,9.526021e-03_rb,9.107615e-03_rb,8.693714e-03_rb,8.284096e-03_rb, &
       7.878558e-03_rb,7.476910e-03_rb,7.078974e-03_rb,6.684586e-03_rb,6.293589e-03_rb, &
       5.905839e-03_rb,5.521200e-03_rb,5.139543e-03_rb/)
      absice2(:,3) = (/ &
! band 3
       1.065397e-01_rb,8.005726e-02_rb,6.546428e-02_rb,5.589131e-02_rb,4.898681e-02_rb, &
       4.369932e-02_rb,3.947901e-02_rb,3.600676e-02_rb,3.308299e-02_rb,3.057561e-02_rb, &
       2.839325e-02_rb,2.647040e-02_rb,2.475872e-02_rb,2.322164e-02_rb,2.183091e-02_rb, &
       2.056430e-02_rb,1.940407e-02_rb,1.833586e-02_rb,1.734787e-02_rb,1.643034e-02_rb, &
       1.557512e-02_rb,1.477530e-02_rb,1.402501e-02_rb,1.331924e-02_rb,1.265364e-02_rb, &
       1.202445e-02_rb,1.142838e-02_rb,1.086257e-02_rb,1.032445e-02_rb,9.811791e-03_rb, &
       9.322587e-03_rb,8.855053e-03_rb,8.407591e-03_rb,7.978763e-03_rb,7.567273e-03_rb, &
       7.171949e-03_rb,6.791728e-03_rb,6.425642e-03_rb,6.072809e-03_rb,5.732424e-03_rb, &
       5.403748e-03_rb,5.086103e-03_rb,4.778865e-03_rb/)
      absice2(:,4) = (/ &
! band 4
       1.804566e-01_rb,1.168987e-01_rb,8.680442e-02_rb,6.910060e-02_rb,5.738174e-02_rb, &
       4.902332e-02_rb,4.274585e-02_rb,3.784923e-02_rb,3.391734e-02_rb,3.068690e-02_rb, &
       2.798301e-02_rb,2.568480e-02_rb,2.370600e-02_rb,2.198337e-02_rb,2.046940e-02_rb, &   
       1.912777e-02_rb,1.793016e-02_rb,1.685420e-02_rb,1.588193e-02_rb,1.499882e-02_rb, &   
       1.419293e-02_rb,1.345440e-02_rb,1.277496e-02_rb,1.214769e-02_rb,1.156669e-02_rb, &   
       1.102694e-02_rb,1.052412e-02_rb,1.005451e-02_rb,9.614854e-03_rb,9.202335e-03_rb, &   
       8.814470e-03_rb,8.449077e-03_rb,8.104223e-03_rb,7.778195e-03_rb,7.469466e-03_rb, &   
       7.176671e-03_rb,6.898588e-03_rb,6.634117e-03_rb,6.382264e-03_rb,6.142134e-03_rb, &   
       5.912913e-03_rb,5.693862e-03_rb,5.484308e-03_rb/)                                    
      absice2(:,5) = (/ &                                                                   
! band 5                                                                                    
       2.131806e-01_rb,1.311372e-01_rb,9.407171e-02_rb,7.299442e-02_rb,5.941273e-02_rb, &   
       4.994043e-02_rb,4.296242e-02_rb,3.761113e-02_rb,3.337910e-02_rb,2.994978e-02_rb, &   
       2.711556e-02_rb,2.473461e-02_rb,2.270681e-02_rb,2.095943e-02_rb,1.943839e-02_rb, &   
       1.810267e-02_rb,1.692057e-02_rb,1.586719e-02_rb,1.492275e-02_rb,1.407132e-02_rb, &   
       1.329989e-02_rb,1.259780e-02_rb,1.195618e-02_rb,1.136761e-02_rb,1.082583e-02_rb, &   
       1.032552e-02_rb,9.862158e-03_rb,9.431827e-03_rb,9.031157e-03_rb,8.657217e-03_rb, &   
       8.307449e-03_rb,7.979609e-03_rb,7.671724e-03_rb,7.382048e-03_rb,7.109032e-03_rb, &   
       6.851298e-03_rb,6.607615e-03_rb,6.376881e-03_rb,6.158105e-03_rb,5.950394e-03_rb, &   
       5.752942e-03_rb,5.565019e-03_rb,5.385963e-03_rb/)                                    
      absice2(:,6) = (/ &                                                                   
! band 6                                                                                    
       1.546177e-01_rb,1.039251e-01_rb,7.910347e-02_rb,6.412429e-02_rb,5.399997e-02_rb, &   
       4.664937e-02_rb,4.104237e-02_rb,3.660781e-02_rb,3.300218e-02_rb,3.000586e-02_rb, &   
       2.747148e-02_rb,2.529633e-02_rb,2.340647e-02_rb,2.174723e-02_rb,2.027731e-02_rb, &   
       1.896487e-02_rb,1.778492e-02_rb,1.671761e-02_rb,1.574692e-02_rb,1.485978e-02_rb, &   
       1.404543e-02_rb,1.329489e-02_rb,1.260066e-02_rb,1.195636e-02_rb,1.135657e-02_rb, &   
       1.079664e-02_rb,1.027257e-02_rb,9.780871e-03_rb,9.318505e-03_rb,8.882815e-03_rb, &   
       8.471458e-03_rb,8.082364e-03_rb,7.713696e-03_rb,7.363817e-03_rb,7.031264e-03_rb, &   
       6.714725e-03_rb,6.413021e-03_rb,6.125086e-03_rb,5.849958e-03_rb,5.586764e-03_rb, &   
       5.334707e-03_rb,5.093066e-03_rb,4.861179e-03_rb/)                                    
      absice2(:,7) = (/ &                                                                   
! band 7                                                                                    
       7.583404e-02_rb,6.181558e-02_rb,5.312027e-02_rb,4.696039e-02_rb,4.225986e-02_rb, &   
       3.849735e-02_rb,3.538340e-02_rb,3.274182e-02_rb,3.045798e-02_rb,2.845343e-02_rb, &   
       2.667231e-02_rb,2.507353e-02_rb,2.362606e-02_rb,2.230595e-02_rb,2.109435e-02_rb, &   
       1.997617e-02_rb,1.893916e-02_rb,1.797328e-02_rb,1.707016e-02_rb,1.622279e-02_rb, &   
       1.542523e-02_rb,1.467241e-02_rb,1.395997e-02_rb,1.328414e-02_rb,1.264164e-02_rb, &   
       1.202958e-02_rb,1.144544e-02_rb,1.088697e-02_rb,1.035218e-02_rb,9.839297e-03_rb, &   
       9.346733e-03_rb,8.873057e-03_rb,8.416980e-03_rb,7.977335e-03_rb,7.553066e-03_rb, &   
       7.143210e-03_rb,6.746888e-03_rb,6.363297e-03_rb,5.991700e-03_rb,5.631422e-03_rb, &   
       5.281840e-03_rb,4.942378e-03_rb,4.612505e-03_rb/)                                    
      absice2(:,8) = (/ &                                                                   
! band 8                                                                                    
       9.022185e-02_rb,6.922700e-02_rb,5.710674e-02_rb,4.898377e-02_rb,4.305946e-02_rb, &   
       3.849553e-02_rb,3.484183e-02_rb,3.183220e-02_rb,2.929794e-02_rb,2.712627e-02_rb, &   
       2.523856e-02_rb,2.357810e-02_rb,2.210286e-02_rb,2.078089e-02_rb,1.958747e-02_rb, &   
       1.850310e-02_rb,1.751218e-02_rb,1.660205e-02_rb,1.576232e-02_rb,1.498440e-02_rb, &   
       1.426107e-02_rb,1.358624e-02_rb,1.295474e-02_rb,1.236212e-02_rb,1.180456e-02_rb, &   
       1.127874e-02_rb,1.078175e-02_rb,1.031106e-02_rb,9.864433e-03_rb,9.439878e-03_rb, &   
       9.035637e-03_rb,8.650140e-03_rb,8.281981e-03_rb,7.929895e-03_rb,7.592746e-03_rb, &   
       7.269505e-03_rb,6.959238e-03_rb,6.661100e-03_rb,6.374317e-03_rb,6.098185e-03_rb, &   
       5.832059e-03_rb,5.575347e-03_rb,5.327504e-03_rb/)                                    
      absice2(:,9) = (/ &
! band 9
       1.294087e-01_rb,8.788217e-02_rb,6.728288e-02_rb,5.479720e-02_rb,4.635049e-02_rb, &
       4.022253e-02_rb,3.555576e-02_rb,3.187259e-02_rb,2.888498e-02_rb,2.640843e-02_rb, &
       2.431904e-02_rb,2.253038e-02_rb,2.098024e-02_rb,1.962267e-02_rb,1.842293e-02_rb, &
       1.735426e-02_rb,1.639571e-02_rb,1.553060e-02_rb,1.474552e-02_rb,1.402953e-02_rb, &
       1.337363e-02_rb,1.277033e-02_rb,1.221336e-02_rb,1.169741e-02_rb,1.121797e-02_rb, &
       1.077117e-02_rb,1.035369e-02_rb,9.962643e-03_rb,9.595509e-03_rb,9.250088e-03_rb, &
       8.924447e-03_rb,8.616876e-03_rb,8.325862e-03_rb,8.050057e-03_rb,7.788258e-03_rb, &
       7.539388e-03_rb,7.302478e-03_rb,7.076656e-03_rb,6.861134e-03_rb,6.655197e-03_rb, &
       6.458197e-03_rb,6.269543e-03_rb,6.088697e-03_rb/)
      absice2(:,10) = (/ &
! band 10
       1.593628e-01_rb,1.014552e-01_rb,7.458955e-02_rb,5.903571e-02_rb,4.887582e-02_rb, &
       4.171159e-02_rb,3.638480e-02_rb,3.226692e-02_rb,2.898717e-02_rb,2.631256e-02_rb, &
       2.408925e-02_rb,2.221156e-02_rb,2.060448e-02_rb,1.921325e-02_rb,1.799699e-02_rb, &
       1.692456e-02_rb,1.597177e-02_rb,1.511961e-02_rb,1.435289e-02_rb,1.365933e-02_rb, &
       1.302890e-02_rb,1.245334e-02_rb,1.192576e-02_rb,1.144037e-02_rb,1.099230e-02_rb, &   
       1.057739e-02_rb,1.019208e-02_rb,9.833302e-03_rb,9.498395e-03_rb,9.185047e-03_rb, &   
       8.891237e-03_rb,8.615185e-03_rb,8.355325e-03_rb,8.110267e-03_rb,7.878778e-03_rb, &   
       7.659759e-03_rb,7.452224e-03_rb,7.255291e-03_rb,7.068166e-03_rb,6.890130e-03_rb, &   
       6.720536e-03_rb,6.558794e-03_rb,6.404371e-03_rb/)                                    
      absice2(:,11) = (/ &                                                                  
! band 11                                                                                   
       1.656227e-01_rb,1.032129e-01_rb,7.487359e-02_rb,5.871431e-02_rb,4.828355e-02_rb, &   
       4.099989e-02_rb,3.562924e-02_rb,3.150755e-02_rb,2.824593e-02_rb,2.560156e-02_rb, &   
       2.341503e-02_rb,2.157740e-02_rb,2.001169e-02_rb,1.866199e-02_rb,1.748669e-02_rb, &   
       1.645421e-02_rb,1.554015e-02_rb,1.472535e-02_rb,1.399457e-02_rb,1.333553e-02_rb, &   
       1.273821e-02_rb,1.219440e-02_rb,1.169725e-02_rb,1.124104e-02_rb,1.082096e-02_rb, &   
       1.043290e-02_rb,1.007336e-02_rb,9.739338e-03_rb,9.428223e-03_rb,9.137756e-03_rb, &   
       8.865964e-03_rb,8.611115e-03_rb,8.371686e-03_rb,8.146330e-03_rb,7.933852e-03_rb, &   
       7.733187e-03_rb,7.543386e-03_rb,7.363597e-03_rb,7.193056e-03_rb,7.031072e-03_rb, &   
       6.877024e-03_rb,6.730348e-03_rb,6.590531e-03_rb/)                                    
      absice2(:,12) = (/ &                                                                  
! band 12                                                                                   
       9.194591e-02_rb,6.446867e-02_rb,4.962034e-02_rb,4.042061e-02_rb,3.418456e-02_rb, &   
       2.968856e-02_rb,2.629900e-02_rb,2.365572e-02_rb,2.153915e-02_rb,1.980791e-02_rb, &   
       1.836689e-02_rb,1.714979e-02_rb,1.610900e-02_rb,1.520946e-02_rb,1.442476e-02_rb, &   
       1.373468e-02_rb,1.312345e-02_rb,1.257858e-02_rb,1.209010e-02_rb,1.164990e-02_rb, &   
       1.125136e-02_rb,1.088901e-02_rb,1.055827e-02_rb,1.025531e-02_rb,9.976896e-03_rb, &   
       9.720255e-03_rb,9.483022e-03_rb,9.263160e-03_rb,9.058902e-03_rb,8.868710e-03_rb, &   
       8.691240e-03_rb,8.525312e-03_rb,8.369886e-03_rb,8.224042e-03_rb,8.086961e-03_rb, &   
       7.957917e-03_rb,7.836258e-03_rb,7.721400e-03_rb,7.612821e-03_rb,7.510045e-03_rb, &   
       7.412648e-03_rb,7.320242e-03_rb,7.232476e-03_rb/)           
      absice2(:,13) = (/ &
! band 13                                                                                   
       1.437021e-01_rb,8.872535e-02_rb,6.392420e-02_rb,4.991833e-02_rb,4.096790e-02_rb, &   
       3.477881e-02_rb,3.025782e-02_rb,2.681909e-02_rb,2.412102e-02_rb,2.195132e-02_rb, &   
       2.017124e-02_rb,1.868641e-02_rb,1.743044e-02_rb,1.635529e-02_rb,1.542540e-02_rb, &   
       1.461388e-02_rb,1.390003e-02_rb,1.326766e-02_rb,1.270395e-02_rb,1.219860e-02_rb, &   
       1.174326e-02_rb,1.133107e-02_rb,1.095637e-02_rb,1.061442e-02_rb,1.030126e-02_rb, &   
       1.001352e-02_rb,9.748340e-03_rb,9.503256e-03_rb,9.276155e-03_rb,9.065205e-03_rb, &   
       8.868808e-03_rb,8.685571e-03_rb,8.514268e-03_rb,8.353820e-03_rb,8.203272e-03_rb, &   
       8.061776e-03_rb,7.928578e-03_rb,7.803001e-03_rb,7.684443e-03_rb,7.572358e-03_rb, &   
       7.466258e-03_rb,7.365701e-03_rb,7.270286e-03_rb/)                                    
      absice2(:,14) = (/ &                                                                  
! band 14                                                                                   
       1.288870e-01_rb,8.160295e-02_rb,5.964745e-02_rb,4.703790e-02_rb,3.888637e-02_rb, &   
       3.320115e-02_rb,2.902017e-02_rb,2.582259e-02_rb,2.330224e-02_rb,2.126754e-02_rb, &   
       1.959258e-02_rb,1.819130e-02_rb,1.700289e-02_rb,1.598320e-02_rb,1.509942e-02_rb, &   
       1.432666e-02_rb,1.364572e-02_rb,1.304156e-02_rb,1.250220e-02_rb,1.201803e-02_rb, &   
       1.158123e-02_rb,1.118537e-02_rb,1.082513e-02_rb,1.049605e-02_rb,1.019440e-02_rb, &   
       9.916989e-03_rb,9.661116e-03_rb,9.424457e-03_rb,9.205005e-03_rb,9.001022e-03_rb, &   
       8.810992e-03_rb,8.633588e-03_rb,8.467646e-03_rb,8.312137e-03_rb,8.166151e-03_rb, &   
       8.028878e-03_rb,7.899597e-03_rb,7.777663e-03_rb,7.662498e-03_rb,7.553581e-03_rb, &   
       7.450444e-03_rb,7.352662e-03_rb,7.259851e-03_rb/)                                    
      absice2(:,15) = (/ &                                                                  
! band 15                                                                                   
       8.254229e-02_rb,5.808787e-02_rb,4.492166e-02_rb,3.675028e-02_rb,3.119623e-02_rb, &   
       2.718045e-02_rb,2.414450e-02_rb,2.177073e-02_rb,1.986526e-02_rb,1.830306e-02_rb, &   
       1.699991e-02_rb,1.589698e-02_rb,1.495199e-02_rb,1.413374e-02_rb,1.341870e-02_rb, &   
       1.278883e-02_rb,1.223002e-02_rb,1.173114e-02_rb,1.128322e-02_rb,1.087900e-02_rb, &   
       1.051254e-02_rb,1.017890e-02_rb,9.873991e-03_rb,9.594347e-03_rb,9.337044e-03_rb, &   
       9.099589e-03_rb,8.879842e-03_rb,8.675960e-03_rb,8.486341e-03_rb,8.309594e-03_rb, &   
       8.144500e-03_rb,7.989986e-03_rb,7.845109e-03_rb,7.709031e-03_rb,7.581007e-03_rb, &   
       7.460376e-03_rb,7.346544e-03_rb,7.238978e-03_rb,7.137201e-03_rb,7.040780e-03_rb, &   
       6.949325e-03_rb,6.862483e-03_rb,6.779931e-03_rb/)                                    
      absice2(:,16) = (/ &                                                                  
! band 16                                                                                   
       1.382062e-01_rb,8.643227e-02_rb,6.282935e-02_rb,4.934783e-02_rb,4.063891e-02_rb, &   
       3.455591e-02_rb,3.007059e-02_rb,2.662897e-02_rb,2.390631e-02_rb,2.169972e-02_rb, &   
       1.987596e-02_rb,1.834393e-02_rb,1.703924e-02_rb,1.591513e-02_rb,1.493679e-02_rb, &   
       1.407780e-02_rb,1.331775e-02_rb,1.264061e-02_rb,1.203364e-02_rb,1.148655e-02_rb, &   
       1.099099e-02_rb,1.054006e-02_rb,1.012807e-02_rb,9.750215e-03_rb,9.402477e-03_rb, &   
       9.081428e-03_rb,8.784143e-03_rb,8.508107e-03_rb,8.251146e-03_rb,8.011373e-03_rb, &   
       7.787140e-03_rb,7.577002e-03_rb,7.379687e-03_rb,7.194071e-03_rb,7.019158e-03_rb, &   
       6.854061e-03_rb,6.697986e-03_rb,6.550224e-03_rb,6.410138e-03_rb,6.277153e-03_rb, &   
       6.150751e-03_rb,6.030462e-03_rb,5.915860e-03_rb/)                                    
                                                                                           
! ICEFLAG = 3; Fu parameterization. Particle size 5 - 140 micron in 
! increments of 3 microns.
! units = m2/g
! Hexagonal Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice3(:,1) = (/ &
! band 1
       3.110649e-03_rb,4.666352e-02_rb,6.606447e-02_rb,6.531678e-02_rb,6.012598e-02_rb, &
       5.437494e-02_rb,4.906411e-02_rb,4.441146e-02_rb,4.040585e-02_rb,3.697334e-02_rb, &
       3.403027e-02_rb,3.149979e-02_rb,2.931596e-02_rb,2.742365e-02_rb,2.577721e-02_rb, &
       2.433888e-02_rb,2.307732e-02_rb,2.196644e-02_rb,2.098437e-02_rb,2.011264e-02_rb, &
       1.933561e-02_rb,1.863992e-02_rb,1.801407e-02_rb,1.744812e-02_rb,1.693346e-02_rb, &
       1.646252e-02_rb,1.602866e-02_rb,1.562600e-02_rb,1.524933e-02_rb,1.489399e-02_rb, &
       1.455580e-02_rb,1.423098e-02_rb,1.391612e-02_rb,1.360812e-02_rb,1.330413e-02_rb, &
       1.300156e-02_rb,1.269801e-02_rb,1.239127e-02_rb,1.207928e-02_rb,1.176014e-02_rb, &
       1.143204e-02_rb,1.109334e-02_rb,1.074243e-02_rb,1.037786e-02_rb,9.998198e-03_rb, &
       9.602126e-03_rb/)
      absice3(:,2) = (/ &
! band 2
       3.984966e-04_rb,1.681097e-02_rb,2.627680e-02_rb,2.767465e-02_rb,2.700722e-02_rb, &
       2.579180e-02_rb,2.448677e-02_rb,2.323890e-02_rb,2.209096e-02_rb,2.104882e-02_rb, &
       2.010547e-02_rb,1.925003e-02_rb,1.847128e-02_rb,1.775883e-02_rb,1.710358e-02_rb, &
       1.649769e-02_rb,1.593449e-02_rb,1.540829e-02_rb,1.491429e-02_rb,1.444837e-02_rb, &
       1.400704e-02_rb,1.358729e-02_rb,1.318654e-02_rb,1.280258e-02_rb,1.243346e-02_rb, &
       1.207750e-02_rb,1.173325e-02_rb,1.139941e-02_rb,1.107487e-02_rb,1.075861e-02_rb, &
       1.044975e-02_rb,1.014753e-02_rb,9.851229e-03_rb,9.560240e-03_rb,9.274003e-03_rb, &
       8.992020e-03_rb,8.713845e-03_rb,8.439074e-03_rb,8.167346e-03_rb,7.898331e-03_rb, &
       7.631734e-03_rb,7.367286e-03_rb,7.104742e-03_rb,6.843882e-03_rb,6.584504e-03_rb, &
       6.326424e-03_rb/)
      absice3(:,3) = (/ &
! band 3
       6.933163e-02_rb,8.540475e-02_rb,7.701816e-02_rb,6.771158e-02_rb,5.986953e-02_rb, &
       5.348120e-02_rb,4.824962e-02_rb,4.390563e-02_rb,4.024411e-02_rb,3.711404e-02_rb, &
       3.440426e-02_rb,3.203200e-02_rb,2.993478e-02_rb,2.806474e-02_rb,2.638464e-02_rb, &
       2.486516e-02_rb,2.348288e-02_rb,2.221890e-02_rb,2.105780e-02_rb,1.998687e-02_rb, &
       1.899552e-02_rb,1.807490e-02_rb,1.721750e-02_rb,1.641693e-02_rb,1.566773e-02_rb, &
       1.496515e-02_rb,1.430509e-02_rb,1.368398e-02_rb,1.309865e-02_rb,1.254634e-02_rb, &
       1.202456e-02_rb,1.153114e-02_rb,1.106409e-02_rb,1.062166e-02_rb,1.020224e-02_rb, &
       9.804381e-03_rb,9.426771e-03_rb,9.068205e-03_rb,8.727578e-03_rb,8.403876e-03_rb, &
       8.096160e-03_rb,7.803564e-03_rb,7.525281e-03_rb,7.260560e-03_rb,7.008697e-03_rb, &
       6.769036e-03_rb/)
      absice3(:,4) = (/ &
! band 4
       1.765735e-01_rb,1.382700e-01_rb,1.095129e-01_rb,8.987475e-02_rb,7.591185e-02_rb, &
       6.554169e-02_rb,5.755500e-02_rb,5.122083e-02_rb,4.607610e-02_rb,4.181475e-02_rb, &
       3.822697e-02_rb,3.516432e-02_rb,3.251897e-02_rb,3.021073e-02_rb,2.817876e-02_rb, &
       2.637607e-02_rb,2.476582e-02_rb,2.331871e-02_rb,2.201113e-02_rb,2.082388e-02_rb, &
       1.974115e-02_rb,1.874983e-02_rb,1.783894e-02_rb,1.699922e-02_rb,1.622280e-02_rb, &
       1.550296e-02_rb,1.483390e-02_rb,1.421064e-02_rb,1.362880e-02_rb,1.308460e-02_rb, &
       1.257468e-02_rb,1.209611e-02_rb,1.164628e-02_rb,1.122287e-02_rb,1.082381e-02_rb, &
       1.044725e-02_rb,1.009154e-02_rb,9.755166e-03_rb,9.436783e-03_rb,9.135163e-03_rb, &
       8.849193e-03_rb,8.577856e-03_rb,8.320225e-03_rb,8.075451e-03_rb,7.842755e-03_rb, &
       7.621418e-03_rb/)
      absice3(:,5) = (/ &
! band 5
       2.339673e-01_rb,1.692124e-01_rb,1.291656e-01_rb,1.033837e-01_rb,8.562949e-02_rb, &
       7.273526e-02_rb,6.298262e-02_rb,5.537015e-02_rb,4.927787e-02_rb,4.430246e-02_rb, &
       4.017061e-02_rb,3.669072e-02_rb,3.372455e-02_rb,3.116995e-02_rb,2.894977e-02_rb, &
       2.700471e-02_rb,2.528842e-02_rb,2.376420e-02_rb,2.240256e-02_rb,2.117959e-02_rb, &
       2.007567e-02_rb,1.907456e-02_rb,1.816271e-02_rb,1.732874e-02_rb,1.656300e-02_rb, &
       1.585725e-02_rb,1.520445e-02_rb,1.459852e-02_rb,1.403419e-02_rb,1.350689e-02_rb, &
       1.301260e-02_rb,1.254781e-02_rb,1.210941e-02_rb,1.169468e-02_rb,1.130118e-02_rb, &
       1.092675e-02_rb,1.056945e-02_rb,1.022757e-02_rb,9.899560e-03_rb,9.584021e-03_rb, &
       9.279705e-03_rb,8.985479e-03_rb,8.700322e-03_rb,8.423306e-03_rb,8.153590e-03_rb, &
       7.890412e-03_rb/)                                                             
      absice3(:,6) = (/ &                                                            
! band 6                                                                             
       1.145369e-01_rb,1.174566e-01_rb,9.917866e-02_rb,8.332990e-02_rb,7.104263e-02_rb, &
       6.153370e-02_rb,5.405472e-02_rb,4.806281e-02_rb,4.317918e-02_rb,3.913795e-02_rb, &
       3.574916e-02_rb,3.287437e-02_rb,3.041067e-02_rb,2.828017e-02_rb,2.642292e-02_rb, &
       2.479206e-02_rb,2.335051e-02_rb,2.206851e-02_rb,2.092195e-02_rb,1.989108e-02_rb, &
       1.895958e-02_rb,1.811385e-02_rb,1.734245e-02_rb,1.663573e-02_rb,1.598545e-02_rb, &
       1.538456e-02_rb,1.482700e-02_rb,1.430750e-02_rb,1.382150e-02_rb,1.336499e-02_rb, &
       1.293447e-02_rb,1.252685e-02_rb,1.213939e-02_rb,1.176968e-02_rb,1.141555e-02_rb, &
       1.107508e-02_rb,1.074655e-02_rb,1.042839e-02_rb,1.011923e-02_rb,9.817799e-03_rb, &
       9.522962e-03_rb,9.233688e-03_rb,8.949041e-03_rb,8.668171e-03_rb,8.390301e-03_rb, &
       8.114723e-03_rb/)                                                             
      absice3(:,7) = (/ &                                                            
! band 7                                                                             
       1.222345e-02_rb,5.344230e-02_rb,5.523465e-02_rb,5.128759e-02_rb,4.676925e-02_rb, &
       4.266150e-02_rb,3.910561e-02_rb,3.605479e-02_rb,3.342843e-02_rb,3.115052e-02_rb, &
       2.915776e-02_rb,2.739935e-02_rb,2.583499e-02_rb,2.443266e-02_rb,2.316681e-02_rb, &
       2.201687e-02_rb,2.096619e-02_rb,2.000112e-02_rb,1.911044e-02_rb,1.828481e-02_rb, &
       1.751641e-02_rb,1.679866e-02_rb,1.612598e-02_rb,1.549360e-02_rb,1.489742e-02_rb, &
       1.433392e-02_rb,1.380002e-02_rb,1.329305e-02_rb,1.281068e-02_rb,1.235084e-02_rb, &
       1.191172e-02_rb,1.149171e-02_rb,1.108936e-02_rb,1.070341e-02_rb,1.033271e-02_rb, &
       9.976220e-03_rb,9.633021e-03_rb,9.302273e-03_rb,8.983216e-03_rb,8.675161e-03_rb, &
       8.377478e-03_rb,8.089595e-03_rb,7.810986e-03_rb,7.541170e-03_rb,7.279706e-03_rb, &
       7.026186e-03_rb/)                                                             
      absice3(:,8) = (/ &                                                            
! band 8                                                                             
       6.711058e-02_rb,6.918198e-02_rb,6.127484e-02_rb,5.411944e-02_rb,4.836902e-02_rb, &
       4.375293e-02_rb,3.998077e-02_rb,3.683587e-02_rb,3.416508e-02_rb,3.186003e-02_rb, &
       2.984290e-02_rb,2.805671e-02_rb,2.645895e-02_rb,2.501733e-02_rb,2.370689e-02_rb, &
       2.250808e-02_rb,2.140532e-02_rb,2.038609e-02_rb,1.944018e-02_rb,1.855918e-02_rb, &
       1.773609e-02_rb,1.696504e-02_rb,1.624106e-02_rb,1.555990e-02_rb,1.491793e-02_rb, &
       1.431197e-02_rb,1.373928e-02_rb,1.319743e-02_rb,1.268430e-02_rb,1.219799e-02_rb, &
       1.173682e-02_rb,1.129925e-02_rb,1.088393e-02_rb,1.048961e-02_rb,1.011516e-02_rb, &
       9.759543e-03_rb,9.421813e-03_rb,9.101089e-03_rb,8.796559e-03_rb,8.507464e-03_rb, &
       8.233098e-03_rb,7.972798e-03_rb,7.725942e-03_rb,7.491940e-03_rb,7.270238e-03_rb, &
       7.060305e-03_rb/)                                                             
      absice3(:,9) = (/ &                                                            
! band 9                                                                             
       1.236780e-01_rb,9.222386e-02_rb,7.383997e-02_rb,6.204072e-02_rb,5.381029e-02_rb, &
       4.770678e-02_rb,4.296928e-02_rb,3.916131e-02_rb,3.601540e-02_rb,3.335878e-02_rb, &
       3.107493e-02_rb,2.908247e-02_rb,2.732282e-02_rb,2.575276e-02_rb,2.433968e-02_rb, &
       2.305852e-02_rb,2.188966e-02_rb,2.081757e-02_rb,1.982974e-02_rb,1.891599e-02_rb, &
       1.806794e-02_rb,1.727865e-02_rb,1.654227e-02_rb,1.585387e-02_rb,1.520924e-02_rb, &
       1.460476e-02_rb,1.403730e-02_rb,1.350416e-02_rb,1.300293e-02_rb,1.253153e-02_rb, &
       1.208808e-02_rb,1.167094e-02_rb,1.127862e-02_rb,1.090979e-02_rb,1.056323e-02_rb, &
       1.023786e-02_rb,9.932665e-03_rb,9.646744e-03_rb,9.379250e-03_rb,9.129409e-03_rb, &
       8.896500e-03_rb,8.679856e-03_rb,8.478852e-03_rb,8.292904e-03_rb,8.121463e-03_rb, &
       7.964013e-03_rb/)                                                             
      absice3(:,10) = (/ &                                                           
! band 10                                                                            
       1.655966e-01_rb,1.134205e-01_rb,8.714344e-02_rb,7.129241e-02_rb,6.063739e-02_rb, &
       5.294203e-02_rb,4.709309e-02_rb,4.247476e-02_rb,3.871892e-02_rb,3.559206e-02_rb, &
       3.293893e-02_rb,3.065226e-02_rb,2.865558e-02_rb,2.689288e-02_rb,2.532221e-02_rb, &
       2.391150e-02_rb,2.263582e-02_rb,2.147549e-02_rb,2.041476e-02_rb,1.944089e-02_rb, &
       1.854342e-02_rb,1.771371e-02_rb,1.694456e-02_rb,1.622989e-02_rb,1.556456e-02_rb, &
       1.494415e-02_rb,1.436491e-02_rb,1.382354e-02_rb,1.331719e-02_rb,1.284339e-02_rb, &
       1.239992e-02_rb,1.198486e-02_rb,1.159647e-02_rb,1.123323e-02_rb,1.089375e-02_rb, &
       1.057679e-02_rb,1.028124e-02_rb,1.000607e-02_rb,9.750376e-03_rb,9.513303e-03_rb, &
       9.294082e-03_rb,9.092003e-03_rb,8.906412e-03_rb,8.736702e-03_rb,8.582314e-03_rb, &
       8.442725e-03_rb/)                                                             
      absice3(:,11) = (/ &                                                           
! band 11                                                                            
       1.775615e-01_rb,1.180046e-01_rb,8.929607e-02_rb,7.233500e-02_rb,6.108333e-02_rb, &
       5.303642e-02_rb,4.696927e-02_rb,4.221206e-02_rb,3.836768e-02_rb,3.518576e-02_rb, &
       3.250063e-02_rb,3.019825e-02_rb,2.819758e-02_rb,2.643943e-02_rb,2.487953e-02_rb, &
       2.348414e-02_rb,2.222705e-02_rb,2.108762e-02_rb,2.004936e-02_rb,1.909892e-02_rb, &
       1.822539e-02_rb,1.741975e-02_rb,1.667449e-02_rb,1.598330e-02_rb,1.534084e-02_rb, &
       1.474253e-02_rb,1.418446e-02_rb,1.366325e-02_rb,1.317597e-02_rb,1.272004e-02_rb, &
       1.229321e-02_rb,1.189350e-02_rb,1.151915e-02_rb,1.116859e-02_rb,1.084042e-02_rb, &
       1.053338e-02_rb,1.024636e-02_rb,9.978326e-03_rb,9.728357e-03_rb,9.495613e-03_rb, &
       9.279327e-03_rb,9.078798e-03_rb,8.893383e-03_rb,8.722488e-03_rb,8.565568e-03_rb, &
       8.422115e-03_rb/)                                                             
      absice3(:,12) = (/ &                                                           
! band 12                                                                            
       9.465447e-02_rb,6.432047e-02_rb,5.060973e-02_rb,4.267283e-02_rb,3.741843e-02_rb, &
       3.363096e-02_rb,3.073531e-02_rb,2.842405e-02_rb,2.651789e-02_rb,2.490518e-02_rb, &
       2.351273e-02_rb,2.229056e-02_rb,2.120335e-02_rb,2.022541e-02_rb,1.933763e-02_rb, &
       1.852546e-02_rb,1.777763e-02_rb,1.708528e-02_rb,1.644134e-02_rb,1.584009e-02_rb, &
       1.527684e-02_rb,1.474774e-02_rb,1.424955e-02_rb,1.377957e-02_rb,1.333549e-02_rb, &
       1.291534e-02_rb,1.251743e-02_rb,1.214029e-02_rb,1.178265e-02_rb,1.144337e-02_rb, &
       1.112148e-02_rb,1.081609e-02_rb,1.052642e-02_rb,1.025178e-02_rb,9.991540e-03_rb, &
       9.745130e-03_rb,9.512038e-03_rb,9.291797e-03_rb,9.083980e-03_rb,8.888195e-03_rb, &
       8.704081e-03_rb,8.531306e-03_rb,8.369560e-03_rb,8.218558e-03_rb,8.078032e-03_rb, &
       7.947730e-03_rb/)                                                             
      absice3(:,13) = (/ &                                                           
! band 13                                                                            
       1.560311e-01_rb,9.961097e-02_rb,7.502949e-02_rb,6.115022e-02_rb,5.214952e-02_rb, &
       4.578149e-02_rb,4.099731e-02_rb,3.724174e-02_rb,3.419343e-02_rb,3.165356e-02_rb, &
       2.949251e-02_rb,2.762222e-02_rb,2.598073e-02_rb,2.452322e-02_rb,2.321642e-02_rb, &
       2.203516e-02_rb,2.096002e-02_rb,1.997579e-02_rb,1.907036e-02_rb,1.823401e-02_rb, &
       1.745879e-02_rb,1.673819e-02_rb,1.606678e-02_rb,1.544003e-02_rb,1.485411e-02_rb, &
       1.430574e-02_rb,1.379215e-02_rb,1.331092e-02_rb,1.285996e-02_rb,1.243746e-02_rb, &
       1.204183e-02_rb,1.167164e-02_rb,1.132567e-02_rb,1.100281e-02_rb,1.070207e-02_rb, &
       1.042258e-02_rb,1.016352e-02_rb,9.924197e-03_rb,9.703953e-03_rb,9.502199e-03_rb, &
       9.318400e-03_rb,9.152066e-03_rb,9.002749e-03_rb,8.870038e-03_rb,8.753555e-03_rb, &
       8.652951e-03_rb/)                                                             
      absice3(:,14) = (/ &                                                           
! band 14                                                                            
       1.559547e-01_rb,9.896700e-02_rb,7.441231e-02_rb,6.061469e-02_rb,5.168730e-02_rb, &
       4.537821e-02_rb,4.064106e-02_rb,3.692367e-02_rb,3.390714e-02_rb,3.139438e-02_rb, &
       2.925702e-02_rb,2.740783e-02_rb,2.578547e-02_rb,2.434552e-02_rb,2.305506e-02_rb, &
       2.188910e-02_rb,2.082842e-02_rb,1.985789e-02_rb,1.896553e-02_rb,1.814165e-02_rb, &
       1.737839e-02_rb,1.666927e-02_rb,1.600891e-02_rb,1.539279e-02_rb,1.481712e-02_rb, &
       1.427865e-02_rb,1.377463e-02_rb,1.330266e-02_rb,1.286068e-02_rb,1.244689e-02_rb, &
       1.205973e-02_rb,1.169780e-02_rb,1.135989e-02_rb,1.104492e-02_rb,1.075192e-02_rb, &
       1.048004e-02_rb,1.022850e-02_rb,9.996611e-03_rb,9.783753e-03_rb,9.589361e-03_rb, &
       9.412924e-03_rb,9.253977e-03_rb,9.112098e-03_rb,8.986903e-03_rb,8.878039e-03_rb, &
       8.785184e-03_rb/)                                                             
      absice3(:,15) = (/ &                                                           
! band 15                                                                            
       1.102926e-01_rb,7.176622e-02_rb,5.530316e-02_rb,4.606056e-02_rb,4.006116e-02_rb, &
       3.579628e-02_rb,3.256909e-02_rb,3.001360e-02_rb,2.791920e-02_rb,2.615617e-02_rb, &
       2.464023e-02_rb,2.331426e-02_rb,2.213817e-02_rb,2.108301e-02_rb,2.012733e-02_rb, &
       1.925493e-02_rb,1.845331e-02_rb,1.771269e-02_rb,1.702531e-02_rb,1.638493e-02_rb, &
       1.578648e-02_rb,1.522579e-02_rb,1.469940e-02_rb,1.420442e-02_rb,1.373841e-02_rb, &
       1.329931e-02_rb,1.288535e-02_rb,1.249502e-02_rb,1.212700e-02_rb,1.178015e-02_rb, &
       1.145348e-02_rb,1.114612e-02_rb,1.085730e-02_rb,1.058633e-02_rb,1.033263e-02_rb, &
       1.009564e-02_rb,9.874895e-03_rb,9.669960e-03_rb,9.480449e-03_rb,9.306014e-03_rb, &
       9.146339e-03_rb,9.001138e-03_rb,8.870154e-03_rb,8.753148e-03_rb,8.649907e-03_rb, &
       8.560232e-03_rb/)                                                             
      absice3(:,16) = (/ &
! band 16                                                                            
       1.688344e-01_rb,1.077072e-01_rb,7.994467e-02_rb,6.403862e-02_rb,5.369850e-02_rb, &
       4.641582e-02_rb,4.099331e-02_rb,3.678724e-02_rb,3.342069e-02_rb,3.065831e-02_rb, &
       2.834557e-02_rb,2.637680e-02_rb,2.467733e-02_rb,2.319286e-02_rb,2.188299e-02_rb, &
       2.071701e-02_rb,1.967121e-02_rb,1.872692e-02_rb,1.786931e-02_rb,1.708641e-02_rb, &
       1.636846e-02_rb,1.570743e-02_rb,1.509665e-02_rb,1.453052e-02_rb,1.400433e-02_rb, &
       1.351407e-02_rb,1.305631e-02_rb,1.262810e-02_rb,1.222688e-02_rb,1.185044e-02_rb, &
       1.149683e-02_rb,1.116436e-02_rb,1.085153e-02_rb,1.055701e-02_rb,1.027961e-02_rb, &
       1.001831e-02_rb,9.772141e-03_rb,9.540280e-03_rb,9.321966e-03_rb,9.116517e-03_rb, &
       8.923315e-03_rb,8.741803e-03_rb,8.571472e-03_rb,8.411860e-03_rb,8.262543e-03_rb, &
       8.123136e-03_rb/)                       

! For LIQFLAG = 0.                                                                   
      absliq0 = 0.0903614_rb                                                         
                                                                                     
! For LIQFLAG = 1.  In each band, the absorption                                     
! coefficients are listed for a range of effective radii from 2.5                    
! to 59.5 microns in increments of 1.0 micron.                                       
      absliq1(:, 1) = (/ &                                                           
! band  1                                                                            
       1.64047e-03_rb, 6.90533e-02_rb, 7.72017e-02_rb, 7.78054e-02_rb, 7.69523e-02_rb, &
       7.58058e-02_rb, 7.46400e-02_rb, 7.35123e-02_rb, 7.24162e-02_rb, 7.13225e-02_rb, &
       6.99145e-02_rb, 6.66409e-02_rb, 6.36582e-02_rb, 6.09425e-02_rb, 5.84593e-02_rb, &
       5.61743e-02_rb, 5.40571e-02_rb, 5.20812e-02_rb, 5.02245e-02_rb, 4.84680e-02_rb, &
       4.67959e-02_rb, 4.51944e-02_rb, 4.36516e-02_rb, 4.21570e-02_rb, 4.07015e-02_rb, &
       3.92766e-02_rb, 3.78747e-02_rb, 3.64886e-02_rb, 3.53632e-02_rb, 3.41992e-02_rb, &
       3.31016e-02_rb, 3.20643e-02_rb, 3.10817e-02_rb, 3.01490e-02_rb, 2.92620e-02_rb, &
       2.84171e-02_rb, 2.76108e-02_rb, 2.68404e-02_rb, 2.61031e-02_rb, 2.53966e-02_rb, &
       2.47189e-02_rb, 2.40678e-02_rb, 2.34418e-02_rb, 2.28392e-02_rb, 2.22586e-02_rb, &
       2.16986e-02_rb, 2.11580e-02_rb, 2.06356e-02_rb, 2.01305e-02_rb, 1.96417e-02_rb, &
       1.91682e-02_rb, 1.87094e-02_rb, 1.82643e-02_rb, 1.78324e-02_rb, 1.74129e-02_rb, &
       1.70052e-02_rb, 1.66088e-02_rb, 1.62231e-02_rb/)                              
      absliq1(:, 2) = (/ &                                                           
! band  2                                                                            
       2.19486e-01_rb, 1.80687e-01_rb, 1.59150e-01_rb, 1.44731e-01_rb, 1.33703e-01_rb, &
       1.24355e-01_rb, 1.15756e-01_rb, 1.07318e-01_rb, 9.86119e-02_rb, 8.92739e-02_rb, &
       8.34911e-02_rb, 7.70773e-02_rb, 7.15240e-02_rb, 6.66615e-02_rb, 6.23641e-02_rb, &
       5.85359e-02_rb, 5.51020e-02_rb, 5.20032e-02_rb, 4.91916e-02_rb, 4.66283e-02_rb, &
       4.42813e-02_rb, 4.21236e-02_rb, 4.01330e-02_rb, 3.82905e-02_rb, 3.65797e-02_rb, &
       3.49869e-02_rb, 3.35002e-02_rb, 3.21090e-02_rb, 3.08957e-02_rb, 2.97601e-02_rb, &
       2.86966e-02_rb, 2.76984e-02_rb, 2.67599e-02_rb, 2.58758e-02_rb, 2.50416e-02_rb, &
       2.42532e-02_rb, 2.35070e-02_rb, 2.27997e-02_rb, 2.21284e-02_rb, 2.14904e-02_rb, &
       2.08834e-02_rb, 2.03051e-02_rb, 1.97536e-02_rb, 1.92271e-02_rb, 1.87239e-02_rb, &
       1.82425e-02_rb, 1.77816e-02_rb, 1.73399e-02_rb, 1.69162e-02_rb, 1.65094e-02_rb, &
       1.61187e-02_rb, 1.57430e-02_rb, 1.53815e-02_rb, 1.50334e-02_rb, 1.46981e-02_rb, &
       1.43748e-02_rb, 1.40628e-02_rb, 1.37617e-02_rb/)                              
      absliq1(:, 3) = (/ &                                                           
! band  3                                                                            
       2.95174e-01_rb, 2.34765e-01_rb, 1.98038e-01_rb, 1.72114e-01_rb, 1.52083e-01_rb, &
       1.35654e-01_rb, 1.21613e-01_rb, 1.09252e-01_rb, 9.81263e-02_rb, 8.79448e-02_rb, &
       8.12566e-02_rb, 7.44563e-02_rb, 6.86374e-02_rb, 6.36042e-02_rb, 5.92094e-02_rb, &
       5.53402e-02_rb, 5.19087e-02_rb, 4.88455e-02_rb, 4.60951e-02_rb, 4.36124e-02_rb, &
       4.13607e-02_rb, 3.93096e-02_rb, 3.74338e-02_rb, 3.57119e-02_rb, 3.41261e-02_rb, &
       3.26610e-02_rb, 3.13036e-02_rb, 3.00425e-02_rb, 2.88497e-02_rb, 2.78077e-02_rb, &
       2.68317e-02_rb, 2.59158e-02_rb, 2.50545e-02_rb, 2.42430e-02_rb, 2.34772e-02_rb, &
       2.27533e-02_rb, 2.20679e-02_rb, 2.14181e-02_rb, 2.08011e-02_rb, 2.02145e-02_rb, &
       1.96561e-02_rb, 1.91239e-02_rb, 1.86161e-02_rb, 1.81311e-02_rb, 1.76673e-02_rb, &
       1.72234e-02_rb, 1.67981e-02_rb, 1.63903e-02_rb, 1.59989e-02_rb, 1.56230e-02_rb, &
       1.52615e-02_rb, 1.49138e-02_rb, 1.45791e-02_rb, 1.42565e-02_rb, 1.39455e-02_rb, &
       1.36455e-02_rb, 1.33559e-02_rb, 1.30761e-02_rb/)                              
      absliq1(:, 4) = (/ &                                                           
! band  4                                                                            
       3.00925e-01_rb, 2.36949e-01_rb, 1.96947e-01_rb, 1.68692e-01_rb, 1.47190e-01_rb, &
       1.29986e-01_rb, 1.15719e-01_rb, 1.03568e-01_rb, 9.30028e-02_rb, 8.36658e-02_rb, &
       7.71075e-02_rb, 7.07002e-02_rb, 6.52284e-02_rb, 6.05024e-02_rb, 5.63801e-02_rb, &
       5.27534e-02_rb, 4.95384e-02_rb, 4.66690e-02_rb, 4.40925e-02_rb, 4.17664e-02_rb, &
       3.96559e-02_rb, 3.77326e-02_rb, 3.59727e-02_rb, 3.43561e-02_rb, 3.28662e-02_rb, &
       3.14885e-02_rb, 3.02110e-02_rb, 2.90231e-02_rb, 2.78948e-02_rb, 2.69109e-02_rb, &
       2.59884e-02_rb, 2.51217e-02_rb, 2.43058e-02_rb, 2.35364e-02_rb, 2.28096e-02_rb, &
       2.21218e-02_rb, 2.14700e-02_rb, 2.08515e-02_rb, 2.02636e-02_rb, 1.97041e-02_rb, &
       1.91711e-02_rb, 1.86625e-02_rb, 1.81769e-02_rb, 1.77126e-02_rb, 1.72683e-02_rb, &
       1.68426e-02_rb, 1.64344e-02_rb, 1.60427e-02_rb, 1.56664e-02_rb, 1.53046e-02_rb, &
       1.49565e-02_rb, 1.46214e-02_rb, 1.42985e-02_rb, 1.39871e-02_rb, 1.36866e-02_rb, &
       1.33965e-02_rb, 1.31162e-02_rb, 1.28453e-02_rb/)                              
      absliq1(:, 5) = (/ &                                                           
! band  5                                                                            
       2.64691e-01_rb, 2.12018e-01_rb, 1.78009e-01_rb, 1.53539e-01_rb, 1.34721e-01_rb, &
       1.19580e-01_rb, 1.06996e-01_rb, 9.62772e-02_rb, 8.69710e-02_rb, 7.87670e-02_rb, &
       7.29272e-02_rb, 6.70920e-02_rb, 6.20977e-02_rb, 5.77732e-02_rb, 5.39910e-02_rb, &
       5.06538e-02_rb, 4.76866e-02_rb, 4.50301e-02_rb, 4.26374e-02_rb, 4.04704e-02_rb, &
       3.84981e-02_rb, 3.66948e-02_rb, 3.50394e-02_rb, 3.35141e-02_rb, 3.21038e-02_rb, &
       3.07957e-02_rb, 2.95788e-02_rb, 2.84438e-02_rb, 2.73790e-02_rb, 2.64390e-02_rb, &
       2.55565e-02_rb, 2.47263e-02_rb, 2.39437e-02_rb, 2.32047e-02_rb, 2.25056e-02_rb, &
       2.18433e-02_rb, 2.12149e-02_rb, 2.06177e-02_rb, 2.00495e-02_rb, 1.95081e-02_rb, &
       1.89917e-02_rb, 1.84984e-02_rb, 1.80269e-02_rb, 1.75755e-02_rb, 1.71431e-02_rb, &
       1.67283e-02_rb, 1.63303e-02_rb, 1.59478e-02_rb, 1.55801e-02_rb, 1.52262e-02_rb, &
       1.48853e-02_rb, 1.45568e-02_rb, 1.42400e-02_rb, 1.39342e-02_rb, 1.36388e-02_rb, &
       1.33533e-02_rb, 1.30773e-02_rb, 1.28102e-02_rb/)                              
      absliq1(:, 6) = (/ &                                                           
! band  6                                                                            
       8.81182e-02_rb, 1.06745e-01_rb, 9.79753e-02_rb, 8.99625e-02_rb, 8.35200e-02_rb, &
       7.81899e-02_rb, 7.35939e-02_rb, 6.94696e-02_rb, 6.56266e-02_rb, 6.19148e-02_rb, &
       5.83355e-02_rb, 5.49306e-02_rb, 5.19642e-02_rb, 4.93325e-02_rb, 4.69659e-02_rb, &
       4.48148e-02_rb, 4.28431e-02_rb, 4.10231e-02_rb, 3.93332e-02_rb, 3.77563e-02_rb, &
       3.62785e-02_rb, 3.48882e-02_rb, 3.35758e-02_rb, 3.23333e-02_rb, 3.11536e-02_rb, &
       3.00310e-02_rb, 2.89601e-02_rb, 2.79365e-02_rb, 2.70502e-02_rb, 2.62618e-02_rb, &
       2.55025e-02_rb, 2.47728e-02_rb, 2.40726e-02_rb, 2.34013e-02_rb, 2.27583e-02_rb, &
       2.21422e-02_rb, 2.15522e-02_rb, 2.09869e-02_rb, 2.04453e-02_rb, 1.99260e-02_rb, &
       1.94280e-02_rb, 1.89501e-02_rb, 1.84913e-02_rb, 1.80506e-02_rb, 1.76270e-02_rb, &
       1.72196e-02_rb, 1.68276e-02_rb, 1.64500e-02_rb, 1.60863e-02_rb, 1.57357e-02_rb, &
       1.53975e-02_rb, 1.50710e-02_rb, 1.47558e-02_rb, 1.44511e-02_rb, 1.41566e-02_rb, &
       1.38717e-02_rb, 1.35960e-02_rb, 1.33290e-02_rb/)                              
      absliq1(:, 7) = (/ &                                                           
! band  7                                                                            
       4.32174e-02_rb, 7.36078e-02_rb, 6.98340e-02_rb, 6.65231e-02_rb, 6.41948e-02_rb, &
       6.23551e-02_rb, 6.06638e-02_rb, 5.88680e-02_rb, 5.67124e-02_rb, 5.38629e-02_rb, &
       4.99579e-02_rb, 4.86289e-02_rb, 4.70120e-02_rb, 4.52854e-02_rb, 4.35466e-02_rb, &
       4.18480e-02_rb, 4.02169e-02_rb, 3.86658e-02_rb, 3.71992e-02_rb, 3.58168e-02_rb, &
       3.45155e-02_rb, 3.32912e-02_rb, 3.21390e-02_rb, 3.10538e-02_rb, 3.00307e-02_rb, &
       2.90651e-02_rb, 2.81524e-02_rb, 2.72885e-02_rb, 2.62821e-02_rb, 2.55744e-02_rb, &
       2.48799e-02_rb, 2.42029e-02_rb, 2.35460e-02_rb, 2.29108e-02_rb, 2.22981e-02_rb, &
       2.17079e-02_rb, 2.11402e-02_rb, 2.05945e-02_rb, 2.00701e-02_rb, 1.95663e-02_rb, &
       1.90824e-02_rb, 1.86174e-02_rb, 1.81706e-02_rb, 1.77411e-02_rb, 1.73281e-02_rb, &
       1.69307e-02_rb, 1.65483e-02_rb, 1.61801e-02_rb, 1.58254e-02_rb, 1.54835e-02_rb, &
       1.51538e-02_rb, 1.48358e-02_rb, 1.45288e-02_rb, 1.42322e-02_rb, 1.39457e-02_rb, &
       1.36687e-02_rb, 1.34008e-02_rb, 1.31416e-02_rb/)                              
      absliq1(:, 8) = (/ &                                                           
! band  8                                                                            
       1.41881e-01_rb, 7.15419e-02_rb, 6.30335e-02_rb, 6.11132e-02_rb, 6.01931e-02_rb, &
       5.92420e-02_rb, 5.78968e-02_rb, 5.58876e-02_rb, 5.28923e-02_rb, 4.84462e-02_rb, &
       4.60839e-02_rb, 4.56013e-02_rb, 4.45410e-02_rb, 4.31866e-02_rb, 4.17026e-02_rb, &
       4.01850e-02_rb, 3.86892e-02_rb, 3.72461e-02_rb, 3.58722e-02_rb, 3.45749e-02_rb, &
       3.33564e-02_rb, 3.22155e-02_rb, 3.11494e-02_rb, 3.01541e-02_rb, 2.92253e-02_rb, &
       2.83584e-02_rb, 2.75488e-02_rb, 2.67925e-02_rb, 2.57692e-02_rb, 2.50704e-02_rb, &
       2.43918e-02_rb, 2.37350e-02_rb, 2.31005e-02_rb, 2.24888e-02_rb, 2.18996e-02_rb, &
       2.13325e-02_rb, 2.07870e-02_rb, 2.02623e-02_rb, 1.97577e-02_rb, 1.92724e-02_rb, &
       1.88056e-02_rb, 1.83564e-02_rb, 1.79241e-02_rb, 1.75079e-02_rb, 1.71070e-02_rb, &
       1.67207e-02_rb, 1.63482e-02_rb, 1.59890e-02_rb, 1.56424e-02_rb, 1.53077e-02_rb, &
       1.49845e-02_rb, 1.46722e-02_rb, 1.43702e-02_rb, 1.40782e-02_rb, 1.37955e-02_rb, &
       1.35219e-02_rb, 1.32569e-02_rb, 1.30000e-02_rb/)                              
      absliq1(:, 9) = (/ &                                                           
! band  9                                                                            
       6.72726e-02_rb, 6.61013e-02_rb, 6.47866e-02_rb, 6.33780e-02_rb, 6.18985e-02_rb, &
       6.03335e-02_rb, 5.86136e-02_rb, 5.65876e-02_rb, 5.39839e-02_rb, 5.03536e-02_rb, &
       4.71608e-02_rb, 4.63630e-02_rb, 4.50313e-02_rb, 4.34526e-02_rb, 4.17876e-02_rb, &
       4.01261e-02_rb, 3.85171e-02_rb, 3.69860e-02_rb, 3.55442e-02_rb, 3.41954e-02_rb, &
       3.29384e-02_rb, 3.17693e-02_rb, 3.06832e-02_rb, 2.96745e-02_rb, 2.87374e-02_rb, &
       2.78662e-02_rb, 2.70557e-02_rb, 2.63008e-02_rb, 2.52450e-02_rb, 2.45424e-02_rb, &
       2.38656e-02_rb, 2.32144e-02_rb, 2.25885e-02_rb, 2.19873e-02_rb, 2.14099e-02_rb, &
       2.08554e-02_rb, 2.03230e-02_rb, 1.98116e-02_rb, 1.93203e-02_rb, 1.88482e-02_rb, &
       1.83944e-02_rb, 1.79578e-02_rb, 1.75378e-02_rb, 1.71335e-02_rb, 1.67440e-02_rb, &
       1.63687e-02_rb, 1.60069e-02_rb, 1.56579e-02_rb, 1.53210e-02_rb, 1.49958e-02_rb, &
       1.46815e-02_rb, 1.43778e-02_rb, 1.40841e-02_rb, 1.37999e-02_rb, 1.35249e-02_rb, &
       1.32585e-02_rb, 1.30004e-02_rb, 1.27502e-02_rb/)                              
      absliq1(:,10) = (/ &                                                           
! band 10                                                                            
       7.97040e-02_rb, 7.63844e-02_rb, 7.36499e-02_rb, 7.13525e-02_rb, 6.93043e-02_rb, &
       6.72807e-02_rb, 6.50227e-02_rb, 6.22395e-02_rb, 5.86093e-02_rb, 5.37815e-02_rb, &
       5.14682e-02_rb, 4.97214e-02_rb, 4.77392e-02_rb, 4.56961e-02_rb, 4.36858e-02_rb, &
       4.17569e-02_rb, 3.99328e-02_rb, 3.82224e-02_rb, 3.66265e-02_rb, 3.51416e-02_rb, &
       3.37617e-02_rb, 3.24798e-02_rb, 3.12887e-02_rb, 3.01812e-02_rb, 2.91505e-02_rb, &
       2.81900e-02_rb, 2.72939e-02_rb, 2.64568e-02_rb, 2.54165e-02_rb, 2.46832e-02_rb, &
       2.39783e-02_rb, 2.33017e-02_rb, 2.26531e-02_rb, 2.20314e-02_rb, 2.14359e-02_rb, &
       2.08653e-02_rb, 2.03187e-02_rb, 1.97947e-02_rb, 1.92924e-02_rb, 1.88106e-02_rb, &
       1.83483e-02_rb, 1.79043e-02_rb, 1.74778e-02_rb, 1.70678e-02_rb, 1.66735e-02_rb, &
       1.62941e-02_rb, 1.59286e-02_rb, 1.55766e-02_rb, 1.52371e-02_rb, 1.49097e-02_rb, &
       1.45937e-02_rb, 1.42885e-02_rb, 1.39936e-02_rb, 1.37085e-02_rb, 1.34327e-02_rb, &
       1.31659e-02_rb, 1.29075e-02_rb, 1.26571e-02_rb/)                              
      absliq1(:,11) = (/ &                                                           
! band 11                                                                            
       1.49438e-01_rb, 1.33535e-01_rb, 1.21542e-01_rb, 1.11743e-01_rb, 1.03263e-01_rb, &
       9.55774e-02_rb, 8.83382e-02_rb, 8.12943e-02_rb, 7.42533e-02_rb, 6.70609e-02_rb, &
       6.38761e-02_rb, 5.97788e-02_rb, 5.59841e-02_rb, 5.25318e-02_rb, 4.94132e-02_rb, &
       4.66014e-02_rb, 4.40644e-02_rb, 4.17706e-02_rb, 3.96910e-02_rb, 3.77998e-02_rb, &
       3.60742e-02_rb, 3.44947e-02_rb, 3.30442e-02_rb, 3.17079e-02_rb, 3.04730e-02_rb, &
       2.93283e-02_rb, 2.82642e-02_rb, 2.72720e-02_rb, 2.61789e-02_rb, 2.53277e-02_rb, &
       2.45237e-02_rb, 2.37635e-02_rb, 2.30438e-02_rb, 2.23615e-02_rb, 2.17140e-02_rb, &
       2.10987e-02_rb, 2.05133e-02_rb, 1.99557e-02_rb, 1.94241e-02_rb, 1.89166e-02_rb, &
       1.84317e-02_rb, 1.79679e-02_rb, 1.75238e-02_rb, 1.70983e-02_rb, 1.66901e-02_rb, &
       1.62983e-02_rb, 1.59219e-02_rb, 1.55599e-02_rb, 1.52115e-02_rb, 1.48761e-02_rb, &
       1.45528e-02_rb, 1.42411e-02_rb, 1.39402e-02_rb, 1.36497e-02_rb, 1.33690e-02_rb, &
       1.30976e-02_rb, 1.28351e-02_rb, 1.25810e-02_rb/)                              
      absliq1(:,12) = (/ &                                                           
! band 12                                                                            
       3.71985e-02_rb, 3.88586e-02_rb, 3.99070e-02_rb, 4.04351e-02_rb, 4.04610e-02_rb, &
       3.99834e-02_rb, 3.89953e-02_rb, 3.74886e-02_rb, 3.54551e-02_rb, 3.28870e-02_rb, &
       3.32576e-02_rb, 3.22444e-02_rb, 3.12384e-02_rb, 3.02584e-02_rb, 2.93146e-02_rb, &
       2.84120e-02_rb, 2.75525e-02_rb, 2.67361e-02_rb, 2.59618e-02_rb, 2.52280e-02_rb, &
       2.45327e-02_rb, 2.38736e-02_rb, 2.32487e-02_rb, 2.26558e-02_rb, 2.20929e-02_rb, &
       2.15579e-02_rb, 2.10491e-02_rb, 2.05648e-02_rb, 1.99749e-02_rb, 1.95704e-02_rb, &
       1.91731e-02_rb, 1.87839e-02_rb, 1.84032e-02_rb, 1.80315e-02_rb, 1.76689e-02_rb, &
       1.73155e-02_rb, 1.69712e-02_rb, 1.66362e-02_rb, 1.63101e-02_rb, 1.59928e-02_rb, &
       1.56842e-02_rb, 1.53840e-02_rb, 1.50920e-02_rb, 1.48080e-02_rb, 1.45318e-02_rb, &
       1.42631e-02_rb, 1.40016e-02_rb, 1.37472e-02_rb, 1.34996e-02_rb, 1.32586e-02_rb, &
       1.30239e-02_rb, 1.27954e-02_rb, 1.25728e-02_rb, 1.23559e-02_rb, 1.21445e-02_rb, &
       1.19385e-02_rb, 1.17376e-02_rb, 1.15417e-02_rb/)                              

      absliq1(:,13) = (/ &                                                           
! band 13                                                                            
       3.11868e-02_rb, 4.48357e-02_rb, 4.90224e-02_rb, 4.96406e-02_rb, 4.86806e-02_rb, &
       4.69610e-02_rb, 4.48630e-02_rb, 4.25795e-02_rb, 4.02138e-02_rb, 3.78236e-02_rb, &
       3.74266e-02_rb, 3.60384e-02_rb, 3.47074e-02_rb, 3.34434e-02_rb, 3.22499e-02_rb, &
       3.11264e-02_rb, 3.00704e-02_rb, 2.90784e-02_rb, 2.81463e-02_rb, 2.72702e-02_rb, &
       2.64460e-02_rb, 2.56698e-02_rb, 2.49381e-02_rb, 2.42475e-02_rb, 2.35948e-02_rb, &
       2.29774e-02_rb, 2.23925e-02_rb, 2.18379e-02_rb, 2.11793e-02_rb, 2.07076e-02_rb, &
       2.02470e-02_rb, 1.97981e-02_rb, 1.93613e-02_rb, 1.89367e-02_rb, 1.85243e-02_rb, &
       1.81240e-02_rb, 1.77356e-02_rb, 1.73588e-02_rb, 1.69935e-02_rb, 1.66392e-02_rb, &
       1.62956e-02_rb, 1.59624e-02_rb, 1.56393e-02_rb, 1.53259e-02_rb, 1.50219e-02_rb, &
       1.47268e-02_rb, 1.44404e-02_rb, 1.41624e-02_rb, 1.38925e-02_rb, 1.36302e-02_rb, &
       1.33755e-02_rb, 1.31278e-02_rb, 1.28871e-02_rb, 1.26530e-02_rb, 1.24253e-02_rb, &
       1.22038e-02_rb, 1.19881e-02_rb, 1.17782e-02_rb/)                              
      absliq1(:,14) = (/ &                                                           
! band 14                                                                            
       1.58988e-02_rb, 3.50652e-02_rb, 4.00851e-02_rb, 4.07270e-02_rb, 3.98101e-02_rb, &
       3.83306e-02_rb, 3.66829e-02_rb, 3.50327e-02_rb, 3.34497e-02_rb, 3.19609e-02_rb, &
       3.13712e-02_rb, 3.03348e-02_rb, 2.93415e-02_rb, 2.83973e-02_rb, 2.75037e-02_rb, &
       2.66604e-02_rb, 2.58654e-02_rb, 2.51161e-02_rb, 2.44100e-02_rb, 2.37440e-02_rb, &
       2.31154e-02_rb, 2.25215e-02_rb, 2.19599e-02_rb, 2.14282e-02_rb, 2.09242e-02_rb, &
       2.04459e-02_rb, 1.99915e-02_rb, 1.95594e-02_rb, 1.90254e-02_rb, 1.86598e-02_rb, &
       1.82996e-02_rb, 1.79455e-02_rb, 1.75983e-02_rb, 1.72584e-02_rb, 1.69260e-02_rb, &
       1.66013e-02_rb, 1.62843e-02_rb, 1.59752e-02_rb, 1.56737e-02_rb, 1.53799e-02_rb, &
       1.50936e-02_rb, 1.48146e-02_rb, 1.45429e-02_rb, 1.42782e-02_rb, 1.40203e-02_rb, &
       1.37691e-02_rb, 1.35243e-02_rb, 1.32858e-02_rb, 1.30534e-02_rb, 1.28270e-02_rb, &
       1.26062e-02_rb, 1.23909e-02_rb, 1.21810e-02_rb, 1.19763e-02_rb, 1.17766e-02_rb, &
       1.15817e-02_rb, 1.13915e-02_rb, 1.12058e-02_rb/)                              
      absliq1(:,15) = (/ &                                                           
! band 15                                                                            
       5.02079e-03_rb, 2.17615e-02_rb, 2.55449e-02_rb, 2.59484e-02_rb, 2.53650e-02_rb, &
       2.45281e-02_rb, 2.36843e-02_rb, 2.29159e-02_rb, 2.22451e-02_rb, 2.16716e-02_rb, &
       2.11451e-02_rb, 2.05817e-02_rb, 2.00454e-02_rb, 1.95372e-02_rb, 1.90567e-02_rb, &
       1.86028e-02_rb, 1.81742e-02_rb, 1.77693e-02_rb, 1.73866e-02_rb, 1.70244e-02_rb, &
       1.66815e-02_rb, 1.63563e-02_rb, 1.60477e-02_rb, 1.57544e-02_rb, 1.54755e-02_rb, &
       1.52097e-02_rb, 1.49564e-02_rb, 1.47146e-02_rb, 1.43684e-02_rb, 1.41728e-02_rb, &
       1.39762e-02_rb, 1.37797e-02_rb, 1.35838e-02_rb, 1.33891e-02_rb, 1.31961e-02_rb, &
       1.30051e-02_rb, 1.28164e-02_rb, 1.26302e-02_rb, 1.24466e-02_rb, 1.22659e-02_rb, &
       1.20881e-02_rb, 1.19131e-02_rb, 1.17412e-02_rb, 1.15723e-02_rb, 1.14063e-02_rb, &
       1.12434e-02_rb, 1.10834e-02_rb, 1.09264e-02_rb, 1.07722e-02_rb, 1.06210e-02_rb, &
       1.04725e-02_rb, 1.03269e-02_rb, 1.01839e-02_rb, 1.00436e-02_rb, 9.90593e-03_rb, &
       9.77080e-03_rb, 9.63818e-03_rb, 9.50800e-03_rb/)                              
      absliq1(:,16) = (/ &                                                           
! band 16                                                                            
       5.64971e-02_rb, 9.04736e-02_rb, 8.11726e-02_rb, 7.05450e-02_rb, 6.20052e-02_rb, &
       5.54286e-02_rb, 5.03503e-02_rb, 4.63791e-02_rb, 4.32290e-02_rb, 4.06959e-02_rb, &
       3.74690e-02_rb, 3.52964e-02_rb, 3.33799e-02_rb, 3.16774e-02_rb, 3.01550e-02_rb, &
       2.87856e-02_rb, 2.75474e-02_rb, 2.64223e-02_rb, 2.53953e-02_rb, 2.44542e-02_rb, &
       2.35885e-02_rb, 2.27894e-02_rb, 2.20494e-02_rb, 2.13622e-02_rb, 2.07222e-02_rb, &
       2.01246e-02_rb, 1.95654e-02_rb, 1.90408e-02_rb, 1.84398e-02_rb, 1.80021e-02_rb, &
       1.75816e-02_rb, 1.71775e-02_rb, 1.67889e-02_rb, 1.64152e-02_rb, 1.60554e-02_rb, &
       1.57089e-02_rb, 1.53751e-02_rb, 1.50531e-02_rb, 1.47426e-02_rb, 1.44428e-02_rb, &
       1.41532e-02_rb, 1.38734e-02_rb, 1.36028e-02_rb, 1.33410e-02_rb, 1.30875e-02_rb, &
       1.28420e-02_rb, 1.26041e-02_rb, 1.23735e-02_rb, 1.21497e-02_rb, 1.19325e-02_rb, &
       1.17216e-02_rb, 1.15168e-02_rb, 1.13177e-02_rb, 1.11241e-02_rb, 1.09358e-02_rb, &
       1.07525e-02_rb, 1.05741e-02_rb, 1.04003e-02_rb/)                              

!jm not thread safe      hvrclc = '$Revision: 1.8 $'

      ncbands = 1

! This initialization is done in rrtmg_lw_subcol.F90.
!      do lay = 1, nlayers
!         do ig = 1, ngptlw
!            taucmc(ig,lay) = 0.0_rb
!         enddo
!      enddo

! Main layer loop
      do lay = 1, nlayers

        do ig = 1, ngptlw
          cwp = ciwpmc(ig,lay) + clwpmc(ig,lay) + cswpmc(ig,lay)
          if (cldfmc(ig,lay) .ge. cldmin .and.                          &
     &      (cwp .ge. cldmin .or. taucmc(ig,lay) .ge. cldmin)) then


! Ice clouds and water clouds combined.
            if (inflag .eq. 0) then
! Cloud optical depth already defined in taucmc, return to main program
               return

            elseif(inflag .eq. 1) then
                stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'
!               cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
!               taucmc(ig,lay) = abscld1 * cwp

! Separate treatement of ice clouds and water clouds.
            elseif(inflag .ge. 2) then
               radice = reicmc(lay)

! Calculation of absorption coefficients due to ice clouds.
               if ((ciwpmc(ig,lay)+cswpmc(ig,lay)) .eq. 0.0_rb) then
                  abscoice(ig) = 0.0_rb
                  abscosno(ig) = 0.0_rb

               elseif (iceflag .eq. 0) then
                  if (radice .lt. 10.0_rb) stop 'ICE RADIUS TOO SMALL'
                  abscoice(ig) = absice0(1) + absice0(2)/radice
                  abscosno(ig) = 0.0_rb

               elseif (iceflag .eq. 1) then
                  if (radice .lt. 13.0_rb .or. radice .gt. 130._rb) stop&
     &                'ICE RADIUS OUT OF BOUNDS'
                  ncbands = 5
                  ib = icb(ngb(ig))
                  abscoice(ig) = absice1(1,ib) + absice1(2,ib)/radice
                  abscosno(ig) = 0.0_rb

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

               elseif (iceflag .eq. 2) then
                  if (radice .lt. 5.0_rb .or. radice .gt. 131.0_rb) stop&
     &                'ICE RADIUS OUT OF BOUNDS'
                     ncbands = 16
                     factor = (radice - 2._rb)/3._rb
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice(ig) =                                     &
     &                   absice2(index,ib) + fint *                     &
     &                   (absice2(index+1,ib) - (absice2(index,ib)))
                     abscosno(ig) = 0.0_rb

! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

               elseif (iceflag .ge. 3) then
                  if (radice .lt. 5.0_rb .or. radice .gt. 140.0_rb) then
                         write(errmsg,'(a,i5,i5,f8.2,f8.2)' )           &
     &         'ERROR: ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'    &
     &          ,ig, lay, ciwpmc(ig,lay), radice
                         errflg = 1
                         return
                     end if                                                                      
                     ncbands = 16                                                                
                     factor = (radice - 2._rb)/3._rb                                             
                     index = int(factor)                                                         
                     if (index .eq. 46) index = 45                                               
                     fint = factor - float(index)                                                
                     ib = ngb(ig)                                                                
                     abscoice(ig) =                                     &
     &                   absice3(index,ib) + fint *                     &
     &                   (absice3(index+1,ib) - (absice3(index,ib)))                             
                     abscosno(ig) = 0.0_rb                                                       
                                                                                                 
               endif                                                                             
                                                                                                 
!..Incorporate additional effects due to snow.                                                   
               if (cswpmc(ig,lay).gt.0.0_rb .and. iceflag .eq. 5) then                           
                  radsno = resnmc(lay)                                                           
                  if (radsno .lt. 5.0_rb .or. radsno .gt. 140.0_rb) then                         
                         write(errmsg,'(a,i5,i5,f8.2,f8.2)' )           &
     &         'ERROR: SNOW GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'   &                        
     &         ,ig, lay, cswpmc(ig,lay), radsno                                                  
                         errflg = 1
                         return
                     end if                                                                      
                     ncbands = 16                                                                
                     factor = (radsno - 2._rb)/3._rb                                             
                     index = int(factor)                                                         
                     if (index .eq. 46) index = 45                                               
                     fint = factor - float(index)                                                
                     ib = ngb(ig)                                                                
                     abscosno(ig) =                                     &
     &                   absice3(index,ib) + fint *                     &
     &                  (absice3(index+1,ib) - (absice3(index,ib)))                             
               endif                                                                             
                                                                                                 
                                                                                                 
                                                                                                 
! Calculation of absorption coefficients due to water clouds.                                    
               if (clwpmc(ig,lay) .eq. 0.0_rb) then                                              
                  abscoliq(ig) = 0.0_rb                                                          
                                                                                                 
               elseif (liqflag .eq. 0) then                                                      
                   abscoliq(ig) = absliq0                                                        
                                                                                                 
               elseif (liqflag .eq. 1) then                                                      
                  radliq = relqmc(lay)                        
                  if (radliq .lt. 2.5_rb .or. radliq .gt. 60._rb) then
                     write(errmsg,'(a,i5,i5,f8.2,f8.2)' )              &
&                         'ERROR: LIQUID EFFECTIVE SIZE OUT OF BOUNDS' &
&                         ,ig, lay, clwpmc(ig,lay), radliq
                     errflg = 1
                     return
                  end if
                  index = int(radliq - 1.5_rb)                                                   
                  if (index .eq. 0) index = 1                                                    
                  if (index .eq. 58) index = 57                                                  
                  fint = radliq - 1.5_rb - float(index)                                          
                  ib = ngb(ig)                                                                   
                  abscoliq(ig) =                                        &
     &                  absliq1(index,ib) + fint *                      &
     &                  (absliq1(index+1,ib) - (absliq1(index,ib)))                              
               endif                                                                             
                                                                                                 
               taucmc(ig,lay) = ciwpmc(ig,lay) * abscoice(ig) +         &
     &                           clwpmc(ig,lay) * abscoliq(ig) +        &
     &                           cswpmc(ig,lay) * abscosno(ig)                                    
                                                                                                 
            endif                                                                                
         endif                                                                                   
         enddo                                                                                   
      enddo                                                                                      
                                                                                                 
      end subroutine cldprmc                                                                     
                                                                    

!........................................!$
      end module rrtmg_lw                !$
!========================================!$

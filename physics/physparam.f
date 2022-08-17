!> \file physparam.f
!! This file contains module physparam.

!  ==========================================================  !!!!!
!                    module physparam description              !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!     This module defines commonly used control variables/parameters   !
!     in physics related programs.                                     !
!                                                                      !
!     Section 1 contains control variables defined in the form of      !
!     parameter. They are pre-determined choices and not adjustable    !
!     during model's run-time.                                         !
!                                                                      !
!     Section 2 contains control variables defined as module variables.!
!     They are more flexible to be changed during run-time by either   !
!     through input namelist, or through model environment condition.  !
!     They are preassigned here as the default values.                 !
!                                                                      !
!!!!!  ==========================================================  !!!!!

!> \defgroup phy_sparam GFS Physics Parameter Module
!! Those variables are grouped  together in accordance with functionaity
!! and are given brief descriptions and value specifications. There are
!! two types of attributes (parameters vs. save) designated for the
!! control variables. Those with a "parameter" attribute are prescribed
!! with a preferred option value, while the ones with a "save" attribute
!! are given a default value but could be changed at the model's
!! execution-time (usually through an input of name-list file or through
!! run scripts).

!> This module defines commonly used control variables and parameters
!! in physics related programs.
      module physparam            
!
!     implicit   none

!  --- ...  define kind parameters here

!   ** if already exist, use the module containing kind definitions
      use machine

!   ** otherwise, define kind parameter here
!     implicit   none
!     integer, public, parameter :: kind_io4 = 4
!     integer, public, parameter :: kind_io8 = 8
!     integer, public, parameter :: kind_phys= selected_real_kind(13,60) ! the '60' maps to 64-bit real
!      .....

!     implicit   none
!
      public

!==================================================================================
!  Section - 1 -
!     control flags are pre-set as run-time non-adjuztable parameters.
!==================================================================================

! ............................................. !
!> \name  1.1 Control flags for SW radiation
! ............................................. !

!> SW heating rate unit control flag: =1:k/day; =2:k/second.
      integer,parameter :: iswrate = 2

!> SW minor gases effect control flag (CH4 and O2): =0:no; =1:yes.
!!\n =0: minor gases' effects are not included in calculations
!!\n =1: minor gases' effects are included in calculations
      integer,parameter :: iswrgas = 1

!> SW optical property for liquid clouds
!!\n =0:input cld opt depth, ignoring iswcice setting
!!\n =1:cloud optical property scheme based on Hu and Stamnes(1993) \cite hu_and_stamnes_1993 method
!!\n =2:cloud optical property scheme based on Hu and Stamnes(1993) -updated
      integer,save      :: iswcliq = 1

!> SW optical property for ice clouds (only iswcliq>0)
!!\n =1:optical property scheme based on Ebert and Curry (1992)
!!      \cite ebert_and_curry_1992 method
!!\n =2:optical property scheme based on Streamer v3.0
!!      \cite key_2002 method
!!\n =3:optical property scheme based on Fu's method (1996)
!!      \cite fu_1996 method
      integer,save      :: iswcice = 3

!> SW control flag for scattering process approximation
!!\n =1:two-stream delta-eddington    (Joseph et al. 1976
!!                         \cite joseph_et_al_1976)
!!\n =2:two-stream PIFM               (Zdunkowski et al. 1980
!!                         \cite zdunkowski_et_al_1980)
!!\n =3:discrete ordinates (Liou, 1973
!!                         \cite liou_1973)
      integer,parameter :: iswmode = 2

! ............................................. !
!> \name  1.2 Control flags for LW radiation
! ............................................. !

!> LW heating rate unit: =1:k/day; =2:k/second.
      integer,parameter :: ilwrate = 2

!> LW minor gases effect control flag (CH4,N2O,O2,and some CFCs):
!!\n =0: minor gases' effects are not included in calculations
!!\n =1: minor gases' effects are included in calculations
      integer,parameter :: ilwrgas = 1

!> LW optical property scheme for liquid clouds
!!\n =0:input cloud optical properties directly, not computed within
!!\n =1:input cwp,rew, use Hu and Stamnes(1993)
!!      \cite hu_and_stamnes_1993 method
      integer,save      :: ilwcliq = 1

!> LW optical property scheme for ice clouds (only ilwcliq>0)
!!\n =1:optical property scheme based on Ebert and Curry (1992)
!!      \cite ebert_and_curry_1992 method
!!\n =2:optical property scheme based on Streamer v3
!!      \cite key_2002 method
!!\n =3:optical property scheme use Fu's method (1998)
!!      \cite fu_et_al_1998 method
      integer,save      :: ilwcice = 3

!==================================================================================
!  Section - 2 -
!     values of control flags might be re-set in initialization subroutines
!       (may be adjusted at run time based on namelist input or run condition)
!==================================================================================

! ............................................. !
!> \name  2.3 For module radiation_gases
! ............................................. !

!> co2 data source control flag
!!\n =0:prescribed value(380 ppmv)
!!\n =1:yearly global averaged annual mean from observations
!!\n =2:monthly 15 degree horizontal resolution from observations
!!\n Opr GFS/CFS=2; see ICO2 in run scripts
      integer, save :: ico2flg = 0

!> controls external data at initial time and data usage during
!! forecast time
!!\n =-2:as in 0,but superimpose with seasonal climatology cycle
!!\n =-1:use user data,no extrapolation in overtime
!!\n =0:use IC time to select data,no extrapolation in overtime
!!\n =1:use forecast time to select data,extrapolate when necessary
!!\n =yyyy0:use yyyy year of data, no extrapolation
!!\n =yyyy1:use yyyy year of data, extrapolate when necessary
!!\n Opr GFS/CFS=1; see ICTM in run scripts
      integer, save :: ictmflg = 0

!> ozone data source control flag
!!\n =0:use seasonal climatology ozone data
!!\n >0:use prognostic ozone scheme (also depend on other model control
!!      variable at initial time)
      integer, save :: ioznflg = 1

!> external co2 2d monthly obsv data table: co2historicaldata_2004.txt
      character, save :: co2dat_file*26
!> external co2 global annual mean data tb: co2historicaldata_glob.txt
      character, save :: co2gbl_file*26
!> external co2 user defined data table: co2userdata.txt
      character, save :: co2usr_file*26
!> external co2 clim monthly cycle data tb: co2monthlycyc.txt
      character, save :: co2cyc_file*26
      data co2dat_file   / 'co2historicaldata_2004.txt' /   !year is run-time selected
      data co2gbl_file   / 'co2historicaldata_glob.txt' /
      data co2usr_file   / 'co2userdata.txt           ' /
      data co2cyc_file   / 'co2monthlycyc.txt         ' /

! ............................................. !
!>\name  2.4 For module radiation_clouds
! ............................................. !

!> cloud optical property scheme control flag
!!\n =0:use diagnostic cloud scheme for cloud cover and mean optical properties
!!\n =1:use prognostic cloud scheme for cloud cover and cloud properties
      integer, save :: icldflg = 1

!> cloud overlapping control flag for Radiation
!!\n =0:use random cloud overlapping method
!!\n =1:use maximum-random cloud overlapping method
!!\n =2:use maximum cloud overlapping method
!!\n =3:use decorrelation length overlapping method
!!\n =4:use exponential overlapping method
!!\n =5:use exponential-random overlapping method
!!\n Opr GFS/CFS=1; see IOVR in run scripts
      integer, save :: iovr  = 1
!!\n Decorrelation length type for iovr = 4 or 5
!!\n =0:use constant decorrelation length defined by decorr_con (in module physcons)
!!\n =1:use day-of-year and latitude-varying decorrelation length
      integer, save :: idcor   = 1

!> sub-column cloud approx flag in SW radiation
!!\n =0:no McICA approximation in SW radiation
!!\n =1:use McICA with precribed permutation seeds (test mode)
!!\n =2:use McICA with randomly generated permutation seeds
!!\n Opr GFS/CFS=2; see ISUBC_SW in run scripts
      integer, save :: isubcsw = 0
!> sub-column cloud approx flag in LW radiation
!!\n =0:no McICA approximation in LW radiation
!!\n =1:use McICA with prescribed permutation seeds (test mode)
!!\n =2:use McICA with randomly generatedo
!!\n Opr GFS/CFS=2; see ISUBC_LW in run scripts
      integer, save :: isubclw = 0

!> eliminating CRICK control flag
      logical, save :: lcrick  =.false.
!> in-cld condensate control flag
      logical, save :: lcnorm  =.false.
!> precip effect on radiation flag (Ferrier microphysics)
      logical, save :: lnoprec =.false.
!> shallow convetion flag
      logical, save :: lsashal =.false.

! ............................................. !
!> \name  2.6 general purpose
! ............................................. !

!> vertical profile indexing flag
      integer, save :: ivflip  = 1

!> initial permutaion seed for mcica radiation
      integer, save :: ipsd0   = 0
      integer, save :: ipsdlim = 1e8
!
!...................................!
      end module physparam          !
!===================================!

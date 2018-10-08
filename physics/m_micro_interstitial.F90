!> \file m_micro_interstitial.F90
!! This file contains subroutines that prepare data for and from the Morrison-Gettelman microphysics scheme
!! as part of the GFS physics suite.
module m_micro_pre
contains

! \brief Brief description of the subroutine
!
!> \section arg_table_m_micro_pre_init Argument Table
!!
subroutine m_micro_pre_init()
end subroutine m_micro_pre_init

! \brief Brief description of the subroutine
!!
!! \section arg_table_m_micro_pre_run Argument Table
!! | local_name     | standard_name                                         | long_name                                                                                   | units   | rank | type       | kind      | intent | optional |
!! |----------------|-------------------------------------------------------|---------------------------------------------------------------------------------------------|---------|------|------------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                | horizontal loop extent                                                                      | count   |    0 | integer    |           | in     | F        |
!! | levs           | vertical_dimension                                    | number of vertical layers                                                                   | count   |    0 | integer    |           | in     | F        |
!! | do_shoc        | flag_for_shoc                                         | flag for SHOC                                                                               | flag    |    0 | logical    |           | in     | F        |
!! | imfdeepcnv     | flag_for_mass_flux_deep_convection_scheme             | flag for mass-flux deep convection scheme                                                   | flag    |    0 | integer    |           | in     | F        |
!! | imfshalcnv     | flag_for_mass_flux_shallow_convection_scheme          | flag for mass-flux shallow convection scheme                                                | flag    |    0 | integer    |           | in     | F        |
!! | gq0_ice        | ice_water_mixing_ratio_updated_by_physics             | moist (dry+vapor, no condensates) mixing ratio of ice water updated by physics              | kg kg-1 |    2 | real       | kind_phys | in     | F        |
!! | gq0_water      | cloud_condensed_water_mixing_ratio_updated_by_physics | moist (dry+vapor, no condensates) mixing ratio of cloud condensed water updated by physics  | kg kg-1 |    2 | real       | kind_phys | in     | F        |
!! | cld_shoc       | subgrid_scale_cloud_fraction_from_shoc                | subgrid-scale cloud fraction from the SHOC scheme                                           | frac    |    2 | real       | kind_phys | in     | F        |
!! | cnvc           | convective_cloud_cover                                | convective cloud cover                                                                      | frac    |    2 | real       | kind_phys | in     | F        |
!! | cnvw           | convective_cloud_water_mixing_ratio                   | moist convective cloud water mixing ratio                                                   | kg kg-1 |    2 | real       | kind_phys | in     | F        |
!! | tcr            | cloud_phase_transition_threshold_temperature          | threshold temperature below which cloud starts to freeze                                    | K       |    0 | real       | kind_phys | in     | F        |
!! | tcrf           | cloud_phase_transition_denominator                    | denominator in cloud phase transition = 1/(tcr-tf)                                          | K-1     |    0 | real       | kind_phys | in     | F        |
!! | gt0            | air_temperature_updated_by_physics                    | temperature updated by physics                                                              | K       |    2 | real       | kind_phys | in     | F        |
!! | cld_frc_MG     | cloud_fraction_for_MG                                 | cloud fraction used by Morrison-Gettelman MP                                                | frac    |    2 | real       | kind_phys | inout  | F        |
!! | qlcn           | mass_fraction_of_convective_cloud_liquid_water        | mass fraction of convective cloud liquid water                                              | kg kg-1 |    2 | real       | kind_phys | inout  | F        |
!! | qicn           | mass_fraction_of_convective_cloud_ice                 | mass fraction of convective cloud ice water                                                 | kg kg-1 |    2 | real       | kind_phys | inout  | F        |
!! | cf_upi         | convective_cloud_fraction_for_microphysics            | convective cloud fraction for microphysics                                                  | frac    |    2 | real       | kind_phys | inout  | F        |
!! | clw_water      | cloud_liquid_water_mixing_ratio                       | moist cloud water mixing ratio                                                              | kg kg-1 |    2 | real       | kind_phys | out    | F        |
!! | clw_ice        | cloud_ice_mixing_ratio                                | moist cloud ice mixing ratio                                                                | kg kg-1 |    2 | real       | kind_phys | out    | F        |
!! | errmsg         | ccpp_error_message                                    | error message for error handling in CCPP                                                    | none    |    0 | character  | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                       | error flag for error handling in CCPP                                                       | flag    |    0 | integer    |           | out    | F        |
!!
subroutine m_micro_pre_run (im, levs, do_shoc, imfdeepcnv, imfshalcnv, gq0_ice, gq0_water, cld_shoc, cnvc, cnvw, &
  tcr, tcrf, gt0, cld_frc_MG, qlcn, qicn, cf_upi, clw_water, clw_ice, errmsg, errflg )

use machine, only : kind_phys
implicit none

integer, intent(in) :: im, levs, imfdeepcnv, imfshalcnv
logical, intent(in) :: do_shoc
real(kind=kind_phys), intent(in) :: tcr, tcrf

real(kind=kind_phys), intent(in) ::                               &
    gq0_ice(:,:), gq0_water(:,:), cld_shoc(:,:), cnvc(:,:),       &
    cnvw(:,:), gt0(:,:)

real(kind=kind_phys), intent(inout) ::                            &
    cld_frc_MG(:,:), cf_upi(:,:), qlcn(:,:), qicn(:,:)

real(kind=kind_phys), intent(out) :: clw_ice(:,:), clw_water(:,:)

character(len=*), intent(out) :: errmsg
integer,          intent(out) :: errflg

integer :: i, k
real(kind=kind_phys) :: tem

! Initialize CCPP error handling variables
errmsg = ''
errflg = 0

!       Acheng used clw here for other code to run smoothly and minimum change
!       to make the code work. However, the nc and clw should be treated
!       in other procceses too.  August 28/2015; Hope that can be done next
!       year. I believe this will make the physical interaction more reasonable
!       Anning 12/5/2015 changed ntcw hold liquid only
do k=1,levs
  do i=1,im
    clw_ice(i,k) = gq0_ice(i,k)             ! ice
    clw_water(i,k) = gq0_water(i,k)             ! water
  enddo
enddo


if (do_shoc) then
  do k=1,levs
    do i=1,im
      cld_frc_MG(i,k) = cld_shoc(i,k) ! clouds from shoc
    enddo
  enddo
else if ((imfdeepcnv >= 0) .or. (imfshalcnv > 0)) then
  do k=1,levs
    do i=1,im
      cld_frc_MG(i,k) = max(0.0, min(1.0,cld_frc_MG(i,k)+cnvc(i,k)))
                                             ! clouds from t-dt and cnvc
      tem = cnvw(i,k)* max(0.0, MIN(1.0, (TCR-gt0(i,k))*TCRF))
      qlcn(i,k)   = qlcn(i,k)   +  cnvw(i,k) - tem
      qicn(i,k)   = qicn(i,k)   + tem
      cf_upi(i,k) = cf_upi(i,k) + cnvc(i,k)
    enddo
  enddo
end if


end subroutine m_micro_pre_run

! \brief Brief description of the subroutine
!
!> \section arg_table_m_micro_pre_finalize Argument Table
!!
subroutine m_micro_pre_finalize ()
end subroutine m_micro_pre_finalize

end module m_micro_pre

!> This module contains the CCPP-compliant MG microphysics
!! post intersititial codes.
      module m_micro_post

      contains

! \brief Brief description of the subroutine
!
!> \section arg_table_m_micro_post_init Argument Table
!!
      subroutine m_micro_post_init()
      end subroutine m_micro_post_init

! \brief Brief description of the subroutine
!!
!! \section arg_table_m_micro_post_run Argument Table
!! | local_name     | standard_name                                                            | long_name                                                                        | units   | rank | type      | kind      | intent | optional |
!! |----------------|--------------------------------------------------------------------------|----------------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                   | horizontal loop extent                                                           | count   |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                       | number of vertical layers                                                        | count   |    0 | integer   |           | in     | F        |
!! | fprcp          | number_of_frozen_precipitation_species                                   | number of frozen precipitation species                                           | count   |    0 | integer   |           | in     | F        |
!! | mg3_as_mg2     | flag_mg3_as_mg2                                                          | flag for controlling prep for Morrison-Gettelman microphysics                    | flag    |    0 | logical   |           | in     | F        |
!! | qrn            | rain_water_mixing_ratio_updated_by_physics                               | moist (dry+vapor, no condensates) mixing ratio of rain water updated by physics  | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | qsnw           | snow_water_mixing_ratio_updated_by_physics                               | moist (dry+vapor, no condensates) mixing ratio of snow water updated by physics  | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | qgl            | graupel_mixing_ratio_updated_by_physics                                  | moist (dry+vapor, no condensates) mixing ratio of graupel updated by physics     | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                       | error message for error handling in CCPP                                         | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                          | error flag for error handling in CCPP                                            | flag    |    0 | integer   |           | out    | F        |
!!
      subroutine m_micro_post_run(                                         &
     &  im, levs, fprcp, mg3_as_mg2, qrn, qsnw, qgl, errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, fprcp
      logical, intent(in) :: mg3_as_mg2

      real(kind=kind_phys), intent(inout) ::                                 &
&        qrn(:,:), qsnw(:,:), qgl(:,:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k

      real(kind=kind_phys), parameter :: qsmall  = 1.0e-20

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!     do k=1,levs
!     write(1000+me,*)' maxwatnca=',maxval(Stateout%gq0(1:im,k,ntlnc)),' k=',k,' kdt=',kdt
!     enddo
!     write(1000+me,*)' at latitude = ',lat
!     tx1 = 1000.0
!     call moist_bud(im,ix,ix,levs,me,kdt,con_g,tx1,del,rain1
!    &,                    txa, clw(1,1,2), clw(1,1,1)
!    &,           gq0(1,1,1),gq0(1,1,ntcw),gq0(1,1,ntcw+1),' m_micro  ')

!       if (lprnt) write(0,*) ' rain1=',rain1(ipr)*86400.0, &
!    &' rainc=',diag%rainc(ipr)*86400.0                        &
!    &,' cn_prc=',cn_prc(ipr),' cn_snr=',cn_snr(ipr)
!       if(lprnt) write(0,*) ' aftgt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*) ' aftlsgq0=',stateout%gq0(ipr,:,1),' kdt=',kdt
!       if (lprnt) write(0,*)' clw1aft=',stateout%gq0(ipr,:,ntiw),' kdt=',kdt
!       if (ntgl > 0 .and. lprnt)  &
!                  write(0,*)' cgw1aft=',stateout%gq0(ipr,:,ntgl),' kdt=',kdt
!       if (lprnt) write(0,*)' cloudsm=',tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clw2aft=',stateout%gq0(ipr,:,ntcw),' kdt=',kdt
!       if (lprnt) write(0,*)' qrna=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwa=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglba',qgl(ipr,:),' kdt=',kdt

      if (abs(fprcp) == 1 .or. mg3_as_mg2) then
        do k=1,levs
          do i=1,im
            if (abs(qrn(i,k))  < qsmall) qrn(i,k)  = 0.0
            if (abs(qsnw(i,k)) < qsmall) qsnw(i,k) = 0.0
          end do
        end do
      else if (fprcp > 1) then
        do k=1,levs
          do i=1,im
            if (abs(qrn(i,k))  < qsmall) qrn(i,k)  = 0.0
            if (abs(qsnw(i,k)) < qsmall) qsnw(i,k) = 0.0
            if (abs(qgl(i,k))  < qsmall) qgl(i,k)  = 0.0
          end do
        end do
      end if

!       if (lprnt) write(0,*)' cloudsm=',tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clw2aft=',stateout%gq0(ipr,:,ntcw),' kdt=',kdt
!       if (lprnt) write(0,*)' qrna=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwa=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglba',qgl(ipr,:),' kdt=',kdt
!


      end subroutine m_micro_post_run

! \brief Brief description of the subroutine
!
!> \section arg_table_m_micro_post_finalize Argument Table
!!
      subroutine m_micro_post_finalize()
      end subroutine m_micro_post_finalize

      end module m_micro_post

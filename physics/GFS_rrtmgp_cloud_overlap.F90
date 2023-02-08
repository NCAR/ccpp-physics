!> \file GFS_rrtmgp_cloud_overlap.F90
!! 
!> \defgroup GFS_rrtmgp_cloud_overlap GFS_rrtmgp_cloud_overlap.F90
!!
!! \brief This module contains EMC's interface to the different assumptions of vertical cloud 
!! structuce, cloud overlap, used by McICA for cloud sampling in the RRTMGP longwave
!! and shortwave schemes.
!!
module GFS_rrtmgp_cloud_overlap
  use machine,      only: kind_phys
  use radiation_tools,   only: check_error_msg
  use module_radiation_cloud_overlap, only: cmp_dcorr_lgth, get_alpha_exper

  public GFS_rrtmgp_cloud_overlap_run

contains  

!>\defgroup gfs_rrtmgp_cloud_overlap_mod GFS RRTMGP Cloud Overlap Module
!! \section arg_table_GFS_rrtmgp_cloud_overlap_run
!! \htmlinclude GFS_rrtmgp_cloud_overlap_run.html
!!
!> \ingroup GFS_rrtmgp_cloud_overlap
!!
!! This is identical (shares common-code) to RRTMG. The motivation for RRTMGP to have
!! its own scheme is both organizational and philosophical*.
!!
!! *The number of "clouds" being produced by the model physics is often greater than one.
!! rte-rrtmgp can accomodate multiple cloud-types. This module preservers this enhancement
!! in the EMCs coupling to the RRTMGP scheme.
!!
!! \section GFS_rrtmgp_cloud_overlap_run
  subroutine GFS_rrtmgp_cloud_overlap_run(nCol, nLev, yearlen, doSWrad, doLWrad,         &
       julian, lat, p_lev, p_lay, tv_lay, deltaZc, con_pi, con_g, con_rd, con_epsq,      &
       dcorr_con, idcor, iovr, iovr_dcorr, iovr_exp, iovr_exprand, idcor_con,            &
       idcor_hogan, idcor_oreopoulos, cld_frac, cld_cnv_frac, iovr_convcld, top_at_1,    &
       imfdeepcnv, imfdeepcnv_gf, imfdeepcnv_samf, de_lgth, cloud_overlap_param,         &
       cnv_cloud_overlap_param, precip_overlap_param, errmsg, errflg)
    implicit none
    
    ! Inputs   
    integer, intent(in)     :: &
         nCol,                 & ! Number of horizontal grid points
         nLev,                 & ! Number of vertical layers
         yearlen,              & ! Length of current year (365/366) WTF?
         imfdeepcnv,           & !
         imfdeepcnv_gf,        & !
         imfdeepcnv_samf,      & !
         iovr,                 & ! Choice of cloud-overlap method
         iovr_convcld,         & ! Choice of convective cloud-overlap method
         iovr_dcorr,           & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,             & ! Flag for exponential cloud overlap method
         iovr_exprand,         & ! Flag for exponential-random cloud overlap method
         idcor,                & ! Choice of method for decorrelation length computation
         idcor_con,            & ! Flag for decorrelation-length. Use constant value
         idcor_hogan,          & ! Flag for decorrelation-length. (https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.647)
         idcor_oreopoulos        ! Flag for decorrelation-length. (10.5194/acp-12-9097-2012) 
    logical, intent(in)     :: &
         top_at_1,             & ! Vertical ordering flag
    	 doSWrad,              & ! Call SW radiation?
    	 doLWrad                 ! Call LW radiation
    real(kind_phys), intent(in) :: &
         julian,               & ! Julian day 
         con_pi,               & ! Physical constant: pi
         con_g,                & ! Physical constant: gravitational constant
         con_rd,               & ! Physical constant: gas-constant for dry air
         con_epsq,             & ! Physical constant: Minimum value for specific humidity
         dcorr_con               ! Decorrelation-length (used if idcor = idcor_con)
    real(kind_phys), dimension(:), intent(in) :: &
         lat                     ! Latitude             
    real(kind_phys), dimension(:,:), intent(in) :: &         
         tv_lay,               & ! Virtual temperature (K)
         p_lay,                & ! Pressure at model-layers (Pa)
         cld_frac,             & ! Total cloud fraction
         cld_cnv_frac            ! Convective cloud-fraction
    real(kind_phys), dimension(:,:), intent(in) :: &         
         p_lev,                & ! Pressure at model-level interfaces (Pa)
         deltaZc                 ! Layer thickness (from layer-centers)(m)
    
    ! Outputs     
    real(kind_phys), dimension(:),intent(out) :: &
         de_lgth                   ! Decorrelation length
    real(kind_phys), dimension(:,:),intent(out) :: &
         cloud_overlap_param,    & ! Cloud-overlap parameter
         cnv_cloud_overlap_param,& ! Convective cloud-overlap parameter
         precip_overlap_param      ! Precipitation overlap parameter
    character(len=*), intent(out) :: &
         errmsg                    ! Error message
    integer, intent(out) :: &
         errflg                    ! Error flag
    
    ! Local variables
    integer :: iCol,iLay

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. (doSWrad .or. doLWrad)) return

    !
    ! Cloud decorrelation length
    !
    de_lgth(:) = 0.
    if (idcor == idcor_hogan) then
       call cmp_dcorr_lgth(nCol, lat, con_pi, de_lgth)
    endif
    if (idcor == idcor_oreopoulos) then
       call cmp_dcorr_lgth(nCol, lat*(180._kind_phys/con_pi), julian, yearlen, de_lgth)
    endif
    if (idcor == idcor_con) then
       de_lgth(:) = dcorr_con
    endif

    !
    ! Cloud overlap parameter
    !
    if (iovr == iovr_dcorr .or. iovr == iovr_exp .or. iovr == iovr_exprand) then
       call get_alpha_exper(nCol, nLev, iovr, iovr_exprand, deltaZc*0.001, de_lgth, cld_frac, cloud_overlap_param)
    else
       cloud_overlap_param(:,:) = 0.
    endif

    !
    ! Convective cloud overlap parameter
    !
    if (imfdeepcnv == imfdeepcnv_samf .or. imfdeepcnv == imfdeepcnv_gf) then
       if (iovr_convcld == iovr_dcorr .or. iovr_convcld == iovr_exp .or. iovr_convcld == iovr_exprand) then
          call get_alpha_exper(nCol, nLev, iovr_convcld, iovr_exprand, deltaZc*0.001, de_lgth, cld_cnv_frac, cnv_cloud_overlap_param)
       else
          cnv_cloud_overlap_param(:,:) = 0.
       endif
    endif

    ! 
    ! Compute precipitation overlap parameter (Hack. Using same as cloud for now)
    !
    precip_overlap_param = cloud_overlap_param    
    
  end subroutine GFS_rrtmgp_cloud_overlap_run
end module GFS_rrtmgp_cloud_overlap

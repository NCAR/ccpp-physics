module rrtmgp_lw_cloud_sampling
  use machine,                  only: kind_phys
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use mo_optical_props,         only: ty_optical_props_2str
  use rrtmgp_sampling,          only: sampled_mask, draw_samples
  use mersenne_twister,         only: random_setseed, random_number, random_stat  
  use radiation_tools,          only: check_error_msg
  use rrtmgp_lw_gas_optics,     only: lw_gas_props
  use netcdf

  implicit none

contains

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_lw_cloud_sampling_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_sampling_run
!! \htmlinclude rrtmgp_lw_cloud_sampling_run.html
!!
  subroutine rrtmgp_lw_cloud_sampling_run(doLWrad, nCol, nLev, icseed_lw, iovr,iovr_convcld,&
       iovr_max, iovr_maxrand, iovr_rand, iovr_dcorr, iovr_exp, iovr_exprand, isubc_lw,     &
       cld_frac, precip_frac, cloud_overlap_param, precip_overlap_param, cld_cnv_frac,      &
       cnv_cloud_overlap_param, imfdeepcnv, imfdeepcnv_gf, imfdeepcnv_samf,                 &
       lw_optical_props_cloudsByBand, lw_optical_props_cnvcloudsByBand,                     &
       lw_optical_props_precipByBand, lw_optical_props_clouds, lw_optical_props_cnvclouds,  &
       lw_optical_props_precip, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad                             ! Logical flag for shortwave radiation call
    integer, intent(in) :: &
         nCol,                             & ! Number of horizontal gridpoints
         nLev,                             & ! Number of vertical layers
         imfdeepcnv,                       & !
         imfdeepcnv_gf,                    & !
         imfdeepcnv_samf,                  & !
         iovr,                             & ! Choice of cloud-overlap method
         iovr_convcld,                     & ! Choice of convective cloud-overlap
         iovr_max,                         & ! Flag for maximum cloud overlap method
         iovr_maxrand,                     & ! Flag for maximum-random cloud overlap method
         iovr_rand,                        & ! Flag for random cloud overlap method
         iovr_dcorr,                       & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,                         & ! Flag for exponential cloud overlap method
         iovr_exprand,                     & ! Flag for exponential-random cloud overlap method
         isubc_lw 
    integer,intent(in),dimension(:) :: &
         icseed_lw                           ! auxiliary special cloud related array when module 
                                             ! variable isubc_lw=2, it provides permutation seed 
                                             ! for each column profile that are used for generating 
                                             ! random numbers. when isubc_lw /=2, it will not be used.
    real(kind_phys), dimension(:,:),intent(in) :: &
         cld_frac,                         & ! Total cloud fraction by layer
         cld_cnv_frac,                     & ! Convective cloud fraction by layer 
         precip_frac,                      & ! Precipitation fraction by layer
         cloud_overlap_param,              & ! Cloud overlap parameter
         cnv_cloud_overlap_param,          & ! Convective cloud overlap parameter 
         precip_overlap_param                ! Precipitation overlap parameter 
    type(ty_optical_props_2str),intent(in) :: &
         lw_optical_props_cloudsByBand,    & ! RRTMGP DDT: Longwave optical properties in each band (clouds)
         lw_optical_props_cnvcloudsByBand, & ! RRTMGP DDT: Longwave optical properties in each band (convective cloud) 
         lw_optical_props_precipByBand       ! RRTMGP DDT: Longwave optical properties in each band (precipitation)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                         ! CCPP error message
    integer,          intent(out) :: &
         errflg                         ! CCPP error code
    type(ty_optical_props_2str),intent(inout) :: &
         lw_optical_props_clouds,     & ! RRTMGP DDT: Shortwave optical properties by spectral point (clouds)
         lw_optical_props_cnvclouds,  & ! RRTMGP DDT: Shortwave optical properties by spectral point (convective cloud)
         lw_optical_props_precip        ! RRTMGP DDT: Shortwave optical properties by spectral point (precipitation)

    ! Local variables
    integer :: iCol, iLay, iBand
    integer,dimension(ncol) :: ipseed_lw
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(lw_gas_props%get_ngpt(),nLev,ncol) :: rng3D,rng3D2
    real(kind_phys), dimension(lw_gas_props%get_ngpt()*nLev) :: rng2D
    real(kind_phys), dimension(lw_gas_props%get_ngpt()) :: rng1D
    logical, dimension(ncol,nLev,lw_gas_props%get_ngpt()) :: maskMCICA

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. doLWrad) return
    
    ! ####################################################################################    
    ! First sample the clouds...
    ! ####################################################################################
    lw_optical_props_clouds%band2gpt      = lw_gas_props%get_band_lims_gpoint()
    lw_optical_props_clouds%band_lims_wvn = lw_gas_props%get_band_lims_wavenumber()
    do iBand=1,lw_gas_props%get_nband()
       lw_optical_props_clouds%gpt2band(lw_optical_props_clouds%band2gpt(1,iBand):lw_optical_props_clouds%band2gpt(2,iBand)) = iBand
    end do
    
    ! Change random number seed value for each radiation invocation (isubc_lw =1 or 2).
    if(isubc_lw == 1) then      ! advance prescribed permutation seed
       do iCol = 1, ncol
          ipseed_lw(iCol) = lw_gas_props%get_ngpt() + iCol
       enddo
    elseif (isubc_lw == 2) then ! use input array of permutaion seeds
       do iCol = 1, ncol
          ipseed_lw(iCol) = icseed_lw(iCol)
       enddo
    endif
    
    ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
    ! and layers. ([nGpts,nLev,nColumn]-> [nGpts*nLev]*nColumn)
    do iCol=1,ncol
       call random_setseed(ipseed_lw(icol),rng_stat)
       ! Use same rng for each layer
       if (iovr == iovr_max) then
          call random_number(rng1D,rng_stat)
          do iLay=1,nLev
             rng3D(:,iLay,iCol) = rng1D
          enddo
       else
          do iLay=1,nLev
             call random_number(rng1D,rng_stat)
             rng3D(:,iLay,iCol) = rng1D
          enddo
       endif
    enddo

    ! Cloud-overlap.
    ! Maximum-random, random or maximum.
    if (iovr == iovr_maxrand .or. iovr == iovr_rand .or. iovr == iovr_max) then
       call sampled_mask(rng3D, cld_frac, maskMCICA) 
    endif
	!  Exponential decorrelation length overlap
    if (iovr == iovr_dcorr) then
       ! Generate second RNG
       do iCol=1,ncol
          call random_setseed(ipseed_lw(icol),rng_stat)
          call random_number(rng2D,rng_stat)
          rng3D2(:,:,iCol) = reshape(source = rng2D,shape=[lw_gas_props%get_ngpt(),nLev])
       enddo
       call sampled_mask(rng3D, cld_frac, maskMCICA,                    &
                         overlap_param = cloud_overlap_param(:,1:nLev-1),  &
                         randoms2      = rng3D2)
    endif
    ! Exponential or Exponential-random
    if (iovr == iovr_exp .or. iovr == iovr_exprand) then
       call sampled_mask(rng3D, cld_frac, maskMCICA,  &
                         overlap_param = cloud_overlap_param(:,1:nLev-1))    
    endif

    !
    ! Sampling. Map band optical depth to each g-point using McICA
    !
    call check_error_msg('rrtmgp_lw_cloud_sampling_run_draw_samples',&
         draw_samples(maskMCICA, .true.,                             &
                      lw_optical_props_cloudsByBand,                 &
                      lw_optical_props_clouds))

  end subroutine rrtmgp_lw_cloud_sampling_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_lw_cloud_sampling_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_lw_cloud_sampling_finalize()
  end subroutine rrtmgp_lw_cloud_sampling_finalize 

end module rrtmgp_lw_cloud_sampling

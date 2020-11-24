module rrtmgp_sw_cloud_sampling
  use machine,                  only: kind_phys
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use mo_optical_props,         only: ty_optical_props_2str
  use rrtmgp_sampling,          only: sampled_mask, draw_samples
  use mersenne_twister,         only: random_setseed, random_number, random_stat  
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

  implicit none

contains
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_cloud_sampling_init()
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_sampling_init
!! \htmlinclude rrtmgp_sw_cloud_sampling.html
!!
  subroutine rrtmgp_sw_cloud_sampling_init(sw_gas_props, ipsdsw0, errmsg, errflg)
    ! Inputs
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props ! RRTMGP DDT: K-distribution data
    ! Outputs
    integer, intent(out) :: &
         ipsdsw0      ! Initial permutation seed for McICA
    character(len=*), intent(out) :: &
         errmsg       ! Error message
    integer,          intent(out) :: &
         errflg       ! Error flag

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdsw0 = sw_gas_props%get_ngpt()
    
  end subroutine rrtmgp_sw_cloud_sampling_init

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_sampling_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_sampling_run
!! \htmlinclude rrtmgp_sw_cloud_sampling.html
!!
  subroutine rrtmgp_sw_cloud_sampling_run(doSWrad, nCol, nDay, nLev, ipsdsw0, idxday, iovr, &
       iovr_max, iovr_maxrand, iovr_rand, iovr_dcorr, iovr_exp, iovr_exprand, isubc_sw,     &
       icseed_sw, cld_frac, precip_frac, cloud_overlap_param, precip_overlap_param,         &
       sw_gas_props, sw_optical_props_cloudsByBand, sw_optical_props_precipByBand,          &
       sw_optical_props_clouds, sw_optical_props_precip, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doSWrad                            ! Logical flag for shortwave radiation call
    integer, intent(in) :: &
         nCol,                            & ! Number of horizontal gridpoints
         nDay,                            & ! Number of daylit points.
         nLev,                            & ! Number of vertical layers
         ipsdsw0,                         & ! Initial permutation seed for McICA
         iovr,                            & ! Choice of cloud-overlap method                                                                                                                 
         iovr_max,                        & ! Flag for maximum cloud overlap method                                                                                                          
         iovr_maxrand,                    & ! Flag for maximum-random cloud overlap method                                                                                                   
         iovr_rand,                       & ! Flag for random cloud overlap method                                                                                                           
         iovr_dcorr,                      & ! Flag for decorrelation-length cloud overlap method                                                                                             
         iovr_exp,                        & ! Flag for exponential cloud overlap method                                                                                                      
         iovr_exprand,                    & ! Flag for exponential-random cloud overlap method 
         isubc_sw
    integer,intent(in),dimension(ncol) :: &
         idxday                             ! Indices for daylit points.
    integer,intent(in),dimension(ncol) :: &
         icseed_sw                          ! auxiliary special cloud related array when module 
                                            ! variable isubc_sw=2, it provides permutation seed 
                                            ! for each column profile that are used for generating 
                                            ! random numbers. when isubc_sw /=2, it will not be used.
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         cld_frac,                        & ! Total cloud fraction by layer
         precip_frac                        ! Precipitation fraction by layer
    real(kind_phys), dimension(ncol,nLev), intent(in)  :: &
         cloud_overlap_param,             & ! Cloud overlap parameter
         precip_overlap_param               ! Precipitation overlap parameter
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props                       ! RRTMGP DDT: K-distribution data
    type(ty_optical_props_2str),intent(in) :: &
         sw_optical_props_cloudsByBand,   & ! RRTMGP DDT: Shortwave optical properties in each band (clouds)
         sw_optical_props_precipByBand      ! RRTMGP DDT: Shortwave optical properties in each band (precipitation)    

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                             ! Error message
    integer,          intent(out) :: &
         errflg                             ! Error flag
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_clouds,         & ! RRTMGP DDT: Shortwave optical properties at each spectral point (clouds) 
         sw_optical_props_precip            ! RRTMGP DDT: Shortwave optical properties at each spectral point (precipitation) 

    ! Local variables
    integer :: iday,iLay,iGpt
    integer,dimension(nday) :: ipseed_sw
    type(random_stat) :: rng_stat
    real(kind_phys) :: tauloc,asyloc,ssaloc
    real(kind_phys), dimension(sw_gas_props%get_ngpt(),nLev,nday) :: rng3D,rng3D2
    real(kind_phys), dimension(sw_gas_props%get_ngpt()*nLev) :: rng2D
    real(kind_phys), dimension(sw_gas_props%get_ngpt()) :: rng1D
    logical, dimension(nday,nLev,sw_gas_props%get_ngpt()) :: cldfracMCICA,precipfracSAMP

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. doSWrad) return
    if (nDay .gt. 0) then
       ! #################################################################################
       ! First sample the clouds...
       ! #################################################################################

       ! Allocate space RRTMGP DDTs [nday,nLev,nGpt]
       call check_error_msg('rrtmgp_sw_cloud_sampling_run', & 
            sw_optical_props_clouds%alloc_2str(nday, nLev, sw_gas_props))
       sw_optical_props_clouds%tau(:,:,:) = 0._kind_phys
       sw_optical_props_clouds%ssa(:,:,:) = 1._kind_phys
       sw_optical_props_clouds%g(:,:,:)   = 0._kind_phys
 
       ! Change random number seed value for each radiation invocation (isubc_sw =1 or 2).
       if(isubc_sw == 1) then      ! advance prescribed permutation seed
          do iday = 1, nday
             ipseed_sw(iday) = ipsdsw0 + iday
          enddo
       elseif (isubc_sw == 2) then ! use input array of permutaion seeds
          do iday = 1, nday
             ipseed_sw(iday) = icseed_sw(iday)
          enddo
       endif

       ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
       ! and layers. ([nGpts,nLev,nDayumn]-> [nGpts*nLev]*nDayumn)
       do iday=1,nday
          call random_setseed(ipseed_sw(iday),rng_stat)
          ! Use same rng for each layer
          if (iovr == iovr_max) then
             call random_number(rng1D,rng_stat)
             do iLay=1,nLev
                rng3D(:,iLay,iday) = rng1D
             enddo
          else
             do iLay=1,nLev
                call random_number(rng1D,rng_stat)
                rng3D(:,iLay,iday) = rng1D
             enddo
          endif
       enddo

       do iday=1,nday
          call random_setseed(ipseed_sw(iday),rng_stat)
          call random_number(rng2D,rng_stat)
          rng3D(:,:,iday) = reshape(source = rng2D,shape=[sw_gas_props%get_ngpt(),nLev])
       enddo

       ! Cloud overlap.
       ! Maximum-random, random, or maximum cloud overlap
       if (iovr == iovr_maxrand .or. iovr == iovr_max .or. iovr == iovr_rand) then
          call sampled_mask(rng3D, cld_frac(idxday(1:nDay),:), cldfracMCICA)  
       endif
       ! Decorrelation-length overlap
       if (iovr == iovr_dcorr) then
          do iday=1,nday
             call random_setseed(ipseed_sw(iday),rng_stat)
             call random_number(rng2D,rng_stat)
             rng3D2(:,:,iday) = reshape(source = rng2D,shape=[sw_gas_props%get_ngpt(),nLev])
          enddo
          call sampled_mask(rng3D, cld_frac(idxday(1:nDay),:), cldfracMCICA,             &
	                        overlap_param = cloud_overlap_param(idxday(1:nDay),1:nLev-1),&
	                        randoms2      = rng3D2)
       endif 
       ! Exponential or exponential-random cloud overlap
       if (iovr == iovr_exp .or. iovr == iovr_exprand) then
          call sampled_mask(rng3D, cld_frac(idxday(1:nDay),:), cldfracMCICA, &
                            overlap_param = cloud_overlap_param(idxday(1:nDay),1:nLev-1))
       endif

       !
       ! Sampling. Map band optical depth to each g-point using McICA
       !
       call check_error_msg('rrtmgp_sw_cloud_sampling_run_draw_samples', & 
            draw_samples(cldfracMCICA,                      &
                         sw_optical_props_cloudsByBand,     &
                         sw_optical_props_clouds))
         
       ! #################################################################################       
       ! Next sample precipitation (same as clouds for now)
       ! #################################################################################

       ! Allocate space RRTMGP DDTs [nday,nLev,nGpt]
       call check_error_msg('rrtmgp_sw_cloud_sampling_run', &
           sw_optical_props_precip%alloc_2str( nday, nLev, sw_gas_props))
 
       ! Change random number seed value for each radiation invocation (isubc_sw =1 or 2).
       if(isubc_sw == 1) then      ! advance prescribed permutation seed
          do iday = 1, nday
             ipseed_sw(iday) = ipsdsw0 + iday
          enddo
       elseif (isubc_sw == 2) then ! use input array of permutaion seeds
          do iday = 1, nday
             ipseed_sw(iday) = icseed_sw(iday)
          enddo
       endif

       ! No need to call RNG second time for now, just use the same seeds for precip as clouds.
       !! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
       !! and layers. ([nGpts,nLev,nDay]-> [nGpts*nLev]*nDay)
       !do iday=1,nday
       !   call random_setseed(ipseed_sw(iday),rng_stat)
       !   call random_number(rng1D,rng_stat)
       !   rng3D(:,:,iday) = reshape(source = rng1D,shape=[sw_gas_props%get_ngpt(),nLev])
       !enddo

       ! Precipitation overlap
       ! Maximum-random, random or maximum precipitation overlap
       if (iovr == iovr_maxrand .or. iovr == iovr_max .or. iovr == iovr_rand) then
          call sampled_mask(rng3D, precip_frac(idxday(1:nDay),:), precipfracSAMP)       
       endif
       ! Exponential decorrelation length overlap
       if (iovr == iovr_dcorr) then
          !! Generate second RNG
          !do iday=1,nday
          !   call random_setseed(ipseed_sw(iday),rng_stat)
          !   call random_number(rng1D,rng_stat)
          !   rng3D2(:,:,iday) = reshape(source = rng1D,shape=[sw_gas_props%get_ngpt(),nLev])
          !enddo
          call sampled_mask(rng3D, precip_frac(idxday(1:nDay),:), precipfracSAMP,         & 
                            overlap_param = precip_overlap_param(idxday(1:nDay),1:nLev-1),& 
                            randoms2 = rng3D2)
       endif
       if (iovr == iovr_exp .or. iovr == iovr_exprand) then
          call sampled_mask(rng3D, precip_frac(idxday(1:nDay),:),precipfracSAMP, &
                            overlap_param = precip_overlap_param(idxday(1:nDay),1:nLev-1))
       endif
 
       !
       ! Sampling. Map band optical depth to each g-point using McICA
       !
       call check_error_msg('rrtmgp_sw_precip_sampling_run_draw_samples', & 
            draw_samples(precipfracSAMP,                    &
                         sw_optical_props_precipByBand,     &
                         sw_optical_props_precip))                  
         
       ! #################################################################################        
       ! Just add precipitation optics to cloud-optics
       ! #################################################################################        
       do iGpt=1,sw_gas_props%get_ngpt()
          do iday=1,nDay
             do iLay=1,nLev
                tauloc = sw_optical_props_clouds%tau(iday,iLay,iGpt) + &
                         sw_optical_props_precip%tau(iday,iLay,iGpt)
                if (sw_optical_props_precip%tau(iday,iLay,iGpt) > 0) then
                   ssaloc = (sw_optical_props_clouds%tau(iday,iLay,iGpt)  * &
                             sw_optical_props_clouds%ssa(iday,iLay,iGpt)  + &
                             sw_optical_props_precip%tau(iday,iLay,iGpt)  * &
                             sw_optical_props_precip%ssa(iday,iLay,iGpt)) / &
                             tauloc
                   if (ssaloc > 0) then
                      asyloc = (sw_optical_props_clouds%tau(iday,iLay,iGpt) * &
                                sw_optical_props_clouds%ssa(iday,iLay,iGpt) * &
                                sw_optical_props_clouds%g(iday,iLay,iGpt)   + &
                                sw_optical_props_precip%tau(iday,iLay,iGpt) * &
                                sw_optical_props_precip%ssa(iday,iLay,iGpt) * &
                                sw_optical_props_precip%g(iday,iLay,iGpt))  / &
                                (tauloc*ssaloc)
                   else
                      tauloc = sw_optical_props_clouds%tau(iday,iLay,iGpt) 
                      ssaloc = sw_optical_props_clouds%ssa(iday,iLay,iGpt)
                      asyloc = sw_optical_props_clouds%g(iday,iLay,iGpt)            
                   endif
                   sw_optical_props_clouds%tau(iday,iLay,iGpt) = tauloc	
                   sw_optical_props_clouds%ssa(iday,iLay,iGpt) = ssaloc   
                   sw_optical_props_clouds%g(iday,iLay,iGpt)   = asyloc
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine rrtmgp_sw_cloud_sampling_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_sampling_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_sw_cloud_sampling_finalize()
  end subroutine rrtmgp_sw_cloud_sampling_finalize 

end module rrtmgp_sw_cloud_sampling

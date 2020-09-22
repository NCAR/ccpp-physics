module rrtmgp_lw_cloud_sampling
  use machine,                  only: kind_phys
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use physparam,                only: isubclw, iovrlw
  use mo_optical_props,         only: ty_optical_props_1scl
  use rrtmgp_sampling,          only: sampled_mask, draw_samples
  use mersenne_twister,         only: random_setseed, random_number, random_stat  
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

  implicit none

contains

  ! #########################################################################################
  ! SUBROUTINE mcica_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_sampling_init
!! \htmlinclude rrtmgp_lw_cloud_sampling_init.html
!!
  subroutine rrtmgp_lw_cloud_sampling_init(lw_gas_props, ipsdlw0, errmsg, errflg)
    ! Inputs
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props ! RRTMGP DDT: K-distribution data
    ! Outputs
    integer, intent(out) :: &
         ipsdlw0      ! Initial permutation seed for McICA
    character(len=*), intent(out) :: &
         errmsg       ! Error message
    integer,          intent(out) :: &
         errflg       ! Error flag

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdlw0 = lw_gas_props%get_ngpt()

  end subroutine rrtmgp_lw_cloud_sampling_init

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_lw_cloud_sampling_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_sampling_run
!! \htmlinclude rrtmgp_lw_cloud_sampling_run.html
!!
  subroutine rrtmgp_lw_cloud_sampling_run(doLWrad, nCol, nLev, ipsdlw0, icseed_lw,   &
       cld_frac, precip_frac, cloud_overlap_param, precip_overlap_param, lw_gas_props,      &
       lw_optical_props_cloudsByBand, lw_optical_props_precipByBand,                        &
       lw_optical_props_clouds, lw_optical_props_precip, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad                             ! Logical flag for shortwave radiation call
    integer, intent(in) :: &
         nCol,                             & ! Number of horizontal gridpoints
         nLev,                             & ! Number of vertical layers
         ipsdlw0                             ! Initial permutation seed for McICA
    integer,intent(in),dimension(ncol) :: &
         icseed_lw                           ! auxiliary special cloud related array when module 
                                             ! variable isubclw=2, it provides permutation seed 
                                             ! for each column profile that are used for generating 
                                             ! random numbers. when isubclw /=2, it will not be used.
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         cld_frac,                         & ! Total cloud fraction by layer
         precip_frac                         ! Precipitation fraction by layer
    real(kind_phys), dimension(ncol,nLev), intent(in)  :: &
         cloud_overlap_param,              & ! Cloud overlap parameter
         precip_overlap_param                ! Precipitation overlap parameter 
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props                        ! RRTMGP DDT: K-distribution data
    type(ty_optical_props_1scl),intent(in) :: &
         lw_optical_props_cloudsByBand,    & ! RRTMGP DDT: Longwave optical properties in each band (clouds)
         lw_optical_props_precipByBand       ! RRTMGP DDT: Longwave optical properties in each band (precipitation)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                         ! CCPP error message
    integer,          intent(out) :: &
         errflg                         ! CCPP error code
    type(ty_optical_props_1scl),intent(out) :: &
         lw_optical_props_clouds,     & ! RRTMGP DDT: Shortwave optical properties by spectral point (clouds)
         lw_optical_props_precip        ! RRTMGP DDT: Shortwave optical properties by spectral point (precipitation)

    ! Local variables
    integer :: iCol
    integer,dimension(ncol) :: ipseed_lw
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(lw_gas_props%get_ngpt(),nLev,ncol) :: rng3D,rng3D2
    real(kind_phys), dimension(lw_gas_props%get_ngpt()*nLev) :: rng1D
    logical, dimension(ncol,nLev,lw_gas_props%get_ngpt()) :: cldfracMCICA,precipfracSAMP

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! 
    if (iovrlw .ne. 1 .and. iovrlw .ne. 3 .and. iovrlw .ne. 4 .and. iovrlw .ne. 5) then
       errmsg = 'Cloud overlap assumption not supported.'
       errflg = 1
       call check_error_msg('rrtmgp_lw_cloud_sampling',errmsg)
       return
    endif
    
    if (.not. doLWrad) return
    
    ! ####################################################################################    
    ! First sample the clouds...
    ! ####################################################################################

    ! Allocate space RRTMGP DDTs [nCol,nLev,nGpt]
    call check_error_msg('rrtmgp_lw_cloud_sampling_run',&
         lw_optical_props_clouds%alloc_1scl(nCol, nLev, lw_gas_props))
    
    ! Change random number seed value for each radiation invocation (isubclw =1 or 2).
    if(isubclw == 1) then      ! advance prescribed permutation seed
       do iCol = 1, ncol
          ipseed_lw(iCol) = ipsdlw0 + iCol
       enddo
    elseif (isubclw == 2) then ! use input array of permutaion seeds
       do iCol = 1, ncol
          ipseed_lw(iCol) = icseed_lw(iCol)
       enddo
    endif
    
    ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
    ! and layers. ([nGpts,nLev,nColumn]-> [nGpts*nLev]*nColumn)
    do iCol=1,ncol
       call random_setseed(ipseed_lw(icol),rng_stat)
       call random_number(rng1D,rng_stat)
       rng3D(:,:,iCol) = reshape(source = rng1D,shape=[lw_gas_props%get_ngpt(),nLev])
    enddo

    ! Cloud-overlap.
    ! Maximum-random
    if (iovrlw == 1) then
       call sampled_mask(rng3D, cld_frac, cldfracMCICA) 
    endif
	!  Exponential decorrelation length overlap
    if (iovrlw == 3) then
       ! Generate second RNG
       do iCol=1,ncol
          call random_setseed(ipseed_lw(icol),rng_stat)
          call random_number(rng1D,rng_stat)
          rng3D2(:,:,iCol) = reshape(source = rng1D,shape=[lw_gas_props%get_ngpt(),nLev])
       enddo
       call sampled_mask(rng3D, cld_frac, cldfracMCICA,                    &
                         overlap_param = cloud_overlap_param(:,1:nLev-1),  &
                         randoms2      = rng3D2)
    endif
    ! Exponential or Exponential-random
    if (iovrlw == 4 .or. iovrlw == 5) then
       call sampled_mask(rng3D, cld_frac, cldfracMCICA,  &
                         overlap_param = cloud_overlap_param(:,1:nLev-1))    
    endif

    ! Sampling. Map band optical depth to each g-point using McICA
    call check_error_msg('rrtmgp_lw_cloud_sampling_run_draw_samples',&
         draw_samples(cldfracMCICA,                                  &
                      lw_optical_props_cloudsByBand,                 &
                      lw_optical_props_clouds))

    ! ####################################################################################
    ! Next sample the precipitation...
    ! ####################################################################################
    
    ! Allocate space RRTMGP DDTs [nCol,nLev,nGpt]
    call check_error_msg('rrtmgp_lw_cloud_sampling_run',&
         lw_optical_props_precip%alloc_1scl(nCol, nLev, lw_gas_props))
    
    ! Change random number seed value for each radiation invocation (isubclw =1 or 2).
    if(isubclw == 1) then      ! advance prescribed permutation seed
       do iCol = 1, ncol
          ipseed_lw(iCol) = ipsdlw0 + iCol
       enddo
    elseif (isubclw == 2) then ! use input array of permutaion seeds
       do iCol = 1, ncol
          ipseed_lw(iCol) = icseed_lw(iCol)
       enddo
    endif
    
    ! No need to call RNG second time for now, just use the same seeds for precip as clouds.
    !! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
    !! and layers. ([nGpts,nLev,nColumn]-> [nGpts*nLev]*nColumn)
    !do iCol=1,ncol
    !   call random_setseed(ipseed_lw(icol),rng_stat)
    !   call random_number(rng1D,rng_stat)
    !   rng3D(:,:,iCol) = reshape(source = rng1D,shape=[lw_gas_props%get_ngpt(),nLev])
    !enddo

    ! Precipitation overlap.
	! Maximum-random
    if (iovrlw == 1) then
        call sampled_mask(rng3D, precip_frac, precipfracSAMP)      
    endif 
	!  Exponential decorrelation length overlap
    if (iovrlw == 3) then
       ! No need to call RNG second time for now, just use the same seeds for precip as clouds.
       !! Generate second RNG
       !do iCol=1,ncol
       !   call random_setseed(ipseed_lw(icol),rng_stat)
       !   call random_number(rng1D,rng_stat)
       !   rng3D2(:,:,iCol) = reshape(source = rng1D,shape=[lw_gas_props%get_ngpt(),nLev])
       !enddo
       call sampled_mask(rng3D, precip_frac, precipfracSAMP,               &
                         overlap_param = precip_overlap_param(:,1:nLev-1), &
                         randoms2      = rng3D2)
    endif
    ! Exponential or Exponential-random
    if (iovrlw == 4 .or. iovrlw == 5) then
       call sampled_mask(rng3D, precip_frac, precipfracSAMP,               &
                         overlap_param = precip_overlap_param(:,1:nLev-1))
    endif
    
    
    ! Sampling. Map band optical depth to each g-point using McICA
    call check_error_msg('rrtmgp_lw_precip_sampling_run_draw_samples',&
         draw_samples(precipfracSAMP,                                 &
                      lw_optical_props_precipByBand,                  &
                      lw_optical_props_precip))
         
    ! ####################################################################################
    ! Just add precipitation optics to cloud-optics
    ! ####################################################################################
    lw_optical_props_clouds%tau = lw_optical_props_clouds%tau + lw_optical_props_precip%tau

  end subroutine rrtmgp_lw_cloud_sampling_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_lw_cloud_sampling_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_lw_cloud_sampling_finalize()
  end subroutine rrtmgp_lw_cloud_sampling_finalize 

end module rrtmgp_lw_cloud_sampling

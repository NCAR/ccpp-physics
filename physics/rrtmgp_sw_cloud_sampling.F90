module rrtmgp_sw_cloud_sampling
  use machine,                  only: kind_phys
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use physparam,                only: isubcsw, iovrsw
  use mo_optical_props,         only: ty_optical_props_2str
  use mo_cloud_sampling,        only: sampled_mask_max_ran, sampled_mask_exp_ran, draw_samples
  use mersenne_twister,         only: random_setseed, random_number, random_stat  
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

  implicit none

contains

  ! #########################################################################################
  ! SUBROUTINE mcica_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_sampling_init
!! \htmlinclude rrtmgp_sw_cloud_sampling.html
!!
  subroutine rrtmgp_sw_cloud_sampling_init(sw_gas_props, ipsdsw0)
    ! Inputs
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props ! RRTMGP DDT: K-distribution data
    ! Outputs
    integer, intent(out) :: &
         ipsdsw0      ! Initial permutation seed for McICA

    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdsw0 = sw_gas_props%get_ngpt()

  end subroutine rrtmgp_sw_cloud_sampling_init

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_sampling_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_sampling_run
!! \htmlinclude rrtmgp_sw_cloud_sampling.html
!!
  subroutine rrtmgp_sw_cloud_sampling_run(doSWrad, nCol, nDay, nLev, ipsdsw0, idxday,       &
       icseed_sw, cld_frac, sw_gas_props, sw_optical_props_cloudsByBand,                    &
       sw_optical_props_clouds, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doSWrad                        ! Logical flag for shortwave radiation call
    integer, intent(in) :: &
         nCol,                        & ! Number of horizontal gridpoints
         nDay,                        & ! Number of daylit points.
         nLev,                        & ! Number of vertical layers
         ipsdsw0                        ! Initial permutation seed for McICA
    integer,intent(in),dimension(ncol) :: &
         idxday                         ! Indices for daylit points.
    integer,intent(in),dimension(ncol) :: &
         icseed_sw                      ! auxiliary special cloud related array when module 
                                        ! variable isubcsw=2, it provides permutation seed 
                                        ! for each column profile that are used for generating 
                                        ! random numbers. when isubcsw /=2, it will not be used.
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         cld_frac                       ! Total cloud fraction by layer
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props                   ! RRTMGP DDT: K-distribution data
    type(ty_optical_props_2str),intent(in) :: &
         sw_optical_props_cloudsByBand  ! RRTMGP DDT: Shortwave optical properties (cloudy atmosphere) 

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                         ! Error message
    integer,          intent(out) :: &
         errflg                         ! Error code
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_clouds         ! RRTMGP DDT: Shortwave optical properties (cloudy atmosphere) 

    ! Local variables
    integer :: iCol
    integer,dimension(ncol) :: ipseed_sw
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(sw_gas_props%get_ngpt(),nLev,ncol) :: rng3D
    real(kind_phys), dimension(sw_gas_props%get_ngpt()*nLev) :: rng1D
    logical, dimension(ncol,nLev,sw_gas_props%get_ngpt()) :: cldfracMCICA
    real(kind_phys), dimension(ncol,nLev) :: cld_frac_noSamp

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doSWrad) return
    if (nDay .gt. 0) then
       
       ! Allocate space RRTMGP DDTs [nday,nLev,nGpt]
       call check_error_msg('rrtmgp_sw_cloud_sampling_run',sw_optical_props_clouds%alloc_2str(      &
            nday, nLev, sw_gas_props))
 
       ! Change random number seed value for each radiation invocation (isubcsw =1 or 2).
       if(isubcsw == 1) then      ! advance prescribed permutation seed
          do iCol = 1, ncol
             ipseed_sw(iCol) = ipsdsw0 + iCol
          enddo
       elseif (isubcsw == 2) then ! use input array of permutaion seeds
          do iCol = 1, ncol
             ipseed_sw(iCol) = icseed_sw(iCol)
          enddo
       endif

       ! Call McICA to generate subcolumns.
       ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
       ! and layers. ([nGpts,nLev,nColumn]-> [nGpts*nLev]*nColumn)
       do iCol=1,ncol
          call random_setseed(ipseed_sw(icol),rng_stat)
          call random_number(rng1D,rng_stat)
          rng3D(:,:,iCol) = reshape(source = rng1D,shape=[sw_gas_props%get_ngpt(),nLev])
       enddo
       
       ! Call McICA
       select case ( iovrsw )
          ! Maximumn-random 
       case(1)
          call check_error_msg('rrtmgp_sw_cloud_sampling_run',sampled_mask_max_ran(rng3D,cld_frac,cldfracMCICA))       
       end select
       
       ! Map band optical depth to each g-point using McICA
       call check_error_msg('rrtmgp_sw_cloud_sampling_run',draw_samples(&
            cldfracMCICA(idxday(1:nDay),:,:),sw_optical_props_cloudsByBand,sw_optical_props_clouds))
         
    endif

  end subroutine rrtmgp_sw_cloud_sampling_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_sampling_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_sw_cloud_sampling_finalize()
  end subroutine rrtmgp_sw_cloud_sampling_finalize 

end module rrtmgp_sw_cloud_sampling

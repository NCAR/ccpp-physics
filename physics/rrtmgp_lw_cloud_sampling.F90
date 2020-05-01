module rrtmgp_lw_cloud_sampling
  use machine,                  only: kind_phys
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use physparam,                only: isubclw, iovrlw
  use mo_optical_props,         only: ty_optical_props_1scl
  use mo_cloud_sampling,        only: sampled_mask_max_ran, sampled_mask_exp_ran, draw_samples
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
  subroutine rrtmgp_lw_cloud_sampling_init(lw_gas_props, ipsdlw0)
    ! Inputs
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props ! RRTMGP DDT: K-distribution data
    ! Outputs
    integer, intent(out) :: &
         ipsdlw0      ! Initial permutation seed for McICA

    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdlw0 = lw_gas_props%get_ngpt()

  end subroutine rrtmgp_lw_cloud_sampling_init

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_lw_cloud_sampling_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_sampling_run
!! \htmlinclude rrtmgp_lw_cloud_sampling_run.html
!!
  subroutine rrtmgp_lw_cloud_sampling_run(doLWrad, nCol, nLev, ipsdlw0, icseed_lw, cld_frac,&
       lw_gas_props, lw_optical_props_cloudsByBand, lw_optical_props_clouds, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad                        ! Logical flag for shortwave radiation call
    integer, intent(in) :: &
         nCol,                        & ! Number of horizontal gridpoints
         nLev,                        & ! Number of vertical layers
         ipsdlw0                        ! Initial permutation seed for McICA
    integer,intent(in),dimension(ncol) :: &
         icseed_lw                      ! auxiliary special cloud related array when module 
                                        ! variable isubclw=2, it provides permutation seed 
                                        ! for each column profile that are used for generating 
                                        ! random numbers. when isubclw /=2, it will not be used.
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         cld_frac                       ! Total cloud fraction by layer
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props                   ! RRTMGP DDT: K-distribution data
    type(ty_optical_props_1scl),intent(in) :: &
         lw_optical_props_cloudsByBand  ! RRTMGP DDT: Shortwave optical properties (cloudy atmosphere) 

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                         ! CCPP error message
    integer,          intent(out) :: &
         errflg                         ! CCPP error code
    type(ty_optical_props_1scl),intent(out) :: &
         lw_optical_props_clouds        ! RRTMGP DDT: Shortwave optical properties (cloudy atmosphere) 

    ! Local variables
    integer :: iCol
    integer,dimension(ncol) :: ipseed_lw
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(lw_gas_props%get_ngpt(),nLev,ncol) :: rng3D
    real(kind_phys), dimension(lw_gas_props%get_ngpt()*nLev) :: rng1D
    logical, dimension(ncol,nLev,lw_gas_props%get_ngpt()) :: cldfracMCICA
    real(kind_phys), dimension(ncol,nLev) :: cld_frac_noSamp

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return
       
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
    
    ! Call McICA to generate subcolumns.
    ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
    ! and layers. ([nGpts,nLev,nColumn]-> [nGpts*nLev]*nColumn)
    do iCol=1,ncol
       call random_setseed(ipseed_lw(icol),rng_stat)
       call random_number(rng1D,rng_stat)
       rng3D(:,:,iCol) = reshape(source = rng1D,shape=[lw_gas_props%get_ngpt(),nLev])
    enddo

    ! Call McICA
    select case ( iovrlw )
       ! Maximumn-random 
    case(1)
       call check_error_msg('rrtmgp_lw_cloud_sampling_run',sampled_mask_max_ran(rng3D,cld_frac,cldfracMCICA))       
    end select
    
    ! Map band optical depth to each g-point using McICA
    call check_error_msg('rrtmgp_lw_cloud_sampling_run',draw_samples(&
         cldfracMCICA,lw_optical_props_cloudsByBand,lw_optical_props_clouds))

  end subroutine rrtmgp_lw_cloud_sampling_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_lw_cloud_sampling_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_lw_cloud_sampling_finalize()
  end subroutine rrtmgp_lw_cloud_sampling_finalize 

end module rrtmgp_lw_cloud_sampling

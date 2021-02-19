! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2019,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! This module provides a simple implementation of sampling for the
!   Monte Carlo Independent Pixel Approximation (McICA, doi:10.1029/2002jd003322)
! Cloud optical properties, defined by band and assumed homogenous within each cell (column/layer),
!   are randomly sampled to preserve the mean cloud fraction and one of several possible overlap assumptions
! Users supply random numbers with order ngpt,nlay,ncol
!   These are only accessed if cloud_fraction(icol,ilay) > 0 so many values don't need to be filled in
!
! Adapted by Dustin Swales on 8/11/2020 for use in UFS (NOAA-PSL/CU-CIRES)
!
! -------------------------------------------------------------------------------------------------
module rrtmgp_sampling
  use mo_rte_kind,      only: wp, wl
  use mo_optical_props, only: ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  implicit none
  private
  public :: draw_samples, sampled_mask
contains
  ! -------------------------------------------------------------------------------------------------
  !
  ! Apply a T/F sampled cloud mask to cloud optical properties defined by band to produce
  !   McICA-sampled cloud optical properties
  !
  ! -------------------------------------------------------------------------------------------------
  function draw_samples(cloud_mask,do_twostream,clouds,clouds_sampled) result(error_msg)
	! Inputs
    logical, dimension(:,:,:),      intent(in   ) :: cloud_mask     ! Dimensions ncol,nlay,ngpt
    logical,                        intent(in   ) :: do_twostream   ! Do two-stream?
    class(ty_optical_props_arry),   intent(in   ) :: clouds         ! Defined by band
 
 	! Outputs
    class(ty_optical_props_arry),   intent(inout) :: clouds_sampled ! Defined by g-point
    character(len=128)                            :: error_msg

	! Local variables
    integer :: ncol,nlay,nbnd,ngpt
    integer :: imom
    
    error_msg = ""

    ! Array extents
    ncol = clouds%get_ncol()
    nlay = clouds%get_nlay()
    nbnd = clouds%get_nband()
    ngpt = clouds_sampled%get_ngpt()

    ! Optical depth assignment works for 1scl, 2str (also nstr)
    call apply_cloud_mask(ncol,nlay,nbnd,ngpt,clouds_sampled%get_band_lims_gpoint(),cloud_mask,clouds%tau,clouds_sampled%tau)
    !
    ! For 2-stream
    !
    select type(clouds)
    type is (ty_optical_props_2str)
      select type(clouds_sampled)
      type is (ty_optical_props_2str)
         if (do_twostream) then
            call apply_cloud_mask(ncol,nlay,nbnd,ngpt,clouds_sampled%get_band_lims_gpoint(),cloud_mask,clouds%ssa,clouds_sampled%ssa)
            call apply_cloud_mask(ncol,nlay,nbnd,ngpt,clouds_sampled%get_band_lims_gpoint(),cloud_mask,clouds%g,  clouds_sampled%g  )
         endif
      class default
          error_msg = "draw_samples: by-band and sampled cloud properties need to be the same variable type"
      end select
    end select
  end function draw_samples
  ! -------------------------------------------------------------------------------------------------
  !
  ! Generate a McICA-sampled cloud mask
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine sampled_mask(randoms, cloud_frac, cloud_mask, overlap_param, randoms2)
  	! Inputs
    real(wp), dimension(:,:,:),  intent(in )           :: randoms       ! ngpt,nlay,ncol
    real(wp), dimension(:,:),    intent(in )           :: cloud_frac    ! ncol,nlay

	! Outputs
    logical,  dimension(:,:,:),  intent(out)           :: cloud_mask    ! ncol,nlay,ngpt

	! Inputs (optional)
    real(wp), dimension(:,:),    intent(in ), optional :: overlap_param ! ncol,nlay-1
	real(wp), dimension(:,:,:),  intent(in ), optional :: randoms2      ! ngpt,nlay,ncol
    
    ! Local variables
    integer                              :: ncol, nlay, ngpt, icol, ilay, igpt
    integer                              :: cloud_lay_fst, cloud_lay_lst
    real(wp)                             :: rho
    real(wp), dimension(size(randoms,1)) :: local_rands
    logical,  dimension(size(randoms,2)) :: cloud_mask_layer
    logical                              :: l_use_overlap_param = .false.
    logical                              :: l_use_second_rng = .false.
    character(len=128)                   :: error_msg

	! Array dimensions
    ncol = size(randoms, 3)
    nlay = size(randoms, 2)
    ngpt = size(randoms, 1)

    ! Using cloud-overlap parameter (alpha)?
    if (present(overlap_param)) l_use_overlap_param = .true.
    
    ! Using a second RNG?
    if (present(randoms2)) l_use_second_rng = .true.
    
    ! Construct the cloud mask for each column
    do icol = 1, ncol
      cloud_mask_layer(1:nlay) = cloud_frac(icol,1:nlay) > 0._wp
      if(.not. any(cloud_mask_layer)) then
        cloud_mask(icol,1:nlay,1:ngpt) = .false.
        cycle
      end if
      cloud_lay_fst = findloc(cloud_mask_layer, .true., dim=1)
      cloud_lay_lst = findloc(cloud_mask_layer, .true., dim=1, back = .true.)
      cloud_mask(icol,1:cloud_lay_fst,1:ngpt) = .false.

      ilay = cloud_lay_fst
      local_rands(1:ngpt) = randoms(1:ngpt,cloud_lay_fst,icol)
      cloud_mask(icol,ilay,1:ngpt) = local_rands(1:ngpt) > (1._wp - cloud_frac(icol,ilay))
      do ilay = cloud_lay_fst+1, cloud_lay_lst
	    ! ################################################################################
        ! Max-random overlap
        !   new  random deviates if the adjacent layer isn't cloudy
        !   same random deviates if the adjacent layer is    cloudy
	    ! ################################################################################
	    if (.not. l_use_overlap_param) then
           if(cloud_mask_layer(ilay)) then
             if(.not. cloud_mask_layer(ilay-1)) local_rands(1:ngpt) = randoms(1:ngpt,ilay,icol)
             cloud_mask(icol,ilay,1:ngpt) = local_rands(1:ngpt) > (1._wp - cloud_frac(icol,ilay))
           else
             cloud_mask(icol,ilay,1:ngpt) = .false.
           end if
	    end if   ! END COND: Maximum-random overlap
	    ! ################################################################################
        ! Exponential-random overlap
        !   new  random deviates if the adjacent layer isn't cloudy
        !   correlated  deviates if the adjacent layer is    cloudy
	    ! ################################################################################
		if (l_use_overlap_param) then
           if(cloud_mask_layer(ilay)) then
               if(cloud_mask_layer(ilay-1)) then
               ! Create random deviates correlated between this layer and the previous layer
               !    (have to remove mean value before enforcing correlation).
               rho = overlap_param(icol,ilay-1)
               local_rands(1:ngpt) =  rho*(local_rands(1:ngpt)      -0.5_wp) + &
                      sqrt(1._wp-rho*rho)*(randoms(1:ngpt,ilay,icol)-0.5_wp) + 0.5_wp
             else
               local_rands(1:ngpt) = randoms(1:ngpt,ilay,icol)
             end if        
             cloud_mask(icol,ilay,1:ngpt) = local_rands(1:ngpt) > (1._wp - cloud_frac(icol,ilay))
  		  endif   
		endif   ! END COND: Exponential/Exponential-random overlap
	    ! ################################################################################
        ! Exponential-decorrelation overlap
        !   new  random deviates if the adjacent layer isn't cloudy
        !   correlated  deviates if the adjacent layer is    cloudy and decorrelation-length
	    ! ################################################################################		
		if (l_use_overlap_param .and. l_use_second_rng) then
	       where(randoms2(1:nGpt,iLay,iCol) .le. overlap_param(iCol,iLay))	      
	          cloud_mask(iCol,iLay,1:nGpt) = randoms(1:ngpt,iLay-1,iCol) > (1._wp - cloud_frac(iCol,iLay))
           elsewhere
              cloud_mask(iCol,iLay,1:nGpt) = randoms(1:ngpt,iLay,iCol)   > (1._wp - cloud_frac(iCol,iLay))
	       end where	
		endif 	! END COND: Exponential decorrelation-length
      end do    ! END LOOP: Layers
      
      ! Set cloud-mask in layer below clouds to false      
      cloud_mask(icol,cloud_lay_lst+1:nlay, 1:ngpt) = .false.
    end do		! END LOOP: Columns
    
  end subroutine sampled_mask
  ! -------------------------------------------------------------------------------------------------
  !
  ! Apply a true/false cloud mask to a homogeneous field
  !   This could be a kernel
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine apply_cloud_mask(ncol,nlay,nbnd,ngpt,band_lims_gpt,cloud_mask,input_field,sampled_field)
    integer,                                intent(in ) :: ncol,nlay,nbnd,ngpt
    integer,     dimension(2,nbnd),         intent(in ) :: band_lims_gpt
    logical,     dimension(ncol,nlay,ngpt), intent(in ) :: cloud_mask
    real(wp),    dimension(ncol,nlay,nbnd), intent(in ) :: input_field
    real(wp),    dimension(ncol,nlay,ngpt), intent(out) :: sampled_field

    integer :: icol,ilay,ibnd,igpt

    do ibnd = 1, nbnd
      do igpt = band_lims_gpt(1,ibnd), band_lims_gpt(2,ibnd)
        do ilay = 1, nlay
          sampled_field(1:ncol,ilay,igpt) = merge(input_field(1:ncol,ilay,ibnd), 0._wp, cloud_mask(1:ncol,ilay,igpt))
        end do
      end do
    end do
  end subroutine apply_cloud_mask

end module rrtmgp_sampling

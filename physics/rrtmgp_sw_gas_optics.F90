module rrtmgp_sw_gas_optics
  use machine,                only: kind_phys
  use module_radiation_gases, only: NF_VGAS
  use mo_rte_kind,            only: wl
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,  only: ty_gas_concs
  use rrtmgp_aux,             only: check_error_msg
  use mo_optical_props,       only: ty_optical_props_2str
  use mo_compute_bc,          only: compute_bc
  use netcdf

  implicit none

contains

  ! #########################################################################################
  ! SUBROUTINE sw_gas_optics_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_gas_optics_init
!! \htmlinclude rrtmgp_sw_gas_optics.html
!!
  subroutine rrtmgp_sw_gas_optics_init(rrtmgp_root_dir, rrtmgp_sw_file_gas, rrtmgp_nGases,   &
       active_gases_array, mpicomm, mpirank, mpiroot, sw_gas_props, errmsg, errflg)

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,  & ! RTE-RRTMGP root directory
         rrtmgp_sw_file_gas  ! RRTMGP file containing coefficients used to compute gaseous optical properties
    integer, intent(in) :: &
         rrtmgp_nGases       ! Number of trace gases active in RRTMGP
    character(len=*),dimension(rrtmgp_nGases), intent(in) :: &
         active_gases_array  ! Character array containing trace gases to include in RRTMGP
    integer,intent(in) :: &
         mpicomm,          & ! MPI communicator
         mpirank,          & ! Current MPI rank
         mpiroot             ! Master MPI rank
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg              ! CCPP error message
    integer,          intent(out) :: &
         errflg              ! CCPP error code
    type(ty_gas_optics_rrtmgp),intent(out) :: &
         sw_gas_props        ! RRTMGP DDT: shortwave spectral information

    ! Variables that will be passed to gas_optics%load()
    type(ty_gas_concs)  :: &
         gas_concentrations
    integer, dimension(:), allocatable :: &
         kminor_start_lower,              & ! Starting index in the [1, nContributors] vector for a contributor
                                            ! given by \"minor_gases_lower\" (lower atmosphere)  
         kminor_start_upper                 ! Starting index in the [1, nContributors] vector for a contributor 
                                            ! given by \"minor_gases_upper\" (upper atmosphere)   
    integer, dimension(:,:), allocatable :: &
         band2gpt,                        & ! Beginning and ending gpoint for each band    
         minor_limits_gpt_lower,          & ! Beginning and ending gpoint for each minor interval in lower atmosphere 
         minor_limits_gpt_upper             ! Beginning and ending gpoint for each minor interval in upper atmosphere 
    integer, dimension(:,:,:), allocatable :: &
         key_species                        ! Key species pair for each band 
    real(kind_phys) :: &
         press_ref_trop,                  & ! Reference pressure separating the lower and upper atmosphere [Pa]    
         temp_ref_p,                      & ! Standard spectroscopic reference pressure [Pa] 
         temp_ref_t,                      & ! Standard spectroscopic reference temperature [K] 
	 tsi_default,                     & !
	 mg_default,                      & !
	 sb_default                         !
    real(kind_phys), dimension(:), allocatable :: &
         press_ref,                       & ! Pressures for reference atmosphere; press_ref(# reference layers) [Pa] 
         temp_ref,                        & ! Temperatures for reference atmosphere; temp_ref(# reference layers) [K] 
    	 solar_quiet,                     & !
	 solar_facular,                   & !
	 solar_sunspot                      !
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims                          ! Beginning and ending wavenumber [cm -1] for each band                         
    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref,                         & ! Volume mixing ratios for reference atmosphere
         kminor_lower,                    & ! (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to
                                            ! [nTemp x nEta x nContributors] array)
         kminor_upper,                    & ! (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to
                                            ! [nTemp x nEta x nContributors] array)
         rayl_lower,                      & ! Stored coefficients due to rayleigh scattering contribution
         rayl_upper                         ! Stored coefficients due to rayleigh scattering contribution
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor                             ! Stored absorption coefficients due to major absorbing gases
    character(len=32),  dimension(:), allocatable :: &
         gas_names,                       & ! Names of absorbing gases
         gas_minor,                       & ! Name of absorbing minor gas
         identifier_minor,                & ! Unique string identifying minor gas
         minor_gases_lower,               & ! Names of minor absorbing gases in lower atmosphere
         minor_gases_upper,               & ! Names of minor absorbing gases in upper atmosphere
         scaling_gas_lower,               & ! Absorption also depends on the concentration of this gas
         scaling_gas_upper                  ! Absorption also depends on the concentration of this gas
    logical(wl), dimension(:), allocatable :: &
         minor_scales_with_density_lower, & ! Density scaling is applied to minor absorption coefficients
         minor_scales_with_density_upper, & ! Density scaling is applied to minor absorption coefficients
         scale_by_complement_lower,       & ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T)
         scale_by_complement_upper          ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T)
    ! Dimensions
    integer :: &
         ntemps, npress, ngpts, nabsorbers, nextrabsorbers,       &
         nminorabsorbers, nmixingfracs, nlayers, nbnds, npairs,   &
         nminor_absorber_intervals_lower, nminor_absorber_intervals_upper, &
         ncontributors_lower, ncontributors_upper

    ! Local variables
    integer :: status, ncid, dimid, varID, iGas
    integer,dimension(:),allocatable :: temp1, temp2, temp3, temp4
    character(len=264) :: sw_gas_props_file

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the gphysics_nml
    sw_gas_props_file   = trim(rrtmgp_root_dir)//trim(rrtmgp_sw_file_gas)

    ! Read dimensions for k-distribution fields (only on master processor(0))
!    if (mpirank .eq. mpiroot) then
       ! Open file
       status = nf90_open(trim(sw_gas_props_file), NF90_NOWRITE, ncid)

       ! Read dimensions for k-distribution fields
       status = nf90_inq_dimid(ncid, 'temperature', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ntemps)
       status = nf90_inq_dimid(ncid, 'pressure', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=npress)
       status = nf90_inq_dimid(ncid, 'absorber', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nabsorbers)
       status = nf90_inq_dimid(ncid, 'minor_absorber',dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminorabsorbers)
       status = nf90_inq_dimid(ncid, 'absorber_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nextrabsorbers)
       status = nf90_inq_dimid(ncid, 'mixing_fraction', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nmixingfracs)
       status = nf90_inq_dimid(ncid, 'atmos_layer', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nlayers)
       status = nf90_inq_dimid(ncid, 'bnd', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nbnds)
       status = nf90_inq_dimid(ncid, 'gpt', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ngpts)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=npairs)
       status = nf90_inq_dimid(ncid, 'contributors_lower',dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ncontributors_lower)
       status = nf90_inq_dimid(ncid, 'contributors_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ncontributors_upper)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_lower', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminor_absorber_intervals_lower)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminor_absorber_intervals_upper)
 
       ! Allocate space for arrays
       allocate(gas_names(nabsorbers))
       allocate(scaling_gas_lower(nminor_absorber_intervals_lower))
       allocate(scaling_gas_upper(nminor_absorber_intervals_upper))
       allocate(gas_minor(nminorabsorbers))
       allocate(identifier_minor(nminorabsorbers))
       allocate(minor_gases_lower(nminor_absorber_intervals_lower))
       allocate(minor_gases_upper(nminor_absorber_intervals_upper))
       allocate(minor_limits_gpt_lower(npairs,nminor_absorber_intervals_lower))
       allocate(minor_limits_gpt_upper(npairs,nminor_absorber_intervals_upper))
       allocate(band2gpt(2,nbnds))
       allocate(key_species(2,nlayers,nbnds))
       allocate(band_lims(2,nbnds))
       allocate(press_ref(npress))
       allocate(temp_ref(ntemps))
       allocate(vmr_ref(nlayers, nextrabsorbers, ntemps))
       allocate(kminor_lower(ncontributors_lower, nmixingfracs, ntemps))
       allocate(kmajor(ngpts, nmixingfracs,  npress+1, ntemps))
       allocate(kminor_start_lower(nminor_absorber_intervals_lower))
       allocate(kminor_upper(ncontributors_upper, nmixingfracs, ntemps))
       allocate(kminor_start_upper(nminor_absorber_intervals_upper))
       allocate(minor_scales_with_density_lower(nminor_absorber_intervals_lower))
       allocate(minor_scales_with_density_upper(nminor_absorber_intervals_upper))
       allocate(scale_by_complement_lower(nminor_absorber_intervals_lower))
       allocate(scale_by_complement_upper(nminor_absorber_intervals_upper))
       allocate(rayl_upper(ngpts, nmixingfracs, ntemps))
       allocate(rayl_lower(ngpts, nmixingfracs, ntemps))
       allocate(solar_quiet(ngpts))
       allocate(solar_facular(ngpts))	
       allocate(solar_sunspot(ngpts))
       allocate(temp1(nminor_absorber_intervals_lower))
       allocate(temp2(nminor_absorber_intervals_upper))
       allocate(temp3(nminor_absorber_intervals_lower))
       allocate(temp4(nminor_absorber_intervals_upper))

       ! Read in fields from file
       if (mpirank==mpiroot) write (*,*) 'Reading RRTMGP shortwave k-distribution data ... '
       status = nf90_inq_varid(ncid, 'gas_names', varID)
       status = nf90_get_var(  ncid, varID, gas_names)       
       status = nf90_inq_varid(ncid, 'scaling_gas_lower', varID)
       status = nf90_get_var(  ncid, varID, scaling_gas_lower)       
       status = nf90_inq_varid(ncid, 'scaling_gas_upper', varID)
       status = nf90_get_var(  ncid, varID, scaling_gas_upper)       
       status = nf90_inq_varid(ncid, 'gas_minor', varID)
       status = nf90_get_var(  ncid, varID, gas_minor)       
       status = nf90_inq_varid(ncid, 'identifier_minor', varID)
       status = nf90_get_var(  ncid, varID, identifier_minor)       
       status = nf90_inq_varid(ncid, 'minor_gases_lower', varID)
       status = nf90_get_var(  ncid, varID, minor_gases_lower)       
       status = nf90_inq_varid(ncid, 'minor_gases_upper', varID)
       status = nf90_get_var(  ncid, varID, minor_gases_upper)       
       status = nf90_inq_varid(ncid, 'minor_limits_gpt_lower', varID)
       status = nf90_get_var(  ncid, varID, minor_limits_gpt_lower)       
       status = nf90_inq_varid(ncid, 'minor_limits_gpt_upper', varID)
       status = nf90_get_var(  ncid, varID, minor_limits_gpt_upper)       
       status = nf90_inq_varid(ncid, 'bnd_limits_gpt', varID)
       status = nf90_get_var(  ncid, varID, band2gpt)       
       status = nf90_inq_varid(ncid, 'key_species', varID)
       status = nf90_get_var(  ncid, varID, key_species)       
       status = nf90_inq_varid(ncid,'bnd_limits_wavenumber', varID)
       status = nf90_get_var(  ncid, varID, band_lims)       
       status = nf90_inq_varid(ncid, 'press_ref', varID)
       status = nf90_get_var(  ncid, varID, press_ref)       
       status = nf90_inq_varid(ncid, 'temp_ref', varID)
       status = nf90_get_var(  ncid, varID, temp_ref)       
       status = nf90_inq_varid(ncid, 'absorption_coefficient_ref_P', varID)
       status = nf90_get_var(  ncid, varID, temp_ref_p)
       status = nf90_inq_varid(ncid, 'absorption_coefficient_ref_T', varID)
       status = nf90_get_var(  ncid, varID, temp_ref_t) 
       status = nf90_inq_varid(ncid, 'tsi_default', varID)
       status = nf90_get_var(  ncid, varID, tsi_default)
       status = nf90_inq_varid(ncid, 'mg_default', varID)
       status =	nf90_get_var(  ncid, varID, mg_default)
       status = nf90_inq_varid(ncid, 'sb_default', varID)
       status =	nf90_get_var(  ncid, varID, sb_default)
       status = nf90_inq_varid(ncid, 'press_ref_trop', varID)
       status = nf90_get_var(  ncid, varID, press_ref_trop)       
       status = nf90_inq_varid(ncid, 'kminor_lower', varID)
       status = nf90_get_var(  ncid, varID, kminor_lower)       
       status = nf90_inq_varid(ncid, 'kminor_upper', varID)
       status = nf90_get_var(  ncid, varID, kminor_upper)       
       status = nf90_inq_varid(ncid, 'vmr_ref', varID)
       status = nf90_get_var(  ncid, varID, vmr_ref)       
       status = nf90_inq_varid(ncid, 'kmajor', varID)
       status = nf90_get_var(  ncid, varID, kmajor)      
       status = nf90_inq_varid(ncid, 'kminor_start_lower', varID)
       status = nf90_get_var(  ncid, varID, kminor_start_lower)       
       status = nf90_inq_varid(ncid, 'kminor_start_upper', varID)
       status = nf90_get_var(  ncid, varID, kminor_start_upper)       
       status = nf90_inq_varid(ncid, 'solar_source_quiet', varID)
       status = nf90_get_var(  ncid, varID, solar_quiet)
       status = nf90_inq_varid(ncid, 'solar_source_facular', varID)
       status = nf90_get_var(  ncid, varID, solar_facular)
       status = nf90_inq_varid(ncid, 'solar_source_sunspot', varID)
       status = nf90_get_var(  ncid, varID, solar_sunspot)       
       status = nf90_inq_varid(ncid, 'rayl_lower', varID)
       status = nf90_get_var(  ncid, varID, rayl_lower)
       status = nf90_inq_varid(ncid, 'rayl_upper', varID)
       status = nf90_get_var(  ncid, varID, rayl_upper)

       ! Logical fields are read in as integers and then converted to logicals.
       status = nf90_inq_varid(ncid,'minor_scales_with_density_lower', varID)
       status = nf90_get_var(  ncid, varID,temp1)
       minor_scales_with_density_lower(:) = .false.
       where(temp1 .eq. 1) minor_scales_with_density_lower(:) = .true.
       status = nf90_inq_varid(ncid,'minor_scales_with_density_upper', varID)
       status = nf90_get_var(  ncid, varID,temp2)
       minor_scales_with_density_upper(:) = .false.
       where(temp2 .eq. 1) minor_scales_with_density_upper(:) = .true.
       status = nf90_inq_varid(ncid,'scale_by_complement_lower', varID)
       status = nf90_get_var(  ncid, varID,temp3)
       scale_by_complement_lower(:) = .false.
       where(temp3 .eq. 1) scale_by_complement_lower(:) = .true.       
       status = nf90_inq_varid(ncid,'scale_by_complement_upper', varID)
       status = nf90_get_var(  ncid, varID,temp4)
       scale_by_complement_upper(:) = .false.
       where(temp4 .eq. 1) scale_by_complement_upper(:) = .true.
       
       ! Close
       status = nf90_close(ncid)
!    endif


    ! Initialize gas concentrations and gas optics class
    call check_error_msg('sw_gas_optics_init',gas_concentrations%init(active_gases_array))
    call check_error_msg('sw_gas_optics_init',sw_gas_props%load(gas_concentrations, gas_names, &
         key_species, band2gpt, band_lims, press_ref, press_ref_trop, temp_ref, temp_ref_p,    &
         temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper, gas_minor, identifier_minor, &
         minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower,minor_limits_gpt_upper,  &
         minor_scales_with_density_lower, minor_scales_with_density_upper, scaling_gas_lower,  &
         scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,              &
         kminor_start_lower, kminor_start_upper, solar_quiet, solar_facular, solar_sunspot,    &
	 tsi_default, mg_default, sb_default, rayl_lower, rayl_upper)) 

  end subroutine rrtmgp_sw_gas_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_gas_optics_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_gas_optics_run
!! \htmlinclude rrtmgp_sw_gas_optics.html
!!
  subroutine rrtmgp_sw_gas_optics_run(doSWrad, nCol, nLev, nday, idxday, sw_gas_props, p_lay,&
       p_lev, toa_src_sw, t_lay, t_lev, gas_concentrations, solcon, rrtmgp_nGases,           &
       active_gases_array, sw_optical_props_clrsky, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doSWrad                 ! Flag to calculate SW irradiances
    integer,intent(in) :: &
         nDay,                 & ! Number of daylit points.
         nCol,                 & ! Number of horizontal points
         nLev                    ! Number of vertical levels
    integer,intent(in),dimension(ncol) :: &
         idxday                  ! Indices for daylit points.
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props            ! RRTMGP DDT: spectral information for RRTMGP SW radiation scheme
    real(kind_phys), dimension(ncol,nLev), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers (hPa)
         t_lay                   ! Temperature (K)
    real(kind_phys), dimension(ncol,nLev+1), intent(in) :: &
         p_lev,                & ! Pressure @ model layer-interfaces (hPa)
         t_lev                   ! Temperature @ model levels
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations (vmr)
    real(kind_phys), intent(in) :: &
         solcon                  ! Solar constant
    integer, intent(in) :: &
         rrtmgp_nGases           ! Number of trace gases active in RRTMGP
    character(len=*),dimension(rrtmgp_nGases), intent(in) :: &
         active_gases_array      ! Character array containing trace gases to include in RRTMGP

    ! Output
    character(len=*), intent(out) :: &
         errmsg                  ! CCPP error message
    integer,          intent(out) :: &
         errflg                  ! CCPP error code
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_clrsky ! RRTMGP DDT: clear-sky shortwave optical properties, spectral (tau,ssa,g) 
    real(kind_phys), dimension(ncol,sw_gas_props%get_ngpt()), intent(out) :: &
         toa_src_sw              ! TOA incident spectral flux (W/m2)

    ! Local variables
    integer :: ij,iGas
    real(kind_phys), dimension(ncol,nLev) :: vmrTemp
    real(kind_phys), dimension(nday,sw_gas_props%get_ngpt()) :: toa_src_sw_temp
    type(ty_gas_concs) :: &
         gas_concentrations_daylit    ! RRTMGP DDT: trace gas concentrations   (vmr)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doSWrad) return

    if (nDay .gt. 0) then
       ! Allocate space
       call check_error_msg('rrtmgp_sw_gas_optics_run',sw_optical_props_clrsky%alloc_2str(nday, nLev, sw_gas_props))

       ! Initialize gas concentrations and gas optics class
       call check_error_msg('rrtmgp_sw_rte_run',gas_concentrations_daylit%init(active_gases_array))

       ! Subset the gas concentrations, only need daylit points.
       do iGas=1,rrtmgp_nGases
          call check_error_msg('rrtmgp_sw_rte_run',&
               gas_concentrations%get_vmr(trim(active_gases_array(iGas)),vmrTemp))
          call check_error_msg('rrtmgp_sw_rte_run',&
               gas_concentrations_daylit%set_vmr(trim(active_gases_array(iGas)),vmrTemp(idxday(1:nday),:)))
       enddo

       ! Gas-optics
       call check_error_msg('rrtmgp_sw_gas_optics_run',sw_gas_props%gas_optics(&
            p_lay(idxday(1:nday),:),   & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(idxday(1:nday),:),   & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(idxday(1:nday),:),   & ! IN  - Temperature @ layer-centers (K)
            gas_concentrations_daylit, & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
            sw_optical_props_clrsky,   & ! OUT - RRTMGP DDT: Shortwave optical properties, by
                                         !                   spectral point (tau,ssa,g)
            toa_src_sw_temp))            ! OUT - TOA incident shortwave radiation (spectral)
       toa_src_sw(idxday(1:nday),:) = toa_src_sw_temp
       ! Scale incident flux
       do ij=1,nday
          toa_src_sw(idxday(ij),:) = toa_src_sw(idxday(ij),:)*solcon/ &
                                     sum(toa_src_sw(idxday(ij),:))
       enddo
    else
       toa_src_sw(:,:) = 0.
    endif

  end subroutine rrtmgp_sw_gas_optics_run

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_gas_optics_finalize
  ! #########################################################################################
  subroutine rrtmgp_sw_gas_optics_finalize()
  end subroutine rrtmgp_sw_gas_optics_finalize

end module rrtmgp_sw_gas_optics
 

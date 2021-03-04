module GFS_rrtmgp_sw_pre
  use machine, only: &
       kind_phys                   ! Working type
  use module_radiation_astronomy,only: &
       coszmn                      ! Function to compute cos(SZA)
  use module_radiation_surface,  only: &
       NF_ALBD,                  & ! Number of surface albedo categories (4; nir-direct, nir-diffuse, uvvis-direct, uvvis-diffuse)
       setalb                      ! Routine to compute surface albedo
  use surface_perturbation, only: & 
       cdfnor                      ! Routine to compute CDF (used to compute percentiles)
  use mo_gas_optics_rrtmgp,  only: &
       ty_gas_optics_rrtmgp
  use rrtmgp_sw_gas_optics,  only: sw_gas_props
  public GFS_rrtmgp_sw_pre_run,GFS_rrtmgp_sw_pre_init,GFS_rrtmgp_sw_pre_finalize
  
contains

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_init
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_pre_init ()
  end subroutine GFS_rrtmgp_sw_pre_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_sw_pre_run
!! \htmlinclude GFS_rrtmgp_sw_pre.html
!!
  subroutine GFS_rrtmgp_sw_pre_run(me, nCol, nLev, lndp_type, n_var_lndp,lndp_var_list,     &  
       lndp_prt_list, doSWrad, solhr, lon, coslat, sinlat,  snowd, sncovr, snoalb, zorl,    &
       tsfg, tsfa, hprime, alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc, albdvis,  &
       albdnir, albivis, albinir, lsmask, sfc_wts, p_lay, tv_lay, relhum, p_lev,            &
       nday, idxday, coszen, coszdg, sfc_alb_nir_dir, sfc_alb_nir_dif,                      &
       sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, sfc_alb_dif, errmsg, errflg)

    ! Inputs   
    integer, intent(in)    :: &
         me,                & ! Current MPI rank
         nCol,              & ! Number of horizontal grid points
         nLev,              & ! Number of vertical layers
         n_var_lndp,        &  ! Number of surface variables perturbed
         lndp_type             ! Type of land perturbations scheme used
    character(len=3), dimension(n_var_lndp), intent(in) ::  & 
         lndp_var_list
    real(kind_phys), dimension(n_var_lndp), intent(in) ::   &
         lndp_prt_list
    logical,intent(in) :: &
         doSWrad            ! Call RRTMGP SW radiation?
    real(kind_phys), intent(in) :: &
         solhr                 ! Time in hours after 00z at the current timestep
    real(kind_phys), dimension(nCol), intent(in) :: &
         lsmask,            & ! Landmask: sea/land/ice=0/1/2
         lon,               & ! Longitude
         coslat,            & ! Cosine(latitude)
         sinlat,            & ! Sine(latitude)
         snowd,             & ! Water equivalent snow depth (mm)
         sncovr,            & ! Surface snow area fraction (frac)
         snoalb,            & ! Maximum snow albedo (frac)
         zorl,              & ! Surface roughness length (cm)
         tsfg,              & ! Surface ground temperature for radiation (K)
         tsfa,              & ! Lowest model layer air temperature for radiation (K)         
         hprime,            & ! Standard deviation of subgrid orography (m)
         alvsf,             & ! Mean vis albedo with strong cosz dependency (frac)
         alnsf,             & ! Mean nir albedo with strong cosz dependency (frac)
         alvwf,             & ! Mean vis albedo with weak cosz dependency (frac)
         alnwf,             & ! Mean nir albedo with weak cosz dependency (frac)
         facsf,             & ! Fractional coverage with strong cosz dependency (frac)
         facwf,             & ! Fractional coverage with weak cosz dependency (frac)
         fice,              & ! Ice fraction over open water (frac)
         tisfc                ! Sea ice surface skin temperature (K)
    real(kind_phys), dimension(:), intent(in) :: &
         albdvis,           & ! surface albedo from lsm (direct,vis) (frac)
         albdnir,           & ! surface albedo from lsm (direct,nir) (frac)
         albivis,           & ! surface albedo from lsm (diffuse,vis) (frac)
         albinir              ! surface albedo from lsm (diffuse,nir) (frac)

    real(kind_phys), dimension(nCol,n_var_lndp), intent(in) :: &
         sfc_wts              ! Weights for stochastic surface physics perturbation ()    
    real(kind_phys), dimension(nCol,nLev),intent(in) :: &
         p_lay,             & ! Layer pressure
         tv_lay,            & ! Layer virtual-temperature
         relhum               ! Layer relative-humidity
    real(kind_phys), dimension(nCol,nLev+1),intent(in) :: &
         p_lev                ! Pressure @ layer interfaces (Pa)

    ! Outputs
    integer, intent(out)   :: &
         nday                 ! Number of daylit points
    integer, dimension(ncol), intent(out) :: &
         idxday               ! Indices for daylit points
    real(kind_phys), dimension(ncol), intent(inout) :: &
         coszen,            & ! Cosine of SZA
         coszdg,            & ! Cosine of SZA, daytime
         sfc_alb_dif          ! Mean surface diffused (nIR+uvvis) sw albedo
    real(kind_phys), dimension(sw_gas_props%get_nband(),ncol), intent(out) :: &
         sfc_alb_nir_dir,   & ! Surface albedo (direct) 
         sfc_alb_nir_dif,   & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir, & ! Surface albedo (direct)
         sfc_alb_uvvis_dif    ! Surface albedo (diffuse)
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag

    ! Local variables
    integer :: i, j, iCol, iBand, iLay
    real(kind_phys), dimension(ncol, NF_ALBD) :: sfcalb
    real(kind_phys), dimension(ncol) :: alb1d
    real(kind_phys) :: lndp_alb

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (doSWrad) then

       ! ####################################################################################
       ! Compute cosine of zenith angle (only when SW is called)
       ! ####################################################################################
       call coszmn (lon, sinlat, coslat, solhr, nCol, me, coszen, coszdg)

       ! ####################################################################################
       ! For SW gather daylit points
       ! ####################################################################################
       nday   = 0
       idxday = 0
       do i = 1, NCOL
          if (coszen(i) >= 0.0001) then
             nday = nday + 1
             idxday(nday) = i
          endif
       enddo
       
       ! ####################################################################################
       ! Call module_radiation_surface::setalb() to setup surface albedo.
       ! ####################################################################################
       alb1d(:) = 0.
       lndp_alb = -999.
       call setalb (lsmask, snowd, sncovr, snoalb, zorl, coszen, tsfg, tsfa, hprime, alvsf, &
            alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc, albdvis, albdnir, albivis,      &
            albinir, NCOL, alb1d, lndp_alb, sfcalb)
       
       ! Approximate mean surface albedo from vis- and nir-  diffuse values.
       sfc_alb_dif(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
  
       ! Spread across all SW bands
       do iBand=1,sw_gas_props%get_nband()
          sfc_alb_nir_dir(iBand,1:NCOL)   = sfcalb(1:NCOL,1)
          sfc_alb_nir_dif(iBand,1:NCOL)   = sfcalb(1:NCOL,2)
          sfc_alb_uvvis_dir(iBand,1:NCOL) = sfcalb(1:NCOL,3)
          sfc_alb_uvvis_dif(iBand,1:NCOL) = sfcalb(1:NCOL,4)
       enddo
    else
       nday                        = 0
       idxday                      = 0
       sfc_alb_nir_dir(:,1:nCol)   = 0.
       sfc_alb_nir_dif(:,1:nCol)   = 0.
       sfc_alb_uvvis_dir(:,1:nCol) = 0.
       sfc_alb_uvvis_dif(:,1:nCol) = 0.
       sfc_alb_dif(1:nCol)         = 0.
    endif


  end subroutine GFS_rrtmgp_sw_pre_run
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_pre_finalize ()
  end subroutine GFS_rrtmgp_sw_pre_finalize

end module GFS_rrtmgp_sw_pre

!>\file rrtmgp_aerosol_optics.F90
!!

module rrtmgp_aerosol_optics
  use machine,                   only: kind_phys
  use radiation_tools,           only: check_error_msg
  use rrtmgp_sw_gas_optics,      only: sw_gas_props
  use rrtmgp_lw_gas_optics,      only: lw_gas_props
  use module_radiation_aerosols, only: setaer
  use netcdf

  implicit none

  public rrtmgp_aerosol_optics_run

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_aerosol_optics_run()
  ! #########################################################################################

!>\defgroup rrtmgp_aerosol_optics_mod GFS RRTMGP Aerosol Optics Module
!> @{
!! \section arg_table_rrtmgp_aerosol_optics_run
!! \htmlinclude rrtmgp_aerosol_optics_run.html
!!
  subroutine rrtmgp_aerosol_optics_run(doSWrad, doLWrad, nCol, nLev, nDay, idxday, p_lev,   &
       p_lay, p_lk, tv_lay, relhum, lsmask, tracer, aerfld, lon, lat, iaermdl, iaerflg,     &
       top_at_1, con_pi, con_rd, con_g, aerodp, aerlw_tau, aerlw_ssa, aerlw_g, aersw_tau,   &
       aersw_ssa, aersw_g, errmsg, errflg  )

    ! Inputs
    logical, intent(in) :: &
         doSWrad,               & ! Logical flag for shortwave radiation call
         doLWrad,               & ! Logical flag for longwave radiation call 
         top_at_1                 ! Logical flag for vertical grid direcetion
    integer, intent(in) :: &
         nCol,                  & ! Number of horizontal grid points
         nDay,                  & ! Number of daylit points
         nLev,                  & ! Number of vertical layers
         iaermdl,               & ! Aerosol model scheme flag
         iaerflg                  ! Aerosol effects to include
    integer,intent(in),dimension(:) :: &
         idxday                   ! Indices for daylit points.
    real(kind_phys),intent(in) :: &
         con_pi,                & ! Physical constant (pi)
         con_rd,                & ! Physical constant (gas constant for dry-air)
         con_g                    ! Physical constant (gravitational constant)
    real(kind_phys), dimension(:), intent(in) :: &
         lon,                   & ! Longitude
         lat,                   & ! Latitude
         lsmask                   ! Land/sea/sea-ice mask
    real(kind_phys), dimension(:,:),intent(in) :: &
         p_lay,                 & ! Pressure @ layer-centers (Pa)
         tv_lay,                & ! Virtual-temperature @ layer-centers (K)
         relhum,                & ! Relative-humidity @ layer-centers
         p_lk                     ! Exner function @ layer-centers (1)
    real(kind_phys), dimension(:, :,:),intent(in) :: &
         tracer                   ! trace gas concentrations
    real(kind_phys), dimension(:, :,:),intent(in) :: &
         aerfld                   ! aerosol input concentrations
    real(kind_phys), dimension(:,:),intent(in) :: &
         p_lev                    ! Pressure @ layer-interfaces (Pa)

    ! Outputs
    real(kind_phys), dimension(:,:), intent(out) :: &
         aerodp                   ! Vertical integrated optical depth for various aerosol species 
    real(kind_phys), dimension(:,:,:), intent(out) :: &
         aerlw_tau,             & ! Longwave aerosol optical depth
         aerlw_ssa,             & ! Longwave aerosol single scattering albedo
         aerlw_g,               & ! Longwave aerosol asymmetry parameter
         aersw_tau,             & ! Shortwave aerosol optical depth 
         aersw_ssa,             & ! Shortwave aerosol single scattering albedo
         aersw_g                  ! Shortwave aerosol asymmetry parameter
    integer, intent(out) :: &
         errflg                   ! CCPP error flag
    character(len=*), intent(out) :: &
         errmsg                   ! CCPP error message

    ! Local variables
    real(kind_phys), dimension(nCol, nLev, lw_gas_props%get_nband(), 3) :: &
         aerosolslw            !
    real(kind_phys), dimension(nCol, nLev, sw_gas_props%get_nband(), 3) :: &
         aerosolssw, aerosolssw2
    integer :: iBand

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. (doSWrad .or. doLWrad)) return

    ! Call module_radiation_aerosols::setaer(),to setup aerosols property profile
    call setaer(p_lev*0.01, p_lay*0.01, p_lk, tv_lay, relhum, lsmask, tracer, aerfld, lon, lat, nCol, nLev, &
         nLev+1, .true., .true., iaermdl, iaerflg, top_at_1, con_pi, con_rd, con_g, aerosolssw2, aerosolslw, aerodp, errflg, errmsg)

    ! Shortwave
    if (doSWrad .and. (nDay .gt. 0)) then
       ! Store aerosol optical properties
       ! SW. 
       ! For RRTMGP SW the bands are now ordered from [IR(band) -> nIR -> UV], in RRTMG the 
       ! band ordering was [nIR -> UV -> IR(band)]
       aerosolssw(1:nCol,:,1,1)                          = aerosolssw2(1:nCol,:,sw_gas_props%get_nband(),1)
       aerosolssw(1:nCol,:,1,2)                          = aerosolssw2(1:nCol,:,sw_gas_props%get_nband(),2)
       aerosolssw(1:nCol,:,1,3)                          = aerosolssw2(1:nCol,:,sw_gas_props%get_nband(),3)
       aerosolssw(1:nCol,:,2:sw_gas_props%get_nband(),1) = aerosolssw2(1:nCol,:,1:sw_gas_props%get_nband()-1,1)
       aerosolssw(1:nCol,:,2:sw_gas_props%get_nband(),2) = aerosolssw2(1:nCol,:,1:sw_gas_props%get_nband()-1,2)
       aerosolssw(1:nCol,:,2:sw_gas_props%get_nband(),3) = aerosolssw2(1:nCol,:,1:sw_gas_props%get_nband()-1,3)
     
       ! Copy aerosol optical information/
       aersw_tau = aerosolssw(:,:,:,1)
       aersw_ssa = aerosolssw(:,:,:,2)
       aersw_g   = aerosolssw(:,:,:,3)
    endif

    ! Longwave
    if (doLWrad) then
       aerlw_tau = aerosolslw(:,:,:,1)
       aerlw_ssa = aerosolslw(:,:,:,2)
       aerlw_g   = aerosolslw(:,:,:,3)
    endif

  end subroutine rrtmgp_aerosol_optics_run
!> @}  
end module rrtmgp_aerosol_optics

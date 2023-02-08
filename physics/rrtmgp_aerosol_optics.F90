!>\file rrtmgp_aerosol_optics.F90
!!

module rrtmgp_aerosol_optics
  use machine,                   only: kind_phys
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_optical_props,          only: ty_optical_props_2str, ty_optical_props_1scl
  use radiation_tools,           only: check_error_msg
  use rrtmgp_sw_gas_optics,      only: sw_gas_props
  use rrtmgp_lw_gas_optics,      only: lw_gas_props
  use module_radiation_aerosols, only: &
       NF_AESW,                  & ! Number of optical-fields in SW output (3=tau+g+omega)
       NF_AELW,                  & ! Number of optical-fields in LW output (3=tau+g+omega)
       setaer,                   & ! Routine to compute aerosol radiative properties (tau,g,omega)
       NSPC1                       ! Number of species for vertically integrated aerosol optical-depth
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
  subroutine rrtmgp_aerosol_optics_run(doSWrad, doLWrad, nCol, nLev, nTracer, nTracerAer,   &
       nDay, idxday, p_lev, p_lay, p_lk, tv_lay, relhum, lsmask, tracer, aerfld, lon, lat,  &
       iaermdl, iaerflg, top_at_1, con_pi, con_rd, con_g, aerodp, sw_optical_props_aerosol, &
       lw_optical_props_aerosol, errmsg, errflg  )

    ! Inputs
    logical, intent(in) :: &
         doSWrad,               & ! Logical flag for shortwave radiation call
         doLWrad,               & ! Logical flag for longwave radiation call 
         top_at_1                 ! Logical flag for vertical grid direcetion
    integer, intent(in) :: &
         nCol,                  & ! Number of horizontal grid points
         nDay,                  & ! Number of daylit points
         nLev,                  & ! Number of vertical layers
         nTracer,               & ! Number of tracers
         nTracerAer,            & ! Number of aerosol tracers
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
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_aerosol ! RRTMGP DDT: Longwave aerosol optical properties (tau)
    type(ty_optical_props_1scl),intent(inout) :: &
         lw_optical_props_aerosol ! RRTMGP DDT: Longwave aerosol optical properties (tau)
    integer, intent(out) :: &
         errflg                   ! CCPP error flag
    character(len=*), intent(out) :: &
         errmsg                   ! CCPP error message

    ! Local variables
    real(kind_phys), dimension(nCol, nLev, lw_gas_props%get_nband(), NF_AELW) :: &
         aerosolslw            !
    real(kind_phys), dimension(nCol, nLev, sw_gas_props%get_nband(), NF_AESW) :: &
         aerosolssw, aerosolssw2
    integer :: iBand

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doSWrad) return

    ! Call module_radiation_aerosols::setaer(),to setup aerosols property profile
    call setaer(p_lev*0.01, p_lay*0.01, p_lk, tv_lay, relhum, lsmask, tracer, aerfld, lon, lat, nCol, nLev, &
         nLev+1, .true., .true., iaermdl, iaerflg, top_at_1, con_pi, con_rd, con_g, aerosolssw2, aerosolslw, aerodp, errflg, errmsg)

    ! Shortwave
    if (nDay .gt. 0) then
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
       
       ! Allocate RRTMGP DDT: Aerosol optics [nCol,nlev,nBands]
       call check_error_msg('rrtmgp_aerosol_optics_run',sw_optical_props_aerosol%alloc_2str(      &
            nDay, nlev, sw_gas_props%get_band_lims_wavenumber()))
       
       ! Copy aerosol optical information to RRTMGP DDT
       sw_optical_props_aerosol%tau = aerosolssw(idxday(1:nday),:,:,1)
       sw_optical_props_aerosol%ssa = aerosolssw(idxday(1:nday),:,:,2)
       sw_optical_props_aerosol%g   = aerosolssw(idxday(1:nday),:,:,3)
    endif

    ! Longwave
    if (.not. doLWrad) return
    lw_optical_props_aerosol%tau = aerosolslw(:,:,:,1) * (1. - aerosolslw(:,:,:,2))

    lw_optical_props_aerosol%band_lims_wvn = lw_gas_props%get_band_lims_wavenumber()
    do iBand=1,lw_gas_props%get_nband()
       lw_optical_props_aerosol%band2gpt(1:2,iBand) = iBand
       lw_optical_props_aerosol%gpt2band(iBand)     = iBand
    end do

  end subroutine rrtmgp_aerosol_optics_run
!> @}  
end module rrtmgp_aerosol_optics

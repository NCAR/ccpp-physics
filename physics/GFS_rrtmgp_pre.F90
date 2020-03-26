module GFS_rrtmgp_pre
  use physparam
  use machine, only: &
       kind_phys                   ! Working type
  use GFS_typedefs, only:        &
       GFS_statein_type,         & ! Prognostic state data in from dycore
       GFS_stateout_type,        & ! Prognostic state or tendencies return to dycore
       GFS_sfcprop_type,         & ! Surface fields
       GFS_coupling_type,        & ! Fields to/from coupling with other components (e.g. land/ice/ocean/etc.)
       GFS_control_type,         & ! Model control parameters
       GFS_grid_type,            & ! Grid and interpolation related data
       GFS_tbd_type,             & ! To-Be-Determined data that doesn't fit in any one container
       GFS_radtend_type,         & ! Radiation tendencies needed in physics
       GFS_diag_type               ! Fields targetted for diagnostic output
  use physcons, only:            &
       eps   => con_eps,         & ! Rd/Rv
       epsm1 => con_epsm1,       & ! Rd/Rv-1
       fvirt => con_fvirt,       & ! Rv/Rd-1
       rog   => con_rog            ! Rd/g
  use radcons, only: &
       qmin, epsq                  ! Minimum vlaues for varius calculations
  use funcphys, only:            &
       fpvs                        ! Function ot compute sat. vapor pressure over liq.
  use module_radiation_astronomy,only: &
       coszmn                      ! Function to compute cos(SZA)
  use module_radiation_gases,    only: &
       NF_VGAS,                  & ! Number of active gas species
       getgases,                 & ! Routine to setup trace gases
       getozn                      ! Routine to setup ozone
  use module_radiation_aerosols, only: &
       NF_AESW,                  & ! Number of optical-fields in SW output (3=tau+g+omega)
       NF_AELW,                  & ! Number of optical-fields in LW output (3=tau+g+omega)
       setaer,                   & ! Routine to compute aerosol radiative properties (tau,g,omega)
       NSPC1                       ! Number of species for vertically integrated aerosol optical-depth
  use module_radiation_clouds, only: &
       NF_CLDS,                  & ! Number of fields in "clouds" array (e.g. (cloud(1)=lwp,clouds(2)=ReffLiq,...)
       progcld1,                 & ! Zhao/Moorthi's prognostic cloud scheme
       progcld3,                 & ! Zhao/Moorthi's prognostic cloud+pdfcld
       progcld4,                 & ! GFDL cloud scheme
       progcld5,                 & ! Thompson / WSM6 cloud micrphysics scheme
       progclduni                  ! Unified cloud-scheme
  use surface_perturbation, only: & 
       cdfnor                      ! Routine to compute CDF (used to compute percentiles)
  use module_radiation_surface,  only: &
       setemis,                  & ! Routine to compute surface-emissivity
       NF_ALBD,                  & ! Number of surface albedo categories (4; nir-direct, nir-diffuse, uvvis-direct, uvvis-diffuse)
       setalb                      ! Routine to compute surface albedo
  ! RRTMGP types
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use rrtmgp_aux,            only: check_error_msg!, rrtmgp_minP, rrtmgp_minT
  use mo_rrtmgp_constants,   only: grav, avogad
  use mo_rrtmg_lw_cloud_optics

  real(kind_phys), parameter :: &
       amd   = 28.9644_kind_phys,  & ! Molecular weight of dry-air     (g/mol)
       amw   = 18.0154_kind_phys,  & ! Molecular weight of water vapor (g/mol)
       amo3  = 47.9982_kind_phys,  & ! Modelular weight of ozone       (g/mol)
       amdw  = amd/amw,            & ! Molecular weight of dry air / water vapor
       amdo3 = amd/amo3              ! Molecular weight of dry air / ozone

  ! Some common trace gas on/off flags. 
  ! This allows for control over which trace gases are used in RRTMGP radiation scheme via
  ! namelist.
  logical :: &
       isActive_h2o   = .false., & !
       isActive_co2   = .false., & !
       isActive_o3    = .false., & !
       isActive_n2o   = .false., & !
       isActive_ch4   = .false., & !
       isActive_o2    = .false., & !
       isActive_ccl4  = .false., & !
       isActive_cfc11 = .false., & !
       isActive_cfc12 = .false., & !
       isActive_cfc22 = .false.    !
  integer :: iStr_h2o, iStr_co2, iStr_o3, iStr_n2o, iStr_ch4, iStr_o2, iStr_ccl4, &
       iStr_cfc11, iStr_cfc12, iStr_cfc22 

  public GFS_rrtmgp_pre_run,GFS_rrtmgp_pre_init,GFS_rrtmgp_pre_finalize  
contains
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_pre_init
  ! #########################################################################################
!! \section arg_table_GFS_rrtmgp_pre_init
!! \htmlinclude GFS_rrtmgp_pre_init.html
!!
  subroutine GFS_rrtmgp_pre_init(Model, Radtend, active_gases_array, errmsg, errflg)
    ! Inputs
    type(GFS_control_type), intent(inout) :: &
         Model      ! DDT: FV3-GFS model control parameters
    type(GFS_radtend_type), intent(inout) :: &
         Radtend    ! DDT: FV3-GFS radiation tendencies 

    ! Outputs
    character(len=*),dimension(Model%ngases), intent(out) :: &
         active_gases_array  ! Character array containing trace gases to include in RRTMGP
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer, intent(out) :: &  
         errflg     ! Error flag

    ! Local variables
    character(len=1) :: tempstr
    integer :: ij, count
    integer,dimension(Model%ngases,2) :: gasIndices

    ! Initialize
    errmsg = ''
    errflg = 0

    if (len(Model%active_gases) .eq. 0) return
    
    ! Which gases are active? Provided via physics namelist.

    ! Pull out gas names from list...
    ! First grab indices in character array corresponding to start:end of gas name.
    gasIndices(1,1)=1
    count=1
    do ij=1,len(Model%active_gases)
       tempstr=trim(Model%active_gases(ij:ij))
       if (tempstr .eq. '_') then
          gasIndices(count,2)=ij-1
          gasIndices(count+1,1)=ij+1
          count=count+1
       endif
    enddo
    gasIndices(Model%ngases,2)=len(trim(Model%active_gases))
    
    ! Now extract the gas names
    do ij=1,Model%ngases
       active_gases_array(ij) = Model%active_gases(gasIndices(ij,1):gasIndices(ij,2))
    enddo

    ! Which gases are active? (This is purely for flexibility)
    do ij=1,Model%ngases
       if(trim(active_gases_array(ij)) .eq. 'h2o')   then
          isActive_h2o   = .true. 
          istr_h2o       = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'co2')   then
          isActive_co2   = .true.
          istr_co2       = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'o3')    then
          isActive_o3    = .true.
          istr_o3        = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'n2o')   then
          isActive_n2o   = .true.
          istr_n2o       = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'ch4')   then
          isActive_ch4   = .true.
          istr_ch4       = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'o2')    then
          isActive_o2    = .true.
          istr_o2        = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'ccl4')  then
          isActive_ccl4  = .true.
          istr_ccl4      = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'cfc11') then
          isActive_cfc11 = .true.
          istr_cfc11     = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'cfc12') then
          isActive_cfc12 = .true.
          istr_cfc12     = ij
       endif
       if(trim(active_gases_array(ij)) .eq. 'cfc22') then
          isActive_cfc22 = .true.       
          istr_cfc22     = ij
       endif
    enddo

  end subroutine GFS_rrtmgp_pre_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_pre_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_pre_run
!! \htmlinclude GFS_rrtmgp_pre_run.html
!!
  subroutine GFS_rrtmgp_pre_run (Model, Grid, Statein, Coupling, Radtend, Sfcprop, Tbd,  & ! IN
       ncol, lw_gas_props, active_gases_array,                                           & ! IN
       sec_diff_byband, raddt, p_lay, t_lay, p_lev, t_lev, tsfg, tsfa, cld_frac, cld_lwp,& ! OUT
       cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain,          & ! OUT
       tv_lay, relhum, tracer, cldsa, mtopa, mbota, de_lgth,  gas_concentrations,        & ! OUT
       errmsg, errflg)
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! DDT: FV3-GFS model control parameters
    type(GFS_grid_type), intent(in) :: &
         Grid                 ! DDT: FV3-GFS grid and interpolation related data 
    type(GFS_statein_type), intent(in) :: &
         Statein              ! DDT: FV3-GFS prognostic state data in from dycore    
    type(GFS_coupling_type), intent(in) :: &
         Coupling             ! DDT: FV3-GFS fields to/from coupling with other components 
    type(GFS_radtend_type), intent(inout) :: &
         Radtend              ! DDT: FV3-GFS radiation tendencies 
    type(GFS_sfcprop_type), intent(in) :: &
         Sfcprop              ! DDT: FV3-GFS surface fields
    type(GFS_tbd_type), intent(in) :: &
         Tbd                  ! DDT: FV3-GFS data not yet assigned to a defined container
    integer, intent(in)    :: &
         ncol                 ! Number of horizontal grid points
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props         ! RRTMGP DDT: longwave spectral information
    character(len=*),dimension(Model%ngases), intent(in) :: &
         active_gases_array   ! Character array containing trace gases to include in RRTMGP

    ! Outputs
    real(kind_phys), dimension(ncol,Model%levs), intent(out) :: &
         p_lay,             & ! Pressure at model-layer
         t_lay                ! Temperature at model layer
    real(kind_phys), dimension(ncol,Model%levs+1), intent(out) :: &
         p_lev,             & ! Pressure at model-interface
         t_lev                ! Temperature at model-interface
    real(kind_phys), intent(out) :: &
         raddt                ! Radiation time-step
    real(kind_phys), dimension(ncol), intent(out) :: &
         tsfg,              & ! Ground temperature
         tsfa                 ! Skin temperature
    type(ty_gas_concs),intent(out) :: &
         gas_concentrations   ! RRTMGP DDT: gas volumne mixing ratios
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag
    real(kind_phys), dimension(ncol,Model%levs),intent(out) :: &
         cld_frac,          & ! Total cloud fraction
         cld_lwp,           & ! Cloud liquid water path
         cld_reliq,         & ! Cloud liquid effective radius
         cld_iwp,           & ! Cloud ice water path
         cld_reice,         & ! Cloud ice effecive radius
         cld_swp,           & ! Cloud snow water path
         cld_resnow,        & ! Cloud snow effective radius
         cld_rwp,           & ! Cloud rain water path
         cld_rerain           ! Cloud rain effective radius
    real(kind_phys), dimension(ncol,Model%levs),intent(out) :: &
         tv_lay,            & ! Virtual temperatue at model-layers 
         relhum               ! Relative-humidity at model-layers 
    real(kind_phys), dimension(ncol, Model%levs, 2:Model%ntrac),intent(out) :: &
         tracer               ! Array containing trace gases
    integer,dimension(ncol,3),intent(out) :: &
         mbota,             & ! Vertical indices for cloud tops
         mtopa                ! Vertical indices for cloud bases
    real(kind_phys), dimension(ncol,5), intent(out) :: &
         cldsa                ! Fraction of clouds for low, middle, high, total and BL 
    real(kind_phys), dimension(ncol), intent(out)  :: &
         de_lgth              ! Decorrelation length
    real(kind_phys), dimension(lw_gas_props%get_nband(),ncol),intent(out) :: &
         sec_diff_byband

    ! Local variables
    integer :: i, j, iCol, iBand, iSFC, iTOA, iLay
    logical :: top_at_1
    real(kind_phys),dimension(NCOL,Model%levs) :: vmr_o3, vmr_h2o, coldry, tem0, colamt
    real(kind_phys) :: es, qs, tem1, tem2
    real(kind_phys), dimension(ncol, NF_ALBD) :: sfcalb
    real(kind_phys), dimension(ncol, Model%levs) :: qs_lay, q_lay, deltaZ, deltaP, o3_lay
    real(kind_phys), dimension(ncol, Model%levs, NF_VGAS) :: gas_vmr
    real(kind_phys), dimension(ncol, Model%levs, NF_CLDS) :: clouds
    real(kind_phys), dimension(ncol) ::  precipitableH2o

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. (Model%lsswr .or. Model%lslwr)) return
    
    ! #######################################################################################
    ! What is vertical ordering?
    ! #######################################################################################
    top_at_1 = (Statein%prsi(1,1) .lt.  Statein%prsi(1, Model%levs))
    if (top_at_1) then 
       iSFC = Model%levs
       iTOA = 1
    else
       iSFC = 1
       iTOA = Model%levs
    endif

    ! #######################################################################################
    ! Compute some fields needed by RRTMGP
    ! #######################################################################################
    
    ! Water-vapor mixing-ratio
    q_lay(1:ncol,:)  = Statein%qgrs(1:NCOL,:,1)
    where(q_lay .lt. 1.e-6) q_lay = 1.e-6
    
    ! Pressure at layer-interface
    p_lev(1:NCOL,:) = Statein%prsi(1:NCOL,:)

    ! Pressure at layer-center
    p_lay(1:NCOL,:)   = Statein%prsl(1:NCOL,:)

    ! Temperature at layer-center
    t_lay(1:NCOL,:) = Statein%tgrs(1:NCOL,:)

    ! Temperature at layer-interfaces
    if (top_at_1) then
       t_lev(1:NCOL,1)      = t_lay(1:NCOL,iTOA)
       t_lev(1:NCOL,2:iSFC) = (t_lay(1:NCOL,2:iSFC)+t_lay(1:NCOL,1:iSFC-1))/2._kind_phys
       t_lev(1:NCOL,iSFC+1) = Sfcprop%tsfc(1:NCOL)
    else
       t_lev(1:NCOL,1)      = Sfcprop%tsfc(1:NCOL)
       t_lev(1:NCOL,2:iTOA) = (t_lay(1:NCOL,2:iTOA)+t_lay(1:NCOL,1:iTOA-1))/2._kind_phys
       t_lev(1:NCOL,iTOA+1) = t_lay(1:NCOL,iTOA)
    endif
    
    ! Compute layer pressure thicknes
    deltaP = abs(p_lev(:,2:model%levs+1)-p_lev(:,1:model%levs))

    ! Compute a bunch of thermodynamic fields needed by the macrophysics schemes. Relative humidity, 
    ! saturation mixing-ratio, vapor mixing-ratio, virtual temperature, layer thickness,...
    do iCol=1,NCOL
       do iLay=1,Model%levs
          es                = min( p_lay(iCol,iLay),  fpvs( t_lay(iCol,iLay) ) )  ! fpvs and prsl in pa
          qs                = max( QMIN, eps * es / (p_lay(iCol,iLay) + epsm1*es) )
          relhum(iCol,iLay) = max( 0._kind_phys, min( 1._kind_phys, max(QMIN, q_lay(iCol,iLay))/qs ) )
          qs_lay(iCol,iLay) = qs
          tv_lay(iCol,iLay) = t_lay(iCol,iLay) * (1._kind_phys + fvirt*q_lay(iCol,iLay)) 
          deltaZ(iCol,iLay) = (rog*0.001) * abs(log(p_lev(iCol,iLay)) - log(p_lev(iCol,iLay+1))) * tv_lay(iCol,iLay)
       enddo
    enddo

    ! #######################################################################################
    ! Get layer ozone mass mixing ratio 
    ! #######################################################################################
    ! First recast remaining all tracers (except sphum) forcing them all to be positive
    do j = 2, model%NTRAC
       tracer(1:NCOL,:,j) = Statein%qgrs(1:NCOL,:,j)
       where(tracer(:,:,j) .lt. 0.0) tracer(:,:,j) = 0._kind_phys
    enddo

    if (Model%ntoz > 0) then 
       do iLay=1,Model%levs
          do iCol=1,NCOL
             o3_lay(iCol,iLay) = max( QMIN, tracer(iCol,iLay,Model%ntoz) )
          enddo
       enddo
    ! OR Use climatological ozone data
    else                               
       call getozn (Statein%prslk(1:NCOL,:), Grid%xlat, NCOL, Model%levs, o3_lay)
    endif

    ! #######################################################################################
    ! Set gas concentrations for RRTMGP
    ! #######################################################################################
    ! Call getgases(), to set up non-prognostic gas volume mixing ratios (gas_vmr).
    call getgases (p_lev/100., Grid%xlon, Grid%xlat, NCOL, Model%levs, gas_vmr)

    ! Compute volume mixing-ratios for ozone (mmr) and specific-humidity.
    vmr_h2o = merge((q_lay/(1-q_lay))*amdw, 0., q_lay  .ne. 1.)
    vmr_o3  = merge(o3_lay*amdo3,           0., o3_lay .gt. 0.)
    
    ! Initialize and opulate RRTMGP DDT w/ gas-concentrations
    call check_error_msg('sw_gas_optics_init',gas_concentrations%init(active_gases_array))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr(active_gases_array(iStr_o2),  gas_vmr(:,:,4)))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr(active_gases_array(iStr_co2), gas_vmr(:,:,1)))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr(active_gases_array(iStr_ch4), gas_vmr(:,:,3)))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr(active_gases_array(iStr_n2o), gas_vmr(:,:,2)))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr(active_gases_array(iStr_h2o), vmr_h2o))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr(active_gases_array(iStr_o3),  vmr_o3))

    ! #######################################################################################
    ! Compute diffusivity angle adjustments for each longwave band
    ! *NOTE* Legacy RRTMGP code 
    ! #######################################################################################
    ! Conpute diffusivity angle adjustments.
    ! First need to compute precipitable water in each column
    tem0   = (1._kind_phys - vmr_h2o)*amd + vmr_h2o*amw
    coldry = ( 1.0e-20 * 1.0e3 *avogad)*(deltap*.01) / (100.*grav*tem0*(1._kind_phys + vmr_h2o))
    colamt = max(0._kind_phys, coldry*vmr_h2o)
    do iCol=1,nCol
       tem1   = 0._kind_phys
       tem2   = 0._kind_phys   
       do iLay=1,Model%levs
          tem1 = tem1 + coldry(iCol,iLay)+colamt(iCol,iLay)
          tem2 = tem2 + colamt(iCol,iLay)
       enddo
       precipitableH2o(iCol) = p_lev(iCol,iSFC)*0.01*(10._kind_phys*tem2 / (amdw*tem1*grav))
    enddo

    ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  the function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.    
    do iCol=1,nCol
       do iBand = 1, lw_gas_props%get_nband()
          if (iBand==1 .or. iBand==4 .or. iBand==10) then
             sec_diff_byband(iBand,iCol) = diffusivityB1410
          else
             sec_diff_byband(iBand,iCol) = min( diffusivityHigh, max(diffusivityLow, &
                  a0(iBand)+a1(iBand)*exp(a2(iBand)*precipitableH2o(iCol))))
          endif
       enddo
    enddo

    ! #######################################################################################
    ! Radiation time step (output) (Is this really needed?) (Used by some diangostics)
    ! #######################################################################################
    raddt = min(Model%fhswr, Model%fhlwr)

    ! #######################################################################################
    ! Setup surface ground temperature and ground/air skin temperature if required.
    ! #######################################################################################
    tsfg(1:NCOL) = Sfcprop%tsfc(1:NCOL)
    tsfa(1:NCOL) = Sfcprop%tsfc(1:NCOL)

    ! #######################################################################################
    ! Cloud microphysics
    ! #######################################################################################
    call cloud_microphysics(Model, Tbd, Grid, Sfcprop, ncol, tracer, p_lay, t_lay,  p_lev,  &
         tv_lay, relhum, qs_lay, q_lay, deltaZ, deltaP, clouds, cldsa, mbota, mtopa, de_lgth)

    ! Copy output cloud fields
    cld_frac   = clouds(:,:,1)
    cld_lwp    = clouds(:,:,2)
    cld_reliq  = clouds(:,:,3)
    cld_iwp    = clouds(:,:,4)
    cld_reice  = clouds(:,:,5)   
    cld_rwp    = clouds(:,:,6)  
    cld_rerain = clouds(:,:,7)  
    cld_swp    = clouds(:,:,8)  
    cld_resnow = clouds(:,:,9)  
    
  end subroutine GFS_rrtmgp_pre_run
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_pre_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_pre_finalize ()
  end subroutine GFS_rrtmgp_pre_finalize

  ! #########################################################################################
  ! Subroutine cloud_microphysics()
  ! #########################################################################################
  subroutine cloud_microphysics(Model, Tbd, Grid, Sfcprop, ncol, tracer, p_lay, t_lay, p_lev,&
       tv_lay, relhum, qs_lay, q_lay, deltaZ, deltaP, clouds, cldsa, mbota, mtopa, de_lgth)

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! DDT: FV3-GFS model control parameters
    type(GFS_tbd_type), intent(in) :: &
         Tbd                  ! DDT: FV3-GFS data not yet assigned to a defined container
    type(GFS_grid_type), intent(in) :: &
         Grid                 ! DDT: FV3-GFS grid and interpolation related data 
    type(GFS_sfcprop_type), intent(in) :: &
         Sfcprop              ! DDT: FV3-GFS surface fields
    integer, intent(in) :: &
         ncol                 ! Number of horizontal gridpoints
    real(kind_phys), dimension(ncol, Model%levs, 2:Model%ntrac),intent(in) :: &
         tracer               ! Cloud condensate amount in layer by type ()
    real(kind_phys), dimension(ncol,Model%levs), intent(in) :: &
         p_lay,             & ! Pressure @ model layer centers (Pa)
         t_lay,             & ! Temperature @ layer centers (K)
         tv_lay,            & ! Virtual temperature @ layer centers (K)
         relhum,            & ! Relative humidity @ layer centers(1)
         qs_lay,            & ! Saturation specific humidity @ layer center   (kg/kg)
         q_lay,             & ! Specific humidity @ layer centers(kg/kg)
         deltaZ,            & ! Layer thickness (km)
         deltaP               ! Layer thickness (Pa)
    real(kind_phys), dimension(ncol,Model%levs+1), intent(in) :: &
         p_lev                ! Pressure @ model layer interface (Pa)

    ! Outputs
    real(kind_phys), dimension(ncol, Model%levs, NF_CLDS),intent(out) :: &
         clouds               ! Cloud properties (NCOL,Model%levs,NF_CLDS)
    integer,dimension(ncol,3), intent(out) :: &
         mbota,             & ! Vertical indices for low, mid, hi cloud bases (NCOL,3)
         mtopa                ! Vertical indices for low, mid, hi cloud tops (NCOL,3)
    real(kind_phys), dimension(ncol), intent(out)  ::&
         de_lgth              ! Clouds decorrelation length (km)
    real(kind_phys), dimension(ncol, 5), intent(out) :: &
         cldsa                ! Fraction of clouds for low, mid, hi, tot, bl (NCOL,5)

    ! Local variables
    real(kind_phys), dimension(ncol, Model%levs, Model%ncnd) :: cld_condensate
    integer :: i,k
    real(kind_phys), parameter :: xrc3 = 100.
    real(kind_phys), dimension(ncol, Model%levs) :: delta_q, cnv_w, cnv_c, effr_l, &
         effr_i, effr_r, effr_s, cldcov

    ! #######################################################################################
    !  Obtain cloud information for radiation calculations
    !    (clouds,cldsa,mtopa,mbota)
    !   for  prognostic cloud:
    !    - For Zhao/Moorthi's prognostic cloud scheme,
    !      call module_radiation_clouds::progcld1()
    !    - For Zhao/Moorthi's prognostic cloud+pdfcld,
    !      call module_radiation_clouds::progcld3()
    !      call module_radiation_clouds::progclduni() for unified cloud and ncld=2
    ! #######################################################################################
    cld_condensate = 0.0_kind_phys
    if (Model%ncnd == 1) then                                                                    ! Zhao_Carr_Sundqvist
       cld_condensate(1:NCOL,1:Model%levs,1) = tracer(1:NCOL,1:Model%levs,Model%ntcw)            ! -liquid water/ice
    elseif (Model%ncnd == 2) then                                                                ! MG
       cld_condensate(1:NCOL,1:Model%levs,1) = tracer(1:NCOL,1:Model%levs,Model%ntcw)            ! -liquid water
       cld_condensate(1:NCOL,1:Model%levs,2) = tracer(1:NCOL,1:Model%levs,Model%ntiw)            ! -ice water
    elseif (Model%ncnd == 4) then                                                                ! MG2
       cld_condensate(1:NCOL,1:Model%levs,1) = tracer(1:NCOL,1:Model%levs,Model%ntcw)            ! -liquid water
       cld_condensate(1:NCOL,1:Model%levs,2) = tracer(1:NCOL,1:Model%levs,Model%ntiw)            ! -ice water
       cld_condensate(1:NCOL,1:Model%levs,3) = tracer(1:NCOL,1:Model%levs,Model%ntrw)            ! -rain water
       cld_condensate(1:NCOL,1:Model%levs,4) = tracer(1:NCOL,1:Model%levs,Model%ntsw)            ! -snow water
    elseif (Model%ncnd == 5) then                                                                ! GFDL MP, Thompson, MG3
       cld_condensate(1:NCOL,1:Model%levs,1) = tracer(1:NCOL,1:Model%levs,Model%ntcw)            ! -liquid water
       cld_condensate(1:NCOL,1:Model%levs,2) = tracer(1:NCOL,1:Model%levs,Model%ntiw)            ! -ice water
       cld_condensate(1:NCOL,1:Model%levs,3) = tracer(1:NCOL,1:Model%levs,Model%ntrw)            ! -rain water
       cld_condensate(1:NCOL,1:Model%levs,4) = tracer(1:NCOL,1:Model%levs,Model%ntsw) + &        ! -snow + grapuel
                                               tracer(1:NCOL,1:Model%levs,Model%ntgl) 
    endif
    where(cld_condensate < epsq) cld_condensate = 0.0
    
    ! For GFDL microphysics scheme...
    if (Model%imp_physics == 11 ) then
       if (.not. Model%lgfdlmprad) then
          cld_condensate(:,:,1) =                         tracer(:,1:Model%levs,Model%ntcw)
          cld_condensate(:,:,1) = cld_condensate(:,:,1) + tracer(:,1:Model%levs,Model%ntrw)
          cld_condensate(:,:,1) = cld_condensate(:,:,1) + tracer(:,1:Model%levs,Model%ntiw)
          cld_condensate(:,:,1) = cld_condensate(:,:,1) + tracer(:,1:Model%levs,Model%ntsw)
          cld_condensate(:,:,1) = cld_condensate(:,:,1) + tracer(:,1:Model%levs,Model%ntgl)
       endif
       do k=1,Model%levs
          do i=1,NCOL
             if (cld_condensate(i,k,1) < EPSQ ) cld_condensate(i,k,1) = 0.0
          enddo
       enddo
    endif

    if (Model%uni_cld) then
       if (Model%effr_in) then
          cldcov(:,:)  = Tbd%phy_f3d(:,:,Model%indcld)
          effr_l(:,:)  = Tbd%phy_f3d(:,:,2)
          effr_i(:,:)  = Tbd%phy_f3d(:,:,3)
          effr_r(:,:)  = Tbd%phy_f3d(:,:,4)
          effr_s(:,:)  = Tbd%phy_f3d(:,:,5)
       else
          do k=1,model%levs
             do i=1,ncol
                cldcov(i,k) = Tbd%phy_f3d(i,k,Model%indcld)
             enddo
          enddo
       endif
    elseif (Model%imp_physics == Model%imp_physics_gfdl) then                          ! GFDL MP
       cldcov(1:NCOL,1:Model%levs) = tracer(1:NCOL,1:Model%levs,Model%ntclamt)
       if (Model%effr_in) then
          effr_l(:,:)  = Tbd%phy_f3d(:,:,1)
          effr_i(:,:)  = Tbd%phy_f3d(:,:,2)
          effr_r(:,:)  = Tbd%phy_f3d(:,:,3)
          effr_s(:,:)  = Tbd%phy_f3d(:,:,4)
       endif
    else                                                           ! neither of the other two cases
       cldcov = 0.0
    endif


    ! Add suspended convective cloud water to grid-scale cloud water
    ! only for cloud fraction & radiation computation it is to enhance 
    ! cloudiness due to suspended convec cloud water for zhao/moorthi's 
    ! (imp_phys=99) & ferrier's (imp_phys=5) microphysics schemes
    if ((Model%num_p3d == 4) .and. (Model%npdf3d == 3)) then       ! same as Model%imp_physics = 99
       delta_q(1:ncol,1:Model%levs) = Tbd%phy_f3d(1:ncol,Model%levs:1:-1,5)
       cnv_w  (1:ncol,1:Model%levs) = Tbd%phy_f3d(1:ncol,Model%levs:1:-1,6)
       cnv_c  (1:ncol,1:Model%levs) = Tbd%phy_f3d(1:ncol,Model%levs:1:-1,7)
    elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then ! same as MOdel%imp_physics=98
       delta_q(1:ncol,1:Model%levs) = 0.0
       cnv_w  (1:ncol,1:Model%levs) = Tbd%phy_f3d(1:ncol,1:Model%levs,Model%num_p3d+1)
       cnv_c  (1:ncol,1:Model%levs) = 0.0
    else                                                           ! all the rest
       delta_q(1:ncol,1:Model%levs) = 0.0
       cnv_w  (1:ncol,1:Model%levs) = 0.0
       cnv_c  (1:ncol,1:Model%levs) = 0.0
    endif

    ! For zhao/moorthi's prognostic cloud scheme, add in convective cloud water to liquid-cloud water
    if (Model%imp_physics == 99) then
       cld_condensate(1:NCOL,1:Model%levs,1) = cld_condensate(1:NCOL,1:Model%levs,1) + cnv_w(1:NCOL,1:Model%levs)
    endif
    
    ! For MG prognostic cloud scheme, add in convective cloud water to liquid-and-ice-cloud condensate
    if (Model%imp_physics == 10) then
       cld_condensate(1:NCOL,1:Model%levs,1) = cld_condensate(1:NCOL,1:Model%levs,1) + cnv_w(1:NCOL,1:Model%levs) + cld_condensate(1:NCOL,1:Model%levs,2)
    endif

    ! #######################################################################################
    ! MICROPHYSICS
    ! #######################################################################################
    ! *) zhao/moorthi's prognostic cloud scheme or unified cloud and/or with MG microphysics
    if (Model%imp_physics == 99 .or. Model%imp_physics == 10) then           
       if (Model%uni_cld .and. Model%ncld >= 2) then
          call progclduni(           &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               cld_condensate,       & ! IN  - Cloud condensate amount (Model%ncnd types)     ()
               Model%ncnd,           & ! IN  - Number of cloud condensate types               ()
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               deltaZ,               & ! IN  - Layer thickness                                (km)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               cldcov,               & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
               effr_l,               & ! IN  - Liquid-water effective radius                  (microns)
               effr_i,               & ! IN  - Ice-water effective radius                     (microns)
               effr_r,               & ! IN  - Rain-water effective radius                    (microns)
               effr_s,               & ! IN  - Snow-water effective radius                    (microns)
               Model%effr_in,        & ! IN  - Logical, if .true. use input effective radii
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       else
          call progcld1 (            &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               q_lay,                & ! IN  - Specific humidity at layer center              (kg/kg)
               qs_lay,               & ! IN  - Saturation specific humidity at layer center   (kg/kg)
               relhum,               & ! IN  - Relative humidity at layer center              (1)
               cld_condensate(:,:,1),& ! IN  - Cloud condensate amount                        ()
                                       !       (Zhao: liq+convective; MG: liq+ice+convective) 
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               deltaZ,               & ! IN  - Layer thickness                                (km)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               Model%uni_cld,        & ! IN  - True for cloud fraction from shoc
               Model%lmfshal,        & ! IN  - True for mass flux shallow convection
               Model%lmfdeep2,       & ! IN  - True for mass flux deep convection
               cldcov,               & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
               effr_l,               & ! IN  - Liquid-water effective radius                  (microns)
               effr_i,               & ! IN  - Ice-water effective radius                     (microns)
               effr_r,               & ! IN  - Rain-water effective radius                    (microns)
               effr_s,               & ! IN  - Snow-water effective radius                    (microns)
               Model%effr_in,        & ! IN  - Logical, if .true. use input effective radii
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       endif
       ! *) zhao/moorthi's prognostic cloud+pdfcld
    elseif(Model%imp_physics == 98) then
       call progcld3 (               &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               q_lay,                & ! IN  - Specific humidity at layer center              (kg/kg)
               qs_lay,               & ! IN  - Saturation specific humidity at layer center   (kg/kg)
               relhum,               & ! IN  - Relative humidity at layer center              (1)
               cld_condensate(:,:,1),& ! IN  - Cloud condensate amount (only h20)             ()
               cnv_w,                & ! IN  - Layer convective cloud condensate
               cnv_c,                & ! IN  - Layer convective cloud cover
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               deltaZ,               & ! IN  - Layer thickness                                (km)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               delta_q,              & ! IN  - Total water distribution width
               Model%sup,            & ! IN  - ??? Supersaturation?
               Model%kdt,            & ! IN  - ??? 
               Model%me,             & ! IN  - ??? NOT USED IN PROGCLD3()
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       ! *) GFDL cloud scheme
    elseif (Model%imp_physics == 11) then 
       if (.not.Model%lgfdlmprad) then
          call progcld4 (            &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               q_lay,                & ! IN  - Specific humidity at layer center              (kg/kg)
               qs_lay,               & ! IN  - Saturation specific humidity at layer center   (kg/kg)
               relhum,               & ! IN  - Relative humidity at layer center              (1)
               cld_condensate(:,:,1),& ! IN  - Cloud condensate amount (only h20)             ()
               cnv_w,                & ! IN  - Layer convective cloud condensate
               cnv_c,                & ! IN  - Layer convective cloud cover
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               cldcov,               & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
               deltaZ,               & ! IN  - Layer thickness                                (km)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       else
          call progclduni(           &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               cld_condensate,       & ! IN  - Cloud condensate amount (Model%ncnd types)     ()
               Model%ncnd,           & ! IN  - Number of cloud condensate types               ()
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               deltaZ,               & ! IN  - Layer thickness                                (km)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               cldcov,               & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
               effr_l,               & ! IN  - Liquid-water effective radius                  (microns)
               effr_i,               & ! IN  - Ice-water effective radius                     (microns)
               effr_r,               & ! IN  - Rain-water effective radius                    (microns)
               effr_s,               & ! IN  - Snow-water effective radius                    (microns)
               Model%effr_in,        & ! IN  - Logical, if .true. use input effective radii
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       endif
       ! *) Thompson / WSM6 cloud micrphysics scheme
    elseif(Model%imp_physics == 8 .or. Model%imp_physics == 6) then
 
       call progcld5 ( & ! IN
            p_lay/100.,               & ! IN  - Pressure at model layer centers                (mb)
            p_lev/100.,               & ! IN  - Pressure at model interfaces                   (mb)
            t_lay,                    & ! IN  - Temperature at layer centers                   (K)
            q_lay,                    & ! IN  - Specific humidity at layer center              (kg/kg)
            qs_lay,                   & ! IN  - Saturation specific humidity at layer center   (kg/kg)
            relhum,                   & ! IN  - Relative humidity at layer center              (1)
            tracer,                   & ! IN  - Cloud condensate amount in layer by type       ()
            Grid%xlat,                & ! IN  - Latitude                                       (radians)
            Grid%xlon,                & ! IN  - Longitude                                      (radians)
            Sfcprop%slmsk,            & ! IN  - Land/Sea mask                                  ()
            deltaZ,                   & ! IN  - Layer thickness                                (km)
            deltaP/100.,              & ! IN  - Layer thickness                                (hPa)
            Model%ntrac-1,            & ! IN  - Number of tracers
            Model%ntcw-1,             & ! IN  - Tracer index for cloud condensate (or liquid water)
            Model%ntiw-1,             & ! IN  - Tracer index for ice 
            Model%ntrw-1,             & ! IN  - Tracer index for rain 
            Model%ntsw-1,             & ! IN  - Tracer index for snow 
            Model%ntgl-1,             & ! IN  - Tracer index for groupel
            NCOL,                     & ! IN  - Number of horizontal gridpoints
            MODEL%LEVS,               & ! IN  - Number of model layers
            MODEL%LEVS+1,             & ! IN  - Number of model levels
            Model%uni_cld,            & ! IN  - True for cloud fraction from shoc
            Model%lmfshal,            & ! IN  - True for mass flux shallow convection
            Model%lmfdeep2,           & ! IN  - True for mass flux deep convection
            cldcov(:,1:Model%levs),   & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
            Tbd%phy_f3d(:,:,1),       & ! IN  - Liquid-water effective radius                  (microns)
            Tbd%phy_f3d(:,:,2),       & ! IN  - Ice-water effective radius                     (microns)
            Tbd%phy_f3d(:,:,3),       & ! IN  - LSnow-water effective radius                   (microns)
            clouds,                   & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
            cldsa,                    & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
            mtopa,                    & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
            mbota,                    & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
            de_lgth)                    ! OUT - clouds decorrelation length (km)
    endif ! end if_imp_physics
  end subroutine cloud_microphysics
  !
end module GFS_rrtmgp_pre

module GFS_rrtmgp_pre
  use machine, only: &
       kind_phys                   ! Working type
  use funcphys, only:            &
       fpvs                        ! Function ot compute sat. vapor pressure over liq.
  use module_radiation_gases,    only: &
       NF_VGAS,                  & ! Number of active gas species
       getgases,                 & ! Routine to setup trace gases
       getozn                      ! Routine to setup ozone
  ! RRTMGP types
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use rrtmgp_aux,            only: check_error_msg

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
  subroutine GFS_rrtmgp_pre_init(nGases, active_gases, active_gases_array, errmsg, errflg)
    ! Inputs
	integer, intent(in) :: &
	     nGases     ! Number of active gases in RRTMGP
	character(len=*), intent(in) :: &
	     active_gases ! List of active gases from namelist.     
    ! Outputs
    character(len=*),dimension(nGases), intent(out) :: &
         active_gases_array  ! Character array containing trace gases to include in RRTMGP
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer, intent(out) :: &  
         errflg     ! Error flag

    ! Local variables
    character(len=1) :: tempstr
    integer :: ij, count
    integer,dimension(nGases,2) :: gasIndices

    ! Initialize
    errmsg = ''
    errflg = 0

    if (len(active_gases) .eq. 0) return
    
    ! Which gases are active? Provided via physics namelist.

    ! Pull out gas names from list...
    ! First grab indices in character array corresponding to start:end of gas name.
    gasIndices(1,1)=1
    count=1
    do ij=1,len(active_gases)
       tempstr=trim(active_gases(ij:ij))
       if (tempstr .eq. '_') then
          gasIndices(count,2)=ij-1
          gasIndices(count+1,1)=ij+1
          count=count+1
       endif
    enddo
    gasIndices(nGases,2)=len(trim(active_gases))
    
    ! Now extract the gas names
    do ij=1,nGases
       active_gases_array(ij) = active_gases(gasIndices(ij,1):gasIndices(ij,2))
    enddo

    ! Which gases are active? (This is purely for flexibility)
    do ij=1,nGases
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
  subroutine GFS_rrtmgp_pre_run(nCol, nLev, nGases, nTracers, i_o3, lsswr, lslwr, fhswr, &
       fhlwr, xlat, xlon,  prsl, tgrs, prslk, prsi, qgrs, tsfc, active_gases_array,      &
       con_eps, con_epsm1, con_fvirt, con_epsqs,                                         &
       raddt, p_lay, t_lay, p_lev, t_lev, tsfg, tsfa, tv_lay, relhum, tracer,            &
       gas_concentrations,  errmsg, errflg)
    
    ! Inputs   
    integer, intent(in)    :: &
         nCol,              & ! Number of horizontal grid points
         nLev,              & ! Number of vertical layers
         nGases,            & ! Number of active gases in RRTMGP.
         nTracers,          & ! Number of tracers from model. 
         i_o3                 ! Index into tracer array for ozone
    logical, intent(in) :: &
    	 lsswr,             & ! Call SW radiation?
    	 lslwr                ! Call LW radiation
    character(len=*),dimension(nGases), intent(in) :: &
         active_gases_array   ! Character array containing trace gases to include in RRTMGP
    real(kind_phys), intent(in) :: &
         fhswr,             & ! Frequency of SW radiation call.
         fhlwr                ! Frequency of LW radiation call.
    real(kind_phys), intent(in) :: &
         con_eps,           & ! Physical constant: Epsilon (Rd/Rv)
         con_epsm1,         & ! Physical constant: Epsilon (Rd/Rv) minus one
         con_fvirt,         & ! Physical constant: Inverse of epsilon minus one
         con_epsqs            ! Physical constant: Minimum saturation mixing-ratio (kg/kg)
    real(kind_phys), dimension(nCol), intent(in) :: & 
    	 xlon,              & ! Longitude
    	 xlat,              & ! Latitude
    	 tsfc                 ! Surface skin temperature (K)
    real(kind_phys), dimension(nCol,nLev), intent(in) :: & 
         prsl,              & ! Pressure at model-layer centers (Pa)
         tgrs,              & ! Temperature at model-layer centers (K)
         prslk                ! Exner function at model layer centers (1)
    real(kind_phys), dimension(nCol,nLev+1) :: & 
         prsi                 ! Pressure at model-interfaces (Pa)
    real(kind_phys), dimension(nCol,nLev,nTracers) :: & 
         qgrs                 ! Tracer concentrations (kg/kg)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag    
    real(kind_phys), intent(out) :: &
         raddt                ! Radiation time-step
    real(kind_phys), dimension(ncol), intent(out) :: &
         tsfg,              & ! Ground temperature
         tsfa                 ! Skin temperature    
    real(kind_phys), dimension(nCol,nLev), intent(out) :: &
         p_lay,             & ! Pressure at model-layer
         t_lay,             & ! Temperature at model layer
         tv_lay,            & ! Virtual temperature at model-layers 
         relhum               ! Relative-humidity at model-layers          
    real(kind_phys), dimension(nCol,nLev+1), intent(out) :: &
         p_lev,             & ! Pressure at model-interface
         t_lev                ! Temperature at model-interface
    real(kind_phys), dimension(nCol, nLev, nTracers),intent(out) :: &
         tracer               ! Array containing trace gases
    type(ty_gas_concs),intent(out) :: &
         gas_concentrations   ! RRTMGP DDT: gas volumne mixing ratios
         
    ! Local variables
    integer :: i, j, iCol, iBand, iSFC, iTOA, iLay
    logical :: top_at_1
    real(kind_phys),dimension(nCol,nLev) :: vmr_o3, vmr_h2o
    real(kind_phys) :: es, qs, tem1, tem2
    real(kind_phys), dimension(nCol,nLev) :: o3_lay, q_lay
    real(kind_phys), dimension(nCol,nLev, NF_VGAS) :: gas_vmr

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. (lsswr .or. lslwr)) return
        
    ! #######################################################################################
    ! What is vertical ordering?
    ! #######################################################################################
    top_at_1 = (prsi(1,1) .lt.  prsi(1, nLev))
    if (top_at_1) then 
       iSFC = nLev
       iTOA = 1
    else
       iSFC = 1
       iTOA = nLev
    endif

    ! #######################################################################################
    ! Compute some fields needed by RRTMGP
    ! #######################################################################################
    
    ! Water-vapor mixing-ratio
    q_lay(1:ncol,:)  = qgrs(1:NCOL,:,1)
    where(q_lay .lt. 1.e-6) q_lay = 1.e-6
    
    ! Pressure at layer-interface
    p_lev(1:NCOL,:) = prsi(1:NCOL,:)

    ! Pressure at layer-center
    p_lay(1:NCOL,:)   = prsl(1:NCOL,:)

    ! Temperature at layer-center
    t_lay(1:NCOL,:) = tgrs(1:NCOL,:)

    ! Temperature at layer-interfaces
    if (top_at_1) then
       t_lev(1:NCOL,1)      = t_lay(1:NCOL,iTOA)
       t_lev(1:NCOL,2:iSFC) = (t_lay(1:NCOL,2:iSFC)+t_lay(1:NCOL,1:iSFC-1))/2._kind_phys
       t_lev(1:NCOL,iSFC+1) = tsfc(1:NCOL)
    else
       t_lev(1:NCOL,1)      = tsfc(1:NCOL)
       t_lev(1:NCOL,2:iTOA) = (t_lay(1:NCOL,2:iTOA)+t_lay(1:NCOL,1:iTOA-1))/2._kind_phys
       t_lev(1:NCOL,iTOA+1) = t_lay(1:NCOL,iTOA)
    endif

    ! Compute a bunch of thermodynamic fields needed by the cloud microphysics schemes. 
    ! Relative humidity, saturation mixing-ratio, vapor mixing-ratio, virtual temperature, 
    ! layer thickness,...
    do iCol=1,NCOL
       do iLay=1,nLev
          es                = min( p_lay(iCol,iLay),  fpvs( t_lay(iCol,iLay) ) )  ! fpvs and prsl in pa
          qs                = max( con_epsqs, con_eps * es / (p_lay(iCol,iLay) + con_epsm1*es) )
          relhum(iCol,iLay) = max( 0._kind_phys, min( 1._kind_phys, max(con_epsqs, q_lay(iCol,iLay))/qs ) )
          tv_lay(iCol,iLay) = t_lay(iCol,iLay) * (1._kind_phys + con_fvirt*q_lay(iCol,iLay)) 
       enddo
    enddo

    ! #######################################################################################
    ! Get layer ozone mass mixing ratio 
    ! #######################################################################################
    ! First recast remaining all tracers (except sphum) forcing them all to be positive
    do j = 2, nTracers
       tracer(1:NCOL,:,j) = qgrs(1:NCOL,:,j)
       where(tracer(:,:,j) .lt. 0.0) tracer(:,:,j) = 0._kind_phys
    enddo

    if (i_o3 > 0) then 
       do iLay=1,nlev
          do iCol=1,NCOL
             o3_lay(iCol,iLay) = max( con_epsqs, tracer(iCol,iLay,i_o3) )
          enddo
       enddo
    ! OR Use climatological ozone data
    else                               
       call getozn (prslk(1:NCOL,:), xlat, nCol, nLev, o3_lay)
    endif

    ! #######################################################################################
    ! Set gas concentrations for RRTMGP
    ! #######################################################################################
    ! Call getgases(), to set up non-prognostic gas volume mixing ratios (gas_vmr).
    call getgases (p_lev/100., xlon, xlat, nCol, nLev, gas_vmr)

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
    ! Radiation time step (output) (Is this really needed?) (Used by some diagnostics)
    ! #######################################################################################
    raddt = min(fhswr, fhlwr)

    ! #######################################################################################
    ! Setup surface ground temperature and ground/air skin temperature if required.
    ! #######################################################################################
    tsfg(1:NCOL) = tsfc(1:NCOL)
    tsfa(1:NCOL) = tsfc(1:NCOL)

  end subroutine GFS_rrtmgp_pre_run
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_pre_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_pre_finalize ()
  end subroutine GFS_rrtmgp_pre_finalize

end module GFS_rrtmgp_pre

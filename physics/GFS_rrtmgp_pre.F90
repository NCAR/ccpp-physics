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
  use mo_gas_concentrations, only: ty_gas_concs
  use rrtmgp_aux,            only: check_error_msg

  real(kind_phys), parameter :: &
       amd   = 28.9644_kind_phys,  & ! Molecular weight of dry-air     (g/mol)
       amw   = 18.0154_kind_phys,  & ! Molecular weight of water vapor (g/mol)
       amo3  = 47.9982_kind_phys,  & ! Modelular weight of ozone       (g/mol)
       amdw  = amd/amw,            & ! Molecular weight of dry air / water vapor
       amdo3 = amd/amo3              ! Molecular weight of dry air / ozone

  ! Save trace gas indices.
  integer :: iStr_h2o, iStr_co2, iStr_o3, iStr_n2o, iStr_ch4, iStr_o2, iStr_ccl4, &
       iStr_cfc11, iStr_cfc12, iStr_cfc22 
    character(len=32),dimension(:),allocatable :: &
         active_gases_array 

  public GFS_rrtmgp_pre_run,GFS_rrtmgp_pre_init,GFS_rrtmgp_pre_finalize  
contains
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_pre_init
  ! #########################################################################################
!! \section arg_table_GFS_rrtmgp_pre_init
!! \htmlinclude GFS_rrtmgp_pre_init.html
!!
  subroutine GFS_rrtmgp_pre_init(nGases, active_gases, errmsg, errflg)
    ! Inputs
    integer, intent(in) :: &
         nGases       ! Number of active gases in RRTMGP
    character(len=*), intent(in) :: &
         active_gases ! List of active gases from namelist.     
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg             ! Error message
    integer, intent(out) :: &  
         errflg             ! Error flag

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
    allocate(active_gases_array(nGases))
    do ij=1,nGases
       active_gases_array(ij) = active_gases(gasIndices(ij,1):gasIndices(ij,2))
       if(trim(active_gases_array(ij)) .eq. 'h2o')   istr_h2o       = ij
       if(trim(active_gases_array(ij)) .eq. 'co2')   istr_co2       = ij
       if(trim(active_gases_array(ij)) .eq. 'o3')    istr_o3        = ij
       if(trim(active_gases_array(ij)) .eq. 'n2o')   istr_n2o       = ij
       if(trim(active_gases_array(ij)) .eq. 'ch4')   istr_ch4       = ij
       if(trim(active_gases_array(ij)) .eq. 'o2')    istr_o2        = ij
       if(trim(active_gases_array(ij)) .eq. 'ccl4')  istr_ccl4      = ij
       if(trim(active_gases_array(ij)) .eq. 'cfc11') istr_cfc11     = ij
       if(trim(active_gases_array(ij)) .eq. 'cfc12') istr_cfc12     = ij
       if(trim(active_gases_array(ij)) .eq. 'cfc22') istr_cfc22     = ij
    enddo

  end subroutine GFS_rrtmgp_pre_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_pre_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_pre_run
!! \htmlinclude GFS_rrtmgp_pre_run.html
!!
  subroutine GFS_rrtmgp_pre_run(nCol, nLev, nTracers, i_o3, lsswr, lslwr, fhswr, fhlwr,     &
       xlat, xlon,  prsl, tgrs, prslk, prsi, qgrs, tsfc, con_eps, con_epsm1, con_fvirt,     &
       con_epsqs, minGPpres, minGPtemp, raddt, p_lay, t_lay, p_lev, t_lev, tsfg, tsfa,      & 
       qs_lay, q_lay, tv_lay, relhum, tracer, gas_concentrations, errmsg, errflg)
    
    ! Inputs   
    integer, intent(in)    :: &
         nCol,              & ! Number of horizontal grid points
         nLev,              & ! Number of vertical layers
         nTracers,          & ! Number of tracers from model. 
         i_o3                 ! Index into tracer array for ozone
    logical, intent(in) :: &
    	 lsswr,             & ! Call SW radiation?
    	 lslwr                ! Call LW radiation
    real(kind_phys), intent(in) :: &
         minGPtemp,         & ! Minimum temperature allowed in RRTMGP.
         minGPpres,         & ! Minimum pressure allowed in RRTMGP.
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
    real(kind_phys), intent(inout) :: &
         raddt                ! Radiation time-step
    real(kind_phys), dimension(ncol), intent(inout) :: &
         tsfg,              & ! Ground temperature
         tsfa                 ! Skin temperature    
    real(kind_phys), dimension(nCol,nLev), intent(inout) :: &
         p_lay,             & ! Pressure at model-layer
         t_lay,             & ! Temperature at model layer
         q_lay,             & ! Water-vapor mixing ratio (kg/kg)
         tv_lay,            & ! Virtual temperature at model-layers 
         relhum,            & ! Relative-humidity at model-layers   
         qs_lay               ! Saturation vapor pressure at model-layers
    real(kind_phys), dimension(nCol,nLev+1), intent(inout) :: &
         p_lev,             & ! Pressure at model-interface
         t_lev                ! Temperature at model-interface
    real(kind_phys), dimension(nCol, nLev, nTracers),intent(inout) :: &
         tracer               ! Array containing trace gases
    type(ty_gas_concs),intent(inout) :: &
         gas_concentrations   ! RRTMGP DDT: gas volumne mixing ratios
         
    ! Local variables
    integer :: i, j, iCol, iBand, iSFC, iTOA, iLay
    logical :: top_at_1
    real(kind_phys),dimension(nCol,nLev) :: vmr_o3, vmr_h2o
    real(kind_phys) :: es, tem1, tem2
    real(kind_phys), dimension(nCol,nLev) :: o3_lay, tem2da, tem2db
    real(kind_phys), dimension(nCol,nLev, NF_VGAS) :: gas_vmr
    character(len=32), dimension(gas_concentrations%get_num_gases()) :: active_gases

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

    ! Bound temperature at layer centers.
    do iCol=1,NCOL
       do iLay=1,nLev
          if (t_lay(iCol,iLay) .le. minGPtemp) then
             t_lay = minGPtemp + epsilon(minGPtemp)
          endif
       enddo
    enddo

    ! Temperature at layer-interfaces          
    if (top_at_1) then
       tem2da(1:nCol,2:iSFC) = log(p_lay(1:nCol,2:iSFC))
       tem2db(1:nCol,2:iSFC) = log(p_lev(1:nCol,2:iSFC)) 
       do iCol = 1, nCol
           tem2da(iCol,1)    = log(p_lay(iCol,1) )
           tem2db(iCol,1)    = log(max(minGPpres, p_lev(iCol,1)) )
           tem2db(iCol,iSFC) = log(p_lev(iCol,iSFC) )    
       enddo
       !
       t_lev(1:NCOL,1)      = t_lay(1:NCOL,iTOA)
       do iLay = 2, iSFC
          do iCol = 1, nCol
            t_lev(iCol,iLay) = t_lay(iCol,iLay) + (t_lay(iCol,iLay-1) - t_lay(iCol,iLay))&
                     * (tem2db(iCol,iLay)   - tem2da(iCol,iLay))                   &
                     / (tem2da(iCol,iLay-1) - tem2da(iCol,iLay))
           enddo
        enddo
       t_lev(1:NCOL,iSFC+1) = tsfc(1:NCOL)
    else
       tem2da(1:nCol,2:iTOA) = log(p_lay(1:nCol,2:iTOA))
       tem2db(1:nCol,2:iTOA) = log(p_lev(1:nCol,2:iTOA))     
       do iCol = 1, nCol
           tem2da(iCol,1)    = log(p_lay(iCol,1))
           tem2db(iCol,1)    = log(p_lev(iCol,1))    
           tem2db(iCol,iTOA) = log(max(minGPpres, p_lev(iCol,iTOA)) )
       enddo    
       !
       t_lev(1:NCOL,1)      = tsfc(1:NCOL)
       do iLay = 1, iTOA-1
          do iCol = 1, nCol
            t_lev(iCol,iLay+1) = t_lay(iCol,iLay) + (t_lay(iCol,iLay+1) - t_lay(iCol,iLay))&
                     * (tem2db(iCol,iLay+1) - tem2da(iCol,iLay))                   &
                     / (tem2da(iCol,iLay+1) - tem2da(iCol,iLay))
           enddo
        enddo
       t_lev(1:NCOL,iTOA+1) = t_lay(1:NCOL,iTOA)
    endif

    ! Compute a bunch of thermodynamic fields needed by the cloud microphysics schemes. 
    ! Relative humidity, saturation mixing-ratio, vapor mixing-ratio, virtual temperature, 
    ! layer thickness,...
    do iCol=1,NCOL
       do iLay=1,nLev
          es                = min( p_lay(iCol,iLay),  fpvs( t_lay(iCol,iLay) ) )  ! fpvs and prsl in pa
          qs_lay(iCol,iLay) = max( con_epsqs, con_eps * es / (p_lay(iCol,iLay) + con_epsm1*es) )
          relhum(iCol,iLay) = max( 0._kind_phys, min( 1._kind_phys, max(con_epsqs, q_lay(iCol,iLay))/qs_lay(iCol,iLay) ) )
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
    
    ! Populate RRTMGP DDT w/ gas-concentrations
    gas_concentrations%gas_name(:)                = active_gases_array(:)
    gas_concentrations%concs(istr_o2)%conc(:,:)   = gas_vmr(:,:,4)
    gas_concentrations%concs(istr_co2)%conc(:,:)  = gas_vmr(:,:,1)
    gas_concentrations%concs(istr_ch4)%conc(:,:)  = gas_vmr(:,:,3)
    gas_concentrations%concs(istr_n2o)%conc(:,:)  = gas_vmr(:,:,2)
    gas_concentrations%concs(istr_h2o)%conc(:,:)  = vmr_h2o(:,:)
    gas_concentrations%concs(istr_o3)%conc(:,:)   = vmr_o3(:,:)

    ! #######################################################################################
    ! Radiation time step (output) (Is this really needed?) (Used by some diagnostics)
    ! #######################################################################################
    raddt = min(fhswr, fhlwr)

    ! #######################################################################################
    ! Setup surface ground temperature and ground/air skin temperature if required.
    ! #######################################################################################
    tsfg(1:NCOL) = tsfc(1:NCOL)
    tsfa(1:NCOL) = t_lay(1:NCOL,iSFC)

  end subroutine GFS_rrtmgp_pre_run
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_pre_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_pre_finalize ()
  end subroutine GFS_rrtmgp_pre_finalize

end module GFS_rrtmgp_pre

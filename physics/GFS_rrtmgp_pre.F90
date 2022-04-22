module GFS_rrtmgp_pre
  use machine, only: &
       kind_phys                   ! Working type
  use funcphys, only:            &
       fpvs                        ! Function ot compute sat. vapor pressure over liq.
  use module_radiation_astronomy, only: &
       coszmn 
  use module_radiation_gases,    only: &
       NF_VGAS,                  & ! Number of active gas species
       getgases,                 & ! Routine to setup trace gases
       getozn                      ! Routine to setup ozone
  ! RRTMGP types
  use mo_gas_concentrations, only: ty_gas_concs
  use radiation_tools,       only: check_error_msg,cmp_tlev

  real(kind_phys), parameter :: &
       amd   = 28.9644_kind_phys,  & ! Molecular weight of dry-air     (g/mol)
       amw   = 18.0154_kind_phys,  & ! Molecular weight of water vapor (g/mol)
       amo3  = 47.9982_kind_phys,  & ! Modelular weight of ozone       (g/mol)
       amdw  = amd/amw,            & ! Molecular weight of dry air / water vapor
       amdo3 = amd/amo3              ! Molecular weight of dry air / ozone

  ! Save trace gas indices.
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
         nGases       ! Number of active gases in RRTMGP
    character(len=*), intent(in) :: &
         active_gases ! List of active gases from namelist
    character(len=*), dimension(:), intent(out) :: &
         active_gases_array ! List of active gases from namelist as array

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
  subroutine GFS_rrtmgp_pre_run(me, nCol, nLev, nTracers, i_o3, lsswr, lslwr, fhswr, fhlwr, &
       xlat, xlon,  prsl, tgrs, prslk, prsi, qgrs, tsfc, coslat, sinlat, con_g, con_rd,     &
       con_eps, con_epsm1, con_fvirt, con_epsqs, solhr, minGPpres, maxGPpres, minGPtemp,    &
       maxGPtemp, raddt, p_lay, t_lay, p_lev, t_lev, tsfg, tsfa, qs_lay, q_lay, tv_lay,     &
       relhum, tracer, deltaZ, deltaZc, deltaP, active_gases_array, gas_concentrations,     &
       tsfc_radtime, coszen, coszdg, top_at_1, iSFC, iTOA, errmsg, errflg)
    
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
         maxGPtemp,         & ! Maximum ...
         minGPpres,         & ! Minimum pressure allowed in RRTMGP.
         maxGPpres,         & ! Maximum pressure allowed in RRTMGP. 
         fhswr,             & ! Frequency of SW radiation call.
         fhlwr                ! Frequency of LW radiation call.
    real(kind_phys), intent(in) :: &
         con_g,             & ! Physical constant: gravitational constant
         con_rd,            & ! Physical constant: gas-constant for dry air
         con_eps,           & ! Physical constant: Epsilon (Rd/Rv)
         con_epsm1,         & ! Physical constant: Epsilon (Rd/Rv) minus one
         con_fvirt,         & ! Physical constant: Inverse of epsilon minus one
         con_epsqs,         & ! Physical constant: Minimum saturation mixing-ratio (kg/kg)
         solhr                ! Time in hours after 00z at the current timestep 
    real(kind_phys), dimension(:), intent(in) :: & 
    	 xlon,              & ! Longitude
    	 xlat,              & ! Latitude
    	 tsfc,              & ! Surface skin temperature (K)
         coslat,            & ! Cosine(latitude)
         sinlat               ! Sine(latitude) 
    real(kind_phys), dimension(:,:), intent(in) :: & 
         prsl,              & ! Pressure at model-layer centers (Pa)
         tgrs,              & ! Temperature at model-layer centers (K)
         prslk,             & ! Exner function at model layer centers (1)
         prsi                 ! Pressure at model-interfaces (Pa)
    real(kind_phys), dimension(:,:,:), intent(in) :: & 
         qgrs                 ! Tracer concentrations (kg/kg)
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array   ! List of active gases from namelist as array

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg,            & ! Error flag
         iSFC,              & ! Vertical index for surface
         iTOA                 ! Vertical index for TOA
    logical, intent(out) :: &
         top_at_1             ! Vertical ordering flag
    real(kind_phys), intent(inout) :: &
         raddt                ! Radiation time-step
    real(kind_phys), dimension(:), intent(inout) :: &
         tsfg,              & ! Ground temperature
         tsfa,              & ! Skin temperature    
         tsfc_radtime,      & ! Surface temperature at radiation timestep
         coszen,            & ! Cosine of SZA
         coszdg               ! Cosine of SZA, daytime
    real(kind_phys), dimension(:,:), intent(inout) :: &
         p_lay,             & ! Pressure at model-layer
         t_lay,             & ! Temperature at model layer
         q_lay,             & ! Water-vapor mixing ratio (kg/kg)
         tv_lay,            & ! Virtual temperature at model-layers 
         relhum,            & ! Relative-humidity at model-layers   
         qs_lay,            & ! Saturation vapor pressure at model-layers
         deltaZ,            & ! Layer thickness (m)
         deltaZc,           & ! Layer thickness (m) (between layer centers)
         deltaP,            & ! Layer thickness (Pa)
         p_lev,             & ! Pressure at model-interface
         t_lev                ! Temperature at model-interface
    real(kind_phys), dimension(:,:,:),intent(inout) :: &
         tracer               ! Array containing trace gases
    type(ty_gas_concs), intent(inout) :: &
         gas_concentrations   ! RRTMGP DDT: gas volumne mixing ratios

    ! Local variables
    integer :: i, j, iCol, iBand, iLay, iLev, iSFC_ilev
    real(kind_phys),dimension(nCol,nLev) :: vmr_o3, vmr_h2o
    real(kind_phys) :: es, tem1, tem2, pfac
    real(kind_phys), dimension(nLev+1) :: hgtb
    real(kind_phys), dimension(nLev)   :: hgtc
    real(kind_phys), dimension(nCol,nLev) :: o3_lay
    real(kind_phys), dimension(nCol,nLev, NF_VGAS) :: gas_vmr
    real(kind_phys) :: con_rdog

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
       iSFC_ilev = iSFC + 1
    else
       iSFC = 1
       iTOA = nLev
       iSFC_ilev = 1
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
    p_lay(1:NCOL,:) = prsl(1:NCOL,:)

    ! Temperature at layer-center
    t_lay(1:NCOL,:) = tgrs(1:NCOL,:)

    ! Bound temperature/pressure at layer centers.
    do iLay=1,nLev
       do iCol=1,NCOL
          if (t_lay(iCol,iLay) .le. minGPtemp) then
             t_lay(iCol,iLay) = minGPtemp + epsilon(minGPtemp)
          endif
          if (p_lay(iCol,iLay) .le. minGPpres) then
             p_lay(iCol,iLay) = minGPpres + epsilon(minGPpres)
          endif
          if (t_lay(iCol,iLay) .ge. maxGPtemp) then
             t_lay(iCol,iLay) = maxGPtemp - epsilon(maxGPtemp)
          endif
          if (p_lay(iCol,iLay) .ge. maxGPpres) then
             p_lay(iCol,iLay) = maxGPpres - epsilon(maxGPpres)
          endif
       enddo
    enddo

    ! Temperature at layer-interfaces          
    call cmp_tlev(nCol,nLev,minGPpres,p_lay,t_lay,p_lev,tsfc,t_lev)
    do iLev=1,nLev+1
       do iCol=1,nCol
          if (t_lev(iCol,iLev) .le. minGPtemp) t_lev(iCol,iLev) = minGPtemp + epsilon(minGPtemp)
          if (t_lev(iCol,iLev) .ge. maxGPtemp) t_lev(iCol,iLev) = maxGPtemp - epsilon(maxGPtemp)
       enddo
    enddo

    ! Save surface temperature at radiation time-step, used for LW flux adjustment betwen
    ! radiation calls.
    tsfc_radtime = tsfc

    ! Compute a bunch of thermodynamic fields needed by the cloud microphysics schemes. 
    ! Relative humidity, saturation mixing-ratio, vapor mixing-ratio, virtual temperature, 
    ! layer thickness,...
    do iLay=1,nLev
       do iCol=1,NCOL
          es                = min( p_lay(iCol,iLay),  fpvs( t_lay(iCol,iLay) ) )  ! fpvs and prsl in pa
          qs_lay(iCol,iLay) = max( con_epsqs, con_eps * es / (p_lay(iCol,iLay) + con_epsm1*es) )
          relhum(iCol,iLay) = max( 0._kind_phys, min( 1._kind_phys, max(con_epsqs, q_lay(iCol,iLay))/qs_lay(iCol,iLay) ) )
          tv_lay(iCol,iLay) = t_lay(iCol,iLay) * (1._kind_phys + con_fvirt*q_lay(iCol,iLay)) 
       enddo
    enddo

    !
    ! Compute layer-thickness between layer boundaries (deltaZ) and layer centers (deltaZc)
    !
    deltaP = abs(p_lev(:,2:nLev+1)-p_lev(:,1:nLev))
    con_rdog = con_rd/con_g
    do iCol=1,nCol 
       if (top_at_1) then
          ! Layer thickness (m)
          do iLay=1,nLev
             deltaZ(iCol,iLay) = con_rdog * abs(log(p_lev(iCol,iLay+1)) - log(p_lev(iCol,iLay))) * tv_lay(iCol,iLay)
          enddo
          ! Height at layer boundaries
          hgtb(nLev+1) = 0._kind_phys
          do iLay=nLev,1,-1
             hgtb(iLay)= hgtb(iLay+1) + deltaZ(iCol,iLay)
          enddo
          ! Height at layer centers
          do iLay = nLev, 1, -1
             pfac = abs(log(p_lev(iCol,iLay+1)) - log(p_lay(iCol,iLay))) /  &
                  abs(log(p_lev(iCol,iLay+1)) - log(p_lev(iCol,iLay)))
             hgtc(iLay) = hgtb(iLay+1) + pfac * (hgtb(iLay) - hgtb(iLay+1))
          enddo
          ! Layer thickness between centers
          do iLay = nLev-1, 1, -1
             deltaZc(iCol,iLay) = hgtc(iLay) - hgtc(iLay+1)
          enddo
          deltaZc(iCol,nLev) = hgtc(nLev) - hgtb(nLev+1)
       else
          ! Layer thickness (m)
          do iLay=nLev,1,-1
             deltaZ(iCol,iLay) = con_rdog * abs(log(p_lev(iCol,iLay))  - log(p_lev(iCol,iLay+1))) * tv_lay(iCol,iLay)
          enddo
          ! Height at layer boundaries
          hgtb(1) = 0._kind_phys
          do iLay=1,nLev
             hgtb(iLay+1)= hgtb(iLay) + deltaZ(iCol,iLay)
          enddo
          ! Height at layer centers
          do iLay = 1, nLev
             pfac = abs(log(p_lev(iCol,iLay)) - log(p_lay(iCol,iLay)  )) /  &
                  abs(log(p_lev(iCol,iLay)) - log(p_lev(iCol,iLay+1)))
             hgtc(iLay) = hgtb(iLay) + pfac * (hgtb(iLay+1) - hgtb(iLay))
          enddo
          ! Layer thickness between centers
          do iLay = 2, nLev
             deltaZc(iCol,iLay) = hgtc(iLay) - hgtc(iLay-1)
          enddo
          deltaZc(iCol,1) = hgtc(1) - hgtb(1)
       endif
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
    gas_concentrations%ncol                       = nCol
    gas_concentrations%nlay                       = nLev
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
    tsfg(1:NCOL) = t_lev(1:NCOL,iSFC_ilev)
    tsfa(1:NCOL) = t_lay(1:NCOL,iSFC)

    ! #######################################################################################
    ! Compute cosine of zenith angle (only when SW is called)
    ! #######################################################################################
    if (lsswr) then
       call coszmn (xlon, sinlat, coslat, solhr, nCol, me, coszen, coszdg)
    endif

  end subroutine GFS_rrtmgp_pre_run
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_pre_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_pre_finalize ()
  end subroutine GFS_rrtmgp_pre_finalize
end module GFS_rrtmgp_pre

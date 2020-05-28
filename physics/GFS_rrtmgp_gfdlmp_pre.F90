! ########################################################################################
! This module contains the interface between the GFDL macrophysics and the RRTMGP radiation
! schemes. Only compatable with Model%imp_physics = Model%imp_physics_gfdl
! ########################################################################################
module GFS_rrtmgp_gfdlmp_pre
  use machine,      only: kind_phys
  use GFS_typedefs, only: GFS_control_type, GFS_tbd_type 
  use physcons,     only: con_ttp, & ! Temperature at h2o 3pt (K)
                          con_rd,  & ! Gas constant for dry air (J/KgK)
                          con_pi,  & ! PI
                          con_g      ! Gravity (m/s2)
  use physparam,    only: lcnorm,lcrick
  use rrtmgp_aux,   only: check_error_msg
  
  ! Parameters
  real(kind_phys), parameter :: &
       reice_min = 10.0,     & ! Minimum ice size allowed by scheme
       reice_max = 150.0,    & ! Maximum ice size allowed by scheme
       epsq      = 1.0e-12,  & ! Tiny value
       cllimit   = 0.001,    & ! Lowest cloud fraction in GFDL MP scheme
       gfac      = 1.0e5/con_g
       
contains  
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_gfdlmp_pre_init()
  end subroutine GFS_rrtmgp_gfdlmp_pre_init

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_GFS_rrtmgp_gfdlmp_pre_run
!! \htmlinclude GFS_rrtmgp_gfdlmp_pre_run.html
!!  
  subroutine GFS_rrtmgp_gfdlmp_pre_run(Model, Tbd, nCol, nLev, p_lev, tracer,            &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,   &
       cld_rerain, precip_frac, errmsg, errflg)
    implicit none
    
    ! Inputs  
    type(GFS_control_type), intent(in) :: &
         Model                ! DDT: FV3-GFS model control parameters      
    type(GFS_tbd_type), intent(in) :: &
         Tbd                  ! DDT: FV3-GFS data not yet assigned to a defined container         
    integer, intent(in) :: &
         nCol,              & ! Number of horizontal grid-points
         nLev                 ! Number of vertical-layers
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &         
         p_lev                ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(nCol, nLev, Model%ntrac),intent(in) :: &
         tracer               ! Cloud condensate amount in layer by type ()         
    
    ! Outputs     
    real(kind_phys), dimension(nCol,nLev),intent(out) :: &
         cld_frac,          & ! Total cloud fraction
         cld_lwp,           & ! Cloud liquid water path
         cld_reliq,         & ! Cloud liquid effective radius
         cld_iwp,           & ! Cloud ice water path
         cld_reice,         & ! Cloud ice effecive radius
         cld_swp,           & ! Cloud snow water path
         cld_resnow,        & ! Cloud snow effective radius
         cld_rwp,           & ! Cloud rain water path
         cld_rerain,        & ! Cloud rain effective radius       
         precip_frac          ! Precipitation fraction      
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag
    
    ! Local variables
    real(kind_phys) :: tem1
    real(kind_phys), dimension(nCol, nLev, min(4,Model%ncnd)) :: cld_condensate
    integer :: i,k,l,ncndl
    real(kind_phys), dimension(nCol,nLev) :: deltaP
    
    if (.not. (Model%lsswr .or. Model%lslwr)) return
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    ! Test inputs
    if (Model%ncnd .ne. 5) then
       errmsg = 'Incorrect number of cloud condensates provided'
       errflg = 1
       call check_error_msg('GFS_rrtmgp_gfdlmp_pre_run',errmsg)
       return
    endif
    ! 
    if (lcrick) then
       errmsg = 'Namelist option lcrick is not supported.'
       errflg = 1
       call check_error_msg('GFS_rrtmgp_gfdlmp_pre_run',errmsg)
       return
    endif
    ! 
    if (lcnorm) then
       errmsg = 'Namelist option lcnorm is not supported.'
       errflg = 1
       call check_error_msg('GFS_rrtmgp_gfdlmp_pre_run',errmsg)
       return
    endif
    !
    if (.not. Model%lgfdlmprad) then
       errmsg = 'Namelist option gfdlmprad=F is not supported.'
       errflg = 1
       call check_error_msg('GFS_rrtmgp_gfdlmp_pre_run',errmsg)
       return
    endif
    !
    if(.not. Model%effr_in) then
       errmsg = 'Namelist option effr_in=F is not supported.'
       errflg = 1
       call check_error_msg('GFS_rrtmgp_gfdlmp_pre_run',errmsg)
       return
    endif    
    
    ! Initialize outputs
    cld_lwp(:,:)    = 0.0
    cld_reliq(:,:)  = 0.0
    cld_iwp(:,:)    = 0.0
    cld_reice(:,:)  = 0.0
    cld_rwp(:,:)    = 0.0
    cld_rerain(:,:) = 0.0
    cld_swp(:,:)    = 0.0
    cld_resnow(:,:) = 0.0
    
    ! ####################################################################################
    ! Pull out cloud information for GFDL MP scheme.
    ! ####################################################################################
    ! Condensate
    cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,Model%ntcw)     ! -liquid water
    cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,Model%ntiw)     ! -ice water
    cld_condensate(1:nCol,1:nLev,3) = tracer(1:nCol,1:nLev,Model%ntrw)     ! -rain water
    cld_condensate(1:nCol,1:nLev,4) = tracer(1:nCol,1:nLev,Model%ntsw) + & ! -snow + grapuel
                                      tracer(1:nCol,1:nLev,Model%ntgl) 
                           
    ! Since we combine the snow and grapuel, define local variable for number of condensate types.
    ncndl = min(4,Model%ncnd)

    ! Set really tiny suspended particle amounts to clear
    do l=1,ncndl
       do k=1,nLev
          do i=1,nCol   
             if (cld_condensate(i,k,l) < epsq) cld_condensate(i,k,l) = 0.0
          enddo
       enddo
    enddo
    
    ! Cloud-fraction
    cld_frac(1:nCol,1:nLev)   = tracer(1:nCol,1:nLev,Model%ntclamt)
    
    ! Precipitation fraction (Hack. For now use cloud-fraction)
    precip_frac(1:nCol,1:nLev) = tracer(1:nCol,1:nLev,Model%ntclamt)
    
    ! Condensate and effective size
    deltaP = abs(p_lev(:,2:nLev+1)-p_lev(:,1:nLev))/100.  
    do k = 1, nLev
       do i = 1, nCol
          ! Compute liquid/ice condensate path from mixing ratios (kg/kg)->(g/m2)   
          if (cld_frac(i,k) .ge. cllimit) then         
             tem1          = gfac * deltaP(i,k)
             cld_lwp(i,k)  = cld_condensate(i,k,1) * tem1
             cld_iwp(i,k)  = cld_condensate(i,k,2) * tem1
             cld_rwp(i,k)  = cld_condensate(i,k,3) * tem1
             cld_swp(i,k)  = cld_condensate(i,k,4) * tem1  
          endif
          ! Use radii provided from the macrophysics        
          cld_reliq(i,k)  = Tbd%phy_f3d(i,k,1)
          cld_reice(i,k)  = max(reice_min, min(reice_max,Tbd%phy_f3d(i,k,2)))
          cld_rerain(i,k) = Tbd%phy_f3d(i,k,3)
          cld_resnow(i,k) = Tbd%phy_f3d(i,k,4)            
       enddo
    enddo
    
  end subroutine GFS_rrtmgp_gfdlmp_pre_run
  
  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_gfdlmp_pre_finalize()
  end subroutine GFS_rrtmgp_gfdlmp_pre_finalize
  
end module GFS_rrtmgp_gfdlmp_pre

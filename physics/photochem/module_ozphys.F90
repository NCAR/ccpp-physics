! #########################################################################################
!> \section arg_table_module_ozphys Argument table                               
!! \htmlinclude module_ozphys.html                                               
!!
!
!> The operational GFS currently parameterizes ozone production and destruction based on 
!! monthly mean coefficients (\c global_o3prdlos.f77) provided by Naval Research Laboratory
!! through CHEM2D chemistry model (McCormack et al. (2006) \cite mccormack_et_al_2006).
!!
!! There are two implementations of this parameterization within this module.
!! run_o3prog_2006 - Relies on either two/four mean monthly coefficients. This is explained
!!                   in (https://doi.org/10.5194/acp-6-4943-2006. See Eq.(4)).
!! run_o3prog_2015 - Relies on six mean monthly coefficients, specifically for NRL 
!!                   parameterization and climatological T and O3 are in location 5 and 6 of
!!                   the coefficient array.
!! 
!! Both of these rely on the scheme being setup correctly by invoking the load(), setup(), 
!! and update() procedures prior to calling the run() procedure.
!!
!! load_o3prog()   - Read in data and load into type ty_ozphys (called once from host)
!! setup_o3prog()  - Create spatial interpolation indices      (called once, after model grid is known)
!! update_o3prog() - Update ozone concentration in time        (call in physics loop, before run())
!!                   *CAVEAT* Since the radiation is often run at a lower temporal resolution
!!                            than the rest of the physics, update_o3prog() needs to be
!!                            called before the radiation, which is called before the physics.
!!                            For example, within the physics loop:
!!                                update_o3prog() -> radiation() -> run_o3prog() -> physics....
!!
!! Additionally, there is the functionality to not use interactive ozone, instead reverting
!! to ozone climatology. In this case, analagous to when using prognostic ozone, there are
!! update() and run() procedures that need to be called before the radiation.
!! For example, within the physics loop:
!!     update_o3clim() -> run_o3clim() -> radiation() -> physics...
!!
!!\author   June 2015 - Shrinivas Moorthi
!!\modified Sep  2023 - Dustin Swales
!!
! #########################################################################################
module module_ozphys
  use machine,  only : kind_phys
  use funcphys, only : fpkapx
  implicit none

  public ty_ozphys

! ######################################################################################### 
!> \section arg_table_ty_ozphys Argument Table 
!! \htmlinclude ty_ozphys.html
!!
!> Derived type containing data and procedures needed by ozone photochemistry parameterization
!! *Note* All data field are ordered from surface-to-toa.
!!
! #########################################################################################
  type ty_ozphys
     ! Prognostic ozone.
     integer                      :: nlat          !< Number of latitudes.
     integer                      :: nlev          !< Number of vertical layers.
     integer                      :: ntime         !< Number of times.
     integer                      :: ncf           !< Number of coefficients.
     real(kind_phys), allocatable :: lat(:)        !< Latitude.
     real(kind_phys), allocatable :: pres(:)       !< Pressure levels.
     real(kind_phys), allocatable :: po3(:)        !< Natural log pressure of levels.
     real(kind_phys), allocatable :: time(:)       !< Time.
     real(kind_phys), allocatable :: data(:,:,:,:) !< Ozone forcing data (raw)
     ! Climotological ozone.
     integer                      :: nlatc         !< Number of latitudes.
     integer                      :: nlevc         !< Number of vertical layers.
     integer                      :: ntimec        !< Number of times.
     real(kind_phys)              :: blatc         !< Parameter for ozone climotology
     real(kind_phys)              :: dphiozc       !< Parameter for ozone climotology
     real(kind_phys), allocatable :: pkstr(:)      !<
     real(kind_phys), allocatable :: pstr(:)       !<
     real(kind_phys), allocatable :: datac(:,:,:)  !< Ozone climotological data
     integer                      :: k1oz          !< Lower interpolation index in datac(dim=3), time dim
     integer                      :: k2oz          !< Upper interpolation index in datac(dim=3), time dim 
     real(kind_phys)              :: facoz         !< Parameter for ozone climotology
     contains
       procedure, public :: load_o3prog
       procedure, public :: setup_o3prog
       procedure, public :: update_o3prog
       procedure, public :: run_o3prog_2015
       procedure, public :: run_o3prog_2006
       !
       procedure, public :: load_o3clim
       procedure, public :: update_o3clim
       procedure, public :: run_o3clim
  end type ty_ozphys
  
contains
  ! #########################################################################################
  ! Procedure (type-bound) for loading data for prognostic ozone.
  ! #########################################################################################
  function load_o3prog(this, file, fileID) result (err_message)
    class(ty_ozphys), intent(inout) :: this
    integer,          intent(in)    :: fileID
    character(len=*), intent(in)    :: file
    character(len=128)              :: err_message
    integer :: i1, i2, i3, ierr
    real(kind=4), dimension(:), allocatable :: lat4, pres4, time4, tempin
    real(kind=4) :: blatc4

    ! initialize error message
    err_message = ""

    ! Get dimensions from data file
    open(unit=fileID,file=trim(file), form='unformatted', convert='big_endian', iostat=ierr, iomsg=err_message)
    if (ierr /= 0 ) return
    read (fileID, iostat=ierr, iomsg=err_message) this%ncf, this%nlat, this%nlev, this%ntime
    if (ierr /= 0 ) return
    rewind(fileID)
    
    allocate (this%lat(this%nlat))
    allocate (this%pres(this%nlev))
    allocate (this%po3(this%nlev))
    allocate (this%time(this%ntime+1))
    allocate (this%data(this%nlat,this%nlev,this%ncf,this%ntime))
    
    allocate(lat4(this%nlat), pres4(this%nlev), time4(this%ntime+1))
    read (fileID, iostat=ierr, iomsg=err_message) this%ncf, this%nlat, this%nlev, this%ntime, lat4, pres4, time4
    if (ierr /= 0 ) return
    
    ! Store 
    this%pres(:) = pres4(:)
    this%po3(:)  = log(100.0*this%pres(:)) ! from mb to ln(Pa)
    this%lat(:)  = lat4(:)
    this%time(:) = time4(:)
    deallocate(lat4, pres4, time4)
    
    allocate(tempin(this%nlat))
    do i1=1,this%ntime
       do i2=1,this%ncf
          do i3=1,this%nlev
             read(fileID, iostat=ierr, iomsg=err_message) tempin
             if (ierr /= 0 ) return
             this%data(:,i3,i2,i1) = tempin(:)
          enddo
       enddo
    enddo
    deallocate(tempin)
    close(fileID)

  end function load_o3prog

  ! #########################################################################################
  ! Procedure (type-bound) for setting up interpolation indices between data-grid and 
  ! model-grid. 
  ! Called once during initialization
  ! #########################################################################################
  subroutine setup_o3prog(this, lat, idx1, idx2, idxh)
    class(ty_ozphys), intent(in)  :: this
    real(kind_phys),  intent(in)  :: lat(:)
    integer,          intent(out) :: idx1(:), idx2(:)
    real(kind_phys),  intent(out) :: idxh(:)
    integer :: i,j

    do j=1,size(lat)
       idx2(j) = this%nlat + 1
       do i=1,this%nlat
          if (lat(j) < this%lat(i)) then
             idx2(j) = i
             exit
          endif
       enddo
       idx1(j) = max(idx2(j)-1,1)
       idx2(j) = min(idx2(j),this%nlat)
       if (idx2(j) .ne. idx1(j)) then
          idxh(j) = (lat(j) - this%lat(idx1(j))) / (this%lat(idx2(j)) - this%lat(idx1(j)))
       else
          idxh(j) = 1.0
       endif
    enddo

  end subroutine setup_o3prog

  ! #########################################################################################
  ! Procedure (type-bound) for updating data used in prognostic ozone scheme.
  ! #########################################################################################
  subroutine update_o3prog(this, idx1, idx2, idxh, rjday, idxt1, idxt2, ozpl)
    class(ty_ozphys), intent(in)  :: this
    integer,          intent(in)  :: idx1(:), idx2(:)
    real(kind_phys),  intent(in)  :: idxh(:)
    real(kind_phys),  intent(in)  :: rjday
    integer,          intent(in)  :: idxt1, idxt2
    real(kind_phys),  intent(out) :: ozpl(:,:,:)
    integer :: nc, l, j, j1, j2
    real(kind_phys) :: tem, tx1, tx2

    tx1 = (this%time(idxt2) - rjday) / (this%time(idxt2) - this%time(idxt1))
    tx2 = 1.0 - tx1
 
    do nc=1,this%ncf
       do l=1,this%nlev
          do j=1,size(ozpl(:,1,1))
             j1  = idx1(j)
             j2  = idx2(j)
             tem = 1.0 - idxh(j)
             ozpl(j,l,nc) = tx1*(tem*this%data(j1,l,nc,idxt1)+idxh(j)*this%data(j2,l,nc,idxt1)) &
                  + tx2*(tem*this%data(j1,l,nc,idxt2)+idxh(j)*this%data(j2,l,nc,idxt2))
          enddo
       enddo
    enddo

  end subroutine update_o3prog

  ! #########################################################################################
  ! Procedure (type-bound) for NRL prognostic ozone (2015).
  ! #########################################################################################
  subroutine run_o3prog_2015(this, con_1ovg, dt, p, t, dp, ozpl, oz, do_diag, do3_dt_prd, &
       do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz)
    class(ty_ozphys), intent(in) :: this
    real(kind_phys),  intent(in) :: &
         con_1ovg       ! Physical constant: One divided by gravitational acceleration (m-1 s2)
    real(kind_phys),  intent(in) :: &
         dt             ! Model timestep (sec)
    real(kind_phys),  intent(in), dimension(:,:) :: &
         p,           & ! Model Pressure (Pa)
         t,           & ! Model temperature (K)
         dp             ! Model layer thickness (Pa)
    real(kind_phys), intent(in), dimension(:,:,:) :: &
         ozpl           ! Ozone forcing data
    real(kind_phys), intent(inout), dimension(:,:) :: &
         oz             ! Ozone concentration updated by physics
    logical, intent(in) :: do_diag
    real(kind_phys), intent(inout), dimension(:,:), optional :: &
         do3_dt_prd,  & ! Physics tendency: production and loss effect
         do3_dt_ozmx, & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp, & ! Physics tendency: temperature effect
         do3_dt_ohoz    ! Physics tendency: overhead ozone effect

    integer :: k, kmax, kmin, iLev, iCol, nCol, nLev, iCf
    logical, dimension(size(p,1)) :: flg
    real(kind_phys) :: pmax, pmin, tem, temp
    real(kind_phys), dimension(size(p,1)) :: wk1, wk2, wk3, ozib
    real(kind_phys), dimension(size(p,1),this%ncf) :: prod
    real(kind_phys), dimension(size(p,1),size(p,2)) :: ozi
    real(kind_phys), dimension(size(p,1),size(p,2)+1) :: colo3, coloz

    ! Dimensions
    nCol = size(p,1)
    nLev = size(p,2)

    ! Temporaries
    ozi = oz

    colo3(:,nLev+1) = 0.0
    coloz(:,nLev+1) = 0.0

    do iLev=nLev,1,-1
       pmin =  1.0e10
       pmax = -1.0e10

       do iCol=1,nCol
          wk1(iCol)    = log(p(iCol,iLev))
          pmin         = min(wk1(iCol), pmin)
          pmax         = max(wk1(iCol), pmax)
          prod(iCol,:) = 0._kind_phys
       enddo
       kmax = 1
       kmin = 1
       do k=1,this%nlev-1
          if (pmin < this%po3(k)) kmax = k
          if (pmax < this%po3(k)) kmin = k
       enddo
       !
       do k=kmin,kmax
          temp = 1.0 / (this%po3(k) - this%po3(k+1))
          do iCol=1,nCol
             flg(iCol) = .false.
             if (wk1(iCol) < this%po3(k) .and. wk1(iCol) >= this%po3(k+1)) then
                flg(iCol) = .true.
                wk2(iCol) = (wk1(iCol) - this%po3(k+1)) * temp
                wk3(iCol) = 1.0 - wk2(iCol)
             endif
          enddo
          do iCf=1,this%ncf
             do iCol=1,nCol
                if (flg(iCol)) then
                   prod(iCol,iCf)  = wk2(iCol) * ozpl(iCol,k,iCf) + wk3(iCol) * ozpl(iCol,k+1,iCf)
                endif
             enddo
          enddo
       enddo

       do iCf=1,this%ncf
          do iCol=1,nCol
             if (wk1(iCol) < this%po3(this%nlev)) then
                prod(iCol,iCf) = ozpl(iCol,this%nlev,iCf)
             endif
             if (wk1(iCol) >= this%po3(1)) then
                prod(iCol,iCf) = ozpl(iCol,1,iCf)
             endif
          enddo
       enddo
       do iCol=1,nCol
          colo3(iCol,iLev) = colo3(iCol,iLev+1) + ozi(iCol,iLev)  * dp(iCol,iLev)*con_1ovg
          coloz(iCol,iLev) = coloz(iCol,iLev+1) + prod(iCol,6) * dp(iCol,iLev)*con_1ovg
          prod(iCol,2)     = min(prod(iCol,2), 0.0)
       enddo
       do iCol=1,nCol
          ozib(iCol) = ozi(iCol,iLev)
          tem        = prod(iCol,1) - prod(iCol,2) * prod(iCol,6) &
                                    + prod(iCol,3) * (t(iCol,iLev) - prod(iCol,5)) &
                                    + prod(iCol,4) * (colo3(iCol,iLev)-coloz(iCol,iLev))
          oz(iCol,iLev) = (ozib(iCol)  + tem*dt) / (1.0 - prod(iCol,2)*dt)
       enddo

       ! Diagnostics (optional)
       if (do_diag) then
          do3_dt_prd(:,iLev)  = prod(:,1) * dt
          do3_dt_ozmx(:,iLev) = prod(:,2) * (oz(:,iLev) - prod(:,6)) * dt
          do3_dt_temp(:,iLev) = prod(:,3)*(t(:,iLev)-prod(:,5))*dt
          do3_dt_ohoz(:,iLev) = prod(:,4) * (colo3(:,iLev)-coloz(:,iLev))*dt
       endif
    enddo

    return
  end subroutine run_o3prog_2015

  ! #########################################################################################
  ! Procedure (type-bound) for NRL prognostic ozone (2006).
  ! #########################################################################################
  subroutine run_o3prog_2006(this, con_1ovg, dt, p, t, dp, ozpl, oz, do_diag, do3_dt_prd, &
       do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz)
    class(ty_ozphys), intent(in) :: this
    real(kind_phys),  intent(in) :: &
         con_1ovg       ! Physical constant: One divided by gravitational acceleration (m-1 s2)
    real(kind_phys),  intent(in) :: &
         dt             ! Model timestep (sec)
    real(kind_phys),  intent(in), dimension(:,:) :: &
         p,           & ! Model Pressure (Pa)
         t,           & ! Model temperature (K)
         dp             ! Model layer thickness (Pa)
    real(kind_phys), intent(in), dimension(:,:,:) :: &
         ozpl           ! Ozone forcing data
    real(kind_phys), intent(inout), dimension(:,:) :: &
         oz             ! Ozone concentration updated by physics
    logical, intent(in) :: do_diag
    real(kind_phys), intent(inout), dimension(:,:) :: &
         do3_dt_prd,  & ! Physics tendency: production and loss effect
         do3_dt_ozmx, & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp, & ! Physics tendency: temperature effect
         do3_dt_ohoz    ! Physics tendency: overhead ozone effect

    ! Locals
    integer :: k, kmax, kmin, iLev, iCol, nCol, nLev, iCf
    logical, dimension(size(p,1)) :: flg
    real(kind_phys) :: pmax, pmin, tem, temp
    real(kind_phys), dimension(size(p,1)) :: wk1, wk2, wk3, ozib
    real(kind_phys), dimension(size(p,1),this%ncf) :: prod
    real(kind_phys), dimension(size(p,1),size(p,2)) :: ozi
    real(kind_phys), dimension(size(p,1),size(p,2)+1) :: colo3, coloz

    ! Dimensions
    nCol = size(p,1)
    nLev = size(p,2)

    ! Temporaries
    ozi = oz

    !> - Calculate vertical integrated column ozone values.
    if (this%ncf > 2) then
       colo3(:,nLev+1) = 0.0
       do iLev=nLev,1,-1
          do iCol=1,nCol
             colo3(iCol,iLev) = colo3(iCol,iLev+1) + ozi(iCol,iLev) * dp(iCol,iLev) * con_1ovg
          enddo
       enddo
    endif

    !> - Apply vertically linear interpolation to the ozone coefficients.
    do iLev=1,nLev
       pmin =  1.0e10
       pmax = -1.0e10

       do iCol=1,nCol
          wk1(iCol)    = log(p(iCol,iLev))
          pmin         = min(wk1(iCol), pmin)
          pmax         = max(wk1(iCol), pmax)
          prod(iCol,:) = 0._kind_phys
       enddo
       kmax = 1
       kmin = 1
       do k=1,this%nlev-1
          if (pmin < this%po3(k)) kmax = k
          if (pmax < this%po3(k)) kmin = k
       enddo

       do k=kmin,kmax
          temp = 1.0 / (this%po3(k) - this%po3(k+1))
          do iCol=1,nCol
             flg(iCol) = .false.
             if (wk1(iCol) < this%po3(k) .and. wk1(iCol) >= this%po3(k+1)) then
                flg(iCol) = .true.
                wk2(iCol) = (wk1(iCol) - this%po3(k+1)) * temp
                wk3(iCol) = 1.0 - wk2(iCol)
             endif
          enddo
          do iCf=1,this%ncf
             do iCol=1,nCol
                if (flg(iCol)) then
                   prod(iCol,iCf)  = wk2(iCol) * ozpl(iCol,k,iCf) + wk3(iCol) * ozpl(iCol,k+1,iCf)
                endif
             enddo
          enddo
       enddo

       do iCf=1,this%ncf
          do iCol=1,nCol
             if (wk1(iCol) < this%po3(this%nlev)) then
                prod(iCol,iCf) = ozpl(iCol,this%nlev,iCf)
             endif
             if (wk1(iCol) >= this%po3(1)) then
                prod(iCol,iCf) = ozpl(iCol,1,iCf)
             endif
          enddo
       enddo

       if (this%ncf == 2) then
          do iCol=1,nCol
             ozib(iCol)    = ozi(iCol,iLev)
             oz(iCol,iLev) = (ozib(iCol) + prod(iCol,1)*dt) / (1.0 + prod(iCol,2)*dt)
          enddo
       endif

       if (this%ncf == 4) then
          do iCol=1,nCol
             ozib(iCol)    = ozi(iCol,iLev)
             tem           = prod(iCol,1) + prod(iCol,3)*t(iCol,iLev) + prod(iCol,4)*colo3(iCol,iLev+1)
             oz(iCol,iLev) = (ozib(iCol)  + tem*dt) / (1.0 + prod(iCol,2)*dt)
          enddo
       endif

       ! Diagnostics (optional)
       if (do_diag) then
          do3_dt_prd(:,iLev)  = prod(:,1)*dt
          do3_dt_ozmx(:,iLev) = (oz(:,iLev) - ozib(:))
          do3_dt_temp(:,iLev) = prod(:,3) * t(:,iLev) * dt
          do3_dt_ohoz(:,iLev) = prod(:,4) * colo3(:,iLev) * dt
       endif
    enddo

    return
  end subroutine run_o3prog_2006

  ! #########################################################################################
  ! Procedure (type-bound) for NRL updating climotological ozone.
  ! #########################################################################################
  subroutine run_o3clim(this, lat, prslk, con_pi, oz)
    class(ty_ozphys), intent(in) :: this
    real(kind_phys),  intent(in) :: &
         con_pi  ! Physics constant: Pi
    real(kind_phys),  intent(in), dimension(:)   :: &
         lat     ! Grid latitude
    real(kind_phys),  intent(in), dimension(:,:) :: &
         prslk   ! Exner function
    real(kind_phys),  intent(out), dimension(:,:) :: &
         oz      ! Ozone concentration updated by climotology

    integer :: nCol, iCol, nLev, iLev, j, j1, j2, l, ll
    real(kind_phys) :: elte, deglat, tem, tem1, tem2, tem3, tem4, temp
    real(kind_phys), allocatable :: o3i(:,:),wk1(:)
    logical :: top_at_1

    nCol = size(prslk(:,1))
    nLev = size(prslk(1,:))
    allocate(o3i(nCol, this%nlevc),wk1(nCol))

    ! What is vertical ordering?
    top_at_1 = (prslk(1,1) .lt.  prslk(1, nLev))

    elte = this%blatc + (this%nlatc-1)*this%dphiozc

    do iCol = 1, nCol
       deglat = lat(iCol) * 180.0 / con_pi
       if (deglat > this%blatc .and. deglat < elte) then
          tem1 = (deglat - this%blatc) / this%dphiozc + 1
          j1   = tem1
          j2   = j1 + 1
          tem1 = tem1 - j1
       elseif (deglat <= this%blatc) then
          j1   = 1
          j2   = 1
          tem1 = 1.0
       elseif (deglat >= elte) then
          j1   = this%nlatc
          j2   = this%nlatc
          tem1 = 1.0
       endif
       
       tem2 = 1.0 - tem1
       do j = 1, this%nlevc
          tem3        = tem2*this%datac(j1,j,this%k1oz) + tem1*this%datac(j2,j,this%k1oz)
          tem4        = tem2*this%datac(j1,j,this%k2oz) + tem1*this%datac(j2,j,this%k2oz)
          o3i(iCol,j) = tem4*this%facoz               + tem3*(1.0 - this%facoz)
       enddo
    enddo

    do iLev = 1, nLev
       ll = iLev
       if (.not. top_at_1) ll = nLev - iLev + 1

       do iCol = 1, nCol
          wk1(iCol) = prslk(iCol,ll)
       enddo

       do j = 1, this%nlevc-1
          temp = 1.0 / (this%pkstr(j+1) - this%pkstr(j))
          do iCol = 1, nCol
             if (wk1(iCol) > this%pkstr(j) .and. wk1(iCol) <= this%pkstr(j+1)) then
                tem       = (this%pkstr(j+1) - wk1(iCol)) * temp
                oz(iCol,ll) = tem * o3i(iCol,j) + (1.0 - tem) * o3i(iCol,j+1)
             endif
          enddo
       enddo
       
       do iCol = 1, nCol
          if (wk1(iCol) > this%pkstr(this%nlevc)) oz(iCol,ll) = o3i(iCol,this%nlevc)
          if (wk1(iCol) < this%pkstr(1))          oz(iCol,ll) = o3i(iCol,1)
       enddo
    enddo

    return
  end subroutine run_o3clim

  ! #########################################################################################
  ! Procedure (type-bound) for loading data for climotological ozone.
  ! #########################################################################################
  function load_o3clim(this, file, fileID) result (err_message)
    class(ty_ozphys), intent(inout) :: this
    integer,          intent(in)    :: fileID
    character(len=*), intent(in)    :: file
    character(len=128)              :: err_message

    ! Locals
    real(kind=4) :: blatc4
    integer :: iLev, iLat, imo, ierr
    real(kind=4), allocatable :: o3clim4(:,:,:), pstr4(:)
    integer, allocatable      :: imond(:), ilatt(:,:)

    ! initialize error message
    err_message = ""

    open(unit=fileID,file=trim(file),form='unformatted',convert='big_endian', iostat=ierr, iomsg=err_message)
    if (ierr /= 0 ) return
    read (fileID,end=101,iostat=ierr,iomsg=err_message) this%nlatc, this%nlevc, this%ntimec, blatc4
    if (ierr /= 0 ) return

101 if (this%nlevc  < 10 .or. this%nlevc > 100) then
       rewind (fileID)
       this%nlevc = 17
       this%nlatc = 18
       this%blatc = -85.0
    else
       this%blatc = blatc4
    endif
    this%nlat    = 2
    this%nlev    = 1
    this%ntimec  = 1
    this%ncf     = 0
    this%dphiozc = -(this%blatc+this%blatc)/(this%nlatc-1)

    allocate (o3clim4(this%nlatc,this%nlevc,12), pstr4(this%nlevc), imond(12), ilatt(this%nlatc,12) )

    allocate (this%pkstr(this%nlevc), this%pstr(this%nlevc), this%datac(this%nlatc,this%nlevc,12))
    if ( this%nlevc == 17 ) then ! For the operational ozone climatology
       do iLev = 1, this%nlevc
          read (fileID,15,iostat=ierr,iomsg=err_message) pstr4(iLev)
          if (ierr /= 0 ) return
15        format(f10.3)
       enddo

       do imo = 1, 12
          do iLat = 1, this%nlatc
             read (fileID,16,iostat=ierr,iomsg=err_message) imond(imo), ilatt(iLat,imo), (o3clim4(iLat,iLev,imo),iLev=1,10)
             if (ierr /= 0 ) return
16           format(i2,i4,10f6.2)
             read (fileID,20,iostat=ierr,iomsg=err_message) (o3clim4(iLat,iLev,imo),iLev=11,this%nlevc)
             if (ierr /= 0 ) return
20           format(6x,10f6.2)
          enddo
       enddo
    else ! For newer ozone climatology
       read (fileID)
       do iLev = 1, this%nlevc
          read (fileID) pstr4(iLev)
       enddo
       
       do imo = 1, 12
          do iLev = 1, this%nlevc
              read (fileID,iostat=ierr,iomsg=err_message) (o3clim4(iLat,iLev,imo),iLat=1,this%nlatc)
              if (ierr /= 0 ) return
           enddo
        enddo
     endif   ! end if_this%nlevc_block

     do imo = 1, 12
        do iLev = 1, this%nlevc
           do iLat = 1, this%nlatc
              this%datac(iLat,iLev,imo) = o3clim4(iLat,iLev,imo) * 1.655e-6
           enddo
        enddo
     enddo
     
     do iLev = 1, this%nlevc
        this%pstr(iLev)  = pstr4(iLev)
        this%pkstr(iLev) = fpkapx(this%pstr(iLev)*100.0)
     enddo
     
   end function load_o3clim

   ! #########################################################################################
   ! Procedure (type-bound) for updating temporal interpolation index when using climotological
   ! ozone
   ! #########################################################################################
   subroutine update_o3clim(this, imon, iday, ihour, loz1st)
     class(ty_ozphys), intent(inout) :: this
     integer, intent(in) :: imon, iday, ihour
     logical, intent(in) :: loz1st

     integer ::  midmon=15, midm=15, midp=45, id
     integer, parameter, dimension(13) :: mdays = (/31,28,31,30,31,30,31,31,30,31,30,31,30/)
     logical :: change

     midmon = mdays(imon)/2 + 1
     change = loz1st .or. ( (iday==midmon) .and. (ihour==0) )
    
     if ( change ) then
        if ( iday < midmon ) then
           this%k1oz = mod(imon+10, 12) + 1
           midm = mdays(this%k1oz)/2 + 1
           this%k2oz = imon
           midp = mdays(this%k1oz) + midmon
        else
           this%k1oz = imon
           midm = midmon
           this%k2oz = mod(imon, 12) + 1
           midp = mdays(this%k2oz)/2 + 1 + mdays(this%k1oz)
        endif
     endif
    
     if (iday < midmon) then
        id = iday + mdays(this%k1oz)
     else
        id = iday
     endif
    
     this%facoz = float(id - midm) / float(midp - midm)

   end subroutine update_o3clim

 end module module_ozphys

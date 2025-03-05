! #########################################################################################
!> \section arg_table_module_h2ophys Argument table                               
!! \htmlinclude module_h2ophys.html                                               
!!
! #########################################################################################
module module_h2ophys
  use machine,  only : kind_phys
  implicit none

  public ty_h2ophys

! ######################################################################################### 
!> \section arg_table_ty_h2ophys Argument Table 
!! \htmlinclude ty_h2ophys.html
!!
!> Derived type containing data and procedures needed by h2o photochemistry parameterization
!! *Note* All data field are ordered from surface-to-toa.
!!
! #########################################################################################
  type ty_h2ophys
     integer                      :: nlat          !< Number of latitudes.
     integer                      :: nlev          !< Number of vertical layers.
     integer                      :: ntime         !< Number of times.
     integer                      :: ncf           !< Number of coefficients.
     real(kind_phys), allocatable :: lat(:)        !< Latitude.
     real(kind_phys), allocatable :: pres(:)       !< Pressure levels.
     real(kind_phys), allocatable :: ph2o(:)       !< Natural log pressure of levels.
     real(kind_phys), allocatable :: time(:)       !< Time.
     real(kind_phys), allocatable :: data(:,:,:,:) !< H20 forcing data (raw)
     contains
       procedure, public :: load
       procedure, public :: setup
       procedure, public :: update
       procedure, public :: run
  end type ty_h2ophys
  
contains
  ! #########################################################################################
  ! Procedure (type-bound) for loading data.
  ! #########################################################################################
  function load(this, file, fileID) result (err_message)
    class(ty_h2ophys), intent(inout) :: this
    integer,           intent(in)    :: fileID
    character(len=*),  intent(in)    :: file
    character(len=128)               :: err_message
    integer :: i1, i2, i3, ierr
    real(kind=4), dimension(:), allocatable :: lat4, pres4, time4, tempin

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
    allocate (this%ph2o(this%nlev))
    allocate (this%time(this%ntime+1))
    allocate (this%data(this%nlat,this%nlev,this%ncf,this%ntime))

    allocate(lat4(this%nlat), pres4(this%nlev), time4(this%ntime+1))
    read (fileID, iostat=ierr, iomsg=err_message) this%ncf, this%nlat, this%nlev, this%ntime, lat4, pres4, time4
    if (ierr /= 0 ) return

    ! Store
    this%pres(:) = pres4(:)
    this%ph2o(:) = log(100.0*this%pres(:)) ! from mb to ln(Pa)
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

  end function load

  ! #########################################################################################
  ! Procedure (type-bound) for setting up interpolation indices between data-grid and 
  ! model-grid. 
  ! #########################################################################################
  subroutine setup(this, lat, idx1, idx2, idxh)
    class(ty_h2ophys), intent(in)  :: this
    real(kind_phys),   intent(in)  :: lat(:)
    integer,           intent(out) :: idx1(:), idx2(:)
    real(kind_phys),   intent(out) :: idxh(:)
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

  end subroutine setup

  ! #########################################################################################
  ! Procedure (type-bound) for updating data.
  ! #########################################################################################
  subroutine update(this, idx1, idx2, idxh, rjday, idxt1, idxt2, h2opl)
    class(ty_h2ophys), intent(in)  :: this
    integer,           intent(in)  :: idx1(:), idx2(:)
    real(kind_phys),   intent(in)  :: idxh(:)
    real(kind_phys),   intent(in)  :: rjday
    integer,           intent(in)  :: idxt1, idxt2
    real(kind_phys),   intent(out) :: h2opl(:,:,:)
    integer :: nc, l, j, j1, j2
    real(kind_phys) :: tem, tx1, tx2

    tx1 = (this%time(idxt2) - rjday) / (this%time(idxt2) - this%time(idxt1))
    tx2 = 1.0 - tx1

    do nc=1,this%ncf
       do l=1,this%nlev
          do j=1,size(h2opl(:,1,1))
             j1  = idx1(j)
             j2  = idx2(j)
             tem = 1.0 - idxh(j)
             h2opl(j,l,nc) = tx1*(tem*this%data(j1,l,nc,idxt1)+idxh(j)*this%data(j2,l,nc,idxt1)) &
                  + tx2*(tem*this%data(j1,l,nc,idxt2)+idxh(j)*this%data(j2,l,nc,idxt2))
          enddo
       enddo
    enddo
    
  end subroutine update

  ! #########################################################################################
  ! Procedure (type-bound) for NRL stratospheric h2o photochemistry physics.
  ! #########################################################################################
  subroutine run(this, dt, p, h2opltc, h2o, dqv_dt_prd, dqv_dt_qv)
    class(ty_h2ophys), intent(in) :: this
    real(kind_phys),   intent(in) :: &
         dt        ! Model timestep (sec)
    real(kind_phys),   intent(in), dimension(:,:) :: &
         p         ! Model Pressure (Pa)
    real(kind_phys), intent(in), dimension(:,:,:) :: &
         h2opltc   ! h2o forcing data
    real(kind_phys), intent(inout), dimension(:,:) :: &
         h2o       ! h2o concentration (updated)
    real(kind_phys), intent(inout), dimension(:, :), optional :: &
         dqv_dt_prd, & ! Net production/loss effect
         dqv_dt_qv     ! water vapor effect
    
    integer :: nCol, nLev, iCol, iLev, iCf, kmax, kmin, k
    logical, dimension(size(p,1)) :: flg
    real(kind_phys) :: pmax, pmin, temp
    real(kind_phys), dimension(size(p,1)) :: wk1, wk2, wk3, h2oib
    real(kind_phys), dimension(size(p,1),this%ncf) :: pltc
    real(kind_phys), parameter :: prsmax=10000.0, pmaxl=log(prsmax)
    
    ! Dimensions
    nCol = size(p,1)
    nLev = size(p,2)
    
    do iLev=1,nLev
       pmin =  1.0e10
       pmax = -1.0e10
       do iCol=1,nCol
          wk1(iCol)    = log(p(iCol,iLev))
          pmin         = min(wk1(iCol), pmin)
          pmax         = max(wk1(iCol), pmax)
          pltc(iCol,:) = 0._kind_phys
       enddo
       if (pmin < pmaxl) then
          kmax = 1
          kmin = 1
          do k=1,this%nlev-1
             if (pmin < this%ph2o(k)) kmax = k
             if (pmax < this%ph2o(k)) kmin = k
          enddo

          do k=kmin,kmax
             temp = 1.0 / (this%ph2o(k) - this%ph2o(k+1))
             do iCol=1,nCol
                flg(iCol) = .false.
                if (wk1(iCol) < this%ph2o(k) .and. wk1(iCol) >= this%ph2o(k+1)) then
                   flg(iCol) = .true.
                   wk2(iCol) = (wk1(iCol) - this%ph2o(k+1)) * temp
                   wk3(iCol) = 1.0 - wk2(iCol)
                endif
             enddo
             do iCf=1,this%ncf
                do iCol=1,nCol
                   if (flg(iCol)) then
                      pltc(iCol,iCf)  = wk2(iCol) * h2opltc(iCol,k,iCf) + wk3(iCol) * h2opltc(iCol,k+1,iCf)
                   endif
                enddo
             enddo
          enddo
          
          do iCf=1,this%ncf
             do iCol=1,nCol
                if (wk1(iCol) < this%ph2o(this%nlev)) then
                   pltc(iCol,iCf) = h2opltc(iCol,this%nlev,iCf)
                endif
                if (wk1(iCol) >= this%ph2o(1)) then
                   pltc(iCol,iCf) = h2opltc(iCol,1,iCf)
                endif
             enddo
          enddo
       endif

       do iCol=1,nCol
          if (p(iCol,iLev) < prsmax .and. pltc(iCol,2) /= 0.0) then
             h2oib(iCol)  = h2o(iCol,iLev)
             temp      = 1.0 / pltc(iCol,2)
             h2o(iCol,iLev)  = (h2oib(iCol) + (pltc(iCol,1)+pltc(iCol,3)*temp)*dt) / (1.0 + temp*dt)
          endif
       enddo

       if (present(dqv_dt_prd)) dqv_dt_prd(:, iLev) = pltc(:, 1) * dt
       if (present(dqv_dt_qv))  dqv_dt_qv(:, iLev)  = (h2o(:, iLev) - pltc(:, 3)) * dt / pltc(:, 2)
    enddo


    return
  end subroutine run

end module module_h2ophys

! ###########################################################################################
!> \file ozphys_2015.F90
!!
! ###########################################################################################
module ozphys_2015
  use machine,       only: kind_phys, kind_dbl_prec, kind_sngl_prec
  use module_ozphys, only: ty_ozphys
  implicit none
  public ozphys_2015_run
contains

! ###########################################################################################
!>\defgroup GFS_ozphys_2015 GFS Ozone Photochemistry (2015) Module
!! This module contains the CCPP-compliant Ozone 2015 photochemistry scheme.
!> @{
!> The operational GFS currently parameterizes ozone production and destruction based on 
!! monthly mean coefficients ( \c ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77) provided by 
!! Naval Research Laboratory through CHEM2D chemistry model
!! (McCormack et al. (2006) \cite mccormack_et_al_2006).
!! (https://doi.org/10.5194/acp-6-4943-2006)
!!
!> \section genal_ozphys_2015 GFS ozphys_2015_run General Algorithm
!> -  This code assumes that both 2D fields are ordered from bottom to top.
!> -  This code is specifically for NRL parameterization and climatological T and O3 are in 
!     location 5 and 6 of ozpl array
!!\author   June 2015 - Shrinivas Moorthi
!!\modified May  2023 - Dustin Swales
! ###########################################################################################

! ###########################################################################################
! SUBROUTINE ozphys_2015_run
! ###########################################################################################
!! \section arg_table_ozphys_2015_run Argument Table
!! \htmlinclude ozphys_2015_run.html
!!
  subroutine ozphys_2015_run (oz_phys, ozphys, nCol, nLev, dt, oz, tin, prsl, ozpl,     &
       delp, con_1ovg, do3_dt_prd, do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         oz_phys        ! Flag for ozone_physics_2015 scheme.
    type(ty_ozphys),intent(in) :: &
         ozphys
    real(kind_phys),intent(in) :: &
         con_1ovg       ! Physical constant: One divided by gravitational acceleration (m-1 s2)
    integer, intent(in) :: &
         nCol,        & ! Horizontal dimension
         nLev           ! Number of vertical layers
    real(kind_phys), intent(in) :: &
         dt             ! Physics timestep (seconds)
    real(kind_phys), intent(in), dimension(:,:) :: &
         prsl,        & ! Air-pressure (Pa)
         tin,         & ! Temperature of new-state (K)
         delp           ! Difference between mid-layer pressures (Pa)
    real(kind_phys), intent(in), dimension(:,:,:) :: &
         ozpl           ! Ozone forcing data

    ! Outputs (optional)
    real(kind=kind_phys), intent(inout), dimension(:,:), pointer, optional :: &
         do3_dt_prd,  & ! Physics tendency: production and loss effect
         do3_dt_ozmx, & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp, & ! Physics tendency: temperature effect
         do3_dt_ohoz    ! Physics tendency: overhead ozone effect

    ! Outputs
    real(kind=kind_phys), intent(inout), dimension(:,:) :: &
         oz             ! Ozone concentration updated by physics
    character(len=*), intent(out) :: &
         errmsg         ! CCPP error message
    integer, intent(out) :: &
         errflg         ! CCPP error flag

    ! Locals
    integer :: k, kmax, kmin, l, i, j
    logical, dimension(nCol) :: flg
    real(kind_phys) :: pmax, pmin, tem, temp
    real(kind_phys), dimension(nCol) :: wk1, wk2, wk3, ozib
    real(kind_phys), dimension(nCol,ozphys%ncf) :: prod
    real(kind_phys), dimension(nCol,nLev) :: ozi
    real(kind_phys), dimension(nCol,nLev+1) :: colo3, coloz

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Sanity checks
    if (.not.oz_phys) then
       write (errmsg,'(*(a))') 'Logic error: oz_phys_2015 == .false.'
       errflg = 1
       return
    endif

    ! Temporaries
    ozi = oz

    colo3(:,nLev+1) = 0.0
    coloz(:,nLev+1) = 0.0

    do l=nLev,1,-1
       pmin =  1.0e10
       pmax = -1.0e10

       do i=1,nCol
          wk1(i) = log(prsl(i,l))
          pmin   = min(wk1(i), pmin)
          pmax   = max(wk1(i), pmax)
          prod(i,:) = 0.0
       enddo
       kmax = 1
       kmin = 1
       do k=1,ozphys%nlev-1
          if (pmin < ozphys%po3(k)) kmax = k
          if (pmax < ozphys%po3(k)) kmin = k
       enddo
       !
       do k=kmin,kmax
          temp = 1.0 / (ozphys%po3(k) - ozphys%po3(k+1))
          do i=1,nCol
             flg(i) = .false.
             if (wk1(i) < ozphys%po3(k) .and. wk1(i) >= ozphys%po3(k+1)) then
                flg(i) = .true.
                wk2(i) = (wk1(i) - ozphys%po3(k+1)) * temp
                wk3(i) = 1.0 - wk2(i)
             endif
          enddo
          do j=1,ozphys%ncf
             do i=1,nCol
                if (flg(i)) then
                   prod(i,j)  = wk2(i) * ozpl(i,k,j) + wk3(i) * ozpl(i,k+1,j)
                endif
             enddo
          enddo
       enddo

       do j=1,ozphys%ncf
          do i=1,nCol
             if (wk1(i) < ozphys%po3(ozphys%nlev)) then
                prod(i,j) = ozpl(i,ozphys%nlev,j)
             endif
             if (wk1(i) >= ozphys%po3(1)) then
                prod(i,j) = ozpl(i,1,j)
             endif
          enddo
       enddo
       do i=1,nCol
          colo3(i,l) = colo3(i,l+1) + ozi(i,l)  * delp(i,l)*con_1ovg
          coloz(i,l) = coloz(i,l+1) + prod(i,6) * delp(i,l)*con_1ovg
          prod(i,2)  = min(prod(i,2), 0.0)
       enddo
       do i=1,nCol
          ozib(i)  = ozi(i,l)            ! no filling
          tem      = prod(i,1) - prod(i,2) * prod(i,6) + prod(i,3) * (tin(i,l) - prod(i,5)) &
                                                       + prod(i,4) * (colo3(i,l)-coloz(i,l))
          oz(i,l) = (ozib(i)  + tem*dt) / (1.0 - prod(i,2)*dt)
       enddo

       ! Diagnostics (optional)
       if (associated(do3_dt_prd))  do3_dt_prd(:,l)  = (prod(:,1)-prod(:,2)*prod(:,6))*dt
       if (associated(do3_dt_ozmx)) do3_dt_ozmx(:,l) = (oz(:,l) - ozib(:))
       if (associated(do3_dt_temp)) do3_dt_temp(:,l) = prod(:,3)*(tin(:,l)-prod(:,5))*dt
       if (associated(do3_dt_ohoz)) do3_dt_ohoz(:,l) = prod(:,4) * (colo3(:,l)-coloz(:,l))*dt
    enddo

    return
  end subroutine ozphys_2015_run
!> @}
end module ozphys_2015

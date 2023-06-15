! ###########################################################################################
!> \file ozphys_2015.F90
!!
! ###########################################################################################
module ozphys_2015
  use machine, only : kind_phys, kind_dbl_prec, kind_sngl_prec
  implicit none
  public ozphys_2015_init, ozphys_2015_timestep_init, ozphys_2015_run
contains

! ###########################################################################################
!>\defgroup GFS_ozphys_2015 GFS Ozone Photochemistry (2015) Module
!! This module contains the CCPP-compliant Ozone 2015 photochemistry scheme.
!> @{
!> \section arg_table_ozphys_2015_init Argument Table
!! \htmlinclude ozphys_2015_init.html
!!
! ###########################################################################################
  subroutine ozphys_2015_init(oz_phys_2015, nPts, latsozp, oz_lat, dlat, jindx1, jindx2,    &
       ddy, errmsg, errflg)
    ! Inputs
    logical, intent(in) :: &
         oz_phys_2015 ! Control flag for NRL 2015 ozone physics scheme
    integer, intent(in) :: &
         nPts,      & ! Horizontal dimension
         latsozp      ! Number of latitudes in ozone data
    real(kind_phys),  intent(in), dimension(:) :: &
         oz_lat,    & ! Latitudes of ozone data
         dlat         ! Latitudes of grid
    ! Outputs
    integer, intent(out), dimension(:) :: &
         jindx1,    & ! Interpolation index (low) for ozone data
         jindx2       ! Interpolation index (high) for ozone data
    real(kind_phys), intent(out), dimension(:) :: &
         ddy          ! Interpolation high index for ozone data
    character(len=*), intent(out) :: &
         errmsg       ! CCPP error message
    integer, intent(out) :: &
         errflg       ! CCPP error flag

    ! Local
    integer i,j

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    ! Sanity check
    if (.not.oz_phys_2015) then
       write (errmsg,'(*(a))') 'Logic error: oz_phys_2015 == .false.'
       errflg = 1
       return
    endif

    ! Set indices
    do j=1,nPts
       jindx2(j) = latsozp + 1
       do i=1,latsozp
          if (dlat(j) < oz_lat(i)) then
             jindx2(j) = i
             exit
          endif
       enddo
       jindx1(j) = max(jindx2(j)-1,1)
       jindx2(j) = min(jindx2(j),latsozp)
       if (jindx2(j) .ne. jindx1(j)) then
          ddy(j) = (dlat(j) - oz_lat(jindx1(j))) / (oz_lat(jindx2(j)) - oz_lat(jindx1(j)))
       else
          ddy(j) = 1.0
       endif
    enddo
    
  end subroutine ozphys_2015_init

! ###########################################################################################
!> \section arg_table_ozphys_2015_timestep_init Argument Table
!! \htmlinclude ozphys_2015_timestep_init.html
!!
! ###########################################################################################
  subroutine ozphys_2015_timestep_init(nPts, idate, fhour, jindx1, jindx2, latsozp, levozp, &
       oz_coeff, timeoz, ozplin, oz_time, oz_lat, ddy, ozplout, errmsg, errflg)
    ! Inputs
    integer, intent(in) :: &
         nPts,     & ! Horizontal dimension
         latsozp,  & ! Number of latitudes in ozone data
         levozp,   & ! Number of vertical layers in ozone data
         oz_coeff, & ! Number of coefficients in ozone data
         timeoz      ! Number of times in ozone data
    integer, intent(in),dimension(:) :: &
         idate,    & ! Initial date with different size and ordering
         jindx1,   & ! Interpolation index (low) for ozone
         jindx2      ! Interpolation index (high) for ozone
    real(kind_phys), intent(in) :: &
         fhour       ! Forecast hour
    real(kind_phys), intent(in), dimension(:) :: &
         ddy,      & ! Interpolation high index for ozone data
         oz_lat,   & ! Latitudes for ozone data
         oz_time     ! Time for ozone data
    real(kind_phys), intent(in), dimension(:,:,:,:) :: &
         ozplin      ! Ozone data

    ! Outputs
    real(kind_phys), intent(out), dimension(:,:,:) :: &
         ozplout     ! Ozone forcing data
    character(len=*), intent(out) :: &
         errmsg      ! CCPP error message
    integer,          intent(out) :: &
         errflg      ! CCPP error flag

    ! Local
    integer :: idat(8),jdat(8),iday,j,j1,j2,l,nc,n1,n2,jdow,jdoy,&
         jday,w3kindreal,w3kindint
    real(kind_phys) :: tem, tx1, tx2, rjday
    real(8) :: rinc(5)
    real(4) :: rinc4(5)
    !real(kind_dbl_prec) :: rinc(5)
    !real(kind_sngl_prec) :: rinc4(5)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    !
    idat=0
    idat(1)=idate(4)
    idat(2)=idate(2)
    idat(3)=idate(3)
    idat(5)=idate(1)
    rinc=0.
    rinc(2)=fhour
    call w3kind(w3kindreal,w3kindint)
    if(w3kindreal==4) then
       rinc4=rinc
       CALL w3movdat(rinc4,idat,jdat)
    else
       CALL w3movdat(rinc,idat,jdat)
    endif
    !
    jdow = 0
    jdoy = 0
    jday = 0
    call w3doxdat(jdat,jdow,jdoy,jday)
    rjday = jdoy + jdat(5) / 24.
    IF (RJDAY < oz_time(1)) RJDAY = RJDAY + 365.
    !
    n2 = timeoz + 1
    do j=2,timeoz
       if (rjday < oz_time(j)) then
          n2 = j
          exit
       endif
    enddo
    n1 = n2 - 1
    
    tx1 = (oz_time(n2) - rjday) / (oz_time(n2) - oz_time(n1))
    tx2 = 1.0 - tx1
    
    if (n2 > timeoz) n2 = n2 - timeoz
    !
    do nc=1,oz_coeff
       do L=1,levozp
          do J=1,npts
             J1  = jindx1(J)
             J2  = jindx2(J)
             TEM = 1.0 - ddy(J)
             ozplout(j,L,nc) = tx1*(TEM*ozplin(J1,L,nc,n1)+ddy(J)*ozplin(J2,L,nc,n1)) & 
                  + tx2*(TEM*ozplin(J1,L,nc,n2)+ddy(J)*ozplin(J2,L,nc,n2))
          enddo
       enddo
    enddo

    !
    return

  end subroutine ozphys_2015_timestep_init

! ###########################################################################################  
!> The operational GFS currently parameterizes ozone production and
!! destruction based on monthly mean coefficients (
!! \c ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77) provided by Naval
!! Research Laboratory through CHEM2D chemistry model
!! (McCormack et al. (2006) \cite mccormack_et_al_2006).
!! \section arg_table_ozphys_2015_run Argument Table
!! \htmlinclude ozphys_2015_run.html
!!
!> \section genal_ozphys_2015 GFS ozphys_2015_run General Algorithm
!> -  This code assumes that both prsl and po3 are from bottom to top
!!     as are all other variables.
!> -  This code is specifically for NRL parameterization and
!!     climatological T and O3 are in location 5 and 6 of prdout array
!!\author June 2015 - Shrinivas Moorthi
!!\author May 2023  - Dustin Swales
! ###########################################################################################
  subroutine ozphys_2015_run ( im, levs, ko3, dt, oz, tin, po3, prsl, prdout, pl_coeff,     &
       delp, ldiag3d, dtend, dtidx, ntoz, index_of_process_prod_loss,                       &
       index_of_process_ozmix, index_of_process_temp, index_of_process_overhead_ozone,      &
       con_g, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         ldiag3d                         ! Flag to output GFS diagnostic tendencies
    real(kind_phys),intent(in) :: &
         con_g                           ! Physical constant: Gravitational acceleration (ms-2)
    integer, intent(in) :: &
         im,                           & ! Horizontal dimension
         levs,                         & ! Number of vertical layers
         ko3,                          & ! Number of vertical layers in ozone forcing data
         pl_coeff,                     & ! Number of coefficients in ozone forcing data
         ntoz,                         & ! Index for ozone mixing ratio
         index_of_process_prod_loss,   & ! Index for process in diagnostic tendency output
         index_of_process_ozmix,       & ! Index for process in diagnostic tendency output
         index_of_process_temp,        & ! Index for process in diagnostic tendency output
         index_of_process_overhead_ozone ! Index for process in diagnostic tendency output
    integer, intent(in), dimension(:,:) :: &
         dtidx                           ! Bookkeeping indices for GFS diagnostic tendencies
    real(kind_phys), intent(in) :: &
         dt                              ! Physics timestep (seconds)
    real(kind_phys), intent(in), dimension(:) :: &
         po3                             ! Natural log of ozone forcing data pressure levels
    real(kind_phys), intent(in), dimension(:,:) :: &
         prsl,                         & ! Air-pressure (Pa)
         tin,                          & ! Temperature of new-state (K)
         delp                            ! Difference between mid-layer pressures (Pa)
    real(kind_phys), intent(in), dimension(:,:,:) :: &
         prdout                          ! Ozone forcing data

    ! In/Outs
    real(kind=kind_phys), intent(inout), dimension(:,:,:) :: &
         dtend                           ! Diagnostic tendencies for state variables

    ! Outputs
    real(kind=kind_phys), intent(inout), dimension(:,:) :: &
         oz                              ! Ozone concentration updated by physics
    character(len=*), intent(out) :: &
         errmsg                          ! CCPP error message
    integer, intent(out) :: &
         errflg                          ! CCPP error flag

    ! Locals
    integer :: k, kmax, kmin, l, i, j
    integer, dimension(4) :: idtend
    logical, dimension(im) :: flg
    real :: gravi
    real(kind_phys) :: pmax, pmin, tem, temp
    real(kind_phys), dimension(im) :: wk1, wk2, wk3, ozib
    real(kind_phys), dimension(im,pl_coeff) :: prod
    real(kind_phys), dimension(im,levs) :: ozi
    real(kind_phys), dimension(im,levs+1) :: colo3, coloz

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Are UFS diagnostic tendencies requested? If so, set up bookeeping indices...
    if(ldiag3d) then
       idtend(1) = dtidx(100+ntoz,index_of_process_prod_loss)          ! was ozp1
       idtend(2) = dtidx(100+ntoz,index_of_process_ozmix)              ! was ozp2
       idtend(3) = dtidx(100+ntoz,index_of_process_temp)               ! was ozp3
       idtend(4) = dtidx(100+ntoz,index_of_process_overhead_ozone)     ! was ozp4
    else
       idtend=0
    endif

    ! Temporaries
    ozi = oz
    gravi=1.0/con_g

    colo3(:,levs+1) = 0.0
    coloz(:,levs+1) = 0.0

    do l=levs,1,-1
       pmin =  1.0e10
       pmax = -1.0e10

       do i=1,im
          wk1(i) = log(prsl(i,l))
          pmin   = min(wk1(i), pmin)
          pmax   = max(wk1(i), pmax)
          prod(i,:) = 0.0
       enddo
       kmax = 1
       kmin = 1
       do k=1,ko3-1
          if (pmin < po3(k)) kmax = k
          if (pmax < po3(k)) kmin = k
       enddo
       !
       do k=kmin,kmax
          temp = 1.0 / (po3(k) - po3(k+1))
          do i=1,im
             flg(i) = .false.
             if (wk1(i) < po3(k) .and. wk1(i) >= po3(k+1)) then
                flg(i) = .true.
                wk2(i) = (wk1(i) - po3(k+1)) * temp
                wk3(i) = 1.0 - wk2(i)
             endif
          enddo
          do j=1,pl_coeff
             do i=1,im
                if (flg(i)) then
                   prod(i,j)  = wk2(i) * prdout(i,k,j) + wk3(i) * prdout(i,k+1,j)
                endif
             enddo
          enddo
       enddo

       do j=1,pl_coeff
          do i=1,im
             if (wk1(i) < po3(ko3)) then
                prod(i,j) = prdout(i,ko3,j)
             endif
             if (wk1(i) >= po3(1)) then
                prod(i,j) = prdout(i,1,j)
             endif
          enddo
       enddo
       do i=1,im
          colo3(i,l) = colo3(i,l+1) + ozi(i,l)  * delp(i,l)*gravi
          coloz(i,l) = coloz(i,l+1) + prod(i,6) * delp(i,l)*gravi
          prod(i,2)  = min(prod(i,2), 0.0)
       enddo
       do i=1,im
          ozib(i)  = ozi(i,l)            ! no filling
          tem      = prod(i,1) - prod(i,2) * prod(i,6) + prod(i,3) * (tin(i,l) - prod(i,5)) &
                                                       + prod(i,4) * (colo3(i,l)-coloz(i,l))
          oz(i,l) = (ozib(i)  + tem*dt) / (1.0 - prod(i,2)*dt)
       enddo
       if(idtend(1)>=1) then
          dtend(:,l,idtend(1)) = dtend(:,l,idtend(1)) + (prod(:,1)-prod(:,2)*prod(:,6))*dt
       endif
       if(idtend(2)>=1) then
          dtend(:,l,idtend(2)) = dtend(:,l,idtend(2)) + (oz(:,l) - ozib(:))
       endif
       if(idtend(3)>=1) then
          dtend(:,l,idtend(3)) = dtend(:,l,idtend(3)) + prod(:,3)*(tin(:,l)-prod(:,5))*dt
       endif
       if(idtend(4)>=1) then
          dtend(:,l,idtend(4)) = dtend(:,l,idtend(4)) + prod(:,4) * (colo3(:,l)-coloz(:,l))*dt
       endif
    enddo

    return
  end subroutine ozphys_2015_run
!> @}
end module ozphys_2015

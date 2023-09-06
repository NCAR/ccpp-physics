! ###########################################################################################
!> \file ozphys_time_vary.F90
!!
! ###########################################################################################
module ozphys_time_vary
  use machine, only : kind_phys, kind_dbl_prec, kind_sngl_prec
  implicit none
  public ozphys_time_vary_init, ozphys_time_vary_timestep_init
contains

! ###########################################################################################
!>\defgroup GFS Ozone Data Module
!! This module updates the ozone data used by physics.
!> @{
!> \section arg_table_ozphys_time_vary_init Argument Table
!! \htmlinclude ozphys_time_vary_init.html
!!
! ###########################################################################################
  subroutine ozphys_time_vary_init(nPts, latsozp, oz_lat, dlat, jindx1, jindx2, ddy,        &
       errmsg, errflg)
    ! Inputs
    integer, intent(in) :: &
         nPts,      & ! Horizontal dimension
         latsozp      ! Number of latitudes in ozone data
    real(kind_phys),  intent(in), dimension(:) :: &
         oz_lat,    & ! Latitudes of ozone data
         dlat         ! Latitudes of grid
    ! Outputs
    integer, intent(inout), dimension(:) :: &
         jindx1,    & ! Interpolation index (low) for ozone data
         jindx2       ! Interpolation index (high) for ozone data
    real(kind_phys), intent(inout), dimension(:) :: &
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
    
  end subroutine ozphys_time_vary_init

! ###########################################################################################
!> \section arg_table_ozphys_time_vary_timestep_init Argument Table
!! \htmlinclude ozphys_time_vary_timestep_init.html
!!
! ###########################################################################################
  subroutine ozphys_time_vary_timestep_init(nPts, idate, fhour, jindx1, jindx2, latsozp,    &
       levozp, oz_coeff, timeoz, ozplin, oz_time, oz_lat, ddy, oz_data, errmsg, errflg)
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
    real(kind_phys), intent(inout), dimension(:,:,:) :: &
         oz_data     ! Ozone forcing data
    character(len=*), intent(out) :: &
         errmsg      ! CCPP error message
    integer,          intent(out) :: &
         errflg      ! CCPP error flag

    ! Local
    integer :: idat(8),jdat(8),iday,j,j1,j2,l,nc,n1,n2,jdow,jdoy,&
         jday,w3kindreal,w3kindint
    real(kind_phys) :: tem, tx1, tx2, rjday
    real(kind_dbl_prec) :: rinc(5)
    real(kind_sngl_prec) :: rinc4(5)

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
             oz_data(j,L,nc) = tx1*(TEM*ozplin(J1,L,nc,n1)+ddy(J)*ozplin(J2,L,nc,n1)) & 
                  + tx2*(TEM*ozplin(J1,L,nc,n2)+ddy(J)*ozplin(J2,L,nc,n2))
          enddo
       enddo
    enddo

    return

  end subroutine ozphys_time_vary_timestep_init
!> @}
end module ozphys_time_vary

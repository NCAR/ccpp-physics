!>\file GFS_rad_time_vary.mpas.F90
!!  Contains code related to GFS radiation suite setup (radiation part of time_vary_step)
module GFS_rad_time_vary
  implicit none

  private

  public GFS_rad_time_vary_timestep_init

contains

!> This module contains code related to GFS radiation setup.

!> \section arg_table_GFS_rad_time_vary_timestep_init Argument Table
!! \htmlinclude GFS_rad_time_vary_timestep_init.html
!!
  subroutine GFS_rad_time_vary_timestep_init (lrseeds, rseeds, lslwr, lsswr, isubc_lw, &
       isubc_sw, icsdsw, icsdlw, sec, kdt,  ipsd0, ipsdlim, errmsg, errflg)
    use mersenne_twister, only: random_setseed, random_index, random_stat
    use machine,          only: kind_phys
    use radcons,          only: con_100
    implicit none

    ! Interface variables
    logical,                intent(in)    :: lrseeds
    integer,                intent(in), optional :: rseeds(:,:)
    integer,                intent(in)    :: isubc_lw, isubc_sw, kdt
    integer,                intent(in)    :: ipsd0, ipsdlim
    logical,                intent(in)    :: lslwr, lsswr
    integer,                intent(inout), optional :: icsdsw(:), icsdlw(:)
    real(kind_phys),        intent(in)    :: sec
    character(len=*),       intent(out)   :: errmsg
    integer,                intent(out)   :: errflg

    ! Local variables
    type (random_stat) :: stat
    integer :: ix, j, i, ipseed, ixx
    integer, allocatable, dimension(:) :: numrdm

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (lsswr .or. lslwr) then
       !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
       if ((isubc_lw==2) .or. (isubc_sw==2)) then
          !NRL If random seeds supplied by NEPTUNE
          if(lrseeds) then
             do ix=1,size(icsdsw)
                icsdsw(ix) = rseeds(ix,1)
                icsdlw(ix) = rseeds(ix,2)
             enddo
          else
             allocate(numrdm(size(icsdlw)*2))
             ipseed = mod(nint(con_100*sqrt(sec)), ipsdlim) + 1 + ipsd0
             call random_setseed (ipseed, stat)
             call random_index (ipsdlim, numrdm, stat)

             ixx = 1
             do ix=1,size(icsdsw)*2,2
                icsdsw(ixx) = numrdm(ix)
                icsdlw(ixx) = numrdm(ix+1)
                ixx = ixx + 1
             enddo
             deallocate(numrdm)
          end if ! lrseeds
       endif     ! isubc_lw and isubc_sw
    endif        ! lsswr or lslwr

  end subroutine GFS_rad_time_vary_timestep_init

end module GFS_rad_time_vary

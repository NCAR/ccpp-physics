!> \file GFS_zhao_carr_pre.f90
!! This file contains the subroutines that calculates physics/diagnotics variables 
!! before calling microphysics scheme:

      module GFS_zhao_carr_pre
      contains 

!> \defgroup GFS_zhao_carr_pre GFS Zhao-Carr Scheme pre
!! @{
!!\section arg_table_GFS_zhao_carr_pre_init Argument Table
!!
      subroutine GFS_zhao_carr_pre_init     
      end subroutine GFS_zhao_carr_pre_init


!!\section arg_table_GFS_zhao_carr_pre_run Argument Table
!!| local var name | longname                                               |description                                               | units       | rank |  type   |   kind    | intent | optional |
!!|----------------|--------------------------------------------------------|----------------------------------------------------------|-------------|------|---------|-----------|--------|----------|
!!|   im           | horizontal_loop_extent                                 | horizontal loop extent, start at 1                       | index       | 0    | integer |           | in     |  F       |
!!|   ix           | horizontal_dimension                                   | horizontal dimension                                     | index       | 0    | integer |           | in     |  F       |
!!|   levs         | vertical_dimension                                     | vertical layer dimension                                 | index       | 0    | integer |           | in     |  F       |
!!|   cwm          | cloud_condensed_water_specific_humidity                | cloud condensed water specific humidity                  | kg kg-1     | 2    | real    | kind_phys | in     |  F       |
!!|   clw1         | cloud_ice_specific_humidity                            | cloud ice specific humidity                              | kg kg-1     | 2    | real    | kind_phys | out    |  F       |
!!
      subroutine GFS_zhao_carr_pre_run (im, ix, levs, cwm, clw1 )  
     
!
      use machine,               only: kind_phys

      implicit none
!    
!     declare variables.
!
      integer,intent(in)   :: im, ix, levs
      integer              :: i,k
      real(kind=kind_phys),dimension(ix,levs), intent(in)  ::  cwm 
      real(kind=kind_phys),dimension(ix,levs), intent(out) ::  clw1

            do i = 1, im
              do k = 1, levs
                   clw1(i,k) = cwm(i,k)                     !Stateout%gq0(:,:,Model%ntcw)
              enddo
            enddo


      end subroutine GFS_zhao_carr_pre_run

!!\section arg_table_GFS_zhao_carr_pre_finalize Argument Table
!!
      subroutine GFS_zhao_carr_pre_finalize
      end subroutine GFS_zhao_carr_pre_finalize
!! @}
      end module GFS_zhao_carr_pre


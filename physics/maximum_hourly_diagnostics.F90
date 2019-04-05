module maximum_hourly_diagnostics

   use machine, only: kind_phys

   implicit none

   private

   public maximum_hourly_diagnostics_init, maximum_hourly_diagnostics_run, maximum_hourly_diagnostics_finalize

   ! DH* TODO - THIS CAME FROM PHYSCONS.F90 BUT IS IT BETTER PLACED IN HERE?
   real(kind=kind_phys), parameter ::PQ0=379.90516E0, A2A=17.2693882, A3=273.16, A4=35.86, RHmin=1.0E-6
   ! *DH

contains

   subroutine maximum_hourly_diagnostics_init()
   end subroutine maximum_hourly_diagnostics_init

   subroutine maximum_hourly_diagnostics_finalize()
   end subroutine maximum_hourly_diagnostics_finalize

#if 0
!> \section arg_table_maximum_hourly_diagnostics_run Argument Table
!! | local_name           | standard_name                                                      | long_name                                                          | units      | rank | type      | kind      | intent | optional |
!! |----------------------|--------------------------------------------------------------------|--------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im                   | horizontal_loop_extent                                             | horizontal loop extent                                             | count      |    0 | integer   |           | in     | F        |
!! | levs                 | vertical_dimension                                                 | number of vertical levels                                          | count      |    0 | integer   |           | in     | F        |
!! | kdt                  | index_of_time_step                                                 | current forecast iteration                                         | index      |    0 | integer   |           | in     | F        |
!! | nsteps_per_reset     | number_of_time_steps_per_maximum_hourly_time_interval              | number_of_time_steps_per_maximum_hourly_time_interval              | count      |    0 | integer   |           | in     | F        |
!! | lradar               | flag_for_radar_reflectivity                                        | flag for radar reflectivity                                        | flag       |    0 | logical   |           | in     | F        |
!! | imp_physics          | flag_for_microphysics_scheme                                       | choice of microphysics scheme                                      | flag       |    0 | integer   |           | in     | F        |
!! | imp_physics_gfdl     | flag_for_gfdl_microphysics_scheme                                  | choice of GFDL microphysics scheme                                 | flag       |    0 | integer   |           | in     | F        |
!! | imp_physics_thompson | flag_for_thompson_microphysics_scheme                              | choice of Thompson microphysics scheme                             | flag       |    0 | integer   |           | in     | F        |
!! | con_g                | gravitational_acceleration                                         | gravitational acceleration                                         | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | phil                 | geopotential                                                       | geopotential at model layer centers                                | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | gt0                  | air_temperature_updated_by_physics                                 | temperature updated by physics                                     | K          |    2 | real      | kind_phys | in     | F        |
!! | refl_10cm            | radar_reflectivity_10cm                                            | instantaneous refl_10cm                                            | dBZ        |    2 | real      | kind_phys | in     | F        |
!! | refdmax              | maximum_reflectivity_at_1km_agl_over_maximum_hourly_time_interval  | maximum reflectivity at 1km agl over maximum hourly time interval  | dBZ        |    1 | real      | kind_phys | inout  | F        |
!! | refdmax263k          | maximum_reflectivity_at_minus10c_over_maximum_hourly_time_interval | maximum reflectivity at minus10c over maximum hourly time interval | dBZ        |    1 | real      | kind_phys | inout  | F        |
!! | u10m                 | x_wind_at_10m                                                      | 10 meter u wind speed                                              | m s-1      |    1 | real      | kind_phys | in     | F        |
!! | v10m                 | y_wind_at_10m                                                      | 10 meter v wind speed                                              | m s-1      |    1 | real      | kind_phys | in     | F        |
!! | u10max               | maximum_u_wind_at_10m_over_maximum_hourly_time_interval            | maximum u wind at 10m over maximum hourly time interval            | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | v10max               | maximum_v_wind_at_10m_over_maximum_hourly_time_interval            | maximum v wind at 10m over maximum hourly time interval            | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | spd10max             | maximum_wind_at_10m_over_maximum_hourly_time_interval              | maximum wind at 10m over maximum hourly time interval              | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | pgr                  | surface_air_pressure                                               | surface pressure                                                   | Pa         |    1 | real      | kind_phys | in     | F        |
!! | t2m                  | temperature_at_2m                                                  | 2 meter temperature                                                | K          |    1 | real      | kind_phys | in     | F        |
!! | q2m                  | specific_humidity_at_2m                                            | 2 meter specific humidity                                          | kg kg-1    |    1 | real      | kind_phys | in     | F        |
!! | t02max               | maximum_temperature_at_2m_over_maximum_hourly_time_interval        | maximum temperature at 2m over maximum hourly time interval        | K          |    1 | real      | kind_phys | inout  | F        |
!! | t02min               | minimum_temperature_at_2m_over_maximum_hourly_time_interval        | minumum temperature at 2m over maximum hourly time interval        | K          |    1 | real      | kind_phys | inout  | F        |
!! | rh02max              | maximum_relative_humidity_at_2m_over_maximum_hourly_time_interval  | maximum relative humidity at 2m over maximum hourly time interval  | %          |    1 | real      | kind_phys | inout  | F        |
!! | rh02min              | minimum_relative_humidity_at_2m_over_maximum_hourly_time_interval  | minumum relative humidity at 2m over maximum hourly time interval  | %          |    1 | real      | kind_phys | inout  | F        |
!! | errmsg               | ccpp_error_message                                                 | error message for error handling in CCPP                           | none       |    0 | character | len=*     | out    | F        |
!! | errflg               | ccpp_error_flag                                                    | error flag for error handling in CCPP                              | flag       |    0 | integer   |           | out    | F        |
!!
#endif
   subroutine maximum_hourly_diagnostics_run(im, levs, kdt, nsteps_per_reset, lradar, imp_physics, &
                                             imp_physics_gfdl, imp_physics_thompson, con_g, phil,  &
                                             gt0, refl_10cm, refdmax, refdmax263k, u10m, v10m,     &
                                             u10max, v10max, spd10max, pgr, t2m, q2m, t02max,      &
                                             t02min, rh02max, rh02min, errmsg, errflg)

       ! Interface variables
       integer, intent(in) :: im, levs, kdt, nsteps_per_reset
       logical, intent(in) :: lradar
       integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson
       real(kind_phys), intent(in   ) :: con_g
       real(kind_phys), intent(in   ) :: phil(im,levs)
       real(kind_phys), intent(in   ) :: gt0(im,levs)
       real(kind_phys), intent(in   ) :: refl_10cm(im,levs)
       real(kind_phys), intent(inout) :: refdmax(im)
       real(kind_phys), intent(inout) :: refdmax263k(im)
       real(kind_phys), intent(in   ) :: u10m(im)
       real(kind_phys), intent(in   ) :: v10m(im)
       real(kind_phys), intent(inout) :: u10max(im)
       real(kind_phys), intent(inout) :: v10max(im)
       real(kind_phys), intent(inout) :: spd10max(im)
       real(kind_phys), intent(in   ) :: pgr(im)
       real(kind_phys), intent(in   ) :: t2m(im)
       real(kind_phys), intent(in   ) :: q2m(im)
       real(kind_phys), intent(inout) :: t02max(im)
       real(kind_phys), intent(inout) :: t02min(im)
       real(kind_phys), intent(inout) :: rh02max(im)
       real(kind_phys), intent(inout) :: rh02min(im)
       character(len=*), intent(out)  :: errmsg
       integer, intent(out)           :: errflg

       ! Local variables
       real(kind_phys), dimension(:), allocatable :: refd, refd263k
       real(kind_phys) :: tem, pshltr, QCQ, rh02
       integer :: kdtminus1, i

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       kdtminus1 = kdt-1

!Calculate hourly max 1-km agl and -10C reflectivity
       if (lradar .and. (imp_physics == imp_physics_gfdl .or. imp_physics == imp_physics_thompson)) then
          allocate(refd(im))
          allocate(refd263k(im))
          call max_fields(phil,refl_10cm,con_g,im,levs,refd,gt0,refd263k)
          if(mod(kdtminus1,nsteps_per_reset)==0)then
             do i=1,im
               refdmax(i) = -35.
               refdmax263k(i) = -35.
             enddo
          endif
          do i=1,im
             !if(mod(kdtminus1,nsteps_per_reset)==0)then
             !  refdmax(I) = -35.
             !  refdmax263k(I) = -35.
             !endif
             refdmax(i) = max(refdmax(i),refd(i))
             refdmax263k(i) = max(refdmax263k(i),refd263k(i))
          enddo
          deallocate (refd) 
          deallocate (refd263k)
       endif
!
       if(mod(kdtminus1,nsteps_per_reset)==0)then
          do i=1,im
             spd10max(i) = -999.
             u10max(i)   = -999.
             v10max(i)   = -999.
             t02max(i)   = -999.
             t02min(i)   = 999.
             rh02max(i)  = -999.
             rh02min(i)  = 999.
          enddo
       endif
       do i=1,im
! find max hourly wind speed then decompose
          tem = sqrt(u10m(i)*u10m(i) + v10m(i)*v10m(i))
          !if(mod(kdtminus1,nsteps_per_reset)==0)then
          !   spd10max(i) = -999.
          !   u10max(i)   = -999.
          !   v10max(i)   = -999.
          !   t02max(i)   = -999.
          !   t02min(i)   = 999.
          !   rh02max(i)  = -999.
          !   rh02min(i)  = 999.
          !endif
          if (tem > spd10max(i)) then
             spd10max(i) = tem
             u10max(i)   = u10m(i)
             v10max(i)   = v10m(i)
          endif
          pshltr=pgr(i)*exp(-0.068283/gt0(i,1))
          QCQ=PQ0/pshltr*EXP(A2A*(t2m(i)-A3)/(t2m(i)-A4))
          rh02=q2m(i)/QCQ
          IF (rh02.GT.1.0) THEN
             rh02=1.0
          ENDIF
          IF (rh02.LT.RHmin) THEN !use smaller RH limit for stratosphere
             rh02=RHmin
          ENDIF
          rh02max(i)=max(rh02max(i),rh02)
          rh02min(i)=min(rh02min(i),rh02)
          t02max(i)=max(t02max(i),t2m(i))  !<--- hourly max 2m t
          t02min(i)=min(t02min(i),t2m(i))  !<--- hourly min 2m t
       enddo

   end subroutine maximum_hourly_diagnostics_run

   subroutine max_fields(phil,ref3D,grav,im,levs,refd,tk,refd263k)
      integer, intent(in)               :: im,levs
      real (kind=kind_phys), intent(in) :: grav
      real (kind=kind_phys), intent(in),dimension(im,levs)  :: phil,ref3D,tk
      integer               :: i,k,ll,ipt,kpt
      real :: dbz1avg,zmidp1,zmidloc,refl,fact
      real, dimension(im,levs) :: z
      real, dimension(im) :: zintsfc
      real, dimension(im), intent(inout) :: refd,refd263k
      REAL :: dbz1(2),dbzk,dbzk1
      logical :: counter
      do i=1,im
         do k=1,levs
            z(i,k)=phil(i,k)/grav
         enddo
      enddo
      do i=1,im
         refd(I) = -35.
         vloop:  do k=1,levs-1
            if ( (z(i,k+1)) .ge. 1000.     &
             .and.(z(i,k))   .le. 1000.)  then
               zmidp1=z(i,k+1)
               zmidLOC=z(i,k)
               dbz1(1)=ref3d(i,k+1)   !- dBZ (not Z) values
               dbz1(2)=ref3d(i,k) !- dBZ values
               exit vloop
            endif
         enddo vloop

!!! Initial curefl value without reduction above freezing level
!
!         curefl=0.
!         if (cprate(i,j)>0.) then
!           cuprate=rdtphs*cprate(i,j)
!           curefl=cu_a*cuprate**cu_b
!         endif
         do ll=1,2
           refl=0.
           if (dbz1(ll)>-35.) refl=10.**(0.1*dbz1(ll))
!           dbz1(l)=curefl+refl    !- in Z units
             dbz1(ll)=refl
         enddo
!-- Vertical interpolation of Z (units of mm**6/m**3)
         fact=(1000.-zmidloc)/(zmidloc-zmidp1)
         dbz1avg=dbz1(2)+(dbz1(2)-dbz1(1))*fact
!-- Convert to dBZ (10*logZ) as the last step
         if (dbz1avg>0.01) then
           dbz1avg=10.*alog10(dbz1avg)
         else
           dbz1avg=-35.
         endif
         refd(I)=max(refd(I),dbz1avg)
      enddo

!-- refl at -10C
      do i=1,im
         dbz1(1) = -35.
         dbz1(2) = -35.
         vloopm10:  do k=1,levs-1
            if (tk(i,k+1) .le. 263.15 .and. tk(i,k) .ge. 263.15)  then
               dbz1(1)=ref3d(i,k+1)   !- dBZ (not Z) values
               dbz1(2)=ref3d(i,k) !- dBZ values
               exit vloopm10
            endif
         enddo vloopm10

         do ll=1,2
           refl=0.
           if (dbz1(ll)>-35.) refl=10.**(0.1*dbz1(ll))
!           dbz1(l)=curefl+refl    !- in Z units
             dbz1(ll)=refl
         enddo
!-- Take max of bounding reflectivity values 
         dbz1avg=maxval(dbz1)
!-- Convert to dBZ (10*logZ) as the last step
         if (dbz1avg>0.01) then
           dbz1avg=10.*alog10(dbz1avg)
         else
           dbz1avg=-35.
         endif
         refd263K(I)=dbz1avg
      enddo
   end subroutine max_fields

end module maximum_hourly_diagnostics
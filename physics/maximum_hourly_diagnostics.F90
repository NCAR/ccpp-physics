module maximum_hourly_diagnostics

   use machine, only: kind_phys

   implicit none

   private

   public maximum_hourly_diagnostics_init, maximum_hourly_diagnostics_run, maximum_hourly_diagnostics_finalize

   ! DH* TODO - cleanup use of constants
   real(kind=kind_phys), parameter ::PQ0=379.90516E0, A2A=17.2693882, A3=273.16, A4=35.86, RHmin=1.0E-6
   ! *DH

contains

   subroutine maximum_hourly_diagnostics_init()
   end subroutine maximum_hourly_diagnostics_init

   subroutine maximum_hourly_diagnostics_finalize()
   end subroutine maximum_hourly_diagnostics_finalize

#if 0
!> \section arg_table_maximum_hourly_diagnostics_run Argument Table
!! \htmlinclude maximum_hourly_diagnostics_run.html
!!
#endif
   subroutine maximum_hourly_diagnostics_run(im, levs, reset, lradar, imp_physics,                 &
                                             imp_physics_gfdl, imp_physics_thompson,               &
                                             imp_physics_fer_hires,imp_physics_nssl,con_g, phil,   &
                                             gt0, refl_10cm, refdmax, refdmax263k, u10m, v10m,     &
                                             u10max, v10max, spd10max, pgr, t2m, q2m, t02max,      &
                                             t02min, rh02max, rh02min, errmsg, errflg)

       ! Interface variables
       integer, intent(in) :: im, levs
       logical, intent(in) :: reset, lradar
       integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_fer_hires, imp_physics_nssl
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
       integer :: i

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

!Calculate hourly max 1-km agl and -10C reflectivity
       if (lradar .and. (imp_physics == imp_physics_gfdl .or. &
          imp_physics == imp_physics_thompson .or. imp_physics == imp_physics_fer_hires &
          .or. imp_physics == imp_physics_nssl)) then

          allocate(refd(im))
          allocate(refd263k(im))
          call max_fields(phil,refl_10cm,con_g,im,levs,refd,gt0,refd263k)
          if (reset) then
             do i=1,im
               refdmax(i) = -35.
               refdmax263k(i) = -35.
             enddo
          endif
          do i=1,im
             refdmax(i) = max(refdmax(i),refd(i))
             refdmax263k(i) = max(refdmax263k(i),refd263k(i))
          enddo
          deallocate (refd) 
          deallocate (refd263k)
       endif
!
       if (reset) then
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
          if (tem > spd10max(i)) then
             spd10max(i) = tem
             u10max(i)   = u10m(i)
             v10max(i)   = v10m(i)
          endif
          pshltr=pgr(i)*exp(-0.068283/gt0(i,1))
          QCQ=PQ0/pshltr*EXP(A2A*(t2m(i)-A3)/(t2m(i)-A4))
          rh02=q2m(i)/QCQ
          IF (rh02 > 1.0) THEN
             rh02 = 1.0
          ENDIF
          IF (rh02 < RHmin) THEN !use smaller RH limit for stratosphere
             rh02 = RHmin
          ENDIF
          rh02max(i) = max(rh02max(i),rh02)
          rh02min(i) = min(rh02min(i),rh02)
          t02max(i)  = max(t02max(i),t2m(i))  !<--- hourly max 2m t
          t02min(i)  = min(t02min(i),t2m(i))  !<--- hourly min 2m t
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
            if ( z(i,k+1) >= 1000. .and. z(i,k) <= 1000.)  then
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

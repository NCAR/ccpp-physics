!>\file maximum_hourly_diagnostics.F90
!!

module maximum_hourly_diagnostics

   use machine, only: kind_phys

   implicit none

   private

   public  maximum_hourly_diagnostics_run

   ! DH* TODO - cleanup use of constants
   real(kind=kind_phys), parameter ::PQ0=379.90516E0, A2A=17.2693882, A3=273.16, A4=35.86, RHmin=1.0E-6
   ! *DH

   ! Conversion from flashes per five minutes to flashes per minute.
   real(kind=kind_phys), parameter :: scaling_factor = 0.2

contains

#if 0
!> \section arg_table_maximum_hourly_diagnostics_run Argument Table
!! \htmlinclude maximum_hourly_diagnostics_run.html
!!
#endif
   subroutine maximum_hourly_diagnostics_run(im, levs, reset, lradar, imp_physics,                 &
                                             imp_physics_gfdl, imp_physics_thompson,               &
                                             imp_physics_fer_hires, imp_physics_nssl,              &
                                             con_g, phil,                                          &
                                             gt0, refl_10cm, refdmax, refdmax263k, u10m, v10m,     &
                                             u10max, v10max, spd10max, pgr, t2m, q2m, t02max,      &
                                             t02min, rh02max, rh02min, dtp, rain, pratemax,        &
                                             lightning_threat, ltg1_max,ltg2_max,ltg3_max,         &
                                             wgrs, prsi, qgraupel, qsnowwat, qicewat, tgrs, con_rd,&
                                             prsl, kdt, errmsg, errflg)

       ! Interface variables
       integer, intent(in) :: im, levs, kdt
       logical, intent(in) :: reset, lradar, lightning_threat
       integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_fer_hires, &
                              imp_physics_nssl
       real(kind_phys), intent(in   ) :: con_g
       real(kind_phys), intent(in   ) :: con_rd
       real(kind_phys), intent(in   ) :: phil(:,:)
       real(kind_phys), intent(in   ) :: gt0(:,:)
       real(kind_phys), intent(in   ) :: refl_10cm(:,:)
       real(kind_phys), intent(inout) :: refdmax(:)
       real(kind_phys), intent(inout) :: refdmax263k(:)
       real(kind_phys), intent(in   ) :: u10m(:)
       real(kind_phys), intent(in   ) :: v10m(:)
       real(kind_phys), intent(inout) :: u10max(:)
       real(kind_phys), intent(inout) :: v10max(:)
       real(kind_phys), intent(inout) :: spd10max(:)
       real(kind_phys), intent(in   ) :: pgr(:)
       real(kind_phys), intent(in   ) :: t2m(:)
       real(kind_phys), intent(in   ) :: q2m(:)
       real(kind_phys), intent(inout) :: t02max(:)
       real(kind_phys), intent(inout) :: t02min(:)
       real(kind_phys), intent(inout) :: rh02max(:)
       real(kind_phys), intent(inout) :: rh02min(:)
       real(kind_phys), intent(in   ) :: dtp
       real(kind_phys), intent(in   ) :: rain(:)
       real(kind_phys), intent(in   ) :: tgrs(:,:)
       real(kind_phys), intent(in   ) :: prsl(:,:)
       real(kind_phys), intent(inout) :: pratemax(:)

       real(kind_phys), intent(in), dimension(:,:) :: prsi, qgraupel, qsnowwat, qicewat
       real(kind_phys), intent(in), dimension(:,:), optional :: wgrs
       real(kind_phys), intent(inout), dimension(:), optional :: ltg1_max, ltg2_max, ltg3_max
       character(len=*), intent(out)  :: errmsg
       integer, intent(out)           :: errflg

       ! Local variables
       real(kind_phys), dimension(:), allocatable :: refd, refd263k
       real(kind_phys) :: tem, pshltr, QCQ, rh02, dP, Q
       integer :: i

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

!Lightning threat indices
       if (lightning_threat) then
         call lightning_threat_indices
       endif

!Calculate hourly max 1-km agl and -10C reflectivity
       if (lradar .and. (imp_physics == imp_physics_gfdl .or. &
          imp_physics == imp_physics_thompson .or.    &
          imp_physics == imp_physics_fer_hires .or.   &
          imp_physics == imp_physics_nssl )) then
          allocate(refd(im))
          allocate(refd263k(im))
          call max_fields(phil,refl_10cm,con_g,im,levs,refd,gt0,refd263k)
          if (reset) then
             IF ( imp_physics == imp_physics_nssl ) THEN ! ERM: might not need this as a separate assignment
              do i=1,im
                refdmax(i) = 0.
                refdmax263k(i) = 0.
              enddo
            ELSE
              do i=1,im
                refdmax(i) = -35.
                refdmax263k(i) = -35.
              enddo
            ENDIF
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
             pratemax(i) = 0.
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
          pratemax(i) = max(pratemax(i),(3.6E6/dtp)*rain(i))
       enddo

     contains

       subroutine lightning_threat_indices
         implicit none
         REAL(kind_phys), PARAMETER    :: clim1=1.50
         REAL(kind_phys), PARAMETER    :: clim2=0.40*1.22
         REAL(kind_phys), PARAMETER    :: clim3=0.02*1.22
         !  coef1 and coef2 are modified from the values given
         !  in McCaul et al. 
         !  coef1 is x 1000 x 1.22
         !  coef2 is x 1.22
         !  are these tuning factors, scale factors??
         !  McCaul et al. used a 2-km WRF simulation
         REAL(kind_phys), PARAMETER    :: coef1=0.042*1000.*1.22
         REAL(kind_phys), PARAMETER    :: coef2=0.20*1.22
         
         REAL(kind_phys) :: totice_colint(im), ltg1, ltg2, high_ltg1, high_wgrs, high_graupel, rho
         LOGICAL :: ltg1_calc(im)
         integer :: k, i, count

         count = 0
         high_ltg1 = 0
         high_wgrs = 0
         high_graupel = 0

          totice_colint = 0
          ltg1_calc = .false.
          do k=1,levs-1
             do i=1,im
                dP = prsi(i,k) - prsi(i,k+1)
                Q = qgraupel(i,k) + qsnowwat(i,k) + qicewat(i,k)
                rho = prsl(i,k) / (con_rd * tgrs(i,k))
                totice_colint(i) = totice_colint(i) + Q * rho * dP / con_g

                IF ( .not.ltg1_calc(i) ) THEN
                   IF ( 0.5*(tgrs(i,k+1) + tgrs(i,k)) < 258.15 ) THEN
                      count = count + 1
                      ltg1_calc(i) = .true.
                      
                      ltg1 = coef1*wgrs(i,k)* &
                           (( qgraupel(i,k+1) + qgraupel(i,k) )*0.5 )
                      if(ltg1 > high_ltg1) then
                        high_ltg1 = ltg1
                        high_graupel = qgraupel(i,k)
                        high_wgrs = wgrs(i,k)
                      endif
                      
                      IF ( ltg1 .LT. clim1 ) ltg1 = 0.

                      ! Scale to flashes per minue
                      ltg1 = ltg1 * scaling_factor

                      IF ( ltg1 .GT. ltg1_max(i) ) THEN
                         ltg1_max(i) = ltg1
                      ENDIF
                   ENDIF
                ENDIF
             enddo
          enddo

          do i=1,im
             ltg2 = coef2 * totice_colint(i)

             IF ( ltg2 .LT. clim2 ) ltg2 = 0.

             ! Scale to flashes per minute
             ltg2 = ltg2 * scaling_factor
             
             IF ( ltg2 .GT. ltg2_max(i) ) THEN
                ltg2_max(i) = ltg2
             ENDIF

             ! This calculation is already in flashes per minute.
             ltg3_max(i) = 0.95 * ltg1_max(i) + 0.05 * ltg2_max(i)

             ! Thus, we must scale clim3. The compiler will optimize this away.
             IF ( ltg3_max(i) .LT. clim3 * scaling_factor ) ltg3_max(i) = 0.
          enddo

       end subroutine lightning_threat_indices

   end subroutine maximum_hourly_diagnostics_run

   subroutine max_fields(phil,ref3D,grav,im,levs,refd,tk,refd263k)
      integer, intent(in)               :: im,levs
      real (kind=kind_phys), intent(in) :: grav
      real (kind=kind_phys), intent(in),dimension(:,:)  :: phil,ref3D,tk
      integer               :: i,k,ll,ipt,kpt
      real :: dbz1avg,zmidp1,zmidloc,refl,fact
      real, dimension(im,levs) :: z
      real, dimension(im) :: zintsfc
      real, dimension(:), intent(inout) :: refd,refd263k
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

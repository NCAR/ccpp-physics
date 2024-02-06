!>\file  module_plumerise.F90
!! This file is the fire plume rise driver.

 module module_plumerise

  use machine , only : kind_phys
!  real(kind=kind_phys),parameter :: p1000mb = 100000.  ! p at 1000mb (pascals)
!- Implementing the fire radiative power (FRP) methodology for biomass burning
!- emissions and convective energy estimation.
!- Saulo Freitas, Gabriel Pereira (INPE/UFJS, Brazil)
!- Ravan Ahmadov, Georg Grell (NOAA, USA)
!- The flag "plumerise_flag" defines the method:
!-    =1 => original method
!-    =2 => FRP based

CONTAINS
subroutine ebu_driver (      flam_frac,ebu_in,ebu,                   &
                             theta_phy,q_vap,                        &   ! RAR: moist is replaced with q_vap, SRB: t_phy is repalced by theta_phy
                             rho_phy,vvel,u_phy,v_phy,pi_phy,        &   ! SRB: p_phy is replaced by pi_phy
                             wind_phy,                               &   ! SRB: added wind_phy
                             z_at_w,z,g,con_cp,con_rd,               &   ! scale_fire_emiss is part of config_flags
                             frp_inst, k_min, k_max,                 &   ! RAR:
                             wind_eff_opt,                           &
                             kpbl_thetav,                            &   ! SRB: added kpbl_thetav
                             ids,ide, jds,jde, kds,kde,              &
                             ims,ime, jms,jme, kms,kme,              &
                             its,ite, jts,jte, kts,kte, errmsg, errflg,curr_secs, &
                             xlat, xlong , uspdavg2, hpbl_thetav2, mpiid)

  use rrfs_smoke_config
  !use plume_data_mod
  USE module_zero_plumegen_coms
  USE module_smoke_plumerise
  IMPLICIT NONE

   REAL(kind_phys), PARAMETER :: frp_threshold= 1.e+7   ! Minimum FRP (Watts) to distribute smoke in PBL 
   
   REAL(kind_phys), PARAMETER :: zpbl_threshold   = 2000.    ! SRB: Minimum PBL depth to have plume rise 
   REAL(kind_phys), PARAMETER :: uspd_threshold   = 5.       ! SRB: Wind speed averaged across PBL depth to control smoke release levels 
   REAL(kind_phys), PARAMETER :: frp_threshold500 = 500.e+6  ! SRB: Minimum FRP (Watts) to have plume rise 
    
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme), INTENT(IN ) :: frp_inst         ! RAR: FRP array

   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  xlat,xlong ! SRB

   real(kind_phys), DIMENSION(ims:ime, jms:jme),  INTENT(IN)  ::     kpbl_thetav ! SRB  

   character(*), intent(inout) :: errmsg
   integer, intent(inout) :: errflg
   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   real(kind_phys) :: curr_secs
   INTEGER,      INTENT(IN   ) :: wind_eff_opt
   real(kind=kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT ) ::  ebu
   real(kind=kind_phys), INTENT(IN )  :: g, con_cp, con_rd
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme ), INTENT(IN )  :: ebu_in
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme ), INTENT(OUT ) :: flam_frac
   real(kind=kind_phys), DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   z,z_at_w,vvel,u_phy,v_phy,rho_phy,pi_phy,q_vap,theta_phy,wind_phy                     ! RAR, SRB

! Local variables...
      INTEGER :: nv, i, j, k,  kp1, kp2
      INTEGER :: icall=0
      INTEGER, DIMENSION(ims:ime, jms:jme), INTENT (OUT) :: k_min, k_max      ! Min and max ver. levels for BB injection spread
      REAL, DIMENSION(ims:ime, jms:jme), INTENT (OUT) :: uspdavg2, hpbl_thetav2 ! SRB
      real(kind_phys), dimension (kte) :: u_in ,v_in ,w_in ,theta_in ,pi_in, rho_phyin ,qv_in ,zmid, z_lev, uspd ! SRB
      real(kind=kind_phys) :: dz_plume, cpor, con_rocp, uspdavg ! SRB

! MPI variables
      INTEGER, INTENT(IN) :: mpiid

        cpor    =con_cp/con_rd
        con_rocp=con_rd/con_cp

        if (mod(int(curr_secs),1800) .eq. 0) then
            icall = 0
        endif

        IF ( dbg_opt .and. icall .le. n_dbg_lines) then
           WRITE(1000+mpiid,*) 'module_plumerise: its,ite,jts,jte ', its,ite,jts,jte
           WRITE(1000+mpiid,*) 'module_plumerise: ims,ime,jms,jme ', ims,ime,jms,jme
           WRITE(1000+mpiid,*) 'module_plumerise: maxval(ebu(:,kts,:)) ', maxval(ebu(:,kts,:))
        END IF

! RAR: setting to zero the ebu emissions at the levels k>1, this is necessary when the plumerise is called, so the emissions at k>1 are updated
       !do nv=1,num_ebu
          do j=jts,jte
            do k=kts,kte
               do i=its,ite
                 ebu(i,k,j)=0.
               enddo
            enddo
          enddo
       !enddo

! For now the flammable fraction is constant, based on the namelist. The next
! step to use LU index and meteorology to parameterize it
       do j=jts,jte
        do i=its,ite
           flam_frac(i,j)= 0.
           if (frp_inst(i,j) > frp_threshold) then 
              flam_frac(i,j)= 0.9
           end if
        enddo
       enddo


! RAR: new FRP based approach
! Haiqin: do_plumerise is added to the namelist options
check_pl:  IF (do_plumerise) THEN    ! if the namelist option is set for plumerise
       do j=jts,jte
          do i=its,ite

               do k=kts,kte
                  u_in(k)=  u_phy(i,k,j)
                  v_in(k)=  v_phy(i,k,j)
                  w_in(k)=  vvel(i,k,j)
                  qv_in(k)= q_vap(i,k,j)
                  pi_in(k)= pi_phy(i,k,j)
                  zmid(k)=  z(i,k,j)-z_at_w(i,kts,j)
                  z_lev(k)= z_at_w(i,k,j)-z_at_w(i,kts,j)
                  rho_phyin(k)= rho_phy(i,k,j)
                  theta_in(k)= theta_phy(i,k,j)
                  uspd(k)= wind_phy(i,k,j) ! SRB
               enddo


             IF (dbg_opt .and. (icall .le. n_dbg_lines) .and. (frp_inst(i,j) .ge. frp_threshold) ) then
               WRITE(1000+mpiid,*) 'module_plumerise_before:xlat,xlong,curr_secs,ebu(kts),frp_inst',xlat(i,j), xlong(i,j), int(curr_secs),ebu(i,kts,j),frp_inst(i,j)
               WRITE(1000+mpiid,*) 'module_plumerise_before:xlat,xlong,curr_secs,u(10),v(10),w(10),qv(10)',xlat(i,j), xlong(i,j),int(curr_secs), u_in(10),v_in(10),w_in(kte),qv_in(10)
             END IF

! RAR: the plume rise calculation step:
               CALL plumerise(kte,1,1,1,1,1,1,                      &
                              u_in, v_in, w_in, theta_in ,pi_in,    &
                              rho_phyin, qv_in, zmid, z_lev,        &
                              wind_eff_opt,                         &
                              frp_inst(i,j), k_min(i,j),            & 
                              k_max(i,j), dbg_opt, g, con_cp,       &
                              con_rd, cpor, errmsg, errflg,         &
                              icall, mpiid, xlat(i,j), xlong(i,j), curr_secs )
               if(errflg/=0) return

               kp1= k_min(i,j)
               kp2= k_max(i,j)   
               dz_plume= z_at_w(i,kp2,j) - z_at_w(i,kp1,j)

! SRB: Adding condition for overwriting plumerise levels               
               uspdavg=SUM(uspd(kts:kpbl_thetav(i,j)))/kpbl_thetav(i,j) !Average wind speed within the boundary layer

! SRB: Adding output
               uspdavg2(i,j) = uspdavg
               hpbl_thetav2(i,j) = z_lev(kpbl_thetav(i,j))
                
               IF ((frp_inst(i,j) .gt. frp_threshold) .AND. (frp_inst(i,j) .le. frp_threshold500) .AND. & 
                  (z_lev(kpbl_thetav(i,j)) .gt. zpbl_threshold) .AND. (wind_eff_opt .eq. 1)) THEN
                  kp1=1
                  IF (uspdavg .ge. uspd_threshold) THEN ! Too windy 
                     kp2=kpbl_thetav(i,j)/3 
                  ELSE
                     kp2=kpbl_thetav(i,j)
                  END IF
                  dz_plume= z_at_w(i,kp2,j) - z_at_w(i,kp1,j)
                  do k=kp1,kp2-1     
                     ebu(i,k,j)= ebu_in(i,j)* (z_at_w(i,k+1,j)-z_at_w(i,k,j))/dz_plume
                  enddo
               ELSE
                  do k=kp1,kp2-1     
                     ebu(i,k,j)= flam_frac(i,j)* ebu_in(i,j)* (z_at_w(i,k+1,j)-z_at_w(i,k,j))/dz_plume
                  enddo
                  ebu(i,kts,j)=   (1.-flam_frac(i,j))* ebu_in(i,j)
               END IF
! SRB: End modification 

               IF ( dbg_opt .and. (icall .le. n_dbg_lines)  .and. (frp_inst(i,j) .ge. frp_threshold) ) then
                   WRITE(1000+mpiid,*) 'mod_plumerise_after:xlat,xlong,curr_secs,k_min(i,j), k_max(i,j) ',xlat(i,j),xlong(i,j),int(curr_secs),kp1,kp2
                   WRITE(1000+mpiid,*) 'mod_plumerise_after:xlat,xlong,curr_secs,ebu(kts),frp_inst',xlat(i,j),xlong(i,j),int(curr_secs),ebu(i,kts,j),frp_inst(i,j)
                   WRITE(1000+mpiid,*) 'mod_plumerise_after:xlat,xlong,curr_secs,u(10),v(10),w(10),qv(10)',xlat(i,j),xlong(i,j),int(curr_secs),u_in(10),v_in(10),w_in(kte),qv_in(10)
                   WRITE(1000+mpiid,*) 'mod_plumerise_after:xlat,xlong,curr_secs,uspdavg,kpbl_thetav',xlat(i,j),xlong(i,j),int(curr_secs),uspdavg,kpbl_thetav(i,j)
                 IF ( frp_inst(i,j) .ge. 3.e+9 ) then
                   WRITE(1000+mpiid,*) 'mod_plumerise_after:High FRP at : xlat,xlong,curr_secs,frp_inst',xlat(i,j),xlong(i,j),int(curr_secs),frp_inst(i,j)
                 END IF
                 icall = icall + 1
               END IF
!              endif check_frp
!              icall = icall + 1
            enddo
          enddo

        ENDIF check_pl

end subroutine ebu_driver

END module module_plumerise

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
                             kpbl_thetav, kpbl,                      &   ! SRB: added kpbl_thetav and kpbl
                             curr_secs,xlat, xlong , uspdavg2d,      &
                             hpbl2d,  mpiid, alpha,                  &
                             frp_min, frp_wthreshold,zpbl_lim,uspd_lim,   &
                             ids,ide, jds,jde, kds,kde,              &
                             ims,ime, jms,jme, kms,kme,              &
                             its,ite, jts,jte, kts,kte,              & 
                             errmsg, errflg                          ) 

  use rrfs_smoke_config
  !use plume_data_mod
  USE module_zero_plumegen_coms
  USE module_smoke_plumerise
  IMPLICIT NONE

   REAL(kind_phys), intent(in) :: frp_min, frp_wthreshold, zpbl_lim, uspd_lim
    
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme), INTENT(IN ) :: frp_inst         ! RAR: FRP array

   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  xlat,xlong ! SRB

   integer,         dimension(ims:ime, jms:jme),     intent(in)    :: kpbl, kpbl_thetav

   character(*), intent(inout) :: errmsg
   integer, intent(inout) :: errflg
   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   real(kind_phys) :: curr_secs
   INTEGER,      INTENT(IN   ) :: wind_eff_opt
   REAL(kind_phys), INTENT(IN)    :: alpha !  SRB: Enrainment constant for plumerise scheme
   real(kind=kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT ) ::  ebu
   real(kind=kind_phys), INTENT(IN )  :: g, con_cp, con_rd
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme ), INTENT(IN )  :: ebu_in
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme ), INTENT(OUT ) :: flam_frac
   real(kind=kind_phys), DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   z,z_at_w,vvel,u_phy,v_phy,rho_phy,pi_phy,q_vap,theta_phy,wind_phy                     ! RAR, SRB

! Local variables...
      INTEGER :: nv, i, j, k,  kp1, kp2
      INTEGER :: icall
      INTEGER, DIMENSION(ims:ime, jms:jme), INTENT (OUT) :: k_min, k_max      ! Min and max ver. levels for BB injection spread
      REAL, DIMENSION(ims:ime, jms:jme), INTENT (IN) :: uspdavg2d, hpbl2d ! SRB
      real(kind_phys), dimension (kte) :: u_in ,v_in ,w_in ,theta_in ,pi_in, rho_phyin ,qv_in ,zmid, z_lev, uspd ! SRB
      real(kind=kind_phys) :: dz_plume, cpor, con_rocp ! SRB

! MPI variables
      INTEGER, INTENT(IN) :: mpiid

        cpor    =con_cp/con_rd
        con_rocp=con_rd/con_cp

        if ( dbg_opt .and. (mod(int(curr_secs),1800) .eq. 0) ) then
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
                  !uspd(k)= wind_phy(i,k,j) ! SRB
               enddo

             IF (dbg_opt .and. (icall .le. n_dbg_lines) .and. (frp_inst(i,j) .ge. frp_min) ) then
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
                              icall, mpiid, xlat(i,j), xlong(i,j),  & 
                              curr_secs, alpha, frp_min )
               if(errflg/=0) return

               kp1= k_min(i,j)
               kp2= k_max(i,j)   

! SRB: Adding condition for overwriting plumerise levels               
               !uspdavg=SUM(uspd(kts:kpbl(i)))/kpbl(i) !Average wind speed within the boundary layer

! SRB: Adding output
               !uspdavg2(i,j) = uspdavg
               !hpbl_thetav2(i,j) = z_lev(kpbl(i))
               
               IF (frp_inst(i,j) .le. frp_min) THEN
                  !kp1=1
                  !kp2=2
                  flam_frac(i,j)= 0. 
               ELSE IF ( (frp_inst(i,j) .le. frp_wthreshold) .AND. ( uspdavg2d(i,1) .ge. uspd_lim ) .AND. & 
                       ( hpbl2d(i,1) .gt. zpbl_lim) .AND. (wind_eff_opt .eq. 1)) THEN
                  kp1=2
                  kp2=MAX(3,NINT(real(kpbl(i,j))/3._kind_phys))
                  flam_frac(i,j)=0.85 
               ELSE
                  flam_frac(i,j)=0.9  ! kp1,2 come from the plumerise scheme
               END IF
! SRB: End modification 

               ! RAR: emission distribution
               dz_plume= z_at_w(i,kp2,j) - z_at_w(i,kp1,j)
               do k=kp1,kp2-1
                     ebu(i,k,j)=flam_frac(i,j)*ebu_in(i,j)*(z_at_w(i,k+1,j)-z_at_w(i,k,j))/dz_plume
               enddo
               ebu(i,kts,j)= (1.-flam_frac(i,j))* ebu_in(i,j)

               ! For output diagnostic
               k_min(i,j) = kp1
               k_max(i,j) = kp2

               IF ( dbg_opt .and. (icall .le. n_dbg_lines)  .and. (frp_inst(i,j) .ge. frp_min) ) then
                   WRITE(1000+mpiid,*) 'mod_plumerise_after:xlat,xlong,curr_secs,k_min(i,j), k_max(i,j) ',xlat(i,j),xlong(i,j),int(curr_secs),kp1,kp2
                   WRITE(1000+mpiid,*) 'mod_plumerise_after:xlat,xlong,curr_secs,ebu(kts),frp_inst',xlat(i,j),xlong(i,j),int(curr_secs),ebu(i,kts,j),frp_inst(i,j)
                   WRITE(1000+mpiid,*) 'mod_plumerise_after:xlat,xlong,curr_secs,u(10),v(10),w(10),qv(10)',xlat(i,j),xlong(i,j),int(curr_secs),u_in(10),v_in(10),w_in(kte),qv_in(10)
                   !WRITE(1000+mpiid,*) 'mod_plumerise_after:xlat,xlong,curr_secs,uspdavg,kpbl_thetav,kpbl',xlat(i,j),xlong(i,j),int(curr_secs),uspdavg,kpbl_thetav(i,j),kpbl(i)
                 IF ( frp_inst(i,j) .ge. 2.e+10 ) then
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

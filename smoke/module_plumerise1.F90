!>\file  module_plumerise1.F90
!! This file is the fire plume rise driver.

 module module_plumerise1

  use rrfs_smoke_data
  use machine , only : kind_phys
  real(kind=kind_phys),parameter :: p1000mb = 100000.  ! p at 1000mb (pascals)
!- Implementing the fire radiative power (FRP) methodology for biomass burning
!- emissions and convective energy estimation.
!- Saulo Freitas, Gabriel Pereira (INPE/UFJS, Brazil)
!- Ravan Ahmadov, Georg Grell (NOAA, USA)
!- The flag "plumerise_flag" defines the method:
!-    =1 => original method
!-    =2 => FRP based
!-------------------------------------------------------------------------
!
! use module_zero_plumegen_coms
!  integer, parameter :: nveg_agreg      = 4
! integer, parameter :: tropical_forest = 1
! integer, parameter :: boreal_forest   = 2
! integer, parameter :: savannah        = 3

! integer, parameter :: grassland       = 4
!  real(kind_phys), dimension(nveg_agreg) :: firesize,mean_fct
! character(len=20), parameter :: veg_name(nveg_agreg) = (/ &
!                              'Tropical-Forest', &
!                              'Boreal-Forest  ', &
!                              'Savanna        ', &
!                              'Grassland      ' /)
! character(len=20), parameter :: spc_suf(nveg_agreg) = (/ &
!                              'agtf' , &  ! trop forest
!                              'agef' , &  ! extratrop forest
!                              'agsv' , &  ! savanna
!                              'aggr'   /) ! grassland

CONTAINS
subroutine ebu_driver (      data,flam_frac,ebb_smoke,ebu,           &
                             t_phy,q_vap,                            &   ! RAR: moist is replaced with q_vap
                             rho_phy,vvel,u_phy,v_phy,p_phy,         &
                             z_at_w,z,ktau,g,con_cp,con_rd,          &   ! scale_fire_emiss is part of config_flags
                             plume_frp, k_min, k_max,                &   ! RAR:
                             ids,ide, jds,jde, kds,kde,              &
                             ims,ime, jms,jme, kms,kme,              &
                             its,ite, jts,jte, kts,kte, errmsg, errflg)

  use rrfs_smoke_config
  use plume_data_mod
  USE module_zero_plumegen_coms
  USE module_smoke_plumerise
  IMPLICIT NONE
  type(smoke_data), intent(inout) :: data

   REAL(kind_phys), PARAMETER :: frp_threshold= 1.e+7   ! Minimum FRP (Watts) to have plume rise 
    
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme, 2 ), INTENT(IN ) :: plume_frp         ! RAR: FRP etc. array

!   TYPE(grid_config_rec_type),  INTENT(IN )    :: config_flags
   character(*), intent(inout) :: errmsg
   integer, intent(inout) :: errflg
   INTEGER,      INTENT(IN   ) :: ktau,                                    &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
!   real(kind=kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
!         INTENT(IN ) ::                                   moist
   real(kind=kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT ) ::  ebu

   real(kind=kind_phys), INTENT(IN )  :: g, con_cp, con_rd
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme ), INTENT(IN )  :: ebb_smoke
   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme ), INTENT(OUT ) :: flam_frac

!   real(kind=kind_phys), DIMENSION( ims:ime, 1, jms:jme ),                 &
!         INTENT(IN ) ::                                   ebu_in
!   real(kind=kind_phys), DIMENSION( ims:ime, jms:jme ),                 &
!         INTENT(IN ) ::                                                &
!           mean_fct_agtf,mean_fct_agef,&
!           mean_fct_agsv,mean_fct_aggr,firesize_agtf,firesize_agef,       &
!           firesize_agsv,firesize_aggr

   real(kind=kind_phys), DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   t_phy,z,z_at_w,vvel,u_phy,v_phy,rho_phy,p_phy,q_vap                     ! RAR
  ! real(kind=kind_phys), INTENT(IN   ) ::      dtstep

! Local variables...
      INTEGER :: nv, i, j, k,  kp1, kp2
      INTEGER, DIMENSION(ims:ime, jms:jme), INTENT (OUT) :: k_min, k_max      ! Min and max ver. levels for BB injection spread
      !real(kind_phys), dimension (num_ebu) :: eburn_in
      !real(kind_phys), dimension (kte,num_ebu) :: eburn_out
      real(kind_phys), dimension (kte) :: u_in ,v_in ,w_in ,theta_in ,pi_in, rho_phyin ,qv_in ,zmid, z_lev
      real(kind=kind_phys) :: dz_plume, cpor, con_rocp

      !INTEGER, PARAMETER :: kfire_max=30  
! real(kind_phys), dimension(nveg_agreg) :: firesize,mean_fct
!      real(kind_phys) :: sum, ffirs, ratio
!     real(kind_phys),save,dimension(its:ite,jts:jte) :: ffirs
!      nspecies=num_ebu
!     write(0,*)'plumerise'

! RAR:
!      if (config_flags%biomass_burn_opt == BIOMASSB_SMOKE) then
!          do j=jts,jte:
!             do i=its,ite
!                 ebu(i,kts,j,p_ebu_smoke)= ebb_smoke(i,j)
!                ebu(i,kts,j,p_ebu_no)   = ebu_in(i,1,j,p_ebu_in_no)
!                ebu(i,kts,j,p_ebu_co)   = ebu_in(i,1,j,p_ebu_in_co)
!                ebu(i,kts,j,p_ebu_so2)  = ebu_in(i,1,j,p_ebu_in_so2)
!                ebu(i,kts,j,p_ebu_dms)  = ebu_in(i,1,j,p_ebu_in_dms)
!                ebu(i,kts,j,p_ebu_oc)   = ebu_in(i,1,j,p_ebu_in_oc)
!                ebu(i,kts,j,p_ebu_bc)   = ebu_in(i,1,j,p_ebu_in_bc)
!                ebu(i,kts,j,p_ebu_pm25) = ebu_in(i,1,j,p_ebu_in_pm25)
!                ebu(i,kts,j,p_ebu_pm10) = ebu_in(i,1,j,p_ebu_in_pm10)
!            enddo
!          enddo
        cpor    =con_cp/con_rd
        con_rocp=con_rd/con_cp

        IF ( dbg_opt .and. ktau<2000) then
           WRITE(*,*) 'module_plumerise1: its,ite,jts,jte ', its,ite,jts,jte
           WRITE(*,*) 'module_plumerise1: ims,ime,jms,jme ', ims,ime,jms,jme
          !WRITE(*,*) 'module_plumerise1: p_ebu_smoke,num_ebu: ', p_ebu_smoke,num_ebu
           WRITE(*,*) 'module_plumerise1: maxval(ebu(:,kts,:)) ', maxval(ebu(:,kts,:))
         END IF
      !endif

! RAR: setting to zero the ebu emissions at the levels k>1, this is necessary when the plumerise is called, so the emissions at k>1 are updated
       !do nv=1,num_ebu
          do j=jts,jte
            do k=kts+1,kte
               do i=its,ite
                 ebu(i,k,j)=0.
               enddo
            enddo
          enddo
       !enddo

! For now the flammable fraction is constant, based on the namelist. The next
! step to use LU index and meteorology to parameterize it
!    IF (ktau==2) THEN
       do j=jts,jte
        do i=its,ite
           flam_frac(i,j)= 0.
           if (plume_frp(i,j,1) > frp_threshold) then 
              flam_frac(i,j)= 0.9
           end if
        enddo
       enddo
 !   ENDIF


! RAR: new FRP based approach
!check_pl:  IF (config_flags%plumerise_flag == 2 ) THEN    ! if the namelist option is set for plumerise 
! Haiqin: plumerise_flag is added to the namelist options
!check_pl:  IF (do_plumerise) THEN    ! if the namelist option is set for plumerise 
       do j=jts,jte
          do i=its,ite
              ! k_min(i,j)=0
              ! k_max(i,j)=0

!         check_frp: if (.NOT.do_plumerise) then      ! namelist option
!                        ebu(i,kts,j)= ebb_smoke(i,j)  
!                    else                     

               do k=kts,kte
                  u_in(k)=  u_phy(i,k,j)
                  v_in(k)=  v_phy(i,k,j)
                  w_in(k)=  vvel(i,k,j)
                  qv_in(k)= q_vap(i,k,j)    ! RAR: moist(i,k,j,p_qv)
                 !pi_in(k)= cp*(p_phy(i,k,j)/p1000mb)**rcp
                  pi_in(k)= con_cp*(p_phy(i,k,j)/p1000mb)**con_rocp
                  zmid(k)=  z(i,k,j)-z_at_w(i,kts,j)
                  z_lev(k)= z_at_w(i,k,j)-z_at_w(i,kts,j)
                  rho_phyin(k)= rho_phy(i,k,j)
                  theta_in(k)= t_phy(i,k,j)/pi_in(k)*con_cp
                 !theta_in(k)= t_phy(i,k,j)/pi_in(k)*cp
               enddo

             IF (dbg_opt .and. ktau<2000) then
               WRITE(*,*) 'module_plumerise1: i,j ',i,j
               WRITE(*,*) 'module_plumerise1: plume_frp(i,j,:) ',plume_frp(i,j,:)
               WRITE(*,*) 'module_plumerise1: ebu(i,kts,j) ',ebu(i,kts,j)
               WRITE(*,*) 'module_plumerise1: u_in(10),v_in(10),w_in(kte),qv_in(10),pi_in(10) ',u_in(10),v_in(10),w_in(kte),qv_in(10),pi_in(10)
               WRITE(*,*) 'module_plumerise1: zmid(kte),z_lev(kte),rho_phyin(kte),theta_in(kte) ',zmid(kte),z_lev(kte),rho_phyin(kte),theta_in(kte)
               WRITE(*,*) 'module_plumerise1: t_phy(i,kte,j),pi_in(kte)',t_phy(i,kte,j),pi_in(kte)
             END IF

! RAR: the plume rise calculation step:
               CALL plumerise(data,kte,1,1,1,1,1,1,         &
                              !firesize,mean_fct,                    & 
                              !num_ebu, eburn_in, eburn_out,         &
                              u_in, v_in, w_in, theta_in ,pi_in,    &
                              rho_phyin, qv_in, zmid, z_lev,        &
                              plume_frp(i,j,1), k_min(i,j),         & 
                              k_max(i,j), ktau, dbg_opt, g, con_cp, &
                              con_rd, cpor, errmsg, errflg )
                             !k_max(i,j), ktau, config_flags%debug_chem )
               if(errflg/=0) return

               kp1= k_min(i,j)
               kp2= k_max(i,j)   
               dz_plume= z_at_w(i,kp2,j) - z_at_w(i,kp1,j)

                  do k=kp1,kp2-1     
                     ebu(i,k,j)= flam_frac(i,j)* ebb_smoke(i,j)* (z_at_w(i,k+1,j)-z_at_w(i,k,j))/dz_plume
                  enddo
                  ebu(i,kts,j)=   (1.-flam_frac(i,j))* ebb_smoke(i,j)

               IF ( dbg_opt .and. ktau<2000) then
                   WRITE(*,*) 'module_plumerise1: i,j ',i,j
                   WRITE(*,*) 'module_plumerise1: k_min(i,j), k_max(i,j) ',k_min(i,j), k_max(i,j)  
               END IF
!              endif check_frp  
            enddo
          enddo

!        ENDIF check_pl

end subroutine ebu_driver

END module module_plumerise1

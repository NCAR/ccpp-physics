!>\file  module_add_emiss_burn.F90
!! This file adds the biomass burning emissions to the smoke field.

module module_add_emiss_burn
!RAR: significantly modified for the new BB emissions
  use machine ,        only : kind_phys
  use rrfs_smoke_config
CONTAINS
  subroutine add_emis_burn(dtstep,dz8w,rho_phy,pi,ebb_min,          &
                           chem,julday,gmt,xlat,xlong,              &
                           fire_end_hr, peak_hr,time_int,           &
                           coef_bb_dc, fire_hist, hwp, hwp_prevd,   &
                           swdown,ebb_dcycle, ebu_in, ebu,fire_type,&
                           q_vap, add_fire_moist_flux,              &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte,mpiid          )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: julday,                           &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte

   INTEGER, INTENT(IN) :: mpiid
   real(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem   ! shall we set num_chem=1 here?

   real(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         INTENT(INOUT ) ::                                   ebu, q_vap ! SRB: added q_vap

   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: xlat,xlong, swdown
   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp, peak_hr, fire_end_hr, ebu_in !RAR: Shall we make fire_end integer?
   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(INOUT)  :: coef_bb_dc    ! RAR:
   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp_prevd
 
   real(kind_phys), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: dz8w,rho_phy  !,rel_hum
   real(kind_phys), INTENT(IN) ::  dtstep, gmt
   real(kind_phys), INTENT(IN) ::  time_int, pi, ebb_min       ! RAR: time in seconds since start of simulation
   INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: fire_type
   integer, INTENT(IN) ::  ebb_dcycle     ! RAR: this is going to be namelist dependent, ebb_dcycle=means 
   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: fire_hist
!>--local 
   logical, intent(in)  :: add_fire_moist_flux
   integer :: i,j,k,n,m
   integer :: icall=0
   real(kind_phys) :: conv_rho, conv, dm_smoke, dc_hwp, dc_gp, dc_fn !daero_num_wfa, daero_num_ifa !, lu_sum1_5, lu_sum12_14
   
   INTEGER, PARAMETER :: kfire_max=51    ! max vertical level for BB plume rise
   real(kind_phys), PARAMETER :: ef_h2o=324.22  ! Emission factor for water vapor
   ! Constants for the fire diurnal cycle calculation ! JLS - needs to be
   ! defined below due to intent(in) of pi
   real(kind_phys) :: coef_con
   real(kind_phys) :: timeq, fire_age, age_hr, dt1,dt2,dtm         ! For BB emis. diurnal cycle calculation

! For Gaussian diurnal cycle
   real(kind_phys), PARAMETER :: sc_factor=1.  ! to scale up the wildfire emissions, TBD later
   real(kind_phys), PARAMETER :: rinti=2.1813936e-8, ax2=3400., const2=130., &
                   coef2=10.6712963e-4, cx2=7200., timeq_max=3600.*24.
!>-- Fire parameters: Fores west, Forest east, Shrubland, Savannas, Grassland, Cropland
   real(kind_phys), dimension(1:5), parameter :: avg_fire_dur   = (/8.9, 4.2, 3.3, 3.0, 1.4/)
   real(kind_phys), dimension(1:5), parameter :: sigma_fire_dur = (/8.7, 6.0, 5.5, 5.2, 2.4/)

    timeq= gmt*3600._kind_phys + real(time_int,4)
    timeq= mod(timeq,timeq_max)
    coef_con=1._kind_phys/((2._kind_phys*pi)**0.5)

    if (ebb_dcycle==2) then
    
     do j=jts,jte
       do i=its,ite
        fire_age= time_int/3600._kind_phys + (fire_end_hr(i,j)-1._kind_phys)  !One hour delay is due to the latency of the RAVE files
        fire_age= MAX(0.1_kind_phys,fire_age) ! in hours 

          SELECT CASE ( fire_type(i,j) )   !Ag, urban fires, bare land etc.
          CASE (1)
             ! these fires will have exponentially decreasing diurnal cycle,
             coef_bb_dc(i,j) = coef_con*1._kind_phys/(sigma_fire_dur(5) *fire_age) *                          &
                             exp(- ( log(fire_age) - avg_fire_dur(5))**2 /(2._kind_phys*sigma_fire_dur(5)**2 ))

             IF ( dbg_opt .AND. time_int<5000.) then
               WRITE(6,*) 'i,j,peak_hr(i,j) ',i,j,peak_hr(i,j)
               WRITE(6,*) 'coef_bb_dc(i,j) ',coef_bb_dc(i,j)
             END IF

          CASE (2)    ! Savanna and grassland fires
              coef_bb_dc(i,j) = coef_con*1._kind_phys/(sigma_fire_dur(4) *fire_age) *                          &
                              exp(- ( log(fire_age) - avg_fire_dur(4))**2 /(2._kind_phys*sigma_fire_dur(4)**2 ))

              IF ( dbg_opt .AND. time_int<5000.) then
                WRITE(6,*) 'i,j,peak_hr(i,j) ',i,j,peak_hr(i,j)
                WRITE(6,*) 'coef_bb_dc(i,j) ',coef_bb_dc(i,j)
              END IF



          CASE (3)
             !age_hr= fire_age/3600._kind_phys

             IF (swdown(i,j)<.1 .AND. fire_age> 12. .AND. fire_hist(i,j)>0.75) THEN
                 fire_hist(i,j)= 0.75_kind_phys
             ENDIF
             IF (swdown(i,j)<.1 .AND. fire_age> 24. .AND. fire_hist(i,j)>0.5) THEN
                 fire_hist(i,j)= 0.5_kind_phys
             ENDIF
             IF (swdown(i,j)<.1 .AND. fire_age> 48. .AND. fire_hist(i,j)>0.25) THEN
                 fire_hist(i,j)= 0.25_kind_phys
             ENDIF
   
             ! this is based on hwp, hourly or instantenous TBD
             dc_hwp= hwp(i,j)/ MAX(10._kind_phys,hwp_prevd(i,j))
             dc_hwp= MAX(0._kind_phys,dc_hwp)
             dc_hwp= MIN(20._kind_phys,dc_hwp)
   
             ! RAR: Gaussian profile for wildfires
             dt1= abs(timeq - peak_hr(i,j))
             dt2= timeq_max - peak_hr(i,j) + timeq   ! peak hour is always <86400.
             dtm= MIN(dt1,dt2)
             dc_gp = rinti*( ax2 * exp(- dtm**2/(2._kind_phys*cx2**2) ) + const2 - coef2*timeq )
             dc_gp = MAX(0._kind_phys,dc_gp)
   
             !dc_fn = MIN(dc_hwp/dc_gp,3._kind_phys)
             coef_bb_dc(i,j) = fire_hist(i,j)* dc_hwp

             IF ( dbg_opt .AND. time_int<5000.) then
               WRITE(6,*) 'i,j,fire_hist(i,j),peak_hr(i,j) ', i,j,fire_hist(i,j),peak_hr(i,j)
               WRITE(6,*) 'dc_gp,dc_hwp,dc_fn,coef_bb_dc(i,j) ',dc_gp,dc_hwp,dc_fn,coef_bb_dc(i,j)
             END IF

          CASE DEFAULT
          END SELECT
       enddo
      enddo
     endif

     if (mod(int(time_int),1800) .eq. 0) then
        icall = 0
     endif

     do j=jts,jte
      do i=its,ite
       do k=kts,kfire_max
          if (ebu(i,k,j)<ebb_min) cycle

           if (ebb_dcycle==1) then
            conv= dtstep/(rho_phy(i,k,j)* dz8w(i,k,j))
           elseif (ebb_dcycle==2) then
            conv= sc_factor*coef_bb_dc(i,j)*dtstep/(rho_phy(i,k,j)* dz8w(i,k,j))
           endif
           dm_smoke= conv*ebu(i,k,j)
           chem(i,k,j,p_smoke) = chem(i,k,j,p_smoke) + dm_smoke
           chem(i,k,j,p_smoke) = MIN(MAX(chem(i,k,j,p_smoke),0._kind_phys),5.e+3_kind_phys)

           ! SRB: Modifying Water Vapor content based on Emissions
           if (add_fire_moist_flux) then
             q_vap(i,k,j) = q_vap(i,k,j) + (dm_smoke * ef_h2o * 1.e-9)  ! kg/kg:used 1.e-9 as dm_smoke is in ug/kg
             q_vap(i,k,j) = MIN(MAX(q_vap(i,k,j),0._kind_phys),1.e+3_kind_phys)
           endif

           if ( dbg_opt .and. (k==kts .OR. k==kfire_max) .and. (icall .le. n_dbg_lines) ) then
             WRITE(1000+mpiid,*) 'add_emiss_burn:xlat,xlong,curr_secs,fire_type,fire_hist,peak_hr', xlat(i,j),xlong(i,j),int(time_int),fire_type(i,j),fire_hist(i,j),peak_hr(i,j)
             WRITE(1000+mpiid,*) 'add_emiss_burn:xlat,xlong,curr_secs,coef_bb_dc,ebu',xlat(i,j),xlong(i,j),int(time_int),coef_bb_dc(i,j),ebu(i,k,j)
          endif
       enddo
       icall = icall + 1
      enddo
     enddo

    END subroutine add_emis_burn

END module module_add_emiss_burn


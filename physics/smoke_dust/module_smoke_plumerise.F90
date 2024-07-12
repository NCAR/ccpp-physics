!>\file  module_smoke_plumerise.F90
!! This file contains the fire plume rise module.

!-------------------------------------------------------------------------
!- 12 April 2016
!- Implementing the fire radiative power (FRP) methodology for biomass burning
!- emissions and convective energy estimation.
!- Saulo Freitas, Gabriel Pereira (INPE/UFJS, Brazil)
!- Ravan Ahmadov, Georg Grell (NOAA, USA)
!- The flag "plumerise_flag" defines the method:
!-    =1 => original method
!-    =2 => FRP based
!-------------------------------------------------------------------------
module module_smoke_plumerise

  use machine , only : kind_phys
  !use plume_data_mod,  only : num_frp_plume, p_frp_hr, p_frp_std
                              !tropical_forest, boreal_forest, savannah, grassland,   &
                             ! wind_eff
  USE module_zero_plumegen_coms
  USE rrfs_smoke_config, only : n_dbg_lines

  !real(kind=kind_phys),parameter :: rgas=r_d
  !real(kind=kind_phys),parameter :: cpor=cp/r_d
CONTAINS

! RAR:
    subroutine plumerise(m1,m2,m3,ia,iz,ja,jz,                      &
                         up,vp,wp,theta,pp,dn0,rv,zt_rams,zm_rams,  &
                         wind_eff_opt,                              &
                         frp_inst,k1,k2, dbg_opt, g, cp, rgas,      &
                         cpor,  errmsg, errflg, icall, mpiid, lat, long, curr_secs  )

  implicit none

  LOGICAL, INTENT (IN) :: dbg_opt
  INTEGER, INTENT (IN) :: wind_eff_opt, mpiid
  real(kind_phys),  INTENT(IN) ::  lat,long, curr_secs ! SRB

!  INTEGER, PARAMETER ::  ihr_frp=1, istd_frp=2!, imean_fsize=3, istd_fsize=4       ! RAR:

!  integer, intent(in) :: PLUMERISE_flag
  real(kind=kind_phys) :: frp_inst   ! This is the instantenous FRP, at a given time step
  real(kind=kind_phys) :: g, cp, rgas, cpor

  integer :: ng,m1,m2,m3,ia,iz,ja,jz,ibcon,mynum,i,j,k,imm,ixx,ispc !,nspecies


  INTEGER, INTENT (OUT) :: k1,k2
  character(*), intent(inout) :: errmsg
  integer, intent(inout) :: errflg

!  integer :: ncall = 0
  integer :: kmt
!  real(kind=kind_phys),dimension(m1,nspecies), intent(inout) :: eburn_out
!  real(kind=kind_phys),dimension(nspecies), intent(in) :: eburn_in

  real(kind=kind_phys),    dimension(m1,m2,m3) :: up, vp, wp,theta,pp,dn0,rv 
  real(kind=kind_phys),    dimension(m1)       :: zt_rams,zm_rams
  real(kind=kind_phys)                         :: burnt_area,dzi,FRP   ! RAR:
  real(kind=kind_phys),    dimension(2)        :: ztopmax
  real(kind=kind_phys)                         :: q_smold_kgm2
 
  REAL(kind_phys), PARAMETER :: frp_threshold= 1.e+7   ! Minimum FRP (Watts) to have plume rise 

! From plumerise1.F routine
  integer, parameter :: iveg_ag=1 
!  integer, parameter :: tropical_forest = 1
!  integer, parameter :: boreal_forest   = 2
!  integer, parameter :: savannah        = 3
!  integer, parameter :: grassland       = 4
!  real(kind=kind_phys), dimension(nveg_agreg) :: firesize,mean_fct

  INTEGER ::  wind_eff
  INTEGER, INTENT(IN) :: icall
  type(plumegen_coms), pointer :: coms

! Set wind effect from namelist
  wind_eff = wind_eff_opt

!  integer:: iloop
  !REAL(kind=kind_phys), INTENT (IN)   :: convert_smold_to_flam

  !Fator de conversao de unidades
  !!fcu=1.	  !=> kg [gas/part] /kg [ar]
  !!fcu =1.e+12   !=> ng [gas/part] /kg [ar]
  !!real(kind=kind_phys),parameter :: fcu =1.e+6 !=> mg [gas/part] /kg [ar] 
  !----------------------------------------------------------------------
  !               indexacao para o array "plume(k,i,j)" 
  ! k 
  ! 1   => area media (m^2) dos focos  em biomas floresta dentro do gribox i,j 
  ! 2   => area media (m^2) dos focos  em biomas savana dentro do gribox i,j
  ! 3   => area media (m^2) dos focos  em biomas pastagem dentro do gribox i,j
  ! 4   => desvio padrao da area media (m^2) dos focos : floresta
  ! 5   => desvio padrao da area media (m^2) dos focos : savana
  ! 6   => desvio padrao da area media (m^2) dos focos : pastagem
  ! 7 a 9 =>  sem uso
  !10(=k_CO_smold) => parte da emissao total de CO correspondente a fase smoldering
  !11, 12 e 13 => este array guarda a relacao entre
  !               qCO( flaming, floresta) e a quantidade total emitida 
  !               na fase smoldering, isto e;
  !               qCO( flaming, floresta) =  plume(11,i,j)*plume(10,i,j)
  !               qCO( flaming, savana  ) =  plume(12,i,j)*plume(10,i,j)
  !               qCO( flaming, pastagem) =  plume(13,i,j)*plume(10,i,j)
  !20(=k_PM25_smold),21,22 e 23 o mesmo para PM25               
  !
  !24-n1 =>  sem uso
  !----------------------------------------------------------------------
! print *,' Plumerise_scalar 1',ncall
  coms => get_thread_coms()

IF (frp_inst<frp_threshold) THEN
   k1=1
   k2=2
   !return
END IF
    
! print *,' Plumerise_scalar 2',m1
  j=1
  i=1
! do j = ja,jz           ! loop em j
!   do i = ia,iz         ! loop em i
     
     !- if the max value of flaming is close to zero => there is not emission with
     !- plume rise => cycle

  do k = 1,m1 ! loop over vertical grid
    coms%ucon  (k)=up(k,i,j)       ! u wind
    coms%vcon  (k)=vp(k,i,j)       ! v wind
    !coms%wcon  (k)=wp(k,i,j) 	      ! w wind
    coms%thtcon(k)=theta(k,i,j)      ! pot temperature
    coms%picon (k)=pp(k,i,j)                ! exner function
    !coms%tmpcon(k)=coms%thtcon(k)*coms%picon(k)/cp   ! temperature (K)
    !coms%dncon (k)=dn0(k,i,j) 	      ! dry air density (basic state)
    !coms%prcon (k)=(coms%picon(k)/cp)**cpor*p00 ! pressure (Pa)
    coms%rvcon (k)=rv(k,i,j)        ! water vapor mixing ratio
    coms%zcon  (k)=zt_rams(k)               ! termod-point height
    coms%zzcon (k)=zm_rams(k)               ! W-point height
  enddo
 
!  do ispc=2,nspecies
   !  eburn_out(1,ispc) = eburn_in(ispc)    ! eburn_in is the emissions at the 1st level
!     eburn_out(2:m1,ispc)= 0.  ! RAR: k>1 are used from eburn_out
!  enddo

         !- get envinronmental state (temp, water vapor mix ratio, ...)
         call get_env_condition(coms,1,m1,kmt,wind_eff,g,cp,rgas,cpor,errmsg,errflg)
         if(errflg/=0) return

         !- loop over the four types of aggregate biomes with fires for plumerise version 1
         !- for plumerise version 2, there is exist only one loop
        ! iloop=1
!         IF (PLUMERISE_flag == 1) iloop=nveg_agreg

    !lp_veg:  do iveg_ag=1,iloop
    FRP = max(1000.,frp_inst)

    !- loop over the minimum and maximum heat fluxes/FRP
    lp_minmax: do imm=1,2
        if(imm==1 ) then
          burnt_area = 0.7* 0.0006* FRP   !    0.00021* FRP          ! - 0.5*plume_fre(istd_fsize))
        elseif(imm==2 ) then
          burnt_area = 1.3* 0.0006* FRP   ! RAR: Based on Laura's paper I increased the fire size *3. This should depend on the fuel type and meteorology/HWP
        endif
        burnt_area= max(1.0e4,burnt_area)
        
        IF ( dbg_opt .and. (icall .le. n_dbg_lines) .and. (frp_inst .ge. frp_threshold) ) THEN
            WRITE(1000+mpiid,*) 'inside plumerise: xlat,xlong,curr_secs, m1 ', lat,long, int(curr_secs), m1
            WRITE(1000+mpiid,*) 'inside plumerise: xlat,xlong,curr_secs,imm,FRP,burnt_area ', lat, long, int(curr_secs), imm, FRP,burnt_area
        END IF

       !- get fire properties (burned area, plume radius, heating rates ...)
       call get_fire_properties(coms,imm,iveg_ag,burnt_area,FRP,errmsg,errflg)
       if(errflg/=0) return

       !------  generates the plume rise    ------
       call makeplume (coms,kmt,ztopmax(imm),ixx,imm)

       IF ( dbg_opt .and.  (icall .le. n_dbg_lines) .and. (frp_inst .ge. frp_threshold) ) then
            WRITE(1000+mpiid,*) 'inside plumerise after makeplume:xlat,xlong,curr_secs,imm,kmt,ztopmax(imm) ', lat, long, int(curr_secs), imm,kmt, ztopmax(imm)
       END IF

    enddo lp_minmax

        !- define o dominio vertical onde a emissao flaming ira ser colocada
        call set_flam_vert(ztopmax,k1,k2,nkp,coms%zzcon)     !,W_VMD,VMD)

        !     WRITE(6,*) 'module_chem_plumerise_scalar: eburn_out(:,3) ', eburn_out(:,3)

        !- thickness of the vertical layer between k1 and k2 eta levels (lower and upper bounds for the injection height )
        !dzi= 1./(coms%zzcon(k2)-coms%zzcon(k1))   ! RAR: k2>=k1+1  

        !- emission during flaming phase is evenly distributed between levels k1 and k2
        !do k=k1,k2
        !   do ispc= 2,nspecies
        !     eburn_out(k,ispc)= dzi* eburn_in(ispc)
        !   enddo
        !enddo

    !IF (dbg_opt) then
    !    WRITE(*,*) 'plumerise after set_flam_vert: nkp,k1,k2, ', nkp,k1,k2
    !    WRITE(*,*) 'plumerise after set_flam_vert: dzi ', dzi
       !WRITE(*,*) 'plumerise after set_flam_vert: eburn_in(2) ', eburn_in(2)
       !WRITE(*,*) 'plumerise after set_flam_vert: eburn_out(:,2) ',eburn_out(:,2)
    !END IF

!   enddo lp_veg   ! sub-grid vegetation, currently it's aggregated

end subroutine plumerise
!-------------------------------------------------------------------------

subroutine get_env_condition(coms,k1,k2,kmt,wind_eff,g,cp,rgas,cpor,errmsg,errflg)

!se module_zero_plumegen_coms
!use rconstants
implicit none
type(plumegen_coms), pointer :: coms
real(kind=kind_phys) :: g,cp,rgas,cpor
integer :: k1,k2,k,kcon,klcl,kmt,nk,nkmid,i
real(kind=kind_phys),parameter :: p1000mb = 100000.  ! p at 1000mb (pascals)
real(kind=kind_phys),parameter :: p00=p1000mb
real(kind=kind_phys) :: znz,themax,tlll,plll,rlll,zlll,dzdd,dzlll,tlcl,plcl,dzlcl,dummy
!integer :: n_setgrid = 0 
integer :: wind_eff
character(*), intent(inout) :: errmsg
integer, intent(inout) :: errflg

if(.not.coms%initialized) then
 ! n_setgrid = 1
  call set_grid(coms) ! define vertical grid of plume model
                ! coms%zt(k) =  thermo and water levels
                ! coms%zm(k) =  dynamical levels 
endif

znz=coms%zcon(k2)
errflg=1
do k=nkp,1,-1
  if(coms%zt(k).lt.znz) then
    errflg=0
    exit
  endif
enddo
if(errflg/=0) then
  errmsg=' envir stop 12'
  return
endif
!-srf-mb
kmt=min(k,nkp-1)

nk=k2-k1+1
!call htint(nk, coms%wcon,coms%zzcon,kmt,wpe,coms%zt,errmsg,errflg)
!if(errflg/=0) return
 call htint(nk,  coms%ucon,coms%zcon,kmt,coms%upe,coms%zt,errmsg,errflg)
 if(errflg/=0) return
 call htint(nk,  coms%vcon,coms%zcon,kmt,coms%vpe,coms%zt,errmsg,errflg)
 if(errflg/=0) return
 call htint(nk,coms%thtcon,coms%zcon,kmt,coms%the  ,coms%zt,errmsg,errflg)
 if(errflg/=0) return
 call htint(nk, coms%rvcon,coms%zcon,kmt,coms%qvenv,coms%zt,errmsg,errflg)
 if(errflg/=0) return
do k=1,kmt
  coms%qvenv(k)=max(coms%qvenv(k),1e-8)
enddo

coms%pke(1)=coms%picon(1)
do k=1,kmt
  coms%thve(k)=coms%the(k)*(1.+.61*coms%qvenv(k)) ! virtual pot temperature
enddo
do k=2,kmt
  coms%pke(k)=coms%pke(k-1)-g*2.*(coms%zt(k)-coms%zt(k-1))  & ! exner function
        /(coms%thve(k)+coms%thve(k-1))
enddo
do k=1,kmt
  coms%te(k)  = coms%the(k)*coms%pke(k)/cp         ! temperature (K) 
  coms%pe(k)  = (coms%pke(k)/cp)**cpor*p00    ! pressure (Pa)
  coms%dne(k)= coms%pe(k)/(rgas*coms%te(k)*(1.+.61*coms%qvenv(k))) !  dry air density (kg/m3)
!  print*,'ENV=',coms%qvenv(k)*1000., coms%te(k)-273.15,coms%zt(k)
!-srf-mb
  coms%vel_e(k) = sqrt(coms%upe(k)**2+coms%vpe(k)**2)         !-env wind (m/s)
  !print*,'k,coms%vel_e(k),coms%te(k)=',coms%vel_e(k),coms%te(k)
enddo

!-ewe - env wind effect
if(wind_eff < 1)  coms%vel_e(1:kmt) = 0.

!-use este para gerar o RAMS.out
! ------- print environment state
!print*,'k,coms%zt(k),coms%pe(k),coms%te(k)-273.15,coms%qvenv(k)*1000'
!do k=1,kmt
! write(*,100)  k,coms%zt(k),coms%pe(k),coms%te(k)-273.15,coms%qvenv(k)*1000.
! 100 format(1x,I5,4f20.12)
!enddo
!stop 333


!--------- nao eh necessario este calculo
!do k=1,kmt
!  call thetae(coms%pe(k),coms%te(k),coms%qvenv(k),coms%thee(k))
!enddo


!--------- converte press de Pa para kPa para uso modelo de plumerise
do k=1,kmt
 coms%pe(k) = coms%pe(k)*1.e-3
enddo 

return 
end subroutine get_env_condition

!-------------------------------------------------------------------------

subroutine set_grid(coms)
!use module_zero_plumegen_coms  
implicit none
type(plumegen_coms), pointer :: coms
integer :: k,mzp

coms%dz=100. ! set constant grid spacing of plume grid model(meters)

mzp=nkp
coms%zt(1) = coms%zsurf
coms%zm(1) = coms%zsurf
coms%zt(2) = coms%zt(1) + 0.5*coms%dz
coms%zm(2) = coms%zm(1) + coms%dz
do k=3,mzp
 coms%zt(k) = coms%zt(k-1) + coms%dz ! thermo and water levels
 coms%zm(k) = coms%zm(k-1) + coms%dz ! dynamical levels	
enddo
!print*,coms%zsurf
!Print*,coms%zt(:)
do k = 1,mzp-1
   coms%dzm(k) = 1. / (coms%zt(k+1) - coms%zt(k))
enddo 
coms%dzm(mzp)=coms%dzm(mzp-1)

do k = 2,mzp
   coms%dzt(k) = 1. / (coms%zm(k) - coms%zm(k-1))
enddo
coms%dzt(1) = coms%dzt(2) * coms%dzt(2) / coms%dzt(3)

coms%initialized = .true.
   
!   coms%dzm(1) = 0.5/coms%dz
!   coms%dzm(2:mzp) = 1./coms%dz
return
end subroutine set_grid
!-------------------------------------------------------------------------

  SUBROUTINE set_flam_vert(ztopmax,k1,k2,nkp,zzcon) !,W_VMD,VMD)

    REAL(kind=kind_phys)    , INTENT(IN)  :: ztopmax(2)
    INTEGER , INTENT(OUT) :: k1,k2

    ! plumegen_coms
    INTEGER , INTENT(IN)  :: nkp
    REAL(kind=kind_phys)    , INTENT(IN)  :: zzcon(nkp)

    INTEGER imm,k
    INTEGER, DIMENSION(2)  :: k_lim

    !- version 2
!    REAL(kind=kind_phys)    , INTENT(IN)  :: W_VMD(nkp,2)
!    REAL(kind=kind_phys)    , INTENT(OUT) ::   VMD(nkp,2)
!    real(kind=kind_phys)   w_thresold,xxx
!    integer k_initial,k_final,ko,kk4,kl

    !- version 1
    DO imm=1,2
       ! checar 
       !    do k=1,m1-1
       DO k=1,nkp-1
          IF(zzcon(k) > ztopmax(imm)) EXIT
       ENDDO
       k_lim(imm) = k
    ENDDO
    k1= MIN(MAX(4,k_lim(1)),51)
    k2= MIN(51,k_lim(2))   ! RAR: the model doesn't simulate very high injection heights, so it's safe to assume maximum heigh of 12km AGL for HRRR grid

    IF (k2 <= k1) THEN
       !print*,'1: ztopmax k=',ztopmax(1), k1
       !print*,'2: ztopmax k=',ztopmax(2), k2
       k2= k1+1     ! RAR: I added k1+1
    ENDIF
    
    !- version 2    
    !- vertical mass distribution
    !- 
!    w_thresold = 1.
!    DO imm=1,2

!       VMD(1:nkp,imm)= 0.
!       xxx=0.
!       k_initial= 0
!       k_final  = 0
    
       !- define range of the upper detrainemnt layer
!       do ko=nkp-10,2,-1
     
!        if(w_vmd(ko,imm) < w_thresold) cycle
     
!        if(k_final==0) k_final=ko
     
!        if(w_vmd(ko,imm)-1. > w_vmd(ko-1,imm)) then
!          k_initial=ko
!          exit
!        endif
      
!       enddo
       !- if there is a non zero depth layer, make the mass vertical distribution 
!       if(k_final > 0 .and. k_initial > 0) then 
       
!           k_initial=int((k_final+k_initial)*0.5)
       
           !- parabolic vertical distribution between k_initial and k_final
!           kk4 = k_final-k_initial+2
!           do ko=1,kk4-1
!               kl=ko+k_initial-1
!               VMD(kl,imm) = 6.* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
!           enddo
!           if(sum(VMD(1:NKP,imm)) .ne. 1.) then
!             xxx= ( 1.- sum(VMD(1:NKP,imm)) )/float(k_final-k_initial+1)
!             do ko=k_initial,k_final
!                VMD(ko,imm) = VMD(ko,imm)+ xxx !- values between 0 and 1.
!              enddo
               ! print*,'new mass=',sum(mass)*100.,xxx
               !pause
!           endif
!        endif !k_final > 0 .and. k_initial > 

!    ENDDO
    
  END SUBROUTINE set_flam_vert
!-------------------------------------------------------------------------

subroutine get_fire_properties(coms,imm,iveg_ag,burnt_area,FRP,errmsg,errflg)
!use module_zero_plumegen_coms
implicit none
type(plumegen_coms), pointer :: coms
integer ::  moist,  i,  icount,imm,iveg_ag  !,plumerise_flag
real(kind=kind_phys)::   bfract,  effload,  heat,  hinc ,burnt_area,heat_fluxW,FRP
!real(kind=kind_phys),    dimension(2,4) :: heat_flux
integer, intent(inout) :: errflg
character(*), intent(inout) :: errmsg
INTEGER, parameter :: use_last = 1    ! RAR 10/31/2022: I set to one, checking with Saulo

!real(kind=kind_phys), parameter :: beta = 5.0   !ref.: Wooster et al., 2005
REAL(kind=kind_phys), parameter :: beta = 0.88  !ref.: Paugam et al., 2015

!
coms%area = burnt_area! area of burn, m^2

!IF ( PLUMERISE_flag == 1) THEN
!    !fluxo de calor para o bioma
!    heat_fluxW = heat_flux(imm,iveg_ag) * 1000. ! converte para W/m^2

!ELSEIF ( PLUMERISE_flag == 2) THEN
    ! "beta" factor converts FRP to convective energy
    heat_fluxW = beta*(FRP/coms%area)/0.55 ! in W/m^2
! FIXME: These five lines were not in the known-working version. Delete them?
!    if(coms%area<1e-6) then
!      heat_fluxW = 0
!    else
!      heat_fluxW = beta*(FRP/coms%area)/0.55 ! in W/m^2
!    endif

!ENDIF

coms%mdur = 53        ! duration of burn, minutes
coms%bload = 10.      ! total loading, kg/m**2 
moist = 10       ! fuel moisture, %. average fuel moisture,percent dry
coms%maxtime =coms%mdur+2  ! model time, min
!heat = 21.e6    !- joules per kg of fuel consumed                   
!heat = 15.5e6   !joules/kg - cerrado
heat = 19.3e6    !joules/kg - floresta em alta floresta (mt)
!coms%alpha = 0.1      !- entrainment constant
coms%alpha = 0.05      !- entrainment constant

!-------------------- printout ----------------------------------------

!!WRITE ( * ,  * ) ' SURFACE =', COMS%ZSURF, 'M', '  LCL =', COMS%ZBASE, 'M'  
!
!PRINT*,'======================================================='
!print * , ' FIRE BOUNDARY CONDITION   :'  
!print * , ' DURATION OF BURN, MINUTES =',COMS%MDUR  
!print * , ' AREA OF BURN, HA	      =',COMS%AREA*1.e-4
!print * , ' HEAT FLUX, kW/m^2	      =',heat_fluxW*1.e-3
!print * , ' TOTAL LOADING, KG/M**2    =',COMS%BLOAD  
!print * , ' FUEL MOISTURE, %	      =',MOIST !average fuel moisture,percent dry
!print * , ' MODEL TIME, MIN.	      =',COMS%MAXTIME  
!
!
!
! ******************** fix up inputs *********************************
!
                                             
!IF (MOD (COMS%MAXTIME, 2) .NE.0) COMS%MAXTIME = COMS%MAXTIME+1  !make coms%maxtime even
                                                  
COMS%MAXTIME = COMS%MAXTIME * 60  ! and put in seconds
!
COMS%RSURF = SQRT (COMS%AREA / 3.14159) !- entrainment surface radius (m)

COMS%FMOIST   = MOIST / 100.       !- fuel moisture fraction
!
!
! calculate the energy flux and water content at lboundary.
! fills heating() on a minute basis. could ask for a file at this po
! in the program. whatever is input has to be adjusted to a one
! minute timescale.
!
                        
  DO I = 1, ntime         !- make sure of energy release
    COMS%HEATING (I) = 0.0001  !- avoid possible divide by 0
  enddo  
!                                  
  COMS%TDUR = COMS%MDUR * 60.       !- number of seconds in the burn

  bfract = 1.             !- combustion factor

  EFFLOAD = COMS%BLOAD * BFRACT  !- patchy burning
  
!     spread the burning evenly over the interval
!     except for the first few minutes for stability
  ICOUNT = 1  
!
  if(COMS%MDUR > NTIME) then
    errmsg = 'Increase time duration (ntime) in min - see file "module_zero_plumegen_coms.F90"'
    errflg = 1
    return
  endif

  DO WHILE (ICOUNT.LE.COMS%MDUR)                             
!  COMS%HEATING (ICOUNT) = HEAT * EFFLOAD / COMS%TDUR  ! W/m**2 
!  COMS%HEATING (ICOUNT) = 80000.  * 0.55         ! W/m**2 

   COMS%HEATING (ICOUNT) = heat_fluxW  * 0.55     ! W/m**2 (0.55 converte para energia convectiva)
   ICOUNT = ICOUNT + 1  
  ENDDO  
!     ramp for 5 minutes, RAR: in the current version this is inactive
 IF(use_last /= 1) THEN

    HINC = COMS%HEATING (1) / 4.  
    COMS%HEATING (1) = 0.1  
    COMS%HEATING (2) = HINC  
    COMS%HEATING (3) = 2. * HINC  
    COMS%HEATING (4) = 3. * HINC  
 ELSE
    HINC = COMS%HEATING (1) / 4.   ! RAR: this needs to be revised later
    IF(imm==1) THEN
       !HINC = COMS%HEATING (1) / 4.
       COMS%HEATING (1) = 0.1  
       COMS%HEATING (2) = HINC  
       COMS%HEATING (3) = 2. * HINC  
       COMS%HEATING (4) = 3. * HINC 
    ELSE 
       COMS%HEATING (2) = COMS%HEATING (1)+ HINC  
       COMS%HEATING (3) = COMS%HEATING (2)+ HINC  
       COMS%HEATING (4) = COMS%HEATING (3)+ HINC 
    ENDIF
 ENDIF

return
end subroutine get_fire_properties
!-------------------------------------------------------------------------------
!
SUBROUTINE MAKEPLUME (coms,kmt,ztopmax,ixx,imm)  
!
! *********************************************************************
!
!    EQUATION SOURCE--Kessler Met.Monograph No. 32 V.10 (K)
!    Alan Weinstein, JAS V.27 pp 246-255. (W),
!    Ogura and Takahashi, Monthly Weather Review V.99,pp895-911 (OT)
!    Roger Pielke,Mesoscale Meteorological Modeling,Academic Press,1984
!    Originally developed by: Don Latham (USFS)
!
!
! ************************ VARIABLE ID ********************************
!
!     DT=COMPUTING TIME INCREMENT (SEC)
!     DZ=VERTICAL INCREMENT (M)
!     LBASE=LEVEL ,CLOUD BASE
!
!     CONSTANTS:
!       G = GRAVITATIONAL ACCELERATION 9.80796 (M/SEC/SEC).
!       R = DRY AIR GAS CONSTANT (287.04E6 JOULE/KG/DEG K)
!       CP = SPECIFIC HT. (1004 JOULE/KG/DEG K)
!       HEATCOND = HEAT OF CONDENSATION (2.5E6 JOULE/KG)
!       HEATFUS = HEAT OF FUSION (3.336E5 JOULE/KG)
!       HEATSUBL = HEAT OF SUBLIMATION (2.83396E6 JOULE/KG)
!       EPS = RATIO OF MOL.WT. OF WATER VAPOR TO THAT OF DRY AIR (0.622)
!       DES = DIFFERENCE BETWEEN VAPOR PRESSURE OVER WATER AND ICE (MB)
!       TFREEZE = FREEZING TEMPERATURE (K)
!
!
!     PARCEL VALUES:
!       T = TEMPERATURE (K)
!       TXS = TEMPERATURE EXCESS (K)
!       QH = HYDROMETEOR WATER CONTENT (G/G DRY AIR)
!       QHI = HYDROMETEOR ICE CONTENT (G/G DRY AIR)
!       QC = WATER CONTENT (G/G DRY AIR)
!       QVAP = WATER VAPOR MIXING RATIO (G/G DRY AIR)
!       QSAT = SATURATION MIXING RATIO (G/G DRY AIR)
!       RHO = DRY AIR DENSITY (G/M**3) MASSES = RHO*Q'S IN G/M**3
!       ES = SATURATION VAPOR PRESSURE (kPa)
!
!     ENVIRONMENT VALUES:
!       TE = TEMPERATURE (K)
!       PE = PRESSURE (kPa)
!       QVENV = WATER VAPOR (G/G)
!       RHE = RELATIVE HUMIDITY FRACTION (e/esat)
!       DNE = dry air density (kg/m^3)
!
!     HEAT VALUES:
!       HEATING = HEAT OUTPUT OF FIRE (WATTS/M**2)
!       MDUR = DURATION OF BURN, MINUTES
!
!       W = VERTICAL VELOCITY (M/S)
!       RADIUS=ENTRAINMENT RADIUS (FCN OF Z)
!	RSURF = ENTRAINMENT RADIUS AT GROUND (SIMPLE PLUME, TURNER)
!	ALPHA = ENTRAINMENT CONSTANT
!       MAXTIME = TERMINATION TIME (MIN)
!
!
!**********************************************************************
!**********************************************************************               
!use module_zero_plumegen_coms 
implicit none 
!logical :: endspace  
type(plumegen_coms), pointer :: coms
character (len=10) :: varn
integer ::  izprint, iconv,  itime, k, kk, kkmax, deltak,ilastprint,kmt &
           ,ixx,nrectotal,i_micro,n_sub_step
real(kind=kind_phys) ::  vc, g,  r,  cp,  eps,  &
         tmelt,  heatsubl,  heatfus,  heatcond, tfreeze, &
         ztopmax, wmax, rmaxtime, es, esat, heat,dt_save !ESAT_PR,
character (len=2) :: cixx 
! Set threshold to be the same as dz=100., the constant grid spacing of plume grid model(meters) found in set_grid()
    REAL(kind=kind_phys) :: DELZ_THRESOLD = 100. 

    INTEGER     :: imm

!  real(kind=kind_phys), external:: esat_pr!
!
! ******************* SOME CONSTANTS **********************************
!
!      XNO=10.0E06 median volume diameter raindrop (K table 4)
!      VC = 38.3/(XNO**.125) mean volume fallspeed eqn. (K)
!
parameter (vc = 5.107387)  
parameter (g = 9.80796, r = 287.04, cp = 1004., eps = 0.622,  tmelt = 273.3)
parameter (heatsubl = 2.834e6, heatfus = 3.34e5, heatcond = 2.501e6)
parameter (tfreeze = 269.3)  
!
coms%tstpf = 2.0  !- timestep factor
coms%viscosity = 500.!- coms%viscosity constant (original value: 0.001)

nrectotal=150
!
!*************** PROBLEM SETUP AND INITIAL CONDITIONS *****************
coms%mintime = 1  
ztopmax = 0. 
coms%ztop    = 0. 
   coms%time = 0.  
     coms%dt = 1.
   wmax = 1. 
kkmax   = 10
deltaK  = 20
ilastprint=0
COMS%L       = 1   ! COMS%L initialization

!--- initialization
CALL INITIAL(coms,kmt)  

!--- initial print fields:
izprint  = 0          ! if = 0 => no printout
!if (izprint.ne.0) then
! write(cixx(1:2),'(i2.2)') ixx
! open(2, file = 'debug.'//cixx//'.dat')  
! open(19,file='plumegen9.'//cixx//'.gra',         &
!     form='unformatted',access='direct',status='unknown',  &
!     recl=4*nrectotal)  !PC   
!     recl=1*nrectotal) !sx6 e tupay
! call printout (izprint,nrectotal)
! ilastprint=2
!endif     

! ******************* model evolution ******************************
rmaxtime = float(coms%maxtime)
!
!print * ,' TIME=',coms%time,' RMAXTIME=',rmaxtime
!print*,'======================================================='
 DO WHILE (COMS%TIME.LE.RMAXTIME)  !beginning of time loop

!   do itime=1,120

!-- set model top integration
    coms%nm1 = min(kmt, kkmax + deltak)
!sam 81  format('nm1=',I0,' from kmt=',I0,' kkmax=',I0,' deltak=',I0)
!sam     write(0,81) coms%nm1,kmt,kkmax,deltak
!-- set timestep
    !coms%dt = (coms%zm(2)-coms%zm(1)) / (coms%tstpf * wmax)  
    coms%dt = min(5.,(coms%zm(2)-coms%zm(1)) / (coms%tstpf * wmax))
                                
!-- elapsed time, sec
    coms%time = coms%time+coms%dt 
!-- elapsed time, minutes                                      
    coms%mintime = 1 + int (coms%time) / 60     
    wmax = 1.  !no zeroes allowed.
!************************** BEGIN SPACE LOOP **************************

!-- zerout all model tendencies
    call tend0_plumerise(coms)

!-- bounday conditions (k=1)
    COMS%L=1
    call lbound(coms)

!-- dynamics for the level k>1 
!-- W advection 
!   call vel_advectc_plumerise(COMS%NM1,COMS%WC,COMS%WT,COMS%DNE,COMS%DZM)
    call vel_advectc_plumerise(COMS%NM1,COMS%WC,COMS%WT,COMS%RHO,COMS%DZM)
  
!-- scalars advection 1
    call scl_advectc_plumerise(coms,'SC',COMS%NM1)

!-- scalars advection 2
    !call scl_advectc_plumerise2(coms,'SC',COMS%NM1)

!-- scalars entrainment, adiabatic
    call scl_misc(coms,COMS%NM1)
    
!-- scalars dinamic entrainment
     call  scl_dyn_entrain(COMS%NM1,nkp,coms%wbar,coms%w,coms%adiabat,coms%alpha,coms%radius,coms%tt,coms%t,coms%te,coms%qvt,coms%qv,coms%qvenv,coms%qct,coms%qc,coms%qht,coms%qh,coms%qit,coms%qi,&
                    coms%vel_e,coms%vel_p,coms%vel_t,coms%rad_p,coms%rad_t)

!-- gravity wave damping using Rayleigh friction layer fot COMS%T
    call damp_grav_wave(1,coms%nm1,deltak,coms%dt,coms%zt,coms%zm,coms%w,coms%t,coms%tt,coms%qv,coms%qh,coms%qi,coms%qc,coms%te,coms%pe,coms%qvenv)

!-- microphysics
!   goto 101 ! bypass microphysics
    dt_save=coms%dt
    n_sub_step=3
    coms%dt=coms%dt/float(n_sub_step)

    do i_micro=1,n_sub_step
!-- sedim ?
     call fallpart(coms,COMS%NM1)
!-- microphysics
     coms%L=2
     do while(coms%L<=coms%nm1-1)
     !do L=2,coms%nm1-1
        COMS%WBAR    = 0.5*(coms%W(COMS%L)+coms%W(COMS%L-1))
        ES      = ESAT_PR (COMS%T(COMS%L))            !BLOB SATURATION VAPOR PRESSURE, EM KPA
        COMS%QSAT(COMS%L) = (EPS * ES) / (COMS%PE(COMS%L) - ES)  !BLOB SATURATION LWC G/G DRY AIR
        COMS%EST (COMS%L) = ES  
!sam         if(.not.coms%pe(coms%L)>0 .or. .not. coms%T(coms%L)>200) then
!sam 1304      format('(1304) bad input to rho at L=',I0,' with pe=',F12.5,' T=',F12.5)
!sam           write(0,1304) coms%L,coms%PE(coms%L),coms%T(coms%L)
!sam         endif
        COMS%RHO (COMS%L) = 3483.8 * COMS%PE (COMS%L) / COMS%T (COMS%L) ! AIR PARCEL DENSITY , G/M**3
!srf18jun2005
!	IF (COMS%W(COMS%L) .ge. 0.) COMS%DQSDZ = (COMS%QSAT(COMS%L  ) - COMS%QSAT(COMS%L-1)) / (COMS%ZT(COMS%L  ) -COMS%ZT(COMS%L-1))
!	IF (COMS%W(COMS%L) .lt. 0.) COMS%DQSDZ = (COMS%QSAT(COMS%L+1) - COMS%QSAT(COMS%L  )) / (COMS%ZT(COMS%L+1) -COMS%ZT(COMS%L  ))
        IF (COMS%W(COMS%L) .ge. 0.) then 
           COMS%DQSDZ = (COMS%QSAT(COMS%L+1) - COMS%QSAT(COMS%L-1)) / (COMS%ZT(COMS%L+1 )-COMS%ZT(COMS%L-1))
        ELSE
           COMS%DQSDZ = (COMS%QSAT(COMS%L+1) - COMS%QSAT(COMS%L-1)) / (COMS%ZT(COMS%L+1) -COMS%ZT(COMS%L-1))
        ENDIF 

        call waterbal(coms)
        coms%L=coms%L+1
     enddo
    enddo
    coms%dt=dt_save
!
    101 continue
!
!-- W-viscosity for stability 
    call visc_W(coms,coms%nm1,deltak,kmt)

!-- update scalars
    call update_plumerise(coms,coms%nm1,'S')
    
    call hadvance_plumerise(1,coms%nm1,coms%dt,COMS%WC,COMS%WT,COMS%W,coms%mintime) 

!-- Buoyancy
    call buoyancy_plumerise(COMS%NM1, COMS%T, COMS%TE, COMS%QV, COMS%QVENV, COMS%QH, COMS%QI, COMS%QC, COMS%WT, COMS%SCR1)
 
!-- Entrainment 
    call entrainment(coms,COMS%NM1,COMS%W,COMS%WT,COMS%RADIUS,COMS%ALPHA)

!-- update W
    call update_plumerise(coms,coms%nm1,'W')

    call hadvance_plumerise(2,coms%nm1,coms%dt,COMS%WC,COMS%WT,COMS%W,coms%mintime) 


!-- misc
    do k=2,coms%nm1
!    coms%pe esta em kpa  - esat do rams esta em mbar = 100 Pa = 0.1 kpa
!    es       = 0.1*esat (coms%t(k)) !blob saturation vapor pressure, em kPa
!    rotina do plumegen calcula em kPa
     es       = esat_pr (coms%t(k))  !blob saturation vapor pressure, em kPa
     coms%qsat(k) = (eps * es) / (coms%pe(k) - es)  !blob saturation lwc g/g dry air
     coms%est (k) = es  
     coms%txs (k) = coms%t(k) - coms%te(k)
!sam         if(.not.coms%pe(K)>0 .or. .not. coms%T(K)>200) then
!sam 1305      format('(1305) bad input to rho at K=',I0,' with pe=',F12.5,' T=',F12.5)
!sam           write(0,1305) K,coms%PE(K),coms%T(K)
!sam         endif
     coms%rho (k) = 3483.8 * coms%pe (k) / coms%t (k) ! air parcel density , g/m**3
                                       ! no pressure diff with radius
     if((abs(coms%wc(k))).gt.wmax) wmax = abs(coms%wc(k)) ! keep wmax largest w
    enddo  

! Gravity wave damping using Rayleigh friction layer for W
    call damp_grav_wave(2,coms%nm1,deltak,coms%dt,coms%zt,coms%zm,coms%w,coms%t,coms%tt,coms%qv,coms%qh,coms%qi,coms%qc,coms%te,coms%pe,coms%qvenv)
!---
       !- update radius
       do k=2,coms%nm1
        coms%radius(k) = coms%rad_p(k)
       enddo
      !-- try to find the plume top (above surface height)
       kk = 1
       DO WHILE (coms%w (kk) .GT. 1.)  
          kk = kk + 1  
          coms%ztop =  coms%zm(kk) 
          !print*,'W=',coms%w (kk)
       ENDDO
       !
       coms%ztop_(coms%mintime) = coms%ztop
       ztopmax = MAX (coms%ztop, ztopmax) 
       kkmax   = MAX (kk  , kkmax  ) 
       !print * ,'ztopmax=', coms%mintime,'mn ',coms%ztop_(coms%mintime), ztopmax

       !
       ! if the solution is going to a stationary phase, exit
       IF(coms%mintime > 10) THEN                 
          !   if(coms%mintime > 20) then                     
          !    if( abs(coms%ztop_(coms%mintime)-coms%ztop_(coms%mintime-10)) < COMS%DZ ) exit   
          IF( ABS(coms%ztop_(coms%mintime)-coms%ztop_(coms%mintime-10)) < DELZ_THRESOLD) then 
            !- determine W parameter to determine the VMD
            !do k=2,coms%nm1
            !   W_VMD(k,imm) = coms%w(k)
            !enddo
          EXIT ! finish the integration
         ENDIF  
       ENDIF

   ! if(ilastprint == coms%mintime) then
   !   call printout (izprint,nrectotal)  
   !   ilastprint = coms%mintime+1
   ! endif      
                               

ENDDO   !do next timestep

!print * ,' ztopmax=',ztopmax,'m',coms%mintime,'mn '
!print*,'======================================================='
!
!the last printout
!if (izprint.ne.0) then
! call printout (izprint,nrectotal)  
! close (2)            
! close (19)            
!endif

RETURN  
END SUBROUTINE MAKEPLUME
!-------------------------------------------------------------------------------
!
SUBROUTINE BURN(COMS, EFLUX, WATER)  
!	
!- calculates the energy flux and water content at lboundary
!use module_zero_plumegen_coms                               
implicit none
type(plumegen_coms), pointer :: coms
!real(kind=kind_phys), parameter :: HEAT = 21.E6 !Joules/kg
!real(kind=kind_phys), parameter :: HEAT = 15.5E6 !Joules/kg - cerrado
real(kind=kind_phys), parameter :: HEAT = 19.3E6 !Joules/kg - floresta em Alta Floresta (MT)
real(kind=kind_phys) :: eflux,water
!
! The emission factor for water is 0.5. The water produced, in kg,
! is then  fuel mass*0.5 + (moist/100)*mass per square meter.
! The fire burns for DT out of TDUR seconds, the total amount of
! fuel burned is AREA*COMS%BLOAD*(COMS%DT/TDUR) kg. this amount of fuel is
! considered to be spread over area AREA and so the mass burned per
! unit area is COMS%BLOAD*(COMS%DT/TDUR), and the rate is COMS%BLOAD/TDUR.
!        
IF (COMS%TIME.GT.COMS%TDUR) THEN !is the burn over?   
   EFLUX = 0.000001    !prevent a potential divide by zero
   WATER = 0.  
   RETURN  
ELSE  
!                                                   
   EFLUX = COMS%HEATING (COMS%MINTIME)                          ! Watts/m**2                                                   
!  WATER = EFLUX * (COMS%DT / HEAT) * (0.5 + COMS%FMOIST)       ! kg/m**2 
   WATER = EFLUX * (COMS%DT / HEAT) * (0.5 + COMS%FMOIST) /0.55 ! kg/m**2 
   WATER = WATER * 1000.                              ! g/m**2
!
!        print*,'BURN:',coms%time,EFLUX/1.e+9
ENDIF  
!
RETURN  
END SUBROUTINE BURN
!-------------------------------------------------------------------------------
!
SUBROUTINE LBOUND (coms)  
!
! ********** BOUNDARY CONDITIONS AT ZSURF FOR PLUME AND CLOUD ********
!
! source of equations: J.S. Turner Buoyancy Effects in Fluids
!                      Cambridge U.P. 1973 p.172,
!                      G.A. Briggs Plume Rise, USAtomic Energy Commissio
!                      TID-25075, 1969, P.28
!
! fundamentally a point source below ground. at surface, this produces
! a velocity w(1) and temperature T(1) which vary with time. There is
! also a water load which will first saturate, then remainder go into
! QC(1).
! EFLUX = energy flux at ground,watt/m**2 for the last DT
!
!use module_zero_plumegen_coms  
implicit none
type(plumegen_coms), pointer :: coms
real(kind=kind_phys), parameter :: g = 9.80796, r = 287.04, cp = 1004.6, eps = 0.622,tmelt = 273.3
real(kind=kind_phys), parameter :: tfreeze = 269.3, pi = 3.14159, e1 = 1./3., e2 = 5./3.
real(kind=kind_phys) :: es,  esat, eflux, water,  pres, c1,  c2, f, zv,  denscor, xwater !,ESAT_PR
!  real(kind=kind_phys), external:: esat_pr!

!            
COMS%QH (1) = COMS%QH (2)   !soak up hydrometeors
COMS%QI (1) = COMS%QI (2)              
COMS%QC (1) = 0.       !no cloud here
!
!
   CALL BURN (COMS, EFLUX, WATER)  
!
!  calculate parameters at boundary from a virtual buoyancy point source
!
   PRES = COMS%PE (1) * 1000.   !need pressure in N/m**2
                              
   C1 = 5. / (6. * COMS%ALPHA)  !alpha is entrainment constant

   C2 = 0.9 * COMS%ALPHA  

   F = EFLUX / (PRES * CP * PI)  
                             
   F = G * R * F * COMS%AREA  !buoyancy flux
                 
   ZV = C1 * COMS%RSURF  !virtual boundary height
                                   
   COMS%W (1) = C1 * ( (C2 * F) **E1) / ZV**E1  !boundary velocity
                                         
   DENSCOR = C1 * F / G / (C2 * F) **E1 / ZV**E2   !density correction

   COMS%T (1) = COMS%TE (1) / (1. - DENSCOR)    !temperature of virtual plume at zsurf
   
!
   COMS%WC(1) = COMS%W(1)
    COMS%VEL_P(1) = 0.
    coms%rad_p(1) = coms%rsurf

   !COMS%SC(1) = COMS%SCE(1)+F/1000.*coms%dt  ! gas/particle (g/g)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     match dw/dz,dt/dz at the boundary. F is conserved.
!
   !COMS%WBAR = COMS%W (1) * (1. - 1. / (6. * ZV) )  
   !COMS%ADVW = COMS%WBAR * COMS%W (1) / (3. * ZV)  
   !COMS%ADVT = COMS%WBAR * (5. / (3. * ZV) ) * (DENSCOR / (1. - DENSCOR) )  
   !COMS%ADVC = 0.  
   !COMS%ADVH = 0.  
   !COMS%ADVI = 0.  
   !COMS%ADIABAT = - COMS%WBAR * G / CP  
   COMS%VTH (1) = - 4.  
   COMS%VTI (1) = - 3.  
   COMS%TXS (1) = COMS%T (1) - COMS%TE (1)  

   COMS%VISC (1) = COMS%VISCOSITY  

!sam         if(.not.coms%pe(1)>0 .or. .not. coms%T(1)>200) then
!sam 1306      format('(1306) bad input to rho at 1=',I0,' with pe=',F12.5,' T=',F12.5)
!sam           write(0,1306) 1,coms%PE(1),coms%T(1)
!sam         endif
   COMS%RHO (1) = 3483.8 * COMS%PE (1) / COMS%T (1)   !air density at level 1, g/m**3

   XWATER = WATER / max(1e-20, COMS%W (1) * COMS%DT * COMS%RHO (1) )   !firewater mixing ratio
                                            
   COMS%QV (1) = XWATER + COMS%QVENV (1)  !plus what's already there 


!  COMS%PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
!  ES       = 0.1*ESAT (COMS%T(1)) !blob saturation vapor pressure, em kPa
!  rotina do plumegen ja calcula em kPa
   ES       = ESAT_PR (COMS%T(1))  !blob saturation vapor pressure, em kPa

   COMS%EST  (1)  = ES                                  
   COMS%QSAT (1) = (EPS * ES) / max(1e-20, COMS%PE (1) - ES)   !blob saturation lwc g/g dry air
  
   IF (COMS%QV (1) .gt. COMS%QSAT (1) ) THEN  
       COMS%QC (1) = COMS%QV   (1) - COMS%QSAT (1) + COMS%QC (1)  !remainder goes into cloud drops
       COMS%QV (1) = COMS%QSAT (1)  
   ENDIF  
!
   CALL WATERBAL  (COMS)
!
RETURN  
END SUBROUTINE LBOUND
!-------------------------------------------------------------------------------
!
SUBROUTINE INITIAL (coms,kmt)  
!
! ************* SETS UP INITIAL CONDITIONS FOR THE PROBLEM ************
!use module_zero_plumegen_coms 
implicit none 
type(plumegen_coms), pointer :: coms
real(kind=kind_phys), parameter :: tfreeze = 269.3
integer ::  isub,  k,  n1,  n2,  n3,  lbuoy,  itmp,  isubm1 ,kmt
real(kind=kind_phys) ::     xn1,  xi,  es,  esat!,ESAT_PR
!
COMS%N=kmt
! initialize temperature structure,to the end of equal spaced sounding,
  do k = 1, COMS%N 
  COMS%TXS (k) = 0.0  
    COMS%W (k) = 0.0             
    COMS%T (k) = COMS%TE(k)       !blob set to environment		  
    COMS%WC(k) = 0.0
    COMS%WT(k) = 0.0
    COMS%QV(k) = COMS%QVENV (k)   !blob set to environment             
   COMS%VTH(k) = 0.		!initial rain velocity = 0	                     
   COMS%VTI(k) = 0.		!initial ice  velocity = 0	                     
    COMS%QH(k) = 0.		!no rain			     
    COMS%QI(k) = 0.		!no ice 			     
    COMS%QC(k) = 0.		!no cloud drops	                     
!  COMS%PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
!  ES       = 0.1*ESAT (COMS%T(k)) !blob saturation vapor pressure, em kPa
!  rotina do plumegen calcula em kPa
   ES       = ESAT_PR (COMS%T(k))  !blob saturation vapor pressure, em kPa
   COMS%EST  (k) = ES  
   COMS%QSAT (k) = (.622 * ES) / (COMS%PE (k) - ES) !saturation lwc g/g
!sam         if(.not.coms%pe(k)>0 .or. .not. coms%T(k)>200) then
!sam 1307      format('(1307) bad input to rho at k=',I0,' with pe=',F12.5,' T=',F12.5)
!sam           write(0,1307) k,coms%PE(k),coms%T(k)
!sam         endif
   COMS%RHO  (k) = 3483.8 * COMS%PE (k) / COMS%T (k) 	!dry air density g/m**3    
       COMS%VEL_P(k) = 0.
       coms%rad_p(k) = 0.
  enddo  

! Initialize the entrainment radius, Turner-style plume
  coms%radius(1) = coms%rsurf
  do k=2,COMS%N
     coms%radius(k) = coms%radius(k-1)+(6./5.)*coms%alpha*(coms%zt(k)-coms%zt(k-1))
  enddo
! Initialize the entrainment radius, Turner-style plume
    coms%radius(1) = coms%rsurf
    coms%rad_p(1)  = coms%rsurf
    DO k=2,COMS%N
       coms%radius(k) = coms%radius(k-1)+(6./5.)*coms%alpha*(coms%zt(k)-coms%zt(k-1))
       coms%rad_p(k)  = coms%radius(k)
   ENDDO
    
!  Initialize the viscosity
   COMS%VISC (1) = COMS%VISCOSITY
   do k=2,COMS%N
     !COMS%VISC (k) = COMS%VISCOSITY!max(1.e-3,coms%visc(k-1) - 1.* COMS%VISCOSITY/float(nkp))
     COMS%VISC (k) = max(1.e-3,coms%visc(k-1) - 1.* COMS%VISCOSITY/float(nkp))
   enddo
!--   Initialize gas/concentration
  !DO k =10,20
  !   COMS%SC(k) = 20.
  !ENDDO
  !stop 333

   CALL LBOUND(COMS)

RETURN  
END SUBROUTINE INITIAL
!-------------------------------------------------------------------------------
!
subroutine damp_grav_wave(ifrom,nm1,deltak,dt,zt,zm,w,t,tt,qv,qh,qi,qc,te,pe,qvenv)
implicit none
integer nm1,ifrom,deltak
real(kind=kind_phys) dt
real(kind=kind_phys), dimension(nm1) :: w,t,tt,qv,qh,qi,qc,te,pe,qvenv,dummy,zt,zm

if(ifrom==1) then
 call friction(ifrom,nm1,deltak,dt,zt,zm,t,tt    ,te)
!call friction(ifrom,nm1,dt,zt,zm,qv,coms%qvt,qvenv)
 return
endif 

dummy(:) = 0.
if(ifrom==2) call friction(ifrom,nm1,deltak,dt,zt,zm,w,dummy ,dummy)
!call friction(ifrom,nm1,dt,zt,zm,qi,coms%qit ,dummy)
!call friction(ifrom,nm1,dt,zt,zm,qh,coms%qht ,dummy)
!call friction(ifrom,nm1,dt,zt,zm,qc,coms%qct ,dummy)
return
end subroutine damp_grav_wave
!-------------------------------------------------------------------------------
!
subroutine friction(ifrom,nm1,deltak,dt,zt,zm,var1,vart,var2)
implicit none
real(kind=kind_phys), dimension(nm1) :: var1,var2,vart,zt,zm
integer k,nfpt,kf,nm1,ifrom,deltak
real(kind=kind_phys) zmkf,ztop,distim,c1,c2,dt

!nfpt=50
!kf = nm1 - nfpt
!kf = nm1 - int(deltak/2)
 kf = nm1 - int(deltak)

zmkf = zm(kf) !old: float(kf )*coms%dz
ztop = zm(nm1)
!distim = min(4.*dt,200.)
!distim = 60.
 distim = min(3.*dt,60.)

c1 = 1. / (distim * (ztop - zmkf))
c2 = dt * c1

if(ifrom == 1) then  
  do k = nm1,2,-1
   if (zt(k) .le. zmkf) cycle
   vart(k) = vart(k)   + c1 * (zt(k) - zmkf)*(var2(k) - var1(k))
  enddo
elseif(ifrom == 2) then
  do k = nm1,2,-1
   if (zt(k) .le. zmkf) cycle
   var1(k) =  var1(k) + c2 * (zt(k) - zmkf)*(var2(k) - var1(k))
  enddo
endif
return
end subroutine friction
!-------------------------------------------------------------------------------
!
subroutine vel_advectc_plumerise(m1,wc,wt,rho,dzm)

implicit none
integer :: k,m1
real(kind=kind_phys), dimension(m1) :: wc,wt,flxw,dzm,rho
real(kind=kind_phys), dimension(m1) :: dn0 ! var local
real(kind=kind_phys) :: c1z

!dzm(:)= 1./coms%dz

dn0(1:m1)=rho(1:m1)*1.e-3 ! converte de cgs para mks

flxw(1) = wc(1) * dn0(1) 

do k = 2,m1-1
   flxw(k) = wc(k) * .5 * (dn0(k) + dn0(k+1))
enddo

! Compute advection contribution to W tendency

c1z = .5 

do k = 2,m1-2

   wt(k) = wt(k)  &
      + c1z * dzm(k) / (dn0(k) + dn0(k+1)) *     (   &
	(flxw(k) + flxw(k-1))  * (wc(k) + wc(k-1))   &
      - (flxw(k) + flxw(k+1))  * (wc(k) + wc(k+1))   &
      + (flxw(k+1) - flxw(k-1)) * 2.* wc(k)       )

enddo

return
end subroutine vel_advectc_plumerise
!-------------------------------------------------------------------------------
!
subroutine hadvance_plumerise(iac,m1,dt,wc,wt,wp,mintime)

implicit none
integer :: k,iac
integer :: m1,mintime
real(kind=kind_phys), dimension(m1) :: dummy, wc,wt,wp
real(kind=kind_phys) eps,dt
!     It is here that the Asselin filter is applied.  For the velocities
!     and pressure, this must be done in two stages, the first when
!     IAC=1 and the second when IAC=2.


eps = .2
if(mintime == 1) eps=0.5

!     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.
!
call predict_plumerise(m1,wc,wp,wt,dummy,iac,2.*dt,eps)
!print*,'mintime',mintime,eps
!do k=1,m1
!   print*,'W-HAD',k,wc(k),wp(k),wt(k)
!enddo
return
end subroutine hadvance_plumerise
!-------------------------------------------------------------------------------
!
subroutine predict_plumerise(npts,ac,ap,fa,af,iac,dtlp,epsu)
implicit none
integer :: npts,iac,m
real(kind=kind_phys) :: epsu,dtlp
real(kind=kind_phys), dimension(*) :: ac,ap,fa,af

!     For IAC=3, this routine moves the arrays AC and AP forward by
!     1 time level by adding in the prescribed tendency. It also
!     applies the Asselin filter given by:

!              {AC} = AC + EPS * (AP - 2 * AC + AF)

!     where AP,AC,AF are the past, current and future time levels of A.
!     All IAC=1 does is to perform the {AC} calculation without the AF
!     term present.  IAC=2 completes the calculation of {AC} by adding
!     the AF term only, and advances AC by filling it with input AP
!     values which were already updated in ACOUSTC.
!

if (iac .eq. 1) then
   do m = 1,npts
      ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))
   enddo
   return
elseif (iac .eq. 2) then
   do m = 1,npts
      af(m) = ap(m)
      ap(m) = ac(m) + epsu * af(m)
   enddo
!elseif (iac .eq. 3) then
!   do m = 1,npts
!      af(m) = ap(m) + dtlp * fa(m)
!   enddo
!   if (ngrid .eq. 1 .and. ipara .eq. 0) call cyclic(nzp,nxp,nyp,af,'T')
!   do m = 1,npts
!      ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
!   enddo
endif

do m = 1,npts
  ac(m) = af(m)
enddo
return
end subroutine predict_plumerise
!-------------------------------------------------------------------------------
!
subroutine  buoyancy_plumerise(m1, T, TE, QV, QVENV, QH, QI, QC, WT, scr1)
implicit none
integer :: k,m1
real(kind=kind_phys), parameter :: g = 9.8, eps = 0.622, gama = 0.5 ! mass virtual coeff.
real(kind=kind_phys), dimension(m1) :: T, TE, QV, QVENV, QH, QI, QC, WT, scr1
real(kind=kind_phys) :: TV,TVE,QWTOTL,umgamai
real(kind=kind_phys), parameter :: mu = 0.15 

!- orig
umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
                       ! das pertubacoes nao-hidrostaticas no campo de pressao

!- new                 ! Siesbema et al, 2004
!umgamai = 1./(1.-2.*mu)

do k = 2,m1-1

    TV =   T(k) * (1. + (QV(k)   /EPS))/(1. + QV(k)   )  !blob virtual temp.                                        	   
    TVE = TE(k) * (1. + (QVENV(k)/EPS))/(1. + QVENV(k))  !and environment

    QWTOTL = QH(k) + QI(k) + QC(k)                       ! QWTOTL*G is drag
!- orig
   !scr1(k)= G*( umgamai*(  TV - TVE) / TVE   - QWTOTL) 
    scr1(k)= G*  umgamai*( (TV - TVE) / TVE   - QWTOTL) 

    !if(k .lt. 10)print*,'BT',k,TV,TVE,TVE,QWTOTL
enddo

do k = 2,m1-2
    wt(k) = wt(k)+0.5*(scr1(k)+scr1(k+1))
!   print*,'W-BUO',k,wt(k),scr1(k),scr1(k+1)
enddo

end subroutine  buoyancy_plumerise
!-------------------------------------------------------------------------------
!
subroutine ENTRAINMENT(coms,m1,w,wt,radius,ALPHA)
implicit none
type(plumegen_coms), pointer :: coms
integer :: k,m1
real(kind=kind_phys), dimension(m1) :: w,wt,radius
REAL(kind=kind_phys) DMDTM,WBAR,RADIUS_BAR,umgamai,DYN_ENTR,ALPHA
real(kind=kind_phys), parameter :: mu = 0.15 ,gama = 0.5 ! mass virtual coeff.

!- new - Siesbema et al, 2004
!umgamai = 1./(1.-2.*mu)

!- orig
!umgamai = 1
umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
                       ! das pertubacoes nao-hidrostaticas no campo de pressao

!
!-- ALPHA/RADIUS(COMS%L) = (1/M)DM/COMS%DZ  (W 14a)
  do k=2,m1-1

!-- for W: WBAR is only W(k)
!     WBAR=0.5*(W(k)+W(k-1))           
      WBAR=W(k)          
      RADIUS_BAR = 0.5*(RADIUS(k) + RADIUS(k-1))
! orig
     !DMDTM =           2. * ALPHA * ABS (WBAR) / RADIUS_BAR  != (1/M)DM/COMS%DT
      DMDTM = umgamai * 2. * ALPHA * ABS (WBAR) / RADIUS_BAR  != (1/M)DM/COMS%DT

!--  DMDTM*W(COMS%L) entrainment,
      wt(k) = wt(k)  - DMDTM*ABS (WBAR)
      !print*,'W-ENTR=',k,w(k),- DMDTM*ABS (WBAR)
      
      !if(COMS%VEL_P (k) - COMS%VEL_E (k) > 0.) cycle

       !-   dynamic entrainment
       DYN_ENTR =  (2./3.1416)*0.5*ABS (COMS%VEL_P(k)-COMS%VEL_E(k)+COMS%VEL_P(k-1)-COMS%VEL_E(k-1)) /RADIUS_BAR

       wt(k) = wt(k)  - DYN_ENTR*ABS (WBAR)
       
       !- entraiment acceleration for output only
       !dwdt_entr(k) =  - DMDTM*ABS (WBAR)- DYN_ENTR*ABS (WBAR)
  enddo
end subroutine  ENTRAINMENT
!-------------------------------------------------------------------------------
!
subroutine scl_advectc_plumerise(coms,varn,mzp)
!use module_zero_plumegen_coms
implicit none
type(plumegen_coms), pointer :: coms
integer :: mzp
character(len=*) :: varn
real(kind=kind_phys) :: dtlto2
integer :: k

!  wp => w
!- Advect  scalars
   dtlto2   = .5 * coms%dt
!  coms%vt3dc(1) =      (coms%w(1) + coms%wc(1)) * dtlto2 * coms%dne(1)
   coms%vt3dc(1) =      (coms%w(1) + coms%wc(1)) * dtlto2 * coms%rho(1)*1.e-3!converte de CGS p/ MKS
   coms%vt3df(1) = .5 * (coms%w(1) + coms%wc(1)) * dtlto2 * coms%dzm(1)

   do k = 2,mzp
!     coms%vt3dc(k) =  (coms%w(k) + coms%wc(k)) * dtlto2 *.5 * (coms%dne(k) + coms%dne(k+1))
      coms%vt3dc(k) =  (coms%w(k) + coms%wc(k)) * dtlto2 *.5 * (coms%rho(k) + coms%rho(k+1))*1.e-3
      coms%vt3df(k) =  (coms%w(k) + coms%wc(k)) * dtlto2 *.5 *  coms%dzm(k)
     !print*,'coms%vt3df-coms%vt3dc',k,coms%vt3dc(k),coms%vt3df(k)
   enddo

 
!-srf-24082005
!  do k = 1,mzp-1
  do k = 1,mzp
     coms%vctr1(k) = (coms%zt(k+1) - coms%zm(k)) * coms%dzm(k)
     coms%vctr2(k) = (coms%zm(k)   - coms%zt(k)) * coms%dzm(k)
!    coms%vt3dk(k) = coms%dzt(k) / coms%dne(k)
     coms%vt3dk(k) = coms%dzt(k) /(coms%rho(k)*1.e-3)
     !print*,'Coms%Vt3dk',k,coms%dzt(k) , coms%dne(k)
  enddo

!      scalarp => scalar_tab(coms%n,ngrid)%var_p
!      scalart => scalar_tab(coms%n,ngrid)%var_t

!- temp advection tendency (COMS%TT)
   coms%scr1=COMS%T
   call fa_zc_plumerise(mzp                   &
             	       ,COMS%T	  ,coms%scr1  (1)  &
             	       ,coms%vt3dc (1) ,coms%vt3df (1)  &
             	       ,coms%vt3dg (1) ,coms%vt3dk (1)  &
             	       ,coms%vctr1,coms%vctr2	      )

   call advtndc_plumerise(mzp,COMS%T,coms%scr1(1),COMS%TT,coms%dt)

!- water vapor advection tendency (COMS%QVT)
   coms%scr1=COMS%QV
   call fa_zc_plumerise(mzp                  &
             	       ,COMS%QV	  ,coms%scr1  (1)  &
             	       ,coms%vt3dc (1) ,coms%vt3df (1)  &
             	       ,coms%vt3dg (1) ,coms%vt3dk (1)  &
             	       ,coms%vctr1,coms%vctr2	     )

   call advtndc_plumerise(mzp,COMS%QV,coms%scr1(1),COMS%QVT,coms%dt)

!- liquid advection tendency (COMS%QCT)
   coms%scr1=COMS%QC
   call fa_zc_plumerise(mzp                  &
             	       ,COMS%QC	  ,coms%scr1  (1)  &
             	       ,coms%vt3dc (1) ,coms%vt3df (1)  &
             	       ,coms%vt3dg (1) ,coms%vt3dk (1)  &
             	       ,coms%vctr1,coms%vctr2	     )

   call advtndc_plumerise(mzp,COMS%QC,coms%scr1(1),COMS%QCT,coms%dt)

!- ice advection tendency (COMS%QIT)
   coms%scr1=COMS%QI
   call fa_zc_plumerise(mzp                  &
             	       ,COMS%QI	  ,coms%scr1  (1)  &
             	       ,coms%vt3dc (1) ,coms%vt3df (1)  &
             	       ,coms%vt3dg (1) ,coms%vt3dk (1)  &
             	       ,coms%vctr1,coms%vctr2	     )

   call advtndc_plumerise(mzp,COMS%QI,coms%scr1(1),COMS%QIT,coms%dt)

!- hail/rain advection tendency (COMS%QHT)
!   if(ak1 > 0. .or. ak2 > 0.) then

      coms%scr1=COMS%QH
      call fa_zc_plumerise(mzp                  &
             	          ,COMS%QH	    ,coms%scr1  (1)  &
             	          ,coms%vt3dc (1) ,coms%vt3df (1)  &
             	          ,coms%vt3dg (1) ,coms%vt3dk (1)  &
             	          ,coms%vctr1,coms%vctr2	       )

      call advtndc_plumerise(mzp,COMS%QH,coms%scr1(1),COMS%QHT,coms%dt)
!   endif
    !- horizontal wind advection tendency (COMS%VEL_T)
    coms%scr1=COMS%VEL_P
    call fa_zc_plumerise(mzp		      &
    			,COMS%VEL_P     ,coms%scr1  (1)  &
    			,coms%vt3dc (1) ,coms%vt3df (1)  &
    			,coms%vt3dg (1) ,coms%vt3dk (1)  &
    			,coms%vctr1,coms%vctr2	     )

    call advtndc_plumerise(mzp,COMS%VEL_P,coms%scr1(1),COMS%VEL_T,coms%dt)

    !- vertical radius transport

    coms%scr1=coms%rad_p
    call fa_zc_plumerise(mzp                  &
             	        ,coms%rad_p     ,coms%scr1  (1)  &
             	        ,coms%vt3dc (1) ,coms%vt3df (1)  &
             	        ,coms%vt3dg (1) ,coms%vt3dk (1)  &
             	        ,coms%vctr1,coms%vctr2	     )

    call advtndc_plumerise(mzp,coms%rad_p,coms%scr1(1),coms%rad_t,coms%dt)


   return
!
!- gas/particle advection tendency (COMS%SCT)
!    if(varn == 'SC')return
   coms%scr1=COMS%SC
   call fa_zc_plumerise(mzp		    &
   	     	       ,COMS%SC	 ,coms%scr1  (1)  &
   	     	       ,coms%vt3dc (1) ,coms%vt3df (1)  &
   	     	       ,coms%vt3dg (1) ,coms%vt3dk (1)  &
   	     	       ,coms%vctr1,coms%vctr2	     )
   
   call advtndc_plumerise(mzp,COMS%SC,coms%scr1(1),COMS%SCT,coms%dt)


return
end subroutine scl_advectc_plumerise
!-------------------------------------------------------------------------------
!
subroutine fa_zc_plumerise(m1,scp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2)

implicit none
integer :: m1,k
real(kind=kind_phys) :: dfact
real(kind=kind_phys), dimension(m1) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dk
real(kind=kind_phys), dimension(m1) :: vctr1,vctr2

dfact = .5

! Compute scalar flux VT3DG
      do k = 1,m1-1
         vt3dg(k) = vt3dc(k)                   &
                  * (vctr1(k) * scr1(k)        &
                  +  vctr2(k) * scr1(k+1)      &
                  +  vt3df(k) * (scr1(k) - scr1(k+1)))
      enddo
      
! Modify fluxes to retain positive-definiteness on scalar quantities.
!    If a flux will remove 1/2 quantity during a timestep,
!    reduce to first order flux. This will remain positive-definite
!    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!    both fluxes are evacuating the box.

do k = 1,m1-1
 if (vt3dc(k) .gt. 0.) then
   if (vt3dg(k) * vt3dk(k)    .gt. dfact * scr1(k)) then
	 vt3dg(k) = vt3dc(k) * scr1(k)
   endif
 elseif (vt3dc(k) .lt. 0.) then
   if (-vt3dg(k) * vt3dk(k+1) .gt. dfact * scr1(k+1)) then
	 vt3dg(k) = vt3dc(k) * scr1(k+1)
   endif
 endif

enddo

! Compute flux divergence
do k = 2,m1-1
    scr1(k) = scr1(k)  &
            + vt3dk(k) * ( vt3dg(k-1) - vt3dg(k) &
            + scp  (k) * ( vt3dc(k)   - vt3dc(k-1)))
enddo
return
end subroutine fa_zc_plumerise
!-------------------------------------------------------------------------------
!
subroutine advtndc_plumerise(m1,scp,sca,sct,dtl)
implicit none
integer :: m1,k
real(kind=kind_phys) :: dtl,dtli
real(kind=kind_phys), dimension(m1) :: scp,sca,sct

dtli = 1. / dtl
do k = 2,m1-1
   sct(k) = sct(k) + (sca(k)-scp(k)) * dtli
enddo
return
end subroutine advtndc_plumerise
!-------------------------------------------------------------------------------
!
subroutine tend0_plumerise(coms)
implicit none
type(plumegen_coms), pointer :: coms
 coms%wt(1:coms%nm1)  = 0.
 coms%tt(1:coms%nm1)  = 0.
coms%qvt(1:coms%nm1)  = 0.
coms%qct(1:coms%nm1)  = 0.
coms%qht(1:coms%nm1)  = 0.
coms%qit(1:coms%nm1)  = 0.
coms%vel_t(1:coms%nm1)  = 0.
coms%rad_t(1:coms%nm1)  = 0.
!coms%sct(1:coms%nm1)  = 0.
end subroutine tend0_plumerise

!     ****************************************************************

subroutine scl_misc(coms,m1)
!use module_zero_plumegen_coms
implicit none
type(plumegen_coms), pointer :: coms
real(kind=kind_phys), parameter :: g = 9.81, cp=1004.
integer m1,k
real(kind=kind_phys) dmdtm

 do k=2,m1-1
      COMS%WBAR    = 0.5*(COMS%W(k)+COMS%W(k-1))  
!-- dry adiabat
      COMS%ADIABAT = - COMS%WBAR * G / CP 
!      
!-- entrainment     
      DMDTM = 2. * COMS%ALPHA * ABS (COMS%WBAR) / COMS%RADIUS (k)  != (1/M)DM/COMS%DT
      
!-- tendency temperature = adv + adiab + entrainment
      COMS%TT(k) = COMS%TT(K) + COMS%ADIABAT - DMDTM * ( COMS%T  (k) -    COMS%TE (k) ) 

!-- tendency water vapor = adv  + entrainment
      COMS%QVT(K) = COMS%QVT(K)         - DMDTM * ( COMS%QV (k) - COMS%QVENV (k) )

      COMS%QCT(K) = COMS%QCT(K)	      - DMDTM * ( COMS%QC (k)  )
      COMS%QHT(K) = COMS%QHT(K)	      - DMDTM * ( COMS%QH (k)  )
      COMS%QIT(K) = COMS%QIT(K)	      - DMDTM * ( COMS%QI (k)  )

      !-- tendency horizontal speed = adv  + entrainment
      COMS%VEL_T(K) = COMS%VEL_T(K)     - DMDTM * ( COMS%VEL_P (k) - COMS%VEL_E (k) )

      !-- tendency horizontal speed = adv  + entrainment
      coms%rad_t(K) = coms%rad_t(K)     + 0.5*DMDTM*(6./5.)*COMS%RADIUS (k)
!-- tendency gas/particle = adv  + entrainment
!      COMS%SCT(K) = COMS%SCT(K)         - DMDTM * ( COMS%SC (k) -   COMS%SCE (k) )

enddo
end subroutine scl_misc
!     ****************************************************************

  SUBROUTINE scl_dyn_entrain(m1,nkp,wbar,w,adiabat,alpha,radius,tt,t,te,qvt,qv,qvenv,qct,qc,qht,qh,qit,qi,&
                    vel_e,vel_p,vel_t,rad_p,rad_t)
    implicit none

    INTEGER , INTENT(IN)    :: m1

    ! plumegen_coms
    INTEGER , INTENT(IN)    :: nkp
    REAL(kind=kind_phys)    , INTENT(INOUT) :: wbar 
    REAL(kind=kind_phys)    , INTENT(IN)    :: w(nkp)
    REAL(kind=kind_phys)    , INTENT(INOUT) :: adiabat 
    REAL(kind=kind_phys)    , INTENT(IN)    :: alpha
    REAL(kind=kind_phys)    , INTENT(IN)    :: radius(nkp)
    REAL(kind=kind_phys)    , INTENT(INOUT) :: tt(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: t(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: te(nkp)
    REAL(kind=kind_phys)    , INTENT(INOUT) :: qvt(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: qv(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: qvenv(nkp)
    REAL(kind=kind_phys)    , INTENT(INOUT) :: qct(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: qc(nkp)
    REAL(kind=kind_phys)    , INTENT(INOUT) :: qht(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: qh(nkp)
    REAL(kind=kind_phys)    , INTENT(INOUT) :: qit(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: qi(nkp)

    REAL(kind=kind_phys)    , INTENT(IN)    :: vel_e(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: vel_p(nkp)
    REAL(kind=kind_phys)    , INTENT(INOUT) :: vel_t(nkp)
    REAL(kind=kind_phys)    , INTENT(INOUT) :: rad_T(nkp)
    REAL(kind=kind_phys)    , INTENT(IN)    :: rad_p(nkp)

    real(kind=kind_phys), parameter :: g = 9.81, cp=1004., pi=3.1416

    integer k
    real(kind=kind_phys) dmdtm

    DO k=2,m1-1
      !      
      !-- tendency horizontal radius from dyn entrainment
     	   !rad_t(K) = rad_t(K)   +	(vel_e(k)-vel_p(k)) /pi
     	    rad_t(K) = rad_t(K)   + ABS((vel_e(k)-vel_p(k)))/pi
      
      !-- entrainment	  
     	   !DMDTM = (2./3.1416)  *     (VEL_E (k) - VEL_P (k)) / RADIUS (k)  
     	    DMDTM = (2./3.1416)  *  ABS(VEL_E (k) - VEL_P (k)) / RADIUS (k)  
      
      !-- tendency horizontal speed  from dyn entrainment
     	    VEL_T(K) = VEL_T(K)     - DMDTM * ( VEL_P (k) - VEL_E (k) )
      
      !     if(VEL_P (k) - VEL_E (k) > 0.) cycle
      
      !-- tendency temperature  from dyn entrainment
     	    TT(k) = TT(K)	    - DMDTM * ( T (k) - TE  (k) ) 
      
      !-- tendency water vapor  from dyn entrainment
   	    QVT(K) = QVT(K)	    - DMDTM * ( QV (k) - QVENV (k) )
      
     	    QCT(K) = QCT(K)	    - DMDTM * ( QC (k)  )
     	    QHT(K) = QHT(K)	    - DMDTM * ( QH (k)  )
     	    QIT(K) = QIT(K)	    - DMDTM * ( QI (k)  )
      
      !-- tendency gas/particle  from dyn entrainment
      !	 COMS%SCT(K) = COMS%SCT(K)	 - DMDTM * ( SC (k) - COMS%SCE (k) )
    
    ENDDO
   END SUBROUTINE scl_dyn_entrain

!     ****************************************************************

subroutine visc_W(coms,m1,deltak,kmt)
!use module_zero_plumegen_coms
implicit none
type(plumegen_coms), pointer :: coms
integer m1,k,deltak,kmt,m2
real(kind=kind_phys) dz1t,dz1m,dz2t,dz2m,d2wdz,d2tdz  ,d2qvdz ,d2qhdz ,d2qcdz ,d2qidz ,d2scdz, &
 d2vel_pdz,d2rad_dz
!sam real(kind=kind_phys) :: old_tt
logical, save, volatile :: printed = .false.


!srf--- 17/08/2005
!m2=min(m1+deltak,kmt)
m2=min(m1,kmt)

!do k=2,m1-1
do k=2,m2-1
 DZ1T   = 0.5*(COMS%ZT(K+1)-COMS%ZT(K-1))
 DZ2T   = COMS%VISC (k) / (DZ1T * DZ1T)  
 DZ1M   = 0.5*(COMS%ZM(K+1)-COMS%ZM(K-1))
 DZ2M   = COMS%VISC (k) / (DZ1M * DZ1M)  
 D2WDZ  = (COMS%W  (k + 1) - 2 * COMS%W  (k) + COMS%W  (k - 1) ) * DZ2M  
 D2TDZ  = (COMS%T  (k + 1) - 2 * COMS%T  (k) + COMS%T  (k - 1) ) * DZ2T  
 D2QVDZ = (COMS%QV (k + 1) - 2 * COMS%QV (k) + COMS%QV (k - 1) ) * DZ2T  
 D2QHDZ = (COMS%QH (k + 1) - 2 * COMS%QH (k) + COMS%QH (k - 1) ) * DZ2T 
 D2QCDZ = (COMS%QC (k + 1) - 2 * COMS%QC (k) + COMS%QC (k - 1) ) * DZ2T  
 D2QIDZ = (COMS%QI (k + 1) - 2 * COMS%QI (k) + COMS%QI (k - 1) ) * DZ2T  
 !D2SCDZ = (COMS%SC (k + 1) - 2 * COMS%SC (k) + COMS%SC (k - 1) ) * DZ2T 
 d2vel_pdz=(coms%vel_p  (k + 1) - 2 * coms%vel_p  (k) + coms%vel_p  (k - 1) ) * DZ2T
 d2rad_dz =(coms%rad_p  (k + 1) - 2 * coms%rad_p  (k) + coms%rad_p  (k - 1) ) * DZ2T
 
  COMS%WT(k) =   COMS%WT(k) + D2WDZ 
!sam   old_tt=coms%tt(k)
  COMS%TT(k) =   COMS%TT(k) + D2TDZ                          
!sam   if(.not. coms%tt(k)>-10 .and. .not. printed) then
!sam 1924 format("(1924) visc_W Bad TT at k=",I0," TT=",F12.5," old_TT=",F12.5," d2tdz=",F12.5," visc=",F12.5)
!sam 1925 format("(1925)   T = ",F12.5,",",F12.5,",",F12.5," ZT=",F12.5,",",F12.5)
!sam      write(0,1924) k, COMS%TT(k), old_TT, d2tdz, coms%visc(k)
!sam      write(0,1925) coms%T(k-1),coms%T(k),coms%T(k+1),coms%ZT(k-1),coms%ZT(k+1)
!sam      printed = .true.
!sam   endif
 COMS%QVT(k) =  COMS%QVT(k) + D2QVDZ 
 COMS%QCT(k) =  COMS%QCT(k) + D2QCDZ 
 COMS%QHT(k) =  COMS%QHT(k) + D2QHDZ 
 COMS%QIT(k) =  COMS%QIT(k) + D2QIDZ     
 coms%vel_t(k) =   coms%vel_t(k) + d2vel_pdz
 coms%rad_t(k) =   coms%rad_t(k) + d2rad_dz
 !COMS%SCT(k) =  COMS%SCT(k) + D2SCDZ
 !print*,'W-COMS%VISC=',k,D2WDZ
enddo  

end subroutine visc_W

!     ****************************************************************

subroutine update_plumerise(coms,m1,varn)
!use module_zero_plumegen_coms
implicit none
type(plumegen_coms), pointer :: coms
integer m1,k
character(len=*) :: varn
!sam real(kind_phys) :: old_t
 
if(varn == 'W') then

 do k=2,m1-1
   COMS%W(k) =  COMS%W(k) +  COMS%WT(k) * COMS%DT  
 enddo
 return

else 
do k=2,m1-1
!sam   old_t = coms%t(k)
   COMS%T(k) =  COMS%T(k) +  COMS%TT(k) * COMS%DT  
!sam    if(.not. coms%t(k)>200) then
!sam 1921 format("(1921) update_plumerise Bad T at k=",I0," T=",F12.5," old_T=",F12.5," TT=",F12.5," DT=",F12.5)
!sam      write(0,1921) k, COMS%T(k), old_T, coms%tt(k), coms%dt
!sam    endif

  COMS%QV(k) = COMS%QV(k) + COMS%QVT(k) * COMS%DT  

  COMS%QC(k) = COMS%QC(k) + COMS%QCT(k) * COMS%DT !cloud drops travel with air 
  COMS%QH(k) = COMS%QH(k) + COMS%QHT(k) * COMS%DT  
  COMS%QI(k) = COMS%QI(k) + COMS%QIT(k) * COMS%DT 
! COMS%SC(k) = COMS%SC(k) + COMS%SCT(k) * COMS%DT 

!srf---18jun2005  
  COMS%QV(k) = max(0., COMS%QV(k))
  COMS%QC(k) = max(0., COMS%QC(k))
  COMS%QH(k) = max(0., COMS%QH(k))
  COMS%QI(k) = max(0., COMS%QI(k))
  
  COMS%VEL_P(k) =  COMS%VEL_P(k) + COMS%VEL_T(k) * COMS%DT  
  coms%rad_p(k) =  coms%rad_p(k) + coms%rad_t(k) * COMS%DT  
! COMS%SC(k) = max(0., COMS%SC(k))

 enddo
endif
end subroutine update_plumerise
!-------------------------------------------------------------------------------
!
subroutine fallpart(coms,m1)
!use module_zero_plumegen_coms
implicit none
type(plumegen_coms), pointer :: coms
integer m1,k
real(kind=kind_phys) vtc, dfhz,dfiz,dz1
!srf==================================
!   verificar se o gradiente esta correto 
!  
!srf==================================
!
!     XNO=1.E7  [m**-4] median volume diameter raindrop,Kessler
!     VC = 38.3/(XNO**.125), median volume fallspeed eqn., Kessler
!     for ice, see (OT18), use F0=0.75 per argument there. coms%rho*q
!     values are in g/m**3, velocities in m/s

real(kind=kind_phys), PARAMETER :: VCONST = 5.107387, EPS = 0.622, F0 = 0.75  
real(kind=kind_phys), PARAMETER :: G = 9.81, CP = 1004.
!
do k=2,m1-1
!sam   if(.not. coms%rho(k)>1e-20) then
!sam 33 format('(33) Bad density at k=',I0,' rho=',F12.5,' T=',F12.5,' PE=',F12.5,' test=',I0)
!sam     write(0,33) k,coms%rho(k),coms%T(k),coms%PE(k),coms%testval
!sam   endif
   VTC = VCONST * COMS%RHO (k) **.125   ! median volume fallspeed (KTable4)
                                
!  hydrometeor assembly velocity calculations (K Table4)
!  COMS%VTH(k)=-VTC*COMS%QH(k)**.125  !median volume fallspeed, water            
   COMS%VTH (k) = - 4.	    !small variation with coms%qh
   
   COMS%VHREL = COMS%W (k) + COMS%VTH (k)  !relative to surrounding cloud
 
!  rain ventilation coefficient for evaporation
   COMS%CVH(k) = 1.6 + 0.57E-3 * (ABS (COMS%VHREL) ) **1.5  
!
!  COMS%VTI(k)=-VTC*F0*COMS%QI(k)**.125    !median volume fallspeed,ice             
   COMS%VTI (k) = - 3.                !small variation with coms%qi

   COMS%VIREL = COMS%W (k) + COMS%VTI (k)       !relative to surrounding cloud
!
!  ice ventilation coefficient for sublimation
   COMS%CVI(k) = 1.6 + 0.57E-3 * (ABS (COMS%VIREL) ) **1.5 / F0  
!
!
   IF (COMS%VHREL.GE.0.0) THEN  
    DFHZ=COMS%QH(k)*(COMS%RHO(k  )*COMS%VTH(k  )-COMS%RHO(k-1)*COMS%VTH(k-1))/COMS%RHO(k-1)
   ELSE  
    DFHZ=COMS%QH(k)*(COMS%RHO(k+1)*COMS%VTH(k+1)-COMS%RHO(k  )*COMS%VTH(k  ))/COMS%RHO(k)
   ENDIF  
   !
   !
   IF (COMS%VIREL.GE.0.0) THEN  
    DFIZ=COMS%QI(k)*(COMS%RHO(k  )*COMS%VTI(k  )-COMS%RHO(k-1)*COMS%VTI(k-1))/COMS%RHO(k-1)
   ELSE  
    DFIZ=COMS%QI(k)*(COMS%RHO(k+1)*COMS%VTI(k+1)-COMS%RHO(k  )*COMS%VTI(k  ))/COMS%RHO(k)
   ENDIF
   
   DZ1=COMS%ZM(K)-COMS%ZM(K-1)
   
   coms%qht(k) = coms%qht(k) - DFHZ / DZ1 !hydrometeors don't
   coms%qit(k) = coms%qit(k) - DFIZ / DZ1  !nor does ice? hail, what about

enddo
end subroutine fallpart

! *********************************************************************
SUBROUTINE WATERBAL(coms)
implicit none
type(plumegen_coms), pointer :: coms

!use module_zero_plumegen_coms  
!
                                        
IF (COMS%QC (COMS%L) .LE.1.0E-10) COMS%QC (COMS%L) = 0.  !DEFEAT UNDERFLOW PROBLEM
IF (COMS%QH (COMS%L) .LE.1.0E-10) COMS%QH (COMS%L) = 0.  
IF (COMS%QI (COMS%L) .LE.1.0E-10) COMS%QI (COMS%L) = 0.  
!
CALL EVAPORATE(COMS)    !vapor to cloud,cloud to vapor  
!                             
CALL SUBLIMATE(COMS)    !vapor to ice  
!                            
CALL GLACIATE(COMS)     !rain to ice 
                           
CALL MELT(COMS)         !ice to rain
!         
!if(ak1 > 0. .or. ak2 > 0.) &
CALL CONVERT(COMS) !(auto)conversion and accretion 
!CALL CONVERT2 () !(auto)conversion and accretion 
!

RETURN  
END SUBROUTINE WATERBAL
! *********************************************************************
SUBROUTINE EVAPORATE(coms)
!
!- evaporates cloud,rain and ice to saturation
!
!use module_zero_plumegen_coms  
implicit none
type(plumegen_coms), pointer :: coms
!
!     XNO=10.0E06
!     HERC = 1.93*1.E-6*XN035        !evaporation constant
!
real(kind=kind_phys), PARAMETER :: HERC = 5.44E-4, CP = 1.004, HEATCOND = 2.5E3  
real(kind=kind_phys), PARAMETER :: HEATSUBL = 2834., TMELT = 273., TFREEZE = 269.3

real(kind=kind_phys), PARAMETER :: FRC = HEATCOND / CP, SRC = HEATSUBL / CP

real(kind=kind_phys) :: evhdt, evidt, evrate, evap, sd,	quant, dividend, divisor, devidt

!
!
SD = COMS%QSAT (COMS%L) - COMS%QV (COMS%L)  !vapor deficit
IF (SD.EQ.0.0)  RETURN  
!IF (abs(SD).lt.1.e-7)  RETURN  


EVHDT = 0.  
EVIDT = 0.  
!evrate =0.; evap=0.; sd=0.0; quant=0.0; dividend=0.0; divisor=0.0; devidt=0.0
                                 
EVRATE = ABS (COMS%WBAR * COMS%DQSDZ)   !evaporation rate (Kessler 8.32)
EVAP = EVRATE * COMS%DT            !what we can get in DT
                                  

IF (SD.LE.0.0) THEN  !     condense. SD is negative

   IF (EVAP.GE.ABS (SD) ) THEN    !we get it all
                                  
      COMS%QC (COMS%L) = COMS%QC  (COMS%L) - SD  !deficit,remember?
      COMS%QV (COMS%L) = COMS%QSAT(COMS%L)       !set the vapor to saturation  
      COMS%T  (COMS%L) = COMS%T   (COMS%L) - SD * FRC  !heat gained through condensation
                                !per gram of dry air
      RETURN  

   ELSE  
                                 
      COMS%QC (COMS%L) = COMS%QC (COMS%L) + EVAP         !get what we can in DT 
      COMS%QV (COMS%L) = COMS%QV (COMS%L) - EVAP         !remove it from the vapor
      COMS%T  (COMS%L) = COMS%T  (COMS%L) + EVAP * FRC   !get some heat

      RETURN  

   ENDIF  
!
ELSE                                !SD is positive, need some water
!
! not saturated. saturate if possible. use everything in order
! cloud, rain, ice. SD is positive
                                         
   IF (EVAP.LE.COMS%QC (COMS%L) ) THEN        !enough cloud to last DT  
!
                                         
      IF (SD.LE.EVAP) THEN          !enough time to saturate
                                         
         COMS%QC (COMS%L) = COMS%QC (COMS%L) - SD       !remove cloud                                          
         COMS%QV (COMS%L) = COMS%QSAT (COMS%L)          !saturate
         COMS%T (COMS%L) = COMS%T (COMS%L) - SD * FRC   !cool the parcel                                          
         RETURN  !done
!
                                         
      ELSE   !not enough time
                                        
         SD = SD-EVAP               !use what there is
         COMS%QV (COMS%L) = COMS%QV (COMS%L) + EVAP     !add vapor
         COMS%T (COMS%L) = COMS%T (COMS%L) - EVAP * FRC !lose heat
         COMS%QC (COMS%L) = COMS%QC (COMS%L) - EVAP     !lose cloud
	                            !go on to rain.                                      
      ENDIF     
!
   ELSE                !not enough cloud to last DT
!      
      IF (SD.LE.COMS%QC (COMS%L) ) THEN   !but there is enough to sat
                                          
         COMS%QV (COMS%L) = COMS%QSAT (COMS%L)  !use it
         COMS%QC (COMS%L) = COMS%QC (COMS%L) - SD  
         COMS%T  (COMS%L) = COMS%T (COMS%L) - SD * FRC  
         RETURN  
	                              
      ELSE            !not enough to sat
         SD = SD-COMS%QC (COMS%L)  
         COMS%QV (COMS%L) = COMS%QV (COMS%L) + COMS%QC (COMS%L)  
         COMS%T  (COMS%L) = COMS%T (COMS%L) - COMS%QC (COMS%L) * FRC         
         COMS%QC (COMS%L) = 0.0  !all gone
                                          
      ENDIF       !on to rain                           
   ENDIF          !finished with cloud
!
!  but still not saturated, so try to use some rain
!  this is tricky, because we only have time DT to evaporate. if there
!  is enough rain, we can evaporate it for dt. ice can also sublimate
!  at the same time. there is a compromise here.....use rain first, then
!  ice. saturation may not be possible in one DT time.
!  rain evaporation rate (W12),(OT25),(K Table 4). evaporate rain first
!  sd is still positive or we wouldn't be here.


   IF (COMS%QH (COMS%L) > 1.E-10) THEN

!srf-25082005
!  QUANT = ( COMS%QC (COMS%L)  + COMS%QV (COMS%L) - COMS%QSAT (COMS%L) ) * COMS%RHO (COMS%L)   !g/m**3
   QUANT = ( COMS%QSAT (COMS%L)- COMS%QC (COMS%L) - COMS%QV (COMS%L)   ) * COMS%RHO (COMS%L)   !g/m**3
!
   EVHDT = (COMS%DT * HERC * (QUANT) * (COMS%QH (COMS%L) * COMS%RHO (COMS%L) ) **.65) / COMS%RHO (COMS%L)
!             rain evaporation in time DT
                                         
   IF (EVHDT.LE.COMS%QH (COMS%L) ) THEN           !enough rain to last DT

      IF (SD.LE.EVHDT) THEN  		!enough time to saturate	  
         COMS%QH (COMS%L) = COMS%QH (COMS%L) - SD   	!remove rain	  
         COMS%QV (COMS%L) = COMS%QSAT (COMS%L)  		!saturate	  
         COMS%T (COMS%L) = COMS%T (COMS%L) - SD * FRC  	!cool the parcel		  
	 
	 RETURN  			!done
!                       
      ELSE                               !not enough time
         SD = SD-EVHDT  		 !use what there is
         COMS%QV (COMS%L) = COMS%QV (COMS%L) + EVHDT  	 !add vapor
         COMS%T (COMS%L) = COMS%T (COMS%L) - EVHDT * FRC  	 !lose heat
         COMS%QH (COMS%L) = COMS%QH (COMS%L) - EVHDT  	 !lose rain

      ENDIF  				  !go on to ice.
!                                    
   ELSE  !not enough rain to last DT
!
      IF (SD.LE.COMS%QH (COMS%L) ) THEN             !but there is enough to sat
         COMS%QV (COMS%L) = COMS%QSAT (COMS%L)                !use it
         COMS%QH (COMS%L) = COMS%QH (COMS%L) - SD  
         COMS%T (COMS%L) = COMS%T (COMS%L) - SD * FRC  
         RETURN  
!                            
      ELSE                              !not enough to sat
         SD = SD-COMS%QH (COMS%L)  
         COMS%QV (COMS%L) = COMS%QV (COMS%L) + COMS%QH (COMS%L)  
         COMS%T (COMS%L) = COMS%T (COMS%L) - COMS%QH (COMS%L) * FRC    
         COMS%QH (COMS%L) = 0.0                   !all gone
                                          
      ENDIF                             !on to ice
!
                                          
   ENDIF                                !finished with rain
!
!
!  now for ice
!  equation from (OT); correction factors for units applied
!
   ENDIF
   IF (COMS%QI (COMS%L) .LE.1.E-10) RETURN            !no ice there
!
   DIVIDEND = ( (1.E6 / COMS%RHO (COMS%L) ) **0.475) * (SD / COMS%QSAT (COMS%L) &
            - 1) * (COMS%QI (COMS%L) **0.525) * 1.13
   DIVISOR = 7.E5 + 4.1E6 / (10. * COMS%EST (COMS%L) )  
                                                 
   DEVIDT = - COMS%CVI(COMS%L) * DIVIDEND / DIVISOR   !rate of change
                                                  
   EVIDT = DEVIDT * COMS%DT                      !what we could get
!
! logic here is identical to rain. could get fancy and make subroutine
! but duplication of code is easier. God bless the screen editor.
!
                                         
   IF (EVIDT.LE.COMS%QI (COMS%L) ) THEN             !enough ice to last DT
!
                                         
      IF (SD.LE.EVIDT) THEN  		  !enough time to saturate
         COMS%QI (COMS%L) = COMS%QI (COMS%L) - SD   	  !remove ice
         COMS%QV (COMS%L) = COMS%QSAT (COMS%L)  		  !saturate
         COMS%T (COMS%L) = COMS%T (COMS%L) - SD * SRC  	  !cool the parcel
	 
         RETURN  			  !done
!
                                          
      ELSE                                !not enough time
                                          
         SD = SD-EVIDT  		  !use what there is
         COMS%QV (COMS%L) = COMS%QV (COMS%L) + EVIDT  	  !add vapor
          COMS%T (COMS%L) =  COMS%T (COMS%L) - EVIDT * SRC  	  !lose heat
         COMS%QI (COMS%L) = COMS%QI (COMS%L) - EVIDT  	  !lose ice
                                          
      ENDIF  				  !go on,unsatisfied
!                                          
   ELSE                                   !not enough ice to last DT
!                                         
      IF (SD.LE.COMS%QI (COMS%L) ) THEN             !but there is enough to sat
                                          
         COMS%QV (COMS%L) = COMS%QSAT (COMS%L)                !use it
         COMS%QI (COMS%L) = COMS%QI   (COMS%L) - SD  
          COMS%T (COMS%L) =  COMS%T   (COMS%L) - SD * SRC  
	  
         RETURN  
!
      ELSE                                 !not enough to sat
         SD = SD-COMS%QI (COMS%L)  
         COMS%QV (COMS%L) = COMS%QV (COMS%L) + COMS%QI (COMS%L)  
         COMS%T (COMS%L) = COMS%T (COMS%L) - COMS%QI (COMS%L) * SRC             
         COMS%QI (COMS%L) = 0.0                      !all gone

      ENDIF                                !on to better things
                                           !finished with ice
   ENDIF  
!                                 
ENDIF                                      !finished with the SD decision
!
RETURN  
!
END SUBROUTINE EVAPORATE
!
! *********************************************************************
SUBROUTINE CONVERT (coms)
!
!- ACCRETION AND AUTOCONVERSION
!
implicit none
type(plumegen_coms), pointer :: coms

!use module_zero_plumegen_coms  
!
real(kind=kind_phys),      PARAMETER ::  AK1 = 0.001    !conversion rate constant
real(kind=kind_phys),      PARAMETER ::  AK2 = 0.0052   !collection (accretion) rate
real(kind=kind_phys),      PARAMETER ::  TH  = 0.5      !Kessler threshold
integer,   PARAMETER ::  iconv = 1        !- Kessler conversion (=0)
                                     
!real(kind=kind_phys), parameter :: ANBASE =  50.!*1.e+6 !Berry-number at cloud base #/m^3(maritime)
 real(kind=kind_phys), parameter :: ANBASE =100000.!*1.e+6 !Berry-number at cloud base #/m^3(continental)
!real(kind=kind_phys), parameter :: BDISP = 0.366       !Berry--size dispersion (maritime)
 real(kind=kind_phys), parameter :: BDISP = 0.146       !Berry--size dispersion (continental)
real(kind=kind_phys), parameter :: TFREEZE = 269.3  !ice formation temperature
!
real(kind=kind_phys) ::   accrete, con, q, h, bc1,   bc2,  total


IF (COMS%T (COMS%L)  .LE. TFREEZE) RETURN  !process not allowed above ice
!
IF (COMS%QC (COMS%L) .EQ. 0.     ) RETURN  

ACCRETE = 0.  
CON = 0.  
Q = COMS%RHO (COMS%L) * COMS%QC (COMS%L)  
H = COMS%RHO (COMS%L) * COMS%QH (COMS%L)  
!
!     selection rules
!                         
!            
IF (COMS%QH (COMS%L) .GT. 0.     ) ACCRETE = AK2 * Q * (H**.875)  !accretion, Kessler
!
IF (ICONV.NE.0) THEN   !select Berry or Kessler
!
!old   BC1 = 120.  
!old   BC2 = .0266 * ANBASE * 60.  
!old   CON = BDISP * Q * Q * Q / (BC1 * Q * BDISP + BC2) 	  

   CON = Q*Q*Q*BDISP/(60.*(5.*Q*BDISP+0.0366*ANBASE))
!
ELSE  
!                             
!   CON = AK1 * (Q - TH)   !Kessler autoconversion rate
!      
!   IF (CON.LT.0.0) CON = 0.0   !havent reached threshold
 
   CON = max(0.,AK1 * (Q - TH)) ! versao otimizada
!
ENDIF  
!
!
TOTAL = (CON + ACCRETE) * COMS%DT / COMS%RHO (COMS%L)  

!
IF (TOTAL.LT.COMS%QC (COMS%L) ) THEN  
!
   COMS%QC (COMS%L) = COMS%QC (COMS%L) - TOTAL  
   COMS%QH (COMS%L) = COMS%QH (COMS%L) + TOTAL    !no phase change involved
   RETURN  
!
ELSE  
!              
   COMS%QH (COMS%L) = COMS%QH (COMS%L) + COMS%QC (COMS%L)    !uses all there is
   COMS%QC (COMS%L) = 0.0  
!
ENDIF  
!
RETURN  
!
END SUBROUTINE CONVERT
!
!**********************************************************************
!
SUBROUTINE SUBLIMATE(coms)
!
implicit none
type(plumegen_coms), pointer :: coms

! ********************* VAPOR TO ICE (USE EQUATION OT22)***************
!use module_zero_plumegen_coms  
!
real(kind=kind_phys), PARAMETER :: EPS = 0.622, HEATFUS = 334., HEATSUBL = 2834., CP = 1.004
real(kind=kind_phys), PARAMETER :: SRC = HEATSUBL / CP, FRC = HEATFUS / CP, TMELT = 273.3
real(kind=kind_phys), PARAMETER :: TFREEZE = 269.3

real(kind=kind_phys) ::dtsubh,  dividend,divisor, subl
!
DTSUBH = 0.  
!
!selection criteria for sublimation
IF (COMS%T (COMS%L)  .GT. TFREEZE  ) RETURN  
IF (COMS%QV (COMS%L) .LE. COMS%QSAT (COMS%L) ) RETURN  
!
!     from (OT); correction factors for units applied
!
 DIVIDEND = ( (1.E6 / COMS%RHO (COMS%L) ) **0.475) * (COMS%QV (COMS%L) / COMS%QSAT (COMS%L) &
            - 1) * (COMS%QI (COMS%L) **0.525) * 1.13
 DIVISOR = 7.E5 + 4.1E6 / (10. * COMS%EST (COMS%L) )  
!
                                         
 DTSUBH = ABS (DIVIDEND / DIVISOR)   !sublimation rate
 SUBL = DTSUBH * COMS%DT                  !and amount possible
!
!     again check the possibilities
!
IF (SUBL.LT.COMS%QV (COMS%L) ) THEN  
!
   COMS%QV (COMS%L) = COMS%QV (COMS%L) - SUBL             !lose vapor
   COMS%QI (COMS%L) = COMS%QI (COMS%L) + SUBL             !gain ice
   COMS%T (COMS%L) = COMS%T (COMS%L) + SUBL * SRC         !energy change, warms air

   RETURN  
!
ELSE  
!                                     
   COMS%QI (COMS%L) = COMS%QV (COMS%L)                    !use what there is
   COMS%T  (COMS%L) = COMS%T (COMS%L) + COMS%QV (COMS%L) * SRC      !warm the air
   COMS%QV (COMS%L) = 0.0  
!
ENDIF  
!
RETURN  
END SUBROUTINE SUBLIMATE
!
! *********************************************************************
!
SUBROUTINE GLACIATE  (coms)
!
! *********************** CONVERSION OF RAIN TO ICE *******************
!     uses equation OT 16, simplest. correction from W not applied, but
!     vapor pressure differences are supplied.
!
!use module_zero_plumegen_coms  
!
implicit none
type(plumegen_coms), pointer :: coms
real(kind=kind_phys), PARAMETER :: HEATFUS = 334., CP = 1.004, EPS = 0.622, HEATSUBL = 2834.
real(kind=kind_phys), PARAMETER :: FRC = HEATFUS / CP, FRS = HEATSUBL / CP, TFREEZE =  269.3
real(kind=kind_phys), PARAMETER :: GLCONST = 0.025   !glaciation time constant, 1/sec
real(kind=kind_phys) dfrzh
!
                                    
 DFRZH = 0.    !rate of mass gain in ice
!
!selection rules for glaciation
IF (COMS%QH (COMS%L) .LE. 0.       ) RETURN  
IF (COMS%QV (COMS%L) .LT. COMS%QSAT (COMS%L) ) RETURN                                        
IF (COMS%T  (COMS%L) .GT. TFREEZE  ) RETURN  
!
!      NT=TMELT-COMS%T(COMS%L)
!      IF (NT.GT.50) NT=50
!
                                    
 DFRZH = COMS%DT * GLCONST * COMS%QH (COMS%L)    ! from OT(16)
!
IF (DFRZH.LT.COMS%QH (COMS%L) ) THEN  
!
   COMS%QI (COMS%L) = COMS%QI (COMS%L) + DFRZH  
   COMS%QH (COMS%L) = COMS%QH (COMS%L) - DFRZH  
   COMS%T (COMS%L) = COMS%T (COMS%L) + FRC * DFRZH  !warms air
   
   
   RETURN  
!
ELSE  
!
   COMS%QI (COMS%L) = COMS%QI (COMS%L) + COMS%QH (COMS%L)  
   COMS%T  (COMS%L) = COMS%T  (COMS%L) + FRC * COMS%QH (COMS%L)  
   COMS%QH (COMS%L) = 0.0  

 !print*,'8',coms%l,coms%qi(coms%l), COMS%QH (COMS%L)  
!
ENDIF  
!
RETURN  
!
END SUBROUTINE GLACIATE
!
!
! *********************************************************************
SUBROUTINE MELT(coms)
!
! ******************* MAKES WATER OUT OF ICE **************************
!use module_zero_plumegen_coms  
!                                              
implicit none
type(plumegen_coms), pointer :: coms

real(kind=kind_phys), PARAMETER :: FRC = 332.27, TMELT = 273., F0 = 0.75   !ice velocity factor
real(kind=kind_phys) DTMELT
!                                    
 DTMELT = 0.   !conversion,ice to rain
!
!selection rules
IF (COMS%QI (COMS%L) .LE. 0.0  ) RETURN  
IF (COMS%T (COMS%L)  .LT. TMELT) RETURN  
!
                                                      !OT(23,24)
 DTMELT = COMS%DT * (2.27 / COMS%RHO (COMS%L) ) * COMS%CVI(COMS%L) * (COMS%T (COMS%L) - TMELT) * ( (COMS%RHO(COMS%L)  &
         * COMS%QI (COMS%L) * 1.E-6) **0.525) * (F0** ( - 0.42) )
                                                      !after Mason,1956
!
!     check the possibilities
!
IF (DTMELT.LT.COMS%QI (COMS%L) ) THEN  
!
   COMS%QH (COMS%L) = COMS%QH (COMS%L) + DTMELT  
   COMS%QI (COMS%L) = COMS%QI (COMS%L) - DTMELT  
   COMS%T  (COMS%L) = COMS%T (COMS%L) - FRC * DTMELT     !cools air
   
   RETURN  
!
ELSE  
!
   COMS%QH (COMS%L) = COMS%QH (COMS%L) + COMS%QI (COMS%L)   !get all there is to get
   COMS%T  (COMS%L) = COMS%T (COMS%L) - FRC * COMS%QI (COMS%L)  
   COMS%QI (COMS%L) = 0.0  
!
ENDIF  
!
RETURN  
!
END SUBROUTINE MELT

SUBROUTINE htint (nzz1, vctra, eleva, nzz2, vctrb, elevb, errmsg, errflg)
  IMPLICIT NONE
  INTEGER, INTENT(IN ) :: nzz1
  INTEGER, INTENT(IN ) :: nzz2
  REAL(kind=kind_phys),    INTENT(IN ) :: vctra(nzz1)
  REAL(kind=kind_phys),    INTENT(OUT) :: vctrb(nzz2)
  REAL(kind=kind_phys),    INTENT(IN ) :: eleva(nzz1)
  REAL(kind=kind_phys),    INTENT(IN ) :: elevb(nzz2)
  character(*), intent(inout) :: errmsg
  integer, intent(inout) :: errflg
  INTEGER :: l
  INTEGER :: k
  INTEGER :: kk
  REAL(kind=kind_phys)    :: wt

  l=1

  DO k=1,nzz2
     DO
        IF ( (elevb(k) <  eleva(1)) .OR. &
             ((elevb(k) >= eleva(l)) .AND. (elevb(k) <= eleva(l+1))) ) THEN
           wt       = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
           vctrb(k) = vctra(l)+(vctra(l+1)-vctra(l))*wt
           EXIT
        ELSE IF ( elevb(k) >  eleva(nzz1))  THEN
           wt       = (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
           vctrb(k) = vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
           EXIT
        END IF

        l=l+1
        IF(l == nzz1) THEN
           PRINT *,'htint:nzz1',nzz1
           DO kk=1,l
              PRINT*,'kk,eleva(kk),elevb(kk)',kk,eleva(kk),elevb(kk)
           END DO
           errmsg='htint assertion failure (see print for details)'
           errflg=1
        END IF
     END DO
  END DO
END SUBROUTINE htint
!-----------------------------------------------------------------------------
FUNCTION ESAT_PR (TEM)  
!
! ******* Vapor Pressure  A.L. Buck JAM V.20 p.1527. (1981) ***********
!
real(kind=kind_phys), PARAMETER :: CI1 = 6.1115, CI2 = 22.542, CI3 = 273.48
real(kind=kind_phys), PARAMETER :: CW1 = 6.1121, CW2 = 18.729, CW3 = 257.87, CW4 = 227.3
real(kind=kind_phys), PARAMETER :: TMELT = 273.3

real(kind=kind_phys) ESAT_PR
real(kind=kind_phys) temc , tem,esatm
!
!     formulae from Buck, A.L., JAM 20,1527-1532
!     custom takes esat wrt water always. formula for h2o only
!     good to -40C so:
!
!
TEMC = TEM - TMELT  
IF (TEMC<= - 40.0) then
  ESATM = CI1 * EXP (CI2 * TEMC / (TEMC + CI3) )  !ice, millibars  
  ESAT_PR = ESATM / 10.	!kPa			  

  RETURN  
ENDIF
!
ESATM = CW1 * EXP ( ( (CW2 - (TEMC / CW4) ) * TEMC) / (TEMC + CW3))
                          
ESAT_PR = ESATM / 10.	!kPa			  
RETURN  
END function ESAT_PR
!     ******************************************************************

!      ------------------------------------------------------------------------
END Module module_smoke_plumerise

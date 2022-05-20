!>\file dep_dry_mod.F90
!! This file is for the dry depostion driver.

module dep_dry_mod

  use machine ,        only : kind_phys
  use rrfs_smoke_config, only : epsilc, GOCART_SIMPLE => CHEM_OPT_GOCART, CTRA_OPT_NONE
!  use chem_tracers_mod, only : p_o3,p_dust_1,p_vash_1,p_vash_4,p_vash_10,p_dms,
!  &
!            config_flags => chem_config
  use dep_dry_gocart_mod
  use dep_simple_mod
  use dep_vertmx_mod
! use aero_soa_vbs_mod, only : soa_vbs_depdriver

  implicit none


  private

  public :: dry_dep_driver

contains

    subroutine dry_dep_driver(data,ktau,dtstep,julday,current_month,t_phy,p_phy,      &
               moist,p8w,rmol,alt,gmt,t8w,raincv,                         &
               chem,rho_phy,dz8w,exch_h,hfx,                              &
               ivgtyp,tsk,gsw,vegfra,pbl,ust,znt,z,z_at_w,                &
               xland,xlat,xlong,h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,     &
               anh3,ddep,dep_vel_o3,g,                                    &
               e_co,kemit,snowh,numgas,                                   &
               num_chem,num_moist,                                        &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
!----------------------------------------------------------------------
! USE module_model_constants
! USE module_configure
! USE module_state_description
! USE module_dep_simple
! USE module_initial_chem_namelists,only:p_o3,p_dust_1,p_vash_1,p_vash_4,p_vash_10,p_dms
! USE module_vertmx_wrf
! USE module_chemvars,only:epsilc
! USE module_data_sorgam
! USE module_aerosols_sorgam
! USE module_gocart_settling
! use module_dep_simple
! USE module_gocart_drydep,only: gocart_drydep_driver
! USE module_aerosols_soa_vbs, only: soa_vbs_depdriver
! USE module_mosaic_drydep, only:  mosaic_drydep_driver
! USE module_mixactivate_wrappers, only: mosaic_mixactivate, sorgam_mixactivate
  IMPLICIT NONE
  type(smoke_data), pointer, intent(inout) :: data

   INTEGER,      INTENT(IN   ) :: numgas, current_month,         &
               num_chem,num_moist, julday,                       &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte
   INTEGER,      INTENT(IN   ) ::                               &
                                  ktau
   REAL(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
         INTENT(IN ) ::                                   moist
   REAL(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),         &
         INTENT(INOUT ) ::                                 chem

   INTEGER,      INTENT(IN   ) :: kemit
   REAL(kind_phys), DIMENSION( ims:ime, kms:kemit, jms:jme ),            &
         INTENT(IN ) ::                                                    &
          e_co




   REAL(kind_phys),  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                        alt,    &
                                                        t8w,    &
                                                      dz8w,     &
                                              p8w,z_at_w ,  &
                                              exch_h,rho_phy,z
 REAL(kind_phys),  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(INOUT) ::                                      &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,    &
          INTENT(IN   ) ::                                      &
                                                     ivgtyp
   REAL(kind_phys),  DIMENSION( ims:ime , jms:jme )                   ,    &
          INTENT(INOUT) ::                                      &
                                                     tsk,       &
                                                     gsw,       &
                                                     vegfra,     &
                                                     pbl,       &
                                                     snowh,       &
                                                     raincv,       &
                                                     ust,       &
                                                     hfx,       &
                                                   xland,       &
                                                   xlat,        &
                                                   xlong,       &
                                                    znt,rmol
   REAL(kind_phys), DIMENSION( ims:ime, jms:jme, num_chem ),         &
         INTENT(OUT ) ::                                   ddep
   REAL(kind_phys),  DIMENSION( ims:ime , jms:jme )                   ,    &
          INTENT(OUT) ::                                      &
                                                     dep_vel_o3
       REAL(kind_phys),  DIMENSION( ims:ime , kms:kme , jms:jme ),              &
               INTENT(IN   ) ::                                      &
                                                           p_phy,    &
                                                           t_phy

      REAL(kind_phys),      INTENT(IN   ) ::                               &
                             dtstep,g,gmt

!--- deposition and emissions stuff
! .. Parameters ..
! ..
! .. Local Scalars ..

      REAL(kind_phys) :: cdt, factor

      INTEGER :: idrydep_onoff

!      INTEGER :: chem_conv_tr, chem_opt

!     CHARACTER (4) :: luse_typ,mminlu_loc
! ..
! .. Local Arrays ..
   REAL(kind_phys), DIMENSION( its:ite, jts:jte, num_chem ) ::   ddvel

!  REAL(kind_phys),  DIMENSION( ims:ime , kms:kme , jms:jme ) :: dryrho_phy
   REAL(kind_phys),  DIMENSION( kts:kte ) :: dryrho_1d

! turbulent transport
      real(kind_phys) :: pblst(kts:kte),ekmfull(kts:kte+1),zzfull(kts:kte+1),zz(kts:kte)
      integer :: i,j,k,nv
!
! necessary for aerosols (module dependent)
!
  REAL(kind_phys), DIMENSION( its:ite, jts:jte ) ::   aer_res
  REAL(kind_phys), DIMENSION( its:ite, jts:jte ) ::   aer_res_def
  REAL(kind_phys), DIMENSION( its:ite, jts:jte ) ::   aer_res_zcen

! ..
! .. Intrinsic Functions ..
      INTRINSIC max, min

!      chem_opt     = chem_opt
!      chem_conv_tr = chem_conv_tr

!
! compute dry deposition velocities = ddvel
!
! 28-jun-2005 rce - initialize ddvel=0; call aerosol drydep routine
!           only when drydep_opt == WESELY
!       the wesely_driver routine computes aer_res, and currently
!	    you cannot compute aerosol drydep without it !!
! 08-jul-2005 rce - pass idrydep_onoff to mixactivate routines
!
!  write(6,*)'call dry dep driver'
   dep_vel_o3(:,:)=0.
   ddvel(:,:,:) = 0.0
   idrydep_onoff = 0

!  drydep_select: SELECT CASE(drydep_opt)

!    CASE ( WESELY )
!
! drydep_opt == WESELY means 
!     wesely for gases 
!     other (appropriate) routine for aerosols
!
!      CALL wrf_debug(15,'DOING DRY DEP VELOCITIES WITH WESELY METHOD')

      IF( chem_opt /= GOCART_SIMPLE ) THEN
         call wesely_driver(data,ktau,dtstep,                                  &
              current_month,                                              &
              gmt,julday,t_phy,moist,p8w,t8w,raincv,                      &
              p_phy,chem,rho_phy,dz8w,ddvel,aer_res_def,aer_res_zcen,    &
              ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
              snowh,numgas,                                                    &
              ids,ide, jds,jde, kds,kde,                                 &
              ims,ime, jms,jme, kms,kme,                                 &
              its,ite, jts,jte, kts,kte                                  )
      ENDIF
       IF (( chem_opt == GOCART_SIMPLE ) .or.            &
              ( chem_opt == GOCARTRACM_KPP)  .or.            &
              ( chem_opt == 316)  .or.            &
              ( chem_opt == 317)  .or.            &
!             ( chem_opt == 502)  .or.            &
                (chem_opt == 304          )) then
!
! this does aerosol species (dust,seas, bc,oc) for gocart only
! this does aerosol species (dust,seas, bc,oc,sulf) for gocart only
!,
         call gocart_drydep_driver(numgas,                     &
               moist,p8w,chem,rho_phy,dz8w,ddvel,xland,hfx,    &
               ivgtyp,tsk,pbl,ust,znt,g,                       &
               num_moist,num_chem,                             &
               ids,ide, jds,jde, kds,kde,                      &
               ims,ime, jms,jme, kms,kme,                      &
               its,ite, jts,jte, kts,kte                       )
       ELSE if (chem_opt == 501 ) then
! for caesium .1cm/s
!
       ddvel(:,:,:)=.001
       ELSE if (chem_opt == 108 ) then
!!       call soa_vbs_depdriver (ust,t_phy,                    &
!!               moist,p8w,rmol,znt,pbl,           &
!!               alt,p_phy,chem,rho_phy,dz8w,                    &
!!               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,      &
!!               aer_res,ddvel(:,:,numgas+1:num_chem),           &
!!               num_chem-numgas,                                &
!!               ids,ide, jds,jde, kds,kde,                      &
!!               ims,ime, jms,jme, kms,kme,                      &
!!               its,ite, jts,jte, kts,kte                       )
! limit aerosol ddvels to <= 0.5 m/s
! drydep routines occasionally produce unrealistically-large particle
! diameter leading to unrealistically-large sedimentation velocity 
       ddvel(:,:,numgas+1:num_chem) = min( 0.50, ddvel(:,:,numgas+1:num_chem))
       ELSE
          !Set dry deposition velocity to zero when using the
          !chemistry tracer mode.
          ddvel(:,:,:) = 0.
       END IF
       idrydep_onoff = 1

!
!   Compute dry deposition according to NGAC
!
    cdt = real(dtstep, kind=kind_phys)
    do nv = 1, num_chem
      do j = jts, jte
        do i = its, ite
          factor = 1._kind_phys - exp(-ddvel(i,j,nv)*cdt/dz8w(i,kts,j))
          ddep(i,j,nv) = max(0.0, factor * chem(i,kts,j,nv)) & !ug/m2/s
                         * (p8w(i,kts,j)-p8w(i,kts+1,j))/g/dtstep
        end do
      end do
    end do


!   This will be called later from subgrd_transport_driver.F !!!!!!!!
!
!
      do 100 j=jts,jte
      do 100 i=its,ite
      if(p_dust_1.gt.1)dep_vel_o3(i,j)=ddvel(i,j,p_dust_1)
      pblst=0.
!
!
!-- start with vertical mixing
!
      do k=kts,kte+1
         zzfull(k)=z_at_w(i,k,j)-z_at_w(i,kts,j)
      enddo

      if (chem_conv_tr == CTRA_OPT_NONE) then
        ekmfull = 0.
      else
        ekmfull(kts)=0.
        do k=kts+1,kte
         ekmfull(k)=max(1.e-6,exch_h(i,k,j))
        enddo
        ekmfull(kte+1)=0.
      end if

!!$! UNCOMMENT THIS AND FINE TUNE LEVELS TO YOUR DOMAIN IF YOU WANT TO
!!$! FORCE MIXING TO A CERTAIN DEPTH:
!!$!
!!$! --- Mix the emissions up several layers
!
      do k=kts,kte
         zz(k)=z(i,k,j)-z_at_w(i,kts,j)
      enddo
!   vertical mixing routine (including deposition)
!   need to be careful here with that dumm tracer in spot 1
!   do not need lho,lho2
!   (03-may-2006 rce - calc dryrho_1d and pass it to vertmx)
!
!     if(p_o3.gt.1)dep_vel_o3(i,j)=ddvel(i,j,p_o3)
      do nv=1,num_chem-0
         do k=kts,kte
            pblst(k)=max(epsilc,chem(i,k,j,nv))
            dryrho_1d(k) = 1./alt(i,k,j)
         enddo

         !mix_select: SELECT CASE(chem_opt)
         !CASE DEFAULT
               call vertmx(data,dtstep,pblst,ekmfull,dryrho_1d, &
                           zzfull,zz,ddvel(i,j,nv),kts,kte)

         !END SELECT mix_select

         do k=kts,kte 
            chem(i,k,j,nv)=max(epsilc,pblst(k))
         enddo
      enddo
100   continue

END SUBROUTINE dry_dep_driver

end module dep_dry_mod

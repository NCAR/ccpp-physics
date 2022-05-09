!>\file rrfs_smoke_lsdep_wrapper.F90
!! This file is RRFS-smoke large-scale wet deposition wrapper with CCPP
!! Haiqin.Li@noaa.gov 04/2021

 module rrfs_smoke_lsdep_wrapper

   use machine ,        only : kind_phys
   use rrfs_smoke_config
   use dep_wet_ls_mod
   use dust_data_mod
   use rrfs_smoke_data

   implicit none

   private

   public :: rrfs_smoke_lsdep_wrapper_run

contains

!>\defgroup rrfs_smoke_lsdep_wrapper GSD Chem driver Module  
!> \ingroup gsd_chem_group
!! This is the GSD Chem driver Module
!! \section arg_table_rrfs_smoke_lsdep_wrapper_run Argument Table
!! \htmlinclude rrfs_smoke_lsdep_wrapper_run.html
!!
!>\section rrfs_smoke_lsdep_wrapper GSD Chemistry Scheme General Algorithm
!> @{
    subroutine rrfs_smoke_lsdep_wrapper_run(im, kte, kme, ktau, dt,     &
                   rain_cpl, rainc_cpl, g,                              &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,  &
                   w, dqdt, ntrac,ntsmoke,ntdust,                       &
                   gq0,qgrs,wetdep_ls_opt_in,                    &
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau
    integer,        intent(in) :: ntrac,ntsmoke,ntdust
    real(kind_phys),intent(in) :: dt,g

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    real(kind_phys), dimension(:),     intent(in) :: rain_cpl, rainc_cpl
    real(kind_phys), dimension(:,:), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(:,:), intent(in) :: phl3d, prl3d, tk3d,        &
                us3d, vs3d, spechum, w, dqdt
    real(kind_phys), dimension(:,:,:), intent(inout) :: gq0, qgrs
    integer,           intent(in) :: wetdep_ls_opt_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,       &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, vvel, dqdti

    real(kind_phys), dimension(ims:im, jms:jme) :: rcav, rnav

!>- vapor & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  :: var_rmv

    integer :: ide, ime, ite, kde

    real(kind_phys) :: dtstep
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

    type(smoke_data), pointer :: data

!>-- local variables
    integer :: i, j, jp, k, kp, n

    data=>get_thread_smoke_data()

    errmsg = ''
    errflg = 0

    wetdep_ls_opt     = wetdep_ls_opt_in
    !print*,'hli wetdep_ls_opt',wetdep_ls_opt

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- initialize large-sacle wet depostion
    if (ktau==1) then
     call dep_wet_ls_init(data)
    endif

    ! -- set control flags

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
     rnav(i,1)=max((rain_cpl(i)-rainc_cpl(i))*1000., 0.) ! meter to mm
    enddo

!!!

!>- get ready for chemistry run
    call rrfs_smoke_prep_lsdep(data,ktau,dtstep,                        &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w, dqdt,           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,dqdti,z_at_w,vvel,g,                                        &
        ntsmoke,ntdust,                                                 &
        ntrac,gq0,num_chem, num_moist,                                  &
        ppm2ugkg,moist,chem,                                            &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

     ! -- ls wet deposition
     select case (wetdep_ls_opt)
       case (WDLS_OPT_GSD)
         call wetdep_ls(data,dt,chem,rnav,moist,rho_phy,var_rmv,        &
                        num_moist,num_chem,p_qc,p_qi,dz8w,vvel,         &
                        ids,ide, jds,jde, kds,kde,                      &
                        ims,ime, jms,jme, kms,kme,                      &
                        its,ite, jts,jte, kts,kte)
       case (WDLS_OPT_NGAC)
         call WetRemovalGOCART(data,its,ite, jts,jte, kts,kte, 1,1, dt, &
                               num_chem,var_rmv,chem,p_phy,t_phy,       &
                               rho_phy,dqdti,rcav,rnav, g,              &
                               ims,ime, jms,jme, kms,kme)
         !if (chem_rc_check(localrc, msg="Failure in NGAC wet removal scheme", &
         !  file=__FILE__, line=__LINE__, rc=rc)) return
       case default
         ! -- no further option implemented
    end select


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntsmoke)=ppm2ugkg(p_oc1   ) * max(epsilc,chem(i,k,1,p_oc1))
       gq0(i,k,ntdust )=ppm2ugkg(p_dust_1) * max(epsilc,chem(i,k,1,p_dust_1))
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntsmoke)=gq0(i,k,ntsmoke)
       qgrs(i,k,ntdust )=gq0(i,k,ntdust )
     enddo
    enddo


!
   end subroutine rrfs_smoke_lsdep_wrapper_run
!> @}

  subroutine rrfs_smoke_prep_lsdep(data,ktau,dtstep,                        &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,dqdt,           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                  &
        t8w,dqdti,z_at_w,vvel,g,                                       &
        ntsmoke,ntdust,                                                &
        ntrac,gq0,num_chem, num_moist,                                 &
        ppm2ugkg,moist,chem,                                           &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)
    implicit none
    type(smoke_data), intent(inout) :: data

    !Chem input configuration
    integer, intent(in) :: ktau
    real(kind=kind_phys), intent(in) :: dtstep,g

    !FV3 input variables
    integer, intent(in) :: ntrac,ntsmoke,ntdust
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,dqdt
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg
    
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel, dqdti
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w

    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys) ::  factor,factor2,pu,pl,aln,pwant
    real(kind_phys) ::  xhour,xmin,xlonn,xtime,real_time
    real(kind_phys), DIMENSION (1,1) :: sza,cosszax
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour

    ! -- initialize output arrays
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    vvel           = 0._kind_phys
    dqdti          = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys


    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,1)/g)
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=abs(ph3d(ip,kp+1)-ph3d(ip,kp))/g
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          p8w(i,k,j)=pr3d(ip,kp)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(ip,kkp)
          p_phy(i,k,j)=prl3d(ip,kkp)
          u_phy(i,k,j)=us3d(ip,kkp)
          dqdti(i,k,j)=dqdt(ip,kkp)
          v_phy(i,k,j)=vs3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(ip,kkp)*rri(i,k,j)/g 
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          !--
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=2,kte
        do i=its,ite
          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
        enddo
      enddo
    enddo

    ! -- only used in phtolysis....
    do j=jts,jte
      do i=its,ite
        t8w(i,1,j)=t_phy(i,1,j)
        t8w(i,kte+1,j)=t_phy(i,kte,j)
      enddo
    enddo

 
    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_oc1   )=max(epsilc,gq0(i,k,ntsmoke)/ppm2ugkg(p_oc1))
       chem(i,k,jts,p_dust_1)=max(epsilc,gq0(i,k,ntdust )/ppm2ugkg(p_dust_1))
     enddo
    enddo


  end subroutine rrfs_smoke_prep_lsdep
!> @}
  end module rrfs_smoke_lsdep_wrapper

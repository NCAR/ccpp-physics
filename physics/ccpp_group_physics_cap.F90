
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated cap module for the CCPP physics group
!!
!
module ccpp_group_physics_cap

   use lsm_noah, only: lsm_noah_init
   use GFS_suite_interstitial_phys_reset, only: GFS_suite_interstitial_phys_reset_run
   use GFS_suite_stateout_reset, only: GFS_suite_stateout_reset_run
   use get_prs_fv3, only: get_prs_fv3_run
   use GFS_suite_interstitial_1, only: GFS_suite_interstitial_1_run
   use dcyc2t3, only: dcyc2t3_run
   use GFS_surface_generic_pre, only: GFS_surface_generic_pre_run
   use GFS_suite_interstitial_2, only: GFS_suite_interstitial_2_run
   use sfc_ex_coef, only: sfc_ex_coef_run
   use GFS_surface_loop_control_part1, only: GFS_surface_loop_control_part1_run
   use sfc_nst_pre, only: sfc_nst_pre_run
   use sfc_nst, only: sfc_nst_run
   use sfc_nst_post, only: sfc_nst_post_run
   use lsm_noah, only: lsm_noah_run
   use sfc_sice, only: sfc_sice_run
   use GFS_surface_loop_control_part2, only: GFS_surface_loop_control_part2_run
   use dcyc2t3_post, only: dcyc2t3_post_run
   use sfc_diag, only: sfc_diag_run
   use sfc_diag_post, only: sfc_diag_post_run
   use GFS_surface_generic_post, only: GFS_surface_generic_post_run
   use GFS_PBL_generic_pre, only: GFS_PBL_generic_pre_run
   use hedmf, only: hedmf_run
   use GFS_PBL_generic_post, only: GFS_PBL_generic_post_run
   use gwdps_pre, only: gwdps_pre_run
   use gwdps, only: gwdps_run
   use gwdps_post, only: gwdps_post_run
   use rayleigh_damp, only: rayleigh_damp_run
   use GFS_suite_stateout_update, only: GFS_suite_stateout_update_run
   use ozphys, only: ozphys_run
   use GFS_DCNV_generic_pre, only: GFS_DCNV_generic_pre_run
   use get_phi_fv3, only: get_phi_fv3_run
   use GFS_suite_interstitial_3, only: GFS_suite_interstitial_3_run
   use samfdeepcnv, only: samfdeepcnv_run
   use GFS_DCNV_generic_post, only: GFS_DCNV_generic_post_run
   use gwdc_pre, only: gwdc_pre_run
   use gwdc, only: gwdc_run
   use gwdc_post, only: gwdc_post_run
   use GFS_SCNV_generic_pre, only: GFS_SCNV_generic_pre_run
   use samfshalcnv, only: samfshalcnv_run
   use samfshalcnv_post, only: samfshalcnv_post_run
   use GFS_SCNV_generic_post, only: GFS_SCNV_generic_post_run
   use GFS_suite_interstitial_4, only: GFS_suite_interstitial_4_run
   use cnvc90, only: cnvc90_run
   use GFS_MP_generic_pre, only: GFS_MP_generic_pre_run
   use zhaocarr_gscond, only: zhaocarr_gscond_run
   use zhaocarr_precpd, only: zhaocarr_precpd_run
   use GFS_MP_generic_post, only: GFS_MP_generic_post_run
   use sfc_sice_post, only: sfc_sice_post_run
   use lsm_noah, only: lsm_noah_finalize


   implicit none

   private
   public :: physics_init_cap,physics_run_cap,physics_finalize_cap

   logical, save :: initialized = .false.

   contains

   function physics_init_cap(cdata,GFS_Control) result(ierr)

      use ccpp_types, only: ccpp_t
      use GFS_typedefs, only: GFS_control_type

      implicit none

      integer                     :: ierr

      type(ccpp_t), intent(inout) :: cdata
      type(GFS_control_type), intent(in) :: GFS_Control

      ierr = 0


      if (initialized) return



      call lsm_noah_init(me=GFS_Control%me,isot=GFS_Control%isot,ivegsrc=GFS_Control%ivegsrc,nlunit=GFS_Control%nlunit, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in lsm_noah_init"
        ierr=cdata%errflg
        return
      end if

    


      initialized = .true.


   end function physics_init_cap

   function physics_run_cap(GFS_Control,GFS_Data,con_epsm1,con_hvap,GFS_Interstitial,con_rd,con_rv, &
                  con_g,con_cp,con_cvap,con_t0c,con_eps,con_pi,cdata,con_fvirt,con_cliq) result(ierr)

      use GFS_typedefs, only: GFS_control_type
      use GFS_typedefs, only: GFS_data_type
      use machine, only: kind_phys
      use GFS_typedefs, only: GFS_interstitial_type
      use ccpp_types, only: ccpp_t

      implicit none

      integer                     :: ierr

      type(GFS_control_type), intent(inout) :: GFS_Control
      type(GFS_data_type), intent(inout) :: GFS_Data(:)
      real(kind_phys), intent(in) :: con_epsm1
      real(kind_phys), intent(in) :: con_hvap
      type(GFS_interstitial_type), intent(inout) :: GFS_Interstitial(:)
      real(kind_phys), intent(in) :: con_rd
      real(kind_phys), intent(in) :: con_rv
      real(kind_phys), intent(in) :: con_g
      real(kind_phys), intent(in) :: con_cp
      real(kind_phys), intent(in) :: con_cvap
      real(kind_phys), intent(in) :: con_t0c
      real(kind_phys), intent(in) :: con_eps
      real(kind_phys), intent(in) :: con_pi
      type(ccpp_t), intent(inout) :: cdata
      real(kind_phys), intent(in) :: con_fvirt
      real(kind_phys), intent(in) :: con_cliq

      ierr = 0


      if (.not.initialized) then
        write(cdata%errmsg,'(*(a))') 'physics_run called before physics_init'
        cdata%errflg = 1
        return
      end if



      call GFS_suite_interstitial_phys_reset_run(Interstitial=GFS_Interstitial(cdata%thrd_no),Model=GFS_Control,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_suite_interstitial_phys_reset_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_suite_stateout_reset_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,ntrac=GFS_Control%ntrac, &
                  tgrs=GFS_Data(cdata%blk_no)%Statein%tgrs,ugrs=GFS_Data(cdata%blk_no)%Statein%ugrs, &
                  vgrs=GFS_Data(cdata%blk_no)%Statein%vgrs,qgrs=GFS_Data(cdata%blk_no)%Statein%qgrs, &
                  gt0=GFS_Data(cdata%blk_no)%Stateout%gt0,gu0=GFS_Data(cdata%blk_no)%Stateout%gu0, &
                  gv0=GFS_Data(cdata%blk_no)%Stateout%gv0,gq0=GFS_Data(cdata%blk_no)%Stateout%gq0, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_suite_stateout_reset_run"
        ierr=cdata%errflg
        return
      end if

    
      call get_prs_fv3_run(ix=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,phii=GFS_Data(cdata%blk_no)%Statein%phii, &
                  prsi=GFS_Data(cdata%blk_no)%Statein%prsi,tgrs=GFS_Data(cdata%blk_no)%Statein%tgrs, &
                  qgrs1=GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntqv),del=GFS_Interstitial(cdata%thrd_no)%del, &
                  del_gz=GFS_Interstitial(cdata%thrd_no)%del_gz,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in get_prs_fv3_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_suite_interstitial_1_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,ntrac=GFS_Control%ntrac, &
                  crtrh=GFS_Control%crtrh,dtf=GFS_Control%dtf,dtp=GFS_Control%dtp,slmsk=GFS_Data(cdata%blk_no)%Sfcprop%slmsk, &
                  area=GFS_Data(cdata%blk_no)%Grid%area,dxmin=GFS_Control%dxmin,dxinv=GFS_Control%dxinv, &
                  pgr=GFS_Data(cdata%blk_no)%Statein%pgr,rhbbot=GFS_Interstitial(cdata%thrd_no)%rhcbot, &
                  rhpbl=GFS_Interstitial(cdata%thrd_no)%rhcpbl,rhbtop=GFS_Interstitial(cdata%thrd_no)%rhctop, &
                  frain=GFS_Interstitial(cdata%thrd_no)%frain,islmsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  frland=GFS_Interstitial(cdata%thrd_no)%frland,work1=GFS_Interstitial(cdata%thrd_no)%work1, &
                  work2=GFS_Interstitial(cdata%thrd_no)%work2,psurf=GFS_Data(cdata%blk_no)%Intdiag%psurf, &
                  dudt=GFS_Interstitial(cdata%thrd_no)%dudt,dvdt=GFS_Interstitial(cdata%thrd_no)%dvdt, &
                  dtdt=GFS_Interstitial(cdata%thrd_no)%dtdt,dtdtc=GFS_Interstitial(cdata%thrd_no)%dtdtc, &
                  dqdt=GFS_Interstitial(cdata%thrd_no)%dqdt,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_suite_interstitial_1_run"
        ierr=cdata%errflg
        return
      end if

    
      call dcyc2t3_run(solhr=GFS_Control%solhr,slag=GFS_Control%slag,sdec=GFS_Control%sdec,cdec=GFS_Control%cdec, &
                  sinlat=GFS_Data(cdata%blk_no)%Grid%sinlat,coslat=GFS_Data(cdata%blk_no)%Grid%coslat, &
                  xlon=GFS_Data(cdata%blk_no)%Grid%xlon,coszen=GFS_Data(cdata%blk_no)%Radtend%coszen, &
                  tsea=GFS_Data(cdata%blk_no)%Sfcprop%tsfc,tf=GFS_Data(cdata%blk_no)%Statein%tgrs(:,1), &
                  tsflw=GFS_Data(cdata%blk_no)%Radtend%tsflw,sfcemis=GFS_Data(cdata%blk_no)%Radtend%semis, &
                  sfcdsw=GFS_Data(cdata%blk_no)%Coupling%sfcdsw,sfcnsw=GFS_Data(cdata%blk_no)%Coupling%sfcnsw, &
                  sfcdlw=GFS_Data(cdata%blk_no)%Coupling%sfcdlw,swh=GFS_Data(cdata%blk_no)%Tbd%htswc, &
                  swhc=GFS_Data(cdata%blk_no)%Tbd%htsw0,hlw=GFS_Data(cdata%blk_no)%Tbd%htlwc, &
                  hlwc=GFS_Data(cdata%blk_no)%Tbd%htlw0,sfcnirbmu=GFS_Data(cdata%blk_no)%Coupling%nirbmui, &
                  sfcnirdfu=GFS_Data(cdata%blk_no)%Coupling%nirdfui,sfcvisbmu=GFS_Data(cdata%blk_no)%Coupling%visbmui, &
                  sfcvisdfu=GFS_Data(cdata%blk_no)%Coupling%visdfui,sfcnirbmd=GFS_Data(cdata%blk_no)%Coupling%nirbmdi, &
                  sfcnirdfd=GFS_Data(cdata%blk_no)%Coupling%nirdfdi,sfcvisbmd=GFS_Data(cdata%blk_no)%Coupling%visbmdi, &
                  sfcvisdfd=GFS_Data(cdata%blk_no)%Coupling%visdfdi,ix=GFS_Control%blksz(cdata%blk_no), &
                  im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,deltim=GFS_Control%dtf, &
                  dtdt=GFS_Interstitial(cdata%thrd_no)%dtdt,dtdtc=GFS_Interstitial(cdata%thrd_no)%dtdtc, &
                  adjsfcdsw=GFS_Data(cdata%blk_no)%Intdiag%dswsfci,adjsfcnsw=GFS_Data(cdata%blk_no)%Intdiag%nswsfci, &
                  adjsfcdlw=GFS_Data(cdata%blk_no)%Intdiag%dlwsfci,adjsfculw=GFS_Data(cdata%blk_no)%Intdiag%ulwsfci, &
                  xmu=GFS_Interstitial(cdata%thrd_no)%xmu,xcosz=GFS_Interstitial(cdata%thrd_no)%xcosz, &
                  adjnirbmu=GFS_Interstitial(cdata%thrd_no)%adjnirbmu,adjnirdfu=GFS_Interstitial(cdata%thrd_no)%adjnirdfu, &
                  adjvisbmu=GFS_Interstitial(cdata%thrd_no)%adjvisbmu,adjvisdfu=GFS_Interstitial(cdata%thrd_no)%adjvisdfu, &
                  adjnirbmd=GFS_Interstitial(cdata%thrd_no)%adjnirbmd,adjnirdfd=GFS_Interstitial(cdata%thrd_no)%adjnirdfd, &
                  adjvisbmd=GFS_Interstitial(cdata%thrd_no)%adjvisbmd,adjvisdfd=GFS_Interstitial(cdata%thrd_no)%adjvisdfd, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in dcyc2t3_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_surface_generic_pre_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,vfrac=GFS_Data(cdata%blk_no)%Sfcprop%vfrac, &
                  islmsk=GFS_Interstitial(cdata%thrd_no)%islmsk,isot=GFS_Control%isot,ivegsrc=GFS_Control%ivegsrc, &
                  stype=GFS_Data(cdata%blk_no)%Sfcprop%stype,vtype=GFS_Data(cdata%blk_no)%Sfcprop%vtype, &
                  slope=GFS_Data(cdata%blk_no)%Sfcprop%slope,prsik_1=GFS_Data(cdata%blk_no)%Statein%prsik(:,1), &
                  prslk_1=GFS_Data(cdata%blk_no)%Statein%prslk(:,1),semis=GFS_Data(cdata%blk_no)%Radtend%semis, &
                  adjsfcdlw=GFS_Data(cdata%blk_no)%Intdiag%dlwsfci,tsfc=GFS_Data(cdata%blk_no)%Sfcprop%tsfc, &
                  phil=GFS_Data(cdata%blk_no)%Statein%phil,con_g=con_g,sigmaf=GFS_Interstitial(cdata%thrd_no)%sigmaf, &
                  soiltyp=GFS_Interstitial(cdata%thrd_no)%soiltype,vegtype=GFS_Interstitial(cdata%thrd_no)%vegtype, &
                  slopetyp=GFS_Interstitial(cdata%thrd_no)%slopetype,work3=GFS_Interstitial(cdata%thrd_no)%work3, &
                  gabsbdlw=GFS_Interstitial(cdata%thrd_no)%gabsbdlw,tsurf=GFS_Interstitial(cdata%thrd_no)%tsurf, &
                  zlvl=GFS_Data(cdata%blk_no)%Intdiag%zlvl,do_sppt=GFS_Control%do_sppt,dtdtr=GFS_Data(cdata%blk_no)%Tbd%dtdtr, &
                  drain_cpl=GFS_Data(cdata%blk_no)%Tbd%drain_cpl,dsnow_cpl=GFS_Data(cdata%blk_no)%Tbd%dsnow_cpl, &
                  rain_cpl=GFS_Data(cdata%blk_no)%Coupling%rain_cpl,snow_cpl=GFS_Data(cdata%blk_no)%Coupling%snow_cpl, &
                  do_sfcperts=GFS_Control%do_sfcperts,nsfcpert=GFS_Control%nsfcpert,sfc_wts=GFS_Data(cdata%blk_no)%Coupling%sfc_wts, &
                  pertz0=GFS_Control%pertz0,pertzt=GFS_Control%pertzt,pertshc=GFS_Control%pertshc, &
                  pertlai=GFS_Control%pertlai,pertvegf=GFS_Control%pertvegf,z01d=GFS_Interstitial(cdata%thrd_no)%z01d, &
                  zt1d=GFS_Interstitial(cdata%thrd_no)%zt1d,bexp1d=GFS_Interstitial(cdata%thrd_no)%bexp1d, &
                  xlai1d=GFS_Interstitial(cdata%thrd_no)%xlai1d,vegf1d=GFS_Interstitial(cdata%thrd_no)%vegf1d, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_surface_generic_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_suite_interstitial_2_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,lssav=GFS_Control%lssav, &
                  ldiag3d=GFS_Control%ldiag3d,lsidea=GFS_Control%lsidea,cplflx=GFS_Control%cplflx, &
                  flag_cice=GFS_Interstitial(cdata%thrd_no)%flag_cice,shal_cnv=GFS_Control%shal_cnv, &
                  old_monin=GFS_Control%old_monin,mstrat=GFS_Control%mstrat,do_shoc=GFS_Control%do_shoc, &
                  imfshalcnv=GFS_Control%imfshalcnv,dtf=GFS_Control%dtf,xcosz=GFS_Interstitial(cdata%thrd_no)%xcosz, &
                  adjsfcdsw=GFS_Data(cdata%blk_no)%Intdiag%dswsfci,adjsfcdlw=GFS_Data(cdata%blk_no)%Intdiag%dlwsfci, &
                  pgr=GFS_Data(cdata%blk_no)%Statein%pgr,ulwsfc_cice=GFS_Interstitial(cdata%thrd_no)%ulwsfc_cice, &
                  lwhd=GFS_Data(cdata%blk_no)%Radtend%lwhd,htrsw=GFS_Data(cdata%blk_no)%Radtend%htrsw, &
                  htrlw=GFS_Data(cdata%blk_no)%Radtend%htrlw,xmu=GFS_Interstitial(cdata%thrd_no)%xmu, &
                  ctei_rm=GFS_Control%ctei_rm,work1=GFS_Interstitial(cdata%thrd_no)%work1, &
                  work2=GFS_Interstitial(cdata%thrd_no)%work2,prsi=GFS_Data(cdata%blk_no)%Statein%prsi, &
                  tgrs=GFS_Data(cdata%blk_no)%Statein%tgrs,prsl=GFS_Data(cdata%blk_no)%Statein%prsl, &
                  qgrs_water_vapor=GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntqv), &
                  qgrs_cloud_water=GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntcw), &
                  cp=con_cp,hvap=con_hvap,prslk=GFS_Data(cdata%blk_no)%Statein%prslk,suntim=GFS_Data(cdata%blk_no)%Intdiag%suntim, &
                  adjsfculw=GFS_Data(cdata%blk_no)%Intdiag%ulwsfci,dlwsfc=GFS_Data(cdata%blk_no)%Intdiag%dlwsfc, &
                  ulwsfc=GFS_Data(cdata%blk_no)%Intdiag%ulwsfc,psmean=GFS_Data(cdata%blk_no)%Intdiag%psmean, &
                  dt3dt_lw=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,1),dt3dt_sw=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,2), &
                  dt3dt_pbl=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,3),dt3dt_dcnv=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,4), &
                  dt3dt_scnv=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,5),dt3dt_mp=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,6), &
                  ctei_rml=GFS_Interstitial(cdata%thrd_no)%ctei_rml,ctei_r=GFS_Interstitial(cdata%thrd_no)%ctei_r, &
                  kinver=GFS_Interstitial(cdata%thrd_no)%kinver,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_suite_interstitial_2_run"
        ierr=cdata%errflg
        return
      end if

    
      associate(cnt => cdata%loop_cnt)
      do cnt=1,2


      call sfc_ex_coef_run(im=GFS_Control%blksz(cdata%blk_no),ps=GFS_Data(cdata%blk_no)%Statein%pgr, &
                  u1=GFS_Data(cdata%blk_no)%Statein%ugrs(:,1),v1=GFS_Data(cdata%blk_no)%Statein%vgrs(:,1), &
                  t1=GFS_Data(cdata%blk_no)%Statein%tgrs(:,1),q1=GFS_Data(cdata%blk_no)%Statein%qgrs(:,1,GFS_Control%ntqv), &
                  z1=GFS_Data(cdata%blk_no)%Intdiag%zlvl,snwdph=GFS_Data(cdata%blk_no)%Sfcprop%snowd, &
                  tskin=GFS_Data(cdata%blk_no)%Sfcprop%tsfc,z0rl=GFS_Data(cdata%blk_no)%Sfcprop%zorl, &
                  cm=GFS_Interstitial(cdata%thrd_no)%cd,ch=GFS_Interstitial(cdata%thrd_no)%cdq, &
                  rb=GFS_Interstitial(cdata%thrd_no)%rb,prsl1=GFS_Data(cdata%blk_no)%Statein%prsl(:,1), &
                  prslki=GFS_Interstitial(cdata%thrd_no)%work3,islimsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  stress=GFS_Interstitial(cdata%thrd_no)%stress,fm=GFS_Data(cdata%blk_no)%Sfcprop%ffmm, &
                  fh=GFS_Data(cdata%blk_no)%Sfcprop%ffhh,ustar=GFS_Data(cdata%blk_no)%Sfcprop%uustar, &
                  wind=GFS_Interstitial(cdata%thrd_no)%wind,ddvel=GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,GFS_Control%num_p2d), &
                  fm10=GFS_Interstitial(cdata%thrd_no)%fm10,fh2=GFS_Interstitial(cdata%thrd_no)%fh2, &
                  sigmaf=GFS_Interstitial(cdata%thrd_no)%sigmaf,vegtype=GFS_Interstitial(cdata%thrd_no)%vegtype, &
                  shdmax=GFS_Data(cdata%blk_no)%Sfcprop%shdmax,ivegsrc=GFS_Control%ivegsrc, &
                  z0pert=GFS_Interstitial(cdata%thrd_no)%z01d,ztpert=GFS_Interstitial(cdata%thrd_no)%zt1d, &
                  tsurf=GFS_Interstitial(cdata%thrd_no)%tsurf,flag_iter=GFS_Interstitial(cdata%thrd_no)%flag_iter, &
                  redrag=GFS_Control%redrag,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in sfc_ex_coef_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_surface_loop_control_part1_run(im=GFS_Control%blksz(cdata%blk_no),iter=cdata%loop_cnt,wind=GFS_Interstitial(cdata%thrd_no)%wind, &
                  flag_guess=GFS_Interstitial(cdata%thrd_no)%flag_guess,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_surface_loop_control_part1_run"
        ierr=cdata%errflg
        return
      end if

    
      call sfc_nst_pre_run(im=GFS_Control%blksz(cdata%blk_no),islimsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  oro=GFS_Data(cdata%blk_no)%Sfcprop%oro,oro_uf=GFS_Data(cdata%blk_no)%Sfcprop%oro_uf, &
                  tsfc=GFS_Data(cdata%blk_no)%Sfcprop%tsfc,tsurf=GFS_Interstitial(cdata%thrd_no)%tsurf, &
                  tskin=GFS_Interstitial(cdata%thrd_no)%tseal,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in sfc_nst_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call sfc_nst_run(im=GFS_Control%blksz(cdata%blk_no),ps=GFS_Data(cdata%blk_no)%Statein%pgr, &
                  u1=GFS_Data(cdata%blk_no)%Statein%ugrs(:,1),v1=GFS_Data(cdata%blk_no)%Statein%vgrs(:,1), &
                  t1=GFS_Data(cdata%blk_no)%Statein%tgrs(:,1),q1=GFS_Data(cdata%blk_no)%Statein%qgrs(:,1,GFS_Control%ntqv), &
                  tref=GFS_Data(cdata%blk_no)%Sfcprop%tref,cm=GFS_Interstitial(cdata%thrd_no)%cd, &
                  ch=GFS_Interstitial(cdata%thrd_no)%cdq,prsl1=GFS_Data(cdata%blk_no)%Statein%prsl(:,1), &
                  prslki=GFS_Interstitial(cdata%thrd_no)%work3,islimsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  xlon=GFS_Data(cdata%blk_no)%Grid%xlon,sinlat=GFS_Data(cdata%blk_no)%Grid%sinlat, &
                  stress=GFS_Interstitial(cdata%thrd_no)%stress,sfcemis=GFS_Data(cdata%blk_no)%Radtend%semis, &
                  dlwflx=GFS_Interstitial(cdata%thrd_no)%gabsbdlw,sfcnsw=GFS_Data(cdata%blk_no)%Intdiag%nswsfci, &
                  rain=GFS_Data(cdata%blk_no)%Sfcprop%tprcp,timestep=GFS_Control%dtf,kdt=GFS_Control%kdt, &
                  solhr=GFS_Control%solhr,xcosz=GFS_Interstitial(cdata%thrd_no)%xcosz,ddvel=GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,GFS_Control%num_p2d), &
                  flag_iter=GFS_Interstitial(cdata%thrd_no)%flag_iter,flag_guess=GFS_Interstitial(cdata%thrd_no)%flag_guess, &
                  nstf_name1=GFS_Control%nstf_name(1),nstf_name4=GFS_Control%nstf_name(4), &
                  nstf_name5=GFS_Control%nstf_name(5),lprnt=GFS_Control%lprnt,ipr=GFS_Interstitial(cdata%thrd_no)%ipr, &
                  tskin=GFS_Interstitial(cdata%thrd_no)%tseal,tsurf=GFS_Interstitial(cdata%thrd_no)%tsurf, &
                  xt=GFS_Data(cdata%blk_no)%Sfcprop%xt,xs=GFS_Data(cdata%blk_no)%Sfcprop%xs, &
                  xu=GFS_Data(cdata%blk_no)%Sfcprop%xu,xv=GFS_Data(cdata%blk_no)%Sfcprop%xv, &
                  xz=GFS_Data(cdata%blk_no)%Sfcprop%xz,zm=GFS_Data(cdata%blk_no)%Sfcprop%zm, &
                  xtts=GFS_Data(cdata%blk_no)%Sfcprop%xtts,xzts=GFS_Data(cdata%blk_no)%Sfcprop%xzts, &
                  dt_cool=GFS_Data(cdata%blk_no)%Sfcprop%dt_cool,z_c=GFS_Data(cdata%blk_no)%Sfcprop%z_c, &
                  c_0=GFS_Data(cdata%blk_no)%Sfcprop%c_0,c_d=GFS_Data(cdata%blk_no)%Sfcprop%c_d, &
                  w_0=GFS_Data(cdata%blk_no)%Sfcprop%w_0,w_d=GFS_Data(cdata%blk_no)%Sfcprop%w_d, &
                  d_conv=GFS_Data(cdata%blk_no)%Sfcprop%d_conv,ifd=GFS_Data(cdata%blk_no)%Sfcprop%ifd, &
                  qrain=GFS_Data(cdata%blk_no)%Sfcprop%qrain,qsurf=GFS_Interstitial(cdata%thrd_no)%qss, &
                  gflux=GFS_Interstitial(cdata%thrd_no)%gflx,cmm=GFS_Data(cdata%blk_no)%Intdiag%cmm, &
                  chh=GFS_Data(cdata%blk_no)%Intdiag%chh,evap=GFS_Interstitial(cdata%thrd_no)%evap, &
                  hflx=GFS_Interstitial(cdata%thrd_no)%hflx,ep=GFS_Interstitial(cdata%thrd_no)%ep1d, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in sfc_nst_run"
        ierr=cdata%errflg
        return
      end if

    
      call sfc_nst_post_run(im=GFS_Control%blksz(cdata%blk_no),islimsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  oro=GFS_Data(cdata%blk_no)%Sfcprop%oro,oro_uf=GFS_Data(cdata%blk_no)%Sfcprop%oro_uf, &
                  nstf_name1=GFS_Control%nstf_name(1),nstf_name4=GFS_Control%nstf_name(4), &
                  nstf_name5=GFS_Control%nstf_name(5),xt=GFS_Data(cdata%blk_no)%Sfcprop%xt, &
                  xz=GFS_Data(cdata%blk_no)%Sfcprop%xz,dt_cool=GFS_Data(cdata%blk_no)%Sfcprop%dt_cool, &
                  z_c=GFS_Data(cdata%blk_no)%Sfcprop%z_c,rslimsk=GFS_Data(cdata%blk_no)%Sfcprop%slmsk, &
                  tref=GFS_Data(cdata%blk_no)%Sfcprop%tref,xlon=GFS_Data(cdata%blk_no)%Grid%xlon, &
                  tsurf=GFS_Interstitial(cdata%thrd_no)%tsurf,dtzm=GFS_Interstitial(cdata%thrd_no)%dtzm, &
                  tsfc=GFS_Data(cdata%blk_no)%Sfcprop%tsfc,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in sfc_nst_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call lsm_noah_run(im=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%lsoil,ps=GFS_Data(cdata%blk_no)%Statein%pgr, &
                  u1=GFS_Data(cdata%blk_no)%Statein%ugrs(:,1),v1=GFS_Data(cdata%blk_no)%Statein%vgrs(:,1), &
                  t1=GFS_Data(cdata%blk_no)%Statein%tgrs(:,1),q1=GFS_Data(cdata%blk_no)%Statein%qgrs(:,1,GFS_Control%ntqv), &
                  soiltyp=GFS_Interstitial(cdata%thrd_no)%soiltype,vegtype=GFS_Interstitial(cdata%thrd_no)%vegtype, &
                  sigmaf=GFS_Interstitial(cdata%thrd_no)%sigmaf,sfcemis=GFS_Data(cdata%blk_no)%Radtend%semis, &
                  dlwflx=GFS_Interstitial(cdata%thrd_no)%gabsbdlw,dswsfc=GFS_Data(cdata%blk_no)%Intdiag%dswsfci, &
                  snet=GFS_Data(cdata%blk_no)%Intdiag%nswsfci,delt=GFS_Control%dtf,tg3=GFS_Data(cdata%blk_no)%Sfcprop%tg3, &
                  cm=GFS_Interstitial(cdata%thrd_no)%cd,ch=GFS_Interstitial(cdata%thrd_no)%cdq, &
                  prsl1=GFS_Data(cdata%blk_no)%Statein%prsl(:,1),prslki=GFS_Interstitial(cdata%thrd_no)%work3, &
                  zf=GFS_Data(cdata%blk_no)%Intdiag%zlvl,islimsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  ddvel=GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,GFS_Control%num_p2d),slopetyp=GFS_Interstitial(cdata%thrd_no)%slopetype, &
                  shdmin=GFS_Data(cdata%blk_no)%Sfcprop%shdmin,shdmax=GFS_Data(cdata%blk_no)%Sfcprop%shdmax, &
                  snoalb=GFS_Data(cdata%blk_no)%Sfcprop%snoalb,sfalb=GFS_Data(cdata%blk_no)%Radtend%sfalb, &
                  flag_iter=GFS_Interstitial(cdata%thrd_no)%flag_iter,flag_guess=GFS_Interstitial(cdata%thrd_no)%flag_guess, &
                  isot=GFS_Control%isot,ivegsrc=GFS_Control%ivegsrc,bexppert=GFS_Interstitial(cdata%thrd_no)%bexp1d, &
                  xlaipert=GFS_Interstitial(cdata%thrd_no)%xlai1d,vegfpert=GFS_Interstitial(cdata%thrd_no)%vegf1d, &
                  pertvegf=GFS_Control%pertvegf,weasd=GFS_Data(cdata%blk_no)%Sfcprop%weasd, &
                  snwdph=GFS_Data(cdata%blk_no)%Sfcprop%snowd,tskin=GFS_Data(cdata%blk_no)%Sfcprop%tsfc, &
                  tprcp=GFS_Data(cdata%blk_no)%Sfcprop%tprcp,srflag=GFS_Data(cdata%blk_no)%Sfcprop%srflag, &
                  smc=GFS_Data(cdata%blk_no)%Sfcprop%smc,stc=GFS_Data(cdata%blk_no)%Sfcprop%stc, &
                  slc=GFS_Data(cdata%blk_no)%Sfcprop%slc,canopy=GFS_Data(cdata%blk_no)%Sfcprop%canopy, &
                  trans=GFS_Interstitial(cdata%thrd_no)%trans,tsurf=GFS_Interstitial(cdata%thrd_no)%tsurf, &
                  zorl=GFS_Data(cdata%blk_no)%Sfcprop%zorl,sncovr1=GFS_Data(cdata%blk_no)%Sfcprop%sncovr, &
                  qsurf=GFS_Interstitial(cdata%thrd_no)%qss,gflux=GFS_Interstitial(cdata%thrd_no)%gflx, &
                  drain=GFS_Interstitial(cdata%thrd_no)%drain,evap=GFS_Interstitial(cdata%thrd_no)%evap, &
                  hflx=GFS_Interstitial(cdata%thrd_no)%hflx,ep=GFS_Interstitial(cdata%thrd_no)%ep1d, &
                  runoff=GFS_Interstitial(cdata%thrd_no)%runoff,cmm=GFS_Data(cdata%blk_no)%Intdiag%cmm, &
                  chh=GFS_Data(cdata%blk_no)%Intdiag%chh,evbs=GFS_Interstitial(cdata%thrd_no)%evbs, &
                  evcw=GFS_Interstitial(cdata%thrd_no)%evcw,sbsno=GFS_Interstitial(cdata%thrd_no)%sbsno, &
                  snowc=GFS_Interstitial(cdata%thrd_no)%snowc,stm=GFS_Data(cdata%blk_no)%Intdiag%soilm, &
                  snohf=GFS_Interstitial(cdata%thrd_no)%snohf,smcwlt2=GFS_Data(cdata%blk_no)%Intdiag%smcwlt2, &
                  smcref2=GFS_Data(cdata%blk_no)%Intdiag%smcref2,wet1=GFS_Data(cdata%blk_no)%Sfcprop%wet1, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in lsm_noah_run"
        ierr=cdata%errflg
        return
      end if

    
      call sfc_sice_run(im=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%lsoil,ps=GFS_Data(cdata%blk_no)%Statein%pgr, &
                  u1=GFS_Data(cdata%blk_no)%Statein%ugrs(:,1),v1=GFS_Data(cdata%blk_no)%Statein%vgrs(:,1), &
                  t1=GFS_Data(cdata%blk_no)%Statein%tgrs(:,1),q1=GFS_Data(cdata%blk_no)%Statein%qgrs(:,1,GFS_Control%ntqv), &
                  delt=GFS_Control%dtf,sfcemis=GFS_Data(cdata%blk_no)%Radtend%semis,dlwflx=GFS_Interstitial(cdata%thrd_no)%gabsbdlw, &
                  sfcnsw=GFS_Data(cdata%blk_no)%Intdiag%nswsfci,sfcdsw=GFS_Data(cdata%blk_no)%Intdiag%dswsfci, &
                  srflag=GFS_Data(cdata%blk_no)%Sfcprop%srflag,cm=GFS_Interstitial(cdata%thrd_no)%cd, &
                  ch=GFS_Interstitial(cdata%thrd_no)%cdq,prsl1=GFS_Data(cdata%blk_no)%Statein%prsl(:,1), &
                  prslki=GFS_Interstitial(cdata%thrd_no)%work3,islimsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  ddvel=GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,GFS_Control%num_p2d),flag_iter=GFS_Interstitial(cdata%thrd_no)%flag_iter, &
                  mom4ice=GFS_Control%mom4ice,lsm=GFS_Control%lsm,lprnt=GFS_Control%lprnt, &
                  ipr=GFS_Interstitial(cdata%thrd_no)%ipr,hice=GFS_Data(cdata%blk_no)%Sfcprop%hice, &
                  fice=GFS_Data(cdata%blk_no)%Sfcprop%fice,tice=GFS_Data(cdata%blk_no)%Sfcprop%tisfc, &
                  weasd=GFS_Data(cdata%blk_no)%Sfcprop%weasd,tskin=GFS_Data(cdata%blk_no)%Sfcprop%tsfc, &
                  tprcp=GFS_Data(cdata%blk_no)%Sfcprop%tprcp,stc=GFS_Data(cdata%blk_no)%Sfcprop%stc, &
                  ep=GFS_Interstitial(cdata%thrd_no)%ep1d,snwdph=GFS_Data(cdata%blk_no)%Sfcprop%snowd, &
                  qsurf=GFS_Interstitial(cdata%thrd_no)%qss,snowmt=GFS_Interstitial(cdata%thrd_no)%snowmt, &
                  gflux=GFS_Interstitial(cdata%thrd_no)%gflx,cmm=GFS_Data(cdata%blk_no)%Intdiag%cmm, &
                  chh=GFS_Data(cdata%blk_no)%Intdiag%chh,evap=GFS_Interstitial(cdata%thrd_no)%evap, &
                  hflx=GFS_Interstitial(cdata%thrd_no)%hflx,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in sfc_sice_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_surface_loop_control_part2_run(im=GFS_Control%blksz(cdata%blk_no),iter=cdata%loop_cnt,wind=GFS_Interstitial(cdata%thrd_no)%wind, &
                  flag_guess=GFS_Interstitial(cdata%thrd_no)%flag_guess,flag_iter=GFS_Interstitial(cdata%thrd_no)%flag_iter, &
                  islmsk=GFS_Interstitial(cdata%thrd_no)%islmsk,nstf_name1=GFS_Control%nstf_name(1), &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_surface_loop_control_part2_run"
        ierr=cdata%errflg
        return
      end if

    
      end do
      end associate

      call dcyc2t3_post_run(im=GFS_Control%blksz(cdata%blk_no),adjsfcdsw=GFS_Data(cdata%blk_no)%Intdiag%dswsfci, &
                  adjsfcnsw=GFS_Data(cdata%blk_no)%Intdiag%nswsfci,adjsfcusw=GFS_Data(cdata%blk_no)%Intdiag%uswsfci, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in dcyc2t3_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call sfc_diag_run(im=GFS_Control%blksz(cdata%blk_no),grav=con_g,cp=con_cp,eps=con_eps,epsm1=con_epsm1, &
                  ps=GFS_Data(cdata%blk_no)%Statein%pgr,u1=GFS_Data(cdata%blk_no)%Stateout%gu0(:,1), &
                  v1=GFS_Data(cdata%blk_no)%Stateout%gv0(:,1),t1=GFS_Data(cdata%blk_no)%Stateout%gt0(:,1), &
                  q1=GFS_Data(cdata%blk_no)%Stateout%gq0(:,1,GFS_Control%ntqv),tskin=GFS_Data(cdata%blk_no)%Sfcprop%tsfc, &
                  qsurf=GFS_Interstitial(cdata%thrd_no)%qss,f10m=GFS_Data(cdata%blk_no)%Sfcprop%f10m, &
                  u10m=GFS_Data(cdata%blk_no)%Intdiag%u10m,v10m=GFS_Data(cdata%blk_no)%Intdiag%v10m, &
                  t2m=GFS_Data(cdata%blk_no)%Sfcprop%t2m,q2m=GFS_Data(cdata%blk_no)%Sfcprop%q2m, &
                  prslki=GFS_Interstitial(cdata%thrd_no)%work3,evap=GFS_Interstitial(cdata%thrd_no)%evap, &
                  fm=GFS_Data(cdata%blk_no)%Sfcprop%ffmm,fh=GFS_Data(cdata%blk_no)%Sfcprop%ffhh, &
                  fm10=GFS_Interstitial(cdata%thrd_no)%fm10,fh2=GFS_Interstitial(cdata%thrd_no)%fh2, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in sfc_diag_run"
        ierr=cdata%errflg
        return
      end if

    
      call sfc_diag_post_run(im=GFS_Control%blksz(cdata%blk_no),lssav=GFS_Control%lssav,dtf=GFS_Control%dtf, &
                  con_eps=con_eps,con_epsm1=con_epsm1,pgr=GFS_Data(cdata%blk_no)%Statein%pgr, &
                  t2m=GFS_Data(cdata%blk_no)%Sfcprop%t2m,q2m=GFS_Data(cdata%blk_no)%Sfcprop%q2m, &
                  u10m=GFS_Data(cdata%blk_no)%Intdiag%u10m,v10m=GFS_Data(cdata%blk_no)%Intdiag%v10m, &
                  tmpmin=GFS_Data(cdata%blk_no)%Intdiag%tmpmin,tmpmax=GFS_Data(cdata%blk_no)%Intdiag%tmpmax, &
                  spfhmin=GFS_Data(cdata%blk_no)%Intdiag%spfhmin,spfhmax=GFS_Data(cdata%blk_no)%Intdiag%spfhmax, &
                  wind10mmax=GFS_Data(cdata%blk_no)%Intdiag%wind10mmax,u10mmax=GFS_Data(cdata%blk_no)%Intdiag%u10mmax, &
                  v10mmax=GFS_Data(cdata%blk_no)%Intdiag%v10mmax,dpt2m=GFS_Data(cdata%blk_no)%Intdiag%dpt2m, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in sfc_diag_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_surface_generic_post_run(im=GFS_Control%blksz(cdata%blk_no),cplflx=GFS_Control%cplflx,lssav=GFS_Control%lssav, &
                  islmsk=GFS_Interstitial(cdata%thrd_no)%islmsk,dtf=GFS_Control%dtf,ep1d=GFS_Interstitial(cdata%thrd_no)%ep1d, &
                  gflx=GFS_Interstitial(cdata%thrd_no)%gflx,tgrs_1=GFS_Data(cdata%blk_no)%Statein%tgrs(:,1), &
                  qgrs_1=GFS_Data(cdata%blk_no)%Statein%qgrs(:,1,GFS_Control%ntqv),ugrs_1=GFS_Data(cdata%blk_no)%Statein%ugrs(:,1), &
                  vgrs_1=GFS_Data(cdata%blk_no)%Statein%vgrs(:,1),adjsfcdlw=GFS_Data(cdata%blk_no)%Intdiag%dlwsfci, &
                  adjsfcdsw=GFS_Data(cdata%blk_no)%Intdiag%dswsfci,adjnirbmd=GFS_Interstitial(cdata%thrd_no)%adjnirbmd, &
                  adjnirdfd=GFS_Interstitial(cdata%thrd_no)%adjnirdfd,adjvisbmd=GFS_Interstitial(cdata%thrd_no)%adjvisbmd, &
                  adjvisdfd=GFS_Interstitial(cdata%thrd_no)%adjvisdfd,adjsfculw=GFS_Data(cdata%blk_no)%Intdiag%ulwsfci, &
                  adjnirbmu=GFS_Interstitial(cdata%thrd_no)%adjnirbmu,adjnirdfu=GFS_Interstitial(cdata%thrd_no)%adjnirdfu, &
                  adjvisbmu=GFS_Interstitial(cdata%thrd_no)%adjvisbmu,adjvisdfu=GFS_Interstitial(cdata%thrd_no)%adjvisdfu, &
                  t2m=GFS_Data(cdata%blk_no)%Sfcprop%t2m,q2m=GFS_Data(cdata%blk_no)%Sfcprop%q2m, &
                  u10m=GFS_Data(cdata%blk_no)%Intdiag%u10m,v10m=GFS_Data(cdata%blk_no)%Intdiag%v10m, &
                  tsfc=GFS_Data(cdata%blk_no)%Sfcprop%tsfc,pgr=GFS_Data(cdata%blk_no)%Statein%pgr, &
                  xcosz=GFS_Interstitial(cdata%thrd_no)%xcosz,evbs=GFS_Interstitial(cdata%thrd_no)%evbs, &
                  evcw=GFS_Interstitial(cdata%thrd_no)%evcw,trans=GFS_Interstitial(cdata%thrd_no)%trans, &
                  sbsno=GFS_Interstitial(cdata%thrd_no)%sbsno,snowc=GFS_Interstitial(cdata%thrd_no)%snowc, &
                  snohf=GFS_Interstitial(cdata%thrd_no)%snohf,epi=GFS_Data(cdata%blk_no)%Intdiag%epi, &
                  gfluxi=GFS_Data(cdata%blk_no)%Intdiag%gfluxi,t1=GFS_Data(cdata%blk_no)%Intdiag%t1, &
                  q1=GFS_Data(cdata%blk_no)%Intdiag%q1,u1=GFS_Data(cdata%blk_no)%Intdiag%u1, &
                  v1=GFS_Data(cdata%blk_no)%Intdiag%v1,dlwsfci_cpl=GFS_Data(cdata%blk_no)%Coupling%dlwsfci_cpl, &
                  dswsfci_cpl=GFS_Data(cdata%blk_no)%Coupling%dswsfci_cpl,dlwsfc_cpl=GFS_Data(cdata%blk_no)%Coupling%dlwsfc_cpl, &
                  dswsfc_cpl=GFS_Data(cdata%blk_no)%Coupling%dswsfc_cpl,dnirbmi_cpl=GFS_Data(cdata%blk_no)%Coupling%dnirbmi_cpl, &
                  dnirdfi_cpl=GFS_Data(cdata%blk_no)%Coupling%dnirdfi_cpl,dvisbmi_cpl=GFS_Data(cdata%blk_no)%Coupling%dvisbmi_cpl, &
                  dvisdfi_cpl=GFS_Data(cdata%blk_no)%Coupling%dvisdfi_cpl,dnirbm_cpl=GFS_Data(cdata%blk_no)%Coupling%dnirbm_cpl, &
                  dnirdf_cpl=GFS_Data(cdata%blk_no)%Coupling%dnirdf_cpl,dvisbm_cpl=GFS_Data(cdata%blk_no)%Coupling%dvisbm_cpl, &
                  dvisdf_cpl=GFS_Data(cdata%blk_no)%Coupling%dvisdf_cpl,nlwsfci_cpl=GFS_Data(cdata%blk_no)%Coupling%nlwsfci_cpl, &
                  nlwsfc_cpl=GFS_Data(cdata%blk_no)%Coupling%nlwsfc_cpl,t2mi_cpl=GFS_Data(cdata%blk_no)%Coupling%t2mi_cpl, &
                  q2mi_cpl=GFS_Data(cdata%blk_no)%Coupling%q2mi_cpl,u10mi_cpl=GFS_Data(cdata%blk_no)%Coupling%u10mi_cpl, &
                  v10mi_cpl=GFS_Data(cdata%blk_no)%Coupling%v10mi_cpl,tsfci_cpl=GFS_Data(cdata%blk_no)%Coupling%tsfci_cpl, &
                  psurfi_cpl=GFS_Data(cdata%blk_no)%Coupling%psurfi_cpl,nnirbmi_cpl=GFS_Data(cdata%blk_no)%Coupling%nnirbmi_cpl, &
                  nnirdfi_cpl=GFS_Data(cdata%blk_no)%Coupling%nnirdfi_cpl,nvisbmi_cpl=GFS_Data(cdata%blk_no)%Coupling%nvisbmi_cpl, &
                  nvisdfi_cpl=GFS_Data(cdata%blk_no)%Coupling%nvisdfi_cpl,nswsfci_cpl=GFS_Data(cdata%blk_no)%Coupling%nswsfci_cpl, &
                  nswsfc_cpl=GFS_Data(cdata%blk_no)%Coupling%nswsfc_cpl,nnirbm_cpl=GFS_Data(cdata%blk_no)%Coupling%nnirbm_cpl, &
                  nnirdf_cpl=GFS_Data(cdata%blk_no)%Coupling%nnirdf_cpl,nvisbm_cpl=GFS_Data(cdata%blk_no)%Coupling%nvisbm_cpl, &
                  nvisdf_cpl=GFS_Data(cdata%blk_no)%Coupling%nvisdf_cpl,gflux=GFS_Data(cdata%blk_no)%Intdiag%gflux, &
                  evbsa=GFS_Data(cdata%blk_no)%Intdiag%evbsa,evcwa=GFS_Data(cdata%blk_no)%Intdiag%evcwa, &
                  transa=GFS_Data(cdata%blk_no)%Intdiag%transa,sbsnoa=GFS_Data(cdata%blk_no)%Intdiag%sbsnoa, &
                  snowca=GFS_Data(cdata%blk_no)%Intdiag%snowca,snohfa=GFS_Data(cdata%blk_no)%Intdiag%snohfa, &
                  ep=GFS_Data(cdata%blk_no)%Intdiag%ep,runoff=GFS_Data(cdata%blk_no)%Intdiag%runoff, &
                  srunoff=GFS_Data(cdata%blk_no)%Intdiag%srunoff,runof=GFS_Interstitial(cdata%thrd_no)%runoff, &
                  drain=GFS_Interstitial(cdata%thrd_no)%drain,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_surface_generic_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_PBL_generic_pre_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,nvdiff=GFS_Interstitial(cdata%thrd_no)%nvdiff, &
                  ntrac=GFS_Control%ntrac,ntqv=GFS_Control%ntqv,ntcw=GFS_Control%ntcw,ntiw=GFS_Control%ntiw, &
                  ntrw=GFS_Control%ntrw,ntsw=GFS_Control%ntsw,ntlnc=GFS_Control%ntlnc,ntinc=GFS_Control%ntinc, &
                  ntwa=GFS_Control%ntwa,ntia=GFS_Control%ntia,ntgl=GFS_Control%ntgl,ntoz=GFS_Control%ntoz, &
                  ntke=GFS_Control%ntke,ntkev=GFS_Interstitial(cdata%thrd_no)%ntkev,imp_physics=GFS_Control%imp_physics, &
                  imp_physics_gfdl=GFS_Control%imp_physics_gfdl,imp_physics_thompson=GFS_Control%imp_physics_thompson, &
                  imp_physics_wsm6=GFS_Control%imp_physics_wsm6,ltaerosol=GFS_Control%ltaerosol, &
                  satmedmf=GFS_Control%satmedmf,qgrs=GFS_Data(cdata%blk_no)%Statein%qgrs, &
                  vdftra=GFS_Interstitial(cdata%thrd_no)%vdftra,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_PBL_generic_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call hedmf_run(ix=GFS_Control%blksz(cdata%blk_no),im=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%levs, &
                  ntrac=GFS_Interstitial(cdata%thrd_no)%nvdiff,ntcw=GFS_Control%ntcw,dv=GFS_Interstitial(cdata%thrd_no)%dvdt, &
                  du=GFS_Interstitial(cdata%thrd_no)%dudt,tau=GFS_Interstitial(cdata%thrd_no)%dtdt, &
                  rtg=GFS_Interstitial(cdata%thrd_no)%dvdftra,u1=GFS_Data(cdata%blk_no)%Statein%ugrs, &
                  v1=GFS_Data(cdata%blk_no)%Statein%vgrs,t1=GFS_Data(cdata%blk_no)%Statein%tgrs, &
                  q1=GFS_Interstitial(cdata%thrd_no)%vdftra,swh=GFS_Data(cdata%blk_no)%Tbd%htswc, &
                  hlw=GFS_Data(cdata%blk_no)%Tbd%htlwc,xmu=GFS_Interstitial(cdata%thrd_no)%xmu, &
                  psk=GFS_Data(cdata%blk_no)%Statein%prsik(:,1),rbsoil=GFS_Interstitial(cdata%thrd_no)%rb, &
                  zorl=GFS_Data(cdata%blk_no)%Sfcprop%zorl,u10m=GFS_Data(cdata%blk_no)%Intdiag%u10m, &
                  v10m=GFS_Data(cdata%blk_no)%Intdiag%v10m,fm=GFS_Data(cdata%blk_no)%Sfcprop%ffmm, &
                  fh=GFS_Data(cdata%blk_no)%Sfcprop%ffhh,tsea=GFS_Data(cdata%blk_no)%Sfcprop%tsfc, &
                  heat=GFS_Interstitial(cdata%thrd_no)%hflx,evap=GFS_Interstitial(cdata%thrd_no)%evap, &
                  stress=GFS_Interstitial(cdata%thrd_no)%stress,spd1=GFS_Interstitial(cdata%thrd_no)%wind, &
                  kpbl=GFS_Interstitial(cdata%thrd_no)%kpbl,prsi=GFS_Data(cdata%blk_no)%Statein%prsi, &
                  del=GFS_Interstitial(cdata%thrd_no)%del,prsl=GFS_Data(cdata%blk_no)%Statein%prsl, &
                  prslk=GFS_Data(cdata%blk_no)%Statein%prslk,phii=GFS_Data(cdata%blk_no)%Statein%phii, &
                  phil=GFS_Data(cdata%blk_no)%Statein%phil,delt=GFS_Control%dtp,dspheat=GFS_Control%dspheat, &
                  dusfc=GFS_Interstitial(cdata%thrd_no)%dusfc1,dvsfc=GFS_Interstitial(cdata%thrd_no)%dvsfc1, &
                  dtsfc=GFS_Interstitial(cdata%thrd_no)%dtsfc1,dqsfc=GFS_Interstitial(cdata%thrd_no)%dqsfc1, &
                  hpbl=GFS_Data(cdata%blk_no)%Intdiag%hpbl,hgamt=GFS_Interstitial(cdata%thrd_no)%gamt, &
                  hgamq=GFS_Interstitial(cdata%thrd_no)%gamq,dkt=GFS_Interstitial(cdata%thrd_no)%dkt, &
                  kinver=GFS_Interstitial(cdata%thrd_no)%kinver,xkzm_m=GFS_Control%xkzm_m, &
                  xkzm_h=GFS_Control%xkzm_h,xkzm_s=GFS_Control%xkzm_s,lprnt=GFS_Control%lprnt, &
                  ipr=GFS_Interstitial(cdata%thrd_no)%ipr,xkzminv=GFS_Control%xkzminv,moninq_fac=GFS_Control%moninq_fac, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in hedmf_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_PBL_generic_post_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,nvdiff=GFS_Interstitial(cdata%thrd_no)%nvdiff, &
                  ntrac=GFS_Control%ntrac,ntqv=GFS_Control%ntqv,ntcw=GFS_Control%ntcw,ntiw=GFS_Control%ntiw, &
                  ntrw=GFS_Control%ntrw,ntsw=GFS_Control%ntsw,ntlnc=GFS_Control%ntlnc,ntinc=GFS_Control%ntinc, &
                  ntwa=GFS_Control%ntwa,ntia=GFS_Control%ntia,ntgl=GFS_Control%ntgl,ntoz=GFS_Control%ntoz, &
                  ntke=GFS_Control%ntke,ntkev=GFS_Interstitial(cdata%thrd_no)%ntkev,imp_physics=GFS_Control%imp_physics, &
                  imp_physics_gfdl=GFS_Control%imp_physics_gfdl,imp_physics_thompson=GFS_Control%imp_physics_thompson, &
                  imp_physics_wsm6=GFS_Control%imp_physics_wsm6,ltaerosol=GFS_Control%ltaerosol, &
                  cplflx=GFS_Control%cplflx,lssav=GFS_Control%lssav,ldiag3d=GFS_Control%ldiag3d, &
                  lsidea=GFS_Control%lsidea,hybedmf=GFS_Control%hybedmf,do_shoc=GFS_Control%do_shoc, &
                  satmedmf=GFS_Control%satmedmf,dvdftra=GFS_Interstitial(cdata%thrd_no)%dvdftra, &
                  dusfc1=GFS_Interstitial(cdata%thrd_no)%dusfc1,dvsfc1=GFS_Interstitial(cdata%thrd_no)%dvsfc1, &
                  dtsfc1=GFS_Interstitial(cdata%thrd_no)%dtsfc1,dqsfc1=GFS_Interstitial(cdata%thrd_no)%dqsfc1, &
                  dtf=GFS_Control%dtf,dudt=GFS_Interstitial(cdata%thrd_no)%dudt,dvdt=GFS_Interstitial(cdata%thrd_no)%dvdt, &
                  dtdt=GFS_Interstitial(cdata%thrd_no)%dtdt,htrsw=GFS_Data(cdata%blk_no)%Radtend%htrsw, &
                  htrlw=GFS_Data(cdata%blk_no)%Radtend%htrlw,xmu=GFS_Interstitial(cdata%thrd_no)%xmu, &
                  dqdt=GFS_Interstitial(cdata%thrd_no)%dqdt,dusfc_cpl=GFS_Data(cdata%blk_no)%Coupling%dusfc_cpl, &
                  dvsfc_cpl=GFS_Data(cdata%blk_no)%Coupling%dvsfc_cpl,dtsfc_cpl=GFS_Data(cdata%blk_no)%Coupling%dtsfc_cpl, &
                  dqsfc_cpl=GFS_Data(cdata%blk_no)%Coupling%dqsfc_cpl,dusfci_cpl=GFS_Data(cdata%blk_no)%Coupling%dusfci_cpl, &
                  dvsfci_cpl=GFS_Data(cdata%blk_no)%Coupling%dvsfci_cpl,dtsfci_cpl=GFS_Data(cdata%blk_no)%Coupling%dtsfci_cpl, &
                  dqsfci_cpl=GFS_Data(cdata%blk_no)%Coupling%dqsfci_cpl,dusfc_diag=GFS_Data(cdata%blk_no)%Intdiag%dusfc, &
                  dvsfc_diag=GFS_Data(cdata%blk_no)%Intdiag%dvsfc,dtsfc_diag=GFS_Data(cdata%blk_no)%Intdiag%dtsfc, &
                  dqsfc_diag=GFS_Data(cdata%blk_no)%Intdiag%dqsfc,dusfci_diag=GFS_Data(cdata%blk_no)%Intdiag%dusfci, &
                  dvsfci_diag=GFS_Data(cdata%blk_no)%Intdiag%dvsfci,dtsfci_diag=GFS_Data(cdata%blk_no)%Intdiag%dtsfci, &
                  dqsfci_diag=GFS_Data(cdata%blk_no)%Intdiag%dqsfci,dt3dt=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,3), &
                  du3dt_PBL=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,1),du3dt_OGWD=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,2), &
                  dv3dt_PBL=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,1),dv3dt_OGWD=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,2), &
                  dq3dt=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,1),dq3dt_ozone=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,5), &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_PBL_generic_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call gwdps_pre_run(im=GFS_Control%blksz(cdata%blk_no),nmtvr=GFS_Control%nmtvr,mntvar=GFS_Data(cdata%blk_no)%Sfcprop%hprime, &
                  hprime=GFS_Interstitial(cdata%thrd_no)%hprime1,oc=GFS_Interstitial(cdata%thrd_no)%oc, &
                  oa4=GFS_Interstitial(cdata%thrd_no)%oa4,clx=GFS_Interstitial(cdata%thrd_no)%clx, &
                  theta=GFS_Interstitial(cdata%thrd_no)%theta,sigma=GFS_Interstitial(cdata%thrd_no)%sigma, &
                  gamma=GFS_Interstitial(cdata%thrd_no)%gamma,elvmax=GFS_Interstitial(cdata%thrd_no)%elvmax, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in gwdps_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call gwdps_run(im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%levs, &
                  A=GFS_Interstitial(cdata%thrd_no)%dvdt,B=GFS_Interstitial(cdata%thrd_no)%dudt, &
                  C=GFS_Interstitial(cdata%thrd_no)%dtdt,u1=GFS_Data(cdata%blk_no)%Statein%ugrs, &
                  v1=GFS_Data(cdata%blk_no)%Statein%vgrs,t1=GFS_Data(cdata%blk_no)%Statein%tgrs, &
                  q1=GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntqv),kpbl=GFS_Interstitial(cdata%thrd_no)%kpbl, &
                  prsi=GFS_Data(cdata%blk_no)%Statein%prsi,del=GFS_Interstitial(cdata%thrd_no)%del, &
                  prsl=GFS_Data(cdata%blk_no)%Statein%prsl,prslk=GFS_Data(cdata%blk_no)%Statein%prslk, &
                  phii=GFS_Data(cdata%blk_no)%Statein%phii,phil=GFS_Data(cdata%blk_no)%Statein%phil, &
                  deltim=GFS_Control%dtp,kdt=GFS_Control%kdt,hprime=GFS_Interstitial(cdata%thrd_no)%hprime1, &
                  oc=GFS_Interstitial(cdata%thrd_no)%oc,oa4=GFS_Interstitial(cdata%thrd_no)%oa4, &
                  clx4=GFS_Interstitial(cdata%thrd_no)%clx,theta=GFS_Interstitial(cdata%thrd_no)%theta, &
                  sigma=GFS_Interstitial(cdata%thrd_no)%sigma,gamma=GFS_Interstitial(cdata%thrd_no)%gamma, &
                  elvmax=GFS_Interstitial(cdata%thrd_no)%elvmax,dusfc=GFS_Interstitial(cdata%thrd_no)%dusfcg, &
                  dvsfc=GFS_Interstitial(cdata%thrd_no)%dvsfcg,g=con_g,cp=con_cp,rd=con_rd, &
                  rv=con_rv,imx=GFS_Control%lonr,nmtvr=GFS_Control%nmtvr,cdmbgwd=GFS_Control%cdmbgwd, &
                  me=GFS_Control%me,lprnt=GFS_Control%lprnt,ipr=GFS_Interstitial(cdata%thrd_no)%ipr, &
                  rdxzb=GFS_Data(cdata%blk_no)%Intdiag%zmtnblck,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in gwdps_run"
        ierr=cdata%errflg
        return
      end if

    
      call gwdps_post_run(lssav=GFS_Control%lssav,ldiag3d=GFS_Control%ldiag3d,dtf=GFS_Control%dtf, &
                  dusfcg=GFS_Interstitial(cdata%thrd_no)%dusfcg,dvsfcg=GFS_Interstitial(cdata%thrd_no)%dvsfcg, &
                  dudt=GFS_Interstitial(cdata%thrd_no)%dudt,dvdt=GFS_Interstitial(cdata%thrd_no)%dvdt, &
                  dtdt=GFS_Interstitial(cdata%thrd_no)%dtdt,dugwd=GFS_Data(cdata%blk_no)%Intdiag%dugwd, &
                  dvgwd=GFS_Data(cdata%blk_no)%Intdiag%dvgwd,du3dt=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,2), &
                  dv3dt=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,2),dt3dt=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,2), &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in gwdps_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call rayleigh_damp_run(lsidea=GFS_Control%lsidea,im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no), &
                  km=GFS_Control%levs,A=GFS_Interstitial(cdata%thrd_no)%dvdt,B=GFS_Interstitial(cdata%thrd_no)%dudt, &
                  C=GFS_Interstitial(cdata%thrd_no)%dtdt,u1=GFS_Data(cdata%blk_no)%Statein%ugrs, &
                  v1=GFS_Data(cdata%blk_no)%Statein%vgrs,dt=GFS_Control%dtp,cp=con_cp,levr=GFS_Control%levr, &
                  pgr=GFS_Data(cdata%blk_no)%Statein%pgr,prsl=GFS_Data(cdata%blk_no)%Statein%prsl, &
                  prslrd0=GFS_Control%prslrd0,ral_ts=GFS_Control%ral_ts,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in rayleigh_damp_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_suite_stateout_update_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,ntrac=GFS_Control%ntrac, &
                  dtp=GFS_Control%dtp,tgrs=GFS_Data(cdata%blk_no)%Statein%tgrs,ugrs=GFS_Data(cdata%blk_no)%Statein%ugrs, &
                  vgrs=GFS_Data(cdata%blk_no)%Statein%vgrs,qgrs=GFS_Data(cdata%blk_no)%Statein%qgrs, &
                  dudt=GFS_Interstitial(cdata%thrd_no)%dudt,dvdt=GFS_Interstitial(cdata%thrd_no)%dvdt, &
                  dtdt=GFS_Interstitial(cdata%thrd_no)%dtdt,dqdt=GFS_Interstitial(cdata%thrd_no)%dqdt, &
                  gt0=GFS_Data(cdata%blk_no)%Stateout%gt0,gu0=GFS_Data(cdata%blk_no)%Stateout%gu0, &
                  gv0=GFS_Data(cdata%blk_no)%Stateout%gv0,gq0=GFS_Data(cdata%blk_no)%Stateout%gq0, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_suite_stateout_update_run"
        ierr=cdata%errflg
        return
      end if

    
      call ozphys_run(ix=GFS_Control%blksz(cdata%blk_no),im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs, &
                  ko3=GFS_Interstitial(cdata%thrd_no)%levozp,dt=GFS_Control%dtp,oz=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntoz), &
                  tin=GFS_Data(cdata%blk_no)%Stateout%gt0,po3=GFS_Interstitial(cdata%thrd_no)%oz_pres, &
                  prsl=GFS_Data(cdata%blk_no)%Statein%prsl,prdout=GFS_Data(cdata%blk_no)%Tbd%ozpl, &
                  oz_coeff=GFS_Interstitial(cdata%thrd_no)%oz_coeff,delp=GFS_Interstitial(cdata%thrd_no)%del, &
                  ldiag3d=GFS_Control%ldiag3d,ozp1=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,6), &
                  ozp2=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,7),ozp3=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,8), &
                  ozp4=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,9),con_g=con_g,me=GFS_Control%me, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in ozphys_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_DCNV_generic_pre_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,ldiag3d=GFS_Control%ldiag3d, &
                  cnvgwd=GFS_Control%cnvgwd,lgocart=GFS_Control%lgocart,gu0=GFS_Data(cdata%blk_no)%Stateout%gu0, &
                  gv0=GFS_Data(cdata%blk_no)%Stateout%gv0,gt0=GFS_Data(cdata%blk_no)%Stateout%gt0, &
                  gq0_water_vapor=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv), &
                  save_u=GFS_Interstitial(cdata%thrd_no)%save_u,save_v=GFS_Interstitial(cdata%thrd_no)%save_v, &
                  save_t=GFS_Interstitial(cdata%thrd_no)%save_t,save_qv=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,1), &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_DCNV_generic_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call get_phi_fv3_run(ix=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,gt0=GFS_Data(cdata%blk_no)%Stateout%gt0, &
                  gq01=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv),del_gz=GFS_Interstitial(cdata%thrd_no)%del_gz, &
                  phii=GFS_Data(cdata%blk_no)%Statein%phii,phil=GFS_Data(cdata%blk_no)%Statein%phil, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in get_phi_fv3_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_suite_interstitial_3_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,nn=GFS_Interstitial(cdata%thrd_no)%nn, &
                  cscnv=GFS_Control%cscnv,satmedmf=GFS_Control%satmedmf,trans_trac=GFS_Control%trans_trac, &
                  do_shoc=GFS_Control%do_shoc,ltaerosol=GFS_Control%ltaerosol,ntrac=GFS_Control%ntrac, &
                  ntcw=GFS_Control%ntcw,ntiw=GFS_Control%ntiw,ntclamt=GFS_Control%ntclamt, &
                  ntrw=GFS_Control%ntrw,ntsw=GFS_Control%ntsw,ntrnc=GFS_Control%ntrnc,ntsnc=GFS_Control%ntsnc, &
                  ntgl=GFS_Control%ntgl,ntgnc=GFS_Control%ntgnc,xlat=GFS_Data(cdata%blk_no)%Grid%xlat, &
                  gq0=GFS_Data(cdata%blk_no)%Stateout%gq0,imp_physics=GFS_Control%imp_physics, &
                  imp_physics_mg=GFS_Control%imp_physics_mg,imp_physics_zhao_carr=GFS_Control%imp_physics_zhao_carr, &
                  imp_physics_zhao_carr_pdf=GFS_Control%imp_physics_zhao_carr_pdf,imp_physics_gfdl=GFS_Control%imp_physics_gfdl, &
                  imp_physics_thompson=GFS_Control%imp_physics_thompson,imp_physics_wsm6=GFS_Control%imp_physics_wsm6, &
                  prsi=GFS_Data(cdata%blk_no)%Statein%prsi,prsl=GFS_Data(cdata%blk_no)%Statein%prsl, &
                  prslk=GFS_Data(cdata%blk_no)%Statein%prslk,rhcbot=GFS_Interstitial(cdata%thrd_no)%rhcbot, &
                  rhcpbl=GFS_Interstitial(cdata%thrd_no)%rhcpbl,rhctop=GFS_Interstitial(cdata%thrd_no)%rhctop, &
                  rhcmax=GFS_Control%rhcmax,islmsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  work1=GFS_Interstitial(cdata%thrd_no)%work1,work2=GFS_Interstitial(cdata%thrd_no)%work2, &
                  kpbl=GFS_Interstitial(cdata%thrd_no)%kpbl,clw=GFS_Interstitial(cdata%thrd_no)%clw, &
                  rhc=GFS_Interstitial(cdata%thrd_no)%rhc,save_qc=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,GFS_Control%ntcw), &
                  save_qi=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,GFS_Control%ntiw),errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_suite_interstitial_3_run"
        ierr=cdata%errflg
        return
      end if

    
      call samfdeepcnv_run(im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%levs, &
                  cliq=con_cliq,cp=con_cp,cvap=con_cvap,eps=con_eps,epsm1=con_epsm1,fv=con_fvirt, &
                  grav=con_g,hvap=con_hvap,rd=con_rd,rv=con_rv,t0c=con_t0c,delt=GFS_Control%dtp, &
                  ntk=GFS_Interstitial(cdata%thrd_no)%ntk,ntr=GFS_Interstitial(cdata%thrd_no)%nsamftrac, &
                  delp=GFS_Interstitial(cdata%thrd_no)%del,prslp=GFS_Data(cdata%blk_no)%Statein%prsl, &
                  psp=GFS_Data(cdata%blk_no)%Statein%pgr,phil=GFS_Data(cdata%blk_no)%Statein%phil, &
                  qtr=GFS_Interstitial(cdata%thrd_no)%clw,q1=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv), &
                  t1=GFS_Data(cdata%blk_no)%Stateout%gt0,u1=GFS_Data(cdata%blk_no)%Stateout%gu0, &
                  v1=GFS_Data(cdata%blk_no)%Stateout%gv0,cldwrk=GFS_Interstitial(cdata%thrd_no)%cld1d, &
                  rn=GFS_Interstitial(cdata%thrd_no)%raincd,kbot=GFS_Interstitial(cdata%thrd_no)%kbot, &
                  ktop=GFS_Interstitial(cdata%thrd_no)%ktop,kcnv=GFS_Interstitial(cdata%thrd_no)%kcnv, &
                  islimsk=GFS_Interstitial(cdata%thrd_no)%islmsk,garea=GFS_Data(cdata%blk_no)%Grid%area, &
                  dot=GFS_Data(cdata%blk_no)%Statein%vvl,ncloud=GFS_Control%ncld,ud_mf=GFS_Interstitial(cdata%thrd_no)%ud_mf, &
                  dd_mf=GFS_Interstitial(cdata%thrd_no)%dd_mf,dt_mf=GFS_Interstitial(cdata%thrd_no)%dt_mf, &
                  cnvw=GFS_Interstitial(cdata%thrd_no)%cnvw,cnvc=GFS_Interstitial(cdata%thrd_no)%cnvc, &
                  qlcn=GFS_Interstitial(cdata%thrd_no)%qlcn,qicn=GFS_Interstitial(cdata%thrd_no)%qicn, &
                  w_upi=GFS_Interstitial(cdata%thrd_no)%w_upi,cf_upi=GFS_Interstitial(cdata%thrd_no)%cf_upi, &
                  cnv_mfd=GFS_Interstitial(cdata%thrd_no)%cnv_mfd,cnv_dqldt=GFS_Interstitial(cdata%thrd_no)%cnv_dqldt, &
                  clcn=GFS_Interstitial(cdata%thrd_no)%clcn,cnv_fice=GFS_Interstitial(cdata%thrd_no)%cnv_fice, &
                  cnv_ndrop=GFS_Interstitial(cdata%thrd_no)%cnv_ndrop,cnv_nice=GFS_Interstitial(cdata%thrd_no)%cnv_nice, &
                  mp_phys=GFS_Control%imp_physics,mp_phys_mg=GFS_Control%imp_physics_mg,clam=GFS_Control%clam_deep, &
                  c0s=GFS_Control%c0s_deep,c1=GFS_Control%c1_deep,betal=GFS_Control%betal_deep, &
                  betas=GFS_Control%betas_deep,evfact=GFS_Control%evfact_deep,evfactl=GFS_Control%evfactl_deep, &
                  pgcon=GFS_Control%pgcon_deep,asolfac=GFS_Control%asolfac_deep,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in samfdeepcnv_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_DCNV_generic_post_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,lssav=GFS_Control%lssav, &
                  ldiag3d=GFS_Control%ldiag3d,lgocart=GFS_Control%lgocart,ras=GFS_Control%ras, &
                  cscnv=GFS_Control%cscnv,frain=GFS_Interstitial(cdata%thrd_no)%frain,rain1=GFS_Interstitial(cdata%thrd_no)%raincd, &
                  dtf=GFS_Control%dtf,cld1d=GFS_Interstitial(cdata%thrd_no)%cld1d,save_u=GFS_Interstitial(cdata%thrd_no)%save_u, &
                  save_v=GFS_Interstitial(cdata%thrd_no)%save_v,save_t=GFS_Interstitial(cdata%thrd_no)%save_t, &
                  save_qv=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,1),gu0=GFS_Data(cdata%blk_no)%Stateout%gu0, &
                  gv0=GFS_Data(cdata%blk_no)%Stateout%gv0,gt0=GFS_Data(cdata%blk_no)%Stateout%gt0, &
                  gq0_water_vapor=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv), &
                  ud_mf=GFS_Interstitial(cdata%thrd_no)%ud_mf,dd_mf=GFS_Interstitial(cdata%thrd_no)%dd_mf, &
                  dt_mf=GFS_Interstitial(cdata%thrd_no)%dt_mf,con_g=con_g,clw_ice=GFS_Interstitial(cdata%thrd_no)%clw(:,:,1), &
                  clw_liquid=GFS_Interstitial(cdata%thrd_no)%clw(:,:,2),npdf3d=GFS_Control%npdf3d, &
                  num_p3d=GFS_Control%num_p3d,ncnvcld3d=GFS_Control%ncnvcld3d,rainc=GFS_Data(cdata%blk_no)%Intdiag%rainc, &
                  cldwrk=GFS_Data(cdata%blk_no)%Intdiag%cldwrk,cnvprcp=GFS_Data(cdata%blk_no)%Intdiag%cnvprcp, &
                  cnvprcpb=GFS_Data(cdata%blk_no)%Intdiag%cnvprcpb,dt3dt=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,4), &
                  dq3dt=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,2),du3dt=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,3), &
                  dv3dt=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,3),upd_mf=GFS_Data(cdata%blk_no)%Intdiag%upd_mf, &
                  dwn_mf=GFS_Data(cdata%blk_no)%Intdiag%dwn_mf,det_mf=GFS_Data(cdata%blk_no)%Intdiag%det_mf, &
                  dqdti=GFS_Data(cdata%blk_no)%Coupling%dqdti,cnvqci=GFS_Data(cdata%blk_no)%Coupling%cnvqci, &
                  upd_mfi=GFS_Data(cdata%blk_no)%Coupling%upd_mfi,dwn_mfi=GFS_Data(cdata%blk_no)%Coupling%dwn_mfi, &
                  det_mfi=GFS_Data(cdata%blk_no)%Coupling%det_mfi,cnvw=GFS_Interstitial(cdata%thrd_no)%cnvw, &
                  cnvc=GFS_Interstitial(cdata%thrd_no)%cnvc,cnvw_phy_f3d=GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ncnvw), &
                  cnvc_phy_f3d=GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ncnvc), &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_DCNV_generic_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call gwdc_pre_run(im=GFS_Control%blksz(cdata%blk_no),cgwf=GFS_Control%cgwf,dx=GFS_Data(cdata%blk_no)%Grid%dx, &
                  work1=GFS_Interstitial(cdata%thrd_no)%work1,work2=GFS_Interstitial(cdata%thrd_no)%work2, &
                  dlength=GFS_Interstitial(cdata%thrd_no)%dlength,cldf=GFS_Interstitial(cdata%thrd_no)%cldf, &
                  levs=GFS_Control%levs,kbot=GFS_Interstitial(cdata%thrd_no)%kbot,ktop=GFS_Interstitial(cdata%thrd_no)%ktop, &
                  dtp=GFS_Control%dtp,gt0=GFS_Data(cdata%blk_no)%Stateout%gt0,gt0_init=GFS_Interstitial(cdata%thrd_no)%save_t, &
                  del=GFS_Interstitial(cdata%thrd_no)%del,cumabs=GFS_Interstitial(cdata%thrd_no)%cumabs, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in gwdc_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call gwdc_run(im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%levs, &
                  lat=GFS_Interstitial(cdata%thrd_no)%latidxprnt,u1=GFS_Data(cdata%blk_no)%Statein%ugrs, &
                  v1=GFS_Data(cdata%blk_no)%Statein%vgrs,t1=GFS_Data(cdata%blk_no)%Statein%tgrs, &
                  q1=GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntqv),deltim=GFS_Control%dtp, &
                  pmid1=GFS_Data(cdata%blk_no)%Statein%prsl,pint1=GFS_Data(cdata%blk_no)%Statein%prsi, &
                  dpmid1=GFS_Interstitial(cdata%thrd_no)%del,qmax=GFS_Interstitial(cdata%thrd_no)%cumabs, &
                  ktop=GFS_Interstitial(cdata%thrd_no)%ktop,kbot=GFS_Interstitial(cdata%thrd_no)%kbot, &
                  kcnv=GFS_Interstitial(cdata%thrd_no)%kcnv,cldf=GFS_Interstitial(cdata%thrd_no)%cldf, &
                  grav=con_g,cp=con_cp,rd=con_rd,fv=con_fvirt,pi=con_pi,dlength=GFS_Interstitial(cdata%thrd_no)%dlength, &
                  lprnt=GFS_Control%lprnt,ipr=GFS_Interstitial(cdata%thrd_no)%ipr,fhour=GFS_Control%fhour, &
                  utgwc=GFS_Interstitial(cdata%thrd_no)%gwdcu,vtgwc=GFS_Interstitial(cdata%thrd_no)%gwdcv, &
                  tauctx=GFS_Interstitial(cdata%thrd_no)%dusfcg,taucty=GFS_Interstitial(cdata%thrd_no)%dvsfcg, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in gwdc_run"
        ierr=cdata%errflg
        return
      end if

    
      call gwdc_post_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,lssav=GFS_Control%lssav, &
                  ldiag3d=GFS_Control%ldiag3d,dtf=GFS_Control%dtf,dtp=GFS_Control%dtp,con_cp=con_cp, &
                  tauctx=GFS_Interstitial(cdata%thrd_no)%dusfcg,taucty=GFS_Interstitial(cdata%thrd_no)%dvsfcg, &
                  gwdcu=GFS_Interstitial(cdata%thrd_no)%gwdcu,gwdcv=GFS_Interstitial(cdata%thrd_no)%gwdcv, &
                  dugwd=GFS_Data(cdata%blk_no)%Intdiag%dugwd,dvgwd=GFS_Data(cdata%blk_no)%Intdiag%dvgwd, &
                  du3dt=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,4),dv3dt=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,4), &
                  gu0=GFS_Data(cdata%blk_no)%Stateout%gu0,gv0=GFS_Data(cdata%blk_no)%Stateout%gv0, &
                  gt0=GFS_Data(cdata%blk_no)%Stateout%gt0,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in gwdc_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_SCNV_generic_pre_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,ldiag3d=GFS_Control%ldiag3d, &
                  lgocart=GFS_Control%lgocart,gt0=GFS_Data(cdata%blk_no)%Stateout%gt0,gq0_water_vapor=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv), &
                  save_t=GFS_Interstitial(cdata%thrd_no)%save_t,save_qv=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,1), &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_SCNV_generic_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call samfshalcnv_run(im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%levs, &
                  cliq=con_cliq,cp=con_cp,cvap=con_cvap,eps=con_eps,epsm1=con_epsm1,fv=con_fvirt, &
                  grav=con_g,hvap=con_hvap,rd=con_rd,rv=con_rv,t0c=con_t0c,delt=GFS_Control%dtp, &
                  ntk=GFS_Interstitial(cdata%thrd_no)%ntk,ntr=GFS_Interstitial(cdata%thrd_no)%nsamftrac, &
                  delp=GFS_Interstitial(cdata%thrd_no)%del,prslp=GFS_Data(cdata%blk_no)%Statein%prsl, &
                  psp=GFS_Data(cdata%blk_no)%Statein%pgr,phil=GFS_Data(cdata%blk_no)%Statein%phil, &
                  qtr=GFS_Interstitial(cdata%thrd_no)%clw,q1=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv), &
                  t1=GFS_Data(cdata%blk_no)%Stateout%gt0,u1=GFS_Data(cdata%blk_no)%Stateout%gu0, &
                  v1=GFS_Data(cdata%blk_no)%Stateout%gv0,rn=GFS_Interstitial(cdata%thrd_no)%raincs, &
                  kbot=GFS_Interstitial(cdata%thrd_no)%kbot,ktop=GFS_Interstitial(cdata%thrd_no)%ktop, &
                  kcnv=GFS_Interstitial(cdata%thrd_no)%kcnv,islimsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  garea=GFS_Data(cdata%blk_no)%Grid%area,dot=GFS_Data(cdata%blk_no)%Statein%vvl, &
                  ncloud=GFS_Control%ncld,hpbl=GFS_Data(cdata%blk_no)%Intdiag%hpbl,ud_mf=GFS_Interstitial(cdata%thrd_no)%ud_mf, &
                  dt_mf=GFS_Interstitial(cdata%thrd_no)%dt_mf,cnvw=GFS_Interstitial(cdata%thrd_no)%cnvw, &
                  cnvc=GFS_Interstitial(cdata%thrd_no)%cnvc,clam=GFS_Control%clam_shal,c0s=GFS_Control%c0s_shal, &
                  c1=GFS_Control%c1_shal,pgcon=GFS_Control%pgcon_shal,asolfac=GFS_Control%asolfac_shal, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in samfshalcnv_run"
        ierr=cdata%errflg
        return
      end if

    
      call samfshalcnv_post_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,lssav=GFS_Control%lssav, &
                  shcnvcw=GFS_Control%shcnvcw,frain=GFS_Interstitial(cdata%thrd_no)%frain, &
                  rain1=GFS_Interstitial(cdata%thrd_no)%raincs,npdf3d=GFS_Control%npdf3d, &
                  num_p3d=GFS_Control%num_p3d,ncnvcld3d=GFS_Control%ncnvcld3d,cnvc=GFS_Interstitial(cdata%thrd_no)%cnvc, &
                  cnvw=GFS_Interstitial(cdata%thrd_no)%cnvw,rainc=GFS_Data(cdata%blk_no)%Intdiag%rainc, &
                  cnvprcp=GFS_Data(cdata%blk_no)%Intdiag%cnvprcp,cnvprcpb=GFS_Data(cdata%blk_no)%Intdiag%cnvprcpb, &
                  cnvw_phy_f3d=GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ncnvw), &
                  cnvc_phy_f3d=GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ncnvc), &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in samfshalcnv_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_SCNV_generic_post_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,nn=GFS_Interstitial(cdata%thrd_no)%nn, &
                  lssav=GFS_Control%lssav,ldiag3d=GFS_Control%ldiag3d,lgocart=GFS_Control%lgocart, &
                  frain=GFS_Interstitial(cdata%thrd_no)%frain,gt0=GFS_Data(cdata%blk_no)%Stateout%gt0, &
                  gq0_water_vapor=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv), &
                  save_t=GFS_Interstitial(cdata%thrd_no)%save_t,save_qv=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,1), &
                  dqdti=GFS_Data(cdata%blk_no)%Coupling%dqdti,dt3dt=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,5), &
                  dq3dt=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,3),clw=GFS_Interstitial(cdata%thrd_no)%clw, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_SCNV_generic_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_suite_interstitial_4_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,ltaerosol=GFS_Control%ltaerosol, &
                  lgocart=GFS_Control%lgocart,tracers_total=GFS_Interstitial(cdata%thrd_no)%tracers_total, &
                  ntrac=GFS_Control%ntrac,ntcw=GFS_Control%ntcw,ntiw=GFS_Control%ntiw,ntclamt=GFS_Control%ntclamt, &
                  ntrw=GFS_Control%ntrw,ntsw=GFS_Control%ntsw,ntrnc=GFS_Control%ntrnc,ntsnc=GFS_Control%ntsnc, &
                  ntgl=GFS_Control%ntgl,ntgnc=GFS_Control%ntgnc,ntlnc=GFS_Control%ntlnc,ntinc=GFS_Control%ntinc, &
                  nn=GFS_Interstitial(cdata%thrd_no)%nn,imp_physics=GFS_Control%imp_physics, &
                  imp_physics_gfdl=GFS_Control%imp_physics_gfdl,imp_physics_thompson=GFS_Control%imp_physics_thompson, &
                  imp_physics_zhao_carr=GFS_Control%imp_physics_zhao_carr,imp_physics_zhao_carr_pdf=GFS_Control%imp_physics_zhao_carr_pdf, &
                  dtf=GFS_Control%dtf,save_qc=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,GFS_Control%ntcw), &
                  save_qi=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,GFS_Control%ntiw),con_pi=con_pi, &
                  gq0=GFS_Data(cdata%blk_no)%Stateout%gq0,clw=GFS_Interstitial(cdata%thrd_no)%clw, &
                  dqdti=GFS_Data(cdata%blk_no)%Coupling%dqdti,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_suite_interstitial_4_run"
        ierr=cdata%errflg
        return
      end if

    
      call cnvc90_run(clstp=GFS_Control%clstp,im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no), &
                  rn=GFS_Data(cdata%blk_no)%Intdiag%rainc,kbot=GFS_Interstitial(cdata%thrd_no)%kbot, &
                  ktop=GFS_Interstitial(cdata%thrd_no)%ktop,km=GFS_Control%levs,prsi=GFS_Data(cdata%blk_no)%Statein%prsi, &
                  acv=GFS_Data(cdata%blk_no)%Tbd%acv,acvb=GFS_Data(cdata%blk_no)%Tbd%acvb, &
                  acvt=GFS_Data(cdata%blk_no)%Tbd%acvt,cv=GFS_Data(cdata%blk_no)%Cldprop%cv, &
                  cvb=GFS_Data(cdata%blk_no)%Cldprop%cvb,cvt=GFS_Data(cdata%blk_no)%Cldprop%cvt, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in cnvc90_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_MP_generic_pre_run(im=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs,ldiag3d=GFS_Control%ldiag3d, &
                  do_aw=GFS_Control%do_aw,ntcw=GFS_Control%ntcw,nncl=GFS_Interstitial(cdata%thrd_no)%nncl, &
                  ntrac=GFS_Control%ntrac,gt0=GFS_Data(cdata%blk_no)%Stateout%gt0,gq0=GFS_Data(cdata%blk_no)%Stateout%gq0, &
                  save_t=GFS_Interstitial(cdata%thrd_no)%save_t,save_q=GFS_Interstitial(cdata%thrd_no)%save_q, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_MP_generic_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call zhaocarr_gscond_run(im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%levs, &
                  dt=GFS_Control%dtp,dtf=GFS_Control%dtf,prsl=GFS_Data(cdata%blk_no)%Statein%prsl, &
                  ps=GFS_Data(cdata%blk_no)%Statein%pgr,q=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv), &
                  clw1=GFS_Interstitial(cdata%thrd_no)%clw(:,:,1),clw2=GFS_Interstitial(cdata%thrd_no)%clw(:,:,2), &
                  cwm=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntcw),t=GFS_Data(cdata%blk_no)%Stateout%gt0, &
                  tp=GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,1),qp=GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,2), &
                  psp=GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,1),tp1=GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,3), &
                  qp1=GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,4),psp1=GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,2), &
                  u=GFS_Interstitial(cdata%thrd_no)%rhc,lprnt=GFS_Control%lprnt,ipr=GFS_Interstitial(cdata%thrd_no)%ipr, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in zhaocarr_gscond_run"
        ierr=cdata%errflg
        return
      end if

    
      call zhaocarr_precpd_run(im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%levs, &
                  dt=GFS_Control%dtp,del=GFS_Interstitial(cdata%thrd_no)%del,prsl=GFS_Data(cdata%blk_no)%Statein%prsl, &
                  q=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv),cwm=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntcw), &
                  t=GFS_Data(cdata%blk_no)%Stateout%gt0,rn=GFS_Interstitial(cdata%thrd_no)%prcpmp, &
                  sr=GFS_Data(cdata%blk_no)%Sfcprop%sr,rainp=GFS_Interstitial(cdata%thrd_no)%rainp, &
                  u00k=GFS_Interstitial(cdata%thrd_no)%rhc,psautco=GFS_Control%psautco,prautco=GFS_Control%prautco, &
                  evpco=GFS_Control%evpco,wminco=GFS_Control%wminco,wk1=GFS_Interstitial(cdata%thrd_no)%work1, &
                  lprnt=GFS_Control%lprnt,jpr=GFS_Interstitial(cdata%thrd_no)%ipr,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in zhaocarr_precpd_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_MP_generic_post_run(im=GFS_Control%blksz(cdata%blk_no),ix=GFS_Control%blksz(cdata%blk_no),levs=GFS_Control%levs, &
                  kdt=GFS_Control%kdt,nrcm=GFS_Control%nrcm,ncld=GFS_Control%ncld,nncl=GFS_Interstitial(cdata%thrd_no)%nncl, &
                  ntcw=GFS_Control%ntcw,ntrac=GFS_Control%ntrac,imp_physics=GFS_Control%imp_physics, &
                  imp_physics_gfdl=GFS_Control%imp_physics_gfdl,imp_physics_thompson=GFS_Control%imp_physics_thompson, &
                  cal_pre=GFS_Control%cal_pre,lssav=GFS_Control%lssav,ldiag3d=GFS_Control%ldiag3d, &
                  cplflx=GFS_Control%cplflx,cplchm=GFS_Control%cplchm,con_g=con_g,dtf=GFS_Control%dtf, &
                  frain=GFS_Interstitial(cdata%thrd_no)%frain,rainc=GFS_Data(cdata%blk_no)%Intdiag%rainc, &
                  rain1=GFS_Interstitial(cdata%thrd_no)%prcpmp,rann=GFS_Data(cdata%blk_no)%Tbd%rann, &
                  xlat=GFS_Data(cdata%blk_no)%Grid%xlat,xlon=GFS_Data(cdata%blk_no)%Grid%xlon, &
                  gt0=GFS_Data(cdata%blk_no)%Stateout%gt0,gq0=GFS_Data(cdata%blk_no)%Stateout%gq0, &
                  prsl=GFS_Data(cdata%blk_no)%Statein%prsl,prsi=GFS_Data(cdata%blk_no)%Statein%prsi, &
                  phii=GFS_Data(cdata%blk_no)%Statein%phii,tsfc=GFS_Data(cdata%blk_no)%Sfcprop%tsfc, &
                  ice=GFS_Data(cdata%blk_no)%Intdiag%ice,snow=GFS_Data(cdata%blk_no)%Intdiag%snow, &
                  graupel=GFS_Data(cdata%blk_no)%Intdiag%graupel,save_t=GFS_Interstitial(cdata%thrd_no)%save_t, &
                  save_qv=GFS_Interstitial(cdata%thrd_no)%save_q(:,:,1),ice0=GFS_Interstitial(cdata%thrd_no)%icemp, &
                  snow0=GFS_Interstitial(cdata%thrd_no)%snowmp,graupel0=GFS_Interstitial(cdata%thrd_no)%graupelmp, &
                  del=GFS_Interstitial(cdata%thrd_no)%del,rain=GFS_Data(cdata%blk_no)%Intdiag%rain, &
                  domr_diag=GFS_Data(cdata%blk_no)%Intdiag%tdomr,domzr_diag=GFS_Data(cdata%blk_no)%Intdiag%tdomzr, &
                  domip_diag=GFS_Data(cdata%blk_no)%Intdiag%tdomip,doms_diag=GFS_Data(cdata%blk_no)%Intdiag%tdoms, &
                  tprcp=GFS_Data(cdata%blk_no)%Sfcprop%tprcp,srflag=GFS_Data(cdata%blk_no)%Sfcprop%srflag, &
                  totprcp=GFS_Data(cdata%blk_no)%Intdiag%totprcp,totice=GFS_Data(cdata%blk_no)%Intdiag%totice, &
                  totsnw=GFS_Data(cdata%blk_no)%Intdiag%totsnw,totgrp=GFS_Data(cdata%blk_no)%Intdiag%totgrp, &
                  totprcpb=GFS_Data(cdata%blk_no)%Intdiag%totprcpb,toticeb=GFS_Data(cdata%blk_no)%Intdiag%toticeb, &
                  totsnwb=GFS_Data(cdata%blk_no)%Intdiag%totsnwb,totgrpb=GFS_Data(cdata%blk_no)%Intdiag%totgrpb, &
                  dt3dt=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,6),dq3dt=GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,4), &
                  rain_cpl=GFS_Data(cdata%blk_no)%Coupling%rain_cpl,rainc_cpl=GFS_Data(cdata%blk_no)%Coupling%rainc_cpl, &
                  snow_cpl=GFS_Data(cdata%blk_no)%Coupling%snow_cpl,pwat=GFS_Data(cdata%blk_no)%Intdiag%pwat, &
                  do_sppt=GFS_Control%do_sppt,dtdtr=GFS_Data(cdata%blk_no)%Tbd%dtdtr,dtdtc=GFS_Interstitial(cdata%thrd_no)%dtdtc, &
                  drain_cpl=GFS_Data(cdata%blk_no)%Tbd%drain_cpl,dsnow_cpl=GFS_Data(cdata%blk_no)%Tbd%dsnow_cpl, &
                  lsm=GFS_Control%lsm,lsm_ruc=GFS_Control%lsm_ruc,raincprv=GFS_Data(cdata%blk_no)%Tbd%raincprv, &
                  rainncprv=GFS_Data(cdata%blk_no)%Tbd%rainncprv,iceprv=GFS_Data(cdata%blk_no)%Tbd%iceprv, &
                  snowprv=GFS_Data(cdata%blk_no)%Tbd%snowprv,graupelprv=GFS_Data(cdata%blk_no)%Tbd%graupelprv, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_MP_generic_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call sfc_sice_post_run(im=GFS_Control%blksz(cdata%blk_no),islmsk=GFS_Interstitial(cdata%thrd_no)%islmsk, &
                  tsfc=GFS_Data(cdata%blk_no)%Sfcprop%tsfc,fice=GFS_Data(cdata%blk_no)%Sfcprop%fice, &
                  hice=GFS_Data(cdata%blk_no)%Sfcprop%hice,tisfc=GFS_Data(cdata%blk_no)%Sfcprop%tisfc, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in sfc_sice_post_run"
        ierr=cdata%errflg
        return
      end if

    



   end function physics_run_cap

   function physics_finalize_cap(cdata) result(ierr)

      use ccpp_types, only: ccpp_t

      implicit none

      integer                     :: ierr

      type(ccpp_t), intent(inout) :: cdata

      ierr = 0


      if (.not.initialized) return



      call lsm_noah_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in lsm_noah_finalize"
        ierr=cdata%errflg
        return
      end if

    


      initialized = .false.


   end function physics_finalize_cap

end module ccpp_group_physics_cap

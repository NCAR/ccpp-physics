
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
!! @brief Auto-generated cap module for the CCPP radiation group
!!
!
module ccpp_group_radiation_cap

   use GFS_suite_interstitial_rad_reset, only: GFS_suite_interstitial_rad_reset_run
   use GFS_rrtmg_pre, only: GFS_rrtmg_pre_run
   use rrtmg_sw_pre, only: rrtmg_sw_pre_run
   use rrtmg_sw, only: rrtmg_sw_run
   use rrtmg_sw_post, only: rrtmg_sw_post_run
   use rrtmg_lw_pre, only: rrtmg_lw_pre_run
   use rrtmg_lw, only: rrtmg_lw_run
   use rrtmg_lw_post, only: rrtmg_lw_post_run
   use GFS_rrtmg_post, only: GFS_rrtmg_post_run


   implicit none

   private
   public :: radiation_init_cap,radiation_run_cap,radiation_finalize_cap

   logical, save :: initialized = .false.

   contains

   function radiation_init_cap() result(ierr)

      

      implicit none

      integer                     :: ierr

      

      ierr = 0


      if (initialized) return





      initialized = .true.


   end function radiation_init_cap

   function radiation_run_cap(cdata,GFS_Interstitial,GFS_Data,GFS_Control,LTP) result(ierr)

      use ccpp_types, only: ccpp_t
      use GFS_typedefs, only: GFS_interstitial_type
      use GFS_typedefs, only: GFS_data_type
      use GFS_typedefs, only: GFS_control_type

      implicit none

      integer                     :: ierr

      type(ccpp_t), intent(inout) :: cdata
      type(GFS_interstitial_type), intent(inout) :: GFS_Interstitial(:)
      type(GFS_data_type), intent(inout) :: GFS_Data(:)
      type(GFS_control_type), intent(in) :: GFS_Control
      integer, intent(in) :: LTP

      ierr = 0


      if (.not.initialized) then
        write(cdata%errmsg,'(*(a))') 'radiation_run called before radiation_init'
        cdata%errflg = 1
        return
      end if



      call GFS_suite_interstitial_rad_reset_run(Interstitial=GFS_Interstitial(cdata%thrd_no),errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_suite_interstitial_rad_reset_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_rrtmg_pre_run(Model=GFS_Control,Grid=GFS_Data(cdata%blk_no)%Grid,Sfcprop=GFS_Data(cdata%blk_no)%Sfcprop, &
                  Statein=GFS_Data(cdata%blk_no)%Statein,Tbd=GFS_Data(cdata%blk_no)%Tbd,Cldprop=GFS_Data(cdata%blk_no)%Cldprop, &
                  Coupling=GFS_Data(cdata%blk_no)%Coupling,Radtend=GFS_Data(cdata%blk_no)%Radtend, &
                  lm=GFS_Interstitial(cdata%thrd_no)%lm,im=GFS_Control%blksz(cdata%blk_no), &
                  lmk=GFS_Interstitial(cdata%thrd_no)%lmk,lmp=GFS_Interstitial(cdata%thrd_no)%lmp, &
                  kd=GFS_Interstitial(cdata%thrd_no)%kd,kt=GFS_Interstitial(cdata%thrd_no)%kt, &
                  kb=GFS_Interstitial(cdata%thrd_no)%kb,raddt=GFS_Interstitial(cdata%thrd_no)%raddt, &
                  delp=GFS_Interstitial(cdata%thrd_no)%delr,dz=GFS_Interstitial(cdata%thrd_no)%dzlyr, &
                  plvl=GFS_Interstitial(cdata%thrd_no)%plvl,plyr=GFS_Interstitial(cdata%thrd_no)%plyr, &
                  tlvl=GFS_Interstitial(cdata%thrd_no)%tlvl,tlyr=GFS_Interstitial(cdata%thrd_no)%tlyr, &
                  tsfg=GFS_Interstitial(cdata%thrd_no)%tsfg,tsfa=GFS_Interstitial(cdata%thrd_no)%tsfa, &
                  qlyr=GFS_Interstitial(cdata%thrd_no)%qlyr,olyr=GFS_Interstitial(cdata%thrd_no)%olyr, &
                  gasvmr_co2=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,1),gasvmr_n2o=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,2), &
                  gasvmr_ch4=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,3),gasvmr_o2=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,4), &
                  gasvmr_co=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,5),gasvmr_cfc11=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,6), &
                  gasvmr_cfc12=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,7),gasvmr_cfc22=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,8), &
                  gasvmr_ccl4=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,9),gasvmr_cfc113=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,10), &
                  faersw1=GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,1),faersw2=GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,2), &
                  faersw3=GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,3),faerlw1=GFS_Interstitial(cdata%thrd_no)%faerlw(:,:,:,1), &
                  faerlw2=GFS_Interstitial(cdata%thrd_no)%faerlw(:,:,:,2),faerlw3=GFS_Interstitial(cdata%thrd_no)%faerlw(:,:,:,3), &
                  aerodp=GFS_Interstitial(cdata%thrd_no)%aerodp,clouds1=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,1), &
                  clouds2=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,2),clouds3=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,3), &
                  clouds4=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,4),clouds5=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,5), &
                  clouds6=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,6),clouds7=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,7), &
                  clouds8=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,8),clouds9=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,9), &
                  cldsa=GFS_Interstitial(cdata%thrd_no)%cldsa,mtopa=GFS_Interstitial(cdata%thrd_no)%mtopa, &
                  mbota=GFS_Interstitial(cdata%thrd_no)%mbota,de_lgth=GFS_Interstitial(cdata%thrd_no)%de_lgth, &
                  alb1d=GFS_Interstitial(cdata%thrd_no)%alb1d,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_rrtmg_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call rrtmg_sw_pre_run(Model=GFS_Control,Grid=GFS_Data(cdata%blk_no)%Grid,Sfcprop=GFS_Data(cdata%blk_no)%Sfcprop, &
                  Radtend=GFS_Data(cdata%blk_no)%Radtend,im=GFS_Control%blksz(cdata%blk_no), &
                  nday=GFS_Interstitial(cdata%thrd_no)%nday,idxday=GFS_Interstitial(cdata%thrd_no)%idxday, &
                  tsfg=GFS_Interstitial(cdata%thrd_no)%tsfg,tsfa=GFS_Interstitial(cdata%thrd_no)%tsfa, &
                  sfcalb1=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,1),sfcalb2=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,2), &
                  sfcalb3=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,3),sfcalb4=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,4), &
                  alb1d=GFS_Interstitial(cdata%thrd_no)%alb1d,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in rrtmg_sw_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call rrtmg_sw_run(plyr=GFS_Interstitial(cdata%thrd_no)%plyr,plvl=GFS_Interstitial(cdata%thrd_no)%plvl, &
                  tlyr=GFS_Interstitial(cdata%thrd_no)%tlyr,tlvl=GFS_Interstitial(cdata%thrd_no)%tlvl, &
                  qlyr=GFS_Interstitial(cdata%thrd_no)%qlyr,olyr=GFS_Interstitial(cdata%thrd_no)%olyr, &
                  gasvmr_co2=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,1),gasvmr_n2o=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,2), &
                  gasvmr_ch4=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,3),gasvmr_o2=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,4), &
                  gasvmr_co=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,5),gasvmr_cfc11=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,6), &
                  gasvmr_cfc12=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,7),gasvmr_cfc22=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,8), &
                  gasvmr_ccl4=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,9),icseed=GFS_Data(cdata%blk_no)%Tbd%icsdsw, &
                  aeraod=GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,1),aerssa=GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,2), &
                  aerasy=GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,3),sfcalb_nir_dir=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,1), &
                  sfcalb_nir_dif=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,2),sfcalb_uvis_dir=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,3), &
                  sfcalb_uvis_dif=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,4),dzlyr=GFS_Interstitial(cdata%thrd_no)%dzlyr, &
                  delpin=GFS_Interstitial(cdata%thrd_no)%delr,de_lgth=GFS_Interstitial(cdata%thrd_no)%de_lgth, &
                  cosz=GFS_Data(cdata%blk_no)%Radtend%coszen,solcon=GFS_Control%solcon,nday=GFS_Interstitial(cdata%thrd_no)%nday, &
                  idxday=GFS_Interstitial(cdata%thrd_no)%idxday,npts=GFS_Control%blksz(cdata%blk_no), &
                  nlay=GFS_Interstitial(cdata%thrd_no)%lmk,nlp1=GFS_Interstitial(cdata%thrd_no)%lmp, &
                  lprnt=GFS_Control%lprnt,cld_cf=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,1), &
                  lsswr=GFS_Control%lsswr,hswc=GFS_Data(cdata%blk_no)%Tbd%htswc,topflx=GFS_Data(cdata%blk_no)%Intdiag%topfsw, &
                  sfcflx=GFS_Data(cdata%blk_no)%Radtend%sfcfsw,cldtau=GFS_Interstitial(cdata%thrd_no)%cldtausw, &
                  hsw0=GFS_Data(cdata%blk_no)%Tbd%htsw0,fdncmp=GFS_Interstitial(cdata%thrd_no)%scmpsw, &
                  cld_lwp=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,2),cld_ref_liq=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,3), &
                  cld_iwp=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,4),cld_ref_ice=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,5), &
                  cld_rwp=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,6),cld_ref_rain=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,7), &
                  cld_swp=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,8),cld_ref_snow=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,9), &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in rrtmg_sw_run"
        ierr=cdata%errflg
        return
      end if

    
      call rrtmg_sw_post_run(Model=GFS_Control,Grid=GFS_Data(cdata%blk_no)%Grid,Diag=GFS_Data(cdata%blk_no)%Intdiag, &
                  Radtend=GFS_Data(cdata%blk_no)%Radtend,Coupling=GFS_Data(cdata%blk_no)%Coupling, &
                  im=GFS_Control%blksz(cdata%blk_no),ltp=LTP,nday=GFS_Interstitial(cdata%thrd_no)%nday, &
                  lm=GFS_Interstitial(cdata%thrd_no)%lm,kd=GFS_Interstitial(cdata%thrd_no)%kd, &
                  htswc=GFS_Data(cdata%blk_no)%Tbd%htswc,htsw0=GFS_Data(cdata%blk_no)%Tbd%htsw0, &
                  sfcalb1=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,1),sfcalb2=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,2), &
                  sfcalb3=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,3),sfcalb4=GFS_Interstitial(cdata%thrd_no)%sfcalb(:,4), &
                  scmpsw=GFS_Interstitial(cdata%thrd_no)%scmpsw,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in rrtmg_sw_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call rrtmg_lw_pre_run(Model=GFS_Control,Grid=GFS_Data(cdata%blk_no)%Grid,Sfcprop=GFS_Data(cdata%blk_no)%Sfcprop, &
                  Radtend=GFS_Data(cdata%blk_no)%Radtend,im=GFS_Control%blksz(cdata%blk_no), &
                  tsfg=GFS_Interstitial(cdata%thrd_no)%tsfg,tsfa=GFS_Interstitial(cdata%thrd_no)%tsfa, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in rrtmg_lw_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call rrtmg_lw_run(plyr=GFS_Interstitial(cdata%thrd_no)%plyr,plvl=GFS_Interstitial(cdata%thrd_no)%plvl, &
                  tlyr=GFS_Interstitial(cdata%thrd_no)%tlyr,tlvl=GFS_Interstitial(cdata%thrd_no)%tlvl, &
                  qlyr=GFS_Interstitial(cdata%thrd_no)%qlyr,olyr=GFS_Interstitial(cdata%thrd_no)%olyr, &
                  gasvmr_co2=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,1),gasvmr_n2o=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,2), &
                  gasvmr_ch4=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,3),gasvmr_o2=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,4), &
                  gasvmr_co=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,5),gasvmr_cfc11=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,6), &
                  gasvmr_cfc12=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,7),gasvmr_cfc22=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,8), &
                  gasvmr_ccl4=GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,9),icseed=GFS_Data(cdata%blk_no)%Tbd%icsdlw, &
                  aeraod=GFS_Interstitial(cdata%thrd_no)%faerlw(:,:,:,1),aerssa=GFS_Interstitial(cdata%thrd_no)%faerlw(:,:,:,2), &
                  sfemis=GFS_Data(cdata%blk_no)%Radtend%semis,sfgtmp=GFS_Interstitial(cdata%thrd_no)%tsfg, &
                  dzlyr=GFS_Interstitial(cdata%thrd_no)%dzlyr,delpin=GFS_Interstitial(cdata%thrd_no)%delr, &
                  de_lgth=GFS_Interstitial(cdata%thrd_no)%de_lgth,npts=GFS_Control%blksz(cdata%blk_no), &
                  nlay=GFS_Interstitial(cdata%thrd_no)%lmk,nlp1=GFS_Interstitial(cdata%thrd_no)%lmp, &
                  lprnt=GFS_Control%lprnt,cld_cf=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,1), &
                  lslwr=GFS_Control%lslwr,hlwc=GFS_Data(cdata%blk_no)%Tbd%htlwc,topflx=GFS_Data(cdata%blk_no)%Intdiag%topflw, &
                  sfcflx=GFS_Data(cdata%blk_no)%Radtend%sfcflw,cldtau=GFS_Interstitial(cdata%thrd_no)%cldtaulw, &
                  hlw0=GFS_Data(cdata%blk_no)%Tbd%htlw0,cld_lwp=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,2), &
                  cld_ref_liq=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,3),cld_iwp=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,4), &
                  cld_ref_ice=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,5),cld_rwp=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,6), &
                  cld_ref_rain=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,7),cld_swp=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,8), &
                  cld_ref_snow=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,9),errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in rrtmg_lw_run"
        ierr=cdata%errflg
        return
      end if

    
      call rrtmg_lw_post_run(Model=GFS_Control,Grid=GFS_Data(cdata%blk_no)%Grid,Radtend=GFS_Data(cdata%blk_no)%Radtend, &
                  Coupling=GFS_Data(cdata%blk_no)%Coupling,im=GFS_Control%blksz(cdata%blk_no), &
                  ltp=LTP,lm=GFS_Interstitial(cdata%thrd_no)%lm,kd=GFS_Interstitial(cdata%thrd_no)%kd, &
                  tsfa=GFS_Interstitial(cdata%thrd_no)%tsfa,htlwc=GFS_Data(cdata%blk_no)%Tbd%htlwc, &
                  htlw0=GFS_Data(cdata%blk_no)%Tbd%htlw0,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in rrtmg_lw_post_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_rrtmg_post_run(Model=GFS_Control,Grid=GFS_Data(cdata%blk_no)%Grid,Diag=GFS_Data(cdata%blk_no)%Intdiag, &
                  Radtend=GFS_Data(cdata%blk_no)%Radtend,Statein=GFS_Data(cdata%blk_no)%Statein, &
                  Coupling=GFS_Data(cdata%blk_no)%Coupling,scmpsw=GFS_Interstitial(cdata%thrd_no)%scmpsw, &
                  im=GFS_Control%blksz(cdata%blk_no),lm=GFS_Interstitial(cdata%thrd_no)%lm, &
                  ltp=LTP,kt=GFS_Interstitial(cdata%thrd_no)%kt,kb=GFS_Interstitial(cdata%thrd_no)%kb, &
                  kd=GFS_Interstitial(cdata%thrd_no)%kd,raddt=GFS_Interstitial(cdata%thrd_no)%raddt, &
                  aerodp=GFS_Interstitial(cdata%thrd_no)%aerodp,cldsa=GFS_Interstitial(cdata%thrd_no)%cldsa, &
                  mtopa=GFS_Interstitial(cdata%thrd_no)%mtopa,mbota=GFS_Interstitial(cdata%thrd_no)%mbota, &
                  clouds1=GFS_Interstitial(cdata%thrd_no)%clouds(:,:,1),cldtaulw=GFS_Interstitial(cdata%thrd_no)%cldtaulw, &
                  cldtausw=GFS_Interstitial(cdata%thrd_no)%cldtausw,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_rrtmg_post_run"
        ierr=cdata%errflg
        return
      end if

    



   end function radiation_run_cap

   function radiation_finalize_cap() result(ierr)

      

      implicit none

      integer                     :: ierr

      

      ierr = 0


      if (.not.initialized) return





      initialized = .false.


   end function radiation_finalize_cap

end module ccpp_group_radiation_cap

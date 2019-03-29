
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
!! @brief Auto-generated cap module for the CCPP stochastics group
!!
!
module ccpp_group_stochastics_cap

   use GFS_stochastics, only: GFS_stochastics_run


   implicit none

   private
   public :: stochastics_init_cap,stochastics_run_cap,stochastics_finalize_cap

   logical, save :: initialized = .false.

   contains

   function stochastics_init_cap() result(ierr)

      

      implicit none

      integer                     :: ierr

      

      ierr = 0


      if (initialized) return





      initialized = .true.


   end function stochastics_init_cap

   function stochastics_run_cap(cdata,GFS_Data,GFS_Control) result(ierr)

      use ccpp_types, only: ccpp_t
      use GFS_typedefs, only: GFS_data_type
      use GFS_typedefs, only: GFS_control_type

      implicit none

      integer                     :: ierr

      type(ccpp_t), intent(inout) :: cdata
      type(GFS_data_type), intent(inout) :: GFS_Data(:)
      type(GFS_control_type), intent(in) :: GFS_Control

      ierr = 0


      if (.not.initialized) then
        write(cdata%errmsg,'(*(a))') 'stochastics_run called before stochastics_init'
        cdata%errflg = 1
        return
      end if



      call GFS_stochastics_run(im=GFS_Control%blksz(cdata%blk_no),km=GFS_Control%levs,do_sppt=GFS_Control%do_sppt, &
                  use_zmtnblck=GFS_Control%use_zmtnblck,do_shum=GFS_Control%do_shum,do_skeb=GFS_Control%do_skeb, &
                  zmtnblck=GFS_Data(cdata%blk_no)%Intdiag%zmtnblck,sppt_wts=GFS_Data(cdata%blk_no)%Coupling%sppt_wts, &
                  skebu_wts=GFS_Data(cdata%blk_no)%Coupling%skebu_wts,skebv_wts=GFS_Data(cdata%blk_no)%Coupling%skebv_wts, &
                  shum_wts=GFS_Data(cdata%blk_no)%Coupling%shum_wts,sppt_wts_inv=GFS_Data(cdata%blk_no)%Intdiag%sppt_wts, &
                  skebu_wts_inv=GFS_Data(cdata%blk_no)%Intdiag%skebu_wts,skebv_wts_inv=GFS_Data(cdata%blk_no)%Intdiag%skebv_wts, &
                  shum_wts_inv=GFS_Data(cdata%blk_no)%Intdiag%shum_wts,diss_est=GFS_Data(cdata%blk_no)%Statein%diss_est, &
                  ugrs=GFS_Data(cdata%blk_no)%Statein%ugrs,vgrs=GFS_Data(cdata%blk_no)%Statein%vgrs, &
                  tgrs=GFS_Data(cdata%blk_no)%Statein%tgrs,qgrs=GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntqv), &
                  gu0=GFS_Data(cdata%blk_no)%Stateout%gu0,gv0=GFS_Data(cdata%blk_no)%Stateout%gv0, &
                  gt0=GFS_Data(cdata%blk_no)%Stateout%gt0,gq0=GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv), &
                  dtdtr=GFS_Data(cdata%blk_no)%Tbd%dtdtr,rain=GFS_Data(cdata%blk_no)%Intdiag%rain, &
                  rainc=GFS_Data(cdata%blk_no)%Intdiag%rainc,tprcp=GFS_Data(cdata%blk_no)%Sfcprop%tprcp, &
                  totprcp=GFS_Data(cdata%blk_no)%Intdiag%totprcp,cnvprcp=GFS_Data(cdata%blk_no)%Intdiag%cnvprcp, &
                  totprcpb=GFS_Data(cdata%blk_no)%Intdiag%totprcpb,cnvprcpb=GFS_Data(cdata%blk_no)%Intdiag%cnvprcpb, &
                  cplflx=GFS_Control%cplflx,rain_cpl=GFS_Data(cdata%blk_no)%Coupling%rain_cpl, &
                  snow_cpl=GFS_Data(cdata%blk_no)%Coupling%snow_cpl,drain_cpl=GFS_Data(cdata%blk_no)%Tbd%drain_cpl, &
                  dsnow_cpl=GFS_Data(cdata%blk_no)%Tbd%dsnow_cpl,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_stochastics_run"
        ierr=cdata%errflg
        return
      end if

    



   end function stochastics_run_cap

   function stochastics_finalize_cap() result(ierr)

      

      implicit none

      integer                     :: ierr

      

      ierr = 0


      if (.not.initialized) return





      initialized = .false.


   end function stochastics_finalize_cap

end module ccpp_group_stochastics_cap

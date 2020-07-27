module phys_tend

   use machine, only: kind_phys

   implicit none

   private

   public phys_tend_init, phys_tend_run, phys_tend_finalize

contains

   subroutine phys_tend_init()
   end subroutine phys_tend_init

   subroutine phys_tend_finalize()
   end subroutine phys_tend_finalize

!> \section arg_table_phys_tend_run Argument Table
!! \htmlinclude phys_tend_run.html
!!
   subroutine phys_tend_run(ldiag3d, qdiag3d,                &
       du3dt_pbl, du3dt_orogwd, du3dt_deepcnv, du3dt_congwd,  &
       du3dt_rdamp, du3dt_shalcnv, du3dt_phys,                &
       dv3dt_pbl, dv3dt_orogwd, dv3dt_deepcnv, dv3dt_congwd,  &
       dv3dt_rdamp, dv3dt_shalcnv, dv3dt_phys,                &
       dt3dt_lw, dt3dt_sw, dt3dt_pbl, dt3dt_deepcnv,          &
       dt3dt_shalcnv, dt3dt_mp, dt3dt_orogwd, dt3dt_rdamp,    &
       dt3dt_congwd, dt3dt_phys,                              &
       dq3dt_pbl, dq3dt_deepcnv, dq3dt_shalcnv, dq3dt_mp,     &
       dq3dt_o3pbl, dq3dt_o3prodloss, dq3dt_o3mix,            &
       dq3dt_o3tmp, dq3dt_o3column, dq3dt_phys, dq3dt_o3phys, &
       errmsg, errflg)

       ! Interface variables
       logical, intent(in) :: ldiag3d, qdiag3d
       real(kind=kind_phys), intent(in   ) :: du3dt_pbl(:,:)
       real(kind=kind_phys), intent(in   ) :: du3dt_orogwd(:,:)
       real(kind=kind_phys), intent(in   ) :: du3dt_deepcnv(:,:)
       real(kind=kind_phys), intent(in   ) :: du3dt_congwd(:,:)
       real(kind=kind_phys), intent(in   ) :: du3dt_rdamp(:,:)
       real(kind=kind_phys), intent(in   ) :: du3dt_shalcnv(:,:)
       real(kind=kind_phys), intent(  out) :: du3dt_phys(:,:)
       real(kind=kind_phys), intent(in   ) :: dv3dt_pbl(:,:)
       real(kind=kind_phys), intent(in   ) :: dv3dt_orogwd(:,:)
       real(kind=kind_phys), intent(in   ) :: dv3dt_deepcnv(:,:)
       real(kind=kind_phys), intent(in   ) :: dv3dt_congwd(:,:)
       real(kind=kind_phys), intent(in   ) :: dv3dt_rdamp(:,:)
       real(kind=kind_phys), intent(in   ) :: dv3dt_shalcnv(:,:)
       real(kind=kind_phys), intent(  out) :: dv3dt_phys(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_lw(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_sw(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_pbl(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_deepcnv(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_shalcnv(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_mp(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_orogwd(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_rdamp(:,:)
       real(kind=kind_phys), intent(in   ) :: dt3dt_congwd(:,:)
       real(kind=kind_phys), intent(  out) :: dt3dt_phys(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_pbl(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_deepcnv(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_shalcnv(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_mp(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_o3pbl(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_o3prodloss(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_o3mix(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_o3tmp(:,:)
       real(kind=kind_phys), intent(in   ) :: dq3dt_o3column(:,:)
       real(kind=kind_phys), intent(  out) :: dq3dt_phys(:,:)
       real(kind=kind_phys), intent(  out) :: dq3dt_o3phys(:,:)
       character(len=*), intent(out) :: errmsg
       integer, intent(out)          :: errflg

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (.not.ldiag3d .and. .not.qdiag3d) return

       du3dt_phys   = du3dt_pbl + du3dt_orogwd + du3dt_deepcnv + &
                      du3dt_congwd + du3dt_rdamp + du3dt_shalcnv

       dv3dt_phys   = dv3dt_pbl + dv3dt_orogwd + dv3dt_deepcnv + &
                      dv3dt_congwd + dv3dt_rdamp + dv3dt_shalcnv

       dt3dt_phys   = dt3dt_lw + dt3dt_sw + dt3dt_pbl +          &
                      dt3dt_deepcnv + dt3dt_shalcnv + dt3dt_mp + &
                      dt3dt_orogwd + dt3dt_rdamp + dt3dt_congwd

       dq3dt_phys   = dq3dt_pbl + dq3dt_deepcnv +                &
                      dq3dt_shalcnv + dq3dt_mp

       dq3dt_o3phys = dq3dt_o3pbl + dq3dt_o3prodloss             &
                      + dq3dt_o3mix + dq3dt_o3tmp + dq3dt_o3column

   end subroutine phys_tend_run

end module phys_tend

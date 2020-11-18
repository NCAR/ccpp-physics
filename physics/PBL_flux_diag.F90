module PBL_flux_diag

   use machine, only: kind_phys

   implicit none

   private

   public PBL_flux_diag_init, PBL_flux_diag_run, PBL_flux_diag_finalize

contains

   subroutine PBL_flux_diag_init()
   end subroutine PBL_flux_diag_init

   subroutine PBL_flux_diag_finalize()
   end subroutine PBL_flux_diag_finalize

!> \section arg_table_PBL_flux_diag_run Argument Table
!! \htmlinclude PBL_flux_diag_run.html
!!
   subroutine PBL_flux_diag_run(im, km, dtp, con_rd, con_g, con_fvirt, ldiag3d, qdiag3d, &
       du_pbl, dv_pbl, dt_pbl, dq_pbl, prslk, prsi, t1, q1,                              &
       sgs_vert_flx_th, sgs_vert_flx_q, sgs_vert_flx_u, sgs_vert_flx_v,                  &
       errmsg, errflg)

       ! Interface variables
       integer, intent(in) :: im, km
       real(kind=kind_phys), intent(in) :: dtp, con_rd, con_g, con_fvirt
       logical, intent(in) :: ldiag3d, qdiag3d
       real(kind=kind_phys), intent(in), dimension(:,:) :: du_pbl, dv_pbl, dt_pbl, dq_pbl, prslk, t1, q1
       real(kind=kind_phys), intent(in), dimension(:,:) :: prsi
       real(kind=kind_phys), intent(out), dimension(:,:) :: sgs_vert_flx_th, sgs_vert_flx_q, sgs_vert_flx_u, sgs_vert_flx_v
       
       character(len=*), intent(out) :: errmsg
       integer, intent(out)          :: errflg
       
       integer :: i,k
       real(kind=kind_phys) :: dt_inv, tv
       real(kind=kind_phys), dimension(im,km) :: dth_pbl, dthdt, dqdt, dudt, dvdt, dz
       
       real(kind=kind_phys), parameter :: qmin=1.e-8
       
       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (.not.ldiag3d .and. .not.qdiag3d) return
       
       !convert temperature change to potential temperature change
       dth_pbl = dt_pbl/prslk
       
       !SCM-only: the cumulative change terms are reset every physics timestep, so the actual tendency is just the change divided by the physics timestep
       dt_inv = 1.0/dtp
       dthdt = dth_pbl*dt_inv
       dqdt = dq_pbl*dt_inv
       dudt = du_pbl*dt_inv
       dvdt = dv_pbl*dt_inv
       
       ! Layer thickness (m)                                                                                                                                                              
       do k=1, km
         do i = 1, im
           tv = t1(i,k)*(1.+con_fvirt*max(q1(i,k),qmin))
           dz(i,k) = ((con_rd/con_g)) * abs(log(prsi(i,k+1)) - log(prsi(i,k))) * tv
         enddo
       enddo
       
       !assume that the fluxes are zero at the top of the model, and integrate downward to obtain SGS fluxes from PBL
       do i=1, im
         sgs_vert_flx_th(i,km) = 0.0
         sgs_vert_flx_q(i,km) = 0.0
         sgs_vert_flx_u(i,km) = 0.0
         sgs_vert_flx_v(i,km) = 0.0
       end do
       do k = km-1, 1, -1
         do i = 1, im
           !positive sign since dz is positive
           write(*,*) i,k,dz(i,k)
           sgs_vert_flx_th(i,k) = sgs_vert_flx_th(i,k+1) + dthdt(i,k)*dz(i,k)
           sgs_vert_flx_q(i,k) = sgs_vert_flx_q(i,k+1) + dqdt(i,k)*dz(i,k)
           sgs_vert_flx_u(i,k) = sgs_vert_flx_u(i,k+1) + dudt(i,k)*dz(i,k)
           sgs_vert_flx_v(i,k) = sgs_vert_flx_v(i,k+1) + dvdt(i,k)*dz(i,k)
         enddo
       enddo

  end subroutine PBL_flux_diag_run

end module PBL_flux_diag

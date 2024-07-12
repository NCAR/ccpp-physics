!>  \file ugwpv1_gsldrag_post.F90
!! This file contains
module ugwpv1_gsldrag_post

contains

!>\defgroup ugwpv1_gsldrag_post ugwpv1_gsldrag Scheme Post
!! @{
!> \section arg_table_ugwpv1_gsldrag_post_run Argument Table
!! \htmlinclude ugwpv1_gsldrag_post_run.html
!!
     subroutine ugwpv1_gsldrag_post_run ( im, levs, ldiag_ugwp,       &
         dtf, dudt_gw, dvdt_gw, dtdt_gw,                              &
         tau_ogw, tau_ngw, zobl, zlwb, zogw, dudt_obl, dvdt_obl,      &
         dudt_ofd, dvdt_ofd, dudt_ogw, dvdt_ogw,                      &
         dudt_oss, dvdt_oss, tot_zmtb, tot_zlwb, tot_zogw,            &
         tot_tofd, tot_mtb, tot_ogw, tot_ngw,                         &
         du3dt_mtb,du3dt_ogw, du3dt_tms, du3dt_ngw, dv3dt_ngw,        &
         dudt_ngw, dvdt_ngw, dtdt_ngw,                                &
         ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw,                          &
         dws3dt_ogw, dws3dt_obl, dws3dt_oss, dws3dt_ofd,              &
         ldu3dt_ogw, ldu3dt_obl, ldu3dt_oss, ldu3dt_ofd,              &
         du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol, du_osscol,       &
         dv_osscol, du_ofdcol, dv_ofdcol, du3_ogwcol, dv3_ogwcol,     &
         du3_oblcol, dv3_oblcol, du3_osscol, dv3_osscol, du3_ofdcol,  &
         dv3_ofdcol, dtdt, dudt, dvdt, errmsg, errflg)

        use machine,                only: kind_phys

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: ldiag_ugwp      !< flag for CIRES UGWP Diagnostics

        real(kind=kind_phys), intent(in),    dimension(:)   :: zobl, zlwb, zogw
        real(kind=kind_phys), intent(in),    dimension(:)   :: tau_ogw, tau_ngw
        real(kind=kind_phys), intent(in),    dimension(:),optional   :: du_ofdcol, du_oblcol
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_mtb, tot_ogw, tot_tofd, tot_ngw
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_zmtb, tot_zlwb, tot_zogw
	
        real(kind=kind_phys), intent(in),    dimension(:,:) :: dtdt_gw, dudt_gw, dvdt_gw
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_obl, dvdt_obl, dudt_ogw
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dvdt_ogw, dudt_ofd, dvdt_ofd
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_oss, dvdt_oss
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: du3dt_mtb, du3dt_ogw, du3dt_tms
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: du3dt_ngw, dv3dt_ngw
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_ngw, dvdt_ngw, dtdt_ngw
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: dws3dt_ogw, dws3dt_obl
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: dws3dt_oss, dws3dt_ofd
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: ldu3dt_ogw, ldu3dt_obl
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: ldu3dt_oss, ldu3dt_ofd
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: du_ogwcol, dv_ogwcol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: dv_oblcol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: du_osscol, dv_osscol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: dv_ofdcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_ogwcol, dv3_ogwcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_oblcol, dv3_oblcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_osscol, dv3_osscol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_ofdcol, dv3_ofdcol
	
        real(kind=kind_phys), intent(inout), dimension(:,:) :: dtdt, dudt, dvdt

        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0
!
! post creates the "time-averaged" diagnostics"
!

        if (ldiag_ugwp) then
          tot_zmtb =  tot_zmtb + dtf *zobl
          tot_zlwb =  tot_zlwb + dtf *zlwb
          tot_zogw =  tot_zogw + dtf *zogw
    
          tot_tofd  = tot_tofd + dtf *du_ofdcol
          tot_mtb   = tot_mtb +  dtf *du_oblcol
          tot_ogw   = tot_ogw +  dtf *tau_ogw
          tot_ngw   = tot_ngw +  dtf *tau_ngw
    
          du3dt_mtb = du3dt_mtb + dtf *dudt_obl
          du3dt_tms = du3dt_tms + dtf *dudt_ofd
          du3dt_ogw = du3dt_ogw + dtf *dudt_ogw
          du3dt_ngw = du3dt_ngw + dtf *dudt_gw
          dv3dt_ngw = dv3dt_ngw + dtf *dvdt_gw

          dws3dt_ogw = dws3dt_ogw + dtf *sqrt(dudt_ogw**2+dvdt_ogw**2)
          dws3dt_obl = dws3dt_obl + dtf *sqrt(dudt_obl**2+dvdt_obl**2)

          du3_ogwcol = du3_ogwcol + dtf *du_ogwcol
          dv3_ogwcol = dv3_ogwcol + dtf *dv_ogwcol
          du3_oblcol = du3_oblcol + dtf *du_oblcol
          dv3_oblcol = dv3_oblcol + dtf *dv_oblcol

          dws3dt_oss = dws3dt_oss + dtf *sqrt(dudt_oss**2+dvdt_oss**2)
          dws3dt_ofd = dws3dt_ofd + dtf *sqrt(dudt_ofd**2+dvdt_ofd**2)

          ldu3dt_ogw = ldu3dt_ogw + dtf *dudt_ogw
          ldu3dt_obl = ldu3dt_obl + dtf *dudt_obl
          ldu3dt_oss = ldu3dt_oss + dtf *dudt_oss
          ldu3dt_ofd = ldu3dt_ofd + dtf *dudt_ofd

          du3_osscol = du3_osscol + dtf*du_osscol
          dv3_osscol = dv3_osscol + dtf*dv_osscol
          du3_ofdcol = du3_ofdcol + dtf*du_ofdcol
          dv3_ofdcol = dv3_ofdcol + dtf*dv_ofdcol

          ! Special treatment for non-stationary GWD diagnostics
          ldu3dt_ngw = ldu3dt_ngw + dtf *dudt_ngw
          ldv3dt_ngw = ldv3dt_ngw + dtf *dvdt_ngw
          ldt3dt_ngw = ldt3dt_ngw + dtf *dtdt_ngw
        endif
	
!=====================================================================
! Updates inside the ugwpv1_gsldrag.F90
!
!        dtdt = dtdt + dtdt_gw
!        dudt = dudt + dudt_gw
!        dvdt = dvdt + dvdt_gw
!
!       "post" may  also create the "time-averaged" diagnostics"
!            
!     if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
!        do k=1,levs
!          do i=1,im
!             ldu3dt_ngw(i,k) = ldu3dt_ngw(i,k) + dudt_ngw(i,k)*dtf
!             ldv3dt_ngw(i,k) = ldv3dt_ngw(i,k) + dvdt_ngw(i,k)*dtf
!             ldt3dt_ngw(i,k) = ldt3dt_ngw(i,k) + dtdt_ngw(i,k)*dtf
!	  
!             ldu3dt_ogw(i,k) = ldu3dt_ogw(i,k) + dudt_ogw(i,k)*dtf
!             ldv3dt_ogw(i,k) = ldv3dt_ogw(i,k) + dvdt_ogw(i,k)*dtf
!             ldt3dt_ogw(i,k) = ldt3dt_ogw(i,k) + dtdt_ogw(i,k)*dtf
!          enddo
!        enddo
!      endif
! 
!=====================================================================
      end subroutine ugwpv1_gsldrag_post_run      

!! @}
end module ugwpv1_gsldrag_post

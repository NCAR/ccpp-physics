!>  \file cs_conv_aw_adj.F90
!! This file contains a subroutine to adjusts surface rainrate for conservation for CSAW.  

!>\defgroup mod_cs_conv_aw_adj CSAW adjustment Module
!! This module adjusts surface rainrate for conservation.
!> @{
module cs_conv_aw_adj

   implicit none

   private

   public :: cs_conv_aw_adj_run

   contains

!>\ingroup cs_scheme
!> This subroutine adjusts surface rainrate for conservation.
!> \section arg_table_cs_conv_aw_adj_run Argument Table
!! \htmlinclude cs_conv_aw_adj_run.html
!!
!\section gen_cs_conv_aw_adj_run CPT cs_conv_aw_adj_run General Algorithm
   subroutine cs_conv_aw_adj_run(im, levs, do_cscnv, do_aw, do_shoc, &
                ntrac, ntcw, ntclamt, nncl, con_g, sigmafrac,        &
                gt0, gq0, save_t, save_q, prsi, cldfrac, subcldfrac, &
                prcp, imp_physics, imp_physics_mg, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

! --- interface variables
      integer,                                    intent(in)    :: im, levs
      logical,                                    intent(in)    :: do_cscnv, do_aw, do_shoc
      integer,                                    intent(in)    :: ntrac, ntcw, ntclamt, nncl
      real(kind_phys),                            intent(in)    :: con_g
      real(kind_phys),  dimension(:,:),       intent(inout) :: sigmafrac
      real(kind_phys),  dimension(:,:),       intent(inout) :: gt0
      real(kind_phys),  dimension(:,:,:), intent(inout) :: gq0
      real(kind_phys),  dimension(:,:),       intent(in)    :: save_t
      real(kind_phys),  dimension(:,:,:), intent(in)    :: save_q
      real(kind_phys),  dimension(:,:),     intent(in)    :: prsi
      real(kind_phys),  dimension(:,:),       intent(inout), optional :: cldfrac
      real(kind_phys),  dimension(:,:),       intent(inout), optional :: subcldfrac
      real(kind_phys),  dimension(:),            intent(inout) :: prcp
      integer,                                    intent(in   ) :: imp_physics, imp_physics_mg
      character(len=*),                           intent(  out) :: errmsg
      integer,                                    intent(  out) :: errflg

! --- local variables
      real(kind=kind_phys), dimension(im) :: temrain1
      real(kind=kind_phys)                :: tem1, tem2
      real(kind=kind_phys)                :: onebg
      integer :: i, k, n

! --- initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !if (do_cscnv .and. do_aw) then

      onebg = 1.0_kind_phys/con_g

!  Arakawa-Wu adjustment of large-scale microphysics tendencies:
!  reduce by factor of (1-sigma)
!  these are microphysics increments. We want to keep (1-sigma) of the increment,
!  we will remove sigma*increment from final values
!         fsigma = 0.  ! don't apply any AW correction, in addition comment next line
!         fsigma = sigmafrac

!  adjust sfc rainrate for conservation
!  vertically integrate reduction of water increments, reduce precip by that amount

      temrain1(:) = 0.0
      do k = 1,levs
        do i = 1,im
          tem1        = sigmafrac(i,k)
          gt0(i,k)    = gt0(i,k) - tem1 * (gt0(i,k)-save_t(i,k))
          tem2        = tem1 * (gq0(i,k,1)-save_q(i,k,1))
          gq0(i,k,1)  = gq0(i,k,1) - tem2
          temrain1(i) = temrain1(i) - (prsi(i,k)-prsi(i,k+1)) * tem2 * onebg
        enddo
      enddo
! add convective clouds if shoc is true and not MG microphysics
      if (do_shoc .and. imp_physics /= imp_physics_mg) then
        do k = 1,levs
          do i = 1,im
            subcldfrac(i,k) = min(1.0, subcldfrac(i,k) + sigmafrac(i,k))
          enddo
        enddo
      endif
      !
      do n=ntcw,ntcw+nncl-1
        do k = 1,levs
          do i = 1,im
            tem1        = sigmafrac(i,k) * (gq0(i,k,n)-save_q(i,k,n))
            gq0(i,k,n)  = gq0(i,k,n) - tem1
            temrain1(i) = temrain1(i) - (prsi(i,k)-prsi(i,k+1)) * tem1 * onebg
          enddo
        enddo
      enddo
      !
      do i = 1,im
        prcp(i) = max(prcp(i) - temrain1(i)*0.001, 0.0_kind_phys)
      enddo

      !endif

      return

   end subroutine cs_conv_aw_adj_run


end module cs_conv_aw_adj
!> @}

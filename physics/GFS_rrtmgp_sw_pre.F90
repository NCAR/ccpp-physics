module GFS_rrtmgp_sw_pre
  use machine,               only: kind_phys
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use rrtmgp_sw_gas_optics,  only: sw_gas_props

  public GFS_rrtmgp_sw_pre_run, GFS_rrtmgp_sw_pre_init, GFS_rrtmgp_sw_pre_finalize
contains

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_init
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_pre_init ()
  end subroutine GFS_rrtmgp_sw_pre_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_sw_pre_run
!! \htmlinclude GFS_rrtmgp_sw_pre.html
!!
  subroutine GFS_rrtmgp_sw_pre_run(nCol, doSWrad, coszen, nday, idxday, sfc_alb_nir_dir,    &
       sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, sfc_alb_nir_dir_byband,       &
       sfc_alb_nir_dif_byband, sfc_alb_uvvis_dir_byband, sfc_alb_uvvis_dif_byband, errmsg,  &
       errflg)

    ! Input
    integer, intent(in)    :: &
         nCol                 ! Number of horizontal grid points
    logical,intent(in) :: &
         doSWrad              ! Call RRTMGP SW radiation?
    real(kind_phys), dimension(:), intent(in) :: &
         coszen
    real(kind_phys), dimension(:), intent(in) :: &
         sfc_alb_nir_dir,   & !
         sfc_alb_nir_dif,   & !
         sfc_alb_uvvis_dir, & !
         sfc_alb_uvvis_dif    !

    ! Outputs
    integer, intent(out)   :: &
         nday                 ! Number of daylit points
    integer, dimension(:), intent(out) :: &
         idxday               ! Indices for daylit points
    real(kind_phys), dimension(:,:), intent(out) :: &
         sfc_alb_nir_dir_byband,   & ! Surface albedo (direct)
         sfc_alb_nir_dif_byband,   & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir_byband, & ! Surface albedo (direct)
         sfc_alb_uvvis_dif_byband    ! Surface albedo (diffuse)
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &
         errflg               ! Error flag

    ! Local variables
    integer :: i, iBand

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (doSWrad) then
       ! ####################################################################################
       ! For SW gather daylit points
       ! ####################################################################################
       nday   = 0
       idxday = 0
       do i = 1, nCol
          if (coszen(i) >= 0.0001) then
             nday = nday + 1
             idxday(nday) = i
          endif
       enddo

       ! Spread across all SW bands
       do iBand=1,sw_gas_props%get_nband()
          sfc_alb_nir_dir_byband(iBand,1:nCol)   = sfc_alb_nir_dir(1:nCol)
          sfc_alb_nir_dif_byband(iBand,1:nCol)   = sfc_alb_nir_dif(1:nCol)
          sfc_alb_uvvis_dir_byband(iBand,1:nCol) = sfc_alb_uvvis_dir(1:nCol)
          sfc_alb_uvvis_dif_byband(iBand,1:nCol) = sfc_alb_uvvis_dif(1:nCol)
       enddo
    else
       nday                               = 0
       idxday                             = 0
       sfc_alb_nir_dir_byband(:,1:nCol)   = 0.
       sfc_alb_nir_dif_byband(:,1:nCol)   = 0.
       sfc_alb_uvvis_dir_byband(:,1:nCol) = 0.
       sfc_alb_uvvis_dif_byband(:,1:nCol) = 0.
    endif

  end subroutine GFS_rrtmgp_sw_pre_run

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_pre_finalize ()
  end subroutine GFS_rrtmgp_sw_pre_finalize

end module GFS_rrtmgp_sw_pre

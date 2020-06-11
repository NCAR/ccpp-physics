! ###########################################################################################
! ###########################################################################################
module rrtmgp_lwadj
  use machine,    only: kind_phys
  use rrtmgp_aux, only: check_error_msg
  implicit none

  logical :: &
       linit_mod  = .false. !

  public rrtmgp_lwadj_init, rrtmgp_lwadj_run, rrtmgp_lwadj_finalize
contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lwadj_init
  ! #########################################################################################
  subroutine rrtmgp_lwadj_init()
  end subroutine rrtmgp_lwadj_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lwadj_run
  ! #########################################################################################
  subroutine rrtmgp_lwadj_run(use_LW_jacobian, nCol, nLev, skt, sktp1r, fluxlwUP_jac,       &
      fluxlwDOWN_jac, adjsfculw, adjsfcdlw, errmsg, errflg)

  	! Inputs
  	logical, intent(in) :: &
  	     use_LW_jacobian    ! If true the GP scheme is using the Jacobians of the upward/downward
  	                        ! to adjust the LW surface fluxes between radiation calls.
    integer, intent(in) :: &
         nCol,            & ! Number of horizontal gridpoints
         nLev               ! Number of vertical levels 
    real(kind_phys), dimension(nCol), intent(in) :: &
         skt                ! Surface(skin) temperature (K)
    real(kind_phys), dimension(nCol), intent(inout) :: &
         sktp1r             ! Surface(skin) temperature from previous radiation time step (K)
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &
         fluxlwUP_jac,    & ! Jacobian of upward LW flux (W/m2/K)
         fluxlwDOWN_jac     ! Jacobian of downward LW flux (W/m2/K)

	! Outputs
    character(len=*), intent(out) :: &
         errmsg             ! CCPP error message
    integer, intent(out) :: &
         errflg             ! CCPP error flag
	real(kind_phys), dimension(nCol), intent(out) :: &
		adjsfculw,        & !
		adjsfcdlw           !

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0    

	if (.not. use_LW_jacobian) return    
    
    ! Compute adjustment to the surface flux using Jacobian.
	if(linit_mod) then
       adjsfculw(:) = (skt(:) - sktp1r(:)) * fluxlwUP_jac(:,nLev+1)
       adjsfcdlw(:) = (skt(:) - sktp1r(:)) * fluxlwDOWN_jac(:,nLev+1)
    else
       adjsfculw(:) = 0.
       adjsfcdlw(:) = 0.
       linit_mod    = .true.
    endif
    
    ! Store surface temperature for next iteration
    sktp1r(:) = skt(:)
    
  end subroutine rrtmgp_lwadj_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lwadj_finalize
  ! #########################################################################################
  subroutine rrtmgp_lwadj_finalize()
  end subroutine rrtmgp_lwadj_finalize  

end module rrtmgp_lwadj

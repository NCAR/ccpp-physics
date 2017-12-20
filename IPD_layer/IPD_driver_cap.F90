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
!! @brief Semi-auto-generated cap module for the IPD_driver scheme
!!
!
module IPD_driver_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr
    use            :: ccpp_types,                                      &
                      only: ccpp_t
    use            :: ccpp_fields,                                     &
                      only: ccpp_fields_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error
    use            :: IPD_typedefs,                                    &
                      only: IPD_init_type,                             &
                            IPD_control_type,                          &
                            IPD_data_type,                             &
                            IPD_restart_type,                          &
                            IPD_diag_type
    use            :: IPD_driver,                                      &
                      only: IPD_initialize,                            &
                            IPD_setup_step,                            &
                            IPD_radiation_step,                        &
                            IPD_physics_step1,                         &
                            IPD_physics_step2
    use            :: machine,                                         &
                      only: kind_phys
    use            :: namelist_soilveg,                                &
                      only: salp_data, snupx, max_vegtyp
    implicit none

    private
    public :: ipd_initialize_cap,     &
              ipd_setup_step_cap,     &
              ipd_radiation_step_cap, &
              ipd_physics_step1_cap,  &
              ipd_physics_step2_cap

    contains

    subroutine ipd_initialize_cap(ptr) bind(c)

        type(c_ptr), intent(inout) :: ptr

        integer                          :: ierr
        integer, allocatable             :: dims(:)
        type(ccpp_t),           pointer  :: cdata       => null()
        type(IPD_control_type), pointer  :: IPD_Control => null()
        type(IPD_data_type),    pointer  :: IPD_Data(:) => null()
        type(IPD_diag_type),    pointer  :: IPD_Diag(:) => null()
        type(IPD_restart_type), pointer  :: IPD_Restart => null()
        type(IPD_init_type),    pointer  :: Init_parm   => null()
        type(c_ptr)                      :: tmp
        real(kind=kind_phys),   pointer  :: l_snupx(:)  => null()
        real(kind=kind_phys),   pointer  :: l_salp_data => null()

!$OMP critical
        call c_f_pointer(ptr, cdata)

        call ccpp_fields_get(cdata, 'IPD_Control', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Control')
        end if
        call c_f_pointer(tmp, IPD_Control)

        call ccpp_fields_get(cdata, 'IPD_Data', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Data')
        end if
        call c_f_pointer(tmp, IPD_Data, dims)
        deallocate(dims)

        call ccpp_fields_get(cdata, 'IPD_Diag', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Diag')
        end if
        call c_f_pointer(tmp, IPD_Diag, dims)
        deallocate(dims)

        call ccpp_fields_get(cdata, 'IPD_Restart', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Restart')
        end if
        call c_f_pointer(tmp, IPD_Restart)

        call ccpp_fields_get(cdata, 'Init_parm', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Init_parm')
        end if
        call c_f_pointer(tmp, Init_parm)

        call ccpp_fields_get(cdata, 'salp_data', l_salp_data, ierr)
        call ccpp_fields_get(cdata, 'snupx', l_snupx, ierr)
!$OMP end critical

        call IPD_initialize(IPD_Control=IPD_Control, &
                            IPD_Data=IPD_Data,       &
                            IPD_Diag=IPD_Diag,       &
                            IPD_Restart=IPD_Restart, &
                            IPD_init_parm=Init_parm)

        l_snupx = snupx
        l_salp_data = salp_data
    end subroutine ipd_initialize_cap

    subroutine ipd_setup_step_cap(ptr) bind(c)

        type(c_ptr), intent(inout) :: ptr

        integer                          :: ierr
        integer, allocatable             :: dims(:)
        type(ccpp_t),           pointer  :: cdata       => null()
        type(IPD_control_type), pointer  :: IPD_Control => null()
        type(IPD_data_type),    pointer  :: IPD_Data    => null()
        type(IPD_diag_type),    pointer  :: IPD_Diag(:) => null()
        type(IPD_restart_type), pointer  :: IPD_Restart => null()
        type(c_ptr)                      :: tmp

!$OMP critical
        call c_f_pointer(ptr, cdata)

        call ccpp_fields_get(cdata, 'IPD_Control', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Control')
        end if
        call c_f_pointer(tmp, IPD_Control)

        call ccpp_fields_get(cdata, 'IPD_Data', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Data')
        end if
        call c_f_pointer(tmp, IPD_Data)

        call ccpp_fields_get(cdata, 'IPD_Diag', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Diag')
        end if
        call c_f_pointer(tmp, IPD_Diag, dims)
        deallocate(dims)

        call ccpp_fields_get(cdata, 'IPD_Restart', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Restart')
        end if
        call c_f_pointer(tmp, IPD_Restart)
!$OMP end critical

        call IPD_setup_step(IPD_Control=IPD_Control, &
                            IPD_Data=IPD_Data,       &
                            IPD_Diag=IPD_Diag,       &
                            IPD_Restart=IPD_Restart)

    end subroutine IPD_setup_step_cap

    subroutine ipd_radiation_step_cap(ptr) bind(c)

        type(c_ptr), intent(inout) :: ptr

        integer                          :: ierr
        integer, allocatable             :: dims(:)
        type(ccpp_t),           pointer  :: cdata       => null()
        type(IPD_control_type), pointer  :: IPD_Control => null()
        type(IPD_data_type),    pointer  :: IPD_Data    => null()
        type(IPD_diag_type),    pointer  :: IPD_Diag(:) => null()
        type(IPD_restart_type), pointer  :: IPD_Restart => null()
        type(c_ptr)                      :: tmp
        integer,                external :: omp_get_thread_num

!$OMP critical
        call c_f_pointer(ptr, cdata)

        call ccpp_fields_get(cdata, 'IPD_Control', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Control')
        end if
        call c_f_pointer(tmp, IPD_Control)

        call ccpp_fields_get(cdata, 'IPD_Data', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Data')
        end if
        call c_f_pointer(tmp, IPD_Data)

        call ccpp_fields_get(cdata, 'IPD_Diag', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Diag')
        end if
        call c_f_pointer(tmp, IPD_Diag, dims)
        deallocate(dims)

        call ccpp_fields_get(cdata, 'IPD_Restart', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Restart')
        end if
        call c_f_pointer(tmp, IPD_Restart)
!$OMP end critical

        call IPD_radiation_step(IPD_Control=IPD_Control, &
                                IPD_Data=IPD_Data,       &
                                IPD_Diag=IPD_Diag,       &
                                IPD_Restart=IPD_Restart)

    end subroutine ipd_radiation_step_cap

    subroutine ipd_physics_step1_cap(ptr) bind(c)

        type(c_ptr), intent(inout)       :: ptr

        integer                          :: ierr
        integer, allocatable             :: dims(:)
        type(ccpp_t),           pointer  :: cdata       => null()
        type(IPD_control_type), pointer  :: IPD_Control => null()
        type(IPD_data_type),    pointer  :: IPD_Data    => null()
        type(IPD_diag_type),    pointer  :: IPD_Diag(:) => null()
        type(IPD_restart_type), pointer  :: IPD_Restart => null()
        type(c_ptr)                      :: tmp

!$OMP critical
        call c_f_pointer(ptr, cdata)

        call ccpp_fields_get(cdata, 'IPD_Control', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Control')
        end if
        call c_f_pointer(tmp, IPD_Control)

        call ccpp_fields_get(cdata, 'IPD_Data', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Data')
        end if
        call c_f_pointer(tmp, IPD_Data)

        call ccpp_fields_get(cdata, 'IPD_Diag', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Diag')
        end if
        call c_f_pointer(tmp, IPD_Diag, dims)
        deallocate(dims)

        call ccpp_fields_get(cdata, 'IPD_Restart', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Restart')
        end if
        call c_f_pointer(tmp, IPD_Restart)
!$OMP end critical

        call IPD_physics_step1(IPD_Control=IPD_Control, &
                               IPD_Data=IPD_Data,       &
                               IPD_Diag=IPD_Diag,       &
                               IPD_Restart=IPD_Restart)

    end subroutine ipd_physics_step1_cap

    subroutine ipd_physics_step2_cap(ptr) bind(c)

        type(c_ptr), intent(inout) :: ptr

        integer                          :: ierr
        integer, allocatable             :: dims(:)
        type(ccpp_t),           pointer  :: cdata       => null()
        type(IPD_control_type), pointer  :: IPD_Control => null()
        type(IPD_data_type),    pointer  :: IPD_Data    => null()
        type(IPD_diag_type),    pointer  :: IPD_Diag(:) => null()
        type(IPD_restart_type), pointer  :: IPD_Restart => null()
        type(c_ptr)                      :: tmp

!$OMP critical
        call c_f_pointer(ptr, cdata)

        call ccpp_fields_get(cdata, 'IPD_Control', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Control')
        end if
        call c_f_pointer(tmp, IPD_Control)

        call ccpp_fields_get(cdata, 'IPD_Data', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Data')
        end if
        call c_f_pointer(tmp, IPD_Data)

        call ccpp_fields_get(cdata, 'IPD_Diag', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Diag')
        end if
        call c_f_pointer(tmp, IPD_Diag, dims)
        deallocate(dims)

        call ccpp_fields_get(cdata, 'IPD_Restart', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Restart')
        end if
        call c_f_pointer(tmp, IPD_Restart)
!$OMP end critical

        call IPD_physics_step2(IPD_Control=IPD_Control, &
                               IPD_Data=IPD_Data,       &
                               IPD_Diag=IPD_Diag,       &
                               IPD_Restart=IPD_Restart)

    end subroutine ipd_physics_step2_cap

end module IPD_driver_cap

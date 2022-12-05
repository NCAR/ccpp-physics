       module canopy_utils_mod

        implicit none


       contains

!--------------------------------------------------------------------------

       function IntegrateTrapezoid(x, y)
           !! Calculates the integral of an array y with respect to x
           !using the trapezoid
           !! approximation. Note that the mesh spacing of x does not
           !have to be uniform.
           real, intent(in)  :: x(:)                !! Variable x
           real, intent(in)  :: y(size(x))          !! Function y(x)
           real              :: IntegrateTrapezoid  !! Integral of y(x)·dx
       ! Integrate using the trapezoidal rule
           associate(n => size(x))
             IntegrateTrapezoid = sum((y(1+1:n-0) + y(1+0:n-1))*
     &                            (x(1+1:n-0) - x(1+0:n-1)))/2
           end associate
       end function

! ---------------------------------------------------------------------------

      function interp_linear1_internal(x,y,xout) result(yout)
        !! Interpolates for the y value at the desired x value,
        !! given x and y values around the desired point.

        implicit none

        real, intent(IN)  :: x(2), y(2), xout
        real :: yout
        real :: alph

        if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
            write(*,*) "interp1: xout < x0 or xout > x1 !"
            write(*,*) "xout = ",xout
            write(*,*) "x0   = ",x(1)
            write(*,*) "x1   = ",x(2)
            stop
        end if

        alph = (xout - x(1)) / (x(2) - x(1))
        yout = y(1) + alph*(y(2) - y(1))

        return

       end function interp_linear1_internal 

       end module canopy_utils_mod

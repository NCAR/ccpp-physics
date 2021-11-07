!>\file surface_perturbation.F90
!! This file includes routines used in the percentile matching algorithm for the
!! albedo and vegetation fraction perturbations.

!>\defgroup gfs_sfcpert GFS Surface Perturbation Module
!> This module contains routines used in the percentile matching algorithm for the
!! albedo and vegetation fraction perturbations.
module surface_perturbation

      implicit none

      private

      public cdfnor, ppfbet

   contains

! mg, sfc-perts ***

! the routines below are used in the percentile matching algorithm for the
! albedo and vegetation fraction perturbations 

!>\ingroup gfs_sfcpert
!> This subrtouine calculates the CDF of the standard normal distribution
!! evaluated at z.
      subroutine cdfnor(z,cdfz)
      use machine

      implicit none
      real(kind=kind_phys), intent(out) :: cdfz
      real(kind=kind_phys),intent(in) ::  z
! local vars
      integer              iflag
      real(kind=kind_phys) del,x,cdfx,eps

      eps = 1.0E-5


      ! definition of passed parameters !
      ! z  = value for which the normal CDF is to be computed
      ! eps = the absolute accuracy requirment for the CDF
      ! iflag = error indicator on output 0->no errors, 1->errorflag from
      ! cdfgam, 2->errorflag from cdfgam
      ! cdfz = the CDF of the standard normal distribution evaluated at z

        del = 2.0*eps
        if (z.eq.0.0) then
          cdfz = 0.5
        else
          x = 0.5*z*z
          call cdfgam(x,0.5_kind_phys,del,iflag, cdfx)
          if (iflag.ne.0) return
          if (z.gt.0.0) then
            cdfz = 0.5+0.5*cdfx
          else
            cdfz = 0.5-0.5*cdfx
          endif
        endif

      return
      end

!>\ingroup gfs_sfcpert
      subroutine cdfgam(x,alpha,eps,iflag,cdfx)
      use machine

      implicit none
      real(kind=kind_phys), intent(out) :: cdfx
      real(kind=kind_phys),intent(in) ::  x, alpha, eps
! local vars
      integer              iflag,i,j,k, imax
      logical              LL
      real(kind=kind_phys) dx, dgln, p,u,epsx,pdfl, eta, bl, uflo
      data imax, uflo / 5000, 1.0E-37 /


      ! definition of passed parameters !
      ! x  = value for which the CDF is to be computed
      ! alpha = parameter of gamma function (>0)
      ! eps = the absolute accuracy requirment for the CDF
      ! iflag = error indicator on output 0->no errors, 1->either alpha or eps
      ! is <= oflo, 2->number of terms evaluated in the infinite series exceeds
      ! imax.
      ! cdf = the CDF evaluated at x

        cdfx = 0.0

        if (alpha.le.uflo.or.eps.le.uflo) then
          iflag=1
          return
        endif
        iflag=0

        ! check for special case of x
        if (x.le.0) return

        dx = x
        call dgamln(alpha,dgln)
        pdfl = (alpha-1.0)*log(dx)-dx-dgln
        if (pdfl.lt.log(uflo)) then
          if (x.ge.alpha) cdfx = 1.0
        else
          p = alpha
          u = exp(pdfl)
          LL = .true.
          if (x.ge.p) then
            k = int(p)
            if (p.le.real(k)) k = k-1
            eta = p - real(k)
            call dgamln(eta,dgln)
            bl = (eta-1)*log(dx)-dx-dgln
            LL = bl.gt.log(eps)
          endif
          epsx = eps/x
          if (LL) then
            do i=0,imax
              if (u.le.epsx*(p-x)) return
              u = x*u/p
              cdfx = cdfx+u
              p = p+1.0
            enddo
            iflag = 2
          else
            do j=1,k
              p=p-1.0
              if (u.le.epsx*(x-p)) continue
              cdfx = cdfx+u
              u = p*u/x
            enddo
            cdfx = 1.0-cdfx
          endif
        endif
        return
      end subroutine cdfgam

!>\ingroup gfs_sfcpert
      subroutine dgamln(x,dgamlnout)

      use machine
      implicit none
      real(kind=kind_phys), intent(in) ::  x
      real(kind=kind_phys), intent(out) ::  dgamlnout
! local vars
      integer              i, n
      real(kind=kind_phys) absacc, b1, b2, b3, b4, b5, b6, b7, b8
      real(kind=kind_phys) c, dx, q, r, xmin, xn
      data xmin, absacc / 6.894d0, 1.0E-15 /
      data c / 0.918938533204672741780329736d0 /
      data b1 / 0.833333333333333333333333333d-1 /
      data b2 / - 0.277777777777777777777777778d-2 /
      data b3 / 0.793650793650793650793650794d-3 /
      data b4 / - 0.595238095238095238095238095d-3 /
      data b5 / 0.841750841750841750841750842d-3 /
      data b6 / - 0.191752691752691752691752692d-2 /
      data b7 / 0.641025641025641025641025641d-2 /
      data b8 / - 0.295506535947712418300653595d-1 /

      if (x.le.0.0) stop '*** x<=0.0 in function dgamln ***'
      dx = x
      n = max(0,int(xmin - dx + 1.0d0) )
      xn = dx + n
      r = 1.0d0/xn
      q = r*r
      dgamlnout = r*( b1+q*( b2+q*( b3+q*( b4+q*( b5+q*( b6+q*( b7+q*b8 ) ) ) ) ) ) ) +c + (xn-0.5d0)*log(xn)-xn

      if (n.gt.0) then
        q = 1.0d0
        do i=0, n-1
          q = q*(dx+i)
        enddo
        dgamlnout = dgamlnout-log(q)
      endif

      if (dgamlnout + absacc.eq.dgamlnout) then
        print *,' ********* WARNING FROM FUNCTION DGAMLN *********'
        print *,' REQUIRED ABSOLUTE ACCURACY NOT ATTAINED FOR X = ',x
      endif
      return
      end subroutine dgamln

!>\ingroup gfs_sfcpert
!> This subroutine computes the beta distribution value that
!! matches the percentile from the random pattern.
      subroutine ppfbet(pr,p,q,iflag,x)
      use machine
        implicit none
        real(kind=kind_phys), intent(in) :: pr, p, q
        real(kind=kind_phys), intent(out) :: x
        ! local variables
        integer         iflag, iter, itmax
        real(kind=kind_phys)            tol, a, b, fa, fb, fc, cdf, tol1
        real(kind=kind_phys)            c, d, e, xm, s, u, v, r, eps
        data    itmax, eps / 50, 1.0E-12 /

        ! Compute beta distribution value corresponding to the
        ! probability and distribution parameters a,b.
        !
        ! pr - a probability value in the interval [0,1]
        ! p  - the first parameter of the beta(p,q) distribution
        ! q  - the second parameter of the beta(p,q) distribution
        ! iflag - erro indicator in output, 0-no errors, 1,2-error flags
        !         from subroutine cdfbet, 3- pr<0 or pr>1, 4-p<=0 or
        !         q<=0, 5-tol<1.E-8, 6-the cdfs at the endpoints have
        !         the same sign and no value of x is defined, 7-maximum
        !         iterations exceeded and current value of x returned

        tol = 1.0E-5


        iflag = 0
        if (pr.lt.0.0.or.pr.gt.1.) then
          iflag = 3
          return
        endif
        if(min(p,q).le.0.) then
          iflag =4
          return
        endif
        if (tol.lt.1.0E-8) then
          iflag = 5
          return
        endif
        a = 0.
        b = 1.
        fa = -pr
        fb = 1.-pr
        if (fb*fa.gt.0.0) then
          iflag = 6
          return
        endif

        fc = fb
        do iter =1,itmax
          if (fb*fc.gt.0.) then
            c=a
            fc=fa
            d = b-a
            e=d
          endif
          if (abs(fc).lt.abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
          endif

          tol1 = 2.*eps*abs(b)+0.5*tol
          xm = 0.5*(c-b)
          if (abs(xm).le.tol1.or.fb.eq.0.0) then
            x=b
            return
          endif
          if (abs(e).ge.tol1.and.abs(fa).gt.abs(fb)) then
            s = fb/fa
            if (a.eq.c) then
              u = 2.0*xm*s
              v = 1.0-s
            else
              v = fa/fc
              r = fb/fc
              u = s*(2.0*xm*v*(v-r)-(b-a)*(r-1.0))
              v = (v-1.0)*(r-1.0)*(s-1.0)
            endif
            if (u.gt.0.0) v = -v
            u = abs(u)
            if (2.0*u.lt.min(3.0*xm*v-ABS(tol1*v),ABS(e*v))) then
              e = d
              d = u/v
            else
              d = xm
              e = d
            endif

          else

            d=xm
            e=d
          endif

          a = b
          fa = fb
          if (abs(d).gt.tol1) then
            b = b+d
          else
            b = b+sign(tol1,xm)
          endif
          call cdfbet(b,p,q,eps,iflag,cdf)
          if (iflag.ne.0) return
          fb = cdf-pr
        enddo
        x = b

        return
      end subroutine ppfbet

!>\ingroup gfs_sfcpert
!> This subroutine computes the value of the cumulative beta distribution
!! at a single point x, given the distribution parameters p,q.
      subroutine cdfbet(x,p,q,eps,iflag,cdfx)
      use machine

        ! Computes the value of the cumulative beta distribution at a
        ! single point x, given the distribution parameters p,q.
        !
        ! x - value at which the CDF is to be computed
        ! p - first parameter of the beta function
        ! q - second parameter of the beta function
        ! eps - desired absolute accuracy

        implicit none
        real(kind=kind_phys), intent(in) :: x, p, q, eps
        real(kind=kind_phys), intent(out) :: cdfx
        ! local vars
        integer         iflag, jmax, j
        logical         LL
        real(kind=kind_phys)            dp, dq, gamln, yxeps, w, uflo
        real(kind=kind_phys)            xy, yx, pq, qp, pdfl, u, r, v
        real(kind=kind_phys)            tmp
        data jmax, w, uflo / 5000, 20.0, 1.0E-30 /

        cdfx = 0.0

        if (p.le.uflo.or.q.le.uflo.or.eps.le.uflo) then
          iflag = 1
        endif
        iflag = 0

        if (x.le.0.0) return
        if (x.ge.1.0) then
           cdfx=1.0
        else
           LL = (p+w).ge.(p+q+2.0*w)*x
           if (LL) then
              xy = x
              yx = 1.-xy
              pq = p
              qp = q
           else
              yx = x
              xy = 1.-yx
              qp = p
              pq = q
           endif

           call gmln(pq,tmp)
           dp = (pq-1.)*log(xy)-tmp
           call gmln(qp,tmp)
           dq = (qp-1.)*log(yx)-tmp
           call gmln(pq+qp,tmp)
           pdfl = tmp+dp+dq

           if (pdfl.ge.log(uflo)) then
              u = exp(pdfl)*xy/pq
              r = xy/yx
              do while (qp.gt.1.) 
                 if (u.le.eps*(1.-(pq+qp)*xy/(pq+1.))) then
                    if (.not.LL) cdfx = 1.-cdfx
                    return
                 endif
                 cdfx = cdfx+u
                 pq = pq+1.
                 qp = qp-1.
                 u = qp*r*u/pq
              enddo
              v = yx*u
              yxeps = yx*eps
              do j = 0, jmax
                 if (v.le.yxeps) then
                    if (.not.LL) cdfx = 1.-cdfx
                    return
                 endif
                 cdfx = cdfx + v
                 pq = pq+1.
                 v = (pq+qp-1.)*xy*v/pq
              enddo
              iflag = 2
           endif
           if (.not.LL) cdfx = 1.-cdfx
        endif

      end subroutine cdfbet

!>\ingroup gfs_sfcpert
!> This subroutine computes the natural logarithm of the gamma distribution.
!! Users can set the absolute accuracy and corresponding xmin.
      subroutine gmln(x,y)
      use machine
      ! Computes the natural logarithm of the gamma distribution. Users
      ! can set the absolute accuracy and corresponding xmin.

      implicit none
      real(kind=kind_phys), intent(in)  ::  x
      real(kind=kind_phys), intent(out) ::  y
! local vars
      integer              i, n
      real(kind=kind_phys) absacc, b1, b2, b3, b4, b5, b6, b7, b8
      real(kind=kind_phys) c, dx, q, r, xmin, xn
!      data xmin, absacc / 6.894d0, 1.0E-15 /
      data xmin, absacc / 1.357d0, 1.0E-3 /
      data c / 0.918938533204672741780329736d0 /
      data b1 / 0.833333333333333333333333333d-1 /
      data b2 / - 0.277777777777777777777777778d-2 /
      data b3 / 0.793650793650793650793650794d-3 /
      data b4 / - 0.595238095238095238095238095d-3 /
      data b5 / 0.841750841750841750841750842d-3 /
      data b6 / - 0.191752691752691752691752692d-2 /
      data b7 / 0.641025641025641025641025641d-2 /
      data b8 / - 0.295506535947712418300653595d-1 /

      if (x.le.0.0) stop '*** x<=0.0 in function gamln ***'
      dx = x
      n = max(0,int(xmin - dx + 1.0d0) )
      xn = dx + n
      r = 1.0d0/xn
      q = r*r
      y = r*( b1+q*( b2+q*( b3+q*( b4+q*( b5+q*( b6+q*( b7+q*b8 )       &
     & )) ) ) ) ) +c + (xn-0.5d0)*log(xn)-xn

      if (n.gt.0) then
        q = 1.0d0
        do i=0, n-1
          q = q*(dx+i)
        enddo
        y = y-log(q)
      endif

      if (y + absacc.eq.y) then
        print *,' ********* WARNING FROM FUNCTION GAMLN *********'
        print *,' REQUIRED ABSOLUTE ACCURACY NOT ATTAINED FOR X = ',x
      endif
      return
      end subroutine gmln

! *** mg, sfc perts
end module surface_perturbation

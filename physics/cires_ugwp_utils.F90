!  
      subroutine um_flow(nz, klow, ktop, up, vp, tp, qp, dp, zpm, zpi, &
                         pmid, pint, bn2, uhm, vhm, bn2hm, rhohm)
!
      use ugwp_common, only : bnv2min, grav, gocp, fv, rdi
      implicit none
!
!  mass-averaged variables   between klow-ktop
!
      integer, intent(in)                ::  nz, klow, ktop
      real, dimension(nz),   intent(in)  ::  up, vp, tp, qp, dp, zpm, pmid
      real, dimension(nz+1), intent(in)  ::  pint, zpi
      real, dimension(nz),   intent(out) ::  bn2

      real ::  vtj, rhok, bnv2, rdz   
      real ::  vtkp, vtk, dzp, rhm,dphm

      real, intent(out) :: uhm, vhm, bn2hm, rhohm

      integer :: k
!
           dphm  = 0.0 !pint(k+1)-pint(k))
 
           uhm   =  0.0 ! dphm*u1(k)
           vhm   =  0.0 ! dphm*v1(k) 
           rhm   =  0.0 !
           bn2hm =  0.0 !
!
           do k=klow, ktop
            vtj  = tp(k)  * (1.+fv*qp(k))
            vtk  = vtj
            vtkp = tp(k+1)  * (1.+fv*qp(k+1))
            rhok = rdi * pmid(k) / vtj                           ! density kg/m**3
            rdz  = 1.0 / (zpm(k+1)-zpm(k))
! dry
!            bnv2 = grav * (rdz * ( tp(k+1)-tp(k)) +grcp) /tp(k)
!
! wet
!
            bnv2 = grav * (rdz * ( vtkp- vtk) +gocp) /vtk
!              if (bnv2 < 0) print *, k, bnv2,  ' bnv2 < 0 ', klow, ktop
            bnv2 = max(bnv2, bnv2min )
            dzp    = pint(k+1)-pint(k)

            dphm   = dphm + dzp
            uhm    = uhm  + up(k)*dzp 
            vhm    = vhm  + vp(k)*dzp
            rhm    = rhm  + rhok*dzp          
            bn2hm  = bn2hm + bnv2 * dzp
            bn2(k) = bnv2
           enddo

           uhm   = uhm/dphm
           vhm   = vhm/dphm
           rhm   = rhm/dphm
          bn2hm  = bn2hm/dphm
          rhohm  = rhm/dphm
!	  
!          print *, ' MF-BV ', bn2hm, bn2(ktop), bn2(klow)
!
      end subroutine um_flow
!
!
      subroutine mflow_tauz(levs, up, vp, tp, qp, dp, zpm, zpi, &
        pmid, pint, rho, ui, vi, ti, bn2i, bvi, rhoi)

        use ugwp_common, only : bnv2min, grav, gocp, fv, rdi

        implicit none
        
        integer :: levs
        real, dimension(levs)   ::  up, vp, tp, qp, dp, zpm, pmid
        real, dimension(levs+1) ::  pint, rho, zpi
        real, dimension(levs)   ::  zdelpi, zdelpm
        real                    ::  zul, bvl
        real, dimension(levs+1) ::  ui, vi, bn2i, bvi, rhoi, ti, qi

        real ::  vtj, rhok, bnv2, rdz   
        real ::  vtkp, vtk, dzp
        real ::  vtji
        integer :: k
!
! get interface values from surf to top
!	
       do k=2,levs                                      
         vi(k) = 0.5 *(vp(k-1) + vp(k)) 
         ui(k) = 0.5 *(up(k-1) + up(k))
         ti(k) = 0.5 *(tp(k-1) + tp(k))
         qi(k) = 0.5 *(qp(k-1) + qp(k))
       enddo

                     k=1                                      
        ti(k) = tp(k) 
        ui(k) = up(k)
        vi(k) = vp(k) 
        qi(k) = qp(k)
                     k= levs
        ti(k+1) = tp(k) 
        ui(k+1) = up(k)
        vi(k+1) = vp(k) 
        qi(k+1)=qp(k)

       do k=1,levs-1 
            vtj  = tp(k)   * (1.+fv*qp(k))
            vtji = ti(k)   *  (1.+fv*qi(k))
            rho(k)   = rdi * pmid(k) / vtj                           ! density kg/m**3
            rhoi(k) =  rdi * pint(k) / vtji
            vtk  = vtj
            vtkp = tp(k+1)  * (1.+fv*qp(k+1))
            rdz  = 1.   / ( zpm(k+1)-zpm(k))
            bnv2 = grav * (rdz * ( vtkp- vtk) +gocp) /vtji
            bn2i(k) = max(bnv2, bnv2min )
            bvi(k) = sqrt( bn2i(k) )         
            vtk = vtkp
       enddo
       k = levs
            vtj  = tp(k)   ! * (1.+fv*qp(k))
            vtji = ti(k)   !*  (1.+fv*qi(k))
            rho(k)   = rdi * pmid(k) / vtj                           
            rhoi(k) =  rdi * pint(k) / vtji               
            bn2i(k) =  bn2i(k-1)
            bvi(k) = sqrt( bn2i(k) )
       k = levs+1
            rhoi(k) =  rdi * pint(k) / ti(k)
            bn2i(k) =  bn2i(k-1)
            bvi(k) = sqrt( bn2i(k) )
!       do k=1,levs
!          write(6, 121) k, zpm(k)*1.e-3, zpi(k)*1.e-3, bvi(k), rho(k), rhoi(k)
!       enddo
 121   format(i5, 2x, 3(2x, F10.3), 2(2x, E10.3))      

      end subroutine mflow_tauz

!     
      subroutine get_unit_vector(u, v, u_n, v_n, mag)
      implicit  none
      real, intent(in) :: u, v
      real, intent(out) :: u_n, v_n, mag
!

      mag = sqrt(u*u + v*v)

      if (mag > 0.0) then
         u_n = u/mag
         v_n = v/mag
      else
         u_n = 0.
         v_n = 0.
      end if

      end subroutine get_unit_vector
!

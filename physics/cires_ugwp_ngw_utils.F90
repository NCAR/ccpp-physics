module cires_ugwp_ngw_utils


contains


      subroutine tau_limb_advance(me, master, im, levs, ddd, curdate,   &
	    j1_tau, j2_tau, ddy_j1tau, ddy_j2tau,  tau_sat,            kdt )  
	    



      use machine, only : kind_phys
      
      use cires_ugwp_module_v1, only : ntau_d1y, ntau_d2t       
      use cires_ugwp_module_v1, only : ugwp_taulat, days_limb,  tau_limb
      
!      use cires_ugwp_module, only : ugwp_qbolat,  days_merra, pmb127, days_y4md, days_y4ddd
!      use cires_ugwp_module, only :  tau_qbo,  stau_qbo,  uqboe, u2 => uzmf_merra  

      implicit none

      integer, intent(in) ::    me, master, im, levs, ddd, curdate, kdt    
      integer, intent(in), dimension(im) :: j1_tau, j2_tau
      
      real , intent(in),  dimension(im) :: ddy_j1tau, ddy_j2tau  
           
      real, intent(out) ::  tau_sat(im)
      
      integer           :: i, j1, j2, k, it1, it2, iday
      real              :: tem,  tx1, tx2, w1, w2, day2, day1, ddx
      integer           :: yr1, yr2  
!
      integer           ::  iqbo1=1      
!

	 
	 
            it1 = 2
         do iday=1, ntau_d2t
	    if (float(ddd) .lt. days_limb(iday) ) then
	    it2 = iday
	    exit
	    endif
	 enddo
	 it2 = min(it2,ntau_d2t)	 
	 it1 = max(it2-1,1)
	 if (it2 > ntau_d2t ) then
	  print *, ' it1, it2, ntau_d2t ', it1, it2, ntau_d2t
	  stop
	 endif
	 w2 = (float(ddd)-days_limb(it1))/(days_limb(it2)-days_limb(it1))
	 w1 = 1.0-w2
      do i=1, im	 
	 j1 = j1_tau(i)
	 j2 = j2_tau(i)
	 tx1 = tau_limb(j1, it1)*ddy_j1tau(i)+tau_limb(j2, it1)*ddy_j2tau(i)
	 tx2 = tau_limb(j1, it2)*ddy_j1tau(i)+tau_limb(j2, it2)*ddy_j2tau(i)	 
	 tau_sat(i) =  tx1*w1 + w2*tx2 
      enddo
      
         if (me == master ) then	    
	    print*, maxval(tau_limb), minval(tau_limb), ' tau_limb '
	    print*, ntau_d2t
	    print*, days_limb(1) ,  days_limb(ntau_d2t)	, ddd,  ' days-taulimb '
	    print*, 'curdate  ', curdate
	    print*, maxval(tau_sat), minval(tau_sat), ' tau_sat_fv3 '	        	    
	 endif 
      return

      end subroutine tau_limb_advance

end module cires_ugwp_ngw_utils

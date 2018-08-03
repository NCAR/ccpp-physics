      module four_to_grid_mod

      use stochy_ccpp, only: num_parthds_stochy => ompthreads

      implicit none

      contains

      subroutine four_to_grid(syn_gr_a_1,syn_gr_a_2,
     &                        lon_dim_coef,lon_dim_grid,lons_lat,lot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use machine
      implicit none
!!
      real(kind=kind_dbl_prec)     syn_gr_a_1(lon_dim_coef,lot)
      real(kind=kind_dbl_prec)     syn_gr_a_2(lon_dim_grid,lot)
      integer                  lon_dim_coef
      integer                  lon_dim_grid
      integer                  lons_lat
      integer                  lot
!________________________________________________________
#ifdef MKL
      integer*8                plan
#else
      real(kind=kind_dbl_prec) aux1crs(42002)
      real(kind=kind_dbl_prec) scale_ibm
      integer                  ibmsign
      integer                  init
#endif
      integer                  lot_thread
      integer                  num_threads
      integer                  nvar_thread_max
      integer                  nvar_1
      integer                  nvar_2
      integer                  thread
#ifdef MKL
      include                  "fftw/fftw3.f"
      integer                  NULL
#else
      external                 dcrft
      external                 scrft
#endif
!________________________________________________________
      num_threads     = min(num_parthds_stochy,lot)
 
      nvar_thread_max = (lot+num_threads-1)/num_threads
 
      if ( kind_dbl_prec == 8 ) then !------------------------------------
#ifdef MKL
!$omp parallel do shared(syn_gr_a_1,syn_gr_a_2,lons_lat)
!$omp+shared(lon_dim_coef,lon_dim_grid)
!$omp+shared(lot,num_threads,nvar_thread_max)
!$omp+private(thread,nvar_1,nvar_2,lot_thread,plan)
#else
!$omp parallel do shared(syn_gr_a_1,syn_gr_a_2,lons_lat)
!$omp+shared(lon_dim_coef,lon_dim_grid)
!$omp+shared(lot,num_threads,nvar_thread_max)
!$omp+shared(ibmsign,scale_ibm)
!$omp+private(thread,nvar_1,nvar_2,lot_thread,init,aux1crs)
#endif
         do thread=1,num_threads   ! start of thread loop ..............
           nvar_1=(thread-1)*nvar_thread_max + 1
           nvar_2=min(nvar_1+nvar_thread_max-1,lot)

           lot_thread=nvar_2 - nvar_1 + 1

           if (nvar_2 >= nvar_1) then
#ifdef MKL
             !call dfftw_plan_many_dft_c2r(
             !            plan, 1, N,        m,                         &
             !            X,          NULL, 1, dimx,                    &
             !            Y,          NULL, 1, dimy,                    &
             !            fftw_flag)
             call dfftw_plan_many_dft_c2r(                              &
     &                    plan, 1, lons_lat, lot_thread,                &
     &                    syn_gr_a_1, NULL, 1, size(syn_gr_a_1,dim=1),  &
     &                    syn_gr_a_2, NULL, 1, size(syn_gr_a_2,dim=1),  &
     &                    FFTW_ESTIMATE)
             call dfftw_execute(plan)
             call dfftw_destroy_plan(plan)
#else
             init      = 1
             ibmsign   = -1
             scale_ibm = 1.0d0

             call dcrft(init,
     &               syn_gr_a_1(1,nvar_1)   ,lon_dim_coef/2,
     &               syn_gr_a_2(1,nvar_1)   ,lon_dim_grid,
     &               lons_lat,lot_thread,ibmsign,scale_ibm,
     &               aux1crs,22000,
     &               aux1crs(22001),20000)

             init = 0
             call dcrft(init,
     &               syn_gr_a_1(1,nvar_1)   ,lon_dim_coef/2,
     &               syn_gr_a_2(1,nvar_1)   ,lon_dim_grid,
     &               lons_lat,lot_thread,ibmsign,scale_ibm,
     &               aux1crs,22000,
     &               aux1crs(22001),20000)
#endif
           endif

         enddo  ! fin thread loop ......................................
      else !------------------------------------------------------------
#ifdef MKL
!$omp parallel do shared(syn_gr_a_1,syn_gr_a_2,lons_lat)
!$omp+shared(lon_dim_coef,lon_dim_grid)
!$omp+shared(lot,num_threads,nvar_thread_max)
!$omp+private(thread,nvar_1,nvar_2,lot_thread,plan)
#else
!$omp parallel do shared(syn_gr_a_1,syn_gr_a_2,lons_lat)
!$omp+shared(lon_dim_coef,lon_dim_grid)
!$omp+shared(lot,num_threads,nvar_thread_max)
!$omp+shared(ibmsign,scale_ibm)
!$omp+private(thread,nvar_1,nvar_2,lot_thread,init,aux1crs)
#endif
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1 = (thread-1)*nvar_thread_max + 1
            nvar_2 = min(nvar_1+nvar_thread_max-1,lot)

            lot_thread = nvar_2 - nvar_1 + 1

            if (nvar_2 >= nvar_1) then
#ifdef MKL
              !call sfftw_plan_many_dft_c2r(
              !            plan, 1, N,        m,                        &
              !            X,          NULL, 1, dimx,                   &
              !            Y,          NULL, 1, dimy,                   &
              !            fftw_flag)
              call sfftw_plan_many_dft_c2r(                             &
     &                     plan, 1, lons_lat, lot_thread,               &
     &                     syn_gr_a_1, NULL, 1, size(syn_gr_a_1,dim=1), &
     &                     syn_gr_a_2, NULL, 1, size(syn_gr_a_2,dim=1), &
     &                     FFTW_ESTIMATE)
              call sfftw_execute(plan)
              call sfftw_destroy_plan(plan)
#else
              init      = 1
              ibmsign   = -1
              scale_ibm = 1.0d0
              call scrft(init,
     &                syn_gr_a_1(1,nvar_1)   ,lon_dim_coef/2,
     &                syn_gr_a_2(1,nvar_1)   ,lon_dim_grid,
     &                lons_lat,lot_thread,ibmsign,scale_ibm,
     &                aux1crs,22000,
     &                aux1crs(22001),20000,
     &                aux1crs(22001),0)
              init = 0
              call scrft(init,
     &                syn_gr_a_1(1,nvar_1)   ,lon_dim_coef/2,
     &                syn_gr_a_2(1,nvar_1)   ,lon_dim_grid,
     &                lons_lat,lot_thread,ibmsign,scale_ibm,
     &                aux1crs,22000,
     &                aux1crs(22001),20000,
     &                aux1crs(22001),0)
#endif
           endif
         enddo  ! fin thread loop ......................................
      endif !-----------------------------------------------------------
!!
      return
      end
      subroutine grid_to_four(anl_gr_a_2,anl_gr_a_1,
     &                        lon_dim_grid,lon_dim_coef,lons_lat,lot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use machine
      implicit none
!!
      real(kind=kind_dbl_prec)     anl_gr_a_2(lon_dim_grid,lot)
      real(kind=kind_dbl_prec)     anl_gr_a_1(lon_dim_coef,lot)
      integer                  lon_dim_grid
      integer                  lon_dim_coef
      integer                  lons_lat
      integer                  lot
!________________________________________________________
      real(kind=kind_dbl_prec) aux1crs(42002)
      real(kind=kind_dbl_prec)     scale_ibm,rone
      integer                  ibmsign
      integer                  init
      integer                  lot_thread
      integer                  num_threads
      integer                  nvar_thread_max
      integer                  nvar_1,nvar_2
      integer                  thread
!________________________________________________________
#ifdef MKL
      write(0,*) "ERROR in grid_to_four: srcft and drcft ",
     & " must be replaced with MKL's FFTW calls. ABORT."
      call sleep(5)
      stop
#endif
      num_threads=min(num_parthds_stochy,lot)

      nvar_thread_max=(lot+num_threads-1)/num_threads
 
      if ( kind_dbl_prec == 8 ) then !------------------------------------
!$omp parallel do shared(anl_gr_a_1,anl_gr_a_2,lons_lat)
!$omp+shared(lon_dim_coef,lon_dim_grid)
!$omp+shared(lot,num_threads,nvar_thread_max)
!$omp+shared(ibmsign,scale_ibm,rone)
!$omp+private(thread,nvar_1,nvar_2,lot_thread,init,aux1crs)
 
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1 = (thread-1)*nvar_thread_max + 1
            nvar_2 = min(nvar_1+nvar_thread_max-1,lot)

            if (nvar_2 >= nvar_1) then
              lot_thread = nvar_2 - nvar_1 + 1
 
              init      = 1
              ibmsign   = 1
              rone      = 1.0d0
              scale_ibm = rone/lons_lat
              call drcft(init,
     &                anl_gr_a_2(1,nvar_1),   lon_dim_grid,
     &                anl_gr_a_1(1,nvar_1),   lon_dim_coef/2,
     &                lons_lat,lot_thread,ibmsign,scale_ibm,
     &                aux1crs,22000,
     &                aux1crs(22001),20000)
              init = 0
              call drcft(init,
     &                anl_gr_a_2(1,nvar_1),   lon_dim_grid,
     &                anl_gr_a_1(1,nvar_1),   lon_dim_coef/2,
     &                lons_lat,lot_thread,ibmsign,scale_ibm,
     &                aux1crs,22000,
     &                aux1crs(22001),20000)
 
            endif
         enddo  ! fin thread loop ......................................
      else !------------------------------------------------------------
!$omp parallel do shared(anl_gr_a_1,anl_gr_a_2,lons_lat)
!$omp+shared(lon_dim_coef,lon_dim_grid)
!$omp+shared(lot,num_threads,nvar_thread_max)
!$omp+shared(ibmsign,scale_ibm,rone)
!$omp+private(thread,nvar_1,nvar_2,lot_thread,init,aux1crs)
 
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1 = (thread-1)*nvar_thread_max + 1
            nvar_2 = min(nvar_1+nvar_thread_max-1,lot)

            if (nvar_2 >= nvar_1) then
              lot_thread=nvar_2 - nvar_1 + 1
 
              init      = 1
              ibmsign   = 1
              rone      = 1.0d0
              scale_ibm = rone/lons_lat
              call srcft(init,
     &                anl_gr_a_2(1,nvar_1),   lon_dim_grid,
     &                anl_gr_a_1(1,nvar_1),   lon_dim_coef/2,
     &                lons_lat,lot_thread,ibmsign,scale_ibm,
     &                aux1crs,22000,
     &                aux1crs(22001),20000,
     &                aux1crs(22001),0)
              init = 0
              call srcft(init,
     &                anl_gr_a_2(1,nvar_1),   lon_dim_grid,
     &                anl_gr_a_1(1,nvar_1),   lon_dim_coef/2,
     &                lons_lat,lot_thread,ibmsign,scale_ibm,
     &                aux1crs,22000,
     &                aux1crs(22001),20000,
     &                aux1crs(22001),0)
 
            endif
         enddo  ! fin thread loop ......................................
      endif !-----------------------------------------------------------
!!
      return
      end

      end module four_to_grid_mod

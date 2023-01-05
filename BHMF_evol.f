      PROGRAM integ

      implicit none
      integer, parameter :: N=3000
      integer :: i, j, k
      real*8 :: time, Phi_M(N-1), z, M_BH(N), Mmin, Mmax, eps,
     &     Phi_0, M_c, alpha, beta, x,  dM_BH(N-1),
     &     lambda, dt, tau, t_s, delta, MM,
     &     Phi_M_prev(N-1), deriv, lambda1, eps1, dPdlnl, dP,
     &     Number, a, l_c, gammp, consv,
     &     lmin, lmax, M_BH_cen(N-1), xmin=0.1

!     set the grid 
      Mmin = 1.0d5
      Mmax = 1.0d13
      eps  = (Mmax/Mmin)**(1.0d0/dble(N-1))-1.0d0

      print*, eps, 'eps'

      do i=1,N
         M_BH(i)=Mmin*(1.0d0+eps)**(i-1)
      enddo
      do i=1,N-1
         dM_BH(i) = M_BH(i+1)-M_BH(i)
         M_BH_cen(i) = M_BH(i) * (M_BH(2)/M_BH(1))**0.5
      enddo

!     set the initial BHMF (z=6)
      Phi_0 = 1310d0
      M_c   = 6.13d7
      alpha = -1.41d0
      beta  = -2.58d0
      do i=1,N-1
         x = M_BH_cen(i)/M_c
         Phi_M(i) = Phi_0/(x**(-alpha-1.0d0) + x**(-beta-1.0d0))
      enddo
      time = 913.0d0            ! Myr

!     Growth model
      tau = 18.76d0             ! in Myr
      t_s = 45.0d0              ! in Myr
      dt  = tau
      delta = 0.0d0
      MM  = 1.0d8

      dt = tau

!====================!
!     Evolve BHMF    !
!====================!
!     time update
      do k=1,10
         consv = 0.
         do i=1,N-1
            Phi_M_prev(i) = Phi_M(i) ! previous time step
         enddo
         do i=1,N-1
            Phi_M(i) = 0.0d0
            lmin = 10000000.
            lmax = 0.
            do j=1,i
               call EDR(M_BH(i),     M_BH_cen(j),delta,MM,lambda)
               call EDR(M_BH(i+1),   M_BH_cen(j),delta,MM,lambda1)
               lambda  = max(lambda * t_s / dt,  xmin)
               lambda1 = max(lambda1* t_s / dt,  xmin)
               
               ! if (xmin>lambda .and. xmin<lambda1) then 
               !    lambda = xmin
               !    ! print*, 'lambda=xmin:',lambda
               ! endif
               if (lmin>lambda .and. lambda>=xmin) lmin = lambda
               if (lmax<lambda1) lmax = lambda1
               call ERDF_intg(lambda,lambda1,dP)

               ! lambda  = 0.01
               ! lambda1 = 10
               ! call ERDF_intg(lambda,lambda1,dP)
               ! print*,'dP', dP, M_BH(i+1)/M_BH(i), eps1
               ! print*, lambda, lambda1, dlog10(M_BH(j+1)/M_BH(j))
               ! stop

               Phi_M(i) = Phi_M(i) + dP !/ dlog10(M_BH(2)/M_BH(1))
     &                     * Phi_M_prev(j)!*dlog10(M_BH_cen(2)/M_BH_cen(1))
               consv = consv + dP
               ! print*, 'consv', consv
            enddo
!  100        continue
         enddo
         print*, 'lmin',lmin, 'lmax',lmax
         print*, 'consv', consv

         time = time + dt

!     write output
         if(mod(k,1)==0) then
            do i=1,N-1
               write(10,*) time, i, M_BH(i), Phi_M(i),  Phi_M_prev(i)
            enddo
         endif

         Number = 0.0d0
         do i=1,N-1
   !          Number= Number + 0.5d0*(Phi_M(i)+Phi_M(i+1))
   !   &           *dlog10(M_BH(i+1)/M_BH(i))
            Number= Number + Phi_M(i)*dlog10(M_BH(i+1)/M_BH(i))
         enddo
         print*, Number
!--------------------
      enddo
      print*, time, 'time'

      END



!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE EDR(x, x0, d, mm, lambda)
      implicit none 
      real*8 :: x,x0,d, lambda, dldx, mm
      
      if(d==0.0d0) then
         lambda = dlog(x/x0)
      else
         lambda = dlog(x/x0) + ((x/mm)**d - (x0/mm)**d)/d
      endif

      return
      END

      SUBROUTINE ERDF(x,P)
      implicit none
      integer :: i, M
      real*8 :: x, x0, a, P, dx, xmin
      real*8, save  :: Sum = 0.0d0
      
      a=0.20d0
      x0=0.87d0
      xmin = 1.0d-2

      Sum = 2.5480
      P = (x/x0)**a * dexp(-x/x0)
      if(x<xmin) P=0.0d0

      if(Sum==0.0d0) then
         M=1000000000
         do i=1,M
            if(x>10.0d0) exit
            dx = xmin*1.0d-5
            x = xmin + dx*(i-1)
            P = (x/x0)**a * dexp(-x/x0)
            Sum = Sum + P*dx/x
         enddo
         print*, Sum, 'Sum' 
!      stop
      endif
      P = P/Sum
      return
      END

      SUBROUTINE ERDF_intg(x1,x2,P)
      implicit none
      integer :: i, M
      real*8 :: x1, x2, x, x0, a, P, dlnx, xmin
      real*8, save  :: Sum = 0.0d0
      
      a=0.20d0
      x0=0.87d0
      xmin = 1.0d-1

      ! Sum = 1.4074868928125490 

      if(x1<xmin) then 
         x=xmin
      else 
         x = x1
      endif
      if(x2<xmin) then
         ! print*, x2, 'x2 wrong' 
         P = 0
         return
      endif

      M=1000000000
      dlnx = 1.0d-3
      P = 0
      do i=1,M
         if(x>=x2) exit
         P = P + (x/x0)**a * dexp(-x/x0) * dlnx
         x = x * dexp(dlnx)
      enddo

      if(Sum==0.0d0) then
         x = xmin
         M=1000000000
         do i=1,M
            if(x>10.0d0) exit
            Sum = Sum + (x/x0)**a * dexp(-x/x0)*dlnx
            x = x * dexp(dlnx)
         enddo
         print*, Sum, 'Sum' 
      ! stop
      endif
      P = P/Sum
      return
      END

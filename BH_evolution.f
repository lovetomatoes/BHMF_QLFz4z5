      PROGRAM BH_MF

      implicit none
      integer, parameter :: N=1000000, N_on=1
      integer :: i,j,k, i_ERDF
      real*8 :: dt, t_0, t_i, t_f
      real*8 :: M_BH(N), x(N), y(N), v(N), w(N), s(N), u(N),
     &     z1(N), z2(N), z3(N), z4(N), z5(N), z6(N), lambda(N), Mdot, 
     &     M1450, L_bol, L_sun
      real*8 :: alpha, beta, aa
      real*8 :: pi, mu, sigma, eta, delta

!     ERDF: 0 (log normal)
!     ERDF: 1 (Schechter)     
      i_ERDF = 1


      L_sun = 3.839d33
      pi = 2.0d0*dacos(0.0d0)
      call random_number(x)
      call random_number(y)
      call random_number(v)
      call random_number(w)
      call random_number(s)
      call random_number(u)

!     log normal
      mu = -1.0d0
      sigma = 0.45d0
!     M-dep efficiency
      delta = 0.3d0
      eta   = 0.1d0
!     Schechter
      alpha = 0.5d0
      beta  = 0.4d0
      

!     time scale in yrs
      t_0 = 4.0d8
      t_i = 3.0d8
      t_f = 8.0d8
      dt = t_f-t_i
      ! write(10,*) 'ls'
      ! do i = 1,N
      !    call random_gamma(alpha,aa)
      !    lambda(i) = aa * beta
      !    write(10,*) lambda(i)
      ! enddo
      ! stop

!     initial condition
      do i=1,N
         M_BH(i) = 1.0d2
         if(i_ERDF==0) then
            z1(i) = dsqrt(-2.0d0*dlog(x(i)))*dcos(2.0d0*pi*y(i))
            z2(i) = dsqrt(-2.0d0*dlog(y(i)))*dcos(2.0d0*pi*x(i))
            z3(i) = dsqrt(-2.0d0*dlog(v(i)))*dcos(2.0d0*pi*w(i))
            z4(i) = dsqrt(-2.0d0*dlog(w(i)))*dcos(2.0d0*pi*v(i))
            z5(i) = dsqrt(-2.0d0*dlog(s(i)))*dcos(2.0d0*pi*u(i))
            z6(i) = dsqrt(-2.0d0*dlog(u(i)))*dcos(2.0d0*pi*s(i))
         endif
c         write(10,*) i, M_BH(i), x(i), z1(i), z2(i)
      enddo

!     BH growth
      write(16,*) 'loglambda ', 'logMBH ', ' M1450'        

      do i=1,N
         do j=1,N_on
            if(i_ERDF==0) then
               if(j==1) lambda(i) = 10.0d0**(z1(i)*sigma + mu)
               if(j==2) lambda(i) = 10.0d0**(z2(i)*sigma + mu)
               if(j==3) lambda(i) = 10.0d0**(z3(i)*sigma + mu)
               if(j==4) lambda(i) = 10.0d0**(z4(i)*sigma + mu)
               if(j==5) lambda(i) = 10.0d0**(z5(i)*sigma + mu)
               if(j==6) lambda(i) = 10.0d0**(z6(i)*sigma + mu)
            endif
            if(i_ERDF==1) then
               call random_gamma(alpha,aa)
               lambda(i) = aa * beta
            endif

            call Direct_integr(M_BH(i),dt/N_on,eta,delta,lambda(i),t_0)
   
            if(i_ERDF==1) then
               call random_gamma(alpha,aa)
               lambda(i) = aa * beta
            endif
       
         enddo

!         call Analytical_BH(M_BH(i), dt, eta, lambda(i), t_0)
!         call Newton_method(M_BH(i), dt, eta, lambda(i), t_0)

         L_bol = 1.3d38*M_BH(i)*lambda(i)
         M1450 = -25.0d0-2.5d0*(dlog10(L_bol/L_sun)-13.15364d0)
         if(M_BH(i)>1.0d2) then
!         if(M1450>-23.0d0) then
            write(16,*) dlog10(lambda(i)), dlog10(M_BH(i)), M1450        
         endif
      enddo

      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Analytical_BH(M0, dt, eta, lambda, t_0)
      implicit none
      integer :: i
      real*8 :: M0, dt, eta, lambda, t_0, M

      M0 = M0 * dexp((1.0d0-eta)/eta*lambda/t_0*dt)
      return
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Newton_method(M0, dt, eta, lambda, t_0)
      implicit none
      integer :: i
      real*8 :: M0, dt, eta, lambda, t_0, M, fM, dfdM

      M = M0
      do i=1,100
         fM   = eta/(1.0d0-eta)/lambda *dlog(M/M0) - dt/t_0
         dfdM = eta/(1.0d0-eta)/lambda /M
         M = M -fM/dfdM
         if(dabs(fM/dfdM/M)<1.0d-2) exit
!         if(i>40) print*, i, M, 'too much iteration'
         if(i==100) print*, M0, lambda, M, 'not converge...'
      enddo
!      if(M>1.0d5) print*, M, M0, lambda
      M0 = M

      return
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Direct_integr(M0, dt, eta0, delta, lambda, t_0)
      implicit none
      integer :: i
      real*8 :: M0, dt, eta, eta0, lambda, delta, 
     &     t_0, M, Mdot, ddt, eps, t

      eps = 3.0d-2

      t = 0.0d0
      M = M0
      do i=1,10000
         eta = eta0*(1.0d0 + (M/1.0d7)**delta)
!         eta = dmin1(eta, 0.42d0)
!         Mdot = (1.0d0-eta)/eta*lambda/t_0*M
         Mdot = 1.0d0/eta*lambda/t_0*M
         ddt = M/Mdot * eps
         M = M + Mdot*ddt
         t = t + ddt
         if(t>dt) exit
      enddo
!      print*, M, M0, t
      M0 = M

      return
      END



      SUBROUTINE random_gamma(alpha,x)
      implicit none
      real*8 :: alpha, x, g, u

      if(alpha<=0) then
         stop "alpha<=0"
      else if(alpha<1.0d0) then
         call random_gamma_alpha_ge_1(alpha+1.0d0,g)
         call random_number(u)
         x = g*u**(1.0d0/alpha)
      else
         call random_gamma_alpha_ge_1(alpha,x)
      end if
      return
      END

      SUBROUTINE random_gamma_alpha_ge_1(alpha,x)
      implicit none
      real*8 :: alpha, x, u1, u2, y, t

      do
         call random_number(u1)
         y = -log(u1)
         t = (y/exp(y-1))**(alpha-1)
         call random_number(u2)
         if(u2 <= t) then
            x = alpha*y
            exit
         end if
      enddo
      return
      END


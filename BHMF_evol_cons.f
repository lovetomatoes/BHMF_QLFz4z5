      PROGRAM integ
!      SUBROUTINE integ(tau, delta, MM, a, l_c, chisq)

!   gfortran BHMF_evol_cons.f -o BHMF_evol_cons.o; ./BHMF_evol_cons.o

      implicit none
      integer, parameter :: N=1000, N_d=18
      integer :: i, j, k, i_read, i_sw
      real*8 :: time, Phi_M(N-1), z, M_BH(N), Mmin, Mmax, eps, 
     &     M_BH_c(N-1), dlogM_BH(N-1), P_M(N-1), P_M_prev(N-1),
     &     Phi_0, M_c, alpha, beta, x, lambda(N), dP_l(N-1),
     &     dt, tau, t_s, delta, MM,
     &     mag(N), mag_c(N-1), Phi_L(N-1), deriv, dP_dlnl, 
     &     magmin, magmax, lambda_b, f_obs(N-1), XX, f_x, 
     &     L_bol, L_sun, L_x,
     &     Number, a, l_c, gammp,
     &     mag_data(N_d), Phi_L_data(N_d), dPhi_L_data(N_d), dumy, chisq

      i_sw=0

!     set the grid 
      Mmin = 1.0d5
      Mmax = 1.0d11
      eps  = (Mmax/Mmin)**(1.0d0/dble(N-1))-1.0d0

      do i=1,N
         M_BH(i) = Mmin*(1.0d0+eps)**(i-1)
      enddo
      do i=1,N-1
         M_BH_c(i) = dsqrt(M_BH(i+1)*M_BH(i))
         dlogM_BH(i) = dlog10(M_BH(i+1)/M_BH(i))
      enddo

      magmin = -15.0d0
      magmax = -30.0d0
      do i=1,N
         mag(i) = magmin + (magmax-magmin)*(i-1)/(N-1)
      enddo
      do i=1,N-1
         mag_c(i) = 0.5d0*(mag(i)+mag(i+1))
      enddo

!     set the initial BHMF (z=6)
      Phi_0 = 1310d0
      M_c   = 6.13d7
      alpha = -1.41d0
      beta  = -2.58d0
      do i=1,N-1
         x = M_BH_c(i)/M_c
         Phi_M(i) = Phi_0/(x**(-alpha-1.0d0) + x**(-beta-1.0d0))
         P_M(i)   = Phi_M(i)*dlogM_BH(i)
      enddo
      time = 935.0d0            ! Myr

!     Growth model
!      tau = 18.76d0             ! in Myr
      tau = 10.368000984191895d0
      t_s = 45.0d0              ! in Myr
      dt  = tau
      delta = 0.01d0
      MM  = 1.0d8

!     ERDF for the given grid set
      a   = 0.89161005020141604d0
      l_c = 0.301d0

      do i=1,N
         call EDR(M_BH(i), M_BH_c(1), delta, MM, lambda(i))
         lambda(i)   = t_s/2.0d0/dt * lambda(i)
         if(lambda(i)<=1.0d-2) lambda(i)=1.0d-2
      enddo
      do i=1,N-1
         dP_l(i) = (gammp(a,lambda(i+1)/l_c)-gammp(a,lambda(i)/l_c))
     &        /(gammp(a,1.0d2/l_c)-gammp(a,1.0d-2/l_c)) 
!         print*, M_BH(i+1), M_BH_c(1), lambda(i), dP_l(i)
!         write(10,*) M_BH(i+1), M_BH_c(1), lambda(i), dP_l(i)
      enddo


!==========================!
!     Evolve BHMF & QLF    !
!==========================!
!     time update
      do k=1,100

!     BHMF
         do i=1,N-1
            P_M_prev(i) = P_M(i) ! previous time step
         enddo
         do i=1,N-1             ! BH mass aftet dt
            P_M(i)=0.0d0
            do j=1,i            ! sum of the contribution from the previous time
               P_M(i) = P_M(i) + dP_l(i-j+1) * P_M_prev(j)
            enddo
            Phi_M(i) = P_M(i)/dlogM_BH(i)
         enddo

!     QLF
         do i=1,N-1               ! mag
            Phi_L(i) = 0.0d0
            L_sun = 3.828d33
            XX = -(mag_c(i)+21.0d0)/2.5d0 - dlog10(L_sun/1.0d45) ! log(L_bol/L_sun)
            L_bol = L_sun*(10.0d0**XX)
            f_x = 10.96d0 *(1.0d0+(XX/11.93d0)**17.79d0)
            L_x = L_bol/f_x

c            print*, mag_c(i), f_x
            f_obs(i) = dmax1(0.73d0-0.24d0*(dlog10(L_x)-43.75d0), 0.2d0)
            f_obs(i) = dmin1(0.84d0, f_obs(i))
         
            do j=1,N-1          ! BH mass
               lambda_b = 10.0d0**(-(mag_c(i)+21.0d0)/2.5d0)
     &              /1.3d-7/M_BH_c(j)
               deriv = dlog(10.0d0)/2.5d0
               call ERDF(a, lambda_b, l_c, dP_dlnl, i_sw)

c               print*, lambda_b, dP_dlnl
               Phi_L(i) = Phi_L(i) + dP_dlnl*deriv*Phi_M(j)*dlogM_BH(j)
            enddo
            
         enddo
         
!     output BHMF
c         do i=1,N-1
c            write(10,*) M_BH(i), P_M(i), Phi_M(i), mag_c(i), Phi_L(i),
c     &           f_obs(i)
c         enddo
c         write(10,*) ' '

!     check the number conservation
c         Number = 0.0d0
c         do i=1,N-1
c            Number= Number + Phi_M(i) * dlog10(M_BH(i+1)/M_BH(i))
c         enddo
c         print*, Number, '#'

!---------------------------
         time = time + dt
         if(time>1.17d3) exit
c         if(time>1.54d3) exit
      enddo

      do i=1,N-1
         write(10,*) M_BH(i), P_M(i), Phi_M(i), mag_c(i), Phi_L(i),
     &        f_obs(i)
      enddo

!     calculate chisquare
      open(12,file='N20_z5.dat')
      do k=1,N_d
         read(12,*) mag_data(k), dumy, Phi_L_data(k), dPhi_L_data(k)
      enddo
      close(12)

      i_read = 1
      chisq  = 0.0d0
      do i=1,N-1
         if(mag_c(i)<mag_data(i_read)) then
c            print*, i_read, mag_c(i), Phi_L(i), mag_data(i_read), 
c     &           Phi_L_data(i_read), dPhi_L_data(i_read)

            chisq = chisq 
     &           + (Phi_L(i)*(1.0d0-f_obs(i))-Phi_L_data(i_read))**2.0d0
     &           /dPhi_L_data(i_read)**2.0d0
            i_read = i_read + 1
         endif
         if(i_read>N_d) exit
      enddo

      print*, chisq, chisq/N_d, 'chisq'
      print*, time, 'final time'
      return
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

      SUBROUTINE ERDF(a,x,x0,P,i_sw)
      implicit none
      integer :: i, M, i_sw
      real*8 :: x, x0, a, P, dx, xmin
      real*8, save  :: Sum = 0.0d0
      
      if(i_sw==0) Sum=0.0d0

!      a=0.20d0
!      x0=0.87d0
      xmin = 1.0d-2

!      Sum = 2.5480149509670924d0
      P = (x/x0)**a * dexp(-x/x0)
      if(x<xmin) P=0.0d0

      if(Sum==0.0d0) then
         M=1000000000
         do i=1,M
            if(x>10.0d0) exit
            dx = xmin*1.0d-3
            x = xmin + dx*(i-1)
            P = (x/x0)**a * dexp(-x/x0)
            Sum = Sum + P*dx/x
         enddo
         print*, Sum, 'Sum'
         i_sw = 1
!      stop
      endif
      P = P/Sum
      return
      END



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION gammp(a,x)
      real*8 :: a,gammp,x
!     USES gcf,gser
!     Returns the incomplete gamma function P(a, x).
      real*8 ::  gammcf,gamser,gln

      if(x.lt.0.0d0 .or. a.le.0.0d0) then
         print*, 'bad argument in gammap'
      endif
      if(x.lt.a+1.0d0) then
!     Use the series representation.
         call gser(gamser,a,x,gln)
         gammp = gamser
      else 
!     Use the continued fraction representation
         call gcf(gammcf,a,x,gln)
         gammp = 1.0d0 - gammcf
      endif
      return
      END


      FUNCTION gammln(xx)
      real*8 :: gammln,xx
!     Returns the value ln[Γ(xx)] for xx > 0.
      integer :: j
      real*8 :: ser,stp,tmp,x,y,cof(6)
!     Internal arithmetic will be done in double precision, 
!     a nicety that you can omit if five-figure
!     accuracy is good enough.
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     &     24.01409824083091d0,-1.231739572450155d0,
     &     .1208650973866179d-2,
     &     -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
      enddo
      gammln=tmp+dlog(stp*ser/x)
      return
      END


      SUBROUTINE gser(gamser,a,x,gln)
      integer :: ITMAX
      real*8 :: a,gamser,gln,x,EPS
!     USES gammln
!     Returns the incomplete gamma function P(a, x) evaluated by its series representation as
!     gamser. Also returns ln Γ(a) as gln.
      integer :: n
      real*8 :: ap,del,sum,gammln
      ITMAX=100
      EPS=3.d-7
      gln=gammln(a)
      if(x.le.0.)then
         if(x.lt.0.) print*, 'x < 0 in gser'
         gamser=0.
         return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,ITMAX
         ap=ap+1.
         del=del*x/ap
         sum=sum+del
         if(abs(del).lt.abs(sum)*EPS) goto 1
      enddo
      print*, 'a too large, ITMAX too small in gser'
 1        gamser=sum*exp(-x+a*dlog(x)-gln)
      return
      END


      SUBROUTINE gcf(gammcf,a,x,gln)
      integer :: ITMAX
      real*8 :: a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!     USES gammln
!     Returns the incomplete gamma function Q(a, x) evaluated 
!     by its continued fraction representation as gammcf. 
!     Also returns ln Γ(a) as gln.
!     Parameters: ITMAX is the maximum allowed number of iterations; 
!     EPS is the relative accuracy; FPMIN is a number near the smallest 
!     representable floating-point number.
      integer :: i
      real*8 :: an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.0d0-a                  !Set up for evaluating continued fraction by modified
      c=1.0d0/FPMIN                !Lentz’s method (§5.2) with b0 = 0.
      d=1.0d0/b
      h=d
      do i=1,ITMAX           !Iterate to convergence.
         an=-i*(i-a)
         b=b+2.0d0
         d=an*d+b
         if(abs(d).lt.FPMIN)d=FPMIN
         c=b+an/c
         if(abs(c).lt.FPMIN)c=FPMIN
         d=1.0d0/d
         del=d*c
         h=h*del
         if(abs(del-1.0d0).lt.EPS)goto 1
      enddo
      print*, 'a too large, ITMAX too small in gcf'
 1        gammcf=exp(-x+a*dlog(x)-gln)*h ! Put factors in front.
      return
      END

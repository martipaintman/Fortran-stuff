      program EulerMacLaurin
      implicit none
      double precision func,pi,b,a,x,Int,h
      integer n
      pi=4.0d0*datan(1.0d0)
      a=0.0d0
      b=pi/2.0d0

      open(1,file='EulerMacLaurin.dat',status='unknown')
      do n=2,1000
         call trapecis(n,Int,a,b)
         h=(b-a)/dble(n)
         write(1,*) n,Int, Int+h**2*(2.0d0*dcos(a)-2.0d0*dcos(b))/12.0d0
      enddo
      close(1)
      stop
      end

      subroutine trapecis(n,Int,a,b)
         double precision Int,h,x,a,b,func
         integer n
         x=a
         Int=func(a)+func(b)
         h=(b-a)/dble(n)
         do n=1,n-1
            x=x+h
            Int=Int+2.0d0*func(x)
         enddo
         Int=Int*h/2.0d0
      return

      end subroutine

      double precision function func(x)
         double precision x
         func=(dsin(x))**2.0d0
         return
      end function
c No s'utilitza pero es podria utilitzar jaja
      double precision function func2(x)
         double precision x
         func2=2dcos(x)
         return
      end function

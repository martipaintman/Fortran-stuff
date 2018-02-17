c bisseccio per trobar zeros de funcions
      program bisseccio
      implicit none
      double precision xa,xb,xt,epsilon,f,x,error,y,dret
      parameter (epsilon=1.d-10)

      f(x)=dcos(x)-x-y
      open(1,file='inversa1.dat')

      y=-10.d0
      dret=10.d0

      do while (y.lt.dret)

        xa=-100.d0
        xb=100.d0
        error = 1000.d0
c        write(*,*)'q pase nen'

        do while(error.gt.epsilon)

          xt=(xa+xb)/2.d0

          if (f(xa)*f(xt).gt.0.d0) xa=xt

          if (f(xa)*f(xt).lt.0.d0) xb=xt

          error=dabs(xa-xb)

        enddo

c        write(*,*) y,xa
c        write(*,*)'q pase nen'
        write(1,*) y,xa

        y=y+0.2d0

      enddo
      close(1)

      stop
      end

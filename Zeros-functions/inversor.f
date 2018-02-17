c bisseccio per trobar zeros de funcions
      program bisseccio
      implicit none
      double precision xa,xb,xt,epsilon,f,x,error,y,dret
      parameter (epsilon=1.d-10 ,step=1.d-4)

      f(x)=dcos(x)-x-y
      open(1,file='inversa1.dat')
      
c next i define the domain of y for wich ill do the inverse
      y=-10.d0.
      dret=10.d0
      
c this do while will make sure we solve our equation for all the y (according to our step)

      do while (y.lt.dret)
      
c here i define the initial x interval where to start looking for solutions of the equation
        xa=-100.d0
        xb=100.d0
        error = 1000.d0
        
c now this do while efectuates the bisection algorithm and asures an error smaller than 'error'
        do while(error.gt.epsilon)

          xt=(xa+xb)/2.d0

          if (f(xa)*f(xt).gt.0.d0) xa=xt

          if (f(xa)*f(xt).lt.0.d0) xb=xt

          error=dabs(xa-xb)

        enddo
c when i have the final solution for the current y i write to the file both coordinates, y and x
        write(1,*) y,xa
        
c here i increment the y value to prepare the next y for the next iteration
        y=y+step

      enddo
      close(1)

      stop
      end

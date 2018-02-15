c bisseccio per trobar zeros de funcions
      program bisseccio
    

      double precision xa,xb,xt,epsilon,f,x
      integer n
      parameter (epsilon=1.d-10)

      f(x)=dcos(x)-x

      error = 1000.d0


 1    write(*,*) 'introdueix lextrem '
      read(*,*) xa
      write(*,*) 'introdueix l"extrem dret'
      read(*,*) xb

      if(xa.gt.xb) then
        write(*,*) 'interval erroni'
        goto 1
      endif

      if (f(xa)*f(xb).gt.0.d0) then
        write(*,*) 'en aquest interval no hi ha cap zero'
        goto 1
      endif



      do while(error.gt.epsilon)

        xt=(xa+xb)/2.d0

        if (f(xa)*f(xt).gt.0.d0) xa=xt

        if (f(xa)*f(xt).lt.0.d0) xb=xt
        error=dabs(xa-xb)
        write(*,*) xa,xb,error
      enddo

      stop
      end

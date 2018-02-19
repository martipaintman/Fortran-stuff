c Programa funció inversa
      program inversa
      implicit none
      double precision yb,dy,xa,xb,xt,epsilon,error,f,y,x
      integer n
c Definim tota la merda
      epsilon=1.d-10
      error=1000.d0
      n=0
      f(x)=dcos(x)-x-y
      open(1,file='inversa.dat')
c Demanem que ens digui els limits de la y
 1    write(*,*)'introdueix l''extrem esquerre de la y'
      read(*,*)y

      write(*,*)'introdueix l''extrem dret de la y'
      read(*,*)yb
      if(y.gt.yb) then
        write(*,*) 'interval erroni'
        goto 1
      endif

 2    write(*,*)'introdueix l''increment de la y'
      read(*,*)dy
      if(dy.ge.dabs(yb-y))then
        write(*,*)'Increment massa gran'
        goto 2
      endif
c Same with x lol
 3    write(*,*) 'introdueix l''extrem superior'
      read(*,*) xa
      write(*,*) 'introdueix l''extrem inferior'
      read(*,*) xb

      if(xa.gt.xb) then
        write(*,*) 'interval erroni'
        goto 3
      endif
c Obrim un while que actui desde el limit de l'esquerra fins arribar al de la
c dreta amb un increment dy per cada punt.
      do while(y.lt.yb)
c Ara un while on trobi la solucio de la funció per cada y
        do while(error.gt.epsilon)
          n=n+1
          xt=(xa+xb)/2.d0
          if (f(xa)*f(xb).gt.0.d0) then
            write(*,*) 'en aquest interval no hi ha cap zero'
            goto 3
          endif

          if (f(xa)*f(xt).gt.0.d0)then
            xa=xt
          else
            xb=xt
          endif
          error=dabs(xa-xb)
         enddo
        write(*,*) y,xt,error,n
        write(1,*) y,xt,error,n
        y=y+dy
      enddo
      stop
      end

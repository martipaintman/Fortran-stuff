c MARTI PINTO BORRELL
c PRACTICA 5
c

      program prac5
      implicit none
      integer ndat,i,ncaixes,NM,contador
      double precision xgaus,xhisto,histo,errhisto,mitja,xmu,xsigma
      double precision gauss,z,a,b,x,y,txi1,txi2,xn,yn,sig,dt,acumy
      double precision acumy2,vary,t,D
      dimension xn(250),yn(250)

c 1) Genero nums gaussians faig histo i comparo amb dist exacte
      dimension xgaus(int(12d4)),xhisto(120),histo(120),errhisto(120)
      ndat = int(12d4)
      ncaixes = 120
      xmu = 0d0
      xsigma = 1d0
      a = -5d0
      b = 5d0

      open(10,file='P5-18P-res1.dat')
      call subgaussians(ndat,xmu,xsigma,xgaus)
      do i = 1,ndat
         write(*,*) xgaus(i)
      enddo
      call histograma(ndat,xgaus,a,b,ncaixes,xhisto,histo,errhisto)
      do i=1,ncaixes
         write(10,*) xhisto(i),histo(i),errhisto(i)
      enddo
c distribucio exacte
      open(11,file='P5-18P-res2.dat')
      do i=1,ndat
         z = a + ((b-a)/ndat)*i
         write(11,*) z,gauss(z,xmu,xsigma)
      enddo
c 2) Simulacio moviment aleatori 250 molecules
c guardare cada xn i yn en 2 vectors
      open(12,file='P5-18P-res3.dat')
      open(13,file='P5-18P-res4.dat')

      dt = 0.02d0
      sig = dsqrt(2.21d-5*dt)
      contador = 1
      acumy = 0d0
      acumy2 = 0d0
      do i = 1,250
         xn(i) = 0d0
         yn(i) = 0d0
      enddo

      do i = 1,240
         t = i*dt
         acumy2 = 0d0
         acumy = 0d0
         do NM=1,250

            xn(NM) = xn(NM) + sig*xgaus(contador)
c            write(*,*) xgaus(contador),xn(NM),i
            contador = contador + 1
            yn(NM) = yn(NM) + sig*xgaus(contador)
            contador = contador + 1
            acumy = acumy + yn(NM)
            acumy2 = acumy2 + (yn(NM))**2

         enddo
         write(12,*) xn(1),yn(1),xn(2),yn(2),xn(3),yn(3),xn(4)
     +   ,yn(4),xn(5),yn(5)
         vary = acumy2/(250d0) - (acumy/250d0)**2

         write(13,*) t,vary
      enddo
      D = vary
      write(*,*) D
      stop
      end program

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      double precision function gauss(x,xmu,xsigma)
      implicit none
      double precision x,num,denom,xmu,xsigma
      num = dexp((-(x-xmu)**2)/(2d0*xsigma**2))
      denom = dsqrt(2d0*dacos(-1d0)*xsigma**2)
      gauss = num/denom
      return
      end function

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine histograma(ndat,xdata,xa,xb,ncaixes,xhisto,histo,
     +errhisto)
c Genera un histograma normalitzat de ncaixes.
c Input: ndat:Nombre de valors de la variable aleatoria x de la qual
c volem obtenir la pdf (estimada amb histograma), xdata: valors de la
c variable aleatoria, xa i xb: interval de les dades xdata del qual
c volem estimar la pdf, ncaixes: num de caixes que tindra l'histo
c output: xhisto es un vector amb els valors ordenats del punt mig de
c cada caixa, histo es un vector amb l'alçada de cada caixa, i errhisto
c es un altre vector amb l'error associat a cada caixa
      implicit none
      integer ndat,ncaixes,Nk,i,j
      double precision xdata,xa,xb,xhisto,histo,errhisto,w,zone,parenth
      double precision pk
      dimension xdata(ndat), xhisto(ncaixes),histo(ncaixes)
      dimension errhisto(ncaixes),Nk(ncaixes),pk(ncaixes)
c ndat es com N
      w = (xb-xa)/dble(ncaixes)
c      write(*,*) w
      zone = xa + w

      do i=1,ncaixes
         Nk(i)=0
         do j=1,ndat
            if ((xdata(j).lt.zone).and.(xdata(j).gt.(zone-w))) then
               Nk(i) = Nk(i) + 1
            else
               continue
            endif
         enddo
c         write(*,*) Nk(i)
         xhisto(i)=zone-w + w/2d0
c         write(*,*) xhisto(i)
         zone = zone + w
      enddo

      do i=1,ncaixes
c         write(*,*) dble(Nk(i)),w*dble(ndat)
         pk(i)=dble(Nk(i))/(w*dble(ndat))
c         write(*,*) pk(i)
         parenth=(dble(Nk(i))/dble(ndat))*(1-(dble(Nk(i))/dble(ndat)))
         errhisto(i)=(1d0/w)*dsqrt((1d0/dble(ndat))*parenth)
      enddo
      do i = 1,ncaixes
         histo(i) = pk(i)
c         write(*,*) histo(i),pk(i)
      enddo
      return
      end subroutine

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine subgaussians(ndat,xmu,xsigma,xgaus)
c Generador de nombres gaussians mitjançant metode Box-Muller
c input: ndat: quantitat de nums que volem generar,xmu: valor mitja,
c xsigma: sigma de la distribucio
c ouptu : vector xgaus amb els valors generats
      implicit none
      integer ndat,iseed,i
      double precision xmu,xsigma,xgaus,txi1,txi2,part1,part2,pi
      parameter (iseed=16850945)
      dimension xgaus(ndat)
      pi = dacos(-1d0)
      call srand(iseed)
      do i=1,ndat
         txi1 = rand()
         txi2 = rand()
         part1= dsqrt(-2d0*dlog(txi1))
         part2= dcos(2d0*pi*txi2)
         xgaus(i)= part1*part2*xsigma+xmu
      enddo
      return
      end subroutine
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

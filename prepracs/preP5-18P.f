C MARTI PINTO BORRELL
C PREPRAC 5 RANDOM NUMBERS OH YASS

      program preprac5
      implicit none
      integer ndat,i,ncaixes
      double precision xgaus,xhisto,histo,errhisto,mitja,xmu,xsigma
      double precision gauss
      dimension xgaus(int(3d4)),xhisto(150),histo(150),errhisto(150)
      ndat = int(3d4)
      ncaixes = 150
      xmu = 1d0
      xsigma = 0.45d0
      open(10,file='P5provagauss.dat')
      call subgaussians(ndat,1d0,0.45d0,xgaus)
      do i=1,ndat
         write(10,*) xgaus(i)
      enddo
c      write(*,*) maxval(xgaus),minval(xgaus)
      open(11,file='P5provahist.dat')
      call histograma(ndat,xgaus,-3d0,5d0,ncaixes,xhisto,histo,errhisto)
      do i=1,ncaixes
         write(11,*) xhisto(i),histo(i),errhisto(i)
      enddo
c estimacio valor mitja
      mitja = 0d0
      do i=1,ndat
c         write(*,*) xgaus(i),gauss(xgaus(i),xmu,xsigma)
         mitja = mitja + xgaus(i)*gauss(xgaus(i),xmu,xsigma)
      enddo
      write(*,*) mitja

      stop
      end program
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      double precision function gauss(x,xmu,xsigma)
      implicit none
      double precision x,num,denom,xmu,xsigma
      num = dexp((-(x-xmu)**2)/(2d0*xsigma**2))
      denom = dsqrt(2d0*dacos(-1d0)*xsigma**2)
      gauss = num/denom
      return
      end function

      double precision function fun(x)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c 1) Histogram

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
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c 2) Distribucio Gaussiana
      subroutine subgaussians(ndat,xmu,xsigma,xgaus)
c Generador de nombres gaussians mitjançant metode Box-Muller
c input: ndat: quantitat de nums que volem generar,xmu: valor mitja,
c xsigma: sigma de la distribucio
c ouptu : vector xgaus amb els valors generats
      implicit none
      integer ndat,iseed,i
      double precision xmu,xsigma,xgaus,txi1,txi2,part1,part2,pi
      parameter (iseed=16874605)
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

c 3) Acceptacio i rebuig

      subroutine subaccepta(ndat,a,b,M,fun,xnums)
c Genera nombres aleatoris segons la distribucio fun
c input: limits a i b , cota superior M, distribucio fun(external) i
c quantitat de punts desitjats ndat
c output: vector xnums amb el resultat
      implicit none
      external fun
      integer ndat,iseed
      double precision a,b,M,fun,xnums,txi1,txi2,x,p,mitja
      dimension xnums(ndat)
      sac = 0
      iseed = 16874605
      mitja = 0d0
      call srand(iseed)
      do while ((sac).lt.ndat)
         txi1= rand()
         txi2= rand()
         x = (b-a)*txi1+a
         p = M*txi2
         if (fun(x).ge.p) then
            sac = sac + 1
            xnums(sac) = x
            mitja = mitja + x
         else
            continue
         endif
      enddo

      return
      end subroutine

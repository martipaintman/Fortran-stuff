c PRE 6
c MARTI PINTO BORRELL

      program pre6
      implicit none
      external f1,f2,f3,f4,f5,prob1,prob2,pdfc,gauss,g,pmulti
      double precision f1,f2,pi,a,b,txi,xdist,xpunts,prob1,err,sol,f3,g
      double precision emerr,emsol,prob2,xgaus,pdfc,xaccept,f4,f5,gauss
      double precision pmulti
      integer N,iseed,i,nmax
      parameter (iseed=16874605)
      common /dades/ xgaus
      dimension xpunts(1024000),emerr(10),emsol(10),xgaus(1024000)
      dimension xaccept(1024000)
      nmax = 1024000
      pi = dacos(-1d0)

c 1a) Integrar amb montecarlo cru, necesito generar nums aleat uniforms
c 1 I1
      a = -dsqrt(5d0)
      b =  -a
      call uniform(nmax,a,b,xpunts)

      N = int(2d3)
      i = 1
      do while (N.le.nmax)
         call montecarloP6(N,xpunts,f1,prob1,err,sol)
         write(*,*) 5d0*pi/2d0,sol,err
         emerr(i) = err
         emsol(i) = sol
         i = i + 1
         N = N*2
      enddo
      write(*,*) '--------------------------------------------------'
c 1 I2
      a = -2d0*pi
      b =  -a
      call uniform(nmax,a,b,xpunts)

      N = int(2d3)
      i = 1
      do while (N.le.nmax)
         call montecarloP6(N,xpunts,f2,prob2,err,sol)
         write(*,*) -(8d0*pi)/15d0,sol,err
         emerr(i) = err
         emsol(i) = sol
         i = i + 1
         N = N*2
      enddo
      write(*,*) '--------------------------------------------------'
c 1b) generar 1d6 nums gaussians valor mitja 0 i variancia 1
c tocar aqet num i els de dalt per fer mes nums
      nmax = int(1d6)
      call subgaussians(nmax,0d0,1d0,xgaus)

c 1c) generar 1d6 segons distr p(x)... usant subaccepta

      open(10,file='preP6-18Pvis.dat')
      do i = 1,int(1d2)
         a = -2d0*pi + (i-1d0)*(4d0*pi)/1d2
         write(10,*) a,pdfc(a)
      enddo
      a = -2d0*pi
      b =  -a
      call subaccepta(nmax,a,b,0.4d0,pdfc,xaccept)

      open(11,file='preP6-18Paccept.dat')
      do i = 1,nmax
         write(11,*) xaccept(i)
      enddo

c 1d) calcular integrals 3 4 i 5 amb els nous nums aleat generats
c I3
c      a = -(2d0**20)
c      b =  -a

      N = 2500
      i = 1
      do while (N.le.nmax)
         call montecarloP6(N,xgaus,f3,gauss,err,sol)
         write(*,*) N,sol,err
         emerr(i) = err
         emsol(i) = sol
         i = i + 1
         N = N*2
      enddo
      write(*,*) '--------------------------------------------------'

c I4
c      a = -(2d0**20)
c      b =  -a

      N = 2500
      i = 1
      do while (N.le.nmax)
         call montecarloP6(N,xgaus,f4,gauss,err,sol)
         write(*,*) N,sol,err
         emerr(i) = err
         emsol(i) = sol
         i = i + 1
         N = N*2
      enddo
      write(*,*) '--------------------------------------------------'

c I5

c      a = -2d0*pi
c      b =  -a
      N = 2500
      i = 1
      do while (N.le.nmax)
         call montecarloP6(N,xaccept,f5,pdfc,err,sol)
         write(*,*) N,sol,err
         emerr(i) = err
         emsol(i) = sol
         i = i + 1
         N = N*2
      enddo
      write(*,*) '--------------------------------------------------'

c 2) multidimensional
      nmax = 200000
      N = 2500
      i = 1
      do while (N.le.nmax)
         call montemultiP6(N,g,pmulti,err,sol)
         write(*,*) N,sol,err
         emerr(i) = err
         emsol(i) = sol
         i = i + 1
         N = N*2
      enddo
      write(*,*) '--------------------------------------------------'
      stop
      end program
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      double precision function f1(x)
      implicit none
      double precision x
      f1 = dsqrt(5d0-x**2)
      return
      end function

      double precision function f2(x)
      implicit none
      double precision x
      f2 = x*((dcos(x))**2)*(dsin(x)**3)
      return
      end function

      double precision function f3(x)
      implicit none
      double precision x
      f3 = exp((-x**2)/2d0)*(1d0/dsqrt(8*dacos(-1d0)))*x*dsin(x)
      return
      end function

      double precision function f4(x)
      implicit none
      double precision x
      f4 = exp((-x**2)/2d0)*dcos(x)
      return
      end function

      double precision function f5(x)
      implicit none
      double precision x
      f5 = (dsin(x)**6)*x**2
      return
      end function

      double precision function g(x1,x2,x3,x4,x5)
      implicit none
      double precision x1,x2,x3,x4,x5,part1,part2,part3
      part1 = (x3**2)*dcos(x4)*dsin(x5**2)
      part2 = (x3**2)*(1d0+x1)*(1d0-x4)*x2**2
      part3 = (dcos(x4)**2)*x5**2
      g = (part1+part2+part3)*exp(-(x1**2+x2**2+x3**2+x4**2+x5**2))
      return
      end function

      double precision function prob1(x)
      implicit none
      double precision x,a,b
      a = -dsqrt(5d0)
      b =  -a
      prob1 = 1d0/(b-a)
      return
      end function

      double precision function prob2(x)
      implicit none
      double precision x,a,b
      a = -2d0*dacos(-1d0)
      b =  -a
      prob2 = 1d0/(b-a)
      return
      end function

      double precision function gauss(x)
      implicit none
      double precision x,num,denom,xmu,xsigma
      xmu = 0d0
      xsigma = 1d0
      num = dexp((-(x-xmu)**2)/(2d0*xsigma**2))
      denom = dsqrt(2d0*dacos(-1d0)*xsigma**2)
      gauss = num/denom
      return
      end function

      double precision function pdfc(x)
      implicit none
      double precision x,pi
      pi = dacos(-1d0)
      pdfc= (1d0/(2*pi**3-(15d0/16d0)*pi))*(dsin(x)**4)*x**2
      return
      end function

      double precision function pmulti(x1,x2,x3,x4,x5)
      implicit none
      double precision x1,x2,x3,x4,x5,gauss
      pmulti = gauss(x1)*gauss(x2)*gauss(x3)*gauss(x4)*gauss(x5)
      return
      end function

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine uniform(ndat,a,b,xout)
c guarda ndat punts segons distribucio uniforme entre a i b al
c vector xout
      implicit none
      double precision xout,a,b,txi
      integer ndat,i
      dimension xout(ndat)

      call srand(16874605)
      do i = 1,ndat
         txi = rand()
         xout(i) = a + (b-a)*txi
      enddo
      return
      end subroutine
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine montecarloP6(ndat,xpunts,func,prob,sigf,festim)
c Subrutina d'integracio metode montecarlo cru. Input: ndat:Nombre de
c punts aleatoris dels que disposem, xpunts: vector amb els punts
c aleatoris amb els que integrarem, func: funcio que integrarem
c prob : pdf amb la que estem integrant
c output: sigf : valor amb l'error asociat a la integral.
c festim: valor estimat de la integral
      implicit none
      external func,prob
      double precision xpunts,sigf,festim,facum,sigacum1,sigacum2,funca
      double precision func,prob,var
      integer ndat,i
      dimension xpunts(1024000)
c      write(*,*) ndat
c Poso a 0 el valor inicial de cada acumulador
      facum = 0d0
      sigacum1 = 0d0
      sigacum2 = 0d0

      do i = 1,ndat
         funca = func(xpunts(i))/prob(xpunts(i))
         facum = facum + funca
         sigacum1 = sigacum1 + funca**2
      enddo
      sigacum2 = facum
      festim = facum/dble(ndat)
      var = sigacum1/dble(ndat)-(sigacum2/dble(ndat))**2
      sigf = dsqrt(var/dble(ndat))
      return
      end subroutine
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine montemultiP6(ndat,func,prob,sigf,festim)
c Subrutina d'integracio metode montecarlo cru. Input: ndat:Nombre de
c punts aleatoris dels que disposem, xpunts: vector amb els punts
c aleatoris amb els que integrarem, func: funcio que integrarem(multi)
c prob : pdf amb la que estem integrant
c output: sigf : valor amb l'error asociat a la integral.
c festim: valor estimat de la integral
      implicit none
      external func,prob
      double precision sigf,festim,facum,sigacum1,sigacum2,funca
      double precision func,prob,var,xgaus,xpunts1,xpunts2,xpunts3
      double precision xpunts4,xpunts5
      integer ndat,i,j,bonus
      common /dades/ xgaus
      dimension xpunts1(200000),xpunts2(200000),xpunts3(200000)
      dimension xpunts4(200000),xpunts5(200000),xgaus(1024000)
c      write(*,*) ndat
      bonus = 200000
      do i = 1,200000
         xpunts1(i)=xgaus(i)
         xpunts2(i)=xgaus(i+bonus)
         xpunts3(i)=xgaus(i+bonus*2)
         xpunts4(i)=xgaus(i+bonus*3)
         xpunts5(i)=xgaus(i+bonus*4)
      enddo
c Poso a 0 el valor inicial de cada acumulador
      facum = 0d0
      sigacum1 = 0d0
      sigacum2 = 0d0

      do i = 1,ndat
         funca = func(xpunts1(i),xpunts2(i),xpunts3(i),xpunts4(i),
     + xpunts5(i))/prob(xpunts1(i),xpunts2(i),xpunts3(i),xpunts4(i),
     + xpunts5(i))
         facum = facum + funca
         sigacum1 = sigacum1 + funca**2
      enddo
      sigacum2 = facum
      festim = facum/dble(ndat)
      var = sigacum1/dble(ndat)-(sigacum2/dble(ndat))**2
      sigf = dsqrt(var/dble(ndat))
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
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine subaccepta(ndat,a,b,M,fun,xnums)
c Genera nombres aleatoris segons la distribucio fun
c input: limits a i b , cota superior M, distribucio fun(external) i
c quantitat de punts desitjats ndat
c output: vector xnums amb el resultat
      implicit none
      external fun
      integer ndat,iseed,sac
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
c            mitja = mitja + x
         else
            continue
         endif
      enddo

      return
      end subroutine

c MARTI PINTO BORRELL
c PRACTICA 6

      program prac6
      implicit none
      external f1,prob1,pdfc,f2,f3,pmulti
      double precision f1,prob1,pi,a,b,xpunts,err,sol,emerr,emsol,L,f2
      double precision xaccept,pdfc,f3,pmulti
      integer nmax,i,N
      common /dades/ xaccept
      dimension xpunts(1000000),emerr(200),emsol(200),xaccept(1000000)

      nmax = int(1d6)
      pi = dacos(-1d0)
c 1a)
      a = 0d0
      b =  3d0*pi/2d0
      call uniform(nmax,a,b,xpunts)

      open(10,file='P6-18Pres1.dat')
c la N es el numero inicial de punts que usare en la primera aplicacio
c del metode , despres la vaig augmentant per anar aplicant el metode
c amb mes punts
      N = int(1d4)
c la i la faig servir per iterar pels vectors on emmagatzemo la solucio
c i l'error (no els utilitzo pero)
      i = 1
      do while (N.le.nmax)
         call montecarloP6(N,xpunts,f1,prob1,err,sol)
c         write(*,*) sol,err
c emerr i emsol son vectors on emmagatzemo l'error i la solucio
c pero no els faig servir
         emerr(i) = err
         emsol(i) = sol
         write(10,*) N,sol,err

         i = i + 1
c aqui augmento la N per a la seguent iteracio
         N = N + 10000
      enddo
      write(*,*) '--------------------------------------------------'
c1b)
      L = 16d0
      a = -L/2d0
      b =  -a
      call subaccepta(nmax,a+1d-6,b-1d-6,2d0/L,pdfc,xaccept)
c1c)

      open(11,file='P6-18Pres2.dat')
      N = int(1d4)
      i = 1
      do while (N.le.nmax)
         call montecarloP6(N,xaccept,f2,pdfc,err,sol)
c         write(*,*) sol,err

         emerr(i) = err
         emsol(i) = sol
         write(11,*) N,sol,err

         i = i + 1
         N = N + 10000
      enddo
      write(*,*) '--------------------------------------------------'

c 2 fermions

      open(12,file='P6-18Pres3.dat')
      nmax = 250000
      N = int(1d4)
      i = 1
      do while (N.le.nmax)
         call montemultiP6(N,f3,pmulti,err,sol)
c         write(*,*) sol,err

         emerr(i) = err
         emsol(i) = sol
         write(12,*) N,sol,err

         i = i + 1
         N = N + 10000
      enddo
      write(*,*) '--------------------------------------------------'

      stop
      end program


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      double precision function f1(x)
      implicit none
      double precision x
      f1 = (x**3)*dsin(x)**3
      return
      end function

      double precision function f2(x)
      implicit none
      double precision x,pi,pdfc,L
      pi = dacos(-1d0)
      L = 16d0
      f2 = (dsin((8*pi*(x-L/2d0))/L)**2)*pdfc(x)
      return
      end function

      double precision function arg1(x)
      implicit none
      double precision x,pi,L
      pi = dacos(-1d0)
      L = 16d0
      arg1 = dsin((pi*(x-L/2d0))/(L))

      return
      end function

      double precision function arg2(x)
      implicit none
      double precision x,pi,L
      pi = dacos(-1d0)
      L = 16d0
      arg2 = dcos((pi*(x-L/2d0))/(L))
      return
      end function

      double precision function f3(x1,x2,x3)
      implicit none
      double precision x1,x2,x3,prod1,prod2,arg1,arg2

      prod1= arg1(x1)*arg1(x2)*arg1(x3)
      prod2=(arg2(x1)-arg2(x2))*(arg2(x2)-arg2(x3))*(arg2(x1)-arg2(x3))
      f3 = abs(prod1*prod2)**2

      return
      end function

      double precision function prob1(x)
      implicit none
      double precision x,a,b
      a = 0d0
      b =  3d0*dacos(-1d0)/2d0
      prob1 = 1d0/(b-a)
      return
      end function


      double precision function pdfc(x)
      implicit none
      double precision x,pi,L
      L = 16d0
      pi = dacos(-1d0)
      pdfc= (2d0/L)*dsin((pi*(x-(L/2d0)))/L)**2
      return
      end function


      double precision function pmulti(x1,x2,x3)
      implicit none
      double precision x1,x2,x3,pdfc

      pmulti=pdfc(x1)*pdfc(x2)*pdfc(x3)
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
      dimension xpunts(1000000)
c      write(*,*) ndat
c Poso a 0 el valor inicial de cada acumulador
c facum es l acumulador de f,dividit per ndat sera lestimador
c sigacum1 i 2 son acumuladors de la sigma (amb ells primer trobare la
c variancia i despres a partir d 'aquesta la sigma)
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

c el seguent codi es montecarlo multidimensional(3 D) usant common
c amb el qual ens estalviem haver d'entrar la llista de nums aleat
c a utilitzar per integrar com un parametre, pero apart d'axio es molt
c similar a ferho sense common
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
      double precision func,prob,var,xaccept,xpunts1,xpunts2,xpunts3
      integer ndat,i,j,bonus
      common /dades/ xaccept
      dimension xpunts1(250000),xpunts2(250000),xpunts3(250000)
      dimension xaccept(1000000)
c la subrutina es molt similar a montecarloP6 amb la diferencia que ara
c func es una funcio de 3 dimensions, igual que prob, per tant el
c primer que faig es separar la llista de 1 milio de punts aleat en
c 3 llistes de un quart de milio cada una 
c      write(*,*) ndat
      bonus = 250000
      do i = 1,250000
         xpunts1(i)=xaccept(i)
         xpunts2(i)=xaccept(i+bonus)
         xpunts3(i)=xaccept(i+bonus*2)

      enddo
c Poso a 0 el valor inicial de cada acumulador
      facum = 0d0
      sigacum1 = 0d0
      sigacum2 = 0d0

      do i = 1,ndat
         funca = func(xpunts1(i),xpunts2(i),xpunts3(i))/prob(xpunts1(i),
     + xpunts2(i),xpunts3(i))
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

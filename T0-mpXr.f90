!ESTRELA ESTRANHA
!Caso relativ�stico (quarks com massa)
!Programa para calcular a massa e a press�o como fun��es do raio, do centro da estrela at� a a superf�cie, onde p = 0. Neste ponto, o programa para, imprimindo o raio total e a massa total da estrela na tela.
MODULE globals
  !pi, velocidade da luz no v�cuo, constante reduzida de Planck
  REAL, PARAMETER :: pi = 4.*atan(1.0), c = 2.99792458d23, hbar = 6.58211915d-25 ![c] = fm/s, [hbar] = GeV.s 
  !R0 = G*Msol/c^2 // G = 6.67408d-11, Msol = 1.989d30 ![G] = m^3kg^-1s^-2, [Msol] = kg
  REAL, PARAMETER :: R0 = 1.47611 !km
  !M0 = 10^54/(6,242d9*Msol*c**2) !incluindo a convers�o de m^3 para fm^3
  REAL, PARAMETER :: M0 = 8.9616d-4 !GeV 
  !press�o de sacola
  REAL, PARAMETER :: B = (0.155**4)*(5.07**3) !GeV/fm^3
  !densidade bari�nica nuclear
  REAL, PARAMETER :: n0 = 1.2277d-3 !GeV^3 // n0 = 0.16 fm^-3
  !flavour
  INTEGER :: f    
    
END MODULE globals

PROGRAM quarkstar_massivequarks
  USE globals
  IMPLICIT NONE  
  !raio, passo do raio
  REAL :: r, dr = 1d-3 !km
  !massa // mi evita a divis�o por zero na equa��o TOV
  REAL :: m, mi = 1d-3 !Msol
  !press�o, equa��o de estado, densidade de energia
  REAL :: p, pk, e !GeV/fm^3
  !momentum de Fermi 
  REAL :: k !GeV
  !densidade bari�nica, densidade bari�nica central
  REAL :: n, nc !GeV^3
  !massa dos quarks, potencial qu�mico
  REAL :: mf(3), mu(3)
10 format(E14.7, 1x, E14.7, 1x, E14.7, 1x, E14.7) !formato de impress�o
11 format(1x, A, 14x, A, 14x, A, 14x, A)

  !massa dos quarks
  mf(1) = 0.005
  mf(2) = 0.007
  mf(3) = 0.150
 
  !condi��es iniciais
  nc = 5*n0
  k = (pi**2*nc)**(1/3.)
  p = pk(k, mf, 0)
  r = 0d0 
  m = 0d0
  n = nc
  OPEN(10,FILE = 'T0-pmXrn.dat')
  WRITE(10,11) 'r', 'p', 'm', 'n/n0'
  DO WHILE (p > 0d0)
     WRITE(10,10) r, p*1d3, m, n/n0 ![r] = km, [p] = MeV/fm^3, [m] = Msol
     CALL newton(p, k, mf) !m�todo para encontrar ra�zes
     r = r + dr !vari�vel de integra��o
     p = -R0*e(k,mf)*m*dr/r**2.*(1.+p/e(k,mf))*(1.+4.*pi*r**3.*p*M0/mi)/(1.-2.*R0*m/r) + p !integra��o da TOV
     m = 4.*pi*r**2.*e(k,mf)*dr*M0 + m !integra��o da massa
     mi = m
     n = k**3/pi**2 !densidade bari�nica
  END DO
  CLOSE(10)
  PRINT*, r, m
  
END PROGRAM quarkstar_massivequarks

SUBROUTINE newton(p, k0, mf)
  USE globals
  IMPLICIT NONE
  REAL :: pk, dpdk, p
  REAL :: k0, k
  REAL :: mf(3), mu(3)
    
  DO
     k = k0 - pk(k0,mf,p)/dpdk(k0,mf)
     IF (abs(pk(k0,mf,p)) <= 1d-11) EXIT
     k0 = k
 END DO
 
END SUBROUTINE newton

REAL FUNCTION pk(k,mf,p)
  USE globals
  IMPLICIT NONE
  REAL :: k, p
  REAL :: mu(3), mf(3)
  
  pk = 0.
  DO f = 1, 3
     mu(f) = sqrt(mf(f)**2 + k**2)
     pk = 1/(4.*pi**2)*(mu(f)*k*(mu(f)**2-5/2.*mf(f)**2)+3/2.*mf(f)**4*log((mu(f)+k)/mf(f))) + pk
  END DO
  pk = pk/(hbar*c)**3 ![pk]: GeV^4 -> GeV/fm^3
  pk = pk - B - p
  
END FUNCTION pk

REAL FUNCTION dpdk(k, mf)
  USE globals
  IMPLICIT NONE
  REAL :: k, mu(3), mf(3)
  
  dpdk = 0.
  DO f = 1, 3
     mu(f) = sqrt(mf(f)**2 + k**2)
     dpdk = 1/(pi**2)*(mu(f)**3 - 2*mf(f)**2*mu(f) + mf(f)**4/mu(f)) + dpdk
  END DO
  dpdk = dpdk/(hbar*c)**3 ![dpdk] = GeV^3 -> fm^3 
  
END FUNCTION dpdk

REAL FUNCTION e(k, mf)
  USE globals
  IMPLICIT NONE
  REAL :: k, mu(3), mf(3)
  
  e = 0.
  DO f = 1, 3
     mu(f) = sqrt(mf(f)**2 + k**2)
     e = 3/(4.*pi**2)*(k*mu(f)*(k**2 - 1/2.*mf(f)**2) - 1/2*mf(f)**4*log((k+mu(f))/mf(f))) + e
  END DO
  e = e/(hbar*c)**3 ![e]: GeV^4 -> GeV/fm^3
  e = e + B
  
END FUNCTION e

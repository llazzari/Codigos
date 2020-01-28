!ESTRELA ESTRANHA
!Limite ultrarrelativ�stico (quarks sem massa)
!Programa para calcular a massa e a press�o como fun��es do raio, do centro da estrela at� a a superf�cie, onde p = 0. Neste ponto, o programa para, imprimindo o raio total e a massa total da estrela na tela.
PROGRAM quarkstar
  IMPLICIT NONE
  !pi, velocidade da luz no v�cuo, constante reduzida de Planck
  REAL, PARAMETER :: pi = 4.*atan(1.0), c = 2.99792458d23, hbar = 6.58211915d-25 ![c] = fm/s, [hbar] = GeV.s 
  !R0 = G*Msol/c^2 // G = 6.67408d-11, Msol = 1.989d30 ![G] = m^3kg^-1s^-2, [Msol] = kg
  REAL, PARAMETER :: R0 = 1.47611 !km
  !M0 = 10^54/(6,242d9*Msol*c**2) !incluindo a convers�o de m^3 para fm^3
  REAL, PARAMETER :: M0 = 8.9616d-4 !GeV 
  !press�o de sacola
  REAL, PARAMETER :: B = 0.185**4 !GeV^4
  !densidade bari�nica nuclear
  REAL, PARAMETER :: n0 = 1.2277d-3 !GeV^3 // n0 = 0.16 fm^-3
  !raio, passo do raio
  REAL :: r, dr = 1d-3 !km
  !massa // mi evita a divis�o por zero na equa��o TOV
  REAL :: m, mi = 1d-3 !Msol
  !press�o, densidade de energia
  REAL :: p, e  !GeV^4
  !densidade bari�nica, densidade bari�nica central
  REAL :: n, nc !GeV^3
10 format(E14.7, 1x, E14.7, 1x, E14.7, 1x, E14.7) !formato de impress�o
11 format(1x, A, 14x, A, 14x, A, 14x, A)

  !condi��es iniciais
  nc = 8*n0
  n = nc
  p = 3*pi**(2/3.)*n**(4/3.)/4. - B 
  e = 3*p + 4*B 
  r = 0d0 
  m = 0d0
  
  OPEN(10, FILE = 'MQ-mpXrn.dat') !arquivo gerado pelo programa
  WRITE(10,11) 'r', 'p', 'm', 'e'
  WRITE(10,10) r, p*1d3/(hbar*c)**3, m, e*1d3/(hbar*c)**3 ![r] = km, [p] = MeV/fm^3, [m] = Msol
  DO WHILE (p > 0d0)
     r = r + dr !vari�vel de integra��o
     p = -R0*e*m*dr/(r**2.)*(1.+p/e)*(1.+4.*pi*r**3.*p*M0/mi)/(1.-2.*R0*m/r) + p !integra��o da TOV
     m = 4.*pi*r**2.*e*dr*M0/((hbar*c)**3) + m !integra��o da massa
     mi = m
     e = 3*p + 4*B !equa��o de estado
     !n = (4*(p+B)/(3*pi**(2/3.)))**(3/4.) !densidade bari�nica como fun��o da press�o
     WRITE(10,10) r, p*1d3/(hbar*c)**3, m, e*1d3/(hbar*c)**3 ![r] = km, [p] = MeV/fm^3, [m] = Msun
  END DO
  CLOSE(10)
  !raio e massa totais da estrela estranha, impressos na tela
  PRINT*, r, m !km, Msol
  
END PROGRAM quarkstar

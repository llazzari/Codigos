!Módulos nrtype, nrutil e subrotinas splint retiradas do Numerical Recipes for Fortran 90, sem intenção de cópia.
!Outros módulos e subrotinas feitas pelo autor.
!Solução da equação TOV a partir de uma equação de estado tabelada, usando a interpolação spline e o método de Runge-Kutta Quarta-Ordem com passo constante para integração.
!Programa utilizado para resolver a equação TOV considerando uma equação de estado com elétrons.
MODULE nrtype
  !Symbolic names for kind types of 4-, 2-, and 1-byte integers:
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(9)
  !Symbolic names for kind types of single- and double-precision reals:
  INTEGER, PARAMETER :: SP = KIND(1.0d0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  !Symbolic names for kind types of single- and double-precision complex:
  INTEGER, PARAMETER :: SPC = KIND((1.0d0,1.0d0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  !Symbolic name for kind type of default logical:
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  !Frequently used mathematical constants (with precision to spare):
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  !Derived data types for sparse matrices, single and double precision (see use in Chapter B2):
  TYPE sprs2_sp
     INTEGER(I4B) :: n,len
     REAL(SP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  TYPE sprs2_dp
     INTEGER(I4B) :: n,len
     REAL(DP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp
END MODULE nrtype

MODULE nr
  INTERFACE
     FUNCTION locate(xx,x)
       USE nrtype
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
       REAL(SP), INTENT(IN) :: x
       INTEGER(I4B) :: locate
     END FUNCTION locate
  END INTERFACE
  INTERFACE
     FUNCTION splint(xa,ya,y2a,x)
       USE nrtype
       REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: splint
     END FUNCTION splint
  END INTERFACE
END MODULE nr

MODULE nrutil
  USE nrtype
  !Parameters for crossover from serial to parallel algorithms (these are used only within this nrutil module):
  IMPLICIT NONE
  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8 !Each NPAR2 must be ≤ the
  INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2 !corresponding NPAR.
  INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
  INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
  INTEGER(I4B), PARAMETER :: NPAR_POLY=8
  INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE assert_eq
CONTAINS
  SUBROUTINE nrerror(string)
    !Report a message, then die.
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
  FUNCTION assert_eq2(n1,n2,string)
    !Report and die if integers not all equal (used for size checking).
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2
  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3.and. n3== n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn
END MODULE nrutil

MODULE global_variables
  !pi, velocidade da luz no vácuo, constante reduzida de Planck
  REAL, PARAMETER :: pi = 4.*atan(1.0d0)
  !REAL, PARAMETER :: c = 2.99792458d23, hbar = 6.58211915d-25 ![c] = fm/s, [hbar] = GeV.s 
  !R0 = G*Msol/c^2 // G = 6.67408d-11, Msol = 1.989d30 ![G] = m^3kg^-1s^-2, [Msol] = kg
  REAL, PARAMETER :: R0 = 1.47611 !km
  !M0 = Msol/10^54 !Msol em GeV e 10^54 converte de km^3 para fm^3
  REAL, PARAMETER :: M0 = 1.115872166d3
  !densidade bariônica nuclear
  REAL, PARAMETER :: n0 = 1.2277d-3 !GeV^3 // n0 = 0.16 fm^-3
  !masses
  REAL, ALLOCATABLE, DIMENSION(:) :: e, p
  REAL :: p0, pp, e0
  REAL :: m, mi0
  INTEGER :: l

END MODULE global_variables

SUBROUTINE spline(x,y,d2y,n)
  IMPLICIT NONE
  INTEGER :: K, J, n
  REAL :: dr, r
  REAL, DIMENSION(0:n) :: d2y, y, x
  REAL, DIMENSION(0:n) :: a, b, c, d 
  REAL, DIMENSION(n) :: s
  REAL, DIMENSION(n-1,n) :: M
  
  d2y = 0.
  DO K = N-1, 1, -1
     DO J = N-1, 1, -1
        M(K,J) = (x(J+1) - x(J))
        IF (K == J) M(K,J) = 2.*(x(J+1)-x(J-1))    
     END DO
     M(K,N) = 6.*(((y(K+1) - y(K))/(x(K+1) - x(K))) - ((y(K) - y(K-1))/(x(K)-x(K-1))))
  END DO
  M(1,N-1) = 0.
  M(N-1,1) = 0.
    
  CALL fatLU(M,d2y,N-1)

  dr = 0.005
  r = x(0)
  
  DO K = 1, N
     a(K) = (d2y(K) - d2y(K-1))/(6*(x(K) - x(K-1)))
     b(K) = d2y(K)/2.
     c(K) = (x(K) - x(K-1))/6. *(2*d2y(K) + d2y(K-1)) + (y(K) - y(K-1))/(x(K) - x(K-1))
     d(K) = y(K)
     DO WHILE (r <= x(K))
        S(K) = a(K)*(r - x(K))**3. + b(K)*(r - x(K))**2. + c(K)*(r - x(K)) + d(K)
        r = r + dr
     END DO
  END DO
  
END SUBROUTINE spline

SUBROUTINE fatLU(A,x,col)
  IMPLICIT NONE
  INTEGER :: col
  INTEGER :: I, J, K
  REAL :: A(col,col+1)
  REAL, ALLOCATABLE :: LU(:,:)
  REAL :: x(0:col+1), y(0:col+1)
  REAL ::  sum, m

  ALLOCATE (LU(col,col))
  y = 0.
  x = 0.
  DO J = 1, col
     LU(1,J) = A(1,J)
  END DO

  DO I = 1, col-1
     DO J = I + 1, col
        m = a(J,I)/a(I,I)
        LU(J,I) = m
        DO K = 1, col
           a(J,K) = a(J,K) - m*a(I,K)
           IF (J <= K) LU(J,K) = a(J,K)
        END DO
     END DO
  END DO

  y = 0.0
  
  DO J = 1, col
     sum = 0.0
     DO K = 1, J-1
        sum = sum + y(K)*LU(J,K)
     END DO
     y(J) = a(J,col+1) - sum
  END DO
  
  DO J = col, 1, -1
     sum = 0.0
     DO K = J+1, col
        sum = sum + x(K)*LU(J,K)
     END DO
     x(J) = (y(J) - sum)/(LU(J,J))
  END DO

END SUBROUTINE fatLU

FUNCTION splint(xa,ya,y2a,x)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  USE nr, ONLY: locate
  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: xa,ya,y2a
  REAL, INTENT(IN) :: x
  REAL :: splint
  INTEGER :: khi,klo,n
  REAL :: a,b,h
  n=assert_eq(size(xa),size(ya),size(y2a),'splint')
  klo=max(min(locate(xa,x),n-1),1)
  khi=klo+1
  h=xa(khi)-xa(klo)
  if (h == 0.0) call nrerror('bad xa input in splint')
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
END FUNCTION splint

FUNCTION locate(xx,x)
  USE nrtype
  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: xx
  REAL, INTENT(IN) :: x
  INTEGER :: locate
  INTEGER :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
     if (ju-jl <= 1) exit
     jm=(ju+jl)/2
     if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm
     else
        ju=jm
     end if
  end do
  if (x == xx(1)) then
     locate=1
  else if (x == xx(n)) then
     locate=n-1
  else
     locate=jl
  end if
END FUNCTION locate

PROGRAM TOV
  USE global_variables
  USE nr, ONLY: splint
  IMPLICIT NONE
10 format(E16.9, 1x, E16.9) !formato de impressão
11 format(3(A,14x))
12 format(E14.7, 1x, E14.7, 1x, E14.7) !formato de impressão
  INTEGER :: j, stat=1
  REAL, ALLOCATABLE :: d2y(:)

  OPEN (10, FILE='EoS-MIT.dat')
  l = 0
  DO 
     READ(10, *, iostat=stat)
     IF (stat < 0) EXIT
     l = l + 1
  END DO
  l = l - 1 
  
  ALLOCATE(e(0:l),p(0:l),d2y(0:l))

  REWIND(10)

  DO j = 0, l
     READ(10,10) e(j), p(j)
  END DO
  CLOSE(10)

  CALL spline(p, e, d2y, l)

  OPEN(12,FILE='MxecR-MIT.dat')
  WRITE(12,11) 'p0','m','r'
  p0 = p(l)
  e0 = e(l)
  DO WHILE (p0 > 0d0)

     CALL rk4o(d2y)

     p0 = p0 - 1d-2
     e0 = splint(p, e, d2y, p0)
  END DO
  CLOSE(12)
END PROGRAM TOV

SUBROUTINE rk4o(d2y)
  USE global_variables
  IMPLICIT NONE
  REAL, DIMENSION(2) :: y0, y, k1, k2, k3, k4, f
  REAL, PARAMETER :: dr = 1d-3
  REAL :: d2y(0:l)
  REAL :: r
12 format(E14.7, 1x, E14.7, 1x, E14.7) !formato de impressão

  r = 1d-10
  y0(1) = p0
  y0(2) = 0d0
  
  DO WHILE (y0(1) > 0d0)
     CALL funcs(r,y0,f,d2y)
     k1 = f
     CALL funcs(r+dr/2.,y0+dr*k1/2.,f,d2y)
     k2 = f
     CALL funcs(r+dr/2.,y0+dr*k2/2.,f,d2y)
     k3 = f
     CALL funcs(r+dr,y0+dr*k3,f,d2y)
     k4 = f
     y = y0 + dr/6.*(k1 + 2*k2 + 2*k3 + k4)
     y0 = y
     r = r + dr
  END DO

  WRITE(12,12) p0*1d3,y(2),R
  
END SUBROUTINE rk4o

SUBROUTINE funcs(r,y,f,d2y)
  USE global_variables
  USE nr, ONLY: splint
  IMPLICIT NONE
  REAL :: r, y(2), dpdr, dmdr
  REAL :: f(2), ee, d2y(0:l)

  pp = y(1)
  m = y(2)
  ee = splint(p,e,d2y,pp)
  
  dpdr = -(ee + pp)/r*(m+4*pi*r**3*pp/M0)/(r/R0-2d0*m)
  dmdr = 4.*pi*r**2.*ee/M0
  
  f(1) = dpdr
  f(2) = dmdr
  
END SUBROUTINE funcs

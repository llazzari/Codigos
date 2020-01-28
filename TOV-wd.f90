MODULE global
  IMPLICIT NONE
  REAL :: pi = 4.*atan(1.0), me = 9.109d-31, c = 2.99792458d8, hbar = 1.055d-34, mn = 1.67262d-27 ![c] = m/s
  REAL :: G = 6.67408d-11, Msol = 1.989d30 ![G] = m^3kg^-1s^-2
  REAL :: e0 = 8.9798d-4 !e0 = me^4*c^5/(pi^2*hbar^3), [e0] = eV/fm^3
  REAL :: R0 = 1.47611    !R0 = G*Msol/c^2, [R0] = km
  REAL :: M0 = 8.9616d-13 !M0 = 10^54/(6,242d18*Msol*c**2)
  
END MODULE global

PROGRAM pressure
  USE global
  REAL :: p0, r, m ![p0] = eV/fm^3, [r] = km, [m] = x Msol
  REAL :: dr, p, kf, e, x, dp, mi ![e] = eV/fm^3
  
  p0 = 6.242197d-8 !1d19/(1.602d26)
  kf = 0.7d-22 ![kf] = kg*m/s
  x = kf/(me*c)
  dp = 1.2

  DO WHILE (p0 < 6.242197d12)
     p = p0
     r = 0d0 
     dr = 0.1  
     m = 0d0
     mi = 1d-3
     IF (p0 > 1d0) dr = 0.01
     DO WHILE (p > 0d0)
        CALL newton(p,kf,x)
        r = r + dr
        p = -R0*e(kf)*m*dr/r**2.*(1.+p/e(kf))*(1.+4.*pi*r**3.*p*M0/mi)/(1.-2.*R0*m/r) + p
        m = 4.*pi*r**2.*e(kf)*dr*M0 + m
        mi = m
        IF (r > 21000.) STOP
     END DO
     WRITE(20,'(3(E14.6))') p0, m, r
     p0 = p0*dp
  END DO
  
END PROGRAM pressure

REAL FUNCTION e(kf)
  USE global
  REAL :: x, kf, elec, enuc
  
  x = kf/(me*c)
  
  elec = e0/8. * ((2.*x**3. + x)*sqrt(1.+x**2.) - log(x + sqrt(x**2. + 1.)))
  enuc = kf**3.*2.*mn*c**2./(3.*hbar**3.*pi**2*1.602d26) !eV/fm^3
  
  e = enuc + elec
END FUNCTION e

SUBROUTINE newton(p,kf,x0)
  USE global
  REAL :: pkf, dpdx
  REAL :: x, x0, kf, p
  INTEGER :: K, Kmax
  
  K = 1
  Kmax = 200

  IF (p < 6.42197) x0 = 1.
  IF (p > 6.42197) x0 = 250.
  
  DO
     x = x0 - pkf(x0,p)/dpdx(x0)
     IF (abs(pkf(x0,p)) <= 1d-11) EXIT
     x0 = x
     K = K + 1
     IF (K > Kmax) EXIT 
 END DO

 kf = me*c*x
 
END SUBROUTINE newton

REAL FUNCTION pkf(x,p)
  USE global
  REAL :: x, p
  
  pkf = e0/24.*((2.*x**3. - 3.*x)*sqrt(1.+x**2.) + 3.*log(x+sqrt(x**2. + 1.))) - p
END FUNCTION pkf

REAL FUNCTION dpdx(x)
  REAL :: x
  
  dpdx = (2.9872d-4)*x**4./sqrt(x**2. + 1.)
END FUNCTION dpdx

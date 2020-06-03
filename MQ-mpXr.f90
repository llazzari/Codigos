!ESTRELA ESTRANHA
!Limite ultrarrelativístico (quarks sem massa)
!Programa para calcular a massa e a pressão como funções do raio, do centro da estrela até a a superfície, onde p = 0. Neste ponto, o programa para, imprimindo o raio total e a massa total da estrela na tela.
module global_variables
  !pi, velocidade da luz no vácuo, constante reduzida de Planck
  real, parameter :: pi = 4.*atan(1.0d0)
  !R0 = G*Msol/c^2 // G = 6.67408d-11, Msol = 1.989d30 ![G] = m^3kg^-1s^-2, [Msol] = kg, or converts from Msun to km
  real, parameter :: R0 = 1.4765 !km
  !conversion
  real, parameter :: c = 1.3237d-3 !conversion from gev/fm^3 (nu) to km^-2 (gu)
  !bag pressure [GeV/fm^3]
  REAL, parameter :: B = 0.075
  !baryon number saturation density [GeV^3]
  real, parameter :: n0 = 1.2277d-3 !=0.16 fm^-3
  !EoS
  real :: e, p
  !intial pressure, energy densities [GeV/fm^3] and baryon number density [GeV^3]
  real :: p0, e0, nc 
  !mass
  real :: m 
  !ode's [dpdr] = GeV/fm^3/km, [dmdr] = km^-1
  real :: dpdr, dmdr
  !generalized metric functions
  real :: lambda
  !star's radius [km]
  real :: Rf
  !speed of sound squared dp/de
  real, parameter :: cs2 = 1/3.
  
contains
  real function sgn(x)
    implicit none
    real :: x
    
    if (x > 0.) sgn = 1.
    if (x == 0.) sgn = 0.
    if (x < 0.) sgn = -1.
    
  end function sgn
end module global_variables

program TOV
  use global_variables
  implicit none
  !printing formats
12 format(2(e16.9, 1x))
13 format(3(e16.9, 1x))
23 format(3(a,16x))
  
  OPEN(10,FILE='MQ-MRxp0.dat')
  WRITE(10,23) '#p0','M','R'

  nc = 5*n0
  p0 = -B + 3*(pi*nc**2)**(2/3.)*130.32/4.
  e0 = cs2*(e0 + 4*B)
  call rk4o()
  write(*,12) Rf, M
  
end program tov

subroutine rk4o()
  use global_variables
  implicit none
  real, dimension(2) :: y0, y, k1, k2, k3, k4, f
  real, parameter :: dr = 1d-3
  real :: r
13 format(3(e16.9,1x)) 
  
  r = 1d-10
  y0(1) = p0
  y0(2) = 0d0

  do while (y0(1) > 0d0)
     write(10,13) r, y0(1)*1d3, y0(2)/R0
     call funcs(r,y0,f)
     k1 = f
     call funcs(r+dr/2.,y0+dr*k1/2.,f)
     k2 = f
     call funcs(r+dr/2.,y0+dr*k2/2.,f)
     k3 = f
     call funcs(r+dr,y0+dr*k3,f)
     k4 = f
     y = y0 + dr/6.*(k1 + 2*k2 + 2*k3 + k4)
     y0 = y
     r = r + dr
  end do
  M = y(2)/R0
  Rf = R
  
end subroutine rk4o

subroutine funcs(r,y,f)
  use global_variables
  implicit none
  real :: r
  real, dimension(2) :: y, f
  
  p = y(1)
  m = y(2)
  e = 3*p + 4*B

  lambda = -log(1 - 2*m/r)/2 !adimensional
  
  dpdr = -(e+p)*(4*pi*r*p*c + m/r**2)*exp(2*lambda)
  dmdr = 4.*pi*r**2.*e*c
  

  f(1) = dpdr
  f(2) = dmdr
  
end subroutine funcs


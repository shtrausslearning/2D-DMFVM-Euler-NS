!> Computes values of dependent variables (pressure, temperature, speed
!! of sound, specific heat ratio, specific heat coeff. at const. pressure)
!! from conservative variables at all grid points. Additionally, laminar
!! viscosity and heat conductivity coefficients are computed in the case
!! of viscous flow.

subroutine DependentVarsAll

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  implicit none

! local variables
  integer     :: i
  real(rtype) :: gam1, rgas, g1cp, rhoq, s1, s2, s12, rat, cppr

! *****************************************************************************

! R -> Specific Gas Constant
! R = Cp - Cv, gam = Cp/Cv

! Equation of State : p=rho*R*T (1) (Ideal Gas)
! Enthalpy : h = Cp*T (2)  
! H = E + (p/rho) / H = h + |v|^2 / 2

! Obtain P = (gam-1)*rho*(E - (u^2 + v^2)/2) (3)

  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma
  g1cp = gam1*cpgas

! Euler equations
  if (kequs == "E") then
    do i=1,nnodes
      rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i) ! vel
      dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i)) ! Static Pressure (3)
      dv(2,i) = dv(1,i)/(rgas*cv(1,i)) ! Static Temperature (from 1)
      dv(3,i) = Sqrt(g1cp*dv(2,i)) ! Ratio of Specific Heats 
      dv(4,i) = gamma ! Ratio of Specific Heats
      dv(5,i) = cpgas ! Specific Heats Coefficient
    enddo

! Navier-Stokes equations
  else
    s1   = 110.D0
    s2   = 288.16D0
    s12  = 1.D0 + s1/s2
    cppr = cpgas/prlam
    do i=1,nnodes
      rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
      dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i)) ! Static Pressure
      dv(2,i) = dv(1,i)/(rgas*cv(1,i)) ! Static Temperature 
      dv(3,i) = Sqrt(g1cp*dv(2,i)) ! Speed of Sound
      dv(4,i) = gamma ! Ratio of specific heats
      dv(5,i) = cpgas ! Specific Heat Coefficient @constP
      rat     = Sqrt(dv(2,i)/s2)*s12/(1.D0+s1/dv(2,i)) 
      dv(6,i) = refvisc*rat   ! Laminar Viscosity Coef
      dv(7,i) = dv(6,i)*cppr  ! Laminar Heat Conductivity Coef
    enddo
  endif

end subroutine DependentVarsAll

!> Computes values of dependent variables (pressure, temperature, speed
!! of sound, specific heat ratio, specific heat coeff. at const. pressure)
!! from conservative variables at the node i. Additionally, laminar
!! viscosity and heat conductivity coefficients are computed in the case
!! of viscous flow.
!!
!! @param i  node index

subroutine DependentVarsOne( i )

  use ModDataTypes
  use ModPhysics
  implicit none

! parameters
  integer, intent(in) :: i

! local variables
  real(rtype) :: gam1, rgas, g1cp, rhoq, s1, s2, s12, rat

! *****************************************************************************

  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma
  g1cp = gam1*cpgas

! Euler equations

  if (kequs == "E") then
    rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
    dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
    dv(2,i) = dv(1,i)/(rgas*cv(1,i))
    dv(3,i) = Sqrt(g1cp*dv(2,i))
    dv(4,i) = gamma
    dv(5,i) = cpgas

! Navier-Stokes equations

  else
    s1      = 110.D0
    s2      = 288.16D0
    s12     = 1.D0 + s1/s2
    rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
    dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
    dv(2,i) = dv(1,i)/(rgas*cv(1,i))
    dv(3,i) = Sqrt(g1cp*dv(2,i))
    dv(4,i) = gamma
    dv(5,i) = cpgas
    rat     = Sqrt(dv(2,i)/s2)*s12/(1.D0+s1/dv(2,i))
    dv(6,i) = refvisc*rat
    dv(7,i) = dv(6,i)*(cpgas/prlam)
  endif

end subroutine DependentVarsOne

!> Computes values of dependent variables (pressure, temperature, speed
!! of sound, specific heat ratio, specific heat coeff. at const. pressure)
!! from conservative variables at all grid points. Additionally, laminar
!! viscosity and heat conductivity coefficients are computed in the case
!! of viscous flow.
!!
subroutine DependentVarsAll

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  implicit none

! local variables
  integer     :: i
  real(rtype) :: gam1, rgas, g1cp, rhoq, s1, s2, s12, rat, cppr

! *****************************************************************************

  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma
  g1cp = gam1*cpgas

! Euler equations

  if (kequs == "E") then
    do i=1,nnodes
      rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
      dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
      dv(2,i) = dv(1,i)/(rgas*cv(1,i))
      dv(3,i) = Sqrt(g1cp*dv(2,i))
      dv(4,i) = gamma
      dv(5,i) = cpgas
    enddo

! Navier-Stokes equations

  else
    s1   = 110.D0
    s2   = 288.16D0
    s12  = 1.D0 + s1/s2
    cppr = cpgas/prlam
    do i=1,nnodes
      rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
      dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
      dv(2,i) = dv(1,i)/(rgas*cv(1,i))
      dv(3,i) = Sqrt(g1cp*dv(2,i))
      dv(4,i) = gamma
      dv(5,i) = cpgas
      rat     = Sqrt(dv(2,i)/s2)*s12/(1.D0+s1/dv(2,i))
      dv(6,i) = refvisc*rat
      dv(7,i) = dv(6,i)*cppr
    enddo
  endif

end subroutine DependentVarsAll

! =============================================================================

!> Computes values of dependent variables (pressure, temperature, speed
!! of sound, specific heat ratio, specific heat coeff. at const. pressure)
!! from conservative variables at the node i. Additionally, laminar
!! viscosity and heat conductivity coefficients are computed in the case
!! of viscous flow.
!!
!! @param i  node index
!!
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


!> Computes reference values of limited variables (density, u, v, pressure)
!! and of the control volumes. The reference values are used to normalize
!! variables within the limiter functions (Roe's upwind scheme).
!!
subroutine LimiterRefvals

  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none

! local variables
  real(rtype) :: gam1, rgas, temp, rho, cs, mach

! *****************************************************************************
! reference volume (= largest control volume)

  volref = Maxval( vol )

! reference density, velocity and pressure

  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma

! external flow

  if (kflow == "E") then
    limref(1) = rhoinf
    limref(2) = Sqrt(uinf*uinf+vinf*vinf)
    limref(3) = limref(2)
    limref(4) = pinf

! internal flow

  else
    temp      = ttinl*(pout/ptinl)**(gam1/gamma)
    rho       = pout/(rgas*temp)
    cs        = Sqrt(gamma*pout/rho)
    mach      = Sqrt(2.D0*((ttinl/temp)-1.D0)/gam1)
    limref(1) = rho
    limref(2) = mach*cs
    limref(3) = limref(2)
    limref(4) = pout
  endif

end subroutine LimiterRefvals

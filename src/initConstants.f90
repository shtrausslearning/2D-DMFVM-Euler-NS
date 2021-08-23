!> Initializes constants used by the solver.
!!
subroutine InitConstants

  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none

! local variables
  real(rtype) :: gam1, rgas

! *****************************************************************************

  pi  = 4.D0*Atan(1.D0)
  rad = 180.D0/pi

  if (kflow == 'E') then

! - external flow; it is assumed that "gamma" and "cpgas" specified in
!   the input file are valid for the complete far-field boundary

    gam1 = gamma - 1.D0
    rgas = gam1*cpgas/gamma

    alpha   = alpha/rad
    rhoinf  = pinf/(rgas*tinf)
    qinf    = machinf * Sqrt(gam1*cpgas*tinf)
    uinf    = qinf*Cos(alpha)
    vinf    = qinf*Sin(alpha)
    refrho  = rhoinf
    refvel  = qinf
    if (kequs == 'N') then
      refvisc = rhoinf*qinf*cref/renum
    else
      refvisc = 0.D0
    endif

  else

! - internal flow

    betainl = betainl/rad
    betaout = betaout/rad
    if (kequs == 'N') then
      refvisc = refrho*refvel*cref/renum
    else
      refvisc = 0.D0
    endif

  endif

end subroutine InitConstants

!> Applies inflow boundary condition at dummy points (subsonic flow assumed).
!!
!! @param ibegn  indirect pointer to first node of the boundary
!! @param iendn  indirect pointer to last node of the boundary
!!
subroutine BcondInflow( ibegn,iendn )

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModInterfaces, only : DependentVarsOne
  implicit none

! parameters
  integer, intent(in) :: ibegn, iendn

! local variables
  integer     :: ib, ibn, idn, ie
  real(rtype) :: ds, sxn, syn, rrho, u, v, uabs, unorm, cosa, c02, rinv, &
                 dis, cb, cc02, tb, pb, rhob, uabsb, ub, vb, gam1, ggm1, &
                 rgas

! *****************************************************************************

  do ib=ibegn,iendn

    ibn = bnode(1,ib)         ! boundary node
    idn = bnode(2,ib)         ! dummy node
    ie  = bnode(3,ib)         ! edge to dummy node

    ds  = Sqrt(sij(1,ie)**2+sij(2,ie)**2)
    sxn = sij(1,ie)/ds
    syn = sij(2,ie)/ds

    gam1  = dv(4,ibn) - 1.D0
    ggm1  = dv(4,ibn)/gam1
    rgas  = gam1*dv(5,ibn)/dv(4,ibn)
    rrho  = 1.D0/cv(1,ibn)
    u     = cv(2,ibn)*rrho
    v     = cv(3,ibn)*rrho
    uabs  = Sqrt(u*u+v*v)
    unorm = u*sxn + v*syn
    if (uabs < 1.D-20) then
      cosa = 1.D0
    else
      cosa = -unorm/uabs
    endif

    c02  = dv(3,ibn)*dv(3,ibn) + 0.5D0*gam1*uabs*uabs
    rinv = unorm - 2.D0*dv(3,ibn)/gam1
    dis  = (gam1*cosa*cosa+2.D0)*c02/(gam1*rinv*rinv) - 0.5D0*gam1
    if (dis < 0.D0) then
      write(*,1000) ibn,dis
      dis = 1.D-20
    endif
    cb    = -rinv*(gam1/(gam1*cosa*cosa+2.D0))*(1.D0+cosa*Sqrt(dis))
    cc02  = Min(cb*cb/c02, 1.D0)
    tb    = ttinl*cc02
    pb    = ptinl*(tb/ttinl)**ggm1
    rhob  = pb/(rgas*tb)
    uabsb = 2.D0*dv(5,ibn)*(ttinl-tb)
    uabsb = Sqrt(uabsb)
    ub    = uabsb*Cos(betainl)
    vb    = uabsb*Sin(betainl)

    cv(1,idn) = rhob
    cv(2,idn) = rhob*ub
    cv(3,idn) = rhob*vb
    cv(4,idn) = pb/gam1 + 0.5D0*rhob*(ub*ub+vb*vb)

    call DependentVarsOne( idn )

  enddo  ! ib

1000  format(" Warning (BcondInflow): discriminant<0 at boundary node ",I3, &
             ", d= ",1PE13.5)

end subroutine BcondInflow

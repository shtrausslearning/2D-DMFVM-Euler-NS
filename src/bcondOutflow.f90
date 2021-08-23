!> Applies outflow boundary condition to dummy points. Sub- or supersonic
!! outflow is possible.
!!
!! @param ibegn  indirect pointer to first node of the boundary
!! @param iendn  indirect pointer to last node of the boundary
!!
subroutine BcondOutflow( ibegn,iendn )

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModInterfaces, only : DependentVarsOne
  implicit none

! parameters
  integer, intent(in) :: ibegn, iendn

! local variables
  integer     :: ib, ibn, idn, ie
  real(rtype) :: gam1, ds, sxn, syn, rrho, u, v, rrhoc, deltp, &
                 rhob, ub, vb, vnd, q, mach

! *****************************************************************************

  do ib=ibegn,iendn

    ibn = bnode(1,ib)         ! boundary node
    idn = bnode(2,ib)         ! dummy node
    ie  = bnode(3,ib)         ! edge to dummy node

    ds  = Sqrt(sij(1,ie)**2+sij(2,ie)**2)
    sxn = sij(1,ie)/ds
    syn = sij(2,ie)/ds

    gam1  = dv(4,ibn) - 1.D0
    rrho  = 1.D0/cv(1,ibn)
    u     = cv(2,ibn)*rrho
    v     = cv(3,ibn)*rrho
    q     = Sqrt(u*u+v*v)
    mach  = q/dv(3,ibn)

    if (mach < 1.D0) then
      rrhoc = rrho/dv(3,ibn)
      deltp = dv(1,ibn) - pout
      rhob  = cv(1,ibn) - deltp/(dv(3,ibn)*dv(3,ibn))
      ub    = u + sxn*deltp*rrhoc
      vb    = v + syn*deltp*rrhoc

! --- special treatment to prevent "deltp" from changing the sign
!     of velocity components. This may happen for very small u, v.

      vnd = ub*sxn + vb*syn
      if (vnd < 0.D0) then
        ub = Sign(1.D0,u)*Max(Abs(ub),Abs(u))
        vb = Sign(1.D0,v)*Max(Abs(vb),Abs(v))
      endif
      cv(1,idn) = rhob
      cv(2,idn) = rhob*ub
      cv(3,idn) = rhob*vb
      cv(4,idn) = pout/gam1 + 0.5D0*rhob*(ub*ub+vb*vb)
    else
      cv(1,idn) = cv(1,ibn)
      cv(2,idn) = cv(2,ibn)
      cv(3,idn) = cv(3,ibn)
      cv(4,idn) = cv(4,ibn)
    endif

    call DependentVarsOne( idn )

  enddo  ! ib

end subroutine BcondOutflow

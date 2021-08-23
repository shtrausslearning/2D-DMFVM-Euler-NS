!> Computes pressure forces and moments acting on the body. The contributions
!! are summed up by looping over ALL walls. Subroutine also computes lift and
!! drag coefficients.
!!
subroutine Forces

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  implicit none

! local variables
  integer     :: ib, ibf, ibegf, iendf, n1, n2
  real(rtype) :: sx, sy, pwall, cp, xa, ya, dcx, dcy, cx, cy

! *****************************************************************************
! initialize force coefficients

  cx = 0.D0
  cy = 0.D0
  cm = 0.D0

! loop over boundaries searching for walls

  ibegf = 1

  do ib=1,nsegs
    iendf = ibound(1,ib)

    if (btype(ib)>=300 .and. btype(ib)<500) then
      do ibf=ibegf,iendf
        n1    = bface(1,ibf)
        n2    = bface(2,ibf)
        sx    = sbf(1,ibf)
        sy    = sbf(2,ibf)
        pwall = 0.5D0*(dv(1,n1)+dv(1,n2))
        cp    = 2.D0*(pwall-pinf)/(rhoinf*qinf*qinf)
        xa    = (0.5D0*(x(n1)+x(n2))-xref)/cref
        ya    = (0.5D0*(y(n1)+y(n2))-yref)/cref
        dcy   = sy*cp
        dcx   = sx*cp
        cy    = cy + dcy
        cx    = cx + dcx
        cm    = cm + dcx*ya - dcy*xa
      enddo
    endif ! btype

    ibegf = iendf + 1
  enddo   ! ib

! final lift and drag coefficients (pressure forces only!)

  cl = cy*Cos(alpha) - cx*Sin(alpha)
  cd = cy*Sin(alpha) + cx*Cos(alpha)

end subroutine Forces

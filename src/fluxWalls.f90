!> Computes convective fluxes in the normal direction at solid walls.
!!
subroutine FluxWalls

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModNumerics
  implicit none

! local variables
  integer     :: i, j, ib, ibf, ibegf, iendf
  real(rtype) :: sx, sy, pl, pr

! *****************************************************************************

  ibegf = 1

  do ib=1,nsegs
    iendf = ibound(1,ib)
    if (btype(ib)>=300 .and. btype(ib)<500) then
      do ibf=ibegf,iendf
        i        = bface(1,ibf)
        j        = bface(2,ibf)
        sx       = sbf(1,ibf)/12.D0
        sy       = sbf(2,ibf)/12.D0
        pl       = 5.D0*dv(1,i) +      dv(1,j)
        pr       =      dv(1,i) + 5.D0*dv(1,j)
        rhs(2,i) = rhs(2,i) + sx*pl
        rhs(3,i) = rhs(3,i) + sy*pl
        rhs(2,j) = rhs(2,j) + sx*pr
        rhs(3,j) = rhs(3,j) + sy*pr
      enddo
    endif
    ibegf = iendf + 1
  enddo

end subroutine FluxWalls

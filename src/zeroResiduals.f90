!> Zeros out normal component of the residual at symmetry and
!! at no-slip boundaries.
!!
subroutine ZeroResiduals

  use ModGeometry
  use ModNumerics
  use ModPhysics
  implicit none

! local variables
  integer :: i, ib, ibn, idir, ibegn, iendn

! *****************************************************************************

  ibegn = 1

  do ib=1,nsegs
    iendn = ibound(2,ib)

! - symmetry boundary

    if (btype(ib)>=500 .and. btype(ib)<600) then
      if (btype(ib)-500 <  2) idir = 2  ! x=const. line -> x-component
      if (btype(ib)-500 >= 2) idir = 3  ! y=const. line -> y-component
      do ibn=ibegn,iendn
        i           = bnode(1,ibn)
        rhs(idir,i) = 0.D0
      enddo

! - viscous (no-slip) wall

    else if ((btype(ib)>=300 .and. btype(ib)<400) .and. kequs=="N") then
      do ibn=ibegn,iendn
        i        = bnode(1,ibn)
        rhs(2,i) = 0.D0       ! velocity components = 0
        rhs(3,i) = 0.D0
      enddo
    endif

    ibegn = iendn + 1
  enddo

end subroutine ZeroResiduals

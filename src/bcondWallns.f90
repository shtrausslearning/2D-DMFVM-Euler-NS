!> Applies no-slip wall boundary condition. Adiabatic walls
!! only are assumed (velocity components are zeroed out).
!!
!! @param ibegn  indirect pointer to first node of the boundary
!! @param iendn  indirect pointer to last node of the boundary
!!
subroutine BcondWallns( ibegn,iendn )

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModInterfaces, only : DependentVarsOne
  implicit none

! parameters
  integer, intent(in) :: ibegn, iendn

! local variables
  integer :: ib, ibn

! *****************************************************************************

  do ib=ibegn,iendn

    ibn       = bnode(1,ib)   ! boundary node
    cv(2,ibn) = 0.D0
    cv(3,ibn) = 0.D0

    call DependentVarsOne( ibn )

  enddo

end subroutine BcondWallns

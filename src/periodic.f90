!> Adds both parts of variable var(:,:) at all periodic boundaries.
!!
!! @param var  variable to update
!!
subroutine Periodic( var )

  use ModDataTypes
  use ModGeometry
  implicit none

! parameters
  real(rtype) :: var(:,:)

! local variables
  integer :: i, j, ib, ibn, ibegn, iendn, n

! *****************************************************************************

  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if (btype(ib)>=700 .and. btype(ib)<800) then
      do n=1,Ubound(var,1)
        do ibn=ibegn,iendn
          i        = bnode(1,ibn)
          j        = bnode(2,ibn)
          var(n,i) = var(n,i) + var(n,j)
          var(n,j) = var(n,i)
        enddo
      enddo
    endif
    ibegn = iendn + 1
  enddo

end subroutine Periodic

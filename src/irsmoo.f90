!> Applies the central implicit smoothing method to the residuals
!! by using Jacobi iteration.
!!
!! @param ncontr  number of contributions to a node
!! @param rhsold  storage for the initial non-smoothed residual
!! @param rhsit   nodal contributions to the residual during the iteration
!!
subroutine Irsmoo( ncontr,rhsold,rhsit )

  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : Periodic
  implicit none

! parameters
  integer     :: ncontr(:)
  real(rtype) :: rhsold(:,:), rhsit(:,:)

! local variables
  integer     :: itirs, i, j, ie, ib, ibn, ibegn, iendn
  real(rtype) :: den

! *****************************************************************************
! initialize counter (no. of contributions to a node);
! store the non-smoothed residual in "rhsold"

  do i=1,nndint
    ncontr(i)   = 0
    rhsold(1,i) = rhs(1,i)
    rhsold(2,i) = rhs(2,i)
    rhsold(3,i) = rhs(3,i)
    rhsold(4,i) = rhs(4,i)
  enddo

! Jacobi iteration ------------------------------------------------------------

  do itirs=1,nitirs

! - zero out nodal contributions
    do i=1,nndint
    rhsit(1,i) = 0.D0
    rhsit(2,i) = 0.D0
    rhsit(3,i) = 0.D0
    rhsit(4,i) = 0.D0
    enddo

! - loop over edges - first iteration => set counter
    if (itirs == 1) then
    do ie=1,nedint
    i          = edge(1,ie)
    j          = edge(2,ie)
    ncontr(i)  = ncontr(i) + 1
    ncontr(j)  = ncontr(j) + 1
    rhsit(1,i) = rhsit(1,i) + rhs(1,j)
    rhsit(2,i) = rhsit(2,i) + rhs(2,j)
    rhsit(3,i) = rhsit(3,i) + rhs(3,j)
    rhsit(4,i) = rhsit(4,i) + rhs(4,j)

    rhsit(1,j) = rhsit(1,j) + rhs(1,i)
    rhsit(2,j) = rhsit(2,j) + rhs(2,i)
    rhsit(3,j) = rhsit(3,j) + rhs(3,i)
    rhsit(4,j) = rhsit(4,j) + rhs(4,i)
    enddo

! - loop over edges - without counter
    else
    do ie=1,nedint
    i          = edge(1,ie)
    j          = edge(2,ie)
    rhsit(1,i) = rhsit(1,i) + rhs(1,j)
    rhsit(2,i) = rhsit(2,i) + rhs(2,j)
    rhsit(3,i) = rhsit(3,i) + rhs(3,j)
    rhsit(4,i) = rhsit(4,i) + rhs(4,j)

    rhsit(1,j) = rhsit(1,j) + rhs(1,i)
    rhsit(2,j) = rhsit(2,j) + rhs(2,i)
    rhsit(3,j) = rhsit(3,j) + rhs(3,i)
    rhsit(4,j) = rhsit(4,j) + rhs(4,i)
    enddo
    endif    ! itirs > 1


! - new smoothed residual
    do i=1,nndint
    den      = 1.D0/(1.D0+epsirs*Real(ncontr(i)))
    rhs(1,i) = (rhsit(1,i)*epsirs+rhsold(1,i))*den
    rhs(2,i) = (rhsit(2,i)*epsirs+rhsold(2,i))*den
    rhs(3,i) = (rhsit(3,i)*epsirs+rhsold(3,i))*den
    rhs(4,i) = (rhsit(4,i)*epsirs+rhsold(4,i))*den
    enddo

  enddo    ! loop over itirs

end subroutine Irsmoo

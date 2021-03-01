
!> Generates final edge list (interior edges, edges to dummy nodes) using
!! the temporary lists "niedge" and "iedge".
!!
!! @param niedge  pointer from a node to iedge()
!! @param iedge   linked list of edge endpoints:
!!                @li (1,*) = point j of edge (i,j)
!!                @li (2,*) = next point j which is also connected to i;
!!                            if <0 - no further connections
!!                @li (3,*) = pointer to edge() - used in InitMetrics to associate
!!                            face vector sij() with the correct edge
!!
subroutine EdgesFinalize( niedge,iedge )

  use ModGeometry
  use ModNumerics
  use ModInterfaces, only : ErrorMessage
  implicit none

! parameters
  integer :: niedge(:), iedge(:,:)

! local variables
  integer :: errFlag, i, ibn, ie, cedge

! *****************************************************************************

  allocate( edge(2,nedges),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for edge()" )

  ie = 0  ! edge counter

  do i=1,nndint

! - loop over all grid nodes

    cedge = niedge(i)
    if (cedge > 0) then
10    continue

! ----- loop over all edges adjacent to node "i"

        ie             = ie + 1
        edge(1,ie)     = i
        edge(2,ie)     = iedge(1,cedge)
        iedge(3,cedge) = ie                 ! we need it in InitMetrics
        cedge          = iedge(2,cedge)     ! next adjacent edge
        if (cedge < 0) goto 20
      goto 10
20    continue
    endif

  enddo

  if (ie /= nedint) then
    call ErrorMessage( "did not get the correct number of interior edges" )
  endif

! add edges to dummy nodes;
! store 'dummy' edges in "bnode(3,*)"

  do ibn=1,nbnodes
    if (bnode(3,ibn) == -1) then      ! dummy node here (see DummyNodes)
      ie           = ie + 1
      edge(1,ie)   = bnode(1,ibn)     ! boundary node first
      edge(2,ie)   = bnode(2,ibn)     ! dummy node second
      bnode(3,ibn) = ie
    endif
  enddo

  if (ie /= nedges) then
    call ErrorMessage( "did not get the correct number of dummy edges" )
  endif

end subroutine EdgesFinalize

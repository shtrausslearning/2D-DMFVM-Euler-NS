! Generation of dummy nodes at inlet, outlet and far-field boundaries.
subroutine DummyNodes

  use ModGeometry
  use ModInterfaces, only : ErrorMessage
  implicit none

! local variables
  character(chrlen) :: msg
  logical :: flag
  integer :: errFlag, ibegf, iendf, ibegn, iendn, itype
  integer :: i, ib, ibf, idn
  integer, allocatable :: marker(:)  ! node marker

! *****************************************************************************

  allocate( marker(nndint),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate temporary marker" )

! loop over boundary segments

  ibegf = 1
  ibegn = 1
  idn   = 0  ! counter of dummy nodes

  do ib=1,nsegs

    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
    itype = btype(ib)
    flag  = .false.  ! true for inlet/outlet/far-field
    if (itype>=100 .and. itype<200) flag = .true.
    if (itype>=200 .and. itype<300) flag = .true.
    if (itype>=600 .and. itype<700) flag = .true.

    if (itype<700 .or. itype>=800) then   ! NOT a periodic boundary

! --- reset node marker

      do i=1,nndint
        marker(i) = -777
      enddo

! --- loop over faces of boundary "ib" and mark nodes

      do ibf=ibegf,iendf
        marker(bface(1,ibf)) = 1
        marker(bface(2,ibf)) = 1
      enddo

! --- store node indexes in "bnode";
!     count dummy nodes (idn) for inlet/outlet/far-field boundary

      do i=1,nndint
        if (marker(i) == 1) then           ! must be on boundary
          if (ibegn > nbnodes) then        ! check dimension
            call ErrorMessage( "max. no. of boundary nodes exceeded" )
          endif
          if (flag) then                   ! *** inlet/outlet/far-field
            idn            = idn + 1
            bnode(1,ibegn) = i             ! index of boundary node
            bnode(2,ibegn) = nndint + idn  ! index of dummy node
            bnode(3,ibegn) = -1            ! set in "EdgesFinalize"
          else                             ! *** other boundary type
            bnode(1,ibegn) = i             ! index of boundary node
          endif
          ibegn = ibegn + 1                ! count boundary nodes
        endif
      enddo

! --- check no. of boundary nodes

      if (ibegn-1 .ne. iendn) then
        write(msg,1000) itype,iendn,ibegn-1
        call ErrorMessage( msg )
      endif

    endif  ! periodic?

! - reset pointers to faces and to nodes

    ibegf = iendf + 1
    ibegn = iendn + 1

  enddo  ! ib

  deallocate( marker )

! set total number of nodes (add dummy nodes)

  nnodes = nndint + idn

1000 format("no. of nodes for boundary ",I3," is wrong. It should be ",I5, &
            " but it is ",I5)

end subroutine DummyNodes

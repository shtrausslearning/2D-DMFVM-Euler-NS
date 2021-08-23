
subroutine ReadGrid

  use ModFiles
  use ModGeometry
  use ModInterfaces, only : DummyNodes, ErrorMessage
  implicit none

! local variables
  integer :: errFlag, i, ib, ibn, ibf, ibegf, iendf, ibegn, iendn, inn
  integer :: wtid,ierr,ii,nwbfnode,j,k,ic
  
! temporary arrays
  integer, allocatable :: wbfsort(:) ! wall bc from bface
  integer, allocatable :: res(:)       ! output for bsort
  integer, allocatable :: wbfnode(:)

! *****************************************************************************

  open(unit=ifGrid, file="./grids/"//fnGrid, status="old", action="read", iostat=errFlag)
  if (errFlag /= 0) call ErrorMessage( "cannot open grid file" )
  
  write(*,*) '' 
  write(*,*) 'Grid File: '//trim(fnGrid)

  read(ifGrid,"(1X)")
  read(ifGrid,"(1X)")
  read(ifGrid,"(1X)")

! numbers of physical nodes, triangles and boundary segments

  read(ifGrid,*) nndint,ntria,nsegs

! boundary type, no. of boundary faces & nodes, boundary name

  allocate( btype(nsegs),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for btype()" )

  allocate( bname(nsegs),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for bname()" )

  allocate( ibound(2,nsegs),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for ibound()" )

  read(ifGrid,"(1X)")
  do ib=1,nsegs
    read(ifGrid,  *  ) btype(ib),ibound(1,ib),ibound(2,ib)
    read(ifGrid,"(A)") bname(ib)
  enddo

  nbfaces = ibound(1,nsegs)
  nbnodes = ibound(2,nsegs)

! definition of boundary faces / periodic nodes

  allocate( bnode(3,nbnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for bnode()" )

  allocate( bface(2,nbfaces),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for bface()" )

  do ibn=1,nbnodes
    bnode(1,ibn) = -777
    bnode(2,ibn) = -777      ! set in DummyNodes
    bnode(3,ibn) = -777      ! set in EdgesFinalize
  enddo
  do ibf=1,nbfaces
    bface(1,ibf) = -777
    bface(2,ibf) = -777
  enddo

  read(ifGrid,"(1X)")
  ibegf = 1
  ibegn = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
    if (btype(ib)>=700 .and. btype(ib)<800) then   ! periodic nodes
      do ibn=ibegn,iendn
        read(ifGrid,*) bnode(1,ibn),bnode(2,ibn)
      enddo
    else                                           ! boundary faces
      do ibf=ibegf,iendf
        read(ifGrid,*) bface(1,ibf),bface(2,ibf)
      enddo
    endif
    ibegf = iendf + 1
    ibegn = iendn + 1
  enddo

! check boundary faces pointer

  do ibf=1,nbfaces
    if (bface(1,ibf)<0 .or. bface(2,ibf)<0) then
      call ErrorMessage( "array bface() not completely defined" )
    endif
  enddo

! generate dummy nodes

  call DummyNodes

! grid nodes

  allocate( x(nnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for x()" )

  allocate( y(nnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for y()" )

  read(ifGrid,"(1X)")
  do i=1,nndint
    read(ifGrid,*) x(i),y(i)
  enddo

! triangles

  allocate( tria(3,ntria),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for tria()" )

  read(ifGrid,"(1X)")
  do i=1,ntria
    read(ifGrid,*) tria(1,i),tria(2,i),tria(3,i)
  enddo

  close(unit=ifGrid)
  
  allocate( ntype(nndint),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot open ptype file" )
  
 wtid  = 0
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if (btype(ib)>=300 .and. btype(ib)<500) then   ! wall
   do ibf=ibegf,iendf
   wtid = wtid + 1
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo
  
  allocate( wbfsort(2*wtid),stat=ierr )
  if (ierr /= 0) call ErrorMessage( "cannot allocate memory for wbfsort()" )
  
  ii=0
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if (btype(ib)>=300 .and. btype(ib)<500) then   ! wall
   do ibf=ibegf,iendf
   ii=ii+1
   wbfsort(ii) = bface(1,ibf)
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo
  
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if (btype(ib)>=300 .and. btype(ib)<500) then   ! wall
   do ibf=ibegf,iendf
   ii=ii+1
   wbfsort(ii) = bface(2,ibf)
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo

! remove duplicates of wbfsort array
  allocate( res(2*wtid),stat=ierr )
  if (ierr /= 0) call ErrorMessage( "cannot allocate memory for res()" )
  
  k = 1
  res(1) = wbfsort(1)
  outer: do i=2,2*wtid
     do j=1,k
        if (res(j) == wbfsort(i)) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     res(k) = wbfsort(i)
     
  end do outer
  
! number of wall boundary nodes
  nwbfnode = k
    
  allocate( wbfnode(nwbfnode),stat=ierr )
  if (ierr /= 0) call ErrorMessage( "cannot allocate memory for wbfnode()" )
  
  do i=1,k
  wbfnode(i) = res(i)
  enddo
  
  deallocate( res )
  deallocate( wbfsort )
  
  ntype(:) = 0
  
! physical wallnodes
  do i=1,nwbfnode
  ic = wbfnode(i)
  do j=1,nndint
  if( j .eq. ic )then
  ntype(j) = 1 
  endif
  enddo
  enddo
  
  

end subroutine ReadGrid

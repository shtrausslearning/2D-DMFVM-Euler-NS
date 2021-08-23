
    Subroutine Gmsh_read
    use ModDataTypes
    use ModPlotQuant
    use ModFiles
    implicit none
!   Read and Convert Gmsh[msh] to urg grid format
    integer        :: ierr,i,phynod,imoji,et,nwalls,j,nv,k,kk,jj
    integer        :: nsegs,ugrt,status
    character(128) :: moji  
    character(128) :: fname2
    double precision, allocatable :: tcoord(:,:)  ! coordinte array
    integer,allocatable :: segmtype(:)          ! boundary type
    integer,allocatable :: ugr(:,:)             ! ugr bc output array
    
    integer,allocatable :: bcarray(:,:)         ! bc type array elements
    integer,allocatable :: triarray(:,:)        ! cell type array elements
    real(rtype) :: dpp
    integer :: trans
    
!   read input 
    write(*,*) 'converting gmsh format...'
    write(*,*) 'enter .msh input'
    read(*,"(A)") fname2

    call tempopen(fname2,trans)  ! slower method, read whole input file and find transition
                                                    ! between bc edge and triangular elements   
!   ##################################################################################
    ierr=0
!!    open(99, file="./mshFILE"//fname2, status="old", action="read", iostat=ierr)
    open(99, file="./mshFILE/"//fname2,form='formatted')
!!    if(ierr /= 0) call ErrorMessage( "cannot open input file" )
    read(99,*) moji            ! $MeshFormat  
    read(99,*) dpp
    read(99,*) moji            ! $EndMeshFormat
    
!   read physical groups : BC groupies
    read(99,*) moji            ! $PhysicalNames
    read(99,*) nsegs           ! number of physical groups [ includes fluid ]
    
    allocate( segmtype(nsegs),stat=ierr )
    if(ierr /= 0) call ErrorMessage( "cannot allocate segmtype" )
    
    do i=1,nsegs
    read(99,*) imoji,segmtype(i),moji
    enddo

    read(99,*) moji           ! $EndPhysicalNames
    read(99,*) moji            ! $Nodes
    read(99,*) phynod          ! number of nodes in physical domain
   
    ierr=0
    allocate( tcoord(2,phynod),stat=ierr )
    if(ierr /= 0) call ErrorMessage( "cannot allocate coord" )
    
    do i=1,phynod
    read(99,*) imoji,tcoord(1,i),tcoord(2,i),imoji
    enddo
    
    read(99,*) moji            ! $EndNodes
    read(99,*) moji            ! $Elements
    read(99,*) et              ! number of elements [ boundary edge + inner cells/triangles ]
    
    !ierr=0
    !allocate( earray(et,8),stat=ierr )
    !if(ierr /= 0) call ErrorMessage( "cannot allocate earray" ) ! old
    
!   allocate bc element array
    ierr=0
    allocate( bcarray(trans-1,7),stat=ierr )
    if(ierr /= 0) call ErrorMessage( "cannot allocate bcarray" )
    
!   allocate tria element array
    ierr=0
    allocate( triarray(et-trans+1,8),stat=ierr )
    if(ierr /= 0) call ErrorMessage( "cannot allocate triarray" )
    
!   reset
    triarray(:,:) = 0 ; bcarray(:,:) = 0
    !earray(:,:) = 0 ! resets old
    
    
!   Initial Read All Element Types ! Max 8 Columns, 7 for bc faces
    
!   bc edges
    do i=1,trans-1                        
    read(99,*) (bcarray(i,j),j=1,7)
    enddo

!   triangular elements
    do i=1,et-(trans-1)
    read(99,*) (triarray(i,j),j=1,8)
    enddo

    read(99,*) moji           ! $EndElements
    close(99)                 ! close gmsh file
!   ##################################################################################
    
    ierr=0
    allocate( ugr(nsegs-1,3),stat=ierr ) ! allocate excluding fluid triangles
    if(ierr /= 0) call ErrorMessage( "cannot allocate ugr " )

    ugrt = 0              ! ugr bc counter 
    j = 0
    
    do k=1,nsegs-1
    
!   general ugr counter
    ugrt = ugrt+1         ! ugr counter
       j = j + 1          ! bc counter
    
!   find what is next in order of bc faces
    nv = bcarray(j,4)   ! next bc type
    
    j=0
    do i=1,trans-1  
    if( bcarray(i,4) .eq. nv  )then
    j = bcarray(i,1) ! updates every time inside req. values
    endif
    enddo
    
    ugr(ugrt,1) = nv    ! next bc type
    ugr(ugrt,2) = j     ! last element line number
    ugr(ugrt,3) = j     ! only valid if no 700 bc type [ periodic ]
    
    enddo

!   #################################################################
    
!   write ugr file
    open(199,file='./grids/'//fnGrid,form='formatted')
    write(199,*) '#',trim(title)
    write(199,*) '#'
    write(199,*) '# no. of nodes, cells, boundaries:'
    write(199,"(3I8)") phynod,et-(trans-1),nsegs-1
    write(199,*) ' # boundaries: type, last face, last node, name: '
    do k=1,nsegs-1
    write(199,'(3I8)') ugr(k,1),ugr(k,2),ugr(k,3)
    if( ugr(k,1) .eq. 400 )then
    write(199,*) 'wall'
    elseif( ugr(k,1) .eq. 100 )then
    write(199,*) 'inflow'
    elseif( ugr(k,1) .eq. 200 )then
    write(199,*) 'outflow'
    elseif( ugr(k,1) .ge. 500 .and. ugr(k,1) .lt. 600  )then
    write(199,*) 'symmetry'
    elseif( ugr(k,1) .ge. 600 .and. ugr(k,1) .lt. 700 )then
    write(199,*) 'farfield'
    endif
    enddo
    
    write(199,*) '# boundary faces / periodic nodes (node1, node2):'
    
    do i=1,trans-1
    write(199,'(2I8)') bcarray(i,6),bcarray(i,7)
    enddo
    
    write(199,*) '# coordinates (x,y):'
    do i=1,phynod
    write(199,'(2E17.9)') tcoord(1,i),tcoord(2,i)
    enddo    
    
    write(199,*) '# triangles (node1, node2, node3):'
    do i=1,et-(trans-1)   ! should only be 777 nodes 
    write(199,'(3I8)') triarray(i,6),triarray(i,7),triarray(i,8)
    enddo
    close(199)
    

contains

    subroutine Usage2
    write(*,"(/,A,/)") "Usage:"
    write(*,"(A,/)")   "Unstruct2D <input file>"
    stop
    end subroutine Usage2
   
    End Subroutine Gmsh_read
    
    subroutine tempopen(fname2,trans)
    implicit none

    integer,allocatable :: temp(:,:)
    character(128) :: moji  
    integer :: imoji,ierr,moji2
    integer :: numb,i,ii,trans
    character(128), intent(in) :: fname2
    real dpp

    ierr = 0
    open(200,file='./mshFILE/'//fname2,form='formatted')
    if(ierr /= 0) call ErrorMessage( "cannot open input file" )
    
    read(200,*) moji
    read(200,*) dpp
    read(200,*) moji
    read(200,*) moji
    read(200,*) moji2
    
    do i=1,moji2
    read(200,*) imoji
    enddo
    
    read(200,*) moji ! $EndPhysicalNames
    read(200,*) moji ! $Nodes
    
    read(200,*) moji2
    do i=1,moji2
    read(200,*) imoji
    enddo
    
    read(200,*) moji ! $EndNodes
    read(200,*) moji  ! $Elements
    read(200,*) moji2 ! elements
    
    ierr=0
    allocate( temp(2,moji2),stat=ierr )
    if(ierr /= 0) call ErrorMessage( "cannot allocate temp " )
  
    do i=1,moji2
    read(200,*) temp(1,i),imoji,imoji,temp(2,i),imoji,imoji,imoji
    enddo

    read(200,*) moji
    close(200)
    
    do i=1,moji2
    if( temp(2,i) .eq. 777 )then
    ii=temp(1,i)
    goto 10
    endif
    enddo
    
10  continue

    trans = ii

    !write(*,*) trans
    !pause
    
    deallocate( temp )
    
    return
    end subroutine tempopen
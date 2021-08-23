
program Unstruct2D

  use ModDataTypes
  use ModControl
  use ModFiles
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces
  implicit none

! local variables
  character(chrlen) :: fname               ! filename (user input, convergence)
  integer           :: errFlag,ierr             ! error flag
  integer, allocatable     :: niedge(:), & ! temporary edge structures (DO NOT overwrite
                              iedge(:,:)   ! between EdgesInitialize and InitMetrics)
  integer, allocatable     :: iwork(:)     ! integer work space (used by Irsmoo)
  real(rtype), allocatable :: work(:)      ! real work space (used inside Solver)
  real(rtype) :: lst1,lst2,ltt1,ltt2

! *****************************************************************************

! read name of input file from command line

  fname = 'input.dat'
  
! read input parameters
  call ReadParams( fname )
  
  if( gmshin .eq. 1 )then
  call Gmsh_read
  write(*,*) 'grid converted'
  endif
  
  call InitConstants
  call PrintParams

! set no. of equations (rho, rho*u, rho*v, rho*E, ...);
! set no. of dependent variables (p, T, c, gamma, cpgas, ...)

  nconv = 4
  ndepv = 5
  if (kequs == "N") then
    ndepv = ndepv + 2   ! laminar viscosity, heat conduction coeff.
  endif

! read grid dimensions, coordinates, and triangle nodes
  call ReadGrid

  allocate( niedge(nndint),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate temporary edge pointer" )
  allocate( iedge(3,2*ntria),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate temporary edge list" )

  call EdgesInitialize( niedge,iedge )
  call EdgesFinalize( niedge,iedge )
  write(*,1000) nndint,nnodes-nndint,ntria,nedint,nedges,nbfaces,nbnodes

! allocate work space (real)
  allocate( work(2*(nconv*nnodes)),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for work space" )

! compute face vectors & cell volumes and check them;
! set coordinates of dummy nodes = those of boundary nodes;
! correct face vectors at symmetry boundaries;
! compute projections of control volumes

  call InitMetrics( niedge,iedge )
  call InitMetricsBound( niedge,iedge )
  call CheckMetrics( work )
  call FaceVectorsSymm( niedge )
  deallocate( iedge  )
  deallocate( niedge )

  call VolumeProjections
  call AllocateMemory

  allocate( iwork(nndint),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for integer works space" )

! read / initialize flow field

  if (lrest == "Y") then
    write(*,"(A,/)") " Reading initial solution ..."
    call ReadSolution
  else
    write(*,"(A,/)") " Guessing initial solution ..."
    call InitSolution
  endif

  call DependentVarsAll

  if (lrest /= "Y") call BoundaryConditions( work )

! compute limiter reference values
  call LimiterRefvals

    write(fnConv,"(A,A,A)") "convall_",Trim(title),".dat"
    open(ifConv, file=fnConv, form='formatted')
    write(fnConv,"(A,A,A)") "convrho_",Trim(title),".dat"
    open(51, file=fnConv, form='formatted')
    write(fnConv,"(A,A,A)") "convcl_",Trim(title),".dat"
    open(52,file='conv_cl.dat',form='formatted')
    write(fnConv,"(A,A,A)") "convcd_",Trim(title),".dat"
    open(53,file='conv_cd.dat',form='formatted')

    if (kflow == "E") then
    write(*,1015)
    else
    write(*,1025)
    endif

    if (lrest /= "Y") iter = 0
  
    open(99,file='timestamp.dat',form='formatted')
    call CPU_TIME(lst1) 

  do
    iter = iter + 1
    
    open(959,file='feid.dat',form='formatted')
    read(959,*) feidv
    close(959)
    
    call Solver( iwork,work )
    call Convergence
    
    if( mod(iter,outstep) == 0)then
    call CPU_TIME(ltt2) ! outstep total time
    write(99,1006) iter,ltt2
    endif
    
    if (Mod(iter,outstep) == 0)then
!    if( vtkout .eq. 1 )then
!    call output_vtk
!    endif
    call PlotSurfaces
    call WriteSolution
    endif

    if (iter>=maxiter .or. drho<=convtol) exit
    if( feidv .eq. 1 ) exit  ! force exit scenario

  enddo
  
    if( vtkout .eq. 1 )then
    call output_vtk
    endif
  
  call CPU_TIME(lst2) 
  write(99,1006) iter,lst2 ! final timestamp
  close(99)

! -----------------------------------------------------------------------------
! close file for convergence history

  close(ifConv)
  close(51);close(52);close(53)

! Post Calculation
  call output_vtk
  call PlotSurfaces
  call WriteSolution

1000 format(" No. of interior nodes: ",I8,/, &
            " No. of dummy nodes   : ",I8,/, &
            " No. of grid cells    : ",I8,/, &
            " No. of interior edges: ",I8,/, &
            " Total number of edges: ",I8,/, &
            " No. of boundary faces: ",I8,/, &
            " No. of boundary nodes: ",I8,/)
1006  format(I0,1P,2X,F0.3)
1015  format(75("-"),/, &
             " step",5X,"resid",7X,"resmax",5X,"i-res",6X,"cl", &
             10X,"cd",10X,"cm",/,75("-"))
1025  format(64("-"),/, &
             " step",5X,"resid",7X,"resmax",5X,"i-res",3X,"mass flow", &
             3X,"mass ratio",/,64("-"))

! *****************************************************************************

contains

  subroutine Usage
    write(*,"(/,A,/)") "Usage:"
    write(*,"(A,/)")   "Unstruct2D <input file>"
    stop
  end subroutine Usage

end program Unstruct2D

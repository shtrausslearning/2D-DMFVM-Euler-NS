
subroutine ReadSolution

  use ModControl
  use ModFiles
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage
  implicit none

! local variables
  integer :: errFlag, i, n, nnodesDum, nconvDum

! *****************************************************************************

  open(unit=ifRsti, file=fnRsti, status="old", action="read", &
       form="unformatted", iostat=errFlag)
  if (errFlag /= 0) call ErrorMessage( "cannot open solution file" )

! dimensions (for checking purposes)

  read(ifRsti) nnodesDum,nconvDum

  if (nnodesDum /= nnodes)  &
    call ErrorMessage( "no. of nodes differs from the grid file" )

  if (nconvDum /= nconv)  &
    call ErrorMessage( "different number of conservative variables" )

! initial residual, iteration # and solution

  read(ifRsti) drho1,iter
  read(ifRsti) ((cv(n,i), i=1,nnodes), n=1,nconv)

  close(ifRsti)

end subroutine ReadSolution

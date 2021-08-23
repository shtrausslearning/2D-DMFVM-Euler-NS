
subroutine WriteSolution

  use ModControl
  use ModFiles
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage
  implicit none

! local variables
  integer :: errFlag, i, n

! *****************************************************************************

  open(unit=ifRsto, file=fnRsto, status="unknown", action="write", &
       form="unformatted", iostat=errFlag)
  if (errFlag /= 0) call ErrorMessage( "cannot open solution file" )

! dimensions

  write(ifRsto) nnodes,nconv

! initial residual, iteration # and solution

  write(ifRsto) drho1,iter
  write(ifRsto) ((cv(n,i), i=1,nnodes), n=1,nconv)

  close(ifRsto)

end subroutine WriteSolution

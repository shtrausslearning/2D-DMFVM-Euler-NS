!> Variables and parameters related to file I/O.
!!
module ModFiles

  use ModDataTypes
  implicit none

  character(chrlen) :: fnGrid, & !< grid and topology data
                       fnFlow, & !< flow field (+ 5 digit iteration number + .v2d)
                       fnSurf, & !< quantities along wall surface(s) (+ 5 digit iteration number + .v2d)
                       fnConv, & !< convergence history (+ .v2d)
                       fnRsti, & !< restart solution - input
                       fnRsto    !< restart solution - output

  integer, parameter :: ifInp  = 10, & !< user input file (name stored in main.f90)
                        ifGrid = 20, &
                        ifFlow = 30, &
                        ifSurf = 40, &
                        ifConv = 50, &
                        ifRsti = 60, &
                        ifRsto = 70

  integer :: gmshin
  integer :: feidv
  integer :: vtkout
                        
end module ModFiles

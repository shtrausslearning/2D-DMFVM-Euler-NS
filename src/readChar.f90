!> Reads in a single-character user option from file (in ASCII format).
!!
!! @param iunit input-file unit
!! @return option as single letter
!!
function ReadChar( iunit )

  use ModDataTypes
  implicit none

! parameters
  integer, intent(in) :: iunit

! result
  character(1) :: ReadChar

! local variables
  character(chrlen) :: str
  integer :: i

! *****************************************************************************

  read(iunit,"(A)") str

  do i=1,Len_trim(str)
    ReadChar = str(i:i)
    if (ReadChar /= " ") exit  ! first non-empty char should be the option ...
  enddo

end function ReadChar


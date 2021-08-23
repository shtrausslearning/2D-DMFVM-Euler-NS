!> Prints out error message and stops program execution.
!!
!! @param message error message
!!
subroutine ErrorMessage( message )

  implicit none

! parameters
  character(*), intent(in) :: message

! *****************************************************************************

  write(*,"(A,A,/)") " Error: ",Trim( message )
  stop

end subroutine ErrorMessage


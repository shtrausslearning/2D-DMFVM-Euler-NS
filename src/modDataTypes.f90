!> Precision of all real variables; length of strings.
!!
module ModDataTypes

  implicit none

  integer, parameter :: chrlen = 256         !< length of strings
  integer, parameter :: rtype  = kind(1.D0)  !< reals with double precision

end module ModDataTypes

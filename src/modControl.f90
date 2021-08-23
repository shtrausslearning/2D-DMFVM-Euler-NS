!> Variables related to program control.
!!
module ModControl

  use ModDataTypes
  implicit none

  character(1) :: lrest !< use of previous solution for restart ("Y"=yes, "N"=no)

  integer :: maxiter, & !< max. number of iterations
             outstep, & !< number of iterations between solution dumps
             iter       !< actual iteration number

  real(rtype) :: convtol !< convergence criterion (2-norm of density change for
                         !! which the iteration process is stopped)

end module ModControl

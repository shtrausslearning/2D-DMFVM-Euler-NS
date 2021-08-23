!> Variables related to plot output.
!!
module ModPlotQuant

  use ModDataTypes
  implicit none

  integer, parameter :: mxquant =13, & !< total number of plot variables
                        mxqfield=11    !< no. of plot variables in the field (cf and Cp only at the boundaries)

  character(chrlen) :: title           !< title of the simulation case

  real(rtype) :: drho,  & !< change of the density residual (convergence criterion)
                 drho1, & !< initial change of the density residual (used for normalization)
                 cl,    & !< lift coefficient (pressure forces only; external flow)
                 cd,    & !< drag coefficient (pressure forces only; external flow)
                 cm,    & !< pitching moment coefficient wrp. to the reference point
                 mflow, & !< average mass flow rate (internal flow)
                 mfratio  !< ratio of mass flow at outlet to mass flow at inlet

end module ModPlotQuant

!> Variables related to the numerics.
!!
module ModNumerics

  use ModDataTypes
  implicit none

  character(1) :: ktimst,  & !< switch between local (="L") and global (="G") time-stepping
                  lvort,   & !< far-field vortex correction ("Y"=yes, "N"=no)
                  kprecond   !< low Mach-number preconditioning ("Y"=yes, "N"=no)

  integer :: nedges, & !< total number of edges (including edges between boundary and dummy nodes)
             nedint, & !< number of edges excluding those to dummy nodes
             iorder, & !< order of Roe's upwind scheme (1 or 2)
             nitirs, & !< number of Jacobi iterations (implicit residual smoothing)
             nrk,    & !< number of stages (Runge-Kutta scheme); max. = 5
             ldiss(5)  !< dissipation evaluation per stage (0=no, 1=yes)

  integer, allocatable :: edge(:,:) !< edge list (node i, node j)\n
    !! @details
    !! For ie > nedint, edge(*,ie) represents the edge from a boundary node
    !! (i) to a dummy node (used at inlet, outlet and far-field boundaries).

  real(rtype) :: cfl,      & !< CFL-number
                 epsirs,   & !< coefficient of implicit residual smoothing
                 limfac,   & !< limiter coefficient (Roe's upwind scheme)
                 epsentr,  & !< entropy correction coefficient (Roe's upwind scheme)
                 precoeff, & !< preconditioning parameter K (low Mach numbers)
                 ark(5),   & !< stage coefficients
                 betrk(5)    !< dissipation-blending coefficients

  real(rtype) :: volref    !< reference volume\n
                           !! @details
                           !! Parameter is required for the computation of limiter
                           !! functions (higher-order Roe scheme).
  real(rtype) :: limref(4) !< reference values of density, u, v and pressure\n
                           !! @details
                           !! Parameter is required for the computation of limiter
                           !! functions (higher-order Roe scheme).

  real(rtype), allocatable :: cvold(:,:), & !< conservative variables from previous time step
                              diss(:,:),  & !< artificial dissipation
                              rhs(:,:),   & !< residual (right-hand side)
                              lim(:,:),   & !< values of the limiter function (density, u, v, pressure)
                              tstep(:)      !< time steps (without the CFL-number)

  real(rtype), allocatable :: gradx(:,:) !< gradients of density, velocity components,
                              !! pressure and temperature with respect to the x-coordinate
  real(rtype), allocatable :: grady(:,:) !< gradients of density, velocity components,
                              !! pressure and temperature with respect to the y-coordinate

  real(rtype) :: pi, & !< 3.14...
                 rad   !< 180./pi

end module ModNumerics

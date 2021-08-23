

!> Integrates the four basic equations (continuity, momentum and energy) by
!! the explicit, multi-stage (Runge-Kutta) time-stepping scheme.
!!
!! @param iwork  integer work space for temporary variables
!! @param work   real work space for temporary variables
!!
subroutine Solver( iwork,work )

  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : BoundaryConditions, CompTheta, Cons2Prim, &
                            DependentVarsAll, DissipRoe1, DissipRoe1Prec, &
                            DissipRoe2, DissipRoe2Prec, ErrorMessage, &
                            FluxRoe1, FluxRoe2, FluxViscous, Gradients, &
                            GradientsVisc, Irsmoo, Limiter, LimiterInit, &
                            MatrixTimesInverse, Periodic, Prim2Cons, &
                            TimeStep, ZeroResiduals
  implicit none

! parameters
  integer     :: iwork(:)
  real(rtype) :: work(:)

! local variables
  integer     :: i, irk, mp, mp2
  real(rtype) :: blend1, fac, adtv, H, q2, rhop, rhoT, hT, theta, u, v
  real(rtype) :: wvec(5), wpvec(5), pmat(5,5), gmat1(5,5), dmat(5,5), r(5)
  real(rtype), allocatable :: dum1(:,:), dum2(:,:)

! *****************************************************************************
! calculate dimensions for dummy arrays (LimiterInit, Limiter, FluxRoe,
! Irsmoo, BoundaryConditions); check them

  mp  = nconv*nnodes
  mp2 = 2*mp
  if (mp2 > Ubound(work,1)) then
    call ErrorMessage( "insufficient work space in Solver" )
  endif

! store previous solution; set dissipation = 0

  cvold(:,:) = cv(:,:)
  diss(:,:)  = 0.D0

! compute the time step
  call TimeStep

! loop over the Runge-Kutta stages ============================================

  do irk=1,nrk

! - initialize dissipation

    if (irk>1 .and. ldiss(irk)/=0) then
      blend1 = 1.D0 - betrk(irk)
      do i=1,nnodes
        diss(1,i) = blend1*diss(1,i)
        diss(2,i) = blend1*diss(2,i)
        diss(3,i) = blend1*diss(3,i)
        diss(4,i) = blend1*diss(4,i)
      enddo
    endif

! - viscous flux (Navier-Stokes eqs.)

    if (ldiss(irk)/=0 .and. kequs=="N") then
      call GradientsVisc
      call FluxViscous( betrk(irk) )
    endif

    ! limiter and upwind dissipation
    if (ldiss(irk) /= 0) then
      if (iorder < 2) then
          call DissipRoe1( betrk(irk) )
      else
        dum1 = Reshape( work(1:mp)      ,(/4, nnodes/) )
        dum2 = Reshape( work((mp+1):mp2),(/4, nnodes/) )
        if (kequs == "E") call Gradients
        call LimiterInit( dum1,dum2 )
        call Limiter( dum1,dum2 )
          call DissipRoe2( betrk(irk) )
      endif
    endif

    ! convective flux; add upwind dissipation => residual
    if (iorder < 2) then
      call FluxRoe1
    else
      call FluxRoe2
    endif

    call ZeroResiduals   ! - correct residuals at symmetry/no-slip boundaries

    fac = ark(irk)*cfl
    do i=1,nndint
      adtv     = fac*tstep(i)/vol(i)
      rhs(1,i) = adtv*rhs(1,i)
      rhs(2,i) = adtv*rhs(2,i)
      rhs(3,i) = adtv*rhs(3,i)
      rhs(4,i) = adtv*rhs(4,i)
    enddo

! - implicit residual smoothing
    if (epsirs > 0.D0) then
      dum1 = Reshape( work(1:mp)      ,(/4, nnodes/) )
      dum2 = Reshape( work((mp+1):mp2),(/4, nnodes/) )
      call Irsmoo( iwork,dum1,dum2 )
      call ZeroResiduals
    endif

! - update - new solution, new dependent variables

    do i=1,nndint
      cv(1,i) = cvold(1,i) - rhs(1,i)
      cv(2,i) = cvold(2,i) - rhs(2,i)
      cv(3,i) = cvold(3,i) - rhs(3,i)
      cv(4,i) = cvold(4,i) - rhs(4,i)
    enddo

    call DependentVarsAll
    call BoundaryConditions( work )

  enddo ! irk
  
end subroutine Solver

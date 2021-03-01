
subroutine AllocateMemory

  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : ErrorMessage
  implicit none

! local variables
  integer :: errFlag

! *****************************************************************************
! base flow variables

  allocate( cv(nconv,nnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for cv()" )

  allocate( dv(ndepv,nnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for dv()" )

! general numerical variables

  allocate( cvold(nconv,nnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for cvold()" )

  allocate( diss(nconv,nnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for diss()" )

  allocate( rhs(nconv,nnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for rhs()" )

  allocate( tstep(nnodes),stat=errFlag )
  if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for tstep()" )

! numerical variables required for higher-order Roe scheme
! and for viscous flow

  if (iorder > 1) then
    allocate( lim(nconv,nnodes),stat=errFlag )
    if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for lim()" )

    if (kequs=="n" .or. kequs=="N") then
      allocate( gradx(5,nnodes),stat=errFlag )
      if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for gradx()" )
      allocate( grady(5,nnodes),stat=errFlag )
      if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for grady()" )
    else
      allocate( gradx(4,nnodes),stat=errFlag )
      if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for gradx()" )
      allocate( grady(4,nnodes),stat=errFlag )
      if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for grady()" )
    endif
  else
    if (kequs=="n" .or. kequs=="N") then
      allocate( gradx(5,nnodes),stat=errFlag )
      if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for gradx()" )
      allocate( grady(5,nnodes),stat=errFlag )
      if (errFlag /= 0) call ErrorMessage( "cannot allocate memory for grady()" )
    endif
  endif

end subroutine AllocateMemory

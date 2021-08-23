!> Computes mass flow at inlet and mass flow ratio between inlet and outlet
!! boundaries. Contributions are summed up by looping over ALL inlet and
!! outlet boundaries.
!!
subroutine Massflow

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModPlotQuant
  implicit none

! local variables
  integer     :: ibegf, iendf, n1, n2, in, out
  integer     :: ib, ibf
  real(rtype) :: sx, sy, massin, massout, mass

! *****************************************************************************
! initialize mass flow before summing up the contributions

  massin  = 0.D0
  massout = 0.D0
  mflow   = 0.D0
  mfratio = 0.D0

  in  = 0  ! flow into the domain (=1)
  out = 0  ! flow out of the domain (=1)

! loop over boundaries searching for inlet / outlet

  ibegf = 1

  do ib=1,nsegs
    iendf = ibound(1,ib)
    if (btype(ib)>=100 .and. btype(ib)<300) then
      do ibf=ibegf,iendf
        n1   = bface(1,ibf)
        n2   = bface(2,ibf)
        sx   = sbf(1,ibf)
        sy   = sbf(2,ibf)
        mass = 0.5D0*((cv(2,n1)+cv(2,n2))*sx+(cv(3,n1)+cv(3,n2))*sy)
        if (btype(ib)>=100 .and. btype(ib)<200) then
! ------- inflow
          massin  = massin - mass
          in = 1
        else
! ------- outflow
          massout = massout + mass
          out = 1
        endif
      enddo
    endif
    ibegf = iendf + 1
  enddo

! mass flow and ratio

  if (in == 1) then
    mflow   = massin
    mfratio = massout/massin
  endif
  if (in==0 .and. out==1) then
    mflow   = massout
    mfratio = 1.D0
  endif

end subroutine Massflow

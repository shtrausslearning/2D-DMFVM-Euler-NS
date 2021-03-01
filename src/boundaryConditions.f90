

!> Sets boundary conditions at dummy points. In a first loop, b.c.'s at
!! inlet, outlet and far-field are specified. In a second loop, b.c.'s are
!! set for all solid walls.
!!
!! @param work  work space for temporary variables (used by far-field b.c.)
!!
subroutine BoundaryConditions( work )

  use ModDataTypes
  use ModGeometry
  use ModPhysics
  use ModInterfaces, only : BcondFarfield, BcondInflow, BcondOutflow, &
                            BcondWallns, ErrorMessage
  implicit none

! parameters
  real(rtype) :: work(:)

! local variables
  integer :: ib, ibegn, iendn, itype, wdim

! *****************************************************************************
! loop over all boundaries with the exception of walls

  ibegn = 1

  do ib=1,nsegs

    itype = btype(ib)
    iendn = ibound(2,ib)

! - inflow

    if (itype>=100 .and. itype<200) then

      call BcondInflow( ibegn,iendn )

! - outflow

    else if (itype>=200 .and. itype<300) then

      call BcondOutflow( ibegn,iendn )

! - far-field

    else if (itype>=600 .and. itype<700) then

      wdim = Ubound(work,1)
      if ((4*nbnodes) > wdim) then
        call ErrorMessage( "insufficient work space in BoundaryConditions" )
      endif
      call BcondFarfield( ibegn,iendn, &
                          work( 1           :  nbnodes), &
                          work((1+  nbnodes):2*nbnodes), &
                          work((1+2*nbnodes):3*nbnodes), &
                          work((1+3*nbnodes):4*nbnodes) )
    endif

    ibegn = iendn + 1

  enddo ! ib

! solid walls (treated last because they should dominate)

  ibegn = 1

  do ib=1,nsegs

    itype = btype(ib)
    iendn = ibound(2,ib)

! - viscous (no-slip) wall - if Navier-Stokes equations solved

    if (itype>=300 .and. itype<400 .and. kequs=="N") then
      call BcondWallns( ibegn,iendn )
    endif

    ibegn = iendn + 1

  enddo ! ib

end subroutine BoundaryConditions

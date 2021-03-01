!> @file readParams.f90
!!
!! Reading in of user-input parameters.
!
! *****************************************************************************
!
!  (c) J. Blazek, CFD Consulting & Analysis, www.cfd-ca.de
!  Created February 25, 2014
!  Last modification: May 20, 2014
!
! *****************************************************************************
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!
! *****************************************************************************

!> Reads in user-input parameters.
!!
!! @param fname  path and name of the user input file
!!
subroutine ReadParams( fname )

  use ModControl
  use ModFiles
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage, ReadChar
  implicit none

! parameters
  character(*), intent(in) :: fname

! local variables
  character(1) :: ch
  integer :: errFlag, i

! *****************************************************************************

  open(unit=ifInp, file=fname, status="old", action="read", iostat=errFlag)
  if (errFlag /= 0) call ErrorMessage( "cannot open input file" )

  read(ifInp,"(A)") title
  read(ifInp,"(A)") fnGrid
  read(ifInp,"(A)") fnFlow
  read(ifInp,"(A)") fnSurf
  read(ifInp,"(A)") fnConv
  read(ifInp,"(A)") fnRsti
  read(ifInp,"(A)") fnRsto
  read(ifInp,*) cpcoff
  read(ifInp,*) gmshin

  ch = ReadChar( ifInp )
  if (ch=="e" .or. ch=="E") then
    kflow = "E"
  else
    kflow = "I"
  endif
  ch = ReadChar( ifInp )
  if (ch=="e" .or. ch=="E") then
    kequs = "E"
  else
    kequs = "N"
  endif
  read(ifInp,*) gamma
  read(ifInp,*) cpgas
  read(ifInp,*) renum
  read(ifInp,*) refvel
  read(ifInp,*) refrho
  read(ifInp,*) prlam

! physics - external flow
  read(ifInp,*) machinf
  read(ifInp,*) alpha
  read(ifInp,*) pinf
  read(ifInp,*) tinf

! physics - internal flow
  read(ifInp,*) ptinl
  read(ifInp,*) ttinl
  read(ifInp,*) betainl
  read(ifInp,*) pout
  read(ifInp,*) betaout
  read(ifInp,*) p12rat

! geometrical reference values
  read(ifInp,*) xref
  read(ifInp,*) yref
  read(ifInp,*) cref

! iteration control
  read(ifInp,*) maxiter
  read(ifInp,*) outstep
  read(ifInp,*) vtkout
  read(ifInp,*) convtol
  ch = ReadChar( ifInp )
  if (ch=="y" .or. ch=="Y") then
    lrest = "Y"
  else
    lrest = "N"
  endif

! numerical parameters
  read(ifInp,*) cfl
  read(ifInp,*) epsirs
  read(ifInp,*) nitirs
  ch = ReadChar( ifInp )
  if (ch=="l" .or. ch=="L") then
    ktimst = "L"
  else
    ktimst = "G"
  endif
  ch = ReadChar( ifInp )
  if (ch=="y" .or. ch=="Y") then
    kprecond = "Y"
  else
    kprecond = "N"
  endif
  read(ifInp,*) precoeff
  read(ifInp,*) iflux
  read(ifInp,*) iorder
  read(ifInp,*) limfac
  read(ifInp,*) epsentr
  ch = ReadChar( ifInp )
  if (ch=="y" .or. ch=="Y") then
    lvort = "Y"
  else
    lvort = "N"
  endif
  read(ifInp,*) nrk
  read(ifInp,*) (ark  (i), i=1,nrk)
  read(ifInp,*) (betrk(i), i=1,nrk)
  read(ifInp,*) (ldiss(i), i=1,nrk)
  close(unit=ifInp)
  

end subroutine ReadParams

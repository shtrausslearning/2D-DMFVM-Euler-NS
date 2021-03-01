!> @file modPhysics.f90
!!
!! Physical variables, settings and boundary conditions.
!! Note that all quantities are expected in SI-units.
!
! *****************************************************************************
!
!  (c) Jiri Blazek, CFD Consulting & Analysis, www.cfd-ca.de
!  Created February 25, 2014
!  Last modification: April 26, 2014
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

!> Variables related to physics and boundary conditions.
!!
module ModPhysics

  use ModDataTypes
  implicit none

  character(1) :: kequs, & !< equations solved ("E"=Euler, "N"=Navier-Stokes)
                  kflow    !< type of flow ("E"=external, "I"=internal)

! reference values

  real(rtype) :: gamma,  & !< ratio of specific heat coefficients
                 cpgas,  & !< specific heat coefficient at constant pressure
                 prlam,  & !< laminar Prandtl number
                 renum,  & !< Reynolds number
                 refvel, & !< reference velocity (internal flow only; for external flow computed from the far-field boundary)
                 refrho, & !< reference density (internal flow only; for external flow computed from the far-field boundary)
                 refvisc   !< reference dynamic viscosity coefficient (computed from renum, refvel, cref and refrho)

! boundary conditions - external flow

  real(rtype) :: machinf, & !< Mach-number at infinity
                 alpha,   & !< angle of attack
                 pinf,    & !< static pressure at infinity
                 tinf,    & !< static temperature at infinity
                 rhoinf,  & !< density at infinity
                 uinf,    & !< u-component of velocity vector at infinity
                 vinf,    & !< v-component of velocity vector at infinity
                 qinf       !< total velocity (= SQRT(uinf**2+vinf**2))

! boundary conditions - internal flow

  real(rtype) :: ptinl,   & !< total pressure at inlet
                 ttinl,   & !< total temperature at inlet
                 betainl, & !< low angle at inlet (with x-axis, positive in the clock-wise direction)
                 betaout, & !< approximate outlet angle (utilized for the initial guess only)
                 p12rat,  & !< ratio of inlet to outlet static pressure (initial guess only)
                 pout       !< static pressure at outlet

! flow variables

  integer :: nconv, & !< number of conservative variables (cv)
             ndepv    !< number of dependent variables (dv)

  real(rtype), allocatable :: cv(:,:) !< conservative variables\n
              !! @details
              !! cv(1,i) = density\n
              !! cv(2,i) = density * u\n
              !! cv(3,i) = density * v\n
              !! cv(4,i) = density * E

  real(rtype), allocatable :: dv(:,:) !< dependent variables\n
              !! @details
              !! dv(1,i) = static pressure\n
              !! dv(2,i) = static temperature\n
              !! dv(3,i) = speed of sound\n
              !! dv(4,i) = ratio of specific heats\n
              !! dv(5,i) = specific heat coefficient at constant pressure\n
              !! dv(6,i) = laminar viscosity coefficient (if viscous flow)\n
              !! dv(7,i) = laminar heat conductivity coefficient (if viscous flow)
  
  
   integer :: iflux

end module ModPhysics

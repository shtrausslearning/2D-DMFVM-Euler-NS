
subroutine PlotSurfaces

  use ModDataTypes
  use ModControl
  use ModFiles
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : ErrorMessage
  implicit none

! local variables
  character(chrlen) :: fname
  integer     :: errFlag, itype, ibegf, iendf, ibegn, iendn, ibf1, ibf2, &
                 nquant, nsurfs,ierr
  integer     :: i, ib, ibf, ibn, m
  real(rtype) :: rrho, u, v, e, press, temp, c, ptot, ttot, mach, machis, &
                 ptloss, pratio, ptotinf, gam1, ggm1
  real(rtype) :: cf, cp, visc, sx, sy, ds, sxn, syn, grdnx, grdny, grdnn, &
                 dvdnx, dvdny, dvdna, sgn

! *****************************************************************************

  ptotinf = 0.D0
  if (kflow == "E") then
    gam1    = gamma - 1.D0
    ggm1    = gamma/gam1
    ptotinf = pinf*(1.D0+0.5D0*gam1*machinf*machinf)**ggm1
  endif

! open plot file

  write(fnSurf,"(A,A,A)") 'surfp_',Trim(title),".dat"
  open(44, file=fnSurf, status="unknown", action="write", iostat=ierr)
  if (ierr /= 0) call ErrorMessage( "cannot open plot file Surface.dat" )

  do i=1,nndint
  if( ntype(i) .eq. 1 )then ! if wall node
  if( i .le. cpcoff )then

    rrho  = 1.D0/cv(1,i)
    u     = cv(2,i)*rrho
    v     = cv(3,i)*rrho
    e     = cv(4,i)*rrho
    press = dv(1,i)
    c     = dv(3,i)
    cp     = 2.0d0*(pinf-press)/(rhoinf*qinf*qinf)
    
    write(44,1006) x(i),cp

  endif
  endif
  enddo 

  close(44)

1006  format(F0.5,1X,F0.5)
1020  format(1P,20E16.8,1P,20E16.8,1P,20E16.8,1P,20E16.8)

end subroutine PlotSurfaces

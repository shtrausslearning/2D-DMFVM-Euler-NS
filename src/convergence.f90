
subroutine Convergence

  use ModDataTypes
  use ModFiles
  use ModControl
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces, only : Forces, Massflow
  implicit none

! local variables
  integer     :: i, idr
  real(rtype) :: dr, drmax

! *****************************************************************************
! compute the residual

    drho  = 0.D0
    drmax = 0.D0
    do i=1,nndint
    dr   = cv(1,i) - cvold(1,i)
    drho = drho + dr*dr
    if (Abs(dr) >= drmax) then
    drmax = Abs(dr)
    idr   = i
    endif
    enddo

    if (iter == 1) then
    drho1 = Sqrt(drho) + 1.D-32
    drho  = 1.D0
    else
    drho  = Sqrt(drho)/drho1
    endif

! compute forces & moments (external flow)

    if (kflow == "E") then
    call Forces
    else
    call Massflow
    endif

!   print out / store
    if (kflow == "E") then
    if(mod(iter,10)==1)then
    write(ifConv,1000) iter,Log10(drho),drmax,idr,cl,cd,cm
    write(51,1006) iter,log10(drho)
    write(52,1006) iter,cl
    write(53,1006) iter,cd
    write(*,1005) iter,Log10(drho),drmax,idr,cl,cd,cm
    endif
    else
    write(ifConv,1010) iter,Log10(drho),drmax,idr,mflow,mfratio
    write(*,1015) iter,Log10(drho),drmax,idr,mflow,mfratio
    endif
  
    if( isnan(Log10(drho)) )then
    write(*,*) 'drho erro'
    pause
    stop
    endif

1000  format(I5,1P,2X,E12.5,2X,E12.5,0P,I8,1P,3(2X,E12.5))
1005  format(I5,1P,2X,E11.4,2X,E10.4,0P,I8,1P,3(2X,E10.3))
1010  format(I5,1P,2X,E12.5,2X,E12.5,0P,I8,1X,1P,2(2X,E12.5))
1015  format(I5,1P,2X,E11.4,2X,E10.4,0P,I8,1X,1P,2(2X,E10.4))
1006  format(I0,1P,2X,F0.4)

end subroutine Convergence

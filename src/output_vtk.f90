
  subroutine output_vtk
! ##########################################################################################
! WRITING VTK V2.0 FILE FORMAT 
! ##########################################################################################
  use ModDataTypes
  use ModControl
  use ModFiles
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModPlotQuant
  use ModInterfaces
  implicit none

  character(chrlen) :: fname
  integer           :: ierr, nquant, i, m, CELLTYPE
  real(rtype)       :: rrho, u, v, e, press, temp, c, mach, ttot, ptot
  real(rtype)       :: ptloss, pratio, ptotinf, gam1, ggm1, visc, machis
  real(rtype)       :: varout(mxqfield+2),lpress,cp
! ######################################################################################

  ptotinf = 0.D0
  if (kflow == "E") then
    gam1    = gamma - 1.D0
    ggm1    = gamma/gam1
    ptotinf = pinf*(1.D0+0.5D0*gam1*machinf*machinf)**ggm1
  endif

  write(fname,"(A,I7.7,A)") Trim(title),iter,".vtk"
  
  open(30, file="./vtk/"//fname, status="unknown", action="write", iostat=ierr)
  if (ierr /= 0) call ErrorMessage( "cannot open plot file (field)" )

  write(30,1006)
  write(30,1007)
  write(30,1008) 
  write(30,*)""
  write(30,1009)
  write(30,1004) nndint
  
  DO i=1,nndint
    varout(1) = x(i)
    varout(2) = y(i)
    varout(3) = 0.0d0
    WRITE(30,1011) (varout(m), m=1,3)
  enddo
  
  CELLTYPE = 4*ntria
      
  write(30,*) ""
  write(30,1013) ntria,CELLTYPE
      
  do i=1,ntria
    write(30,1040) 3,tria(1,i)-1,tria(2,i)-1,tria(3,i)-1
  enddo
  
  write(30,*) ""
  write(30,1014) ntria
  
  do i=1,ntria
    write(30,1015) 5
  enddo 

   write(30,*) ''
   write(30,1012) nndint
      
    write(30,*) 'SCALARS Mach float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,nndint
     rrho  = 1.D0/cv(1,i)
     u     = cv(2,i)*rrho ; v = cv(3,i)*rrho ; c  = dv(3,i)
     mach  = SQRT(u*u+v*v)/c
     write(30,1020) mach
    enddo  
    write(30,*) 'SCALARS tempMach float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,nndint
     rrho  = 1.D0/cv(1,i)
     u     = cv(2,i)*rrho ; v = cv(3,i)*rrho ; c  = dv(3,i)
     mach  = SQRT(u*u+v*v)/c
     write(30,1020) mach
    enddo  
    write(30,*) 'SCALARS Pressure float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,nndint
    write(30,1020) dv(1,i)/pinf
    enddo
    write(30,*) 'SCALARS tempPressure float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,nndint
    write(30,1020) dv(1,i)/pinf
    enddo

    write(30,*) 'SCALARS Cp float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,nndint
     lpress = dv(1,i)
     cp         = 2.0d0*(pinf-lpress)/(rhoinf*qinf*qinf)
     write(30,1020) cp
    enddo

    close(30)

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  FORMAT(1P,20E14.6)
1040  FORMAT(1I2,3I7)
25    format('VECTORS ', a10, '   float')

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")

1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I8,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1013  FORMAT('CELLS ',I8,I8)
1014  FORMAT('CELL_TYPES ',I8)
1015  FORMAT(I1)
1012  FORMAT('POINT_DATA ',I8)

! ######################################################################################
  end subroutine output_vtk
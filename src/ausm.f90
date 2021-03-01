  
  Subroutine flux_ausmpup1
! ###############################################################################
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : FluxWalls
  implicit none

! local variables
  integer     :: i, j, ie,ni
  integer     :: ibegf,ib,ibf,iendf
  real(rtype) :: ds,nx,ny,rrho,gam1,ggm1,rl,ul,vl,pl,ql2,hl,cl2
  real(rtype) :: rr,ur,vr,pr,qr2,hr,cr2,VnLFT,VnRHT,sMLFT,sMRHT
  real(rtype) :: sMLFT2,sMRHT2,sMhat,cbar,skai,VnLFTa,VnRHTa,VnBARa
  real(rtype) :: G,VnLFTBARa,VnRHTBARa,sm,blft,brht,PT
  real(rtype) :: fd(4),sx,sy
  real(rtype) :: Mbar2,Mbar,Mo2,Mo,faMo,alphaco
  real(rtype) :: shigma,ku,kp,coff,fMpLb,fpmRa,Pu,fppLa,fMmRb,r1over2,mp,m1over2
!###############################################################################

  do ie=1,nedges
  
    i  = edge(1,ie) 
    j  = edge(2,ie)
    ds = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
    nx = sij(1,ie)/ds
    ny = sij(2,ie)/ds

    
!   AUSM CONSTANTS [ VARIABLE ]
    shigma = 1.0d0
    Ku = 0.75d0
    Kp = 0.25d0

!   left state 
    rrho = 1.D0/cv(1,i)
    gam1 = dv(4,i) - 1.0d0
    ggm1 = dv(4,i)/gam1
    
    rl   = cv(1,i)
    ul   = cv(2,i)*rrho
    vl   = cv(3,i)*rrho
    pl   = dv(1,i)
    ql2  = ul*ul + vl*vl
    hl   = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
    cl2  = dv(3,i)**2

!   rights state
    rrho = 1.D0/cv(1,j)
    gam1 = dv(4,j) - 1.0d0
    ggm1 = dv(4,j)/gam1
    
    rr   = cv(1,j)
    ur   = cv(2,j)*rrho
    vr   = cv(3,j)*rrho
    pr   = dv(1,j)
    qr2   = ur*ur + vr*vr
    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
    cr2  = dv(3,j)**2
    
!   c_bar, average speed of sound
    cbar = 0.5d0 * ( sqrt(cl2) + sqrt(cr2) ) !SI

!   %%%%%%%%%%%%%%%%%%%%%%
!   SETTING AUSM VARIABLES
!   %%%%%%%%%%%%%%%%%%%%%%
    
!   Cutoff Mach ( need to set )
    Coff = 1.5d0 ! decide

!   lhs,rhs convariant velocity 
    VnLFT = nx*ul + ny*vl !SI
    VnRHT = nx*ur + ny*vr 
!   lhs,rhs Mach states
    sMLFT = VnLFT/cbar !Ml
    sMRHT = VnRHT/cbar !Mr
    
!   added
    Mbar2 = 0.5d0/(cbar**2) * (VnLFT**2 + VnRHT**2)
    Mbar = DSQRT(Mbar2)
    Mo2 = dmin1(1.0d0,dmax1(Mbar2,(Coff)**2))
     
    Mo = DSQRT(Mo2)
    faMo = Mo*(2.0d0-Mo)
    alphaco = 0.1875d0*(-4.0d0 + 5.0d0*(faMo**2))
    
!    F +/-_p | ALPHA FUNCTION
     IF( DABS(sMLFT) >= 1.0d0 )THEN
       fppLa = 0.5d0 * ( 1.0d0 + dmax1(-1.0d0,dmin1(1.0d0,sMLFT) ))
     ELSE
      fppLa = 0.25d0*(sMLFT+1.0d0)**2*(2.0d0-sMLFT) + alphaco*sMLFT*(sMLFT**2 - 1.0d0)**2
     ENDIF 
     IF( DABS(sMRHT) >= 1.0d0 )THEN
      fpmRa = 0.5d0 * ( 1.0d0 - dmax1(-1.0d0,dmin1(1.0d0,sMRHT)) )
     ELSE
      fpmRa = 0.25d0*(sMRHT-1.0d0)**2*(2.0d0+sMRHT) - alphaco*sMRHT*(sMRHT**2 - 1.0d0)**2
     ENDIF  
    
!    F+/-_M | BETA FUNCTION
     IF( DABS(sMLFT) >= 1.0d0 )THEN
       fMpLb = 0.5d0*(sMLFT+dabs(sMLFT))
     ELSE
       fMpLb = 0.25d0*(sMLFT+1.0d0)**2 + 0.125d0*(sMLFT**2 - 1.0d0)**2
     ENDIF
     IF( DABS(sMRHT) >= 1.0d0 )THEN
       fMmRb = 0.5d0*(sMRHT-dabs(sMRHT))
     ELSE
       fMmRb = - 0.25d0*(sMRHT-1.0d0)**2 - 0.125d0*(sMRHT**2 - 1.0d0)**2
     ENDIF  
     
     R1over2 = 0.5d0*(rl+rr)
     Mp = - Kp/faMo*dmax1((1.0d0-shigma*Mbar2),0.0d0)*(pr-pl)/R1over2/(cbar**2)
     M1over2 = fMpLb + fMmRb + Mp

     Pu = - Ku * fppLa * fpmRa * (rl+rr)*(faMo*cbar)*(VnRHT-VnLFT)
    
!    PRESSURE FLUX FUNCTION
     PT = fppLa*pl + fpmRa*pr + Pu

!    MASS FLUX FUNCTION
     IF( M1over2 > 0.0d0 )THEN
      sm = M1over2 * cbar * rl
     ELSE
      sm = M1over2 * cbar * rr
     ENDIF 
     
!   AUSM+up flux
    fd(1) = (0.5d0*(sm+abs(sm))    +0.5d0*(sm-abs(sm))           )
    fd(2) = (0.5d0*(sm+abs(sm))*ul +0.5d0*(sm-abs(sm))*ur+ PT*nx )
    fd(3) = (0.5d0*(sm+abs(sm))*vl +0.5d0*(sm-abs(sm))*vr+ PT*ny )
    fd(4) = (0.5d0*(sm+abs(sm))*hl +0.5d0*(sm-abs(sm))*hr        )

!   dissipation allocation
    do ni=1,4
     diss(ni,i) = diss(ni,i) + fd(ni)*ds
     diss(ni,j) = diss(ni,j) - fd(ni)*ds
    enddo
    
  enddo
  
  do i=1,nnodes
    rhs(1,i) = diss(1,i)
    rhs(2,i) = diss(2,i)
    rhs(3,i) = diss(3,i)
    rhs(4,i) = diss(4,i)
  enddo
  
! treatment of solid walls

  call FluxWalls

!###############################################################################
  end subroutine flux_ausmpup1
  
  subroutine flux_ausmpup2
! ###############################################################################
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : FluxWalls
  implicit none

! local variables
  integer     :: i, j, ie,ni
  integer     :: ibegf,ib,ibf,iendf
  real(rtype) :: ds,nx,ny,rrho,gam1,ggm1,rl,ul,vl,pl,ql2,hl,cl2
  real(rtype) :: rr,ur,vr,pr,qr2,hr,cr2,VnLFT,VnRHT,sMLFT,sMRHT
  real(rtype) :: sMLFT2,sMRHT2,sMhat,cbar,skai,VnLFTa,VnRHTa,VnBARa
  real(rtype) :: G,VnLFTBARa,VnRHTBARa,sm,blft,brht,PT
  real(rtype) :: fd(4),rx,ry,sx,sy
  real(rtype) :: Mbar2,Mbar,Mo2,Mo,faMo,alphaco
  real(rtype) :: shigma,ku,kp,coff,fMpLb,fpmRa,Pu,fppLa,fMmRb,r1over2,mp,m1over2
!###############################################################################

  do ie=1,nedges
  
    i  = edge(1,ie)
    j  = edge(2,ie)
    ds = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
    nx = sij(1,ie)/ds
    ny = sij(2,ie)/ds
    rx = 0.5D0*(x(j)-x(i))
    ry = 0.5D0*(y(j)-y(i))
    
!   ausm constants
    shigma = 1.0d0
    Ku = 0.75d0
    Kp = 0.25d0

!   left state 
    rrho = 1.D0/cv(1,i)
    gam1 = dv(4,i) - 1.0d0
    ggm1 = dv(4,i)/gam1
    
    rl   = cv(1,i)      + lim(1,i)*(gradx(1,i)*rx+grady(1,i)*ry)
    ul   = cv(2,i)*rrho + lim(2,i)*(gradx(2,i)*rx+grady(2,i)*ry)
    vl   = cv(3,i)*rrho + lim(3,i)*(gradx(3,i)*rx+grady(3,i)*ry)
    pl   = dv(1,i)      + lim(4,i)*(gradx(4,i)*rx+grady(4,i)*ry)
    ql2  = ul*ul + vl*vl
    hl   = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
    !cl2  = dv(3,i)**2
    cl2  = gamma*pl/rl

!   rights state
    rrho = 1.D0/cv(1,j)
    gam1 = dv(4,j) - 1.0d0
    ggm1 = dv(4,j)/gam1
    
    rr   = cv(1,j)      - lim(1,j)*(gradx(1,j)*rx+grady(1,j)*ry)
    ur   = cv(2,j)*rrho - lim(2,j)*(gradx(2,j)*rx+grady(2,j)*ry)
    vr   = cv(3,j)*rrho - lim(3,j)*(gradx(3,j)*rx+grady(3,j)*ry)
    pr   = dv(1,j)      - lim(4,j)*(gradx(4,j)*rx+grady(4,j)*ry)
    qr2   = ur*ur + vr*vr
    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
    !cr2  = dv(3,j)**2
    cr2  = gamma*pr/rr
    
!   c_bar, average speed of sound
    cbar = 0.5d0 * ( sqrt(cl2) + sqrt(cr2) ) !SI

!   %%%%%%%%%%%%%%%%%%%%%%
!   SETTING AUSM VARIABLES
!   %%%%%%%%%%%%%%%%%%%%%%
    
!   Cutoff Mach ( need to set )
    Coff = 0.1d0 ! decide

!   lhs,rhs convariant velocity 
    VnLFT = nx*ul + ny*vl !SI
    VnRHT = nx*ur + ny*vr 
!   lhs,rhs Mach states
    sMLFT = VnLFT/cbar !Ml
    sMRHT = VnRHT/cbar !Mr
    !sMLFT = VnLFT/sqrt(cl2)
    !sMRHT = VnRHT/sqrt(cr2) 
    
!   added
    Mbar2 = 0.5d0/(cbar**2) * (VnLFT**2 + VnRHT**2)
    Mbar = DSQRT(Mbar2)
    Mo2 = dmin1(1.0d0,dmax1(Mbar2,(Coff)**2))
     
    Mo = DSQRT(Mo2)
    faMo = Mo*(2.0d0-Mo)
    alphaco = 0.1875d0*(-4.0d0 + 5.0d0*(faMo**2))
    
!    FOR PRESSURE TERM
     IF( DABS(sMLFT) >= 1.0d0 )THEN
       fppLa = 0.5d0 * ( 1.0d0 + dmax1(-1.0d0,dmin1(1.0d0,sMLFT) ))
     ELSE
       fppLa = 0.25d0*(sMLFT+1.0d0)**2*(2.0d0-sMLFT) + alphaco*sMLFT*(sMLFT**2 - 1.0d0)**2
       !fppLa = 0.25d0*(sMLFT+1.0d0)**2*(2.0d0-sMLFT)
     ENDIF 
     
     IF( DABS(sMRHT) >= 1.0d0 )THEN
      fpmRa = 0.5d0 * ( 1.0d0 - dmax1(-1.0d0,dmin1(1.0d0,sMRHT)) )
     ELSE
      fpmRa = 0.25d0*(sMRHT-1.0d0)**2*(2.0d0+sMRHT) - alphaco*sMRHT*(sMRHT**2 - 1.0d0)**2
      !fpmRa = 0.25d0*(sMRHT-1.0d0)**2*(2.0d0+sMRHT)
     ENDIF  
    
!    FOR INTERFACE MACH NUMBER
     IF( DABS(sMLFT) >= 1.0d0 )THEN
       fMpLb = 0.5d0*(sMLFT+dabs(sMLFT))
     ELSE
       fMpLb = 0.25d0*(sMLFT+1.0d0)**2 + 0.125d0*(sMLFT**2 - 1.0d0)**2
     ENDIF
     IF( DABS(sMRHT) >= 1.0d0 )THEN
       fMmRb = 0.5d0*(sMRHT-dabs(sMRHT))
     ELSE
       fMmRb = - 0.25d0*(sMRHT-1.0d0)**2 - 0.125d0*(sMRHT**2 - 1.0d0)**2
     ENDIF  
     
!    Terms for M1/2 !check
     R1over2 = 0.5d0*(rl+rr)
     Mp = - Kp/faMo*dmax1((1.0d0-shigma*Mbar2),0.0d0)*(pr-pl)/R1over2/(cbar**2)
     M1over2 = fMpLb + fMmRb + Mp
     !M1over2 = fMpLb + fMmRb
     
!    Terms for Pu !check
     Pu = - Ku * fppLa * fpmRa * (rl+rr)*(faMo*cbar)*(VnRHT-VnLFT)
    
!    PRESSURE FLUX FUNCTION
     PT = fppLa*pl + fpmRa*pr + Pu
     !PT = fppLa*pl + fpmRa*pr
     !PT = (pl + pr)*0.5d0
     
!    MASS FLUX FUNCTION
     IF( M1over2 > 0.0d0 )THEN
      sm = M1over2 * cbar * rl
     ELSE
      sm = M1over2 * cbar * rr
     ENDIF 
     
!   AUSM+up flux
    fd(1) = (0.5d0*(sm+abs(sm))    +0.5d0*(sm-abs(sm))           )
    fd(2) = (0.5d0*(sm+abs(sm))*ul +0.5d0*(sm-abs(sm))*ur+ PT*nx )
    fd(3) = (0.5d0*(sm+abs(sm))*vl +0.5d0*(sm-abs(sm))*vr+ PT*ny )
    fd(4) = (0.5d0*(sm+abs(sm))*hl +0.5d0*(sm-abs(sm))*hr        )

!   dissipation distribution
    do ni=1,4
     diss(ni,i) = diss(ni,i) + fd(ni)*ds
     diss(ni,j) = diss(ni,j) - fd(ni)*ds
    enddo 
    
  enddo
  
  do i=1,nnodes
    rhs(1,i) = diss(1,i)
    rhs(2,i) = diss(2,i)
    rhs(3,i) = diss(3,i)
    rhs(4,i) = diss(4,i)
  enddo
  
! treatment of solid walls
  call FluxWalls

!###############################################################################
  end subroutine flux_ausmpup2
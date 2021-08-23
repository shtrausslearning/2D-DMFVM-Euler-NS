! 
!  subroutine flux_slau1
!! ###############################################################################
!  use ModDataTypes
!  use ModGeometry
!  use ModNumerics
!  use ModPhysics
!  use ModInterfaces, only : FluxWalls
!  implicit none
!
!! local variables
!  integer     :: i, j, ie,ni
!  real(rtype) :: ds,nx,ny,rrho,gam1,ggm1,rl,ul,vl,pl,ql2,hl,cl2
!  real(rtype) :: rr,ur,vr,pr,qr2,hr,cr2,VnLFT,VnRHT,sMLFT,sMRHT
!  real(rtype) :: sMLFT2,sMRHT2,sMhat,cbar,skai,VnLFTa,VnRHTa,VnBARa
!  real(rtype) :: G,VnLFTBARa,VnRHTBARa,sm,blft,brht,PT
!  real(rtype) :: fd(4)
!!###############################################################################
!
!  do ie=1,nedges
!  
!    i  = edge(1,ie)
!    j  = edge(2,ie)
!    ds = dsqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
!    nx = sij(1,ie)/ds
!    ny = sij(2,ie)/ds
!
!!   left state 
!    rrho = 1.D0/cv(1,i)
!    gam1 = dv(4,i) - 1.0d0
!    ggm1 = dv(4,i)/gam1
!    
!    rl   = cv(1,i)
!    ul   = cv(2,i)*rrho
!    vl   = cv(3,i)*rrho
!    pl   = dv(1,i)
!    ql2  = ul*ul + vl*vl
!    hl   = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
!    cl2  = dv(3,i)**2
!
!!   rights state
!    rrho = 1.D0/cv(1,j)
!    gam1 = dv(4,j) - 1.0d0
!    ggm1 = dv(4,j)/gam1
!    
!    rr   = cv(1,j)
!    ur   = cv(2,j)*rrho
!    vr   = cv(3,j)*rrho
!    pr   = dv(1,j)
!    qr2   = ur*ur + vr*vr
!    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
!    cr2  = dv(3,j)**2
!    
!!   c_bar, average speed of sound
!    cbar = 0.5d0 * ( dsqrt(cl2) + dsqrt(cr2) ) !SI
!    
!    sMhat = dmin1(1.0d0,dsqrt(0.5d0 *( ( ul*ul + vl*vl ) + (ur*ur + vr*vr )))/cbar)
!
!!   %%%%%%%%%%%%%%%%%%%%%%
!!   SETTING SLAU VARIABLES
!!   %%%%%%%%%%%%%%%%%%%%%%
!
!!   lhs,rhs convariant velocity 
!    VnLFT = nx*ul + ny*vl !SI
!    VnRHT = nx*ur + ny*vr 
!!   lhs,rhs Mach states
!    sMLFT = VnLFT/cbar !Ml
!    sMRHT = VnRHT/cbar !Mr
!    
!    sMLFT2 = (VnLFT)**0.5d0/cbar
!    sMRHT2 = (VnRHT)**0.5d0/cbar 
!    
!    skai = (1.0d0-sMhat)**2.0d0
!    
!    VnLFTa = abs(VnLFT) 
!    VnRHTa = abs(VnRHT)
!!   Vn_bar
!    VnBARa = ( rl*abs(VnLFT) + rr*abs(VnRHT) )/(rl+rr)
!    
!    G = -max(min(sMLFT,0.0d0),-1.0d0)*min(max(VnRHT,0.0d0),1.0d0)
!    G = min(max(0.0d0,G),1.0d0)
!    
!    VnLFTBARa = (1.0d0-G)*VnBARa + G*VnLFTa
!    VnRHTBARa = (1.0d0-G)*VnBARa + G*VnRHTa
!    
!!   MASS FLUX FUNCTION
!    sm = 0.5d0*( rl*(VnLFT + VnLFTBARa) + rr*(VnRHT - VnRHTBARa) - skai*(pl - pr)/cbar )
!    
!!   P+/P- w/ alpha = 0
!    IF( abs( sMLFT ) >= 1.0d0 )THEN
!     blft = 0.5d0 * ( 1.0d0 + max(-1.0d0,min(1.0d0,sMLFT) ))
!    ELSE
!     blft = 0.25d0 * ((sMLFT + 1.0d0)**2) * (2.0d0 - sMLFT)  ! P+|a=0
!    ENDIF  
!      
!    IF( abs( sMRHT ) >= 1.0d0 )THEN
!     brht = 0.5d0 * ( 1.0d0 - max(-1.0d0,min(1.0d0,sMRHT)) )
!    ELSE
!     brht = 0.25d0 * ((sMRHT - 1.0d0)**2) * (2.0d0 + sMRHT)  ! P-|a=0
!    ENDIF    
!    
!!   SLAU1 PRESSURE FLUX FUNCTION
!    !PT = 0.5d0*(pl  + pr) &
!    ! + 0.5d0*(blft  - brht)*(pl - pr) &
!    ! + 0.5d0*(1.0d0 - skai)*(blft + brht - 1.0d0)*(pl + pr)
!   !SLAU2 PRESSURE FLUX FUNCTION
!    PT = 0.5d0*(pl  + pr)          &                             
!      + 0.5d0*(blft  - brht)*(pl - pr)   &                      
!      + SQRT(0.5d0*(ql2+qr2))*(blft+brht-1.0d0)*(0.5d0*(rl+rr))*cbar
!     
!!   SLAU flux
!    fd(1) = (0.5d0*(sm+abs(sm))    +0.5d0*(sm-abs(sm))           )
!    fd(2) = (0.5d0*(sm+abs(sm))*ul +0.5d0*(sm-abs(sm))*ur+ PT*nx )
!    fd(3) = (0.5d0*(sm+abs(sm))*vl +0.5d0*(sm-abs(sm))*vr+ PT*ny )
!    fd(4) = (0.5d0*(sm+abs(sm))*hl +0.5d0*(sm-abs(sm))*hr        )
!
!! - diss replaces
!    do ni=1,4
!     diss(ni,i) = diss(ni,i) + fd(ni)*ds
!     diss(ni,j) = diss(ni,j) - fd(ni)*ds
!    enddo
!    
!  enddo
!  
!  do i=1,nnodes
!    rhs(1,i) = diss(1,i)
!    rhs(2,i) = diss(2,i)
!    rhs(3,i) = diss(3,i)
!    rhs(4,i) = diss(4,i)
!  enddo
!  
!! treatment of solid walls
!  call FluxWalls
!
!!###############################################################################
!  end subroutine flux_slau1
!
!  subroutine flux_slau2
!!###############################################################################
!  use ModDataTypes
!  use ModGeometry
!  use ModNumerics
!  use ModPhysics
!  use ModInterfaces, only : FluxWalls
!  
!! local variables
!  integer     :: i, j, ie,ni
!  real(rtype) :: ds,nx,ny,rrho,gam1,ggm1,rl,ul,vl,pl,ql2,hl,cl2
!  real(rtype) :: rr,ur,vr,pr,qr2,hr,cr2,VnLFT,VnRHT,sMLFT,sMRHT
!  real(rtype) :: sMLFT2,sMRHT2,sMhat,cbar,skai,VnLFTa,VnRHTa,VnBARa
!  real(rtype) :: G,VnLFTBARa,VnRHTBARa,sm,blft,brht,PT
!  real(rtype) :: fd(4),rx,ry
!!###############################################################################
!
!  do ie=1,nedges
!  
!    i  = edge(1,ie)
!    j  = edge(2,ie)
!    ds = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
!    nx = sij(1,ie)/ds
!    ny = sij(2,ie)/ds
!    rx = 0.5D0*(x(j)-x(i))
!    ry = 0.5D0*(y(j)-y(i))
!
!! - left & right state
!
!    rrho = 1.D0/cv(1,i)
!    gam1 = dv(4,i) - 1.D0
!    ggm1 = dv(4,i)/gam1
!    
!    rl   = cv(1,i)      + lim(1,i)*(gradx(1,i)*rx+grady(1,i)*ry)
!    ul   = cv(2,i)*rrho + lim(2,i)*(gradx(2,i)*rx+grady(2,i)*ry)
!    vl   = cv(3,i)*rrho + lim(3,i)*(gradx(3,i)*rx+grady(3,i)*ry)
!    pl   = dv(1,i)      + lim(4,i)*(gradx(4,i)*rx+grady(4,i)*ry)
!    ql2  = ul*ul + vl*vl
!    hl   = ggm1*pl/rl + 0.5D0*(ul*ul+vl*vl)
!    !cl2  = dv(3,i)**2
!    cl2  = gamma*pl/rl
!
!    rrho = 1.D0/cv(1,j)
!    gam1 = dv(4,j) - 1.D0
!    ggm1 = dv(4,j)/gam1
!    
!    rr   = cv(1,j)      - lim(1,j)*(gradx(1,j)*rx+grady(1,j)*ry)
!    ur   = cv(2,j)*rrho - lim(2,j)*(gradx(2,j)*rx+grady(2,j)*ry)
!    vr   = cv(3,j)*rrho - lim(3,j)*(gradx(3,j)*rx+grady(3,j)*ry)
!    pr   = dv(1,j)      - lim(4,j)*(gradx(4,j)*rx+grady(4,j)*ry)
!    qr2   = ur*ur + vr*vr
!    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
!    !cr2  = dv(3,j)**2
!    cr2  = gamma*pr/rr
!
!!   c_bar, average speed of sound
!    cbar = 0.5d0 * ( sqrt(cl2) + sqrt(cr2) ) !SI
!    
!    sMhat = dmin1(1.0d0,sqrt(0.5d0 *( ( ul*ul + vl*vl ) + (ur*ur + vr*vr )))/cbar)
!
!!   %%%%%%%%%%%%%%%%%%%%%%
!!   SETTING SLAU VARIABLES
!!   %%%%%%%%%%%%%%%%%%%%%%
!
!!   lhs,rhs convariant velocity 
!    VnLFT = nx*ul + ny*vl !SI
!    VnRHT = nx*ur + ny*vr 
!!   lhs,rhs Mach states
!    sMLFT = VnLFT/cbar !Ml
!    sMRHT = VnRHT/cbar !Mr
!    
!    sMLFT2 = (VnLFT)**0.5d0/cbar
!    sMRHT2 = (VnRHT)**0.5d0/cbar 
!    
!    skai = (1.0d0-sMhat)**2.0d0
!    
!    VnLFTa = abs(VnLFT) 
!    VnRHTa = abs(VnRHT)
!!   Vn_bar
!    VnBARa = ( rl*abs(VnLFT) + rr*abs(VnRHT) )/(rl+rr)
!    
!    G = -max(min(sMLFT,0.0d0),-1.0d0)*min(max(VnRHT,0.0d0),1.0d0)
!    G = min(max(0.0d0,G),1.0d0)
!    
!    VnLFTBARa = (1.0d0-G)*VnBARa + G*VnLFTa
!    VnRHTBARa = (1.0d0-G)*VnBARa + G*VnRHTa
!    
!!   MASS FLUX FUNCTION
!    sm = 0.5d0*( rl*(VnLFT + VnLFTBARa) + rr*(VnRHT - VnRHTBARa) - skai*(pl - pr)/cbar )
!    
!!   P+/P- w/ alpha = 0
!    IF( abs( sMLFT ) >= 1.0d0 )THEN
!     blft = 0.5d0 * ( 1.0d0 + max(-1.0d0,min(1.0d0,sMLFT) ))
!    ELSE
!     blft = 0.25d0 * ((sMLFT + 1.0d0)**2) * (2.0d0 - sMLFT)  ! P+|a=0
!    ENDIF  
!      
!    IF( abs( sMRHT ) >= 1.0d0 )THEN
!     brht = 0.5d0 * ( 1.0d0 - max(-1.0d0,min(1.0d0,sMRHT)) )
!    ELSE
!     brht = 0.25d0 * ((sMRHT - 1.0d0)**2) * (2.0d0 + sMRHT)  ! P-|a=0
!    ENDIF    
!    
!!!   SLAU1 PRESSURE FLUX FUNCTION
!!    PT = 0.5d0*(pl  + pr) &
!!     + 0.5d0*(blft  - brht)*(pl - pr) &
!!     + 0.5d0*(1.0d0 - skai)*(blft + brht - 1.0d0)*(pl + pr)
!  ! SLAU2 PRESSURE FLUX FUNCTION
!   PT = 0.5d0*(pl+pr)+0.5d0*(blft-brht)*(pl-pr) + SQRT(0.5d0*(ql2+qr2))*(blft+brht-1.0d0)*(0.5d0*(rl+rr))*cbar
!     
!!   SLAU flux
!    fd(1) = (0.5d0*(sm+abs(sm))    +0.5d0*(sm-abs(sm))           )
!    fd(2) = (0.5d0*(sm+abs(sm))*ul +0.5d0*(sm-abs(sm))*ur+ PT*nx )
!    fd(3) = (0.5d0*(sm+abs(sm))*vl +0.5d0*(sm-abs(sm))*vr+ PT*ny )
!    fd(4) = (0.5d0*(sm+abs(sm))*hl +0.5d0*(sm-abs(sm))*hr        )
!
!    do ni=1,4
!     diss(ni,i) = diss(ni,i) + fd(ni)*ds
!     diss(ni,j) = diss(ni,j) - fd(ni)*ds
!    enddo    
!    
!  enddo
!  
!  do i=1,allnod
!    rhs(1,i) = diss(1,i)
!    rhs(2,i) = diss(2,i)
!    rhs(3,i) = diss(3,i)
!    rhs(4,i) = diss(4,i)
!  enddo
!  
!! treatment of solid walls
!  call FluxWalls
!
!!###############################################################################
!  end subroutine flux_slau2
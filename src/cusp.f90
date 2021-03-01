
  subroutine cuspdissip1
!###############################################################################
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : FluxWalls
  implicit none

! local variables
  integer     :: i, j, ie,ni
  real(rtype) :: ds, nx, ny, gam1, rrho, rl, ul, vl, pl, hl
  real(rtype) :: rr, ur, vr, pr, hr, rav, dd, dd1, uav, vav, hav, q2a
  real(rtype) :: c2a, cav, uv, du, h1, h2, h3, h4, h5, encordelta, eabs1
  real(rtype) :: eabs2, eabs4
  real(rtype) :: fd(4)
  real(rtype) :: qsl,ccl,mml,rhl,ccr,qsr,mmr,rhr,mn1p,vn1p,eigp,eigm,termbeta1
  real(rtype) :: betaf, alphaf,c1p2
!###############################################################################

  do ie=1,nedges
  
    i  = edge(1,ie)
    j  = edge(2,ie)
    ds = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
    nx = sij(1,ie)/ds
    ny = sij(2,ie)/ds

!   left state
    rrho = 1.D0/cv(1,i)
    rl   = cv(1,i)
    ul   = cv(2,i)*rrho
    vl   = cv(3,i)*rrho
    pl   = dv(1,i)
    hl   = (pl+cv(4,i))*rrho
    ccl  = dv(3,i) ! C_l
    qsl  = ( cv(2,i)*nx + cv(3,i)*ny )/cv(1,i) ! V_l
    mml  = qsl/ccl ! M_l
    rhl = dv(1,i) + cv(4,i)
    
!   right state
    rrho = 1.D0/cv(1,j)
    rr   = cv(1,j)
    ur   = cv(2,j)*rrho
    vr   = cv(3,j)*rrho
    pr   = dv(1,j)
    hr   = (pr+cv(4,j))*rrho
    ccr  = dv(3,j)
    qsr  = ( cv(2,j)*nx + cv(3,j)*ny )/cv(1,j)
    mmr  = qsr/ccr 
    rhr = dv(1,j) + cv(4,j)

!   Interface Mach, convariant velocity
    c1p2 = 0.5d0*( ccl + ccr )
    vn1p = 0.5d0*( (ul+ur)*nx + (vl+vr)*ny )
    mn1p = vn1p/c1p2
    !mn1p = 0.5d0*( mml + mmr )

!   roe averaging
    rav   = Sqrt(rl*rr)
    gam1  = 0.5D0*(dv(4,i)+dv(4,j)) - 1.D0
    dd    = rav/rl
    dd1   = 1.D0/(1.D0+dd)
    uav   = (ul+dd*ur)*dd1
    vav   = (vl+dd*vr)*dd1
    hav   = (hl+dd*hr)*dd1
    q2a   = 0.5D0*(uav*uav+vav*vav)
    c2a   = gam1*(hav-q2a)
    cav   = Sqrt(c2a)
    uv    = uav*nx + vav*ny           ! convariant Roe 

!   eigenvalues
    eigp = 0.5d0*uv*( gamma + 1.0d0 )/gamma + sqrt( (0.5d0*uv*(gamma - 1.0d0)/gamma)**2 + c2a/gamma )
    eigm = 0.5d0*uv*( gamma + 1.0d0 )/gamma - sqrt( (0.5d0*uv*(gamma - 1.0d0)/gamma)**2 + c2a/gamma )

!   factor beta
    if( mn1p < 1.0d0 .and. mn1p >= 0.0d0 )then
     termbeta1 = (vn1p + eigm)/(vn1p - eigm)
     betaf = max1( 0.0d0, termbeta1 ) 
    elseif(  mn1p > -1.0d0 .and. mn1p < 0.0d0 )then
     termbeta1  = (vn1p + eigp)/(vn1p - eigp)
     betaf = - max1( 0.0d0, termbeta1 )
    elseif( abs(mn1p) >= 1.0d0 )then
     betaf = sign(1.0d0,mn1p) 
    endif 

!   factor alpha
    if( betaf == 0.0d0 )then
     alphaf = abs(vn1p)
    elseif( betaf > 0.0d0 .and. mn1p > 0.0d0 .and. mn1p < 1.0d0 )then
     alphaf = - ( 1.0d0 + betaf )*eigm
    elseif( betaf < 0.0d0 .and. mn1p > -1.0d0 .and. mn1p < 0.0d0 )then
     alphaf = ( 1.0d0 - betaf )*eigp
    elseif( abs( mn1p ) > 1.0d0 )then
     alphaf = 0.0d0
    endif 
    
!   CUSP dissipation term
    fd(1) = alphaf*(cv(1,j)-cv(1,i)) + betaf*( qsr*cv(1,j) - qsl*cv(1,i) )
    fd(2) = alphaf*(cv(2,j)-cv(2,i)) + betaf*( (qsr*cv(2,j)+nx*pr) - (qsl*cv(2,i)+nx*pl) )
    fd(3) = alphaf*(cv(3,j)-cv(3,i)) + betaf*( (qsr*cv(3,j)+ny*pr) - (qsl*cv(3,i)+ny*pl) )
    fd(4) = alphaf*(  rhr - rhl    ) + betaf*( qsr*rhr - qsl*rhl )
    
! - edge contributions to dissipation
    do ni = 1,4
     diss(ni,i) = diss(ni,i) + 0.5d0*fd(ni)*ds
     diss(ni,j) = diss(ni,j) - 0.5d0*fd(ni)*ds
    enddo
    
  enddo
  
!###############################################################################
  end subroutine cuspdissip1
  
  subroutine cuspdissip2
!###############################################################################
  use ModDataTypes
  use ModGeometry
  use ModNumerics
  use ModPhysics
  use ModInterfaces, only : FluxWalls
  implicit none

! local variables
  integer     :: i, j, ie,ni
  real(rtype) :: ds, nx, ny, gam1, rrho, rl, ul, vl, pl, hl
  real(rtype) :: rr, ur, vr, pr, hr, rav, dd, dd1, uav, vav, hav, q2a
  real(rtype) :: c2a, cav, uv, du, h1, h2, h3, h4, h5, encordelta, eabs1
  real(rtype) :: eabs2, eabs4
  real(rtype) :: fd(4)
  real(rtype) :: qsl,ccl,mml,rhl,ccr,qsr,mmr,rhr,mn1p,vn1p,eigp,eigm,termbeta1
  real(rtype) :: betaf, alphaf,rx,ry,ggm1,c1p2
!###############################################################################

  do ie=1,nedges
  
    i  = edge(1,ie)
    j  = edge(2,ie)
    ds = Sqrt(sij(1,ie)*sij(1,ie)+sij(2,ie)*sij(2,ie))
    nx = sij(1,ie)/ds
    ny = sij(2,ie)/ds
    rx = 0.5D0*(x(j)-x(i))
    ry = 0.5D0*(y(j)-y(i))

!   left state
    rrho = 1.D0/cv(1,i)
    gam1 = dv(4,i) - 1.D0
    ggm1 = dv(4,i)/gam1
    rl   = cv(1,i)      + lim(1,i)*(gradx(1,i)*rx+grady(1,i)*ry)
    ul   = cv(2,i)*rrho + lim(2,i)*(gradx(2,i)*rx+grady(2,i)*ry)
    vl   = cv(3,i)*rrho + lim(3,i)*(gradx(3,i)*rx+grady(3,i)*ry)
    pl   = dv(1,i)      + lim(4,i)*(gradx(4,i)*rx+grady(4,i)*ry)
    hl   = ggm1*pl/rl + 0.5D0*(ul*ul+vl*vl)
    ccl  = gamma*pl/rl
    qsl  = ul*nx + vl*ny
    mml  = qsl/ccl ! M_l
    
!   right state
    rrho = 1.D0/cv(1,j)
    gam1 = dv(4,j) - 1.D0
    ggm1 = dv(4,j)/gam1
    rr   = cv(1,j)      - lim(1,j)*(gradx(1,j)*rx+grady(1,j)*ry)
    ur   = cv(2,j)*rrho - lim(2,j)*(gradx(2,j)*rx+grady(2,j)*ry)
    vr   = cv(3,j)*rrho - lim(3,j)*(gradx(3,j)*rx+grady(3,j)*ry)
    pr   = dv(1,j)      - lim(4,j)*(gradx(4,j)*rx+grady(4,j)*ry)
    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
    ccr  = gamma*pr/rr
    qsr  = ur*nx + vr*ny
    mmr  = qsr/ccr 

!   Interface Mach, convariant velocity
    c1p2 = 0.5d0*( ccl + ccr )
    vn1p = 0.5d0*( (ul+ur)*nx + (vl+vr)*ny )
    mn1p = vn1p/c1p2
    !mn1p = 0.5d0*( mml + mmr )

!   roe averaging
    rav   = Sqrt(rl*rr)
    gam1  = 0.5D0*(dv(4,i)+dv(4,j)) - 1.D0
    dd    = rav/rl
    dd1   = 1.D0/(1.D0+dd)
    uav   = (ul+dd*ur)*dd1
    vav   = (vl+dd*vr)*dd1
    hav   = (hl+dd*hr)*dd1
    q2a   = 0.5D0*(uav*uav+vav*vav)
    c2a   = gam1*(hav-q2a)
    cav   = Sqrt(c2a)
    uv    = uav*nx + vav*ny           ! convariant Roe 

!   eigenvalues
    eigp = 0.5d0*uv*( gamma + 1.0d0 )/gamma + sqrt( (0.5d0*uv*(gamma - 1.0d0)/gamma)**2 + c2a/gamma )
    eigm = 0.5d0*uv*( gamma + 1.0d0 )/gamma - sqrt( (0.5d0*uv*(gamma - 1.0d0)/gamma)**2 + c2a/gamma )

!   factor beta
    if( mn1p < 1.0d0 .and. mn1p >= 0.0d0 )then
     termbeta1 = (vn1p + eigm)/(vn1p - eigm)
     betaf = max1( 0.0d0, termbeta1 ) 
    elseif(  mn1p > -1.0d0 .and. mn1p < 0.0d0 )then
     termbeta1  = (vn1p + eigp)/(vn1p - eigp)
     betaf = - max1( 0.0d0, termbeta1 )
    elseif( abs(mn1p) >= 1.0d0 )then
     betaf = sign(1.0d0,mn1p) 
    endif 

!   factor alpha
    if( betaf == 0.0d0 )then
     alphaf = abs(vn1p)
    elseif( betaf > 0.0d0 .and. mn1p > 0.0d0 .and. mn1p < 1.0d0 )then
     alphaf = - ( 1.0d0 + betaf )*eigm
    elseif( betaf < 0.0d0 .and. mn1p > -1.0d0 .and. mn1p < 0.0d0 )then
     alphaf = ( 1.0d0 - betaf )*eigp
    elseif( abs( mn1p ) > 1.0d0 )then
     alphaf = 0.0d0
    endif 
    
!   CUSP dissipation term
    fd(1) = alphaf*(rr    - rl   ) + betaf*(  qsr*rr           -  qsl*rl           )
    fd(2) = alphaf*(rr*ur - rl*ul) + betaf*( (qsr*rr*ur+nx*pr) - (qsl*rl*ul+nx*pl) )
    fd(3) = alphaf*(rr*vr - rl*vl) + betaf*( (qsr*rr*vr+ny*pr) - (qsl*rl*vl+ny*pl) )
    fd(4) = alphaf*(rr*hr - rl*hl) + betaf*(  qsr*rr*hr        -  qsl*rl*hl        )
    
! - edge contributions to dissipation
    do ni = 1,4
     diss(ni,i) = diss(ni,i) + 0.5d0*fd(ni)*ds
     diss(ni,j) = diss(ni,j) - 0.5d0*fd(ni)*ds
    enddo

    
  enddo
  
!###############################################################################
  end subroutine cuspdissip2
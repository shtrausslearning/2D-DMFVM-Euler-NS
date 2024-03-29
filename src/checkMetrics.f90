!> Checks metrics by computing min. and max. volume, and by computing
!! the sum of face vectors for each control volume (should be zero).
!!
!! @param work  work space (used to store the sum of face vectors
!!              of a control volume - fvecSum)
!!
subroutine CheckMetrics( work )

  use ModDataTypes
  use ModGeometry
  use ModNumerics
  implicit none

! parameters
  real(rtype) :: work(:)

! local variables
  integer     :: i, j, ib, ie, ibegf, iendf, ibegn, iendn
  real(rtype) :: volmin, volmax, s, smax
  real(rtype), allocatable :: fvecSum(:,:)

! *****************************************************************************

  fvecSum = Reshape( work,(/2, nndint/) )

! obtain mix. and max. control volume

  volmin = Minval( vol )
  volmax = Maxval( vol )

! compute sum of face vectors for each control volume (fvecSum)

  do i=1,nndint
    fvecSum(1,i) = 0.D0
    fvecSum(2,i) = 0.D0
  enddo

  do ie=1,nedint
    i            = edge(1,ie)
    j            = edge(2,ie)
    fvecSum(1,i) = fvecSum(1,i) + sij(1,ie)
    fvecSum(2,i) = fvecSum(2,i) + sij(2,ie)
    fvecSum(1,j) = fvecSum(1,j) - sij(1,ie)
    fvecSum(2,j) = fvecSum(2,j) - sij(2,ie)
  enddo

  ibegf = 1
  do ib=1,nsegs
    iendf = ibound(1,ib)
    if (btype(ib)<700 .or. btype(ib)>=800) then   ! boundary faces (non-periodic)
      do ie=ibegf,iendf
        i            = bface(1,ie)
        j            = bface(2,ie)
        fvecSum(1,i) = fvecSum(1,i) + 0.5D0*sbf(1,ie)
        fvecSum(2,i) = fvecSum(2,i) + 0.5D0*sbf(2,ie)
        fvecSum(1,j) = fvecSum(1,j) + 0.5D0*sbf(1,ie)
        fvecSum(2,j) = fvecSum(2,j) + 0.5D0*sbf(2,ie)
      enddo
    endif
    ibegf = iendf + 1
  enddo

  ibegn = 1
  do ib=1,nsegs
    iendn = ibound(2,ib)
    if (btype(ib)>=700 .and. btype(ib)<800) then   ! periodic nodes
      do ie=ibegn,iendn
        i            = bnode(1,ie)
        j            = bnode(2,ie)
        fvecSum(1,i) = fvecSum(1,i) + fvecSum(1,j)
        fvecSum(2,i) = fvecSum(2,i) + fvecSum(2,j)
        fvecSum(1,j) = fvecSum(1,i)
        fvecSum(2,j) = fvecSum(2,i)
      enddo
    endif
    ibegn = iendn + 1
  enddo

! compute maximum of the sum of face vectors

  smax = -1.D+32
  do i=1,nndint
    s    = Sqrt(fvecSum(1,i)**2+fvecSum(2,i)**2)
    smax = Max(smax,s)
  enddo

! print info

  write(*,1000) smax,volmin,volmax

1000  format(" max. sum(S) = ",E11.4,/," min. volume = ",E11.4,/, &
             " max. volume = ",E11.4,/)

end subroutine CheckMetrics

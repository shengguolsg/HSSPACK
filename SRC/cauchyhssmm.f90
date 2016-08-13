  subroutine CauchyHssMM( D,F,U,V,N,ni,tol,ncol,B,LDB,time )
    use aux_hss
    use BasicMM
    use CauchyHss_VP2
    use ConstructHssd
! 
    implicit none

! Scalar parameters
    integer :: N,ni,ncol,LDB
    double precision :: tol
! Array parameters
    double precision :: D(*),F(*),U(*),V(*),B(LDB,*)
    double precision :: time(2)
!
! Purpose
! ========
! This routine computes the multiplication of Cauchy-like matrix A with a general
! dense matrix, B. The generators of Cauchy-like matrix are D,F,U and V, 
! A(i,j) = U(i)*V(j)/(D(i)-F(j)). 
!
! D   (input/output) double array, dimension(N), the row generator of A;
!     On exit, destroyed. 
!
! F   (input/output) double array, dimension(N), the col generator of A;
!     On exit, destroyed. 
!
! U   (input/output) double array, dimension(N), the row generator of A;
!     On exit, destroyed. 
!
! V   (input/output) double array, dimension(N), the col generator of A;
!     On exit, destroyed. 
! 
! N   (input) integer, the dimension of A;
!
! ni  (input) integer, the row partition of A;
!
! tol (input) double precisionn, the accuracy parameter for low-rank approximation;
!
! B   (input/output) double precision, dimension(LDB,*)
!     Store the output matrix B = A*B
!
! LDB (input) integer, the leading dimension of B
!
! time (output) double precision, timing information
!      time(1) : HSS construction time
!      time(2) : HSS matrix-matrix multiplication time
!
! ========
! Written by S.G. Li, on March 15th, 2014
! ======================================
!
    integer  :: ltr, lm, n1, nlast, ierr, n0, pnh, tnk, nc, lvl, i
    double precision :: time2, s1
    double precision, allocatable :: WW(:,:), SX(:), H(:), DD(:)
    
    ! HSS tree
    integer, allocatable :: tr(:), m(:)
    type(HTinfo)         :: TRE
    type(hssmm)          :: PH

    ! execution lines
    allocate (m(n/ni+1), tr( 2*(n/ni+1)-1), stat=ierr )
    IF(ierr /= 0 ) Then
       write(*,*) "Allocate failed in cauchy2hss! "
       stop
    END IF
    m   = 0
    ltr = 0
    call npart(n, ni, tr, ltr, m, lm)

!   Tree information
    lvl = ceiling( log2(ltr+1) ) -1
    n1 = 0.5* (ltr +1 )

    allocate(TRE%ttr(lvl,n1), TRE%lentr(lvl), TRE%lfx(ltr), TRE%ch(2,ltr) )
    call GetTreeInfo(TR, LTR, M, M, TRE, lvl)       

!   ******************** HSS Matrix Approximation *********************
!   Allocate space and initialize variables
    tnk = 0
    DO i = 1, lm
       tnk = tnk+ m(i)**2
    END DO

    allocate( PH%D(ltr), PH%U(2,ltr), PH%V(2,ltr), PH%B(2,ltr),PH%ph(4,ltr), PH%pd(2,ltr) )
    allocate( H(8*ni*n), DD(tnk) ) ! test only; can be optimized
    call hssexpmm0(PH, ltr)

!   HSS matrix
    pnh = 1
    N0 = N
    time(1) = comp_time()
!    call cauchy2hss_OMP( D,F,U,V,N0,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh,TRE,info )
    call cauchy2hss( D,F,U,V,N0,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh )
    time2= comp_time()
    time(1) = time2-time(1)
    s1 = (1.0D0*pnh )/(ni*n)
!   print*, 'construction time: ',time(1), s1

!   **************** HSS Matrix-Matrix Multiplication *******************
    Allocate(WW(2*n,ncol), SX(2*lvl*n*ncol))

    nc = floor( dble(N0) /ncol )
    time(2) = comp_time()
    Do i = 1, nc
       call fastHssmm(DD,H,PH,m,B(1:N,(i-1)*ncol+1:i*ncol),tr,ltr,n,ncol,TRE,lvl,WW,SX)
    End Do
    nlast = n - nc*ncol
    IF(nlast .ne. 0) THEN
       i = nc + 1
       call fastHssmm(DD,H,PH,m,B(1:N,(i-1)*ncol+1:n),tr,ltr,n,nlast,TRE,lvl,WW,SX)
    END IF
    time2 = comp_time()
    time(2) = time2 - time(2)
    print*, 'multiplication time: ', time(2), 'ncol ', ncol


    deallocate(TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch, m, tr)
    deallocate(PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd )
    deallocate( H,DD, WW, SX)

  end subroutine CauchyHssMM

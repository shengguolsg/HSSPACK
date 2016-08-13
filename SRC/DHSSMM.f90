SUBROUTINE DHSSMM(A,LDA,M,N,B,LDB,K,tol,prank,ni,ncol,time,info)
      use aux_hss
      use BasicMM
      use ConstructHssd
!
      implicit none
!
! .. Scalar Arguments
      integer  :: LDA, M, N, K, prank, ni, ncol, LDB, info
      double precision :: tol 
!
! .. Array Arguments
      double precision :: A(LDA,*), B(LDB,*),time(2)
!
! Purpose
! =======
!
! This routine assume A is an HSS matrix and computes the multiplication of A and B, 
! where B is a general dense, LDB-by-K, matrix. This HSS construction algorithm use
! SVD or randomized SVD to compute the low-rank approximation. 
!
! A    (input/output) double array, dimension(LDA,*) 
!      A is a dense HSS matrix. For simplicity, we assume LDA=M=N. 
!
! LDA  (input) integer, the leading dimension of A
!
! M    (input) integer, the number of rows of A to multiplicate
!
! N    (input) integer, the number of cols of A to multiplicate
!
! B    (input/output), dimension(LDB,*)
!      B is a general dense matrix, and for simplicity, we assume LDB=N=K. 
!
! LDB  (input) integer, the leading dimension of B
!
! K    (input) integer, the number of cols of B to multiplicate
!
! tol  (input) double precision, accuracy parameter for low-rank approximation
!
! prank (input) integer, the prefixed rank parameter
! 
! ni   (input) integer, the size of block row partition
!
! ncol (input) integer, the size of right-hand
!
! time (output) double precision, the timing information
!
! info (output) integer, the information of execution
!
! ==========
! Written by S.G. Li, on March 17th, 2014
! ======================================================
!
    integer  :: ltr, lm, n1, nlast, ierr, pnh, tnk, nc, lvl,i
    double precision :: time2, time1
    double precision, allocatable :: WW(:,:),SX(:),H(:),DD(:),WORK(:)
    
    ! HSS tree
    integer, allocatable :: tr(:), blk(:)
    type(HTinfo)         :: TRE
    type(hssmm)          :: PH

    ! execution lines
    allocate( blk(n/ni+1), tr( 2*(n/ni+1)-1), stat=ierr )
    IF(ierr /= 0 ) Then
       write(*,*) "Allocate failed in cauchy2hss! "
       stop
    END IF
    blk   = 0
    ltr = 0
    call npart(n, ni, tr, ltr, blk, lm)

!   Tree information
    lvl = ceiling( log2(ltr+1) ) -1
    n1 = 0.5* (ltr +1 )

    allocate( TRE%ttr(lvl,n1), TRE%lentr(lvl), TRE%lfx(ltr), TRE%ch(2,ltr) )
    call GetTreeInfo(TR, LTR, BLK, BLK, TRE, lvl)       

!   ******************** HSS Matrix Approximation *********************
!   Allocate space and initialize variables
    tnk = 0
    DO i = 1, lm
       tnk = tnk+ blk(i)**2
    END DO

    allocate( PH%D(ltr), PH%U(2,ltr), PH%V(2,ltr), PH%B(2,ltr),PH%ph(4,ltr), PH%pd(2,ltr) )
    allocate( H(8*ni*n), DD(tnk),WORK(2*lvl*ni*N) ) ! test only; can be optimized
    call hssexpmm0(PH, ltr)

!   HSS matrix
    pnh = 1
    call cpu_time(time1)
    call mat2hssvd(A, n, tr, ltr, blk, lm, PH, H, DD, tol, prank, Work,pnh )
    call cpu_time(time2)
    time(1) = time2-time1
!    print*, 'construction time: ',time

!   **************** HSS Matrix-Matrix Multiplication *******************
    deallocate( work )
    Allocate( WW(2*LDA,ncol), SX(lvl*n*lda) )

    nc = floor( dble(K) /ncol )
    call cpu_time( time2 )
    DO i = 1, nc
       call fastHssmm(DD,H,PH,blk,B(1:N,(i-1)*ncol+1:i*ncol),tr,ltr,n,ncol,TRE,lvl,WW,SX)
!       call fastHssmm_omp(DD,H,PH,blk,B(1:N,(i-1)*ncol+1:i*ncol),tr,ltr,n,ncol,TRE,lvl )
!       call fastHssmmL( DD,H,PH,blk,B((i-1)*ncol+1:i*ncol,1:N),tr,ltr,lda,ncol,TRE,lvl,WW,SX )
    End Do
    nlast = K - nc * ncol
    IF( nlast .ne. 0 ) THEN
       i = nc+1
       call fastHssmm(DD,H,PH,blk,B((i-1)*ncol+1:N,1:N),tr,ltr,n,nlast,TRE,lvl,WW,SX)
!       call fastHssmm_omp(DD,H,PH,blk,B(1:N,(i-1)*ncol+1:N),tr,ltr,n,nlast,TRE,lvl )
!       call fastHssmmL( DD,H,PH,blk,B((i-1)*ncol+1:N,1:N),tr,ltr,lda,nlast,TRE,lvl,WW,SX )
    END IF
    call cpu_time( time1 )        
    time(2) = time2 - time1
!    print*, 'multiplication time: ', time
   
    deallocate(TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch, blk, tr)
    deallocate(PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd )
    deallocate(H,DD,WW,SX)

  end SUBROUTINE DHSSMM

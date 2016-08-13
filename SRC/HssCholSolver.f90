  SUBROUTINE HssCholSolver( RS,WPCG,A,LDA,b,ni,tol,prank,pp,tau,MaxIt,info )
    use aux_hss
    use mat2hsssvd
! 
    implicit none
!
! .. Scalar Arguments ..
    character :: RS, WPCG
    integer   :: LDA, ni, prank, pp,MaxIt,info
    double precision :: tol, tau
! .. Arrary Arguments ..
    double precision :: A(LDA,*), B(*) 
! 
! Purpose 
! =======
! This routine solves linear equation AX=B, where A is a dense HSS 
! matrix and B may have several columns. An HSS Cholesky factorization
! is performed on matrix A. The details are introduced in our SIMAX 
! paper
! "New efficient and robust HSS Cholesky factorization of SPD matrices " SIMAX, 2011
!
!  N.B.  B only is a vector
! 
! Arguments
! =========
! RS    (input) character 
!       'R': use randomized SVD to do low-rank approximation
!       'S': use SVD to do low-rank approximation
!
! WPCG  (input) character
!       'Y': further use PCG iteration method
!       'N': only use HSS direct solver
! 
! A     (input/output) double precision, dimension(LDA,LDA)
!       An HSS SPD matrix, its Cholesky factorization is to be computed.
!
! LDA   (input) integer, the leading dimension of A
!
! B     (input/output) double precision, dimension(LDA,NCOL)
!       the right-hand-size
! 
! ni    (input) integer, the size of block rows
! 
! ncol  (input) integer, the number of cols of B
!
! tol   (input) double precision, the accuracy parameter for low-rank approximation
!
! prank (input) integer, the prefixed rank parameter
!
! pp    (input) integer, the oversampling parameter
!
! tau   (input) double precision, the tolerance for iteration
!
! MaxIt (input) integer, the allowed maximum iterations
!
! info  (output) integer, the execution information 
! 
! ============
! Written by Shengguo Li, on March 17th, 2014
! ===========================================
!
    integer :: lwork,ierr,ps,ltr,lm,i,N, prec
    logical :: PCGY,PCGN,WNR,WNS
    double precision :: time, time1 
    double precision, allocatable :: R(:),WORK(:), A1(:,:)
!
!   HSS tree
    integer, allocatable :: tr(:), m(:)
    type(hss)            :: RH
    type(ddcell)         :: RT
!
    LOGICAL, EXTERNAL ::  LSAME
!
    N = LDA 
    WNR = LSAME(RS, 'R')
    WNS = LSAME(RS, 'S')
    PCGY = LSAME(WPCG, 'Y')
    PCGN = LSAME(WPCG, 'N')
    
    info = 0
    IF( .not.( WNR .OR. WNS ) ) THEN
       info = -1
    ELSE 
       IF( .not. (PCGY .OR. PCGN ) ) THEN
          info = -2
       END IF
    END IF

    IF( info .ne. 0 ) return

    ! Parameters
    ltr = 0 
    lwork = n*n    ! work space size; for test only; could be improved

    ALLOCATE(m(n/ni+1), tr(2*(n/ni+1)-1),A1(LDA,LDA), stat=ierr )
    IF(ierr .ne. 0 ) THEN
       write(*,*) 'Allocation failed in hsscholsolver'
    END IF
    call npart( n,ni,tr,ltr,m,lm )
    call dlacpy( 'A',LDA,LDA,A,LDA,A1,LDA )


    allocate( RH%dmi(ltr), RH%pp(ltr), RH%p(2,ltr), RH%d(4,ltr),RT%pmi(2,ltr), stat=ierr )
    IF(ierr .ne. 0 ) THEN
       write(*,*) 'Allocation failed in hsscholsolver'
       return
    END IF
    allocate( R(n*n), Work(lwork), stat=ierr ) ! test only; can be optimized
    IF(ierr .ne. 0 ) THEN
       write(*,*) 'Allocation failed in hsscholsolver'
       return
    END IF
    call hssexp0( RH, RT, ltr )

    time = comp_time()
    IF( WNR ) THEN
       ps = prank + pp
       call mat2hsscholsvdr( A1, n, tr, ltr, m, lm, ltr, RH, R, tol, prank, Work, RT, ps )
    ELSE
       call mat2hsscholsvd( A1, n, tr, ltr, m, lm, ltr, RH, R, tol, prank, Work, RT )
    END IF
    time1 = comp_time()
    time1 = time1 - time
    print*, 'HSS Chol factorization time is', time1 

    IF( PCGN ) THEN
       print*,'========== Approx HSS soln (low accuracy) =========='
       call solvecholsvd(RH, R, m, tr, ltr, b, n, Work, lwork )
    ELSE
       print*,'========== PCG  with  HSS preconditioning =========='
       prec = 1
       call cg(A, n, n, b, tau, MaxIt, RH, R, m, tr, ltr, Work, lwork, prec )
    END IF
!
    deallocate( m, tr, A1, R, work )
    deallocate(RH%dmi, RH%pp, RH%p, RH%d,RT%pmi )
!    
!
  END SUBROUTINE HssCholSolver

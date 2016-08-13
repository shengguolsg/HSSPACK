  subroutine SvdHssUpdat( D,F,Z,V,N0,TOL,ni )
     use aux_hss
     use basicMM
     use CauchyHssvd_VPS
!
! ..Scalar Arguments..
     INTEGER, INTENT(IN) :: N0, ni
     DOUBLE PRECISION, INTENT(IN) :: TOL
!
! ..Array Arguments..
     DOUBLE PRECISION, INTENT(IN)    :: F(*)
     DOUBLE PRECISION, INTENT(OUT)   :: D(*)
     DOUBLE PRECISION, INTENT(INOUT) :: Z(*), V(N0,*)
! 
! Purpose
! =========
!  This routine computes the singular values of an updated matrix Ab and updates 
!  its right singular vector matrix. The right singular vector matrix of 
!  M =[diag(F); z^T] is an Cauchy-like matrix, which is approximated by an HSS 
!  matrix to accelerate matrix-matrix multiplication. D are the singular values 
!  of M.  
!
! This routine only updates the right singular vector matrix. It uses structured
! low-rank approximation method to construct HSS matrices. 
!
! ..Parameters..
! ==============
! 
! D  (OUTPUT)  DOUBLE PRECISION Array, DIMENSION(N). 
!    It would contain the updated singular values and is not defined originally.
!
! F  (INPUT) DOUBLE PRECISION Array, DIMENSION(N)
!    the previous singular values.
!
! Z  (INOUT) DOUBLE PRECISION Array, DIMENSION(N)
!    the appended row vector, will be modified by Lowner's Theorem. 
!
! V  (INOUT)  DOUBLE PRECISION Array, DIMENSION(N,N)
!    the original right singular vector matrix.
!    It would be the right singular vector matrix of the updated matrix.
!
! N  (INPUT)  INTEGER, the leading dimension of A
! 
! tol (INPUT) DOUBLE PRECISION
!     Precision parameter for low-rank approximation
!
! ni (INPUT)  INTEGER, row block size
! 
! Futher Details
! See "A Fast and Stable Algorithm for Computing the SVD of an updated matrix", S.G. Li, M. Gu, 
!  Submitted to SISC, 2012
!
! ===============
! Written by Shengguo Li, On Oct. 2nd, 2012
! Modified on April 25th, 2013
! 1) Add some comments and uses the new version of HSS construction algorithm for Cauchy-like
!    matrices.
! ===============
     integer   info, nc, nlast
     logical   hflag
     double precision, allocatable :: DIFL(:), DIFR(:),V1(:,:)

!    HSS tree
     integer, allocatable :: tr(:), m(:), iwork(:)
     type(HTinfo)          :: TRE
       
     real(8), allocatable :: WW(:,:),SX(:),H(:),DD(:),U(:),WORK(:)
     real(8)       :: s1, time, time1,timerrlu,nf1
     integer       :: i, lvl, pnh,tnk,ncol,hlef, ierr, n
     integer       :: ltr, lm, n1, nflops, nswap,K,IV,IW
     type(hssmm)   :: PH

! ********************   Perform Deflation  ********************
     allocate( WORK(N0*N0+N0),IWORK(N0), stat=ierr )
     IF(ierr /= 0 ) Then
        write(*,*) "Allocate failed in svdhssupdat1! "
        stop
     END IF

     IV = 1
     IW = IV + N0*N0
     call mdlasd20( N0,K,F,Z,V,N0,D,WORK(IV),N0,IWORK,WORK(IW),HFLAG,INFO )
!     call testorth(V,N0,N0)

! ********************   Construct HSS tree  ********************
     N = K
     ltr   = 0
     allocate (m(n/ni+1), tr( 2*(n/ni+1)-1 ) )
     call npart(n, ni, tr, ltr, m, lm)

!    Tree information
     lvl = ceiling( log2(ltr+1) ) -1
     n1 = 0.5* (ltr +1 )

     allocate(TRE%ttr(lvl,n1),TRE%lentr(lvl),TRE%lfx(ltr),TRE%ch(2,ltr), stat=ierr )
     IF(ierr /= 0 ) Then
        write(*,*) "Allocate failed in svdhssupcauchy1 ! "
        stop
     END IF

     call GetTreeInfo( TR, LTR, M, M, TRE, lvl )

!   ******************** HSS Matrix Approximation *********************
!    Allocate space and initialize variables
     tnk = 0
     DO i = 1, lm
        tnk = tnk+ m(i)**2
     END DO

     allocate( PH%D(ltr),PH%U(2,ltr),PH%V(2,ltr),PH%B(2,ltr),PH%ph(4,ltr),PH%pd(2,ltr) )
     allocate( H(8*ni*n), DD(tnk),U(n), stat=ierr )
     IF(ierr /= 0 ) Then
        write(*,*) "Allocate failed in svdhsscauchy2 ! "
        stop
     END IF
     call hssexpmm0( PH,ltr )
     
!    HSS matrix
     allocate( DIFL(N), DIFR(N) )
     call mdlasd31( K,D,Z,F,DIFL,DIFR,U,info ) ! U is alpha, D is the new svals, V is the updated Z, F is old svals
     if( info /= 0 )  write(*,*) 'Error in solving secular equations'
!     call testOrthCauchy( D,F,U,Z,DIFL,DIFR,N )
!
     pnh = 1
     hlef = max(5,lm/3)
     nswap = -1   ! signal that this is for updating SVD problem
     call cpu_time(time)
     call cauchy2hssvd( D,F,U,Z,DIFL,DIFR,N,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh,timerrlu,nflops,nswap )
     call cpu_time(time1)
     time = time1-time
     s1  = (1.0D0*pnh )/(ni*n)
     nf1 = (1.0D0*nflops)/(n*n)
     print*, 'Construction time: ',time,'rrlu time is ',timerrlu, 'nflops ', nf1, 'swap ',nswap

!   **************** HSS Matrix-Matrix Multiplication *******************
     ncol = 200
     deallocate( WORK,IWORK )
!
     Allocate(WW(2*N,ncol),SX(2*lvl*N*ncol), stat=ierr )
     IF(ierr /= 0 ) Then
        write(*,*) "Allocate failed in svdhssupcauchy multiplication ! "
        stop
     END IF

     call cpu_time(time)
     nc = floor( dble(N0) /ncol )
     nlast = N0 - nc*ncol
     IF(N .LT. N0 )  THEN
        Allocate( V1(N,N0) )
        call dlacpy('A',N,N0,V,N0,V1,N)
        Do i = 1, nc
           call fastHssmm(DD,H,PH,m,V1(1:N,(i-1)*ncol+1:i*ncol),tr,ltr,n,ncol,TRE,lvl,WW,SX)
        End Do

        IF(nlast .ne. 0) THEN
           i = nc + 1
           call fastHssmm(DD,H,PH,m,V1(1:N,(i-1)*ncol+1:N0),tr,ltr,n,nlast,TRE,lvl,WW,SX)
        END IF
        call dlacpy('A',N,N0,V1,N,V,N0)
        deallocate( V1 )
     ELSE
        Do i = 1, nc
           call fastHssmm(DD,H,PH,m,V(1:N,(i-1)*ncol+1:i*ncol),tr,ltr,n,ncol,TRE,lvl,WW,SX)
        End Do
        IF(nlast .ne. 0) THEN
           i = nc + 1
           call fastHssmm(DD,H,PH,m,V(1:N,(i-1)*ncol+1:n),tr,ltr,n,nlast,TRE,lvl,WW,SX)
        END IF
     END IF
     call cpu_time(time1)
     time = time1-time
     print*, 'Multiplication time: ',time

     deallocate(TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch)
     deallocate(PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd)
     deallocate(H,DD,m,tr,DIFL,DIFR,WW,SX)

   end subroutine SvdHssUpdat

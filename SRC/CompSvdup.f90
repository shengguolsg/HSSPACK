module CompSvdup
! This procedure is random svd + column majority. 
  implicit none

contains

!!!!!!!!! HssMV
   subroutine SvdHssup(A, N, V, RS, ni )
     use aux_hss
     use ConstructHssd
!
! .. Scalar Arguments ..
     INTEGER  :: N, ni
     Character(len=1) :: RS
! .. Array Arguments ..
     DOUBLE PRECISION :: A(N,*), V(N,*)
! 
! Purpose
! =========
!  This routine firstly constructs an HSS approximation of matrix A, and 
!  uses this HSS approximation to multiply with a matrix V. It is firstly
!  used for updating SVD problem. 
!
! ..Parameters..
! ==============
! 
! A  (INOUT)  DOUBLE PRECISION Array, DIMENSION(N,N)
!    An HSS matrix with off-diagonal low-rank structure.
!    It will be destroyed after calling
!
! N  (INPUT)  INTEGER, the leading dimension of A
! 
! V  (INOUT)  DOUBLE PRECISION Array, DIMENSION(N,N)
!    The return matrix, V = A*V.
!
! RS (INPUT)   CHARACTER, LENGTH =1
!    'R' : using random sampling algorithm to compress the matrix
!    'S' : using SVD to compress the matrix
!
! ni (INPUT)  INTEGER, row block size
! 
! ===============
! Written by Shengguo Li, On Aug. 24th, 2012
! Modified on April 25th, 2013
! 1) Changed its name from SvdHssup to HssMV. 
! ===============

!  .. External Functions
      LOGICAL     LSAME
      EXTERNAL    LSAME

      LOGICAL  :: CR

      ! HSS tree
      integer, allocatable :: tr(:), m(:)
      type(HTinfo)          :: TRE

      real(8), allocatable :: H(:), Work(:), DD(:)
      real(8), allocatable :: WW(:,:), SX(:)
      real(8)       :: tol, s1, time, time1
      integer       :: pp, i, tnk, lvl, ncol, pnh, prank
      integer       :: ltr, lm, n1, nc, nlast
      type(hssmm)   :: PH

! ********************   Construct HSS tree  ********************
      tol   = 1.0d-15   ! tolerance      
      ltr   = 0
      prank = 50    ! prank never used
      pp = 2*ni+10
      allocate (m(n/ni+1), tr( 2*(n/ni+1)-1 ))
      call npart(n, ni, tr, ltr, m, lm)

!     Tree information
      lvl = ceiling( log2(ltr+1) ) -1
      n1 = 0.5* (ltr +1 )

      allocate(TRE%ttr(lvl,n1), TRE%lentr(lvl), TRE%lfx(ltr), TRE%ch(2,ltr) )
      call GetTreeInfo(TR, LTR, M, M, TRE, lvl)       

!   ******************** HSS Matrix Approximation *********************
!      Allocate space and initialize variables
       tnk = 0
       DO i = 1, lm
          tnk = tnk+ m(i)**2
       END DO

       allocate (PH%D(ltr), PH%U(2,ltr), PH%V(2,ltr), PH%B(2,ltr),PH%ph(4,ltr), PH%pd(2,ltr))
       allocate (H(8*ni*n), DD(tnk), WORK(2*lvl*ni*n) ) ! H may not big enough
       call hssexpmm0(PH, ltr)
       
!      HSS matrix
       pnh = 1
       CR = LSAME(RS, 'R')
       call cpu_time(time)
       IF( CR ) Then
          tol = 1.0D-13
          call mat2hssvdr(A, n, tr, ltr, m, lm, PH, H, DD, tol, prank, Work, pp,lvl,pnh )
       ELSE
          call mat2hssvd(A, n, tr, ltr, m, lm, PH, H, DD, tol, prank, Work, pnh )
       END IF
       call cpu_time(time1)
       time = time1-time
       s1 = (1.0D0*pnh )/(ni*n)
       print*, 'construction time: ',time, s1

!   **************** HSS Matrix-Matrix Multiplication *******************
       deallocate( work)
       ncol = 200
       Allocate(WW(2*n,ncol), SX(2*lvl*n*ncol))

       nc = floor( dble(n) /ncol )
       call cpu_time(time)
       Do i = 1, nc
          call fastHssmm(DD,H,PH,m,V(1:N,(i-1)*ncol+1:i*ncol),tr,ltr,n,ncol,TRE,lvl,WW,SX)
       End Do
       nlast = n - nc*ncol
       IF(nlast .ne. 0) THEN
          i = nc + 1
          call fastHssmm(DD,H,PH,m,V(1:N,(i-1)*ncol+1:n),tr,ltr,n,nlast,TRE,lvl,WW,SX)
       END IF
       call cpu_time(time1)
       time = time1-time
       print*, 'Multiplication time: ',time

       deallocate(TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch)
       deallocate(PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd)
       deallocate(H,DD,m,tr)

     end subroutine SvdHssup

!!!!!!!
  subroutine SvdupdatAll(U,D,V,va,U1,D1,V1,LDA,CDA,tol,ni)
    use aux_hss
    use basicMM
    use CauchyHssvd_VPS
!
! .. Scalar Arguments ..
    integer LDA, CDA,ni
    double precision tol
! .. Array Arguments ..
    double precision :: D(*), D1(*), va(*)
    double precision :: U(LDA,CDA), V(CDA,CDA), U1(LDA+1,CDA),V1(CDA,CDA)
!
! Purpose
! =========
! Assume matrix A is m-by-n and its SVD is A=U*diag(D)*V, U is m-by-n, D is n-by-n, and V is n-by-n.
! We compute the SVD of Ab=[A; va^T], Ab = U1*diag(D1)*V1, m >= n and va is a one-by-n vector. 
! U1 is mp-by-n, D1 is n-by-n and V1 is n-by-n.  
! 
! .. Parameters ..
! U  (in) double precision array, dimension(m,m), the left singular vector matrix of A;
!
! D  (in) double precision array, dimension(n), the singular values of A, in increasing order;
!
! V  (in) double precision array, dimension(n,n), the right singular vector matrix of A;
!
! va (in) double precision array, dimension(n), the updated vector; 
!
! U1 (out) double precision array, dimension(m+1,m+1), the left singular vector matrix of Ab;
!
! D1 (out) double precision array, dimension(n), the singular values of Ab, in increasing order;
! 
! V1 (out) double precision array, dimension(n,n), the right singular vector matrix of Ab;
!
! LDA (in) integer, the row dimension of A;
! 
! CDA (in) integer, the column dimension of A;
! 
! tol (in) double precision, accuracy parameter for HSS approximation;
!
! ni  (in) integer, the blksizes of nodes at the bottom level of HSS tree
!
! ====================
! Written by Shengguo Li, on Nov. 21st, Changsha, 2012
! ====================
       double precision Z(CDA),ZZ(CDA),ALPHA_L(CDA),ALPHA_R(CDA),DIFL(CDA),DIFR(CDA)

!      HSS tree
       integer, allocatable :: tr(:), m(:)
       type(HTinfo)          :: TRE
       
       real(8), allocatable :: WW(:,:),SX(:),H(:),DD(:)
       real(8)       :: s1,time,time1,timerrlu,nf1,temp
       integer       :: i,j,lvl,pnh,tnk,ncol,N2,nc,nlast
       integer       :: ltr,lm,info,n1,nflops,nswap,N,ldap
       type(hssmm)   :: PH

!*****************************************************
!   Write z as the combination of V, z = V*a         *
!*****************************************************
       N = CDA
       LDAP = LDA+1
       call dgemv('N',N,N,1d0,V,N,va,1,0d0,Z,1)  ! appended vector Z

!****************************************************
!    Computing the singular values of Ab            *
!****************************************************
       call dlasrt('I',N,D,info)   ! increasingly order D
       N2 = N/2
       Do J =1, N2                       ! reverse Z
          temp = Z(J)
          Z(J) = Z(N-J+1)
          Z(N-J+1) = temp
       END DO

       call mdlasd32(N,D1,Z,D,DIFL,DIFR,ALPHA_L,ALPHA_R,info)
       
! ***************************************************
!             Construct HSS tree                    *
! ***************************************************
       ltr   = 0
       allocate (m(n/ni+1), tr( 2*(n/ni+1)-1 ))
       call npart(n, ni, tr, ltr, m, lm)

!     Tree information
      lvl = ceiling( log2(ltr+1) ) -1
      n1 = 0.5* (ltr +1 )
      
      allocate(TRE%ttr(lvl,n1), TRE%lentr(lvl), TRE%lfx(ltr), TRE%ch(2,ltr) )
      call GetTreeInfo(TR, LTR, M, M, TRE, lvl)       

! ***************************************************
!             HSS Matrix Approximation              *
! ***************************************************
!      Allocate space and initialize variables
       tnk = 0
       DO i = 1, lm
          tnk = tnk+ m(i)**2
       END DO

       allocate (PH%D(ltr), PH%U(2,ltr), PH%V(2,ltr), PH%B(2,ltr),PH%ph(4,ltr), PH%pd(2,ltr))
       allocate (H(8*ni*n), DD(tnk) ) ! test only; can be optimized
       call hssexpmm0(PH, ltr)

! ****************************************************
!             Right singular vector matrix           *
! ****************************************************
       pnh = 1
       ZZ(1:N) = Z(1:N)
       call cpu_time(time)
       call cauchy2hssvd(D,D1,ZZ,ALPHA_R,DIFL,DIFR,N,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh, &
                         timerrlu,nflops,nswap)
       call cpu_time(time1)
       time = time1-time
       s1  = (1.0D0*pnh )/(ni*n)
       nf1 = (1.0D0*nflops)/(n*n)
!       write(*,*), 'Right const time: ',time,'rrlu time ',timerrlu, 'nflops ', nf1, 'swap ', nswap, s1

!   ************ HSS Matrix-Matrix Multiplication ************
        ncol = 200
        Allocate(WW(2*n,ncol), SX(2*lvl*n*ncol))
        
        nc = floor( dble(n) /ncol )
        V1(1:N,1:N) = V(1:N,1:N)
        call cpu_time(time)
        Do i = 1, nc
           call fastHssmm(DD,H,PH,m,V1(1:N,(i-1)*ncol+1:i*ncol),tr,ltr,n,ncol,TRE,lvl,WW,SX)
        End Do
        nlast = N - nc* ncol
        if( nlast .ne. 0 ) then
           i = nc+1 
           call fastHssmm(DD,H,PH,m,V1(1:N,(i-1)*ncol+1:N ),tr,ltr,n,nlast,TRE,lvl,WW,SX)
        end if
        call cpu_time(time1)
        time = time1-time
!        write(*,*), 'Right mult time: ',time

! **************************************************
!            Left singular vector matrix           *
! **************************************************
       call hssexpmm0(PH, ltr)
       DD = 0d0
       pnh = 1
       nflops = 0
       nswap = 0
       ZZ = D(1:N)*Z(1:N)             ! define \hat{Z} for left singular vectors 

       ! (LDA+1)th row of U1       
       do i = 1, N
           U1(LDAP,i) = ALPHA_L(N+1-i)
        end do
       call cpu_time(time)
       call cauchy2hssvd(D,D1,ZZ,ALPHA_L,DIFL,DIFR,N,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh, &
                         timerrlu,nflops,nswap)
       call cpu_time(time1)
       time = time1-time
       s1  = (1.0D0*pnh )/(ni*n)
       nf1 = (1.0D0*nflops)/(n*n)
       print*, 'Left const time: ',time,'rrlu time ',timerrlu, 'nflops ', nf1, 'swap ', nswap, s1

! ***********************   Comments        ******************************
!  After calling cauchy2hssvd, the orders of coefficients, alpha and Z,  *
!  would be permuted to any order.                                       *
! ************************************************************************

! ****************  Update the left singular vectors **********
        call cpu_time(time)
        U(1:N,1:N) = transpose( U(1:N,1:N) )
        Do i = 1, nc
           call fastHssmm(DD,H,PH,m,U(1:N,(i-1)*ncol+1:i*ncol),tr,ltr,n,ncol,TRE,lvl,WW,SX)
        End Do
        if( nlast .ne. 0 ) then
           i = nc+1 
           call fastHssmm(DD,H,PH,m,U(1:N,(i-1)*ncol+1:N ),tr,ltr,n,nlast,TRE,lvl,WW,SX)
        end if
        call cpu_time(time1)
        time = time1-time
        print*, 'Left mult time: ',time

        U(1:N,1:N) = transpose( U(1:N,1:N) )
        call dlacpy('A',N,N,U,LDA,U1,LDAP)

!!$        ! (n+2)th-mth columns of U1 
!!$        call dlacpy('A',LDA,(LDA-CDA),U(1,CDA+1),U1(1,CDA+1),LDA+1)
        
! *** reorder the new singular values ***
        call dlasrt('D',N,D1,info)

        deallocate(TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch)
        deallocate(PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd)
        deallocate(H,DD,m,tr,WW,SX)

  end subroutine SvdupdatAll

!!!!!!
  SUBROUTINE MDLASD32( K,D,Z,DSIGMA,DIFL,DIFR,ALPHA_L,ALPHA_R,INFO )
!
    implicit none

!     .. Scalar Arguments ..
      INTEGER            INFO, K
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D(*),Z(*),DSIGMA(*),DIFL(*),DIFR(*),ALPHA_L(*),ALPHA_R(*)
!     ..
!
!  Purpose
!  =======
!
!  MDLASD3 finds the square roots of the roots of the secular equation,
!  as defined by the values in DSIGMA and Z. It makes the appropriate
!  calls to DLASD4.
!
!  DLASD8 is called from DLASD6.
!
!  Arguments
!  =========
!
!  K       (input) INTEGER
!          The number of terms in the rational function to be solved
!          by DLASD4.  K >= 1.
!
!  D       (output) DOUBLE PRECISION array, dimension ( K )
!          On output, D contains the updated singular values.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension ( K )
!          On entry, the first K elements of this array contain the
!          components of the deflation-adjusted updating row vector.
!          On exit, Z is updated.
!
!  DSIGMA  (input/output) DOUBLE PRECISION array, dimension ( K )
!          On entry, the first K elements of this array contain the old
!          roots of the deflated updating problem.  These are the poles
!          of the secular equation.
!          On exit, the elements of DSIGMA may be very slightly altered
!          in value.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension at least 3 * K
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J
      DOUBLE PRECISION   RHO, TEMP,DIFLJ,DIFRJ,DJ,DSIGJ,DSIGJP
      DOUBLE PRECISION, allocatable :: WORK(:)
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASCL, DLASD4, DLASET, XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      Allocate( Work(3*K) )
!      
      IF( K.LT.1 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASD8', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
         DIFL( 1 ) = D( 1 )
         RETURN
      END IF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
      END DO
!
!     Book keeping.
!
      IWK1 = 1
      IWK2 = IWK1 + K
      IWK3 = IWK2 + K
      IWK2I = IWK2 - 1
      IWK3I = IWK3 - 1
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Initialize WORK(IWK3).
!
      CALL DLASET( 'A', K, 1, ONE, ONE, WORK( IWK3 ), K )
!
!     Compute the updated singular values, the arrays DIFL, DIFR,
!     and the updated Z.
!
      DO J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ), &
                     WORK( IWK2 ), INFO )
!
!        If the root finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'DLASD4', -INFO )
            RETURN
         END IF
         WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
         DIFL( J ) = -WORK( J )  
         DIFR( J ) = -WORK( J+1 )
         DO I = 1, J - 1
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )*      &
                             WORK( IWK2I+I ) / ( DSIGMA( I )- &
                             DSIGMA( J ) ) / ( DSIGMA( I )+   &
                             DSIGMA( J ) )
         END DO
         DO I = J + 1, K
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )*       &
                             WORK( IWK2I+I ) / ( DSIGMA( I )- &
                             DSIGMA( J ) ) / ( DSIGMA( I )+   &
                             DSIGMA( J ) )
         END DO
      END DO
!
!     Compute updated Z.
!
      DO I = 1, K
         Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
      END DO
!
!     Compute the left and right scalars
!
      DO J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         IF( J.LT.K ) THEN
            DIFRJ = -DIFR( J )
            DSIGJP = -DSIGMA( J+1 )
         END IF
         WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
         DO I = 1, J - 1
            WORK( I ) = Z( I ) / ( ( DSIGMA( I )+ DSIGJ )-DIFLJ ) / ( DSIGMA( I )+DJ )
         END DO
         DO I = J + 1, K
            WORK( I ) = Z( I ) / ( ( DSIGMA( I )+DSIGJP )+DIFRJ ) / ( DSIGMA( I )+DJ )
         END DO
         TEMP = DNRM2( K, WORK, 1 )
         ALPHA_R( J ) = one / TEMP
      END DO

      WORK(K+1)= -1.0D0
      DO J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         IF( J.LT.K ) THEN
            DIFRJ = -DIFR( J )
            DSIGJP = -DSIGMA( J+1 )
         END IF
         WORK( J ) = -DSIGMA(J)*Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
         DO I = 1, J - 1
            WORK( I ) = DSIGMA(I)*Z( I ) / ( ( DSIGMA( I )+ DSIGJ )-DIFLJ ) / ( DSIGMA( I )+DJ )
         END DO
         DO I = J + 1, K
            WORK( I ) = DSIGMA(I)*Z( I ) / ( ( DSIGMA( I )+DSIGJP )+DIFRJ ) / ( DSIGMA( I )+DJ )
         END DO
         TEMP = DNRM2( K+1, WORK, 1 )
         ALPHA_L( J ) = one / TEMP
      END DO

      deallocate( WORK )

!     End of MDLASD32
!
    END SUBROUTINE MDLASD32

!!!!!!!!
  subroutine fastHssmm(DD,H,PH,M,X,TR,LTR,LDX,NCOL,TRE,LVL,WW,SX)
    use aux_hss
!
!  .. Scalar Arguments ..
       INTEGER           ::  LTR, LDX, NCOL, LVL

!  .. Array Arguments ..
       DOUBLE PRECISION  ::  DD(*), H(*), X(LDX,*), WW(2*LDX,*), SX(*)
       INTEGER           ::  M(*), TR(*)
       TYPE(HSSMM)       ::  PH  
       TYPE(HTinfo)      ::  TRE
!
!  Purpose    
!  ========
!   Compute HSS matrix with a general matrix. This HSS matrix is in post ordering. 
!   
! .. Argumrents ..
! DD    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators D's.
! 
! H    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators, U, V and B's. 
! 
! PH   (INPUT)  TYPE HSSMM, stores the information of all these generators.
! 
! M    (INPUT)  INTEGER array, DIMENSION(:), stores the sizes of diagonal blocks.
! 
! X    (INPUT/OUTPUT) DOUBLE PRECISION array, DIMENSION(:,:), the input matix.
!
! TR    (INPUT) INTEGER array, DIMENSION(ltr), The HSS tree in post ordering.
!
! LTR  (INPUT)  INTEGER, the length of HSS Tree, LTR.
! 
! LDX  (INPUT) INTEGER, the row dimension of X
!
! NCOL (INPUT) INTEGER, the column dimension of X
!
! TRE  (INPUT) TYPE HTinfo, contains some information about HSS tree
! 
! LVL  (INPUT) INTEGER, the total levels of HSS tree
! 
! ====================
! Written by Shengguo Li, on Aug 9th, 2012
!
! Modified on Aug. 10th, 2012.
!   1) ww may be changed to 1D array, and the indexing would be much more complicated.
!   2) Some indexing l and lx may be done duing processing. Modify it later and make it work first.  
! ====================

! .. Local Parameters ..
       INTEGER    :: n1, i, it, lt, nj, flag, LL, pww, ldw, &
                     pdx, strid_m, j, k, leafi, ni, pnx, nstart1, nstart2

! .. Array Parameters ..
       INTEGER     :: lx(ltr), l(ltr), px(ltr)
!
! .. Executable Statements ..
       n1 = 0.5*(ltr+1)
       ldw = 2*ldx

!*********************************************************
!                      leaf nodes                        *
!*********************************************************
! lx(i): row partition of x
       lt = 1
       it = 1
       lx = 0
       l  = 0
       n1 = TRE%lentr(lvl)
       do i = 1, n1
          lx(i) = lt
          lt = lt+m(i)

          ni = TRE%ttr(lvl, i)
          l(ni) = it
          it = it+PH%V(1,ni)
       end do

       pnx = 1
       do i = 1, n1
          ni = TRE%ttr(lvl,i)  ! can change it in column major
          call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0,H( PH%ph(2,ni) ),PH%V(1,ni), & 
               X(lx(i),1), ldx, 0d0,ww(pnx,1), ldw)
          ! ww = V{i}'*Xi
          pnx = pnx + PH%V(1,ni)
       end do
!
       pdx = 1
       do i = 1, n1
          ni = TRE%ttr(lvl, i)
          
          if (ni .eq. TRE%ch(1,tr(ni) ) ) then 
             nj = TRE%ch(2, tr(ni) )    ! adjacent of ni
          else
             nj = TRE%ch(1, tr(ni) )
          end if
          
          ! xt{ni}
          call dgemm('N','N',PH%B(1,ni),ncol,PH%B(2,ni),1d0,H( PH%ph(3,ni) ),&
               PH%B(1,ni),ww(l(nj),1),ldw, 0d0, SX(pdx),PH%B(1,ni)) 
          px(ni) = pdx
          pdx = pdx + PH%B(1,ni)*ncol 
       end do

!***************************************************************
!                    Forward Traversal                         *
!***************************************************************
       flag = 0
       pww = 0
       DO LL = lvl-1, 1, -1
          n1 = TRE%lentr(LL)
          
          ! merge and determine the row dim
          lt = 1
          it = 1
          DO i =1, n1
             ni = TRE%ttr(ll, i)
             IF ( TRE%ch(1,ni) .ne. 0  ) THEN
                lx(i)  = lt
                lt = lt+PH%V( 1, TRE%ch(1,ni) )+PH%V( 1, TRE%ch(2,ni) )
             ELSE
                if( flag .eq. 0) then
                   flag  = 1
                   leafi = i   ! leafi is a leaf node, above is full binary tree
                end if
             END IF
             
             l(ni) = it    ! l(ni) records the starting row of ww(ni)
             it = it+PH%V(1,ni)
          END DO !(i)

          ! update ww
          Nstart1 = max( pww * ldx, 1)
          Nstart2 = max( (1-pww) * ldx, 1)
          IF ( flag .eq. 0 ) THEN
             DO i = 1, n1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0,H( PH%ph(2,ni) ), & 
                     PH%V(1,ni),ww(Nstart1+lx(i)-1,1),ldw,0d0,ww(Nstart2+l(ni)-1,1),ldw)
             END DO
          ELSE
             DO i = 1, leafi-1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0,H( PH%ph(2,ni) ), & 
                     PH%V(1,ni),ww(Nstart1+lx(i)-1,1),ldw,0d0,ww(Nstart2+l(ni)-1,1),ldw)
             END DO
             
             DO i = leafi, n1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0, H( PH%ph(2,ni) ), & 
                     PH%V(1,ni), X(TRE%lfx(ni),1), ldx, 0d0,ww(Nstart2+l(ni)-1,1), ldw )
             END DO
             flag = 0
          END IF !(flag)          

          do i = 1, n1
             ni = TRE%ttr(LL, i)
          
             if (ni .eq. TRE%ch(1,tr(ni) ) ) then 
                nj = TRE%ch(2, tr(ni) )    ! adjacent of ni
             else
                nj = TRE%ch(1, tr(ni) )
             end if
          
             ! xt{ni}
             call dgemm('N','N',PH%B(1,ni),ncol,PH%B(2,ni),1d0,H( PH%ph(3,ni) ), &
                  PH%B(1,ni),ww(Nstart2+l(nj)-1,1), ldw, 0d0,SX(pdx), PH%B(1,ni) )      
             px(ni) = pdx
             pdx = pdx + PH%B(1,ni)*ncol 
          end do
          pww = 1 - pww

       END DO ! (LL)

!*********************************************************
!                 Backward Traversal                     *
!*********************************************************

       ! First level
       LL = 1
       lt = 1
       DO i = 1, 2
          ni = TRE%ttr(LL, i)
          call dgemm('N','N',PH%U(1,ni),ncol,PH%U(2,ni),1d0, H( PH%ph(1,ni) ), & 
               PH%U(1,ni),SX( px(ni) ), PH%U(2,ni),0d0, ww(lt,1),ldw)
          
          IF ( TRE%ch(1,ni) .eq. 0 ) THEN   ! add it to X
             call dgemm('N','N',PH%D(ni),ncol,PH%D(ni), 1d0, DD( PH%pd(1,ni) ), PH%D(ni), & 
                  X(TRE%lfx(ni), 1),ldx, 1d0, ww(lt,1), ldw )
             call dlacpy('A', PH%D(ni),ncol,ww(lt,1),ldw,X(TRE%lfx(ni),1), ldx)
          END IF
          lt = lt + PH%U(1,ni)

       END DO
       
       ! Traverse downward
       DO LL = 2, lvl
          n1 = TRE%lentr(LL)

          ! update the data of children
          lt =1
          DO i = 1, n1
             ni = TRE%ttr(LL,i)

             strid_m = PH%B(1,ni)
             it = px(ni)-1
             DO k = 1, ncol
                DO j = 1, strid_m         ! use 2D array there is a huge stride.
                   SX(it+j+(k-1)*strid_m ) = SX(it+j+(k-1)*strid_m ) + ww(lt+j-1,k)
                END DO
             END DO ! matrix-matrix addition
             lt = lt+strid_m
          END DO
!
          ! extend
          lt = 1
          DO i = 1, n1
             ni = TRE%ttr(LL, i)
             call dgemm('N','N',PH%U(1,ni),ncol,PH%U(2,ni),1d0, H( PH%ph(1,ni) ),& 
                  PH%U(1,ni),SX( px(ni) ), PH%U(2,ni),0d0, ww(lt,1),ldw )
             
             IF( TRE%ch(1,ni) .eq. 0 ) THEN   ! add it to X
                call dgemm('N','N',PH%D(ni),ncol,PH%D(ni), 1d0, DD( PH%pd(1,ni) ), PH%D(ni), & 
                     X(TRE%lfx(ni), 1),ldx, 1d0, ww(lt,1), ldw )
                call dlacpy('A', PH%D(ni),ncol,ww(lt,1),ldw,X(TRE%lfx(ni),1), ldx)
             END IF             
             lt = lt + PH%U(1,ni)
          END DO

       END DO ! (LL)

     end subroutine fastHssmm

!!!!!!!!
  subroutine fastHssmmL(DD,H,PH,M,X,TR,LTR,LDX,NCOL,TRE,LVL,WW,SX)
    use aux_hss
!
!  .. Scalar Arguments ..
       INTEGER           ::  LTR, LDX, NCOL, LVL

!  .. Array Arguments ..
       DOUBLE PRECISION  ::  DD(*), H(*), X(LDX,*), WW(LDX,*), SX(*)
       INTEGER           ::  M(*), TR(*)
       TYPE(HSSMM)       ::  PH  
       TYPE(HTinfo)      ::  TRE
!
!  Purpose    
!  ========
!   Compute the multiplication of a general matrix with an HSS matrix. This HSS matrix is in post ordering. 
!   This algorithm is a sequential code and based on the assumption that row and column partitions are same.
!   
! .. Argumrents ..
! DD    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators D's.
! 
! H    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators, U, V and B's. 
! 
! PH   (INPUT)  TYPE HSSMM, stores the information of all these generators.
! 
! M    (INPUT)  INTEGER array, DIMENSION(:), stores the sizes of diagonal blocks.
! 
! X    (INPUT/OUTPUT) DOUBLE PRECISION array, DIMENSION(:,:), the input matix.
!
! TR    (INPUT) INTEGER array, DIMENSION(ltr), The HSS tree in post ordering.
!
! LTR  (INPUT)  INTEGER, the length of HSS Tree, LTR.
! 
! LDX  (INPUT) INTEGER, the row dimension of X
!
! NCOL (INPUT) INTEGER, the column dimension of X
!
! TRE  (INPUT) TYPE HTinfo, contains some information about HSS tree
! 
! LVL  (INPUT) INTEGER, the total levels of HSS tree
! 
! ====================
! Written by Shengguo Li, on Nov 20th, 2012
!
! ====================

! .. Local Parameters ..
       INTEGER    :: n1, i, it, lt, nj, flag, LL, pww, ldw, &
                     pdx, strid_m, j, k, leafi, ni, pnx, nstart1, nstart2

! .. Array Parameters ..
       INTEGER    :: lx(ltr), l(ltr), px(ltr)

! .. Executable Statements ..
       n1 = 0.5*(ltr+1)
       ldw = ldx

!*********************************************************
!                      leaf nodes                        *
!*********************************************************
! lx(i): row partitions of x
       lt = 1
       it = 1
       lx = 0
       l  = 0
       n1 = TRE%lentr(lvl)
       do i = 1, n1
          lx(i) = lt
          lt = lt+m(i)

          ni = TRE%ttr(lvl,i)
          l(ni) = it
          it = it+PH%U(2,ni)
       end do

       pnx = 1
       do i = 1, n1
          ni = TRE%ttr(lvl,i)  ! can change it to column major
          call dgemm('N','N',ldx,PH%U(2,ni),m(i),1d0,X(1,lx(i)),ldx,H( PH%ph(1,ni) ),PH%U(1,ni), & 
               0d0,ww(1,pnx),ldw)          ! ww = Xi * U(i)
          pnx = pnx + PH%U(2,ni)
       end do

       pdx = 1
       do i  = 1, n1
          ni = TRE%ttr(lvl, i)
          
          if (ni .eq. TRE%ch(1,tr(ni) ) ) then 
             nj = TRE%ch(2, tr(ni) )    ! adjacent of ni
          else
             nj = TRE%ch(1, tr(ni) )
          end if
          
          ! xt{ni}: f_i in note
          call dgemm('N','N',ldw,PH%B(2,nj),PH%B(1,nj),1d0,ww(1,l(nj)),ldw,H( PH%ph(3,nj) ),&
               PH%B(1,nj),0d0, SX(pdx),ldx) 
          px(ni) = pdx
          pdx = pdx + PH%B(2,nj)*ldx
       end do

!***************************************************************
!                    Forward Traversal                         *
!***************************************************************
       flag = 0
       pww = 0
       DO LL = lvl-1, 1, -1
          n1 = TRE%lentr(LL)
          
          ! merge and determine the column dim
          lt = 1
          it = 1
          DO i =1, n1
             ni = TRE%ttr(ll, i)
             IF ( TRE%ch(1,ni) .ne. 0  ) THEN
                lx(i)  = lt
                lt = lt+PH%U( 2, TRE%ch(1,ni) )+PH%U( 2, TRE%ch(2,ni) )
             ELSE
                if( flag .eq. 0) then
                   flag  = 1
                   leafi = i   ! leafi is a leaf node, above is full binary tree
                end if
             END IF
             
             l(ni) = it    ! l(ni) records the starting column of ww(ni)
             it = it+PH%U(2,ni)
          END DO !(i)

          ! update ww i.e. g_i in our note
          Nstart1 = max( pww * ncol, 1)
          Nstart2 = max( (1-pww) * ncol, 1)
          IF ( flag .eq. 0 ) THEN          !there is no leaf node except bottom level
             DO i = 1, n1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',ldx,PH%U(2,ni),PH%U(1,ni),1d0,ww(1,Nstart1+lx(i)-1),ldw,H( PH%ph(1,ni) ), & 
                     PH%U(1,ni),0d0,ww(1,Nstart2+l(ni)-1),ldx)
             END DO
          ELSE
             DO i = 1, leafi-1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',ldx,PH%U(2,ni),PH%U(1,ni),1d0,ww(1,Nstart1+lx(i)-1),ldw,H( PH%ph(1,ni) ), & 
                     PH%U(1,ni),0d0,ww(1,Nstart2+l(ni)-1),ldx)
             END DO
             
             DO i = leafi, n1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',ldx,PH%U(2,ni),PH%U(1,ni),1d0,X(1,TRE%lfx(ni)),ldx,H( PH%ph(1,ni) ), & 
                     PH%U(1,ni), 0d0,ww(1,Nstart2+l(ni)-1), ldw )
             END DO
             flag = 0
          END IF !(flag)          

          do i = 1, n1
             ni = TRE%ttr(LL, i)
          
             if (ni .eq. TRE%ch(1,tr(ni) ) ) then 
                nj = TRE%ch(2, tr(ni) )    ! adjacent of ni
             else
                nj = TRE%ch(1, tr(ni) )
             end if
          
             ! xt{ni}
             call dgemm('N','N',ldx,PH%B(2,nj),PH%B(1,nj),1d0,ww(1,Nstart2+l(nj)-1),ldw,H( PH%ph(3,nj) ), &
                  PH%B(1,nj), 0d0,SX(pdx), ldx )      
             px(ni) = pdx
             pdx = pdx + PH%B(2,nj)*ldx
          end do
          pww = 1 - pww

       END DO ! (LL)

!*********************************************************
!                 Backward Traversal                     *
!*********************************************************
       ! First level
       LL = 1
       lt = 1
       DO i = 1, 2
          ni = TRE%ttr(LL, i)
          call dgemm('N','N',ldx,PH%V(2,ni),PH%V(1,ni),1d0,SX( px(ni) ),ldx, H( PH%ph(2,ni) ), & 
               PH%V(1,ni),0d0, ww(1,lt),ldw)
          
          IF ( TRE%ch(1,ni) .eq. 0 ) THEN   ! add it to X
             call dgemm('N','N',ldx,PH%D(ni),PH%D(ni),1d0,X(1,TRE%lfx(ni)),ldx,DD( PH%pd(1,ni) ),PH%D(ni), & 
                  1d0, ww(1,lt), ldw )
             call dlacpy('A', ldx,PH%D(ni),ww(1,lt),ldw,X(1,TRE%lfx(ni)), ldx)
          END IF
          lt = lt + PH%V(2,ni)
       END DO
       
       ! Traverse downward
       DO LL = 2, lvl
          n1 = TRE%lentr(LL)

          ! update the data of children
          lt =1
          DO i = 1, n1
             ni = TRE%ttr(LL,i)

             strid_m = PH%V(1,ni)
             it = px(ni)-1
             DO k = 1, strid_m
                DO j = 1, ldx         ! use 2D array there is a huge stride.
                   SX(it+j+(k-1)*ldx ) = SX(it+j+(k-1)*ldx ) + ww(j,lt+k-1)
                END DO
             END DO ! matrix-matrix addition
             lt = lt+strid_m
          END DO
!
          ! extend
          lt = 1
          DO i = 1, n1
             ni = TRE%ttr(LL, i)
             call dgemm('N','N',ldx,PH%V(2,ni),PH%V(1,ni),1d0, SX( px(ni) ), ldx, H( PH%ph(2,ni) ),& 
                  PH%V(1,ni),0d0, ww(1,lt),ldw )
             
             IF( TRE%ch(1,ni) .eq. 0 ) THEN   ! add it to X
                call dgemm('N','N',ldx,PH%D(ni),PH%D(ni),1d0,X(1,TRE%lfx(ni)),ldx,DD( PH%pd(1,ni) ),PH%D(ni), & 
                     1d0, ww(1,lt), ldw )
                call dlacpy('A', ldx,PH%D(ni),ww(1,lt),ldw,X(1,TRE%lfx(ni)), ldx)
             END IF             
             lt = lt + PH%V(2,ni)
          END DO

       END DO ! (LL)

     end subroutine fastHssmmL


  end module CompSvdup

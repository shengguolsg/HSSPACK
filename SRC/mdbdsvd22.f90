  SUBROUTINE MDBDSVD22( M,N,K,SB,D,ZZ,U,LDU,VT,LDVT,DSIGMA,U2,LDU2, &
                VT2,LDVT2,IDXP,IDXQ,HFLAG,INFO )
!
!  -- HSSPACK computational routine (version 0.0.1) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Def. Tech., China.
!     January 2013
!
        IMPLICIT NONE
!     .. Scalar Arguments ..
      LOGICAL            HFLAG
      INTEGER            INFO,K,SB,LDU,LDU2,LDVT,LDVT2,N,M
!     ..
!     .. Array Arguments ..
      INTEGER            IDXP( * ), IDXQ( * )
      DOUBLE PRECISION   D( * ), DSIGMA( * ), U( LDU,* ), &
                         U2( LDU2,* ), VT( LDVT,* ), VT2( LDVT2,* ), &
                         ZZ( M,* )
!     ..
!  Purpose
!  =======
!
!  MDBDSVD22 deflates a broken arrow matrix. This routine is similar to 
!  MDBDSVD2. It first merges two parts of singular values into an ordered
!  array, and the left corresponding rows are also permuted. The permutation
!  information is stored in IDXQ as inputs. Then it deflates
!  the first broken arrow matrix. There are two ways in which deflation 
!  can occur:  when two or more singular values are close together or if 
!  there is a tiny entry in the Z vector.  For each such occurrence the 
!  order of the related secular equation problem is reduced by one.
!
!  The differences with MDBDSVD2 are that (1) the left and right singular 
!  vector matrices do not have sparse structure to use; (2) the matrices to 
!  deflate are always square. Therefore, it is much simpler than MDBDSVD2. 
!  Some index arrays are not required. For each call, there are two matrix 
!  copy operations. 
! 
!  The main task of this routine is to deflate singular values and move 
!  the deflated singular values and corresponding singular vectors to
!  the back part. MDBDSVD22 is called from MDBDSVD1 after calling MDBDSVD2.
!
!  Arguments
!  =========
!
!  N     (input) INTEGER
!         The dimension of this heavy broken arrow matrix.  N >= 1.
!
!  K      (output) INTEGER
!         Contains the dimension of the non-deflated matrix,
!         This is the order of the related secular equation. 1 <= K <=N.
!
!  SB     (input) INTEGER
!         The number of remaining appended rows, SB >= 1. By calling MDBDSVD
!         SB would gradually reduce to zero. 
!  
!  D      (input/output) DOUBLE PRECISION array, dimension(N)
!         On entry D contains the singular values of the two parts
!         to be combined, and D(SB) = ZERO.  
!         On exit D contains the trailing (N-K) updated singular values 
!         (those which were deflated) sorted into increasing order.
!
!  ZZ      (input/output) DOUBLE PRECISION array, dimension(N,*)
!         Z stores the appended rows which are stored in column form. 
!         On exit Z(SB:SB+K-1,SB) contains the updating row vector in the 
!         secular equation.
!
!  U      (input/output) DOUBLE PRECISION array, dimension(LDU,N)
!         On entry U contains the left singular vectors of original problem,
!         and it usually is a dense matrix. 
!         On exit U is permuted and the last N-K columns are the deflated singular
!         vectors. 
!
!  LDU    (input) INTEGER
!         The leading dimension of the array U.  LDU >= N.
!
!  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)
!         On entry VT**T contains the right singular vectors of the original problem,
!         and it usually is dense. On exit VT is permuted and the last M-K columns of
!         VT**T are the deflated singular vectors.
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT.  LDVT >= M.
!
!  DSIGMA (output) DOUBLE PRECISION array, dimension (N)
!         Contains a copy of the diagonal elements( K-1 singular values
!         and one zero) in the secular equation.
!
!  U2     (output) DOUBLE PRECISION array, dimension(LDU2,N)
!         Contains a copy of the first K-1 left singular vectors which
!         will be used by MDBDSVD4 in a matrix multiply (DGEMM) or 
!         MDLASD8 to solve for the new left singular vectors. 
!         U2 is arranged into two blocks. The first block contains a column 
!         with some nenzeros at positions from NL+1 to NL+SB; 
!         the second block is usually a dense block. 
!
!  LDU2   (input) INTEGER
!         The leading dimension of the array U2.  LDU2 >= N.
!
!  VT2    (output) DOUBLE PRECISION array, dimension(LDVT2,N)
!         VT2**T contains a copy of the first K right singular vectors
!         which will be used by MDBDSVD4 in a matrix multiply (DGEMM) or MDLASD8 to
!         solve for the new right singular vectors. VT2 is arranged into
!         two blocks. The first block contains a row that corresponds
!         to the special 0 diagonal element in SIGMA; the second block
!         is usually dense. 
!
!  LDVT2  (input) INTEGER
!         The leading dimension of the array VT2.  LDVT2 >= M.
!
!  IDXP   (workspace) INTEGER array dimension(N)
!         This will contain the permutation used to place deflated
!         values of D at the end of the array. On output IDXP(SB+1:K)
!         points to the nondeflated D-values and IDXP(K+1:N)
!         points to the deflated singular values.
!
!  IDXQ  (input/output) INTEGER array, dimension( N )
!         This contains the permutation which will reintegrate the
!         subproblem just solved back into sorted order, i.e.
!         D( IDXQ( I = 1, N ) ) will be in ascending order.
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!  Based on LAPACK routine DLASD2
!
!  Written by S.-G. Li, on July. 5th, 2013
!  =======================================
!
!     .. Parameters ..
      INTEGER            BSMLZ, SB1
      DOUBLE PRECISION   ZERO, ONE, TWO, EIGHT
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0, &
                          EIGHT = 8.0D+0, BSMLZ = 500 )
!     ..
!     .. Local Scalars ..
      INTEGER            CT, I, IDXI, IDXJ, IDXJP, J, JP, JPREV, K2
      DOUBLE PRECISION   C, EPS, HLFTOL, S, TAU, TOL, Z1
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   Z(N)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DLASET, DROT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MAXVAL
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      SB1  = SB + 1
      HFLAG = .TRUE. 
!
      IF( N.LT.1 ) THEN
         INFO = -1
      END IF
!
      IF( LDU.LT.N ) THEN
         INFO = -6
      ELSE IF( LDVT.LT.N ) THEN
         INFO = -8
      ELSE IF( LDU2.LT.N ) THEN
         INFO = -11
      ELSE IF( LDVT2.LT.N ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSVD22', -INFO )
         RETURN
      END IF
!
!     Assign the values of DSIGMA
!
!      CALL DLASET( 'A',N-SB,1,ZERO,ZERO,U2(1,SB), LDU2 )
      DO 60 I = SB1, N
         DSIGMA( I ) = D( IDXQ(I)+SB )  ! DSIGMA is in increasing order
         U2(I,SB) = ZZ( IDXQ(I)+SB,SB )
         IDXQ(I) = IDXQ(I) + SB
   60 CONTINUE
!
!     Store DSIGMA back to D
      DO I = SB1, N
         D( I ) = DSIGMA( I )  ! D is in increasing order starting from SB1
         ZZ( I,SB ) = U2( I,SB )
      END DO
!
!     Calculate the allowable deflation tolerance
!
      Z( SB1:N ) = U2( SB1:N,SB )
      EPS = DLAMCH( 'Epsilon' )
      TOL = MAXVAL( ABS( Z(SB1:N) ) )
      TOL = EIGHT*EPS*MAX( ABS( D( N ) ), TOL )
!
!     There are 2 kinds of deflation -- first a value in the z-vector
!     is small, second two (or more) singular values are very close
!     together (their difference is small).
!
!     If the value in the z-vector is small, we simply permute the
!     array so that the corresponding singular value is moved to the
!     end.
!
!     If two values in the D-vector are close, we perform a two-sided
!     rotation designed to make one of the corresponding z-vector
!     entries zero, and then permute the array so that the deflated
!     singular value is moved to the end.
!
!     If there are multiple singular values then the problem deflates.
!     Here the number of equal singular values are found.  As each equal
!     singular value is found, an elementary reflector is computed to
!     rotate the corresponding singular subspace so that the
!     corresponding components of Z are zero in this new basis.
!
      K = SB  
      K2 = N + 1
      DO 80 J = SB1, N
         IF( ABS( Z( J ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
            K2 = K2 - 1
            IDXP( K2 ) = J   ! point to current column
            IF( J.EQ.N ) &
               GO TO 120
         ELSE
            JPREV = J
            GO TO 90  ! If Z(J) is not small, check whether there are multiple svals.
         END IF
   80 CONTINUE
   90 CONTINUE
      J = JPREV
  100 CONTINUE
      J = J + 1
      IF( J.GT.N ) &
         GO TO 110
      IF( ABS( Z( J ) ).LE.TOL ) THEN  ! Next Z(J) is small
!
!        Deflate due to small z component.
!
         K2 = K2 - 1
         IDXP( K2 ) = J
      ELSE
!
!        Check if singular values are close enough to allow deflation.
!
         IF( ABS( D( J )-D( JPREV ) ).LE.TOL ) THEN 
!
!           Deflation is possible.
!
            S = Z( JPREV )
            C = Z( J )
!
!           Find sqrt(a**2+b**2) without overflow or
!           destructive underflow.
!
            TAU = DLAPY2( C,S )
            C = C / TAU
            S = -S / TAU
            Z( J ) = TAU
            Z( JPREV ) = ZERO
!
!           Apply back the Givens rotation to the other added rows, and 
!           the left and right singular vector matrices. The added rows are
!           stored in ZZ as block column.
!
!            CALL DROT( SB-1,ZZ( JPREV,1 ),M,ZZ( J,1 ),M,C,S )
!            CALL DROT( SB-1,U2( JPREV,1 ),LDU2,U2( J,1 ),LDU2,C,S )
            IDXJP = IDXQ( JPREV )
            IDXJ  = IDXQ( J )
            CALL DROT( SB-1,ZZ( IDXJP,1 ),M,ZZ( IDXJ,1 ),M,C,S )
            CALL DROT( N, U( 1,IDXJP ), 1, U( 1,IDXJ ), 1, C, S )
            CALL DROT( M, VT( IDXJP,1 ), LDVT, VT( IDXJ,1 ), LDVT,C,S )
            K2 = K2 - 1
            IDXP( K2 ) = JPREV
            JPREV = J
         ELSE   ! not 
            K = K + 1     ! the size of secular equation increase by one
            U2( K,SB ) = Z( JPREV )   ! store Z in the SB-th column of U2
            DSIGMA( K ) = D( JPREV )  ! store current D( J )
            IDXP( K ) = JPREV         ! point the kth non-deflated sval
            JPREV = J                 ! move on to next J
         END IF
      END IF
      GO TO 100   ! loop
  110 CONTINUE
!
!     Record the last singular value.
!
      K = K + 1
      U2( K, SB ) = Z( JPREV )
      DSIGMA( K ) = D( JPREV )
      IDXP( K ) = JPREV
!
  120 CONTINUE
!     
!     Sort the singular values and corresponding singular vectors into
!     DSIGMA, U2, and VT2 respectively. 
!     Modify the columns (rows) of U2 (VT2) starting from the second one.
!     The first row and column are modifed separately below. 
!
      DO 200 J = SB1, N
         JP = IDXP( J )
         DSIGMA( J ) = D( JP )
         IDXJ = IDXQ( IDXP(J) )  ! The original column in U
         CALL DCOPY( SB-1, ZZ(IDXJ,1), M, U2(J,1), LDU2 )
!        Copy N elements of the first N columns of U to U2
         CALL DCOPY( N, U( 1,IDXJ ), 1, U2( 1,J ), 1 )
!        Copy M elements of the first N rows of VT to VT2
         CALL DCOPY( M, VT( IDXJ,1 ), LDVT, VT2( J,1 ), LDVT2 )
200      CONTINUE
         CALL DLACPY( 'A',N-SB,SB,U2(SB1,1),LDU2,ZZ(SB1,1),M )
!
!     Determine DSIGMA(SB), DSIGMA(SB1) and Z(SB)
!
      DSIGMA( SB ) = ZERO
      HLFTOL = TOL / TWO
      IF( ABS( DSIGMA( SB1 ) ).LE.HLFTOL ) &
         DSIGMA( SB1 ) = HLFTOL
      Z1 = ZZ( SB,SB )
      IF( ABS( Z1 ).LE.TOL ) THEN
         ZZ( SB,SB ) = TOL   ! the SB-th column of ZZ stores Z
      ELSE
         ZZ( SB,SB ) = Z1
      END IF
!
!     Determine the SB-th column of U2, and the SB-th row of VT2
!
!      CALL DLASET( 'A', N-SB, 1, ZERO, ZERO, U2(1,SB), LDU2 )
      CALL DCOPY( N, U( 1,SB ), 1, U2( 1,SB ), 1 )  ! SB-th column of U2
      CALL DCOPY( M, VT( SB,1 ), LDVT, VT2( SB,1 ), LDVT2 ) ! SB-th row of VT2
!
      IF( K .GT. BSMLZ )  HFLAG = .FALSE. 
!
      IF( .NOT. HFLAG ) THEN
!
!        Use HSS techniques
         CALL DCOPY( N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 )
         CALL DLACPY( 'A', N, N-SB, U2( 1,SB1 ), LDU2, U( 1,SB1 ), &
                     LDU )
         CALL DLACPY( 'A', N-SB, M, VT2( SB1,1 ), LDVT2, VT( SB1,1 ), &
                     LDVT )
      ELSE IF( N.GT.K ) THEN
!
!     The deflated singular values and their corresponding vectors go
!     into the back of D, U, and V respectively.
!
         CALL DCOPY( N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 )
         CALL DLACPY( 'A', N, N-K, U2( 1,K+1 ), LDU2, U( 1,K+1 ), &
                     LDU )
         CALL DLACPY( 'A', N-K, M, VT2( K+1,1 ), LDVT2, VT( K+1,1 ), &
                     LDVT )
      END IF
!
      RETURN
!
!     End of MDBDSVD22
!
      END

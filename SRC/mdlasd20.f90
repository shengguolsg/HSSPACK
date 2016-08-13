  SUBROUTINE MDLASD20( N,K,D,Z,VT,LDVT,DSIGMA,VT2,LDVT2,IDXP,WW,HFLAG,INFO )
!
!  -- HSSPACK computational routine (version 1.0.0) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Defense Tech., China.
!     September 2013
!
        IMPLICIT NONE
!     ..Scalar Arguments..
      INTEGER            INFO, K, N, LDVT, LDVT2
      DOUBLE PRECISION   ALPHA, BETA
      LOGICAL            HFLAG
!     ..
!     ..Array Arguments..
      INTEGER            IDXP(*)
      DOUBLE PRECISION   D(*),DSIGMA(*),VT(LDVT, *),VT2(LDVT2, *),Z(*),WW(*)
!     ..
!  Purpose
!  =======
!
!  MDLASD20 deflates a blunt broken arrow matrix in the form of M = [diag(D); Z^T], 
!  of size (N+1)-by-N. It also permutes the right singular vector 
!  matrix correspondingly. There are two ways in which deflation can occur:  
!  when two or more singular values are close together or if there is a tiny entry in the
!  Z vector.  For each such occurrence the order of the related secular
!  equation problem is reduced by one.
!
!  MDLASD21 is designed for updating SVD problem, and is called from SvdHssUpdatAll.
!  A related routine MDLASD20 only permutes the right singular vector matrix, which is 
!  called from SvdHssUpdat.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The dimension of D, the column dimension of M. We always assume
!         M is column full rank. 
!
!  K      (output) INTEGER
!         Contains the dimension of the non-deflated matrix,
!         This is the order of the related secular equation. 1 <= K <=N.
!
!  D      (input/output) DOUBLE PRECISION array, dimension(N)
!         On entry D contains the singular values of the two submatrices
!         to be combined.  On exit D contains the trailing (N-K) updated
!         singular values (those which were deflated) sorted into
!         increasing order.
!
!  Z      (output) DOUBLE PRECISION array, dimension(N)
!         On exit Z contains the updating row vector in the secular
!         equation.
!
!  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,N)
!         On entry VT**T contains the right singular vectors.
!         On exit VT**T contains the trailing (N-K) updated right singular
!         vectors (those which were deflated) in its last N-K columns.
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT.  LDVT >= N.
!
!  DSIGMA (output) DOUBLE PRECISION array, dimension (N)
!         Contains a copy of the diagonal elements (K-1 singular values
!         and one zero) in the secular equation.
!
!  VT2    (output) DOUBLE PRECISION array, dimension(LDVT2,N)
!         VT2**T contains a copy of the first K right singular vectors
!         which will be used by SvdHssUpdatAll in an HSS matrix multiply to
!         solve for the new right singular vectors. 
!
!  LDVT2  (input) INTEGER
!         The leading dimension of the array VT2.  LDVT2 >= N.
!
!  IDXP   (workspace) INTEGER array dimension(N)
!         This will contain the permutation used to place deflated
!         values of D at the end of the array. On output IDXP(2:K)
!         points to the nondeflated D-values and IDXP(K+1:N)
!         points to the deflated singular values.
!
!  WW     (workspace) DOUBLE PRECISION array, dimension(N)
!         This is used to permute the appended vector Z. 
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!  Based on dlasd2 in LAPACK 3.4.0
!
!  Written by S.-G. Li, on April 25th, 2013
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, EIGHT
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0, &
                         EIGHT = 8.0D+0 )
!     ..
!     .. Local Arrays ..
!     ..
!     .. Local Scalars ..
      INTEGER            CT, I, IDXI, IDXJ, IDXJP, J, JP, JPREV, K2, &
                         BSMLZ
      DOUBLE PRECISION   C, EPS, HLFTOL, S, TAU, TOL, Z1
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DLASET, DROT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      BSMLZ = 10
!
      IF( N.LT.1 ) THEN
         INFO = -1
      ELSE IF( LDVT.LT.N ) THEN
         INFO = -8
      ELSE IF( LDVT2.LT.N ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDLASD20', -INFO )
         RETURN
      END IF
!
!     Calculate the allowable deflation tolerance
!
      EPS = DLAMCH( 'Epsilon' )
      TOL = MAX( ABS( ALPHA ), ABS( BETA ) )
      TOL = MAX( ABS( D(N) ), TOL )
      TOL = EIGHT*EPS*MAX( ABS( D( 1 ) ), TOL )
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
      K = 0
      K2 = N + 1
      DO 80 J = 1, N
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
            TAU = DLAPY2( C, S )
            C = C / TAU
            S = -S / TAU
            Z( J ) = TAU
            Z( JPREV ) = ZERO
!
!           Apply back the Givens rotation to the left and right
!           singular vector matrices.
!
!            CALL DROT( N, U( 1, JPREV ), 1, U( 1, J ), 1, C, S )
            CALL DROT( N, VT( JPREV,1 ), LDVT, VT( J,1 ), LDVT, C, &
                       S )
            K2 = K2 - 1
            IDXP( K2 ) = JPREV
            JPREV = J
         ELSE   ! not 
            K = K + 1     ! the size of secular equation increase by one
            WW( K ) = Z( JPREV )   ! store current Z(J)
            DSIGMA( K ) = D( JPREV )  ! store current D(J)
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
      WW( K ) = Z( JPREV )
      DSIGMA( K ) = D( JPREV )
      IDXP( K ) = JPREV
!
  120 CONTINUE
!
!     Move the rest of the updating row to Z. Here we assume D are
!     in increasing order.
      DO J = 1, K
         Z(J) = WW(J)
      END DO
!
!     Sort the singular values and corresponding singular vectors into
!     DSIGMA, U2, and VT2 respectively.  The singular values/vectors
!     which were not deflated go into the first K slots of DSIGMA, U2,
!     and VT2 respectively, while those which were deflated go into the
!     last N - K slots, except that the first column/row will be treated
!     separately.
!
      DO J = 1, N
         JP = IDXP( J )
         WW( J ) = D( JP )
!         CALL DCOPY( N, U( 1, JP ), 1, U2( 1, J ), 1 )
         CALL DCOPY( N, VT( JP,1 ), LDVT, VT2( J,1 ), LDVT2 )
      END DO
!
!     Copy the old singular values to DSIGMA, assuming D are in increasing order
      DO J = 1, K
         D(J) = WW(J)
      END DO
!
!     Determine DSIGMA(K) and Z(K)
!
      HLFTOL = TOL / TWO
      IF( ABS( DSIGMA( 1 ) ).LE.HLFTOL ) &
         DSIGMA( K ) = HLFTOL
      IF( ABS( Z(1) ).LE.TOL ) THEN
         Z( 1 ) = TOL
      END IF
!
      IF( K .GT. BSMLZ )  HFLAG = .FALSE. 
!      write(*,*) 'We entered.', N-K
!
      IF( .NOT. HFLAG ) THEN
!
!     The deflated singular values and all singular vectors go
!     into the back of D, U, and V respectively.
!
         CALL DCOPY( N-K, WW( K+1 ), 1, D( K+1 ), 1 )
!        CALL DLACPY( 'A', N, N, U2( 1, 1 ), LDU2, U( 1, 1 ), &
!                  LDU )
         CALL DLACPY( 'A', N, N, VT2( 1,1 ), LDVT2, VT( 1,1 ), &
                     LDVT )
!
      ELSE IF( N.GT.K ) THEN
!
!     The deflated singular values and their corresponding vectors go
!     into the back of D, U, and V respectively.
!
         CALL DCOPY( N-K, WW( K+1 ), 1, D( K+1 ), 1 )
!         CALL DLACPY( 'A', N, N-K, U2( 1, K+1 ), LDU2, U( 1, K+1 ), &
!                      LDU )
         CALL DLACPY( 'A',N-K,N,VT2( K+1,1 ),LDVT2,VT( K+1,1 ), &
                      LDVT )
      END IF
!
!
      RETURN
!
!     End of MDLASD20
!
      END 

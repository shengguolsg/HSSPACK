!
      SUBROUTINE MDSBED4( K,N,SB,D,Q,LDQ,RHO,DLAMDA,Q2,INDXQ,W,S,Z,Z2,INFO )
! 
        IMPLICIT NONE
!       ..
!       .. Scalar Arguments ..
        INTEGER            INFO, K, LDQ, N, SB
        DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
        INTEGER            INDXQ( * )
        DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ,* ), Q2( * ), &
                           S( * ), W( * ), Z( N,* ), Z2( N,* )
!       ..
!
! Purpose
! ========
!
! MDSBED4 finds the roots of the secular equation, as defined by the
! values in D, W, and RHO, between 1 and K.  It makes the
! appropriate calls to DLAED4 and then updates the eigenvectors by
! multiplying the matrix of eigenvectors of the pair of eigensystems
! being combined by the matrix of eigenvectors of the K-by-K system
! which is solved here.
!
! It is similar to MDSBED3, and the difference is that this routine does not
! explore any sparse structure in Q. The eigenvectors are simply updated via
! matrix-matrix multiplication. The other appended vectors are also updated. 
!
! MDSBED4 is called from MDSBED1, and is designed for MDSBED22.  INDXQ is not
! used in this version now, and it can be removed. 
!
!  Arguments
! ===========
!
! \param[in] K
!          K is INTEGER
!          The number of terms in the rational function to be solved by
!          DLAED4.  K >= 0. 
!
! \param[in] N
!          N is INTEGER
!          The number of rows and columns in the Q matrix.
!          N >= K ( deflation may result in N>K ).
!
! \param[in] N1
!          N1 is INTEGER
!          The location of the last eigenvalue in the leading submatrix.
!          min(1,N) <= N1 <= N/2.
!
! \param[in] SB
!          N1 is INTEGER
!          The number of columns of Z, and it relates to the semi-bandwith of
!          the original symmetric band matrix.
!
! \param[out] D
!          D is DOUBLE PRECISION array, dimension ( N )
!          D(I) contains the updated eigenvalues for
!          1 <= I <= K.
!
! \param[out] Q
!          Q is DOUBLE PRECISION array, dimension ( LDQ,N )
!          Initially the first K columns are used as workspace.
!          On output the columns 1 to K contain
!          the updated eigenvectors.
!
! \param[in] LDQ
!          LDQ is INTEGER
!          The leading dimension of the array Q.  LDQ >= max(1,N).
!
! \param[in] RHO
!          RHO is DOUBLE PRECISION
!          The value of the parameter in the rank one update equation.
!          RHO >= 0 required.
!
! \param[in,out] DLAMDA
!          DLAMDA is DOUBLE PRECISION array, dimension (K)
!          The first K elements of this array contain the old roots
!          of the deflated updating problem.  These are the poles
!          of the secular equation. May be changed on output by
!          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
!          Cray-2, or Cray C-90, as described above.
!
! \param[in] Q2
!          Q2 is DOUBLE PRECISION array, dimension (LDQ2, N)
!          The first K columns of this matrix contain the non-deflated
!          eigenvectors for the split problem.
!
! \param[in] INDXQ
!          INDXQ is INTEGER array, dimension (N)
!          The permutation used to arrange the columns of the deflated
!          Q matrix into three groups (see DLAED2).
!          The rows of the eigenvectors found by DLAED4 must be likewise
!          permuted before the matrix multiply can take place.
!
! \param[in] CTOT
!          CTOT is INTEGER array, dimension (4)
!          A count of the total number of the various types of columns
!          in Q, as described in INDXQ.  The fourth column type is any
!          column which has been deflated.
!
! \param[in,out] W
!          W is DOUBLE PRECISION array, dimension (K)
!          The first K elements of this array contain the components
!          of the deflation-adjusted updating vector. Destroyed on
!          output.
!
! \param[out] S
!          S is DOUBLE PRECISION array, dimension MAX( (N1 + 1)*K, K*K ) 
!          Will contain the eigenvectors of the repaired matrix which
!          will be multiplied by the previously accumulated eigenvectors
!          to update the system.
!
! \param[in] Z
!          Z is DOUBLE PRECISION array, dimension (N,SB)
!          The appended block row of after deflation
!
! \param[out] Z2
!          Z2 is DOUBLE PRECISION array, dimension (N,SB)
!          The returned appended block row of after updation.
!
! \param[out] INFO
!          INFO is INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an eigenvalue did not converge
!
!  Authors:
!  ========
!
! \author Nati. Univ. of Def. Tech. 
! \date October 2014 
! 
! Contributors 
!  ============
! 
!  Shengguo Li, NUDT, China, 2014
! 
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, II, J
      DOUBLE PRECISION   TEMP
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLAED4, DLASET, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.K ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDSBED3', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.0 )  RETURN
!
!     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract ( Cray XMP, Cray YMP, Cray C 90 and Cray 2 ).
!     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DLAMDA(I) if it is 1; this makes the subsequent
!     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DLAMDA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DLAMDA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO I = 1, K
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
      END DO
!
      DO J = 1, K
         CALL DLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) &
            GO TO 120
      END DO
!
      IF( K.EQ.1 )  &
         GO TO 110
      IF( K.EQ.2 ) THEN
         DO J = 1, K
            W( 1 ) = Q( 1, J )
            W( 2 ) = Q( 2, J )
!            II = INDXQ( 1 )
            II = 1
            Q( 1, J ) = W( II )
!            II = INDXQ( 2 )
            II = 2
            Q( 2, J ) = W( II )
         END DO
         GO TO 110
      END IF
!
!     Compute updated W.
!
      CALL DCOPY( K, W, 1, S, 1 )
!
!     Initialize W(I) = Q(I,I)
!
      CALL DCOPY( K, Q, LDQ+1, W, 1 )
      DO J = 1, K
         DO I = 1, J - 1
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
         END DO
         DO I = J + 1, K
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
         END DO
      END DO
      DO I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I ) )
      END DO
!
!     Compute eigenvectors of the modified rank-1 modification.
!
      DO J = 1, K
         DO I = 1, K
            S( I ) = W( I ) / Q( I, J )
         END DO
         TEMP = DNRM2( K, S, 1 )
         DO I = 1, K
!            II = INDXQ( I )
            II = I
            Q( I,J ) = S( II ) / TEMP
         END DO
      END DO
!
!     Compute the updated eigenvectors.
!
  110 CONTINUE
!
!     Update the remaining added vectors
      IF( SB .GT. 1 ) THEN
         CALL DGEMM( 'T','N', K, SB, K, ONE, Q, LDQ, Z, N, ZERO, Z2, N  )
         IF( K .LT. N ) THEN
            CALL DLACPY( 'A',N-K,SB,Z(K+1,1),N,Z2(K+1,1),N )
         END IF
      END IF
!
      CALL DLACPY( 'A', K, K, Q, LDQ, S, K )
      CALL DGEMM( 'N','N',N,K,K,ONE,Q2,N,S,K,ZERO,Q,LDQ )
!     
!
  120 CONTINUE
      RETURN
!
!     End of MDSBED4
!
    END SUBROUTINE MDSBED4

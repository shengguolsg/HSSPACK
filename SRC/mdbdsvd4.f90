    SUBROUTINE MDBDSVD4( M,N,SB,K,D,Q,LDQ, DSIGMA, U, LDU, U2, &
                       LDU2,VT,LDVT,VT2,LDVT2,ZZ,LDZZ,INFO )
!
      Use basicmm
      IMPLICIT NONE
!
!  -- HSSPACK computational routine (version 1.0.0) --
!  -- HSSPACK is a software package provided by 
!     Nati. Univ. of Def. Tech., China, January 2013
!
!     .. Scalar Arguments ..
      INTEGER            M,N,SB,INFO,K,LDQ,LDU,LDU2,LDVT,LDVT2,LDZZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DSIGMA(*), Q( LDQ,*), U( LDU,*), &
                         U2( LDU2,* ), VT( LDVT,* ), VT2( LDVT2,*),&
                         ZZ( LDZZ,* )
!     ..
!
!  Purpose
!  =======
!  MDBDSVD4 finds all the square roots of the roots of the secular
!  equation, as defined by the values in D and Z.  It makes the
!  appropriate calls to DLASD4 and then updates the singular
!  vectors by matrix multiplications.
!
!  It is similar to MDBDSVD3, and the difference is that 
!  this routine does not consider the sparse structure of U and VT 
!  since both are essentially dense. The singular vectors
!  are updated via plain matrix-matrix multiplication. 
! 
!  MDBDSVD4 is called from MDBDSVD1, and is designed for MDBDSVD22. 
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The size of row dimension of the merged subproblem
! 
!  SB     (input) INTEGER
!         The number of current added rows
! 
!  K      (input) INTEGER
!         The size of the secular equation, 1 =< K = < N.
!
!  D      (output) DOUBLE PRECISION array, dimension(K)
!         On exit the square roots of the roots of the secular equation,
!         in ascending order.
!
!  Q      (workspace) DOUBLE PRECISION array,
!                     dimension at least (LDQ,K).
!
!  LDQ    (input) INTEGER
!         The leading dimension of the array Q.  LDQ >= K.
!
!  DSIGMA (input) DOUBLE PRECISION array, dimension(N)
!         The entries from SB+1 to SB+K are the old roots of the deflated 
!         updating problem. These are the poles of the secular equation.
!
!  U      (output) DOUBLE PRECISION array, dimension (LDU, N)
!         The last N - K columns of this matrix contain the deflated
!         left singular vectors, and the columns from SB to (SB+K-1) contain 
!         the non-deflated left singular vectors. 
!
!  LDU    (input) INTEGER
!         The leading dimension of the array U.  LDU >= N.
!
!  U2     (input/output) DOUBLE PRECISION array, dimension (LDU2, N)
!         The K columns from SB to (SB+K-1) contain the non-deflated
!         left singular vectors for the split problem.
!
!  LDU2   (input) INTEGER
!         The leading dimension of the array U2.  LDU2 >= N.
!
!  VT     (output) DOUBLE PRECISION array, dimension (LDVT, M)
!         The last M - K columns of VT**T contain the deflated
!         right singular vectors, and the first K columns of VT**T 
!         contains the non-deflated right singular vectors. 
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT.  LDVT >= N.
!
!  VT2    (input/output) DOUBLE PRECISION array, dimension (LDVT2, N)
!         The first columns from SB to SB+K-1 of VT2**T contain the non-deflated
!         right singular vectors for the split problem.
!
!  LDVT2  (input) INTEGER
!         The leading dimension of the array VT2.  LDVT2 >= N.
!
!  ZZ     (input/output) DOUBLE PRECISION array, dimension(LDZZ,SB)
!         The elements from SB to SB+K-1 of SB-th column contain the components
!         of the deflation-adjusted updating row vector.
!
!  LDZZ   (input) INTEGER
!         The leading dimension of block vector ZZ, M
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!  Based on LAPACK routine DLASD3
!
!  Written by S.-G. Li, on July. 6th, 2013
!  =======================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, NEGONE
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0, NEGONE = -1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            CTEMP,I,J,JC,KTEMP,SBM1,JB,IB,SB1
      DOUBLE PRECISION   RHO,TEMP
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION Z(K)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLASCL, MDLASD4, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      SBM1 = SB -1 
      SB1  = SB + 1
!
      IF( K.LT.1 ) THEN
         INFO = -1
      END IF
!
      IF( LDQ.LT.K ) THEN
         INFO = -4
      ELSE IF( LDU.LT.K ) THEN
         INFO = -7
      ELSE IF( LDU2.LT.K ) THEN
         INFO = -9
      ELSE IF( LDVT.LT.K ) THEN
         INFO = -11
      ELSE IF( LDVT2.LT.K ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSVD4', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      Z(1:K) = ZZ(SB:SBM1+K,SB)
      IF( K.EQ.1 ) THEN
         D( SB ) = ABS( Z( 1 ) )
         CALL DCOPY( M, VT2( SB, 1 ), LDVT2, VT( SB, 1 ), LDVT )
         IF( Z( 1 ).GT.ZERO ) THEN
            CALL DCOPY( N, U2( 1,SB ), 1, U( 1,SB ), 1 )
         ELSE
            DO 10 I = 1, N
               U( I,SB ) = -U2( I,SB )
   10       CONTINUE
         END IF
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
!     2*DSIGMA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO 20 I = SB, K+SBM1
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   20 CONTINUE
!
!     Keep a copy of Z.
!
      CALL DCOPY( K, Z, 1, Q, 1 )
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Find the new singular values. The columns from SB to SBM1+K of 
!     U are used as workspace. So does the columns of VT. It would be 
!     K-by-K matrices. Since the first SB rows of VT are useful, the 
!     row of VT starts from the SB-th row, and hope the first SB row 
!     has been copied to VT2. So does U, its first SB-1 columns are 
!     not referenced. 
!     Ignore the shift SB, U(i,j) stores DSIGMA(i)-D(j), and VT(i,j)
!     stores DSIGMA(i)+D(j).
!
      DO 30 J = SB, K+SBM1
         CALL DLASD4( K, J-SBM1, DSIGMA(SB), Z, U( 1, J ), RHO, D( J ), &
                     VT( SB, J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) THEN
            write(*,*) 'DLASD4 exits in d4', info, 'J', J, SB
            RETURN
         END IF
   30 CONTINUE
!
!     Compute updated Z.
!
      DO 60 I = 1, K
         IB = I + SBM1
         Z( I ) = U( I, K+SBM1 )*VT( IB, K+SBM1 )
         DO 40 J = 1, I - 1
            JB = J + SBM1
            Z( I ) = Z( I )*( U( I,JB )*VT( IB,JB ) / &
                    ( DSIGMA( IB )-DSIGMA( JB ) ) / &
                    ( DSIGMA( IB )+DSIGMA( JB ) ) )
   40    CONTINUE
         DO 50 J = I, K - 1
            JB = J + SBM1
            Z( I ) = Z( I )*( U( I,JB )*VT( IB,JB ) / &
                    ( DSIGMA( IB )-DSIGMA( JB+1 ) ) / &
                    ( DSIGMA( IB )+DSIGMA( JB+1 ) ) )
   50    CONTINUE
         Z( I ) = SIGN( SQRT( ABS( Z( I ) ) ), Q( I, 1 ) )
   60 CONTINUE
!
!     Compute left singular vectors of the modified diagonal matrix,
!     and store related information for the right singular vectors.
!     U and VT are first used as workspace
!
      DO 90 I = 1, K
         IB = I + SBM1
         VT( SB,IB ) = Z( 1 ) / U( 1,IB ) / VT( SB,IB )
         U( 1,IB ) = NEGONE
         DO 70 J = 2, K
            JB = J+SBM1
            VT( JB,IB ) = Z( J ) / U( J,IB ) / VT( JB,IB )
            U( J,IB ) = DSIGMA( JB )*VT( JB,IB )
   70    CONTINUE
         TEMP = DNRM2( K, U( 1,IB ), 1 )
         Q( 1, I ) = U( 1,IB ) / TEMP
         DO 80 J = 2, K
            Q( J, I ) = U( J,IB ) / TEMP
   80    CONTINUE
   90 CONTINUE
!
!     Update the left singular vector matrix.
!
      CALL DGEMM( 'N','N',N,K,K,ONE,U2(1,SB),LDU2,Q,LDQ,&
           ZERO,U(1,SB),LDU )
!
!     Generate the right singular vectors.
!
      DO 120 I = 1, K
         IB = I + SBM1
         TEMP = DNRM2( K, VT( SB, IB ), 1 )
         Q( I, 1 ) = VT( SB, IB ) / TEMP
         DO 110 J = 2, K
            JC = J + SBM1
            Q( I,J ) = VT( JC,IB ) / TEMP
  110    CONTINUE
  120 CONTINUE
!
!     Update the remaining added vectors, use the first SB-1 columns of
!     U2 as workspace
      IF( SB .GT. 1 ) THEN
         CALL DLACPY('A',K,SB-1,ZZ(SB,1),LDZZ,U2(SB,1),LDU2 )
         CALL DGEMM( 'N','N',K,SB-1,K,ONE,Q,LDQ,U2(SB,1),LDU2,&
              ZERO,ZZ(SB,1),LDZZ )
      END IF
!
!     Update the right singular vector matrix.
!
      CALL DGEMM( 'N','N',K,M,K,ONE,Q,LDQ,VT2(SB,1),LDVT2, & 
           ZERO,VT(SB,1),LDVT )
!

      RETURN
!
!     End of MDBDSVD4
!
      END

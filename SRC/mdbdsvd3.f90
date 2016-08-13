    SUBROUTINE MDBDSVD3( NL,NR,SB,SQRE,K,D,Q,LDQ,DSIGMA,U,LDU,U2, &
           LDU2,VT,LDVT,VT2,LDVT2,IDXC,CTOT,ZZ,LDZZ,INFO )
!
      Use basicmm
      IMPLICIT NONE
!
!  -- HSSPACK computational routine (version 0.0.1) --
!  -- HSSPACK is a software package provided by 
!     Nati. Univ. of Def. Tech., China, January 2013
!
!     .. Scalar Arguments ..
      INTEGER            INFO,K,LDQ,LDU,LDU2,LDVT,LDVT2,NL,NR,&
                         SQRE,SB,LDZZ
!     ..
!     .. Array Arguments ..
      INTEGER            CTOT( * ), IDXC( * )
      DOUBLE PRECISION   D( * ), DSIGMA(*), Q( LDQ,*), U( LDU,*), &
                         U2( LDU2,* ), VT( LDVT,* ), VT2( LDVT2,*), &
                         ZZ( LDZZ,* )
!  ..
!  Purpose
!  =======
!  MDBDSVD3 finds all the square roots of the roots of the secular
!  equation, as defined by the values in D and Z.  It makes the
!  appropriate calls to DLASD4 and then updates the singular
!  vectors by matrix multiplication.
!
!  It is similar to DLASD3, but there are some differences:
!  1) The first SB columns of U2 are special, each of which has 
!     SB non-zeros. In DLASQ3 only the first column is special 
!     which only has one nonzero. 
!  2) ZZ is a block column which has SB columns, and after dealing with
!     the last columns should be multiplied by the right singular vector
!     matrix. 
!
!  MDBDSVD3 is called from MDBDSVD1, and designed for MDBDSVD2. 
!
!  Arguments
!  =========
!  NL     (input) INTEGER
!         The row dimension of the upper block.  NL >= 1.
!
!  NR     (input) INTEGER
!         The row dimension of the lower block.  NR >= 1.
!
!  SQRE   (input) INTEGER
!         = 0 : the lower block is an NR-by-NR square matrix.
!         = SB: the lower block is an NR-by-(NR+SB) rectangular matrix.
!
!         The heavy broken arrow matrix has N = NL + NR + SB rows and
!         M = N + SQRE >= N columns.
!
! SB      (input) INTEGER
!         The number of nonzeros in the SB-th column of U2
!
!  K      (input) INTEGER
!         The size of the secular equation, 1 =< K = < N.
!
!  D      (output) DOUBLE PRECISION array, dimension(K+SB-1)
!         On exit the square roots of the roots of the secular equation,
!         in ascending order, and the updated values are stored at the 
!         places from SB to SB+K-1. The first SB-1 entries are zeros. 
!
!  Q      (workspace) DOUBLE PRECISION array,
!                     dimension at least (LDQ,K).
!
!  LDQ    (input) INTEGER
!         The leading dimension of the array Q.  LDQ >= K.
!
!  DSIGMA (input) DOUBLE PRECISION array, dimension(K+SB-1)
!         The elements of this array from SB to SB+K-1 contain the old roots
!         of the deflated updating problem.  These are the poles
!         of the secular equation. The SB-th value of DSIGMA is zero, and so
!         does the first SB-1 entries of DSIGMA.
!
!  U      (output) DOUBLE PRECISION array, dimension (LDU, N)
!         The last N - K columns of this matrix contain the deflated left 
!         singular vectors, and the first K columns contain the non-deflated
!         left singular vectors. The columns from SB to SB+K-1 will be updated
!         via matrix-matrix multiplcations. 
!
!  LDU    (input) INTEGER
!         The leading dimension of the array U.  LDU >= N.
!
!  U2     (input/output) DOUBLE PRECISION array, dimension (LDU2, N)
!         The first K columns of this matrix contain the non-deflated
!         left singular vectors for the split problem. Each column of the 
!         first SB only has SB nonzero entries from NL+1 to NL+SB. The 
!         columns from SB+1 to K are divided into three types: 1, 2 and 3. 
!
!  LDU2   (input) INTEGER
!         The leading dimension of the array U2.  LDU2 >= N.
!
!  VT     (output) DOUBLE PRECISION array, dimension (LDVT, M)
!         The last M - K columns of VT**T contain the deflated
!         right singular vectors, and the columns from SB1 to
!         SB+K-1 of VT**T contains the non-deflated right singular 
!         vectors. 
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT.  LDVT >= N.
!
!  VT2    (input/output) DOUBLE PRECISION array, dimension (LDVT2, N)
!         The first K rows of VT2 contain the non-deflated
!         right singular vectors for the split problem. The first SB
!         rows either only have nonzeros at the first NL+SB columns or 
!         are dense rows. The rows from SB+1 to K are divided into three
!         types: 1, 2, and 3. 
!
!  LDVT2  (input) INTEGER
!         The leading dimension of the array VT2.  LDVT2 >= N.
!
!  IDXC   (input) INTEGER array, dimension ( N )
!         The permutation used to arrange the columns of U (and rows of
!         VT) into three groups:  the first group contains non-zero
!         entries only at and above (or before) NL+SB; the second
!         contains non-zero entries only at and below (or after) NL+SB1;
!         and the third is dense. The first SB columns of U and the first
!         SB rows of VT are treated separately, however.
!
!         The rows of the singular vectors found by DLASD4 must be likewise 
!         permuted before the matrix multiplies can take place.
!
!  CTOT   (input) INTEGER array, dimension ( 4 )
!         A count of the total number of the various types of columns
!         in U (or rows in VT), as described in IDXC. The fourth column
!         type is any column which has been deflated.
!
!  ZZ     (input) DOUBLE PRECISION array, dimension( M,SB )
!         The last column contain the components of the deflation-adjusted 
!         updating row vector. The elements of ZZ(SB:SB+K-1,SB) will 
!         be copied to a new array Z with dimension( K ). After the computation
!         of right singular vector matrix, rows from SB to SB+K-1 of the first
!         SB-1 columns of ZZ should be updated. 
!         
!  LDZZ   (input) INTEGER, leading dimension of ZZ
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!  Based on LAPACK routine dlasd3
!
!  The difference is that each of the first SB columns of U2 has SB nonzeros, SB > 1
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, NEGONE
      PARAMETER          ( ONE=1.0D+0, ZERO=0.0D+0, NEGONE=-1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            CTEMP,SBM1,I,J,JC,SB1,KTEMP,M,N,NLP1,NLPB,&
                         NRP1,NLPB1,IB,JB
      DOUBLE PRECISION   RHO, TEMP, Z(K)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLASCL, DLASD4, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( NL.LT.1 ) THEN
         INFO = -1
      ELSE IF( NR.LT.1 ) THEN
         INFO = -2
      ELSE IF( ( SQRE.NE.SB ) .AND. ( SQRE.NE.0 ) ) THEN
         INFO = -3
      ELSE IF( SB.LT.2 ) THEN
         INFO = -4
      END IF
!
      N = NL + NR + SB
      M = N + SQRE
      NLP1 = NL + 1
      NLPB = NL + SB
      NLPB1 = NLPB+1
      SBM1 = SB - 1
      SB1  = SB + 1
!
      IF( ( K.LT.1 ) .OR. ( K.GT.N ) ) THEN
         INFO = -5
      ELSE IF( LDQ.LT.K ) THEN
         INFO = -8
      ELSE IF( LDU.LT.N ) THEN
         INFO = -11
      ELSE IF( LDU2.LT.N ) THEN
         INFO = -13
      ELSE IF( LDVT.LT.M ) THEN
         INFO = -15
      ELSE IF( LDVT2.LT.M ) THEN
         INFO = -17
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSVD3', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      Z( 1:K ) = ZZ( SB:SBM1+K,SB )
      IF( K.EQ.1 ) THEN
         D( SB ) = ABS( Z( 1 ) )
         CALL DCOPY( M, VT2( SB, 1 ), LDVT2, VT( SB, 1 ), LDVT )
         IF( Z( 1 ).GT.ZERO ) THEN
            CALL DCOPY( N, U2( 1, SB ), 1, U( 1, SB ), 1 )
         ELSE
            DO I = 1, N
               U( I, SB ) = -U2( I, SB )
            END DO
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
      DO I = SB, K+SBM1
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
      END DO
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
         CALL DLASD4( K, J-SBM1, DSIGMA(SB), Z, U( 1,J ), RHO, D(J), &
                     VT( SB,J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) THEN
            write(*,*) 'MDBDSVD3 exits, error'
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
            JB = J+SBM1
            JC = IDXC( JB )-SBM1   ! ? ?
            Q( J, I ) = U( JC,IB ) / TEMP
   80    CONTINUE
   90 CONTINUE
!
!     Update the left singular vector matrix.
!
      IF( K.EQ.2 ) THEN
         CALL DGEMM( 'N','N',N,K,K,ONE,U2(1,SB),LDU2,Q,LDQ,ZERO,U(1,SB), &
                    LDU )
         GO TO 100
      END IF
      IF( CTOT( 1 ).GT.0 ) THEN     ! TYPE 1: start from (SB1)-th column
         CALL DGEMM( 'N', 'N', NL, K, CTOT( 1 ), ONE, U2(1,SB1),LDU2, &
                    Q( 2, 1 ), LDQ, ZERO, U( 1, SB ), LDU )
         IF( CTOT( 3 ).GT.0 ) THEN  ! TYPE 3 is updated with TYPE 1
            KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
            CALL DGEMM( 'N','N',NL,K,CTOT(3),ONE,U2(1,KTEMP+SBM1),&
                      LDU2,Q(KTEMP,1),LDQ,ONE,U(1,SB),LDU )
         END IF
      ELSE IF( CTOT( 3 ).GT.0 ) THEN  ! TYPE 3 is updated lonely
         KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
         CALL DGEMM( 'N','N',NL,K,CTOT(3),ONE,U2(1,KTEMP+SBM1), &
                   LDU2, Q( KTEMP,1 ), LDQ, ZERO, U( 1,SB ), LDU )
      ELSE  ! TYPE 1 and 3 do not exist
         CALL DLACPY( 'F', NL, K, U2, LDU2, U, LDU )
      END IF
!     The rows from NLP1 to NL+SB
      CALL DGEMM('N','N',SB,K,1,ONE,U2(NLP1,SB),LDU2,Q(1,1),LDQ,ZERO, &
           U(NLP1,SB), LDU )
      KTEMP = 2 + CTOT( 1 )
      CTEMP = CTOT( 2 ) + CTOT( 3 )
      CALL DGEMM( 'N','N',NR,K,CTEMP,ONE,U2(NLPB1,KTEMP+SBM1),LDU2, &
                 Q(KTEMP,1),LDQ,ZERO,U(NLPB1,SB),LDU )
!
!     Generate the right singular vectors.
!
  100 CONTINUE
      DO I = 1, K
         IB = I + SBM1
         TEMP = DNRM2( K, VT( SB, IB ), 1 )
         Q( I, 1 ) = VT( SB, IB ) / TEMP
         DO J = 2, K
            JC = IDXC( J+SBM1 )   ! ? ? 
            Q( I, J ) = VT( JC, IB ) / TEMP
         END DO
      END DO
!
!     Update the remaining added vectors, use the first SB-1 columns of
!     U2 as workspace
      CALL DLACPY('A',K,SB-1,ZZ(SB,1),LDZZ,U2(SB,1),LDU2 )
      CALL DGEMM( 'N','N',K,SB-1,K,ONE,Q,LDQ,U2(SB,1),LDU2,&
           ZERO,ZZ(SB,1),LDZZ )
!
!     Update the right singular vector matrix
!
      IF( K.EQ.2 ) THEN
         CALL DGEMM( 'N', 'N', K, M, K, ONE, Q, LDQ, VT2(SB,1), LDVT2, ZERO, &
                    VT(SB,1), LDVT )
         RETURN
      END IF
      KTEMP = 1 + CTOT( 1 )  ! The first row and TYPE 1
      CALL DGEMM( 'N','N',K,NLPB,KTEMP,ONE,Q( 1,1 ),LDQ, &
                  VT2(SB,1), LDVT2, ZERO, VT( SB,1 ), LDVT )
      KTEMP = SB + CTOT( 1 ) + CTOT( 2 ) + 1
      IF( KTEMP.LE.LDVT2 ) &  ! The front part of TYPE 3
        CALL DGEMM( 'N', 'N', K, NLPB, CTOT(3), ONE, Q( 1,KTEMP-SBM1), &
                   LDQ, VT2( KTEMP,1 ), LDVT2, ONE, VT( SB,1 ), &
                   LDVT )
!
      KTEMP = CTOT( 1 ) + 1
      NRP1 = NR + SQRE
      IF( KTEMP.GT.1 ) THEN
         DO 130 I = 1, K
            Q( I, KTEMP ) = Q( I, 1 )    ! move the first column of Q to the KTEMP-th column
  130    CONTINUE
         DO 140 I = NLPB1, M
            VT2( KTEMP+SBM1,I ) = VT2( SB,I ) ! move the back part of SB-th row of VT2 to the KTEMP-th row
  140    CONTINUE
      END IF
      CTEMP = 1 + CTOT( 2 ) + CTOT( 3 )  ! TYPE 2 and 3
      CALL DGEMM( 'N', 'N', K, NRP1, CTEMP, ONE, Q( 1,KTEMP ), LDQ, &
                VT2( KTEMP+SBM1,NLPB1 ),LDVT2,ZERO,VT( SB,NLPB1 ),LDVT )
!
      RETURN
!
!     End of MDBDSVD3
!
      END

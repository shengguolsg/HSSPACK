!  Definition:
!  ===========
!
!       SUBROUTINE MDLAED3( K,N,N1,DLAMDA,RHO,Z,D,DIFL,DIFR,ALPHA,Q,LDQ,WORK,INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            K, N, N1, LDQ, INFO
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   DLAMDA( * ), Z( * ), D( * ), DIFL(*), Q( LDQ, * ), 
!      $                   WORK( * ), ALPHA( * ), DIFR(*)
!       ..
!  
!
! \par Purpose:
!  =============
!
! MDLAED3 finds the roots of the secular equation, as defined by the
! values in DLAMDA, Z, and RHO, between 1 and K, i.e. the eigensystems of
! M = diag(DLAMDA) + RHO * Z * Z^T. 
! It makes the appropriate calls to DLAED4 and get the updated eigenvalues
! D, the distance between D and DLAMDA, DIFL, and updated Z. 
! Then, it calls DHSSEVS to construct an HSS approximation to the eigenvectors
! of M, and updates the eigenvectors via HSS matrix multiplication. 
! The differences between MDLAED3 and DLAED3 is construct the eigenvectors 
! explicitly or not. 
!
! This code makes very mild assumptions about floating point
! arithmetic. It will work on machines with a guard digit in
! add/subtract, or on those binary machines without guard digits
! which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
! It could conceivably fail on hexadecimal or decimal machines
! without guard digits, but we know of none.
!
!  Arguments:
!  ==========
!
! \param[in] K
!          K is INTEGER
!          The number of terms in the rational function to be solved by
!          DLAED4.  K >= 0.
!
! \param[in] N
!          N is INTEGER
!          The number of rows and columns in the Q matrix.
!          N >= K (deflation may result in N>K).
!
! \param[in] N1     
!          N1 is INTEGER
!          The location of the last eigenvalue in the leading submatrix.
!          min(1,N) <= N1 <= N/2. to see later
!
! \param[in,out] DLAMDA
!          DLAMDA is DOUBLE PRECISION array, dimension (K)
!          The first K elements of this array contain the old roots
!          of the deflated updating problem.  These are the poles
!          of the secular equation. May be changed on output by
!          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
!          Cray-2, or Cray C-90, as described above.
!
! \param[in] RHO
!          RHO is DOUBLE PRECISION
!          The value of the parameter in the rank one update equation.
!          RHO >= 0 required.
!
! \param[in,out] Z
!          Z is DOUBLE PRECISION array, dimension (K)
!          The first K elements of this array contain the components
!          of the deflation-adjusted updating vector. Destroyed on
!          output.
!
! \param[out] D
!          D is DOUBLE PRECISION array, dimension (N)
!          D(I) contains the updated eigenvalues for
!          1 <= I <= K.
!
! \param[in,out] DIFL
!          DIFL is DOUBLE PRECISION array, dimension (K)
!          The first K elements of this array contain the distances
!          between D(i) and DLAMDA(i), for 1 <= i <=K, which is positive
!          and DIFL(i) = D(i) - DLAMDA(i).
!
! \param[in,out] DIFR
!          DIFR is DOUBLE PRECISION array, dimension (K)
!          The first K elements of this array contain the distances
!          between D(i) and DLAMDA(i+1), for 1 <= i <=K-1, which is negative
!          and DIFR(i) = D(i) - DLAMDA(i+1).
!
! \param[in,out] ALPHA
!          ALPHA is DOUBLE PRECISION array, dimension (K)
!          The first K elements of this array contain recipral of the 2-norm
!          of each column of the eigenvector matrix of M.
!
! \param[out] Q
!          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!          Initially the first K columns contain undeflated eigenvectors.
!          On output the columns 1 to K contain the updated eigenvectors.
!
! \param[in] LDQ
!          LDQ is INTEGER
!          The leading dimension of the array Q.  LDQ >= max(1,N).
!
! \param[in, out] WORK
!          WORK is DOUBLE PRECISION array, dimension ( 2*K )
!          Used to compute DIFL, DIFR, and updated Z during solving secular 
!          equation.
!
! \param[out] INFO
!          INFO is INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an eigenvalue did not converge
!
!  Further details:
!  ===============
!
!  Modified from DLAED3 in LAPACK-3.4.2, on April 16th, 2013
!
!  =====================================================================
      SUBROUTINE MDLAED33( K,N,D,RHO,DLAMDA,Z,DIFL,DIFR,ALPHA,WORK,&
           INFO)
!
      USE CauchyHssEig
      IMPLICIT NONE
!  -- HSSPACK computational routine (version 1.0.0) --
!  -- HSSPACK is a software package provided by Nat. Univ. of Defense Tech. --
!     September 2013
!
!     .. Scalar Arguments ..
       INTEGER            K, N, N1, LDQ, INFO
       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
       DOUBLE PRECISION   DLAMDA(*), Z(*), D(*), DIFL(*), &
                          WORK( * ), ALPHA( * ), DIFR(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, II, J, N12, N2, N23, Ni
      DOUBLE PRECISION   TEMP, DIFLJ, DLAMJ, DIFRJ, DLAMJP
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLAED4, DLASET, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT, ABS, MINVAL, MAXVAL
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
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDLAED33', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.0 )   RETURN
!
!     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
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
      DO 10 I = 1, K
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
   10 CONTINUE
!
      CALL DLASET( 'A', K, 1, ONE, ONE, WORK( K+1 ), K )
!
      DO 20 J = 1, K
         CALL DLAED4( K, J, DLAMDA, Z, WORK(1), RHO, D( J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
         IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'DLAED4', -INFO )
            RETURN
         END IF
         DIFL(J) = -WORK(J)
         DIFR(J) = -WORK(J+1)
         WORK(K+J) = WORK(K+J)*WORK(J)
         DO 40 I = 1, J - 1
            WORK(K+I)= WORK(K+I)*( WORK(I) / ( DLAMDA(I)-DLAMDA(J) ) )
 40      CONTINUE
         DO 50 I = J + 1, K
            WORK(K+I)= WORK(K+I)*( WORK(I) / ( DLAMDA(I)-DLAMDA(J) ) )
 50      CONTINUE
!
 20   CONTINUE
!
!     Compute updated Z.
      DO 60 I = 1, K
         Z( I ) = SIGN( SQRT( ABS( WORK( K+I ) ) ), Z( I ) )
 60   CONTINUE
!
!     Compute scalars for eigenvector matrix
      DO J = 1, K
         DIFLJ = DIFL( J )
         DLAMJ = -DLAMDA( J )
         IF( J.LT.K ) THEN
            DIFRJ  = -DIFR( J )
            DLAMJP = -DLAMDA( J+1 )
         END IF
         WORK( J ) = Z( J ) / DIFLJ 
         DO I = 1, J - 1
            WORK( I ) = Z( I ) / ( ( DLAMDA( I )+ DLAMJ )-DIFLJ )
         END DO
         DO I = J + 1, K
            WORK( I ) = Z( I ) / ( ( DLAMDA( I )+DLAMJP )+DIFRJ ) 
         END DO
         TEMP = DNRM2( K, WORK, 1 )
         ALPHA( J ) = ONE / TEMP
      END DO
!
!     Calling DHSSEVC to update eigenvectors
!
!      
      RETURN
!
!     End of MDLAED33
!
      END

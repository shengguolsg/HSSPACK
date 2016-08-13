      SUBROUTINE MDLASD31( K,D,Z,DSIGMA,DIFL,DIFR,ALPHA,INFO )
!
        implicit none

!     .. Scalar Arguments ..
      INTEGER            INFO, K
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D(*),Z(*),DSIGMA(*),DIFL(*),DIFR(*),ALPHA(*)
!     ..
!
!  Purpose
!  =======
!
!  MDLASD31 finds the square roots of the roots of the secular equation,
!  as defined by the values in DSIGMA and Z. It makes the appropriate
!  calls to DLASD4.
!
!  MDLASD31 is called from computing the right singular vector matrix of
!  an broken arrow matrix.
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
!  Provided by Shengguo Li, based on dlasd8 written Ming Gu and Huan Ren, 
!  Computer Science Division, University of California at Berkeley, USA
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
!     Update VF and VL.
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
         ALPHA( J ) = one / TEMP
      END DO

      deallocate( WORK )

!     End of MDLASD31
!
    END SUBROUTINE MDLASD31

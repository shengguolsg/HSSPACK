    SUBROUTINE MDBDSVD8( K,SB,D,ZZ,DSIGMA,DIFL,DIFR,APHAU,APHAV,Z, &
                        UT,LDUT,VT,LDVT,WORK,INFO,nod )
      USE Cauchyhssvd_VPS
      IMPLICIT NONE
!
!  -- HSSPACK routine (version 0.0.1) --
!  -- HSSPACK is a software package provided by National Univ. of Defense Tech.    --
!  -- Novermber 2012                                                  --
!
!     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: K,LDUT,LDVT,SB,nod
      INTEGER, INTENT(OUT) :: INFO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION, INTENT(OUT) ::  D(*), DIFL(*),DIFR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: APHAV(*), Z(*),ZZ(LDVT,*),APHAU(*),&
                        DSIGMA(*), UT(LDUT,K), VT(K,LDVT), WORK(*)
!  ..
!  Purpose
!  =======
!
!  MDBDSVD8 finds the square roots of the roots of the secular equation,
!  as defined by the values in DSIGMA and Z. It makes the appropriate
!  calls to MDLASD4, and stores, for each  element in D, the distance
!  to its two nearest poles (elements in DSIGMA). It also updates
!  the left and right singular vectors explicitly. 
!
!  It computes HSS approximations to both left and right singular vector
!  matrices of M = UH * W* VHT, where these matrices are all K-by-K. 
!  We form the HSS approximations to UH and VHT without forming them
!  explicitly. 
!
!  MDLASD8 is called from MDLASD1.
!
!  Arguments
!  =========
!
!  K       (input) INTEGER
!          The number of terms in the rational function to be solved
!          by DLASD4.  1 <= K <= LDUT.
!
!  D       (output) DOUBLE PRECISION array, dimension ( K )
!          On output, D contains the updated singular values.
!
!  ZZ      (input/output) DOUBLE PRECISION array, dimension ( LDVT,SB )
!          On entry, the elements of the last column contain the
!          components of the deflation-adjusted updating row vector. The
!          elements from SB+1 to N of the first SB-1 columns of Z will be 
!          multiplied by the right singular vector matrix of that broken arrow
!          matrix. On exit, Z is updated. 
!
!  DSIGMA  (input/output) DOUBLE PRECISION array, dimension ( K )
!          On entry, the first K elements of this array contain the old
!          roots of the deflated updating problem.  These are the poles
!          of the secular equation.
!          On exit, the elements of DSIGMA may be very slightly altered
!          in value.
!
!  DIFL    (output) DOUBLE PRECISION array, dimension ( K )
!          On exit, DIFL(I) = D(I) - DSIGMA(I).
!
!  DIFR    (output) DOUBLE PRECISION array, dimension ( K )
!          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K) is not
!          defined and will not be referenced.
!
!  APHAU  (output) DOUBLE PRECISION array, dimension ( K )
!          The column normalization factors of U
!
!  APHAV  (output) DOUBLE PRECISION array, dimension ( K )
!          The row normalization factors of VT
!
!  Z      (input) DOUBLE PRECISION array, dimension ( K ) 
!         A temp copy of Z.
!
!  UT     (input/output) DOUBLE PRECISION array, dimension ( K, LDUT).
!          Here UT is the transpose of first K columns of U. 
!
! LDUT     (input) INTEGER, the number of columns of UT
!
!  VT      (input/output) DOUBLE PRECISION array, dimension (K, LDVT)
!          The first K rows of right singular vector matrix.
!
! LDVT     (input) INTEGER, the number of columns of VT
! 
!  WORK    (workspace) DOUBLE PRECISION array, dimension at least 3 * K +(N^2)
!          It is used to construct HSS approximations and also used during HSS
!          matrix multiplication. To store one HSS approximation it needs at least
!          MAX(LDUT, LDVT)*Mi*8, Mi is the row dimension of the matrix at bottom of HSS tree. 
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!  Based on LAPACK routine, DLASD8
!  Written by S.-G. Li, in Changsha, Dec. 1st, 2012
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J, Ni, NDD, &
                         IDD, IH, ILEFT
      DOUBLE PRECISION   DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, RHO, TEMP
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLASCL, DLASD4, DLASET, XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DDOT, DLAMC3, DNRM2
      EXTERNAL           DDOT, DLAMC3, DNRM2
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
      IF( ( K.LT.1 ) ) THEN
         INFO = -1
      ELSE IF( LDUT.LT.K ) THEN
         INFO = -11
      ELSE IF( LDVT.LT.K ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSVD8', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
!!$      IF( K.EQ.1 ) THEN
!!$         D( 1 ) = ABS( Z( 1 ) )
!!$         DIFL( 1 ) = D( 1 )
!!$!         
!!$         DIFL( 2 ) = ONE
!!$         DIFR( 1 ) = ONE
!!$!
!!$         RETURN
!!$      END IF
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
      DO 10 I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   10 CONTINUE
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
      Z(1:K) = ZZ(SB:SB+K-1,SB)
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
      DO 40 J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ), &
                     WORK( IWK2 ), INFO )
!
!        If the root finder fails, the computation is terminated.
         IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'DLASD4', -INFO )
            write(*,*) 'DLASQ4 brokes at J=', J
            RETURN
         END IF
         WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
         DIFL( J ) = -WORK( J )
         DIFR( J ) = -WORK( J+1 )
         DO 20 I = 1, J - 1
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* &
                            WORK( IWK2I+I ) / ( DSIGMA( I )- &
                            DSIGMA( J ) ) / ( DSIGMA( I )+ &
                            DSIGMA( J ) )
   20    CONTINUE
         DO 30 I = J + 1, K
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* &
                            WORK( IWK2I+I ) / ( DSIGMA( I )- &
                            DSIGMA( J ) ) / ( DSIGMA( I )+ &
                            DSIGMA( J ) )
   30    CONTINUE
   40 CONTINUE
!
!     Compute updated Z.
!
      DO 50 I = 1, K
         Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
   50 CONTINUE
      ZZ( SB:SB+K-1,SB ) = Z( 1:K )
!
!     Compute scalars for right and left singular vectors
!
      DO 80 J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         IF( J.LT.K ) THEN
            DIFRJ = -DIFR( J )
            DSIGJP = -DSIGMA( J+1 )
         END IF
         WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
         DO 60 I = 1, J - 1
            WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJ )-DIFLJ ) &
                       / ( DSIGMA( I )+DJ )
   60    CONTINUE
         DO 70 I = J + 1, K
            WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJP )+DIFRJ ) &
                       / ( DSIGMA( I )+DJ )
   70    CONTINUE
         TEMP = DNRM2( K, WORK, 1 )
         APHAV( J ) = ONE / TEMP
   80 CONTINUE
!     
      DO J = 1, K
         DIFLJ = DIFL( J )
         DJ    = D( J )
         DSIGJ = -DSIGMA( J )
         IF( J.LT.K ) THEN
            DIFRJ  = -DIFR( J )
            DSIGJP = -DSIGMA( J+1 )
         END IF
         WORK( J ) = -DSIGMA(J)*Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
         DO I = 1, J - 1
            WORK( I ) = DSIGMA(I)*Z( I ) / ( ( DSIGMA( I )+ DSIGJ )-DIFLJ ) / ( DSIGMA( I )+DJ )
         END DO
         DO I = J + 1, K
            WORK( I ) = DSIGMA(I)*Z( I ) / ( ( DSIGMA( I )+DSIGJP )+DIFRJ ) / ( DSIGMA( I )+DJ )
         END DO
         WORK(1)= -1.0D+0
         TEMP = DNRM2( K, WORK, 1 )
         APHAU( J ) = ONE / TEMP
      END DO
!
! ********* Update the singular vectors ******* 
!
      Ni  = 100    ! block size for nodes at bottom level
      CALL DHSSVCS_BD( K, SB, Ni, D, ZZ, DSIGMA, DIFL, DIFR, APHAU, APHAV, Z, &
                   UT,LDUT, VT, LDVT, WORK, INFO )
!
      RETURN
!
!     End of MDBDSVD8
!
      END


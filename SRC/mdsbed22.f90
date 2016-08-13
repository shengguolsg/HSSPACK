! 
     SUBROUTINE MDSBED22( K, N, N1, SB, D, Q, LDQ, RHO, Z, Z2, DLAMDA, W, & 
                          Q2, INDXQ, INDXP, INFO ) 
!  
        IMPLICIT NONE
!       ..
!       .. Scalar Arguments .. 
        INTEGER            INFO, K, LDQ, N, N1, SB
        DOUBLE PRECISION   RHO 
!       .. 
!       .. Array Arguments .. 
        INTEGER            INDXQ( * ), INDXP( * )
        DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ,* ), Q2( * ), & 
                           W( * ), Z( N,* ), Z2( N,* )
!       .. 
!  Purpose 
!  ========
! 
! MDSBED22 merges the two sets of eigenvalues together into a single
! sorted set.  Then it tries to deflate the size of the problem.
! 
! Different than MDSBED2, this routine does not explore any sparse structure,
! and the parameter INDXQ already contains the permutation index such that
! D( INDXQ(I) ) is the I-th smallest entry. It saves many list arrarys.
! Z and Z2 are using in a ping-pong form. The vectors are firstly stored in
! Z and on exit they are moved to Z2.
! 
! 
!  Arguments: 
!  ========== 
!
! \param[out] K 
!         K is INTEGER 
!         The number of non-deflated eigenvalues, and the order of the 
!         related secular equation. 0 <= K <=N. 
! 
! \param[in] N 
!         N is INTEGER 
!         The dimension of current problem which is merged from two 
!         submatrices.  N >= 0. 
! 
! \param[in] N1 
!         N1 is INTEGER 
!         The location of the last eigenvalue in the leading sub-matrix. 
!         min(1,N) <= N1 <= N/2. 
! 
! \param[in] SB
!         SB is INTEGER 
!         The rank of current low-rank permutations, and it equals to the column
!         dimension of Z. 
! 
! \param[in,out] D 
!          D is DOUBLE PRECISION array, dimension (N) 
!         On entry, D contains the eigenvalues of the two submatrices to 
!         be combined. 
!         On exit, D contains the trailing (N-K) updated eigenvalues 
!         (those which were deflated) sorted into increasing order. 
! 
! \param[in,out] Q 
!         Q is DOUBLE PRECISION array, dimension (LDQ, N) 
!         On entry, Q contains the eigenvectors of two submatrices in 
!         the two square blocks with corners at (1,1), (N1,N1) 
!         and (N1+1, N1+1), (N,N). 
!         On exit, Q contains the trailing (N-K) updated eigenvectors 
!         (those which were deflated) in its last N-K columns. 
! 
! \param[in] LDQ 
!          LDQ is INTEGER 
!         The leading dimension of the array Q.  LDQ >= max(1,N). 
! 
! \param[in,out] RHO 
!         RHO is DOUBLE PRECISION 
!         On entry, the off-diagonal element associated with the rank-1 
!         cut which originally split the two submatrices which are now 
!         being recombined. 
!         On exit, RHO has been modified to the value required by 
!         DLAED3. 
! 
! \param[inout] Z 
!         Z is DOUBLE PRECISION array, dimension (N,SB) 
!         On entry, Z contains the updating vectors (combination of the last KD
!         rows of the first sub-eigenvector matrix and combination of the first 
!         KD rows of the second sub-eigenvector matrix). The column dimension is
!         reduced by one after each call of this routine.  
!         On exit, the contents of Z have been destroyed by the updating 
!         process. 
! 
! \param[inout] Z2
!         Z2 is DOUBLE PRECISION array, dimension (N,SB) 
!         On entry, Z2 is used to contain the permuted vectors, which is used
!         as a copy of Z. On exit, the contents of Z2 are the updated and permutated rows
!         of matrix Z. 
! 
! \param[out] DLAMDA 
!         DLAMDA is DOUBLE PRECISION array, dimension (N) 
!         A copy of the first K eigenvalues which will be used by 
!         DLAED3 to form the secular equation. 
! 
! \param[out] W 
!         W is DOUBLE PRECISION array, dimension (N) 
!         The first K values of the final deflation-altered z-vector 
!         which will be passed to DLAED3. 
! 
! \param[out] Q2 
!         Q2 is DOUBLE PRECISION array, dimension (N1**2+(N-N1)**2) 
!         A copy of the first K eigenvectors which will be used by 
!         DLAED3 in a matrix multiply (DGEMM) to solve for the new 
!         eigenvectors. 
! 
! \param[in,out] INDXQ 
!         INDXQ is INTEGER array, dimension (N) 
!         The permutation which separately sorts the two sub-problems 
!         in D into ascending order.  INDXQ(I) is the index of column that
!         corresponds to the I-th smallest eigenvallue.!         
! 
! \param[out] INDXP 
!         INDXP is INTEGER array, dimension (N) 
!         The permutation used to place deflated values of D at the end 
!         of the array.  INDXP(1:K) points to the nondeflated D-values 
!         and INDXP(K+1:N) points to the deflated eigenvalues. The non-deflated
!         are in ascending order and the deflated ones are in descending order.
! 
! \param[out] INFO 
!          INFO is INTEGER 
!          = 0:  successful exit. 
!          < 0:  if INFO = -i, the i-th argument had an illegal value. 
! 
!  Authors: 
!  ======== 
!
! \author Nat. Univ. of Def. Tech. 
! \date October 2014 
! 
! \par Contributors: 
!  ================== 
! 
! Shengguo Li, NUDT, China, 2014
! 
! ===================================================================== 
!     ..
!     .. Parameters .. 
      DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT 
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0,  &
                         TWO = 2.0D0, EIGHT = 8.0D0 ) 
!     .. 
!     .. Local Scalars .. 
      INTEGER            I, IMAX, IQ1, IQ2, J, JMAX, JS, K2, N1P1, &
                         N2, NJ, PJ, KI, INDXI, JP
      DOUBLE PRECISION   C, EPS, S, T, TAU, TOL 
!     .. 
!     .. External Functions .. 
      INTEGER            IDAMAX 
      DOUBLE PRECISION   DLAMCH, DLAPY2 
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2 
!     .. 
!     .. External Subroutines .. 
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA 
!     .. 
!     .. Intrinsic Functions .. 
      INTRINSIC          ABS, MAX, MIN, SQRT 
!     .. 
!     .. Executable Statements .. 
! 
!     Test the input parameters. 
! 
      INFO = 0 
! 
      IF( N.LT.0 ) THEN 
         INFO = -2 
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN 
         INFO = -7
!      ELSE IF( MIN( 1, ( N/2 ) ).GT.N1 .OR. ( N/2 ).LT.N1 ) THEN 
      ELSE IF( MIN( 1, ( N/2 ) ).GT.N1 ) THEN 
         INFO = -3 
      END IF 
      IF( INFO.NE.0 ) THEN 
         CALL XERBLA( 'MDSBED22', -INFO ) 
         RETURN 
      END IF 
! 
!     Quick return if possible 
! 
      IF( N.EQ.0 ) &
         RETURN 
! 
      KI = 1       ! Column index for Z-matrix
      N2 = N - N1 
      N1P1 = N1 + 1 
! 
      IF( RHO.LT.ZERO ) THEN 
         CALL DSCAL( N2, MONE, Z( N1P1,1 ), 1 )   ! Scale a vector by a constant
      END IF 
! 
!     Normalize z so that norm(z) = 1.  Since z is the concatenation of 
!     two normalized vectors, norm2(z) = sqrt(2). 
! 
      T = ONE / SQRT( TWO ) 
      CALL DSCAL( N, T, Z(1,1), 1 )   ! only deflate the first column
! 
!     RHO = ABS( norm(z)**2 * RHO ) 
! 
      RHO = ABS( TWO*RHO ) 
! 
!     Sort the eigenvalues into increasing order
! 
      KI = 1
      DO I = 1, N 
         INDXI = INDXQ( I )   ! the original index for the I-the smallest eigenvalue
         DLAMDA( I ) = D( INDXI )
      END DO
!
!     Calculate the allowable deflation tolerance 
! 
      IMAX = IDAMAX( N, Z( 1,KI ), 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS  = DLAMCH( 'Epsilon' ) 
      TOL  = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX,KI ) ) ) 
! 
!     If the rank-1 modifier is small enough, no more needs to be done 
!     except to reorganize Q so that its columns correspond with the 
!     elements in D. 
! 
      IF( RHO*ABS( Z( IMAX,KI ) ).LE.TOL ) THEN 
         K = 0 
         IQ2 = 1 
         DO J = 1, N 
            I = INDXQ( J ) 
            CALL DCOPY( N, Q( 1, I ), 1, Q2( IQ2 ), 1 ) 
            DLAMDA( J ) = D( I ) 
            IQ2 = IQ2 + N 
         END DO
         CALL DLACPY( 'A', N, N, Q2, N, Q, LDQ ) 
         CALL DCOPY( N, DLAMDA, 1, D, 1 ) 
         GO TO 190       ! the size of secular equation is zero
      END IF 
! 
!     If there are multiple eigenvalues then the problem deflates.  Here 
!     the number of equal eigenvalues are found.  As each equal 
!     eigenvalue is found, an elementary reflector is computed to rotate 
!     the corresponding eigensubspace so that the corresponding 
!     components of Z are zero in this new basis. 
! 
! 
      K = 0               ! increase toward N
      K2 = N + 1          ! decrease toward one, such that K+K2=N+1
      DO J = 1, N 
         NJ = INDXQ( J )   ! the Jth smallest
         IF( RHO*ABS( Z( NJ,KI ) ).LE.TOL ) THEN 
! 
!           Deflate due to small z component. 
! 
            K2 = K2 - 1 
            INDXP( K2 ) = NJ    ! records the deflated ones, which moved to the back part
            IF( J.EQ.N ) &
               GO TO 100 
         ELSE 
            PJ = NJ    ! PJ means Previous J
            GO TO 80 
         END IF 
      END DO
   80 CONTINUE 
      J = J + 1 
      NJ = INDXQ( J )           ! next larger one
      IF( J.GT.N ) GO TO 100 
      IF( RHO*ABS( Z( NJ, KI) ).LE.TOL ) THEN  ! next z is small enough
! 
!        Deflate due to small z component. 
! 
         K2 = K2 - 1 
         INDXP( K2 ) = NJ   ! deflated eigenvalues
      ELSE 
! 
!        Check if eigenvalues are close enough to allow deflation. 
! 
         S = Z( PJ,KI ) 
         C = Z( NJ,KI ) 
! 
!        Find sqrt(a**2+b**2) without overflow or 
!        destructive underflow. 
! 
         TAU = DLAPY2( C, S ) 
         T = D( NJ ) - D( PJ ) 
         C = C / TAU 
         S = -S / TAU 
         IF( ABS( T*C*S ).LE.TOL ) THEN  
! 
!           Deflation is possible. 
! 
            Z( NJ,KI ) = TAU 
            Z( PJ,KI ) = ZERO 
!
!           Apply the Givens rotation to the Q matrix and the other rows in Z-matrix.
            CALL DROT( N, Q( 1, PJ ), 1, Q( 1, NJ ), 1, C, S )   ! Is it faster to use AQ ?
            IF( SB .GT. 1 ) THEN
               CALL DROT( SB-1, Z(PJ,2),N, Z(NJ,2), N, C, S )    ! Here should be wrong ? G or G^T
            END IF

            T = D( PJ )*C**2 + D( NJ )*S**2 
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2 
            D( PJ ) = T 
            K2 = K2 - 1 
            I = 1 
   90       CONTINUE 
            IF( K2+I.LE.N ) THEN   ! order the deflated eigenvalues to descending order
               IF( D( PJ ).LT.D( INDXP( K2+I ) ) ) THEN 
                  INDXP( K2+I-1 ) = INDXP( K2+I ) 
                  INDXP( K2+I ) = PJ 
                  I = I + 1 
                  GO TO 90 
               ELSE 
                  INDXP( K2+I-1 ) = PJ 
               END IF 
            ELSE 
               INDXP( K2+I-1 ) = PJ 
            END IF 
            PJ = NJ 
         ELSE 
            K = K + 1 
            DLAMDA( K ) = D( PJ )  ! the non-deflated eigenvalues
            W( K ) = Z( PJ,KI )    ! W stores the non-deflated Z in the right order.
            INDXP( K ) = PJ 
            PJ = NJ 
         END IF 
      END IF 
      GO TO 80 
  100 CONTINUE 
! 
!     Record the last eigenvalue. 
! 
      K = K + 1 
      DLAMDA( K ) = D( PJ )   ! DLAMDA is in ascending order now
      W( K ) = Z( PJ,KI ) 
      INDXP( K ) = PJ 
! 
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA 
!     and Q2 respectively.  The eigenvalues/vectors which were not 
!     deflated go into the first K slots of DLAMDA and Q2 respectively, 
!     while those which were deflated go into the last N - K slots. 
! 
      I = 1 
      IQ1 = 1 
      DO J = 1, N
         JP = INDXP( I )    ! The permuted columns, first K nondeflated, last N-K deflated.
         JS = INDXQ( JP )   ! The original column in Q
         CALL DCOPY( N, Q( 1,JP ), 1, Q2( IQ1 ), 1 )    ! Is it JS or JP ??
         IF( SB .GT. 1 ) THEN
            CALL DCOPY( SB-1, Z(JP,2), N, Z2(I,1), N )  ! start from the first column
         END IF
         Z( I,KI ) = D( JP ) 
         I = I + 1 
         IQ1 = IQ1 + N
      END DO
!
!     The deflated eigenvalues and their corresponding vectors go back 
!     into the last N - K slots of D and Q respectively. 
! 
      IQ2 = 1 + N*K  ! The starting position for the deflated eigenvectors
      IF( K.LT.N ) THEN 
         CALL DLACPY( 'A',N,N-K,Q2(IQ2),N, Q(1,K+1), LDQ )   ! copy back to Q from Q2
         CALL DCOPY( N-K, Z( K+1,KI ), 1, D( K+1 ), 1 )      ! copy back to D
      END IF
! 
! 
  190 CONTINUE 
      RETURN 
! 
!     End of MDSBED22
! 
    END SUBROUTINE MDSBED22

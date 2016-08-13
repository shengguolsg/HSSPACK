! 
       SUBROUTINE MDSBED2( K, N, N1, KD, D, Q, LDQ, INDXQ, RHO, Z, Z2, DLAMDA, W, & 
                          Q2, INDX, INDXC, INDXP, COLTYP, HFLAG, BSMLZ, INFO ) 
!  
        IMPLICIT NONE
!       ..
!       .. Scalar Arguments .. 
        LOGICAL            HFLAG
        INTEGER            INFO, K, LDQ, N, N1, KD, BSMLZ
        DOUBLE PRECISION   RHO 
!       .. 
!       .. Array Arguments .. 
        INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ), & 
                           INDXQ( * ) 
        DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ), & 
                           W( * ), Z( N,* ), Z2( N,* )
!       .. 
!  Purpose 
!  ========
! 
! MDSBED2 merges the two sets of eigenvalues together into a single 
! sorted set.  Then it tries to deflate the size of the problem. 
! There are two ways in which deflation can occur:  when two or more 
! eigenvalues are close together or if there is a tiny entry in the 
! Z vector.  For each such occurrence the order of the related secular 
! equation problem is reduced by one. 
! Different than DLAED2, this routine updates the other vectors in Z. 
! This routine merges these two submatrices and deflate the first column
! of Z, and the other columns of Z are deflated in a similar routine 
! named MDSBED22 which is simpler than this one since it does not need to
! explore any sparse structrue. Z and Z2 are using in a ping-pong form. 
! 
!  Arguments: 
!  ========== 
!
! \param[out] K 
!          K is INTEGER 
!         The number of non-deflated eigenvalues, and the order of the 
!         related secular equation. 0 <= K <=N. 
! 
! \param[in] N 
!          N is INTEGER 
!         The dimension of current problem which is merged from two 
!         submatrices.  N >= 0. 
! 
! \param[in] N1 
!          N1 is INTEGER 
!         The location of the last eigenvalue in the leading sub-matrix. 
!         min(1,N) <= N1 <= N/2. 
! 
! \param[in] KD
!          KD is INTEGER 
!         The semi-bandwith of the symmetric band matrix
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
! \param[in,out] INDXQ 
!         INDXQ is INTEGER array, dimension (N) 
!         The permutation which separately sorts the two sub-problems 
!         in D into ascending order.  Note that elements in the second 
!         half of this permutation must first have N1 added to their 
!         values. Destroyed on exit. 
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
!          Z is DOUBLE PRECISION array, dimension (N,KD) 
!         On entry, Z contains the updating vectors (combination of the last KD
!         row of the first sub-eigenvector matrix and combination of the first 
!         KD row of the second sub-eigenvector matrix). 
!         On exit, the contents of Z have been destroyed by the updating 
!         process. 
! 
! \param[inout] Z2
!          Z2 is DOUBLE PRECISION array, dimension (N,KD) 
!         On entry, Z2 is used to contain the permuted vectors as a copy of Z. 
!         On exit, the contents of Z2 are the updated and permutated of matrix 
!         Z. 
! 
! \param[out] DLAMDA 
!          DLAMDA is DOUBLE PRECISION array, dimension (N) 
!         A copy of the first K eigenvalues which will be used by 
!         DLAED3 to form the secular equation. 
! 
! \param[out] W 
!          W is DOUBLE PRECISION array, dimension (N) 
!         The first K values of the final deflation-altered z-vector 
!         which will be passed to DLAED3. 
! 
! \param[out] Q2 
!          Q2 is DOUBLE PRECISION array, dimension (N1**2+(N-N1)**2) 
!         A copy of the first K eigenvectors which will be used by 
!         DLAED3 in a matrix multiply (DGEMM) to solve for the new 
!         eigenvectors. 
! 
! \param[out] INDX 
!          INDX is INTEGER array, dimension (N) 
!         The permutation used to sort the contents of DLAMDA into 
!         ascending order. 
! 
! \param[out] INDXC 
!          INDXC is INTEGER array, dimension (N) 
!         The permutation used to arrange the columns of the deflated 
!         Q matrix into three groups:  the first group contains non-zero 
!         elements only at and above N1, the second contains 
!         non-zero elements only below N1, and the third is dense. 
! 
! \param[out] INDXP 
!          INDXP is INTEGER array, dimension (N) 
!         The permutation used to place deflated values of D at the end 
!         of the array.  INDXP(1:K) points to the nondeflated D-values 
!         and INDXP(K+1:N) points to the deflated eigenvalues. 
! 
! \param[out] COLTYP 
!          COLTYP is INTEGER array, dimension (N) 
!         During execution, a label which will indicate which of the 
!         following types a column in the Q2 matrix is: 
!         1 : non-zero in the upper half only; 
!         2 : dense; 
!         3 : non-zero in the lower half only; 
!         4 : deflated. 
!         On exit, COLTYP(i) is the number of columns of type i, 
!         for i=1 to 4 only. 
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
!     .. Local Arrays .. 
      INTEGER            CTOT( 4 ), PSM( 4 ) 
!     .. 
!     .. Local Scalars .. 
      INTEGER            CT, I, IMAX, IQ1, IQ2, J, JMAX, JS, K2, N1P1, &
                         N2, NJ, PJ, KI
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
      ELSE IF( MIN( 1, ( N/2 ) ).GT.N1 .OR. ( N/2 ).LT.N1 ) THEN 
         INFO = -3 
      END IF 
      IF( INFO.NE.0 ) THEN 
         CALL XERBLA( 'MDSBED2', -INFO ) 
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
      DO I = N1P1, N 
         INDXQ( I ) = INDXQ( I ) + N1 
      END DO
! 
!     re-integrate the deflated parts from the last pass 
! 
      DO I = 1, N 
         DLAMDA( I ) = D( INDXQ( I ) ) 
      END DO
!
!     DLAMDA is broken into two pieces, the first N1 entries are in 
!     ascending order, and the last N2 are in ascending order, and they are
!     merged together by calling DLAMRG. 
!
      CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDXC ) 
      DO I = 1, N 
         INDX( I ) = INDXQ( INDXC( I ) )     ! DLAMDA( INDX(I) ) is the I-th smallest one
      END DO
!
!     INDXQ stores the original permutations of two subproblems separately.
!     INDX  stores the permutations after merging them together. INDX(I)
!           means the index of I-th smallest entry in the original merged
!           array, D=[D1  D2];
! 
!     Calculate the allowable deflation tolerance 
! 
      IMAX = IDAMAX( N, Z(1,KI), 1 ) 
      JMAX = IDAMAX( N, D, 1 ) 
      EPS = DLAMCH( 'Epsilon' ) 
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX,KI ) ) ) 
! 
!     If the rank-1 modifier is small enough, no more needs to be done 
!     except to reorganize Q so that its columns correspond with the 
!     elements in D. 
! 
      IF( RHO*ABS( Z( IMAX,KI ) ).LE.TOL ) THEN 
         K = 0 
         IQ2 = 1 
         DO J = 1, N 
            I = INDX( J ) 
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
      DO I = 1, N1 
         COLTYP( I ) = 1  ! the first N1 are nonzero
      END DO
      DO I = N1P1, N 
         COLTYP( I ) = 3  ! the last N2 are nonzero
      END DO
! 
      K = 0               ! increase toward N
      K2 = N + 1          ! decrease toward one, such that K+K2=N+1
      DO J = 1, N 
         NJ = INDX( J )   ! the Jth smallest
         IF( RHO*ABS( Z( NJ,KI ) ).LE.TOL ) THEN 
! 
!           Deflate due to small z component. 
! 
            K2 = K2 - 1 
            COLTYP( NJ ) = 4    ! deflated 
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
      NJ = INDX( J )           ! next larger one
      IF( J.GT.N ) GO TO 100 
      IF( RHO*ABS( Z( NJ, KI) ).LE.TOL ) THEN  ! next z is small enough
! 
!        Deflate due to small z component. 
! 
         K2 = K2 - 1 
         COLTYP( NJ ) = 4 
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
            IF( COLTYP( NJ ).NE.COLTYP( PJ ) ) &  ! change its type to 2
               COLTYP( NJ ) = 2 
            COLTYP( PJ ) = 4             ! change its type to deflated

!           Apply the Givens rotation to the Q matrix and the other rows in Z-matrix.
            CALL DROT( N, Q( 1, PJ ), 1, Q( 1, NJ ), 1, C, S )    ! Is it faster to use AQ ?
            CALL DROT( KD-1, Z(PJ,2),N, Z(NJ,2), N, C, S )        ! Here should be wrong ? G or G^T ??

            T = D( PJ )*C**2 + D( NJ )*S**2 
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2 
            D( PJ ) = T 
            K2 = K2 - 1 
            I = 1 
   90       CONTINUE 
            IF( K2+I.LE.N ) THEN   ! order the deflated eigenvalues to ascending order
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
      IF( K .LT. BSMLZ ) THEN    ! Sparse structure
!     Count up the total number of the various types of columns, then 
!     form a permutation which positions the four column types into 
!     four uniform groups (although one or more of these groups may be 
!     empty). 
! 
         DO J = 1, 4 
            CTOT( J ) = 0 
         END DO
         DO J = 1, N 
            CT = COLTYP( J ) 
            CTOT( CT ) = CTOT( CT ) + 1 
         END DO
! 
!     PSM(*) = Position in SubMatrix (of types 1 through 4) 
! 
         PSM( 1 ) = 1 
         PSM( 2 ) = 1 + CTOT( 1 ) 
         PSM( 3 ) = PSM( 2 ) + CTOT( 2 ) 
         PSM( 4 ) = PSM( 3 ) + CTOT( 3 ) 
         K = N - CTOT( 4 ) 
! 
!     Fill out the INDXC array so that the permutation which it induces 
!     will place all type-1 columns first, all type-2 columns next, 
!     then all type-3's, and finally all type-4's. 
! 
         DO J = 1, N 
            JS = INDXP( J )           ! JS is the original column
            CT = COLTYP( JS ) 
            INDX( PSM( CT ) ) = JS 
            INDXC( PSM( CT ) ) = J    ! Using INDXP INDXC can point to the original col
            PSM( CT ) = PSM( CT ) + 1 ! each type increases by one
         END DO
! 
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA 
!     and Q2 respectively.  The eigenvalues/vectors which were not 
!     deflated go into the first K slots of DLAMDA and Q2 respectively, 
!     while those which were deflated go into the last N - K slots. 
! 
         I = 1 
         IQ1 = 1 
         IQ2 = 1 + ( CTOT( 1 )+CTOT( 2 ) )*N1 
         DO J = 1, CTOT( 1 ) 
            JS = INDX( I ) 
            CALL DCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 ) 
            CALL DCOPY( KD-1, Z(JS,2), N, Z2(I,1), N )  ! start from the first column
            Z( I,KI ) = D( JS ) 
            I = I + 1 
            IQ1 = IQ1 + N1 
         END DO
! 
         DO J = 1, CTOT( 2 ) 
            JS = INDX( I ) 
            CALL DCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 ) 
            CALL DCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 ) 
            CALL DCOPY( KD-1, Z(JS,2), N, Z2(I,1), N )  ! start from the first column
            Z( I,KI ) = D( JS ) 
            I = I + 1 
            IQ1 = IQ1 + N1 
            IQ2 = IQ2 + N2 
         END DO
!
         DO J = 1, CTOT( 3 ) 
            JS = INDX( I ) 
            CALL DCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 ) 
            CALL DCOPY( KD-1, Z(JS,2), N, Z2(I,1), N )  ! start from the first column
            Z( I,KI ) = D( JS ) 
            I = I + 1 
            IQ2 = IQ2 + N2 
         END DO
!
         IQ1 = IQ2 
         DO J = 1, CTOT( 4 ) 
            JS = INDX( I ) 
            CALL DCOPY( N, Q( 1, JS ), 1, Q2( IQ2 ), 1 ) 
            CALL DCOPY( KD-1, Z(JS,2), N, Z2(I,1), N )  ! start from the first column
            IQ2 = IQ2 + N 
            Z( I,KI ) = D( JS ) 
            I = I + 1 
         END DO
! 
      ELSE               ! not use LAPACK technique and use HSS instead
!
         HFLAG = .FALSE. 
         DO J = 1, N
            JS = INDXP( J )
            INDX( J ) = JS
         END DO
!
!     The first K entries of INDX store the original index of non-deflated 
!     eigenvectors, and the last CTOT(4) entries store the original index 
!     of deflated eigenvectors.
!
!     Copy the first K non-deflated eigenvectors to the first K columns of 
!     Q2, which will be treated as an N-by-K matrix. 
         I = 1
         IQ1 = 1
         IQ2 = 1 + K * N
         DO 210 J = 1, K     ! the K nondeflated eigenvectors
            JS = INDX( I )
            CALL DCOPY( N, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
            CALL DCOPY( KD-1, Z(JS,2), N, Z2(I,1), N )  ! start from the first column
            Z( I,KI ) = D( JS )
            I = I + 1
            IQ1 = IQ1 + N
210      CONTINUE
!
         IQ2 = IQ1
         DO 220 J = 1, N-K  ! the CTOT(4) deflated eigenvectors
            JS = INDX( I )
            CALL DCOPY( N, Q( 1, JS ), 1, Q2( IQ2 ), 1 )
            CALL DCOPY( KD-1, Z(JS,2), N, Z2(I,1), N )  ! start from the first column
            IQ2 = IQ2 + N
            Z( I,KI ) = D( JS )
            I = I + 1
220      CONTINUE

      END IF  ! ( HFLAG )
         
!     The deflated eigenvalues and their corresponding vectors go back 
!     into the last N - K slots of D and Q respectively. 
! 
      IF( .not. HFLAG ) THEN
!
         CALL DLACPY( 'A', N, N-K, Q2(IQ1), N, Q(1,K+1), LDQ )  ! copy the whole matrix Q2 ??
         CALL DCOPY( N-K, Z( K+1,KI ), 1, D( K+1 ), 1 )
!
      ELSE  IF( K.LT.N ) THEN
         CALL DLACPY( 'A',N,CTOT( 4 ),Q2( IQ1 ),N,  &
                      Q( 1,K+1 ), LDQ )     ! copy back to Q from Q2
         CALL DCOPY( N-K, Z( K+1,KI ), 1, D( K+1 ), 1 )  ! copy back to D
      END IF          
! 
!     Copy CTOT into COLTYP for referencing in DLAED3. 
! 
      DO J = 1, 4 
         COLTYP( J ) = CTOT( J ) 
      END DO
! 
  190 CONTINUE 
      RETURN 
! 
!     End of MDSBED2 
! 
    END SUBROUTINE MDSBED2

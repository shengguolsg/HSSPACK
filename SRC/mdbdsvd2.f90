    SUBROUTINE MDBDSVD2( NL,NR,SB,SQRE,K,D,Z,ALPHA,BETA,U,LDU,VT,&
             LDVT,DSIGMA,U2,LDU2,VT2,LDVT2,IDXP,IDX,IDXC,IDXQ, &
             COLTYP, WORK, HFLAG, INFO )
!
      use basicmm
      IMPLICIT NONE
!
!  -- HSSPACK computational routine ( version 1.0.0 ) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Def. Tech., China
!     January 2013
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT(IN)   :: LDU,LDU2,LDVT,LDVT2,NL,NR,SQRE,SB
      INTEGER, INTENT(INOUT):: INFO, K
      LOGICAL ::               HFLAG
!     ..
!     .. Array Arguments ..
      INTEGER            COLTYP(*), IDX(*), IDXC(*), IDXP(*), &
                         IDXQ( * )
      DOUBLE PRECISION   D(*),DSIGMA(*),U( LDU,* ), U2( LDU2,* ),&
                         VT( LDVT,* ), VT2( LDVT2,* ), WORK( * ), &
                         Z(NL+NR+SB+SQRE,*),BETA(SB,SB),ALPHA(SB,SB)
!  ..
!  Purpose
!  =======
!  MDBDSVD2 merges the two sets of singular values together into a single
!  sorted set. Then it tries to deflate the size of the problem for the first
!  broken arrow matrix.
!
!  There are two ways in which deflation can occur:  when two or more
!  singular values are close together or if there is a tiny entry in the
!  Z vector.  For each such occurrence the order of the related secular
!  equation problem is reduced by one.
!
!  MDBDSVD2 is called from MDBDSVD1, and a similar routine is MDBDSVD22 which 
!  deflates a broken arrow matrix, and is called after MDBDSVD2. Different from
!  MDBDSVD2, MDBDSVD22 is much simpler and does not use the sparse structure of
!  singular vectors. Note that ALPHA and BETA are now upper and lower triangular
!  matrices respectively. 
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
!         = 0:  the lower block is an NR-by-NR square matrix.
!         = SB: the lower block is an NR-by-(NR+SB) rectangular matrix.
!
!         The upper band triangular matrix has N = NL+NR+SB rows and
!         M = N + SQRE >= N columns.
!
! SB    (input) INTEGER 
!         the semi-bandwith of this upper banded matrix, and the row dimension 
!         of added block row
!
!  K      (output) INTEGER
!         It equals to N - ND, where ND is the number of deflated singular values.
!         It also equals to SB-1 plus the order of the related secular equation, 
!         1 <= K <=N. The size of secular equation is K-SB+1
!
!  D      (input/output) DOUBLE PRECISION array, dimension( N )
!         On entry D contains the singular values of the two submatrices
!         to be combined.  On exit D contains the trailing (N-K) updated
!         singular values (those which were deflated) sorted into
!         increasing order.
!
!  Z      (output) DOUBLE PRECISION array, dimension( M,SB )
!         On exit Z contains the updating row vector in the secular
!         equation.
!
!  ALPHA  (input) DOUBLE PRECISION, dimension( SB,SB )
!         Contains the main diagonal elements associated with the added 
!         block row.
!
!  BETA   (input) DOUBLE PRECISION, dimension( SB,SB )
!         Contains the off-diagonal elements associated with the added
!         block row.
!
!  U      (input/output) DOUBLE PRECISION array, dimension( LDU,N )
!         On entry U contains the left singular vectors of two
!         submatrices in the two square blocks with corners at (1,1),
!         (NL, NL), and (NL+SB+1, NL+SB+1), (N,N).
!         On exit U contains the trailing (N-K) updated left singular
!         vectors (those which were deflated) in its last N-K columns.
!         U is updated only when the second kind deflation happens. 
!
!  LDU    (input) INTEGER
!         The leading dimension of the array U.  LDU >= N.
!
!  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)
!         On entry VT**T contains the right singular vectors of two
!         submatrices in the two square blocks with corners at (1,1),
!         (NL+SB, NL+SB), and (NL+SB+1, NL+SB+1), (M,M).
!         On exit VT**T contains the trailing (N-K) updated right singular
!         vectors (those which were deflated) in its last N-K columns.
!         In case SQRE = SB, the last SB rows of VT spans the right null
!         space.
!         VT is updated only when the second kind deflation happens. 
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT,  LDVT >= M.
!
!  DSIGMA (output) DOUBLE PRECISION array, dimension( N )
!         Contains a copy of the diagonal elements (K-SB singular values
!         and SB zeros) in the secular equation.
!
!  U2     (output) DOUBLE PRECISION array, dimension( LDU2,N )
!         Contains a copy of the first K-SB left singular vectors which
!         will be used by MDBDSVD3 in a matrix multiply (DGEMM) to solve
!         for the new left singular vectors. U2 is arranged into four
!         blocks. The first block contains a block column with SB vectors 
!         with nonzeros at rows from NL+1 to NL+SB, and zero everywhere 
!         else; the second block contains non-zero
!         entries only at and above NL; the third contains non-zero
!         entries only below NL+SB; and the fourth is dense, which is caused 
!         when deflating two multiple roots. 
!         The first SB-1 columns of U2 are not referenced. 
!
!  LDU2   (input) INTEGER
!         The leading dimension of the array U2.  LDU2 >= N.
!
!  VT2    (output) DOUBLE PRECISION array, dimension( LDVT2,N )
!         VT2**T contains a copy of the first K right singular vectors
!         which will be used by DLASD3 in a matrix multiply (DGEMM) to
!         solve for the new right singular vectors. VT2 is arranged into
!         three blocks. The first block contains SB rows that corresponds
!         to the special 0 diagonal elements in SIGMA; the second block
!         contains non-zeros only at and before NL+SB; the third block
!         contains non-zeros only at and after NL+SB+1.
!
!  LDVT2  (input) INTEGER
!         The leading dimension of the array VT2.  LDVT2 >= M.
!
!  IDXP   (workspace) INTEGER array dimension( N )
!         This will contain the permutation used to place deflated
!         values of D at the end of the array. On output IDXP( SB+1:K )
!         points to the nondeflated D-values and IDXP( K+1:N )
!         points to the deflated singular values.
!
!  IDX    (workspace) INTEGER array dimension( N )
!         This will contain the permutation used to sort the contents of
!         D into ascending order.
!
!  IDXC   (output) INTEGER array dimension( N )
!         This will contain the permutation used to arrange the columns
!         of the deflated U matrix into three groups:  the first group
!         contains non-zero entries only at and above NL, the second
!         contains non-zero entries only below NL+SB+1, and the third is
!         dense.
!
!  IDXQ   (input/output) INTEGER array dimension( N )
!         This contains the permutation which separately sorts the two
!         sub-problems in D into ascending order.  Note that entries in
!         the first half of this permutation must first be moved SB
!         positions backward; and entries in the second half
!         must first have NL+SB added to their values.
!
!  COLTYP (workspace/output) INTEGER array dimension( N )
!         As workspace, this will contain a label which will indicate
!         which of the following types a column in the U2 matrix or a
!         row in the VT2 matrix is:
!         1 : non-zero in the upper half only
!         2 : non-zero in the lower half only
!         3 : dense
!         4 : deflated
!
!         On exit, it is an array of dimension 4, with COLTYP(I) being
!         the dimension of the I-th type columns.
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension( * )
!         Compute Z through matrix-matrix multiplication, complete orthogonal 
!         factorization, DTZRZF, and transpose.
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!  Based on LAPACK routine dlasd2
!
!  Written by S.-G. Li, on Jun. 31th, 2013
!  =======================================
!
!     .. Parameters ..
      INTEGER            BSMLZ, DFK, IQ, IVT1
      DOUBLE PRECISION   ZERO, ONE, TWO, EIGHT
      PARAMETER          ( ZERO  = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0, &
                           EIGHT = 8.0D+0, BSMLZ = 500 )
!     ..
!     .. Local Arrays ..
      INTEGER            CTOT( 4 ), PSM( 4 )
!     ..
!     .. Local Scalars ..
      INTEGER            CT, I, IDXI, IDXJ, IDXJP, J, JP, JPREV, K2, M, &
                         N,NLP1,NLP2,NLPD,NRPS,IWF,IWD,IWE,ITAU, &
                         IWK,LWORK,ITZ1,ITZ2,ITAZ,IWKZ,IWE2,SB1,NLPD1
      DOUBLE PRECISION   C, EPS, HLFTOL, S, TAU, TOL
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
!
      IF( NL.LT.1 ) THEN
         INFO = -1
      ELSE IF( NR.LT.1 ) THEN
         INFO = -2
      ELSE IF( ( SQRE.NE.SB ) .AND. ( SQRE.NE.0 ) ) THEN
         INFO = -3
      END IF
!
      N = NL + NR + SB
      M = N + SQRE
      SB1 = SB + 1
      HFLAG = .TRUE. 
!
      IF( LDU.LT.N ) THEN
         INFO = -11
      ELSE IF( LDVT.LT.M ) THEN
         INFO = -13
      ELSE IF( LDU2.LT.N ) THEN
         INFO = -16
      ELSE IF( LDVT2.LT.M ) THEN
         INFO = -18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSVD2', -INFO )
         RETURN
      END IF
!
      NLP1 = NL + 1      ! starting row of added rows
      NLPD = NL + SB     ! dim of VT of left subproblem
      NLPD1 = NLPD + 1   ! starting row of right subproblem
      NRPS  = NR + SQRE  ! col dim of VT of right subproblem
      LWORK = M          ! SB if memory is limitted 
!
!     Generate the first part of the block vector Z; and move the singular
!     values in the first part of D SB positions backward.
!
      IWF = 1               ! front of Z
      IWD = IWF + SB        ! middle of Z
      IWE = IWD + NL        ! back of Z
      IQ  = IWF + SB*M      ! an orthogonal matrix in WORK
      ITAU = IQ + SB*SB     ! TAU used to construct Q
      IWK = ITAU + SB       ! workspace
!
      ! middle --> front
      CALL DGEMM('N','T',SB,SB,SB,ONE,VT(NLP1,NLP1),LDVT,ALPHA,SB, &
                  ZERO, Z(IWF,1), M )
      ! front --> middle
      CALL DGEMM('N','T',NL,SB,SB,ONE,VT(1,NLP1),LDVT,ALPHA,SB, &
                  ZERO, Z(IWD,1), M )
      ! back part
      CALL DGEMM('N','T',NRPS,SB,SB,ONE,VT(NLPD1,NLPD1),LDVT,BETA,SB, &
                  ZERO, Z(IWE,1), M )
      ! move D_1 backward
      DO I = NL, 1, -1
         D( I+SB ) = D( I )
         IDXQ( I+SB ) = IDXQ( I ) + SB  !used to permute Sigma_1 into increasing order
      END DO
      D( 1:SB ) = ZERO
!
!     QR factorization 
      CALL transp2d( Z, M, SB, WORK(IWF) )     ! WORK stores Z in row form
      CALL DGEQRF( SB,M,WORK(IWF),SB,WORK(ITAU),WORK(IWK),LWORK,INFO )
      CALL DLACPY( 'A',SB,SB,WORK(IWF),SB,WORK(IQ),SB )
!     WORK(IQ) stores Q
      CALL DORGQR( SB,SB,SB,WORK(IQ),SB,WORK(ITAU),WORK(IWK),LWORK,INFO )
      CALL transp2d( WORK(IWF),SB,M,Z )
!
!     Initialize some reference arrays
!
      DO I = SB1, NLPD
         COLTYP( I ) = 1
      END DO
      DO I = NLPD1, N
         COLTYP( I ) = 2
      END DO
!
!     Sort the singular values into increasing order
!
      DO I = NLPD1, N
         IDXQ( I ) = IDXQ( I ) + NLPD ! used to permute Sigma_2 into increasing order
      END DO
!
!     DSIGMA, IDXC, and the first SB columns of U2
!     are used as storage space.
!
      CALL DLASET( 'A', N, SB, ZERO, ZERO, U2, LDU2 )
      DO I = SB1, N
         DSIGMA( I ) = D( IDXQ( I ) )     ! Sigma_1 and Sigma_2 store two parts of interested singular 
!                                           values in one array, S_1 and S_2 are in increasing order respectively.
         U2( I,SB ) = Z( IDXQ(I),SB ) ! U2 stores the permuted middle part of Z 
         IDXC( I ) = COLTYP( IDXQ( I ) )  ! labels of columns of U and rows of VT
      END DO
!
      CALL DLAMRG( NL, NR, DSIGMA( SB1 ), 1, 1, IDX( SB1 ) )
!
!     Merge Sigma_1 and Sigma_2 into one ordered array. 
!     DSIGMA( IDX(i)+SB ) would be the i-th smallest singular value of the
!     two previous subproblems. 
!
      DO I = SB1, N
         IDXI = SB + IDX( I )
         D( I ) = DSIGMA( IDXI )       ! store the ordered svals back into the back of D, and
!                                        the whole D is in increasing order. 
         Z( I,SB ) = U2( IDXI,SB )     ! Z is the permuted, middle added rows
         COLTYP( I ) = IDXC( IDXI )    ! labels of cols of current U
      END DO
!
!     Calculate the allowable deflation tolerance
!
      EPS = DLAMCH( 'Epsilon' )
      TOL = MAXVAL( ABS( Z(SB:N, SB) ) )  ! the last col of Z
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
      K = SB         ! index of nondeflated D-values in IDXP
      K2 = N + 1     ! index of delfated D-values in IDXP
      DO 80 J = SB1, N      ! first few Z(j) are negligible
         IF( ABS( Z( J,SB ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
            K2 = K2 - 1
            IDXP( K2 ) = J   ! point to current column
            COLTYP( J ) = 4  ! label the column before the permutation of deflation
            ! If a column is labelled as deflation, the corresponding other added rows should also 
            ! be moved to back part. 
            IF( J.EQ.N ) &
               GO TO 120
         ELSE
            JPREV = J
            GO TO 90  ! If Z(J) is not small, check whether Z(J+1) is small or there are multiple svals.
         END IF
   80 CONTINUE
   90 CONTINUE
      J = JPREV
  100 CONTINUE
      J = J + 1
      IF( J.GT.N ) &
         GO TO 110
      IF( ABS( Z( J,SB ) ).LE.TOL ) THEN  ! Z(J+1) is small
!
!        Deflate due to small z component. 
!
         K2 = K2 - 1      
         IDXP( K2 ) = J   ! IDXP contains the deflated D-values
         COLTYP( J ) = 4
      ELSE
!
!        Check if singular values are close enough to allow deflation.
!
         IF( ABS( D( J )-D( JPREV ) ).LE.TOL ) THEN 
!
!           Deflation is possible.
!
            S = Z( JPREV,SB )
            C = Z( J,SB )
!
!           Find sqrt(a**2+b**2) without overflow or
!           destructive underflow.
!
            TAU = DLAPY2( C,S )
            C = C / TAU
            S = -S / TAU
            Z( J,SB ) = TAU
            Z( JPREV,SB ) = ZERO
!
!           Apply back the Givens rotation to the other added rows, and 
!           the left and right singular vector matrices. The inserted rows are
!           stored in Z as block column. 
!
            IDXJP = IDXQ( IDX( JPREV )+SB )
            IDXJ = IDXQ( IDX( J )+SB )
            IF( IDXJP.LE.NLPD ) THEN
               IDXJP = IDXJP - SB  ! Since we do not permute matrix U and VT, if that column belongs to U1, 
!                                    its index should be subtracted by SB. (We are assuming U1 has been moved
!                                    back by SB positions.)
            END IF
            IF( IDXJ.LE.NLPD ) THEN
               IDXJ = IDXJ - SB
            END IF
            CALL DROT( N, U( 1,IDXJP ), 1, U( 1,IDXJ ), 1, C, S )
            CALL DROT( M, VT( IDXJP,1 ), LDVT, VT( IDXJ,1 ), LDVT, C, &
                       S )
            CALL DROT( SB-1,Z( IDXJP,1 ),M,Z( IDXJ,1 ),M,C,S )
            IF( COLTYP( J ).NE.COLTYP( JPREV ) ) THEN
               COLTYP( J ) = 3   ! Col(J) is a dense vector
            END IF
            COLTYP( JPREV ) = 4  ! Col(JPREV) is a deflated vector
            K2 = K2 - 1 
            IDXP( K2 ) = JPREV
            JPREV = J
         ELSE   ! not 
            K = K + 1     ! the size of secular equation increase by one
            U2( K,SB ) = Z( JPREV,SB )   ! store current nondeflated Z(J)
            DSIGMA( K ) = D( JPREV )  ! store current nondeflated D(J)
            IDXP( K ) = JPREV         ! point the kth non-deflated sval
            JPREV = J                 ! move on to the next J
         END IF
      END IF
      GO TO 100   ! loop
  110 CONTINUE
!
!     Record the last singular value.
!
      K = K + 1
      U2( K,SB ) = Z( JPREV,SB )
      DSIGMA( K ) = D( JPREV )
      IDXP( K ) = JPREV
!
  120 CONTINUE
!
      IF( K .LT. BSMLZ ) THEN    ! Sparse structure
! 
!     Count up the total number of the various types of columns, then
!     form a permutation which positions the four column types into
!     four groups of uniform structure (although one or more of these
!     groups may be empty ).
!
         DO J = 1, 4
            CTOT( J ) = 0    ! total number of each type of column
         END DO
         DO J = SB1, N
            CT = COLTYP( J )
            CTOT( CT ) = CTOT( CT ) + 1
         END DO
!
!     PSM(*) = Position in SubMatrix ( of types 1 through 4 )
!
         PSM( 1 ) = SB1  ! start from the SB1-th column, first SB
!                          columns are an orthogonal matrix, Q
         PSM( 2 ) = SB1 + CTOT( 1 )       ! starting col for type 2
         PSM( 3 ) = PSM( 2 ) + CTOT( 2 )  ! starting col for type 3
         PSM( 4 ) = PSM( 3 ) + CTOT( 3 )  ! starting col for type 4
         DFK = PSM( 4 )
!
!     Fill out the IDXC array so that the permutation which it induces
!     will place all type-1 columns first, all type-2 columns next,
!     then all type-3's, and finally all type-4's, starting from the
!     SB1-th column. This applies similarly to the rows of VT.
!
         DO J = SB1, N
            JP = IDXP( J )        ! JP is the original column
            CT = COLTYP( JP )
            IDXC( PSM( CT ) ) = J ! Uisng IDXP IDXC can point to the original col
            PSM( CT ) = PSM( CT ) + 1 ! Each type increases by one
         END DO
!
!     Sort the singular values and corresponding singular vectors into
!     DSIGMA, U2, and VT2 respectively.  The singular values/vectors
!     which were not deflated go into the first K slots of DSIGMA, U2,
!     and VT2 respectively, while those which were deflated go into the
!     last N - K slots, except that the first SB columns/rows will be 
!     treated separately.
!
         DO J = SB1, N
            JP = IDXP( J )
            DSIGMA( J ) = D( JP )  ! DSIGMA: first part stores the nondeflated svals in increasing order
!                                    the second part stores the deflated svals in increasing order.
            IDXJ = IDXQ( IDX( IDXP( IDXC( J ) ) )+SB ) ! IDXJ is the original index in U and VT
            CALL DCOPY( SB-1, Z(IDXJ,1), M, U2(J,1), LDU2 )  ! The rows of Z will be permuted
            IF( IDXJ.LE.NLPD ) THEN
               IDXJ = IDXJ - SB
            END IF
            CALL DCOPY( N, U( 1, IDXJ ), 1, U2( 1, J ), 1 )
            CALL DCOPY( M, VT( IDXJ, 1 ), LDVT, VT2( J, 1 ), LDVT2 )
         END DO
         CALL DLACPY( 'A',N-SB,SB,U2(SB1,1),LDU2,Z(SB1,1),M )  ! The first SB columns of U2 store appended columns
!
      ELSE        ! Do not use the sparse structure. 
!
!     Sort the singular values and corresponding singular vectors into
!     DSIGMA, U2, and VT2 respectively. 
!     Modify the columns (rows) of U2 (VT2) starting from the SB-th one.
!     The first SB rows and columns are modifed separately below.
!
         HFLAG = .FALSE.
         DO J = SB1, N
            JP = IDXP( J )
            DSIGMA( J ) = D( JP )
            IDXJ = IDXQ( IDX( IDXP( J ) )+SB )
            CALL DCOPY( SB-1, Z(IDXJ,1), M, U2(J,1), LDU2 )  ! The rows of Z will be permuted
            IF( IDXJ.LE.NLPD ) THEN
               IDXJ = IDXJ - SB
            END IF
!   Copy N elements of the first N columns of U to U2
            CALL DCOPY( N, U( 1, IDXJ ), 1, U2( 1, J ), 1 )
!   Copy M elements of the first N rows of VT to VT2
            CALL DCOPY( M, VT( IDXJ, 1 ), LDVT, VT2( J, 1 ), LDVT2 )
         END DO
         CALL DLACPY( 'A',N-SB,SB,U2(SB1,1),LDU2,Z(SB1,1),M )  ! The first SB columns of U2 store appended columns
      END IF
!
!     Determine DSIGMA( SB ), DSIGMA( SB1 ) and Z( SB )
!
      DSIGMA( SB ) = ZERO
      HLFTOL = TOL / TWO
      IF( ABS( DSIGMA( SB1 ) ).LE.HLFTOL ) &
         DSIGMA( SB1 ) = HLFTOL
!
!     Determine the first SB columns of U2, U2(NLP1:NL+SB, 1:SB) = Q.
!     The first SB columns of U (may not needed)
      CALL DLASET( 'A', N, SB, ZERO, ZERO, U, LDU )
      CALL DLACPY( 'A', SB, SB, WORK(IQ), SB, U(NLP1,1), LDU )
!     The first SB columns of U2
      CALL DLACPY( 'A', SB, SB, WORK(IQ), SB, U2(NLP1,1), LDU2 )
!
!     If M > N, use orthogonal transformation to eliminate the last part of block Z, 
!     and apply the transformation to the right singular vector matrix VT. 
!
      IF( M.GT.N ) THEN
         ITZ1 = IQ    ! the front part of Z, can use IQ only
         ITZ2 = ITZ1 + SB*SB  ! the back part of Z
         ITAZ = ITZ2 + SB*SB  ! TAU for dtzrzf
         IWKZ = ITAZ + SB     ! WORK
         IWE2 = IWF  + SB*N
         LWORK= M*SB
         CALL DLACPY( 'A', SB, SB, WORK(IWF), SB, WORK(ITZ1), SB )
         CALL DLACPY( 'A', SB, SB, WORK(IWE2), SB, WORK(ITZ2), SB )
         CALL DLASET( 'Lower',SB-1,2*SB,ZERO,ZERO, WORK(ITZ1+1), SB )
         CALL DTZRZF( SB,2*SB,WORK(ITZ1),SB,WORK(ITAZ),WORK(IWKZ),LWORK,INFO )
!
         ! update the front part of block vector Z
         CALL DLACPY( 'U', SB, SB, WORK(ITZ1), SB, WORK(IWF), SB )
         CALL DLASET( 'A',SB,SB,ZERO,ZERO,Z(IWE+NR,1),M ) ! set the bottom part zero
         CALL transp2di( SB,SB,WORK(IWF),SB,Z,M,INFO )  ! Z is stored in column form
      ELSE 
         IF( ABS( Z(SB,SB) ).LE.TOL ) THEN
            Z(SB,SB) = TOL
         END IF
      END IF
!
!     The first and last SB rows of VT is if SQRE == 0, VT2(1:SB,:) = VT(1:SB,:); 
!     Otherwise, VT2(1:SB,:) is computed from DORMRZ. 
!
!     The first and last SB rows of VT2 and VT
      IF( M.GT.N ) THEN
         IVT1 = IWKZ
         IWK = IVT1+ 2*SB*M
         LWORK= M*SB
         CALL DLACPY( 'A', SB, M, VT(NLP1,1), LDVT, WORK(IVT1),2*SB )
         CALL DLACPY( 'A', SB, M, VT(N+1,1), LDVT, WORK(IVT1+SB),2*SB )
         CALL DORMRZ( 'L', 'N', 2*SB, M, SB, SB, WORK(ITZ1), SB, &
                  WORK(ITAZ),WORK(IVT1),2*SB,WORK(IWK),LWORK,INFO )
!        first SB rows of VT2 and VT
         CALL DLACPY( 'A', SB, M, WORK(IVT1), 2*SB, VT2(1,1), LDVT2 )
         CALL DLACPY( 'A', SB, M, WORK(IVT1), 2*SB, VT(1,1), LDVT )
!        last SB rows of VT and VT2
         CALL DLACPY( 'A', SB, M, WORK(IVT1+SB), 2*SB, VT(N+1,1), LDVT )
         CALL DLACPY( 'A', SB, M, WORK(IVT1+SB), 2*SB, VT2(N+1,1), LDVT2 )
      ELSE
         CALL DLACPY( 'A', SB, M, VT(NLP1,1), LDVT, VT2(1,1), LDVT2 )
         CALL DLACPY( 'A', SB, M, VT(NLP1,1), LDVT, VT(1,1), LDVT )
      END IF
!
      IF( .NOT. HFLAG ) THEN
!
!     The deflated singular values and all singular vectors go
!     into the back of D, U, and V respectively.
!
         CALL DCOPY( N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 )
         CALL DLACPY( 'A', N, N-SB, U2( 1, SB1 ), LDU2, U( 1, SB1 ), &
                      LDU )
         CALL DLACPY( 'A', N-SB, M, VT2( SB1,1 ), LDVT2, VT( SB1,1 ), &
                      LDVT )
!
      ELSE IF( N.GT.K ) THEN
!
!     The deflated singular values and their corresponding vectors go
!     into the back of D, U, and V respectively.
!
!     Whether store the orthogonal matrix in the first SB columns of U
!     is clear right now. We should store it in U. A new question is do
!     we need to store them in U2 and VT2 ?
!
         CALL DCOPY( N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 )
         CALL DLACPY( 'A', N, N-K, U2( 1, K+1 ), LDU2, U( 1, K+1 ), &
                      LDU )
         CALL DLACPY( 'A', N-K, M, VT2( K+1, 1 ), LDVT2, VT( K+1, 1 ), &
                      LDVT )
      END IF
!
!     Copy CTOT into COLTYP for referencing in DLASD3.
!
      DO J = 1, 4
         COLTYP( J ) = CTOT( J )
      END DO
!
!     Set the upper triangular part of Z to ZERO
      CALL DLASET( 'Upper',SB-1, SB-1,ZERO,ZERO,Z(1,2),M )
!
      RETURN
!
!     End of MDBDSVD2
!
      END

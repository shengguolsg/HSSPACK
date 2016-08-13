      SUBROUTINE MDBDSVD0( N,SQRE,SMBD,A,LDA,U,LDU,S,VT,LDVT,SMLSIZ, &
                           IWORK,WORK,INFO )
        IMPLICIT NONE
!
!  -- HSSPACK computational routine (version 0.0.1) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Defense Tech., China.
!     January 2013
!
!  .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: N, SQRE, SMBD, LDA, LDU, LDVT, SMLSIZ
      INTEGER, INTENT(INOUT) :: INFO
!  ..
!  .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA,* ),U( LDU,* ), S( * ), VT( LDVT,* ), &
                         WORK( * )
!  ..
!
!  Purpose
!  =======
!
!  Using a divide and conquer approach, MDBDSVD0 computes the singular
!  value decomposition (SVD) of a real upper banded N-by-M matrix
!  B with semi-bandwidth SMBD, where M = N + SQRE and SQRE=ZERO or SMBD.
!  The algorithm computes orthogonal matrices U and VT such that
!  B = U * S * VT. The singular values are stored in S. Matrix B is stored
!  in A in sparse storage form. A is a (SMBD+1)-by-M matrix. LDA >= SMBD+1. 
!
!  This routine is similar to DLASD0, which constructs a computation tree
!  and solves the problem by following this tree from bottom to top level
!  by level. In this code, we assume JOBU='A', i.e. the left and right 
!  singular vectors will be stored in U and VT respectively. 
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         On entry, the row dimension of the upper banded matrix.
!         This is also the dimension of the main diagonal array S.
!
!  SQRE   (input) INTEGER
!         Specifies the column dimension of the banded matrix.
!         = 0:    The banded matrix has column dimension M = N;
!         = SMBD: The banded matrix has column dimension M = N+SMBD;
!
!  SMBD   (input) INTEGER
!         Specifies the semi-bandwidth of the upper banded matrix, 
!         and this number does not include the main diagonal entry. 
!         For example, the SMBD of a bidiagonal matrix would be 1. 
! 
!  A      (input/output) DOUBLE PRECISION array, dimension( LDA,M )
!         On entry A contains the original banded matrix.
!         On exit A, if INFO = 0, contains the its left singular vectors.
!
!  U      (input/output) DOUBLE PRECISION array, dimension( LDU,N )
!         Contains the left singular vector matrix of A. If JOBU='O',
!         U will not be referred. 
!
!  LDU    (input) INTEGER
!         On entry, leading dimension of U.
!
!  S      (output) DOUBLE PRECISION array, dimension( N )
!         On exit, S contains the singular values of A.
!
!  VT     (output) DOUBLE PRECISION array, dimension( LDVT,M )
!         On entry, VT contains the right singular vector matrix of A.
! 
!  LDVT   (input) INTEGER, LDVT >= M
!         On entry, leading dimension of VT.
!
!  SMLSIZ (input) INTEGER
!         On entry, maximum size of the subproblems at the
!         bottom of the computation tree. It usually equals to 25, defined 
!         in MDBDSDD. 
!
!  IWORK  (workspace) INTEGER work array.
!         with dimension 8*N
!
!  WORK   (workspace) DOUBLE PRECISION work array.
!         Dimension must be at least (3 * M**2 + 2 * M)
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!   Written by S.-G. Li, NUDT, Jan. 16th, 2013
!  ===========================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK, &
                         J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR, &
                         NL, NLF, NLVL, NR, NRF, SQREI, ICR, ICR2, NLD, LWORK, &
                         NRD
      DOUBLE PRECISION   BlkALPHA(SMBD,SMBD), BlkBETA(SMBD,SMBD)
!                        BlkALPHA: upper triangular in main diagonal; 
!                        BlkBETA: lower triangular in off diagonal part;
!     ..
!     .. External Subroutines ..
      EXTERNAL           MDBDSVD1, DGESVD, MDBDSVDT, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( (N.LT.0) .OR. (N.LT.SMBD) ) THEN
         INFO = -1
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.SMBD ) ) THEN
         INFO = -2
      ELSE IF( SMBD.LT.1 ) THEN
         INFO = -3
      END IF
!
      M = N + SQRE
!
      IF( LDU.LT.N ) THEN
         INFO = -6
      ELSE IF( LDVT.LT.M ) THEN
         INFO = -8
      ELSE IF( SMLSIZ.LT.3 ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSVD0', -INFO )
         RETURN
      END IF
!
!     If the input matrix is too small, call DGESVD to find the SVD.
!
      IF( N.LE.SMLSIZ ) THEN
         LWORK = 5*N
         CALL DGESVD('A','A',N,M,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO )
         RETURN
      END IF
!
!     Set up the computation tree.
!
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
      IDXQ  = NDIMR + N
      IWK   = IDXQ  + N
      CALL MDBDSVDT( N, SMBD, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), &
                   IWORK( NDIMR ), SMLSIZ )
!
!     For the nodes at the bottom level of the tree, solve
!     their subproblems by DGESVD.
!
      NDB1 = ( ND+1 ) / 2  ! the left most leaf node
      NCC = 0
      DO I = NDB1, ND      ! leaf nodes, each of which have two subproblems
!
!     IC : center row of each node
!     NL : number of rows of left  subproblem
!     NR : number of rows of right subproblem
!     NLF: starting row of the left  subproblem
!     NRF: starting row of the right subproblem
!
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )   ! row dimension of left  subproblem
         NR = IWORK( NDIMR+I1 )   ! row dimension of right subproblem
         NLF = IC - NL
         NRF = IC + SMBD
         NLD = NL + SMBD
         LWORK = 5*NL      ! choose a bigger one for efficiency
!
         ! copy the upper banded matrix to U in full storage form
         CALL DSB2FLC( 'U',SMBD+1,SMBD,U(NLF,NLF),LDU,A(SMBD+1,NLF),LDA,&
                       INFO )
         CALL DGESVD( 'A','A',NL,NLD,A(NLF,NLF),LDA,S(NLF),U(NLF,NLF),LDU, &
                    VT(NLF,NLF),LDVT,WORK,LWORK,INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         ITEMP = IDXQ + NLF - 2  ! IDXQ would be the first entry
         DO J = 1, NL
            IWORK( ITEMP+J ) = J
         END DO
         IF( I.EQ.ND ) THEN   ! for right subproblem
            SQREI = SQRE
         ELSE
            SQREI = SMBD
         END IF
         NRD = NR + SQREI
         LWORK = 5*NR         ! choose a bigger one for efficiency 
         CALL DGESVD( 'A','A',NR,M,A(NRF,NRF),LDA,S(NRF),U(NRF,NRF),LDU, &
                    VT(NRF,NRF),LDVT,WORK,LWORK,INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         ITEMP = IDXQ + IC + SMBD-1   ! starting row of right subproblem
         DO J = 1, NR
            IWORK( ITEMP+J-1 ) = J
         END DO 
      END DO  ! ( I )
!
!     Now conquer each subproblem bottom-up.
!
      DO LVL = NLVL, 1, -1
!
!        Find the first node LF and last node LL on the
!        current level LVL.
!
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )  ! row dimension of left
            NR = IWORK( NDIMR+IM1 )  ! row dimension of right
            NLF = IC - NL
            IF( ( SQRE.EQ.0 ) .AND. ( I.EQ.LL ) ) THEN   ! for right problem
               SQREI = SQRE
            ELSE
               SQREI = SMBD
            END IF
            IDXQC = IDXQ + NLF - 1   ! IDXQC would be the first entry
            ICR  = IC+SMBD-1
            ICR2 = ICR+SMBD
            BlkALPHA = A( IC:ICR,IC:ICR )
            BlkBETA  = A( IC:ICR,ICR+1:ICR2 )
            CALL MDBDSVD1( NL, NR, SMBD, SQREI, S(NLF), BlkALPHA, BlkBETA, &
                         U( NLF, NLF ), LDU, VT( NLF, NLF ), LDVT, &
                         IWORK( IDXQC ), IWORK( IWK ), WORK, INFO )
            IF( INFO.NE.0 ) THEN
               RETURN
            END IF
         END DO ! ( I )
      END DO ! ( LVL )
!
      RETURN
!
!     End of MDBDSVD0
!
    END SUBROUTINE MDBDSVD0

      SUBROUTINE MDBDSVDT( N, SMBD, LVL, ND, INODE, NDIML, NDIMR, MSUB )
!
!  -- HSSPACK computational routine (version 0.0.1) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Defense Tech., China.
!     January 2013
!
!     .. Scalar Arguments ..
      INTEGER, INTENT(IN)    :: N, SMBD
      INTEGER, INTENT(INOUT) :: LVL, ND, MSUB 
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT(INOUT) :: INODE( * ), NDIML( * ), NDIMR( * )
!     ..
!
!  Purpose
!  =======
!
!  MDBDSVDT creates a tree of subproblems for banded divide and conquer algorithm
!  of singular value problems.
!
!  Arguments
!  =========
!
!   N      (input) INTEGER
!          On entry, the row dimension of this upper banded matrix, and 
!          we assume it is a square matrix. 
!
!   SMBD   (input) INTEGER
!          The semi-bandwidth of the upper banded matrix, which does not 
!          count the main diagonal element. For  bidiagonal matrix SMBD equals
!          to 1. 
! 
!   LVL    (output) INTEGER
!          On exit, the number of levels on the computation tree.
!
!   ND     (output) INTEGER
!          On exit, the number of nodes on the tree.
!
!   INODE  (output) INTEGER array, dimension ( N )
!          On exit, centers of subproblems. 
!
!   NDIML  (output) INTEGER array, dimension ( N )
!          On exit, row dimensions of left children. The row dimension
!          of left child will usually larger than that of right child 
!          by about SMBD.
!
!   NDIMR  (output) INTEGER array, dimension ( N )
!          On exit, row dimensions of right children.
!
!   MSUB   (input) INTEGER
!          On entry, the maximum row dimension each subproblem at the
!          bottom of the tree can be of. MSUB should be at least five
!          times larger than SMBD. 
!
!  Further Details
!  ===============
!  Based on LAPACK routine DLASDT
!
!  Written by S.-G. Li, on Jan. 16th, 2013
!  =======================================
!
!     .. Parameters ..
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IL, IR, LLST, MAXN, NCRNT, NLVL
      DOUBLE PRECISION   TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, LOG, MAX
!     ..
!     .. Executable Statements ..
!
!     Find the number of levels on the tree.
!
      IF( MSUB .LT. 3*SMBD ) THEN
         MSUB = 3*SMBD
      END IF
      MAXN = MAX( 1, N )
      TEMP = LOG( DBLE( MAXN ) / DBLE( MSUB+1 ) ) / LOG( TWO )
      LVL = INT( TEMP )
!
      I = N / 2
      INODE( 1 ) = I + 1
      NDIML( 1 ) = I
      NDIMR( 1 ) = N - I - SMBD
      IL = 0
      IR = 1
      LLST = 1
      DO NLVL = 1, LVL - 1
!
!        Constructing the tree at (NLVL+1)-st level. The number of
!        nodes created on this level is LLST * 2.
!
         DO I = 0, LLST - 1
            IL = IL + 2
            IR = IR + 2
            NCRNT = LLST + I
            NDIML( IL ) = NDIML( NCRNT ) / 2
            NDIMR( IL ) = NDIML( NCRNT ) - NDIML( IL ) - SMBD
            INODE( IL ) = INODE( NCRNT ) - NDIMR( IL ) - SMBD
            NDIML( IR ) = NDIMR( NCRNT ) / 2
            NDIMR( IR ) = NDIMR( NCRNT ) - NDIML( IR ) - SMBD
            INODE( IR ) = INODE( NCRNT ) + NDIML( IR ) + SMBD
         END DO  ! ( I )
         LLST = LLST*2
      END DO  ! ( NLVL )
      ND = LLST*2 - 1
!
      RETURN
!
!     End of MBDSVDT
!
    END SUBROUTINE MDBDSVDT

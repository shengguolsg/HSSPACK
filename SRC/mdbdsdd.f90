     SUBROUTINE MDBDSDD( JOBZ,N,KD,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,IWORK,INFO)
!
       IMPLICIT NONE
!
!  -- HSSPACK computational routine (version 0.0.1) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Def. Tech., China.
!     January 2013
!
!  .. Scalar Arguments ..
      CHARACTER, INTENT(IN)    :: JOBZ
      INTEGER,   INTENT(IN)    :: N,LDA,LDU,LDVT,LWORK,KD
      INTEGER,   INTENT(INOUT) :: INFO
!  ..
!  .. Array Arguments ..
      DOUBLE PRECISION, INTENT(INOUT) :: A( LDA,* ), S( * ), U( LDU,* ), &
                                         VT( LDVT,* ), WORK( * ), IWORK( * )
!  ..
!  Purpose
!  =======
!
!  MDBDSDD computes the SVD of a general dense N-by-N matrix. 
!  It computes both the singular values and singular vectors, and it
!  uses HSS techniques to compute singular vectors when the secular 
!  equations are large. This algorithm uses the HSS techniques, see 
!  S.-G. Li's paper 'banded DC algorithm'.
!
!  This subroutine first reduces a dense matrix into banded form, and 
!  then it calls mdbdsvd to compute the SVD of an upper banded matrix.
!  The singular vector matrices will also be updated inside mdbdsvd. 
! 
! Arguments
! ==========
! JOBZ (in) CHARACTER*1
!          Specifies options for computing all or part of the matrix U:
!          = 'A':  all N columns of U and all N rows of V**T are
!                  returned in the arrays U and VT;
!          = 'O':  A is overwritten with the first N columns
!                  of U (the left singular vectors, stored
!                  columnwise) if M >= N;
!                  A is overwritten with the first M rows
!                  of V**T (the right singular vectors, stored
!                  rowwise) otherwise.
!          = 'N':  no columns of U or rows of V**T are computed.
! 
! N (in) INTEGER
!        The order of the input matrix A.  N >= 0.
!
! KD (in) INTEGER
!         The semi-bandwith of desired upper banded matrix.
! 
! A (inout) DOUBLE PRECISION array, DIMENSION( LDA,N )
!          A is DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N matrix A.
!          On exit,
!          if JOBZ = 'O',  A is overwritten with the first N columns
!                          of U (the left singular vectors, stored
!                          columnwise) if M >= N;
!                          A is overwritten with the first M rows
!                          of V**T (the right singular vectors, stored
!                          rowwise) otherwise.
!          if JOBZ .ne. 'O', the contents of A are destroyed.
!
! LDA (in) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N)
!
! S (out)  DOUBLE PRECISION array, dimension ( N )
!          The singular values of A, sorted so that S(i) >= S(i+1)
!
! U (out)  DOUBLE PRECISION array, dimension( LDU,UCOL )
!          Store the left singular vector matrix. 
!
! LDU (in) INTEGER
!          The leading dimension of the array U.  LDU >= N
!
! VT (out) DOUBLE PRECISION array, dimension( LDVT,N )
!          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
!          N-by-N orthogonal matrix V**T;
!          if JOBZ = 'S', VT contains the first min( M,N ) rows of
!          V**T( the right singular vectors, stored rowwise );
!          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
! 
! LDVT (in) INTEGER
!          The leading dimension of the array VT.  LDVT >= N.
!
! WORK (inout) DOUBLE PRECISION array, dimension( MAX(1,LWORK) )
!          WORK is DOUBLE PRECISION array, dimension( MAX(1,LWORK) )
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!
! LWORK (in) INTEGER, The dimension of the array WORK.
!          If JOBZ = 'N',
!            LWORK >= 3*min(M,N) + max( max(M,N),7*min(M,N) )
!          If JOBZ = 'O',
!            LWORK >= 3*min(M,N) + 
!                     max( max(M,N),5*min(M,N)*min(M,N)+4*min(M,N) )
!          If JOBZ = 'S' or 'A'
!            LWORK >= 3*min(M,N) +
!                     max( max(M,N),4*min(M,N)*min(M,N)+4*min(M,N) )
!          For good performance, LWORK should generally be larger.
!          If LWORK = -1 but other input arguments are legal, WORK( 1 )
!          returns the optimal LWORK.
!          The required size of WORK is not clear yet. LWORK at least should be 
!          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) since we use DGESVD to
!          the basic problems. 
!
! IWORK (inout) INTEGER array, dimension( ceiling(8 * (N/KD) ) )
!          Used to store the computation tree of divide-and-conquer algorithm.
!
! INFO (out) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if MDBDSVD did not converge, updating process failed.
!
! ========
! Written by S.-G. Li, on Jan. 14th, 2013
!
! =======================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ISCL, SMLSIZ, M
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, SMLNUM
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASCL, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters..
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT. 0 ) THEN
         INFO = -3 
      ELSE IF( LDA.LT.1 ) THEN
         INFO = -4
      ELSE IF( LDU.LT. 1 ) THEN
         INFO = -7
      ELSE IF( LDVT.LT.1 ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSDD', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      SMLSIZ = 25
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, N, A, LDA, IERR )
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, N, A, LDA, IERR )
      END IF
!
!     Compute the problem in a divide-and-conquer way
!
      CALL MDBDSVD( JOBZ,N,KD,A,LDA,S,U,LDU,VT,LDVT,SMLSIZ,WORK,LWORK,IWORK,INFO )
!
!     Undo scaling if necessary
!
      IF( ISCL.EQ.1 ) THEN
         IF( ANRM.GT.BIGNUM )  &
              CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, &
                        IERR )
         IF( ANRM.LT.SMLNUM )  &
              CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, &
                         IERR )
      END IF
!
      RETURN
!
    END SUBROUTINE MDBDSDD

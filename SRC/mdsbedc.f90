!
!  Definition:
!  ===========
!
!       SUBROUTINE MDSBEDC( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,
!                          LWORK, IWORK, LIWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> MDSBEDC computes all the eigenvalues and, optionally, eigenvectors of
!> a real symmetric band matrix A. If eigenvectors are desired, it uses
!> a divide and conquer algorithm.
!>
!> It is based on the LAPACK routine DSBEVD without bidiagonalization. 
!> Now it is designed for upper triangular matrix. 
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER!1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB, N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD + 1.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!>          eigenvectors of the matrix A, with the i-th column of Z
!>          holding the eigenvector associated with W(i).
!>          If JOBZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array,
!>                                         dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          IF N <= 1,                LWORK must be at least 1.
!>          If JOBZ  = 'N' and N > 2, LWORK must be at least 2*N.
!>          If JOBZ  = 'V' and N > 2, LWORK must be at least
!>                         ( 1 + 5*N + 2*N**2 ).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK and IWORK
!>          arrays, returns these values as the first entries of the WORK
!>          and IWORK arrays, and no error message related to LWORK or
!>          LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.
!>          If JOBZ  = 'V' and N > 2, LIWORK must be at least 3 + 5*N.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK and IWORK arrays, and no error message related to
!>          LWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the algorithm failed to converge; 
!>                off-diagonal elements of an intermediate tridiagonal
!>                form did not converge to zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Nati. Univ. of Def. Tech.
!> \author Shengguo Li
!
!> \date October 2014
!
!> \ingroup doubleOTHEReigen
!
!  =====================================================================
      SUBROUTINE MDSBEDC( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, &
                        LWORK, IWORK, LIWORK, INFO )
!
!  -- HSSPACK driver routine (version 1.0.0) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Def. Tech. --
!     Oct 2014
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, INDE, INDWK2, INDWRK, ISCALE, LIWMIN, &
                         LLWRK2, LWMIN, MINSIZE, MaxBand
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, &
                         SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSB
      EXTERNAL           LSAME, DLAMCH, DLANSB
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASCL, DSBTRD, DSCAL, DSTEDC, &
                         DSTERF, XERBLA, DSBEVD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
!
      MINSIZE = 60 
      MaxBand = 16
      INFO = 0
      IF( N.LE.1 ) THEN
         LIWMIN = 1
         LWMIN = 1
      ELSE
         IF( WANTZ ) THEN
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 5*N + 2*N**2
         ELSE
            LIWMIN = 1
            LWMIN = 2*N
         END IF
      END IF
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -9
      END IF
!
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
!
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDSBEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
        RETURN
!
      IF( N.EQ.1 ) THEN
         W( 1 ) = AB( 1, 1 )
         IF( WANTZ ) &
            Z( 1, 1 ) = ONE
         RETURN
      END IF
!
      IF( N.LE. MINSIZE .OR. KD .GT. MaxBand ) THEN
         CALL DSBEVD( JOBZ,UPLO,N,KD,AB,LDAB,W,Z,LDZ,WORK,LWORK,IWORK,LIWORK,INFO )
         RETURN  
      END IF

!
!     Get machine constants.
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
      ANRM = DLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            CALL DLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         ELSE
            CALL DLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         END IF
      END IF
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, call MDSBEVD.
!
      IF( .NOT.WANTZ ) THEN
!
!        Call DSBTRD to reduce symmetric band matrix to tridiagonal form.
         INDE = 1
         INDWRK = INDE + N
         CALL DSBTRD( JOBZ,UPLO,N,KD,AB,LDAB, W, WORK( INDE ), Z, LDZ, &
              WORK( INDWRK ), IINFO )

         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
!
!        For eigenvectors, call MDSBEVD
         INDWRK = 1
         INDWK2 = INDWRK + N*N
         LLWRK2 = LWORK - INDWK2 + 1
         CALL MDSBEVD('I','U',N,KD,AB,LDAB,W, WORK( INDWRK ), N, &
                      WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO )
         CALL DLACPY( 'A', N, N, WORK( INDWRK ), N, Z, LDZ )
!
      END IF  ! (WANTZ)
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF( ISCALE.EQ.1 ) &
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
!
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
!
!     End of MDSBEDC
!
    END SUBROUTINE MDSBEDC

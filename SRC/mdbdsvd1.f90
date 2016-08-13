      SUBROUTINE MDBDSVD1( NL,NR,SB,SQRE,D,ALPHA,BETA,U,LDU,VT,LDVT,&
                      IDXQ,IWORK,WORK,INFO,NLF,A,AB,LDAB,nod )
!
        USE TestRoutine
        use BasicMM
        IMPLICIT NONE
!
!  -- HSSPACK computational routine (version 1.0.0) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Def. Tech., China,
!     January 2013.
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDU, LDVT, NL, NR, SB, SQRE, NLF,LDAB,nod
      DOUBLE PRECISION   ALPHA(SB,SB), BETA(SB,SB)
!     ..
!     .. Array Arguments ..
      INTEGER            IDXQ(*), IWORK(*)
      DOUBLE PRECISION   D(*),U(LDU,*),VT(LDVT,*),WORK(*),A(*),&
           AB(LDAB,*)
!     ..
!
!  Purpose
!  =======
!  MDBDSVD1 computes the SVD of an X-shape banded N-by-M matrix B,
!  where N = NL + NR + SB and M = N + SQRE. MDBDSVD1 is called from 
!  MDBDSVD.
!
!  MDBDSVD1 computes the SVD as follows:
!
!                ( D1(in)    0    0       0 )
!    B = U(in) * (   Z1**T  Aa   Z2**T   Bb ) * VT(in)
!                (   0       0   D2(in)   0 )
!
!      = U(out) * ( D(out) 0) * VT(out)
!
!  where Z**T = (Z1**T Aa Z2**T Bb) = u**T VT**T, and u is a block vector of dimension M
!  with ALPHA and BETA near the center and zeros elsewhere; and the 
!  block entry Bb is empty if SQRE = 0.
!
!  The left singular vectors of the original matrix are stored in U, and
!  the transpose of the right singular vectors are stored in VT, and the
!  singular values are in D.  The algorithm consists of three stages:
!
!     The first stage consists of deflating the size of the problem
!     when there are multiple singular values or when there are zeros in
!     the rows of Z. Z would be first reduced to upper triangular form by
!     using for example QR factorization. For each such occurence the dimension of the
!     secular equation problem is reduced by one.  This stage is
!     performed by the routine MDBDSVD2. After deflating one row of Z, the
!     other rows of Z should be modified accordingly. 
!
!     The second stage consists of calculating the updated singular values. 
!     This is done by finding the square roots of the roots of the secular 
!     equation via the routine DLASD4 (as called by DBDSVD3). This routine 
!     also calculates the singular vectors of the current problem.
!
!     The final stage consists of computing the updated singular vectors
!     directly using the updated singular values.  The singular vectors
!     for the current problem are multiplied with the singular vectors
!     from the overall problem. When the size of a secular equation is 
!     larger than 2500, the HSS matrix techniques would be used. 
!
!  The above three stages should be performed for each row of Z, which can 
!  be done similarly. After deflating the bottom row of Z, we have several
!  updating SVD problems which is deflated by MDBDSVD22 and updated by MDBDSVD4. 
! 
!  Arguments
!  =========
!  NL     (input) INTEGER
!         The row dimension of the upper block.  NL >= 1.
!
!  NR     (input) INTEGER
!         The row dimension of the lower block.  NR >= 1.
!
!  SB     (input) INTEGER
!         The semi-bandwith of this upper banded matrix
!  
!  SQRE   (input) INTEGER
!         = 0: the lower block is an NR-by-NR square matrix.
!         = SB: the lower block is an NR-by-(NR+SB) rectangular matrix.
!
!         The upper banded matrix has row dimension N = NL + NR + SB,
!         and column dimension M = N + SQRE.
!
!  D      (input/output) DOUBLE PRECISION array, dimension (N = NL+NR+SB).
!         On entry D(1:NL) contains the singular values of the
!         upper block; and D(NL+SB+1:N) contains the singular values of
!         the lower block. On exit D(1:N) contains the singular values
!         of B.
!
!  ALPHA  (input/output) DOUBLE PRECISION
!         Contains the diagonal element associated with the inserted block row.
!
!  BETA   (input/output) DOUBLE PRECISION
!         Contains the off-diagonal element associated with the inserted block
!         row.
!
!  U      (input/output) DOUBLE PRECISION array, dimension(LDU,N)
!         On entry U(1:NL, 1:NL) contains the left singular vectors of
!         the upper block; U(NL+SB+1:N, NL+SB+1:N) contains the left singular
!         vectors of the lower block. On exit U contains the left
!         singular vectors of the merged subproblems
!
!  LDU    (input) INTEGER
!         The leading dimension of the array U.  LDU >= max( 1, N ).
!
!  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)
!         where M = N + SQRE.
!         On entry VT(1:NL+SB,1:NL+SB)**T contains the right singular
!         vectors of the upper block; VT(NL+SB+1:M,NL+SB+1:M)**T contains
!         the right singular vectors of the lower block. On exit
!         VT**T contains the right singular vectors of the upper banded matrix.
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT.  LDVT >= max( 1,M ).
!
!  IDXQ  (output) INTEGER array, dimension( N )
!         This contains the permutation which will reintegrate the
!         subproblem just solved back into sorted order, i.e.
!         D( IDXQ( I = 1,N ) ) will be in ascending order.
!
!  IWORK  (workspace) INTEGER array, dimension( 4 * N )
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension( 3*M**2 + 2*M )
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!  Based on LAPACK routine dlasd1
!
!  Written by S.-G. Li, NUDT, on Dec. 29th, 2012
!  =============================================
!
!     ..Parameters..
!
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     ..Local Scalars..
      LOGICAL            HFLAG
      INTEGER            COLTYP, I, IDX, IDXC, IDXP, ISIGMA, IU2, &
                         IVT2, IZ, K, LDQ, LDU2, LDVT2, M, N, N1, N2, &
                         IDIFL, IDIFR, IAPHAU, IAPHAV, ILEFT, J, IZZ, & 
                         IUF,IWK,ISTAZ
      DOUBLE PRECISION   ORGNRM, TEMP, TEMP1
!     ..
!     ..External Subroutines..
      EXTERNAL           DLAMRG,DLASCL,MDBDSVD2,DLASD3,XERBLA,DGER, &
                         MDBDSVD22,MDLASD4,MDBDSVD8
!     ..
!     ..Intrinsic Functions..
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
      ELSE IF( SB.LT.1 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.SB ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSVD1', -INFO )
         RETURN
      END IF
!
      N = NL + NR + SB     ! row dim
      M = N + SQRE         ! col dim
!
!     The following values are for bookkeeping purposes only.  They are
!     integer pointers which indicate the portion of the workspace
!     used by a particular array in MDBDSVD2 and MDBDSVD3.
!     Usually LDU2 << LDU, LDVT2 << LDVT. LDU and LDVT are fixed, equal to
!     the size of original matrix. LDU2 and LDVT2 are for local problems which 
!     may be much smaller.
!
      LDU2 = N
      LDVT2 = M
!     
!     WORK 
      IZ = 1               ! the block inserted row
      ISIGMA = IZ + M*SB   ! new singular values
      IU2 = ISIGMA + N     ! the left singular vectors of current problem
      IVT2 = IU2 + LDU2*N  ! the right singular vectors of current problem
      IWK = IVT2 + LDVT2*M ! double precision workspace
!
!     IWORK
      IDX = 1              ! permutation used to sort D into ascending order 
      IDXC = IDX + N       ! arrange the columns of U into four groups
      COLTYP = IDXC + N    ! type of each column of U2 and each row of VT2
      IDXP = COLTYP + N    ! permutation used to place deflated svals to the end of the array. 
!
!     Scale the X-shape matrix, does it only needs to scale once ? 
!     
      TEMP  = MAXVAL( ABS(ALPHA) )
      TEMP1 = MAXVAL( ABS(BETA)  )
      ORGNRM = MAX( TEMP, TEMP1 )
      D( NL+1:NL+SB ) = ZERO
      DO I = 1, N
         IF( ABS( D( I ) ).GT.ORGNRM ) THEN
            ORGNRM = ABS( D( I ) )
         END IF
      END DO
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      ALPHA = ALPHA / ORGNRM
      BETA  = BETA / ORGNRM
!
! ********************************************************************
!               Merging and deflating the first row                  *
! ********************************************************************
      CALL MDBDSVD2( NL,NR,SB,SQRE,K,D,WORK(IZ),ALPHA,BETA,U,LDU,VT, & 
               LDVT,WORK( ISIGMA ),WORK( IU2 ),LDU2,WORK( IVT2 ),LDVT2, & 
               IWORK(IDXP),IWORK(IDX),IWORK(IDXC),IDXQ,IWORK( COLTYP ), & 
               WORK(IWK),HFLAG,INFO )
      K = K-SB+1   ! K is now the size of secular equation
!
!     Use HFLAG to determine whether using HSS techniques or not. 
!     HFLAG = TRUE: use LAPACK routine; HFLAG = FALSE: use HSS matrices
      IF( HFLAG ) THEN
         LDQ = K
!        Define U2 and VT2, U and VT carefully
         CALL MDBDSVD3( NL,NR,SB,SQRE,K,D,WORK(IWK),LDQ,WORK(ISIGMA), &
              U, LDU,WORK( IU2 ),LDU2,VT,LDVT,WORK( IVT2 ), LDVT2,&
              IWORK( IDXC ),IWORK( COLTYP ),WORK( IZ ),M,INFO )
      ELSE
!
!        After calling dlasd2, the columns (rows) from SB+1 to N of U and VT store
!        the left and right singular vectors to be updated. There may be a better
!        order such as the place of IUF.
!
         IDIFL  = ISIGMA + N
         IDIFR  = IDIFL  + K
         IAPHAU = IDIFR  + K
         IAPHAV = IAPHAU + K
         IZZ    = IAPHAV + K
         IU2    = IZZ    + K
         IVT2   = IU2  + LDU2 *K
         IUF    = IVT2 + LDVT2*K
         ILEFT  = IUF  + LDU2      ! UF is the first column of U2
!
         ISTAZ = IZ+(LDVT2+1)*(SB-1)
         WORK( IZZ:IU2-1 ) = WORK( ISTAZ:(ISTAZ+K-1) )  ! ZZ = Z
!
!        Store the first K cols of U. 
         CALL DLACPY( 'A', LDU2, K, U(1,SB), LDU, WORK(IU2), LDU2 )
         CALL DLACPY( 'A', K, LDVT, VT(SB,1), LDVT, WORK(IVT2), K )   ! Store the rows of VT from SB to N
         CALL DLACPY( 'A', LDU2, 1, U(1,SB), LDU, WORK(IUF), LDU2 )   ! The SB-th col of U
!
         CALL MDBDSVD8( K,SB,D(SB),WORK(IZ),WORK(ISIGMA+SB-1),WORK(IDIFL), &
             WORK(IDIFR), WORK(IAPHAU), WORK(IAPHAV),WORK(IZZ), &
             WORK(IU2),LDU2,WORK(IVT2),LDVT2,WORK(ILEFT),INFO,nod )
         IF( INFO .NE. 0 ) THEN
            write(*,*) 'MDBDSVD8 failed ', -info
            return
         END IF
!
!        Update U and VT
         CALL DLACPY( 'A', K, LDVT2, WORK(IVT2), K, VT(SB,1), LDVT ) ! Update the first K rows of VT
!
!        Add the contribution of first column to the left singular vector matrix
!        IUF is a column vector and IAHPAU is a row vector. 
         CALL DGER( LDU2,K,1.0D0,WORK(IUF),1,WORK(IAPHAU),1,WORK(IU2), &
                    LDU2 )
         CALL DLACPY( 'A',LDU2,K,WORK(IU2),LDU2,U(1,SB),LDU )
!
      END IF ! ( HFLAG )
!
!     Prepare the IDXQ sorting permutation.
!
      N1 = K  ! The first K entries are in increasing order 
      N2 = N-K-SB+1  ! The last N2 entries are in increasing order
      CALL DLAMRG( N1, N2, D(SB), 1, -1, IDXQ(SB) )
!
! ****************************************************************
!                     Updating SVD problems                      *
! ****************************************************************
!
      DO I = SB-1, 1, -1  ! Number of added rows will be decreased by one
         IU2 = ISIGMA + N     ! the left singular vectors of current problem
         IVT2 = IU2 + LDU2*N  ! the right singular vectors of current problem
         IWK = IVT2 + LDVT2*M ! double precision workspace
         CALL MDBDSVD22( M,N,K,I,D,WORK(IZ),U,LDU,VT,LDVT,WORK(ISIGMA),&
              WORK(IU2),LDU2,WORK(IVT2),LDVT2,IWORK(IDXP),&
              IDXQ,HFLAG,INFO )
         K = K-I+1   ! K is now the size of secular equation
!
         IF( HFLAG ) THEN
!           Update the left and right singular vector matrices
            LDQ = K
            CALL MDBDSVD4( M,N,I,K,D,WORK(IWK),LDQ,WORK(ISIGMA), &
                 U, LDU,WORK( IU2 ),LDU2,VT,LDVT,WORK( IVT2 ), &
                 LDVT2,WORK(IZ),M,INFO )
         ELSE  ! use HSS technique
!
!           After calling mdbdsvd22, the columns (rows) from I to N of U and VT store
!           the left and right singular vectors to be updated. There may be a better
!           order such as the place of IUF.
            IDIFL  = ISIGMA + N
            IDIFR  = IDIFL  + K
            IAPHAU = IDIFR  + K
            IAPHAV = IAPHAU + K
            IZZ    = IAPHAV + K
            IU2    = IZZ    + K
            IVT2   = IU2  + LDU2 *K
            IUF    = IVT2 + LDVT2*K
            ILEFT  = IUF  + LDU2      ! UF is the first column of U2
!
            ISTAZ = IZ+(LDVT2+1)*(I-1)
            WORK( IZZ:IU2-1 ) = WORK( ISTAZ: (ISTAZ+K-1) )  ! ZZ = Z
!
!           Store the first K cols of U. 
            CALL DLACPY( 'A', LDU2, K, U(1,I), LDU, WORK(IU2), LDU2 )
            CALL DLACPY( 'A', K, LDVT, VT(I,1), LDVT, WORK(IVT2), K )   ! Store the rows of VT from I to N
            CALL DLACPY( 'A', LDU2, 1, U(1,I), LDU, WORK(IUF), LDU2 )   ! The I-th column of U
!
            CALL MDBDSVD8( K,I,D(I),WORK(IZ),WORK(ISIGMA+I-1),WORK(IDIFL), &
                 WORK(IDIFR), WORK(IAPHAU), WORK(IAPHAV),WORK(IZZ), &
                 WORK(IU2),LDU2,WORK(IVT2),LDVT2,WORK(ILEFT),INFO,nod )
!
!           Update U and VT
            CALL DLACPY( 'A', K, LDVT2, WORK(IVT2), K, VT(I,1), LDVT ) ! Update the first K rows of VT
!
!           Add the contribution of first column to the left singular vector matrix
!           IUF is a column vector and IAHPAU is a row vector. 
            CALL DGER( LDU2,K,1.0D0,WORK(IUF),1,WORK(IAPHAU),1,WORK(IU2), &
                 LDU2 )
            CALL DLACPY( 'A',LDU2,K,WORK(IU2),LDU2,U(1,I),LDU )
!
         END IF ! ( HFLAG )
!
         N1 = K         ! The first K entries are in increasing order 
         N2 = N-K-I+1   ! The last N2 entries are in increasing order
         CALL DLAMRG( N1, N2, D(I), 1, -1, IDXQ(I) )
      END DO
!
      IF( INFO.NE.0 ) THEN
         write(*,*) 'Info equals to ', info
         RETURN
      END IF
!
!     Unscale.
!
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
!
      RETURN
!
!     End of MDBDSVD1
!
    END SUBROUTINE MDBDSVD1

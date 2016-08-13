!
!  Definition:
!  ===========
!
!       SUBROUTINE MDSBEVD( COMPZ, UPLO, N, KD, AB, LDAB, D, Q, LDQ, WORK, 
!                          LWORK, IWORK, LIWORK, INFO )
! 
!       ..
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ, UPLO
!       INTEGER            INFO, KD, LDQ, LDAB, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), D( * ), WORK( * ), Q( LDQ, * )
!       ..
!
!  Purpose:
!  =============
!
!
! MDSBEDC computes all eigenvalues and, optionally, eigenvectors of a
! symmetric banded matrix using the divide and conquer method.
! Different than DSBEVD, it does not use DSYTRD or DSPTRD or DSBTRD to 
! reduce this band matrix to tridiagonal form.
!
! For band matrices with narrow bandwidth, it can be faster than DSBEVD. 
!  Different than DSTEVD, this routine does not check whether there are
!  splittings. 
!
! \endverbatim
!
!  Arguments:
!  ==========
!
! \param[in] COMPZ
! \verbatim
!          COMPZ is CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'I':  Compute eigenvectors of symmetric band matrix also.
!                  Now only this option is implemented. 
! \endverbatim
!
! \param[in] UPLO
! \verbatim
!          UPLO is CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!                  Now only this one is implemented.
!          = 'L':  Lower triangle of A is stored.
! \endverbatim
!
! \param[in] N
! \verbatim
!          N is INTEGER
!          The dimension of the symmetric banded matrix.  N >= 0.
! \endverbatim
!
! \param[in] KD
! \verbatim
!          KD is INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
! \endverbatim
!
! \param[in,out] AB
! \verbatim
!          AB is DOUBLE PRECISION array, dimension (LDAB, N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!          
! \endverbatim
!
! \param[in] LDAB
! \verbatim
!          LDAB is INTEGER
!          The leading dimension of the array AB.  LDAB >= KD + 1.
! \endverbatim
!
! \param[out] D
! \verbatim
!          D is DOUBLE PRECISION array, dimension (N)
!          On exit, if INFO = 0, the eigenvalues in ascending order.
! \endverbatim
!
! \param[in,out] Q
! \verbatim
!          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!          On exit, if INFO = 0, then if COMPZ = 'I', Z contains the
!          orthonormal eigenvectors of the original symmetric band matrix.
!          If  COMPZ = 'N', then Z is not referenced.
! \endverbatim
!
! \param[in] LDQ
! \verbatim
!          LDQ is INTEGER
!          The leading dimension of the array Q.  LDQ >= 1.
!          If eigenvectors are desired, then LDQ >= max(1,N).
! \endverbatim
!
! \param[out] WORK
! \verbatim
!          WORK is DOUBLE PRECISION array,
!                                         dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
! \endverbatim
!
! \param[in] LWORK
! \verbatim
!          LWORK is INTEGER
!          The dimension of the array WORK.
!          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
!          If COMPZ = 'V' and N > 1 then LWORK must be at least
!                         ( 1 + 3*N + 2*N*lg N + 4*N**2 ),
!                         where lg( N ) = smallest integer k such
!                         that 2**k >= N.
!          If COMPZ = 'I' and N > 1 then LWORK must be at least
!                         ( 1 + 4*N + N**2 ).
!          Note that for COMPZ = 'I' or 'V', then if N is less than or
!          equal to the minimum divide size, usually 25, then LWORK need
!          only be max(1,2*(N-1)).
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
! \endverbatim
!
! \param[out] IWORK
! \verbatim
!          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
! \endverbatim
!
! \param[in] LIWORK
! \verbatim
!          LIWORK is INTEGER
!          The dimension of the array IWORK.
!          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
!          If COMPZ = 'V' and N > 1 then LIWORK must be at least
!                         ( 6 + 6*N + 5*N*lg N ).
!          If COMPZ = 'I' and N > 1 then LIWORK must be at least
!                         ( 3 + 5*N ).
!          Note that for COMPZ = 'I' or 'V', then if N is less than or
!          equal to the minimum divide size, usually 25, then LIWORK
!          need only be 1.
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
! \endverbatim
!
! \param[out] INFO
! \verbatim
!          INFO is INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  The algorithm failed to compute an eigenvalue while
!                working on the submatrix lying in rows and columns
!                INFO/(N+1) through mod(INFO,N+1).
! \endverbatim
!
!  Authors:
!  ========
!
! \author Nati. Univ. of Def. Tech.
! \author Shengguo Li
!
! \date October 2014
!
! \ingroup auxOTHERcomputational
!
! \par Contributors:
!  ==================
!
! Shengguo Li, College of Computer Science, NUDT, China, 2014
!
!  =====================================================================
      SUBROUTINE MDSBEVD( COMPZ,UPLO,N,KD,AB,LDAB,D,Q,LDQ,WORK,LWORK, &
                        IWORK,LIWORK,INFO )
!
        USE BasicMM
!
        IMPLICIT NONE
!  -- HSSPACK computational routine (version 1.0.0) --
!  -- HSSPACK is a software package provided by NUDT
!  -- October 2014
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ, UPLO
      INTEGER            INFO, KD, LDQ, LDAB, LIWORK, LWORK, N, SUBPBS, TLVLS, MSD2, &
                         SPM1, IC, ICR, CURLVL, SPM2, SUBMAT, MATSIZ, CURPRB, INODE 
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AB( LDAB,* ), D( * ), WORK( * ), Q( LDQ,* )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, ICOMPZ, II, J, K, LGN, LIWMIN, &
                         LWMIN,M,SMLSIZ,IU,IVT,IWW,IWT,LIWT,&
                         IINFO,LWORK1,LIWORK1,IA,NL,INDXQ
      DOUBLE PRECISION   EPS, P
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLAED0, DLASCL, DLASET, DLASRT, &
                         DSTEQR, DSTERF, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX, MOD, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDQ.LT.1 ) .OR. &
               ( ICOMPZ.GT.0 .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
         SMLSIZ = ILAENV( 9, 'DSTEDC', ' ', 0, 0, 0, 0 )
!        Modify the definition of SMLSIZ when KD is large ??
!
         IF( N.LE.1 .OR. ICOMPZ.EQ.0 ) THEN
            LIWMIN = 1
            LWMIN = 1
         ELSE IF( N.LE.SMLSIZ ) THEN
            LIWMIN = 1
            LWMIN = 2*( N - 1 )
         ELSE
            LGN = INT( LOG( DBLE( N ) )/LOG( TWO ) )
            IF( 2**LGN.LT.N ) &
               LGN = LGN + 1
            IF( 2**LGN.LT.N ) &
               LGN = LGN + 1
            IF( ICOMPZ.EQ.1 ) THEN
               LWMIN = 1 + 3*N + 2*N*LGN + 4*N**2
               LIWMIN = 6 + 6*N + 5*N*LGN
            ELSE IF( ICOMPZ.EQ.2 ) THEN
               LWMIN = 1 + 4*N + N**2
               LIWMIN = 3 + 5*N
            END IF
         END IF
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
!
         IF( LWORK.LT.LWMIN .AND. .NOT. LQUERY ) THEN
            INFO = -8
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT. LQUERY ) THEN
            INFO = -10
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDSBEVD', -INFO )
         RETURN
      ELSE IF (LQUERY) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      IF( N.EQ.1 ) THEN 
         IF( ICOMPZ.NE.0 ) &
            Q( 1, 1 ) = ONE
         RETURN
      END IF
!
      IF( ICOMPZ.EQ.2 ) THEN
         CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      END IF
!
!        Scale.
!
!      ORGNRM = DLANST( 'M', N, D, E )
!      IF( ORGNRM.EQ.ZERO ) &
!           GO TO 50
!
      EPS = DLAMCH( 'Epsilon' )
!
!     Set up the computation tree, and determine the size and placement 
!     of the submatrices, and save in the leading elements of IWORK.
      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
         END DO
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
!     IWORK stores the ending position of each bottom subproblem.
      DO J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
      END DO
!
!     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
!     using rank-KD modifications (cuts). 
!
!     The workspace index for B_i and its SVD
!     Does the initialization of WORK causes problem ??
      IWT = 1
      LIWT = 5*KD
!
!     The following part can be modified to run in parallel. But, it currently
!     runs sequentially. The SVD of B_i are stored sequentially in WORK.  
!
!     Here I do not know whether the following is the right way, it requires much more
!     workspace than the sequential version. 
      SPM1 = SUBPBS - 1
      DO I = 1, SPM1
         IC = IWORK( I )
         ICR = IC + 1

         IU  = IWT         ! first stores the off-diagonal block, and finally
!                            the left singular vectors of the off-diagonal block
         IVT = IU + KD*KD  ! the right singular vectors of the off-diagonal block
         IWW = IVT+ KD*KD  ! the singular values of the off-diagonal block
         IWT = IWW+ KD     ! the left workspace

!        Compute the SVD of off-diagonal block B_i 
         CALL DSB2FLC( 'L',KD,KD,KD,AB(1,ICR),LDAB,WORK(IU),KD,IINFO )
         CALL DGESVD( 'O','A',KD,KD,WORK(IU),KD,WORK(IWW),WORK(IU),KD,WORK(IVT),KD, &
                     WORK(IWT),LIWT,IINFO )

!        Compute A_i = A_i - U_{i-1} \Sigma_{i-1} U_{i-1}^T - V_i \Sigma_i V_i^T
!        The correctness of DSymSplit is checked. 
         CALL DSymSplit( AB,LDAB,KD,WORK(IU),KD,WORK(IVT),KD,WORK(IWW),IC,WORK(IWT),INFO )
!
      END DO

!     For the nodes at the bottom level, solve these subproblems by DSYEVD.
      IC = 1
      IA = IWT
      INDXQ = 4*N+3    ! Do not know why ??
      DO I = 1, SUBPBS
        NL = IWORK(I) - IC +1
        IWT = IA + NL*NL
        CALL DSB2FLC( 'U',NL,NL,KD,AB(1,IC),LDAB,WORK(IA),NL,INFO ) 
        LWORK1 = 2*NL*NL+6*NL+1
        LIWORK1 = 5*NL + 3
        CALL DSYEVD( 'V','U',NL,WORK(IA),NL,D(IC),WORK(IWT),LWORK1,IWORK(SUBPBS+1),LIWORK1,INFO ) 
        IF( INFO.NE.0 ) THEN
           RETURN
        END IF
        CALL DLACPY( 'A',NL,NL,WORK(IA),NL,Q(IC,IC),LDQ )

        K = 1
        DO J = IC, IWORK(I)
           IWORK( INDXQ+J ) = K
           K = K + 1
        END DO  
        
        IC = IWORK(I)+1   ! The starting position of each subproblem
      END DO  ! (I, bottom subproblems)
!
!
!     Successively merge eigensystems of adjacent submatrices into eigensystem 
!     for the corresponding larger matrix.
!
!     while ( SUBPBS > 1 )
!
      CURLVL = 1
   80 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0

!           Determine the starting position of the singular vectors of Bi
               INODE = 2**(CURLVL-1)-1
               IU = (2*KD*KD+KD)*INODE+1
               IVT = IU + KD*KD
               IWW = IVT + KD*KD
            ELSE
               SUBMAT = IWORK( I ) + 1  ! the starting position of current merged prb
               MATSIZ = IWORK( I+2 ) - IWORK( I )  ! the size of the current merged prb
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1

!           Determine the starting position of the singular vectors of Bi
               INODE = INODE+ 2**CURLVL
               IU = (2*KD*KD+KD)*INODE+1
               IVT = IU + KD*KD
               IWW = IVT + KD*KD
            END IF
!
!     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
!     into an eigensystem of size MATSIZ. MDSBED1 is used for solving this problem.
!
            IF( ICOMPZ.EQ.2 ) THEN
!
               CALL MDSBED1( MATSIZ, D(SUBMAT), Q( SUBMAT,SUBMAT ), LDQ, &
                           IWORK( INDXQ+SUBMAT ),MSD2,WORK(IU),KD,WORK(IVT),&
                           WORK(IWW),WORK(IA),IWORK( SUBPBS+1 ),INFO,INODE )
            ELSE
             WRITE(*,*) 'If eigenvetors not require, use DSBEVD instead'
               GO TO 50
            END IF
            IF( INFO.NE.0 )  &
               GO TO 50
            IWORK( I/2+1 ) = IWORK( I+2 )
         END DO
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      END IF
!
!     end while

!     If the problem split any number of times, then the eigenvalues
!     will not be properly ordered.  Here we permute the eigenvalues
!     (and the associated eigenvectors) into ascending order.
!
         IF( M.NE.N ) THEN
            IF( ICOMPZ.EQ.0 ) THEN
!
!              Use Quick Sort
!
               CALL DLASRT( 'I', N, D, INFO )
!
            ELSE
!
!              Use Selection Sort to minimize swaps of eigenvectors
!
               DO II = 2, N
                  I = II - 1
                  K = I
                  P = D( I )
                  DO J = II, N
                     IF( D( J ).LT.P ) THEN
                        K = J
                        P = D( J )
                     END IF
                  END DO
                  IF( K.NE.I ) THEN
                     D( K ) = D( I )
                     D( I ) = P
                     CALL DSWAP( N, Q( 1, I ), 1, Q( 1, K ), 1 )
                  END IF
               END DO
            END IF  ! (ICOMPZ)
         END IF ! (M .ne. N)
!
   50 CONTINUE
      WORK( 1 )  = LWMIN
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of MDSBEVD
!
      END SUBROUTINE MDSBEVD

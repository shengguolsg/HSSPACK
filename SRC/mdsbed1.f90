!  
     SUBROUTINE MDSBED1( N,D,Q,LDQ,INDXQ,CUTPNT,U,KD,VT,W,WORK,IWORK, &
                         INFO, INODE )
       USE BasicMM
!
       IMPLICIT NONE
!       ..
!       .. Scalar Arguments ..
        INTEGER            CUTPNT, INFO, LDQ, N, KD, INODE
!       ..
!       .. Array Arguments ..
        INTEGER            INDXQ( * ), IWORK( * )
        DOUBLE PRECISION   D( * ), Q( LDQ,* ), U(KD,*), VT(KD,*), W(*), WORK( * )
!       ..
!
!  Purpose
!  ========
!
! MDSBED1 computes the updated eigensystem of a diagonal matrix after 
! modification by a rank-KD symmetric matrix.  This routine is used only for 
! the eigenproblem which requires all eigenvalues and eigenvectors of a symmetric 
! band matrix.
!
!    T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
!
!    where Z = Q**T*u, u is a matrix with dimension N-by-KD, at the middle
!    of which contains the singular vectors of a lower triangular matrix
!    and zeros elsewhere. u = [0 U^T VT 0 ]^T, where B_i = U*Sigma*VT. 
!
!    The eigenvectors of the original matrix are stored in Q, and the
!    eigenvalues are in D.  The algorithm consists of a loop, each
!    iteration has three stages:
!
!       The first stage consists of deflating the size of the problem
!       when there are multiple eigenvalues or if there is a zero in
!       the Z vector.  For each such occurence the dimension of the
!       secular equation problem is reduced by one.  This stage is
!       performed by the routine DLAED2.
!
!       The second stage consists of calculating the updated eigenvalues. 
!       This is done by finding the roots of the secular
!       equation via the routine DLAED4 (as called by DLAED3).
!       This routine also calculates the eigenvectors of the current
!       problem.
!
!       The final stage consists of computing the updated eigenvectors
!       directly using the updated eigenvalues.  The eigenvectors for
!       the current problem are multiplied with the eigenvectors from
!       the overall problem.
!
! Note that this routine can also be compiled by using F77 after small modification, 
! it does not use any new features of F90. 
! 
!  Arguments
!  ===========
!
! \param[in] N
!         N is INTEGER
!         The dimension of the current considered symmetric band matrix.  N >= 0.
!
! \param[in,out] D
!         D is DOUBLE PRECISION array, dimension (N)
!         On entry, the eigenvalues of the rank-KD-perturbed matrix.
!         On exit, the eigenvalues of the repaired matrix.
!
! \param[in,out] Q
!         Q is DOUBLE PRECISION array, dimension (LDQ,N)
!         On entry, the eigenvectors of the rank-KD-perturbed matrix.
!         On exit, the eigenvectors of the repaired tridiagonal matrix.
!
! \param[in] LDQ
!         LDQ is INTEGER
!         The leading dimension of the array Q.  LDQ >= max(1,N).
!
! \param[in,out] INDXQ
!         INDXQ is INTEGER array, dimension (N)
!         On entry, the permutation which separately sorts the two
!         subproblems in D into ascending order.
!         On exit, the permutation which will reintegrate the
!         subproblems back into sorted order,
!         i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.
!
! \param[in] CUTPNT
!         CUTPNT is INTEGER
!         The location of the last eigenvalue in the leading sub-matrix.
!         min(1,N) <= CUTPNT <= N/2.
!
! \param[in] U
!          U is DOUBLE PRECISION array, dimension (KD,KD) 
!          The left singular vector matrix of the lower triangular
!          matrix B_i
!
! \param[in] KD
!          KD is an integer, the semi-bandwith of this banded matrix
!
! \param[in] VT
!          VT is DOUBLE PRECISION array, dimension (KD,KD) 
!          The transpose of right singular vectors of the lower triangular
!          matrix B_i
!
! \param[in] W
!          W is DOUBLE PRECISION array, dimension (KD) 
!          The singular values of the lower triangular matrix B_i
!
! \param[out] WORK
!          WORK is DOUBLE PRECISION array, dimension ( 4*N + N**2 )
!
! \param[out] IWORK
!          IWORK is INTEGER array, dimension ( 4*N )
!
! \param[out] INFO
!          INFO is INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an eigenvalue did not converge
!
!  Authors
!  ========
!
! \author Nat. Univ. of Def. Tech. 
! \date October 2014
!
! Contributors
! =============
!
! Shengguo Li, NUDT, China, 2014
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            HFLAG
      INTEGER            COLTYP, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS, &
                         IW, IZ, K, N1, N2, ZPP1, IQPP, IZ2, SB, IQ1, II
      INTEGER, PARAMETER :: BSMLZ = 1500
      DOUBLE PRECISION   :: RHO
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D0, ONE=1.0D0
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAMRG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDQ.LT.MAX( 1,N ) ) THEN
         INFO = -4
      ELSE IF( MIN( 1, N/2 ).GT.CUTPNT .OR. ( N/2 ).LT.CUTPNT ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDSBED1', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN

      HFLAG = .TRUE.   ! If it is ture, use LAPACK routines, else use HSS techniques.
!
!     The following values are integer pointers which indicate
!     the portion of the workspace used by a particular array in DLAED2 and DLAED3.
!     Matrices Z and Z2 are used to store the appended vectors, and they are used in 
!     a ping-pong form in MDSBED22. Z and Z2 must return. 
!
      IZ  = 1                   ! The start position of Z in WORK
      IZ2 = IZ  + N*KD          ! The start position of Z2 in WORK
      IQ1 = IZ2 + N*KD          ! The start position for the last KD rows of Q1
      IQ2 = IQ1 + CUTPNT*KD     ! The start position for the first KD rows of Q2
      IW  = IQ2 + (N-CUTPNT)*KD ! The remaining workspace
!
      INDX   = 1
      INDXC  = INDX + N
      COLTYP = INDXC + N
      INDXP  = COLTYP + N
!
!     Form the Z-matrix which consists of the last KD rows of Q_1 and the
!     first KD rows of Q_2.
      IQPP = CUTPNT-KD+1
      CALL DLACPY( 'A', KD, CUTPNT, Q(IQPP,1), LDQ, WORK(IQ1), KD )
      ZPP1 = CUTPNT + 1
      CALL DLACPY( 'A', KD, N-CUTPNT, Q(ZPP1,ZPP1), LDQ, WORK(IQ2), KD )
!      
!     The first half of Z-matrix
      CALL DGEMM( 'T','N',CUTPNT,KD,KD,ONE,WORK(IQ1),KD,U,KD,ZERO,WORK(IZ),N )
!     The last half of Z-matrix
      CALL DGEMM( 'T','T',N-CUTPNT,KD,KD,ONE,WORK(IQ2),KD,VT,KD,ZERO,WORK(IZ+CUTPNT),N )
!
!     The main loop for KD updates
!
      IDLMDA = IQ1         ! The updated singular values 
      IW = IDLMDA + N      ! A copy of the first column of Z, used for the broken arrow matrix
      IQ2 = IW + N         ! The size of Q2 is N*K ??
!
!     Deflate eigenvalues for the first column
      RHO = W( 1 )
      CALL MDSBED2( K,N,CUTPNT,KD,D,Q,LDQ,INDXQ,RHO,WORK(IZ),WORK(IZ2), &
           WORK( IDLMDA ), WORK( IW ), WORK( IQ2 ), &
           IWORK( INDX ), IWORK( INDXC ), IWORK( INDXP ), &
           IWORK( COLTYP ), HFLAG, BSMLZ, INFO ) 
!
!     After MDSBED2, Z2 stores the permuted appended vectors which is transformed from 
!     Z, and the column dimension is reduced to KD-1 from KD. 
!
      IF( INFO.NE.0 ) GO TO 20
!
!     Solve Secular Equation.
      IF( K.NE.0 ) THEN
         SB = KD - 1 
!
         IF( HFLAG ) THEN
            IS = ( IWORK( COLTYP )+IWORK( COLTYP+1 ) )*CUTPNT + &
                 ( IWORK( COLTYP+1 )+IWORK( COLTYP+2 ) )*( N-CUTPNT ) + IQ2
            ! Use Lapack technique
            CALL MDSBED3( K, N, CUTPNT,SB, D, Q, LDQ, RHO, WORK( IDLMDA ), &
                 WORK( IQ2 ), IWORK( INDXC ), IWORK( COLTYP ), &
                 WORK( IW ), WORK( IS ), WORK(IZ2),WORK(IZ), INFO ) 
!           After MDSBED3 the appended vectors are stored in Z again. 
            IF( INFO.NE.0 ) THEN 
               write(*,*) 'MDSBED3 failed info=',info
               GO TO 20
            END IF
         ELSE
            ! use HSS techniques 
            IS = IQ2 + N*K
            CALL MDSBED5( K,N,SB,D,Q,LDQ,RHO,WORK(IDLMDA),WORK(IQ2),WORK(IW),&
                 WORK(IS),WORK(IZ2),WORK(IZ),INFO )
!           After MDSBED5, the entries in Z2 are transformed to Z.
            IF( INFO .NE. 0 )  THEN
               write(*,*) 'MDSBED5 is failed', INFO
               GO TO 20
            END IF
!
         END IF  ! ( HFLAG )
!
!       Prepare the INDXQ sorting permutation.
         N1 = K
         N2 = N - K
         CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         DO I = 1, N
            INDXQ( I ) = I
         END DO
      END IF

! ****************************************************************
!                     Updating SVD problems                      *
! ****************************************************************
!
      DO II = 2, KD
         RHO = W(II)
         SB = KD -II+1
!
         CALL MDSBED22( K,N,N1,SB,D,Q,LDQ,RHO,WORK(IZ),WORK(IZ2),WORK(IDLMDA), &
              WORK(IW),WORK(IQ2),INDXQ,IWORK(INDXP),INFO )
         IF( INFO .NE. 0 )  THEN
            write(*,*) 'MDSBED22 is failed', INFO, K
            GO TO 20
         END IF
!        After calling MDSBED22, the entries in Z are transformed to Z2
!         
!        Solve Secular Equation.
         IF( K.NE.0 ) THEN
            IS = IQ2 + N*K
!
            IF( K.LT. BSMLZ ) THEN
               CALL MDSBED4( K, N, SB, D, Q, LDQ, RHO, WORK( IDLMDA ), &
                    WORK(IQ2),INDXQ, WORK(IW), WORK(IS), WORK(IZ2), WORK(IZ),INFO )
!              After MDSBED4, the entries in Z2 are transformed to Z.
               IF( INFO .NE. 0 )  THEN
                  write(*,*) 'MDSBED4 is failed', INFO
                  GO TO 20
               END IF
!
            ELSE
               CALL MDSBED5( K, N, SB, D, Q, LDQ, RHO, WORK( IDLMDA ), &
                    WORK(IQ2),WORK(IW), WORK(IS), WORK(IZ2), WORK(IZ),INFO )
!              After MDSBED5, the entries in Z2 are transformed to Z.
               IF( INFO .NE. 0 )  THEN
                  write(*,*) 'MDSBED5 is failed', INFO
                  GO TO 20
               END IF
!
            END IF
!
!       Prepare the INDXQ sorting permutation.
            N1 = K
            N2 = N - K
            CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
         ELSE
            DO I = 1, N
               INDXQ( I ) = I
            END DO
         END IF
         
      END DO ! (II)
!
20    CONTINUE
!
      RETURN
!
!     End of MDSBED1
!
     END SUBROUTINE MDSBED1

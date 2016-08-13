      SUBROUTINE MDBDSVD( N,SB,AB,LDAB,A,S,U,LDU,VT,LDVT,SMLSIZ,WORK,&
          LWORK,IWORK,INFO )
!
        USE TestRoutine
        use basicmm
        IMPLICIT NONE
!
!  -- HSSPACK computational routine (version 0.0.1) --
!  -- HSSPACK is a software package provided by Nati. Univ. of Def. Tech., China.
!     January 2013
!
!  .. Scalar Arguments ..
      INTEGER, INTENT(IN)    :: N,SB,LDAB,LDU,LDVT
      INTEGER, INTENT(INOUT) :: INFO,LWORK,SMLSIZ
!  ..
!  .. Array Arguments ..
      DOUBLE PRECISION, INTENT(INOUT) :: AB(LDAB,*),A(*),S(*),U(LDU,*), &
               VT(LDVT,*),WORK(*)
      INTEGER, INTENT(INOUT) :: IWORK(*)
!  ..
!  Purpose
!  =======
! 
!  MDBDSVD computes the SVD of an N-by-N upper banded triangular matrix AB. 
!  AB is stored in compact form, only its nonzero diagonals are stored, and 
!  its semibandwidth is SB. It computes both the singular values and singular 
!  vectors,  and it uses HSS techniques to compute the singular vectors when 
!  the secular equations are large. 
!
!  See S.-G. Li's paper about structured banded divide-and-conquer 
!  (SBDC) algorithm.
!  
! Argumente
! ==========
!
! N (in) INTEGER
!        The order the matrix AB.  N >= 0.
!
! SB (in) INTEGER
!         The number of superdiagonals of the matrix A if UPLO = 'U',
!         or the number of subdiagonals if UPLO = 'L'.  SB >= 0. 
!         Currently we only consider the upper banded case.
! 
! AB (inout) DOUBLE PRECISION array, DIMENSION( LDAB,N )
!          On entry, an upper banded SB-by-N matrix AB. AB is stored in
!          sparse form, only the nonzero diagonals are stored. 
!
!          On exit, AB is unchanged, and its entries would be copied to 
!          a square matrix. 
!
! LDAB (in) INTEGER
!          The leading dimension of the array A.  LDAB >= SB+1
!
! A (out)  DOUBLE PRECISION array, dimension( * )
!          It is used as workspace to store AB in full form. It can be
!          the original dense matrix which is transformed to banded form 
!          by orthogonal transformations. 
!
! S (out)  DOUBLE PRECISION array, dimension( N )
!          The singular values of AB, sorted in descending order S(i) >= S(i+1).
!
! U (out)  DOUBLE PRECISION array, dimension( LDU,N )
!          U contains the N-by-N orthogonal matrix U;
!          ( the left singular vectors, stored columnwise ); 
!          Another way of writing the codes is to copy the entries of AB to the main 
!          diagonal and super diagonals of U.
!
! LDU (in) INTEGER
!          The leading dimension of the array U.  LDU >= N
!
! VT (out) DOUBLE PRECISION array, dimension( LDVT,N )
!          VT contains the N-by-N orthogonal matrix V**T;
!          ( the right singular vectors, stored rowwise );
!
! LDVT (in) INTEGER
!           The leading dimension of the array VT.  LDVT >= N.
!
! SMLSIZ (in) INTEGER
!          The size of smallest subproblem, it usually equals to 25.
! 
! WORK (inout) DOUBLE PRECISION array, dimension( MAX(1,LWORK) )
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!
! LWORK (in) INTEGER, The dimension of the array WORK.
!          The required size of WORK is not clear yet. 
!
! IWORK  (out) INTEGER work array.
!         with dimension 8*N 
!
! INFO (out) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  the algorithm did not converge.
!
! ==========
! Written by S.-G. Li, on Jun 4th, 2013
! 
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK, &
                         J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR, &
                         NL, NLF, NLVL, NR, NRF, SQREI, ICR, NLD, SQRE, &
                         NRD, WU, WVT, WWK, II, NRFM1, NLFM1, KK,RS,TI, &
                         TN,TM
      DOUBLE PRECISION   BlkALPHA(SB,SB), BlkBETA(SB,SB), P
!     BlkALPHA: upper triangular in main diagonal; 
!     BlkBETA:  lower triangular in off diagonal part;
!     ..
!     .. External Subroutines ..
      EXTERNAL           MDBDSVD1, DGESVD, MDBDSVDT, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      SQRE = 0
!
      IF( (N.LT.0) .OR. (N.LT.SB) ) THEN
         INFO = -2
      ELSE IF( SB.LT.1 ) THEN
         INFO = -3
      END IF
!
      IF( LDU.LT.N ) THEN
         INFO = -9
      ELSE IF( LDVT.LT.N ) THEN
         INFO = -10
      ELSE IF( SMLSIZ.LT.3 ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MDBDSVD0', -INFO )
         RETURN
      END IF
!
!     If the input matrix is too small, call DGESVD to finsd the SVD.
      IF( N.LE.SMLSIZ ) THEN
         LWORK = 5*N
         CALL DSB2FLC( 'U',N,N,SB,AB,LDAB,A,N,INFO )   ! copy sparse to full
         CALL DGESVD( 'A','A',N,M,A,N,S,U,LDU,VT,LDVT,WORK,LWORK,INFO )
         RETURN
      END IF
!
!     Set up the computation tree.
!     A tree needs 3N integer workspace, center, row dims of left and right
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
      IDXQ  = NDIMR + N
      IWK   = IDXQ  + N
      CALL MDBDSVDT( N, SB, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), &
                   IWORK( NDIMR ), SMLSIZ )
!
!     For the nodes at the bottom level of the tree, solve
!     these subproblems by DGESVD.
!
      NDB1 = ( ND+1 ) / 2  ! the left most leaf node
      NCC  = 0
      WU = 1
      DO I = NDB1, ND      ! leaf nodes, each of which have two subproblems
!
!     IC : center row of each node
!     NL : number of rows of the left subproblem
!     NR : number of rows of the right subproblem
!     NLF: starting row of the left subproblem
!     NRF: starting row of the right subproblem
!     
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )   ! row dimension of the left subproblem
         NR = IWORK( NDIMR+I1 )   ! row dimension of the right subproblem
         NLF = IC - NL
         NRF = IC + SB
         NLD = NL + SB
         WVT = WU + NL*NL
         WWK = WVT + NLD*NLD
         LWORK = 5*NLD      ! choose a bigger one for efficiency 
         CALL DSB2FLC( 'U',NL,NLD,SB,AB(1,NLF),LDAB,A,NLD,INFO )
         CALL DGESVD( 'A','A',NL,NLD,A,NLD,S(NLF),WORK(WU),NL,&
                    WORK(WVT),NLD,WORK(WWK),LWORK,INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
!        Use Selection Sort to minimize swaps of singular vectors
         NLFM1 = NLF -1
         DO II = 2, NL
            TI = II - 1
            KK = TI
            P = S( NLFM1+TI )
            DO J = II, NL
               IF( S( J+NLFM1 ).LT.P ) THEN
                  KK = J
                  P = S( J+NLFM1 )
               END IF
            END DO
            IF( KK.NE.TI ) THEN
               S( KK+NLFM1 ) = S( TI+NLFM1 )
               S( TI+NLFM1 ) = P
               CALL DSWAP( NL, WORK( WU+(TI-1)*NL ), 1, WORK( WU+(KK-1)*NL ), 1 )
               CALL DSWAP( NLD, WORK( WVT+TI-1 ), NLD, WORK( WVT+KK-1 ), NLD )
            END IF
         END DO
         CALL DLACPY( 'A',NL,NL,WORK(WU),NL,U(NLF,NLF),LDU )
         CALL DLACPY( 'A',NLD,NLD,WORK(WVT),NLD,VT(NLF,NLF),LDVT )
!
!!$         write(*,*) 'Left node at the bottom level'
!!$         Call TestSVD(NL,NLD,SB,AB(1,NLF),LDAB,A,U(NLF,NLF),LDU,S(NLF), &
!!$              VT(NLF,NLF),LDVT,WORK, info)
!!$         call testorth(VT(NLF,NLF),LDVT,NLD,NLD)
!
         ITEMP = IDXQ + NLF - 2  ! Here minus 2 is because ITEMP will be added by one later.
         DO J = 1, NL
            IWORK( ITEMP+J ) = J ! IDXQ for the left subproblem
         END DO
         IF( I.EQ.ND ) THEN      ! the most right subproblem
            SQREI = SQRE
         ELSE
            SQREI = SB
         END IF
         NRD = NR + SQREI
         WVT = WU + NR*NR
         WWK = WVT + NRD*NRD
         LWORK = 5*NRD            ! choose a bigger one for efficiency
         CALL DSB2FLC( 'U',NR,NRD,SB,AB(1,NRF),LDAB,A,NRD,INFO )
         CALL DGESVD( 'A','A',NR,NRD,A,NRD,S(NRF),WORK(WU),NR, &
                    WORK(WVT),NRD,WORK(WWK),LWORK,INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
!        Use Selection Sort to minimize swaps of singular vectors
         NRFM1 = NRF -1
         DO II = 2, NR
            TI = II - 1
            KK = TI
            P = S( NRFM1+TI )
            DO J = II, NR
               IF( S( J+NRFM1 ).LT.P ) THEN
                  KK = J
                  P = S( J+NRFM1 )
               END IF
            END DO
            IF( KK.NE.TI ) THEN
               S( KK+NRFM1 ) = S( TI+NRFM1 )
               S( TI+NRFM1 ) = P
               CALL DSWAP( NR,WORK( WU+(TI-1)*NR ),1,WORK( WU+(KK-1)*NR ),1 )
               CALL DSWAP( NRD,WORK( WVT+TI-1 ),NRD, WORK( WVT+KK-1 ),NRD )
            END IF
         END DO
         CALL DLACPY( 'A',NR,NR,WORK(WU),NR,U(NRF,NRF),LDU )
         CALL DLACPY( 'A',NRD,NRD,WORK(WVT),NRD,VT(NRF,NRF),LDVT )
!
!!$         write(*,*) 'Right node at the bottom level', NRD, NR
!!$         Call TestSVD( NR,NRD,SB,AB(1,NRF),LDAB,A,U(NRF,NRF),LDU,S(NRF), &
!!$              VT(NRF,NRF),LDVT,WORK,info )
!!$         call testorth( VT(NRF,NRF),LDVT,NRD,NRD )
!
         ITEMP = IDXQ + IC + SB-1   ! starting row of the right subproblem
         DO J = 1, NR
            IWORK( ITEMP+J-1 ) = J  ! IDXQ for the right subproblem
         END DO 
      END DO  ! ( I )
!
!     Now conquer each subproblem bottom-up.
!
!      write(*,*) 'There are totally ', nlvl, ' levels '
      DO LVL = NLVL, 1, -1
!
!        Find the first node LF and last node LL on the
!        current level LVL.
!
!!$         write(*,*) '******************************************' 
!!$         write(*,*) 'Now we are at level ', lvl
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
            NL = IWORK( NDIML+IM1 )  ! row dimension of the left
            NR = IWORK( NDIMR+IM1 )  ! row dimension of the right
            NLF = IC - NL            ! starting row position of left subproblem
            IF( ( SQRE.EQ.0 ) .AND. ( I.EQ.LL ) ) THEN   ! for right problem
               SQREI = SQRE
            ELSE
               SQREI = SB
            END IF
            IDXQC = IDXQ + NLF - 1   ! IDXQC would be the first position of IDXQ for current merged subproblem
            ICR  = IC + SB   ! starting row of right subproblem
            RS   = LDAB-SB
!            write(*,*) 'We are at node ', I, 'of level ', lvl, SB
            CALL DSB2FLC( 'U',SB,SB,SB,AB(1,IC ),LDAB,BlkALPHA,SB,INFO )   ! copy the inserted main diagonal block
            CALL DSB2FLC( 'L',SB,SB,SB,AB(RS,ICR),LDAB,BlkBETA,SB,INFO )   ! copy the inserted super diagonal block
            CALL MDBDSVD1( NL, NR, SB, SQREI, S(NLF), BlkALPHA, BlkBETA, &
                          U( NLF,NLF ), LDU, VT( NLF,NLF ), LDVT, IWORK(IDXQC),&
                          IWORK(IWK),WORK,INFO,NLF,A,AB,LDAB, I )
!
            TN = NL + NR + SB
            TM = TN + SQREI
!!$            CALL TestSVD( TN,TM,SB,AB(1,NLF),LDAB,A,U(NLF,NLF),LDU,S(NLF), &
!!$              VT(NLF,NLF),LDVT,WORK,info )
!
            IF( INFO.NE.0 ) THEN
               RETURN
            END IF
         END DO ! ( I )
      END DO ! ( LVL )
!
!     Use Selection Sort to minimize swaps of singular vectors
!
      DO 60 II = 2, N
         I = II - 1
         KK = I
         P = S( I )
         DO 50 J = II, N
            IF( S( J ).GT.P ) THEN
               KK = J
               P = S( J )
            END IF
   50    CONTINUE
         IF( KK.NE.I ) THEN
            S( KK ) = S( I )
            S( I ) = P
            CALL DSWAP( N, U( 1, I ), 1, U( 1, KK ), 1 )
            CALL DSWAP( N, VT( I, 1 ), LDVT, VT( KK, 1 ), LDVT )
         END IF
   60 CONTINUE
!      
      RETURN
!
!     End of MDBDSVD
!
      END SUBROUTINE MDBDSVD

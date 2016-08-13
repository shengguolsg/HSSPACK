*  Definition:
*  ===========
*
*       SUBROUTINE TORTHSVD( M, N, A, LDA, Dt, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            M, N, LDA
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), Dt( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> TORTHSVD tests the orthogonality of A by computing its SVD. 
*> A = U * S * VT. Compute Max( Abs(S(1) -ONE ), Abs( S(N)-ONE ) ).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*> \endverbatim
*>
*> \param[out] Dt
*> \verbatim
*>          Dt is DOUBLE PRECISION array, dimension ( MIN( M, N ) )
*>          The singular values of A, used for testing orthognoality.
*> \endverbatim
*>
*> \param[in] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension ( LWORK )
*>          Used as workspace.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          INTEGER, the size of WORK.
*>          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)). The larger the better.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if DBDSQR did not converge, INFO specifies how many
*>                superdiagonals of an intermediate bidiagonal form B
*>                did not converge to zero. See the description of WORK
*>                above for details.
*> \endverbatim
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Written by
*>   S. G. Li, Dept. of Math., Nat. Univ. of Defense Tech., Changsha China
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE TORTHSVD( M, N, A, LDA, Dt, WORK, LWORK, INFO )
*
*  -- HSSPACK auxiliary routine (version 0.0.1 ) --
*  -- HSSPACK is a software package provided by National Univ. of Defense Tech.    --
*  -- Novermber 2012                                                      --
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA, LWORK, INFO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), Dt( * ), WORK( * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            MN
      DOUBLE PRECISION   ERR, ONE
      PARAMETER          (ONE = 1.0E+0)
*     ..
*     .. Executable Statements ..
*
      LWORK = 5* MAX(M, N)
      CALL DGESVD('N','N',M,N,A,LDA,Dt,A,M,A,N,WORK,LWORK,INFO)
      IF( INFO .NE. 0 ) THEN
         WRITE(*,*) 'Not complete, INFO is ', info
         RETURN
      END IF
      
      MN = MIN( M, N )
      ERR = MAX( ABS(Dt(1) - ONE), ABS(Dt(MN) - ONE) )
      WRITE(*,*) 'Orthogonality is ', ERR
*
      RETURN
*
*     End of TORTHSVD
*
      END


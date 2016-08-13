*  Definition:
*  ===========
*
*       SUBROUTINE TRANSP2DI( M, N, A, LDA, B, LDB, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            M, N, LDA, LDB, INFO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> TRANSP2D transpose matrix A(1:M,1:N) to matrix B(1:N,1:M). 
*> \endverbatim
*
*  Arguments:
*  ==========
*
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
*> \param[out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (N,LDA)
*>          The transpose of matrix A
*> \endverbatim
*
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
      SUBROUTINE TRANSP2DI( M, N, A, LDA, B, LDB, INFO )
*
*  -- HSSPACK auxiliary routine (version 0.0.1 ) --
*  -- HSSPACK is a software package provided by National Univ. of Defense Tech.    --
*  -- Novermber 2012                                                      --
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA, LDB, INFO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Executable Statements ..
*
      INFO = 0 
      IF( M .GT. LDA ) THEN
         INFO = -1
      ELSE
         IF( N .GT. LDB ) THEN
            INFO = -2
         END IF
      END IF
      IF( INFO .NE. 0 ) THEN
         WRITE(*,*) 'Input error in transp2di'
         RETURN
      END IF

      DO 10 I = 1, M
         DO 20 J = 1, N
            B( J, I ) = A( I, J )
 20      CONTINUE
 10   CONTINUE
*
      RETURN
*
*     End of TRANSP2DI
*
      END

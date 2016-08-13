*  Definition:
*  ===========
*
*       SUBROUTINE TRANSP2D( A, LDA, N, B )
* 
*       .. Scalar Arguments ..
*       INTEGER            LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, N ), B( N, LDA )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> TRANSP2D transpose matrix A to matrix B. 
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
      SUBROUTINE TRANSP2D( A, LDA, N, B )
*
*  -- HSSPACK auxiliary routine (version 0.0.1 ) --
*  -- HSSPACK is a software package provided by National Univ. of Defense Tech.    --
*  -- Novermber 2012                                                      --
*
*     .. Scalar Arguments ..
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, N ), B( N, LDA)
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Executable Statements ..
*
*
      DO 10 I = 1, LDA 
         DO 20 J = 1, N
            B( J, I ) = A( I, J )
 20      CONTINUE
 10   CONTINUE
*
      RETURN
*
*     End of TRANSP2D
*
      END

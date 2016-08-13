module BasicMM
  implicit none

contains

! ********************************************************
!           Basic routines for Cauchy-like matrices      *
! ********************************************************
!!!!!!!!
  SUBROUTINE Cauchy(A,X,Y,Z,N)
! 
! Construction an Cauchy matrix A_{ij} = z_j / (x_j - y_i)
! 
! X column basis
! Y row basis
! Z nominator
!
    REAL(8)  :: A(N,*), X(*), Y(*), Z(*)
    INTEGER  :: N       ! Leading dim of A
!
    INTEGER  :: I,J
    
    DO J = 1, N
       DO I = 1, N
          A(I,J) = Z(J) / ( X(J)-Y(I) )
       END DO
    END DO
    
  END SUBROUTINE Cauchy  

!!!!!!!!
  subroutine Cauchy2(A,D,F,u,v,M,N)
! A(i,j) = u(i)v(j) / (D(i)-F(j))

    integer M,N,i,j
    double precision A(M,*), D(*),F(*),u(*),v(*)

    do j = 1,N
       do i = 1,M
          A(i,j) = u(i)*v(j) / (D(i)-F(j))
       end do
    end do
  end subroutine Cauchy2


! ************************************************************
!         Basic routines for Cauchy-like SVD                 *
! ************************************************************

!!!!!!!
  subroutine Cauchylikesvdm(A,D,F,U,V,DIFL,DIFR,M,N)
!
! This routine is for computing the singular vectors after solving 
! some secular equations. It constructs an Cauchy-like matrix A, 
! defined as A(i,j) = U(i)V(j) / (D(i)^2-F(j)^2). 
! 
! _N.B._ It only works for M == N. 
! 
! Parameters 
! ==========
! D (input) the updated singular values 
! F (input) the old singular values 
! U (input) the normalizing scalars 
! V (input) the modified vector Z
!
! ===========
! Written by S.G. Li, Sep. 10th, 2013
! ===================

    integer M,N,i,j
    double precision DIFLJ,VJ,FJ
    double precision A(M,*),D(*),F(*),U(*),V(*),DIFL(*),DIFR(*)

!     .. External Functions ..
    double precision   DLAMC3
    external           DLAMC3

    DO J = 1, N
         DIFLJ = DIFL( J )
         VJ = V(J)
         FJ = F( J )
         A(J, J ) = U(J)*VJ / DIFLJ / ( D( J )+FJ ) 
         DO I = 1, J - 1
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I+1),-FJ )+DIFR(I) ) / ( D( I )+FJ ) 
         END DO
         DO I = J + 1, M
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I),-FJ )+DIFL(I) ) / ( D( I )+FJ ) 
         END DO
      END DO

    end subroutine Cauchylikesvdm

!!!!!!!
  subroutine CauchylikesvdU(A,D,F,U,V,DIFL,DIFR,M,N)
! A(i,j) = U(i)V(j) / (D(i)^2-F(j)^2)
! only works for M == N
! U is the normalizing scalars and V is Z
! D is the updated singular values and F is the old ones.

    integer M,N,i,j
    double precision DIFLJ,VJ,FJ
    double precision A(M,*),D(*),F(*),U(*),V(*),DIFL(*),DIFR(*)

!     .. External Functions ..
    double precision   DLAMC3
    external           DLAMC3


    DO J = 1, N-1
         DIFLJ = DIFL( J )
         VJ = V(J)
         FJ = F( J )
         A(J, J ) = U(J)*VJ / DIFLJ / ( D( J )+FJ ) 
         DO I = 1, J - 1
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I+1),-FJ )+DIFR(I) ) / ( D( I )+FJ ) 
         END DO
         DO I = J + 1, M
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I),-FJ )+DIFL(I) ) / ( D( I )+FJ ) 
         END DO
      END DO

      DO I = 1, M
         A(I,N) = U(I)
      END DO

    end subroutine CauchylikesvdU

!!!!!!!
  subroutine CauchylikesvdU1(A,D,F,U,V,DIFL,DIFR,M,N)
! A(i,j) = U(i)V(j) / (D(i)^2-F(j)^2)
! only works for M == N
! U is the normalizing scalars and V is Z
! D is the updated singular values and F is the old ones.

    integer M,N,i,j
    integer JX(M)
    double precision DIFLJ,VJ,FJ
    double precision A(M,*),D(*),F(*),U(*),V(*),DIFL(*),DIFR(*)

!     .. External Functions ..
    double precision   DLAMC3
    external           DLAMC3

    DO J = 1, N-1
         DIFLJ = DIFL( J )
         VJ = V(J)
         FJ = F( J )
         A(J, J ) = U(J)*VJ / DIFLJ / ( D( J )+FJ ) 
         DO I = 1, J - 1
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I+1),-FJ )+DIFR(I) ) / ( D( I )+FJ ) 
         END DO
         DO I = J + 1, M
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I),-FJ )+DIFL(I) ) / ( D( I )+FJ ) 
         END DO
      END DO

      DO I = 1, M
         A(I,N) = U(I)
      END DO

!  reorder matrix A
      JX(1:M) = (/ (i, i=M,1,-1) /)
      A(1:M,1:N)= A(JX,1:N)
      A(1:M,1:M)=A(1:M,JX)

    end subroutine CauchylikesvdU1

!!!!!!!
    subroutine testOrthCauchy(D,F,U,V,DIFL,DIFR,N)
!
! .. Scalar Arguments ..
      integer N
! .. Array Arguments ..
      double precision :: D(*),F(*),U(*),V(*),DIFL(*),DIFR(*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of the transpose of 
! right singular vector matrix, VT(i,j) = U_i* V_j / (D_i^2 - F_j^2),
! where D are the old singular values and W are the new singular values, 
! F(i) <= D(i). 
! 
! Parameters
! ===========
! F       (input)  double precision, the old singular values
! D       (input)  double precision, the updated singular values
! U       (input)  double precision, row generator
! V       (input)  double precision, column generator
! DIFL    (input)  double precision
! DIFR    (input)  double precision

      integer lwork,info
      double precision err
      double precision, allocatable :: A(:,:),work(:),Dt(:)

      lwork = N*(N+1)
      allocate( A(N,N),work(lwork),Dt(N) )

      call Cauchylikesvdm(A,D,F,U,V,DIFL,DIFR,N,N) 

      ! *****************************************************
      !           Check the orthogonality of V              *
      ! *****************************************************
      Dt = 0.0D0
      call dgesvd('N','N',N,N,A,N,Dt,A,N,A,N,work,lwork,info)  !svd of V1  V1 would be destroyed
      !     write(*,*) 'The svals of A are ', Dt(1:n)
      err = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
      write(*,*) "Othogonality of exact: ", err

      deallocate( A,work,Dt )

    end subroutine testOrthCauchy

!!!!!!!
    subroutine testOrthLeft(D,D1,ALPHA_L,ZZ,DIFL,DIFR,N)
!
! .. Scalar Arguments ..
      integer N
! .. Array Arguments ..
      double precision :: D(*),D1(*),ALPHA_L(*),ZZ(*),DIFL(*),DIFR(*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of left singular vector matrix.
! 
! .. Parameters ..
! D       (input)  double precision, the old singular values
! D1      (input)  double precision, the updated singular values
! ALPHA_L (input)  double precision, the scalars for the left singular vector matrix
! ZZ      (input)  double precision, the vector of the numerator
! DIFL    (input)  double precision
! DIFR    (input)  double precision

      integer lwork,info
      double precision err
      double precision, allocatable :: A(:,:),work(:),Dt(:)

      lwork = N*(N+1)
      allocate( A(N,N+1),work(lwork),Dt(N) )

      call CauchylikesvdU(A,D1,D,ALPHA_L,ZZ,DIFL,DIFR,N,N+1) 

      ! *****************************************************
      !           Check the orthogonality of V              *
      ! *****************************************************
      Dt = 0.0D0
      call dgesvd('N','N',N,N+1,A,N,Dt,A,N+1,A,N,work,lwork,info)  !svd of V1  V1 would be destroyed
      !     write(*,*) 'The svals of A are ', Dt(1:n)
      err = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
      write(*,*) "Othogonality of exact: ", err

      deallocate( A,work,Dt )

    end subroutine testOrthLeft

!!!!!!!
    subroutine testOrthLeft1(D,D1,ALPHA_L,ZZ,DIFL,DIFR,N)
!
! .. Scalar Arguments ..
      integer N
! .. Array Arguments ..
      double precision :: D(*),D1(*),ALPHA_L(*),ZZ(*),DIFL(*),DIFR(*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of left singular vector matrix.
! 
! .. Parameters ..
! D       (input)  double precision, the old singular values
! D1      (input)  double precision, the updated singular values
! ALPHA_L (input)  double precision, the scalars for the left singular vector matrix
! ZZ      (input)  double precision, the vector of the numerator
! DIFL    (input)  double precision
! DIFR    (input)  double precision

      integer lwork,info
      double precision err
      double precision, allocatable :: A(:,:),work(:),Dt(:)

      lwork = N*(N+1)
      allocate( A(N,N+1),work(lwork),Dt(N) )

      call CauchylikesvdU1(A,D1,D,ALPHA_L,ZZ,DIFL,DIFR,N,N+1) 

      ! *****************************************************
      !           Check the orthogonality of V              *
      ! *****************************************************
      Dt = 0.0D0
      call dgesvd('N','N',N,N+1,A,N,Dt,A,N+1,A,N,work,lwork,info)  !svd of V1  V1 would be destroyed
      !     write(*,*) 'The svals of A are ', Dt(1:n)
      err = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
      write(*,*) "Othogonality of exact: ", err

      deallocate( A,work,Dt )

    end subroutine testOrthLeft1

!!!!!!!
    subroutine testOrth(V,M,N)
!
! .. Scalar Arguments ..
      integer M,N
! .. Array Arguments ..
      double precision :: V(M,*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of matrix V.
! 
      integer lwork,info,MN
      double precision err
      double precision, allocatable :: A(:,:),work(:),Dt(:)

      lwork = N*(N+1)
      If( M < N) then
         MN = M
      else
         MN = N
      end If
      allocate( A(M,N),work(lwork),Dt(MN) )

      ! *****************************************************
      !           Check the orthogonality of V              *
      ! *****************************************************
      Dt = 0.0D0
      A(1:M,1:N) = V(1:M,1:N)
      call dgesvd('N','N',M,N,A,M,Dt,A,M,A,N,work,lwork,info)  !svd of V1  V1 would be destroyed
      err = max( ABS(Dt(1)-1.0D0), ABS(Dt(MN)-1.0D0) )
      write(*,*) "Othogonality of matrix V ", err

      deallocate( A,work,Dt )

    end subroutine testOrth

!!!!!!!
  SUBROUTINE Cauchy_UPSVD_VT( N, D, DSIGMA, VT, Z, INFO )
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D(*), DSIGMA(*), VT(N, *), Z(*)
!     ..
!
!  Purpose
!  =======
!
!  Cauchy_UPSVD_VT finds all the square roots of the roots of the secular
!  equation, as defined by the values in D and Z.  It makes the
!  appropriate calls to DLASD4. It returns the right singular matrix VT.
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Cauchy_SVD1 is called from updating SVD test routines.
!
!  Arguments
!  =========
!
!  N     (input) INTEGER
!         The dimension of the matrix A, which is a square matrix.
!
!  D      (output) DOUBLE PRECISION array, dimension(N)
!         On exit the square roots of the roots of the secular equation,
!         in ascending order.
!
!  DSIGMA (input) DOUBLE PRECISION array, dimension(N)
!         The singular values of previous matrix. These are the poles
!         of the secular equation.
!
!  U      (ouput) DOUBLE PRECISION array, dimension (N, N)
!         It contains the difference of the new singular values and the old ones. 
!         The aim is to be more accurate. 
!         Finally it contains the right singular vectors VT. 
!
!  VT     (workspace) DOUBLE PRECISION array, dimension (N, N)
!         First it contains the sum of new singular values and the old ones.
!         Finally it contains the right singular vectors, which is a Cauchy-like matrix V.
!         Later this will be the output. Make sure whether TRANSPOSE would be faster. 
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         It contains the appended vector and it will be updated after this routine. 
!
!  Q      (workspace) DOUBLE PRECISION array, dimension (N)
!         It is use to copy Z, let the updated Z have the same sign as Z. 
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!
!   Modified from mdlasd3 in LAPACK by Shengguo Li, on Aug. 29th, 2012.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, K
      DOUBLE PRECISION   RHO, TEMP
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION, ALLOCATABLE :: U(:,:), Q(:)
      INTEGER    JX(N)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLASCL, DLASD4
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      K  = N
      allocate( U(N,N), Q(N) )
!
!
!     Quick return if possible
!
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
      END IF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2!DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DSIGMA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO 20 I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   20 CONTINUE
!
!     Keep a copy of Z.
!
      CALL DCOPY( K, Z, 1, Q, 1 )
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Find the new singular values.
!
      DO 30 J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, U( 1, J ), RHO, D( J ), &
                     VT( 1, J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
   30 CONTINUE
!
!     Compute updated Z.
!
      DO 60 I = 1, K
         Z( I ) = U( I, K )*VT( I, K )
         DO 40 J = 1, I - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / &
                    ( DSIGMA( I )-DSIGMA( J ) ) / &
                    ( DSIGMA( I )+DSIGMA( J ) ) )
   40    CONTINUE
         DO 50 J = I, K - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / &
                    ( DSIGMA( I )-DSIGMA( J+1 ) ) /  &
                    ( DSIGMA( I )+DSIGMA( J+1 ) ) )
   50    CONTINUE
         Z( I ) = SIGN( SQRT( ABS( Z( I ) ) ), Q( I ) )
   60 CONTINUE
!
!     Generate the right singular vectors V
!
      DO 90 I = 1, K
         U( 1, I ) = Z( 1 ) / U( 1, I ) / VT( 1, I )
         DO 70 J = 2, K
            U( J, I ) = Z( J ) / U( J, I ) / VT( J, I )
   70    CONTINUE
         TEMP = DNRM2( K, U( 1, I ), 1 )
         TEMP = ONE / TEMP
         DO 80 J = 1, K
            U( J, I ) = U( J, I ) * TEMP
   80    CONTINUE
   90 CONTINUE
!       
      DO  J = 1, K
         DO I = 1, K
            VT( I,J ) = U( J,I )
         END DO
      END DO
!
    deallocate( U, Q )
!
      RETURN
!
!     End of Cauchy_UPSVD_VT
!
  END SUBROUTINE CAUCHY_UPSVD_VT

!!!!!!! Cauchy_SVD1()
  SUBROUTINE Cauchy_SVD1(N, D, DSIGMA, VT, Z, INFO )
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D(*), DSIGMA(*), VT(N, *), Z(*)
!     ..
!
!  Purpose
!  =======
!
!  Cauchy_SVD1 finds all the square roots of the roots of the secular
!  equation, as defined by the values in D and Z.  It makes the
!  appropriate calls to DLASD4
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Cauchy_SVD1 is called from updating SVD test routines.
!
!  Arguments
!  =========
!
!  N     (input) INTEGER
!         The dimension of the matrix A, which is a square matrix.
!
!  D      (output) DOUBLE PRECISION array, dimension(N)
!         On exit the square roots of the roots of the secular equation,
!         in ascending order.
!
!  DSIGMA (input) DOUBLE PRECISION array, dimension(N)
!         The singular values of previous matrix. These are the poles
!         of the secular equation.
!
!  U      (ouput) DOUBLE PRECISION array, dimension (N, N)
!         It contains the difference of the new singular values and the old ones. 
!         The aim is to be more accurate. 
!         Finally it contains the right singular vectors VT. 
!
!  VT     (workspace) DOUBLE PRECISION array, dimension (N, N)
!         First it contains the sum of new singular values and the old ones.
!         Finally it contains the right singular vectors, which is a Cauchy-like matrix V.
!         Later this will be the output. Make sure whether TRANSPOSE would be faster. 
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         It contains the appended vector and it will be updated after this routine. 
!
!  Q      (workspace) DOUBLE PRECISION array, dimension (N)
!         It is use to copy Z, let the updated Z have the same sign as Z. 
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!
!   Modified from mdlasd3 in LAPACK by Shengguo Li, on Aug. 29th, 2012.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, K
      DOUBLE PRECISION   RHO, TEMP
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION, ALLOCATABLE :: U(:,:), Q(:)
      INTEGER    JX(N)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLASCL, DLASD4
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      K  = N
      allocate( U(N,N), Q(N) )
!
!
!     Quick return if possible
!
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
      END IF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2!DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DSIGMA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO 20 I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   20 CONTINUE
!
!     Keep a copy of Z.
!
      CALL DCOPY( K, Z, 1, Q, 1 )
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Find the new singular values.
!
      DO 30 J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, U( 1, J ), RHO, D( J ), &
                     VT( 1, J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
   30 CONTINUE
!
!     Compute updated Z.
!
      DO 60 I = 1, K
         Z( I ) = U( I, K )*VT( I, K )
         DO 40 J = 1, I - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / &
                    ( DSIGMA( I )-DSIGMA( J ) ) / &
                    ( DSIGMA( I )+DSIGMA( J ) ) )
   40    CONTINUE
         DO 50 J = I, K - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / &
                    ( DSIGMA( I )-DSIGMA( J+1 ) ) /  &
                    ( DSIGMA( I )+DSIGMA( J+1 ) ) )
   50    CONTINUE
         Z( I ) = SIGN( SQRT( ABS( Z( I ) ) ), Q( I ) )
   60 CONTINUE
!
!     Generate the right singular vectors V
!
      DO 90 I = 1, K
         VT( 1, I ) = Z( 1 ) / U( 1, I ) / VT( 1, I )
         DO 70 J = 2, K
            VT( J, I ) = Z( J ) / U( J, I ) / VT( J, I )
   70    CONTINUE
         TEMP = DNRM2( K, VT( 1, I ), 1 )
         TEMP = ONE / TEMP
         DO 80 J = 1, K
            VT( J, I ) = VT( J, I ) * TEMP
   80    CONTINUE
   90 CONTINUE
!
!     Generate the right singular vectors VT, which must try 'transpose' later
!
!        VT(1:K,1:K) = transpose(Vt(1:K,1:K))

      DO 120 I = 1, K
         DO 110 J = 1, K
            U( I, J ) = VT( J, I )
  110    CONTINUE
  120 CONTINUE
            
      DO I = 1,N
         JX(I) = N-I+1
      END DO

      DO I = 1, N
         DO J = 1,N
            VT(JX(J),JX(I)) = U(J,I)
        END DO
      END DO

    deallocate( U, Q )
!
      RETURN
!
!     End of Cauchy_SVD1
!
  END SUBROUTINE CAUCHY_SVD1

!!!!!!!!
  SUBROUTINE Cauchy_svd(A,X,Y,Z,N)
! 
! Construction an Cauchy matrix A_{ij} = z_j / (x_j^2 - y_i^2)
! 
! X column basis
! Y row basis
! Z nominator
!
! This procedure is design for updating SVD problem 
!===============
! Written by Shengguo Li, On Aug. 16th, 2012
!===============

    REAL(8)  :: A(N,*), X(*), Y(*), Z(*), temp
    INTEGER  :: N       ! Leading dim of A
!
    INTEGER  :: I,J
    
    DO J = 1, N
       DO I = 1, N
          temp = ( X(J)-Y(I) ) *( X(J) + Y(I) )
          A(I,J) = Z(J) / temp
       END DO
    END DO
    
  END SUBROUTINE Cauchy_svd 

!!!! NormlzM
  subroutine NormLzM(A, N, RC)
!
! .. Scalar Arguments ..
    INTEGER   :: N
    CHARACTER(LEN = 3) :: RC

! .. Array Arguments ..
    DOUBLE PRECISION  :: A(N,N)
!
! Purpose
! ========
! Normalize the column or row of matrix A. 
!
! .. Parameters .. 
! A  (INPUT/OUTPUT)  DOUBLE PRECISION Array, DIMENSION(N,N) 
! 
! N  (INPUT)  INTEGER
!
! RC (INPUT)  CHARACTER STRING, 'ROW' or 'COL'.
!
! =================
! Written by Shengguo Li, On Aug. 16th, 2012
! =================

    INTEGER :: I
    DOUBLE PRECISION  :: TEMP, dnrm2

!  .. External Functions
      LOGICAL     LSAME, GR, GC
      EXTERNAL    LSAME, dnrm2

    GR = LSAME(RC, 'ROW')
    GC = LSAME(RC, 'COL')

    IF( .not.(GR .or. GC) ) THEN
       print *, 'Normalization must specify row or col. '
    END IF
    
    IF( GC ) THEN    
       DO I =1, N
          temp = dnrm2(N,A(1:N, I),1 )
          temp = 1.0D0 / temp
          A(1:N,I) = temp * A(1:N,I)
       END DO
    ELSE
       DO I =1, N
          temp = dnrm2(N,A(I, 1:N),1 )
          temp = 1.0D0 / temp
          A(I, 1:N) = temp * A(I, 1:N)
       END DO
    END IF
       
  end subroutine NormLzM

!!!!!!
      subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************
!
!! R8MAT_PRINT prints a R8MAT.
!
!  Discussion:
!
!    A R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return

end subroutine r8mat_print

!!!!!!!!
   subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************
!
!! R8MAT_PRINT_SOME prints some of a R8MAT.
!
!  Discussion:
!
!    A R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end subroutine r8mat_print_some
   
!!!!!!!!!!!!
   subroutine r8vec_print ( n, a, title )

!*****************************************************************************
!
!! R8VEC_PRINT prints a R8VEC.
!
!  Discussion:
!
!    A R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end subroutine r8vec_print

!!!!!!!!!!!
    subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices to print.
!    The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i_hi
  integer   ( kind = 4 ) i_lo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
  end do

  return
end subroutine r8vec_print_some

!!!!!!!!
   subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end subroutine i4vec_print

!!!!!!!!!
  subroutine i4vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! I4VEC_PRINT_SOME prints "some" of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices to print.
!    The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i_hi
  integer   ( kind = 4 ) i_lo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,i8)' ) i, a(i)
  end do

  return
end subroutine i4vec_print_some

!!!!!!
subroutine init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
          
  CALL SYSTEM_CLOCK(COUNT=clock)
          
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
          
  DEALLOCATE(seed)

end subroutine init_random_seed

!
!!!!!!!!!
   SUBROUTINE DSymSplit( AB,LDAB,SB,U,LDU,VT,LDVT,W,IC,WORK,INFO )
!
      IMPLICIT NONE 
! ..
! .. Scalars arguments..
      INTEGER               LDAB, SB, LDU, LDVT, IC, INFO
! ..
! .. Array arguments ..
      DOUBLE PRECISION      W(*), WORK(*)
      DOUBLE PRECISION      AB(LDAB,*), U(LDU,*), VT(LDVT,*)
!
! Purpose
! =======
!
!   DSymSplit computes the tilling block-diagonal matrix by extracting two small 
!   symmetric blocks, A_i = A_i-V_{i-1}*\Sigma_{i-1}*V_{i-1}^T-U_i*\Sigma_i*U_i^T. 
!
!   Note that AB is a symmetric band matrix and assumed to be stored in upper 
!   triangular form. 
!
! Arguments
! =========
!
!   AB     (inout) double precision array, dimension ( LDAB, N )
!           The banded matrix A in upper banded storage (lower
!           banded storage is not supported).
!
!   LDAB   (in) integer, the leading dimension of AB
!
!   SB     (in) integer, the semi-bandwidth of this banded matrix AB.
!
!   U      (in) double precision array, dimension ( LDU,SB )
!           The left singular vector matrix for B_i.
!
!   LDU    (in) integer, the leading dimension of the array U.
!           LDFULL >= MAX( SB, 1 ).
!
!   VT     (in) double precision array, dimension ( LDVT, SB )
!           The left singular vector matrix of B_i.
!
!   LDVT   (in) integer, the leading dimension of the array VT.
!
!    W     (in) double precision array, dimension(SB)
!            The singualr values of B_i.
!
!   IC     (in) integer, the ending position of current block, A_{i}
!
!   WORK   (inout) double precision array, dimension(2*SB*SB)
!          Double precision workspace, its size should be larger than 2*SB**2
!
!   INFO    (out) integer
!           On exit, info indicates the consistency of the arguments.
!
! ===================
! Author: Shengguo Li, Nati. Univ. of Def. Tech., China
! email:  nudtlsg@gmail.com
!
! Date: Oct. 19, 2014, Version: HSSPACK, V-1.0.0
! 
! ======================================================
!
!  Local variables:
   INTEGER  I, J, ISU, Imat, start,end, IJ, abi, abj
   DOUBLE PRECISION, PARAMETER :: ZERO=0.0D0,ONE=1.0D0  
!
!  External routines
   EXTERNAL DGEMM 
!

   ISU = 1               ! Scaled matrix U or VT
   Imat = ISU + SB*SB    ! The split matrix 
!  Update the bottom right blcok
   start = ISU 
   DO I = 1, SB
     end = start + SB-1
     WORK(start:end) = U(1:SB,I)*W(I)
     start = end + 1   
   END DO

   CALL DGEMM( 'N','T',SB,SB,SB,ONE,WORK(ISU),SB,U,SB,ZERO,WORK(Imat),SB )

!  Compute A_i = A_i - U_i \Sigma_i U_i^T
   DO J = 1, SB
     abj = IC-SB+J
     DO I = 1, J
       IJ = Imat+SB*(J-1)+I-1 
       abi = SB+1+I-J
       AB(abi,abj) = AB(abi,abj) - WORK( IJ )
     END DO
   END DO 
   

!  Update the top right block of next block diagonal block
   start = ISU
   DO I = 1, SB
     end = start + SB-1
     WORK(start:end) = VT(I,1:SB)*W(I)
     start = end + 1   
   END DO

   CALL DGEMM('N','N',SB,SB,SB,ONE,WORK(ISU),SB,VT,SB,ZERO,WORK(Imat),SB )

!  Compute A_{i+1} = A_{i+1} - V_i \Sigma_i V_i^T
   DO J = 1, SB
      abj = IC+J
      DO I = 1, J
        IJ = Imat + SB*(J-1)+I-1
        abi = SB+1+I-J
        AB(abi,abj) = AB(abi,abj) - WORK(IJ)
      END DO
   END DO

   END SUBROUTINE DSymSplit
!
!

end module BasicMM

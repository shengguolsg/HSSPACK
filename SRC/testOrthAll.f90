   subroutine testOrthAll(NRM,V,M,N)
     use BasicMM
     implicit none
!
! .. Scalar Arguments ..
      integer    :: M,N
      character( len= 1 ) :: NRM
! .. Array Arguments ..
      double precision :: V(M,*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of matrix V.
! NRM can be 
!   'F':  Frobenius norm, ||V*V^T -I ||_F, row orthogonality, or
!                         ||V^T*V -I ||_F, col orthogonality;
!   'T':  2-norm, ||V*V^T -I ||_2, row orthogonality, or
!                 ||V^T*V -I ||_2, col orthogonality;
! 
! =========
! Written by S.G. Li, on Sept. 11th, 2013
!
! ==============
!
!     .. External Functions ..
      LOGICAL   LSAME
!
      integer lwork, info, MN, INORM, I, J
      double precision err
      double precision, parameter :: ONE = 1.0E+0, ZERO=0.0E+0
!
      double precision, allocatable :: A(:,:),work(:),Dt(:)
!
      lwork = N*(N+1)
      IF( LSAME( NRM, 'F' ) )  INORM = 1
      IF( LSAME( NRM, 'T' ) )  INORM = 2     

      If( M < N) then   ! Row orth
         MN = M
         allocate( A(M,M),work(lwork),Dt(MN) )
         call dgemm('N','T',M,M,N,ONE,V,M,V,M,ZERO,A,M)
         DO I = 1, M
            A(I,I) = A(I,I)-ONE
         END DO
         
         IF( INORM == 1 ) THEN
            err = ZERO
            DO I = 1, M
               DO J = 1, M
                  err = err + A(J,I)*A(J,I)
               END DO
            END DO
            err = sqrt(err)            
            write(*,*) 'Frobenius norm is ', err
         ELSE
            call testorth(A,M,M)
         END IF         
      else   ! Col orth
         MN = N
         allocate( A(N,N),work(lwork),Dt(MN) )
         call dgemm('T','N',N,N,M,ONE,V,M,V,M,ZERO,A,M)
         DO I = 1, M
            A(I,I) = A(I,I)-ONE
         END DO
         
         IF( INORM == 1 ) THEN
            err = ZERO
            DO I = 1, N
               DO J = 1, N
                  err = err + A(J,I)*A(J,I)
               END DO
            END DO
            err = sqrt(err)
            write(*,*) 'Frobenius norm is ', err
         ELSE
            call testorth(A,N,N)
         END IF         

      end If
      
      deallocate( A,work,Dt )

    end subroutine testOrthAll


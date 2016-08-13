  program testall1
!
    implicit none
! 
! Purpose 
! =========
! This program tests all the matrices in Example 2 of our "Updating SVD Paper". 
! For a N-by-N random matrices, update it several times and compute its SVD.
! In this example, we update twice, four times, 6, 8 and 10 times.
!
! Tests:
!    1) The error of computed singular values compared with compute it directly from the updated matrix
!    2) The orthogonality of computed right singular vector matrix w/o HSS approximation

  INTEGER :: I,J,N,LWORK,ierr,info,N2,ni,mpk,K, lda,Tims,I2,prd
  DOUBLE PRECISION :: time, time1, temp, maxerr,tol
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:),D(:),F(:),V(:,:),Work(:),X(:), &
                              D1(:),AA(:,:),Z(:), Dt(:), A0(:,:),svs(:), serr2(:),Vt(:,:)
! A  : the inital N-by-N matrix;
! D  : the singular values of A;
! Dt : the singular vlaues of V1, used to check the orthogonality;
! W  : the singular values of updated matrix A1; 
! V  : the right singular vector matrix of A;
! V1 : the right singular vector matrix of A1; 
! X  : the appended random vector;
! Z  : the appended vector expressed as a combination of V;  Z=V^T*X
! svs: the singular values updated matrices
! serr2: the maximum errors between computed directly and via updating with HSS
! Work: workspace for calling dgesvd.f;
! 
! ==========
! Modified by S.G. Li, on Sep. 15th, 2013
! =======================================
!
  DOUBLE PRECISION :: ABS
  INTRINSIC        :: ABS
  !
  info = 0
  N = 600
  prd = 2
  tol = 1.0D-16
  Tims = 2
  lda = N+2*Tims
  LWORK = (N/2)*(N+1)
  allocate(A(N+2*Tims,N),D(N),F(N),V(N,N),Work(LWORK),serr2(2*Tims), STAT = ierr )
  IF(ierr /= 0) Then
     write(*,*) "Allocate Error."
     stop
  END IF
  allocate(A0(N+2*Tims,N),AA(N,N),X(N),Z(N),Dt(N),SVS(Tims*N),Vt(N,N),D1(N), STAT = ierr) 
  IF(ierr /= 0) Then
     write(*,*) "Allocate Error."
     stop
  END IF
     
  call random_number(A)
  call dlacpy('A',lda,N,A,lda,A0,lda)
  
  ! compute the singular values of each updated matrix
  DO K =1, Tims
     mpk = N+2*k
     call dgesvd('N','N',mpk,N,A0,lda,svs( 1+N*(k-1) ),A0,N,A0,N,work,lwork,info)  ! svd of updated A
     IF(info /= 0) write(*,*) "Error while solving svals."
     call dlacpy('A',lda,N,A,lda,A0,lda)
  END DO ! K

  call dlacpy('A',N,N,A,lda,AA,N)
  call dgesvd('N','A',N,N,AA,N,D,AA,N,V,N,work,lwork,info)  !svd of A(N,N), D
  IF(info /= 0) Then
     write(*,*) "dgesvd is not wrong. Info= ", info, "I= ", I
     stop
  END IF

  call dlacpy('A',lda,N,A,lda,A0,lda)
  N2 = N/2
  
  DO I =1, 2*Tims        ! updating times

     ! ******************** Computed via HSS matrices ****************
     X(1:N) = A0(N+I,1:N)
     F = D
     call dgemv('N',N,N,1d0,V,N,X,1,0d0,Z,1)   ! appended vector Z
     call dlasrt('I', N, F, info)      ! reverse F
     Do J =1, N2                       ! reverse Z
        temp = Z(J)
        Z(J) = Z(N-J+1)
        Z(N-J+1) = temp
     END DO
  
     ni = 50
     call cpu_time(time)
     call SvdHssUpdat( D,F,Z,V,N,TOL,ni )
     call cpu_time(time1)
     time = time1 - time
!     write(*,*) "hss total time = ", time

! ************** Check the error if I is even *****************
     IF(Mod(I,prd) .eq. 0) Then

        call dlasrt('d', N, D, info)      ! reverse D
        I2 = I/2
        serr2(I2) = maxval( abs( D(1:N)-svs(1+N*(I2-1):I2*N) ) )
        serr2(I2+Tims) = maxval( abs(D(1:N)-svs(1+N*(I2-1):I2*N)) / abs(svs(1+N*(I2-1):I2*N)) )
        
        call dlacpy('A',N,N,V,N,Vt,N)  ! Vt = V
        write(*,*) 'Frobenius Norm of HSS:'
        call testOrthAll('F',V,N,N)

     END IF

  END DO ! loop of updating times

  write(*,*) "The Max error of svals via HSS: "
  write(*,*) serr2(1:Tims/prd+1)
  write(*,*) "The Rel error of svals via HSS: "
  write(*,*) serr2(Tims+1:Tims+Tims/prd+1 )
  
  Deallocate(A, D, V, Work, X, Z, Dt, Vt) 

  end program testall1

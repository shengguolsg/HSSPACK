  program testUpsvdCauchy
    use CompSvdup
    use BasicMM

    implicit none
! 
! This program tests the accuracy of the updating SVD of using structured Cauchy-like
! HSS construction algorithm. This codes would be included as a test routine of 
! HSSPACK. 
!
! ======================
! Written by Shengguo Li, on Dec. 13th, 2012
! 
! Modified by S.G. Li, on Sep. 9th, 2013
! =====================
!
  INTEGER :: I, J, N, LWORK, ierr, info,N2, ni
  INTEGER :: DIM(8)
  DOUBLE PRECISION :: time, time1, temp, tol, maxerr,sa,sb,gap
  DOUBLE PRECISION, ALLOCATABLE :: D(:),V(:,:),V1(:,:),Work(:),X(:),&
                                   Z(:),Dt(:),F(:),A(:,:),FF(:),ZZ(:)
! A  : the inital N-by-N matrix;
! F  : the singular values of A;
! Dt : the singular vlaues of V1, used to check the orthogonality;
! D  : the singular values of updated matrix A1; 
! V  : the right singular vector matrix of A;
! V1 : the right singular vector matrix of A1; 
! X  : the appended random vector;
! Z  : the appended vector expressed as a combination of V;  Z=V^T*X
! Work: workspace for calling dgesvd.f;

  DOUBLE PRECISION :: ABS
  INTRINSIC        :: ABS

    DOUBLE PRECISION ONE, ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
  
  data DIM /299, 1000, 600, 1200, 1400, 1800, 3000, 5000/
  
  DO I =1, 2
     N = DIM(I)
     print *, "The reuslts for Dim = ", N
     LWORK = (N/2)*(N+1)

     allocate(D(N),F(N),FF(N),V(N,N),V1(N,N),Work(LWORK), STAT = ierr )
     IF(ierr /= 0) Then
        write(*,*) "Allocate Error1."
        stop
     END IF
     allocate( A(N,N),X(N),Z(N),Dt(N),ZZ(N), STAT = ierr)
     IF(ierr /= 0) Then
        write(*,*) "Allocate Error2."
        stop
     END IF
     
     call init_random_seed()
     call random_number( X )
     call init_random_seed()
     call random_number( Z )
     Z(10) = 1.0E-16
     Z(200) = 1.0E-17

! *****************************************************
!           Check the orthogonality of V              *
! *****************************************************
     sa = 0.1D0
     sb = 9.0D0
     gap = (sb-sa)/ (2*N)
     DO j = 1, N
        F(j) = sa + (2*j-1)*gap
     END DO

     FF( 1:N ) = F( 1:N )
     ZZ( 1:N ) = Z( 1:N )
     V = zero
     DO j = 1, N
        V( j,j ) = one
     END DO

! ************************************************************
!          Matrix-Matrix Multiplication                      *
! ************************************************************
     info = 0
     call cpu_time(time)
     call Cauchy_UPSVD_VT( N,D,FF,V1,ZZ,INFO )
     call cpu_time(time1)
     time = time1 -time

     A = V1
     call dgesvd('N','N',N,N,A,N,Dt,A,N,V,N,work,lwork,info)  !svd of V1  V1 would be destroyed
     temp = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
     write(*,*) "Othogonality of original : ", temp 

! *********************************************************
!                   HSS method                            *
! *********************************************************
     ni = 50
     tol = 1.0D-15
     call cpu_time(time)
     call SvdHssUpdat( D,F,Z,V,N,TOL,ni )
     call cpu_time(time1)
     time = time1 - time
     write(*,*) "hss total time = ", time 

! *********************************************************
!             Check the 1-norm error of HSS               *
! *********************************************************
     Dt = 0.0D0
     A = V
     call dgesvd('N','N',N,N,A,N,Dt,A,N,A,N,work,lwork,info)  !svd of V
     temp = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0d0) )
     write(*,*) "Othogonality of HSS: ", temp

     Deallocate( D,F,V,V1,Work,A,X,ZZ,Z,Dt,FF )

  END DO ! loop of matrices

end program testUpsvdCauchy

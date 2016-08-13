  program testUpsvdHss
!    use CompSvdup
!    use BasicMM

    implicit none
! 
! This program tests all the matrices in Example 1 of our "Updating SVD Paper". 
! For random matrices with different dimensions, compute all the singular values of the 
! updated matrices by appending one more row at the bottom. 
!
! ======================
! Written by Shengguo Li, on Aug. 25th, 2012
!
! Modified on Aug. 26th, 2012
! =====================
!
  INTEGER :: I, J, N, LWORK, ierr, info,N2, ni
  INTEGER :: DIM(8)
  DOUBLE PRECISION :: time, time1, temp, tol, maxerr,sa,sb,gap
  DOUBLE PRECISION, ALLOCATABLE :: D(:),V(:,:),V1(:,:),Work(:),X(:),&
                                   Z(:),Dt(:),F(:),A(:,:),C(:,:)
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
  INTRINSIC         :: ABS
  Character(len =1) :: RS

    DOUBLE PRECISION ONE, ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
  
  data DIM /1499, 10000, 600, 1200, 1400, 1800, 3000, 5000/
  
  DO I =1, 1
     N = DIM(I)
     print *, "The reuslts for Dim = ", N
     LWORK = (N/2)*(N+1)

     allocate(D(N),F(N),V(N,N),V1(N,N),Work(LWORK), STAT = ierr )
     IF(ierr /= 0) Then
        write(*,*) "Allocate Error1."
        stop
     END IF
     allocate( X(N),Z(N),Dt(N),A(N,N),C(N,N), STAT = ierr)
     IF(ierr /= 0) Then
        write(*,*) "Allocate Error2."
        stop
     END IF
     
!     call init_random_seed()
     call random_number(A)
!     call init_random_seed()
     call random_number(x)

     call cpu_time(time)
     call dgesvd('N','A',N,N,A,N,F,A,N,V,N,work,lwork,info)  !svd of A. 
     call cpu_time(time1)
     time = time1 - time
     write(*,*) "SVD time = ", time

     IF(info /= 0) Then
        write(*,*) "dgesvd is not wrong. Info= ", info, "I= ", I
        stop
     END IF

! *****************************************************
!           Check the orthogonality of V              *
! *****************************************************
     sa = 0.1D0
     sb = 9.0D0
     gap = (sb-sa)/ (2*N)
     DO j = 1, N
        F(j) = sa + (2*j-1)*gap
     END DO
!
     call dgemv('N',N,N,1d0,V,N,X,1,0d0,Z,1)  ! appended vector Z
     call dlasrt('I', N, F, info)      ! order F in increasing order
     N2 = N/2
     Do J =1, N2                       ! reverse Z
        temp = Z(J)
        Z(J) = Z(N-J+1)
        Z(N-J+1) = temp
     END DO

! ************************************************************
!          Matrix-Matrix Multiplication                      *
! ************************************************************
     X = Z
     call Cauchy_SVD1(N,D,F,A,Z,INFO)
     call cpu_time(time)
     call dgemm('N','N',N,N,N,1d0,A,N,V,N,0d0,V1,N )   ! V1 is the CORRECT right singular vector matrix V^T
     call cpu_time(time1)
     time = time1 -time
     write(*,*) "MMM time = ", time

     C = V1
     call dgesvd('N','N',N,N,C,N,Dt,A,N,V,N,work,lwork,info)  !svd of V1  V1 would be destroyed
     temp = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
     write(*,*) "Othogonality of original : ", temp 

! *********************************************************
!                   HSS method                            *
! *********************************************************
     ni = 50
     tol = 1.0D-15
     Z = X
     RS = 'R'
     call cpu_time(time)
     call SvdHssup( A, N, V, RS, ni )
     call cpu_time(time1)
     time = time1 - time
     write(*,*) "hss total time = ", time 

! *********************************************************
!             Check the 1-norm error of HSS               *
! *********************************************************
     Dt = 0.0D0
     A = V
     call dgesvd('N','N',N,N,A,N,Dt,A,N,V,N,work,lwork,info)  !svd of V
     temp = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0d0) )
     write(*,*) "Othogonality of HSS: ", temp 

     V1 = V1 - V
     V1 = ABS(V1)
     maxerr = maxval(V1)
     print*, "The maximum error is ", maxerr

     Deallocate(D,F,V,Work,Z,Dt,X,A)
     
  END DO ! loop of matrices

end program testUpsvdHss

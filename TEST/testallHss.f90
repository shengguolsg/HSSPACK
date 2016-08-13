  program testAllHss
    use ConstructHssd
    use CompSvdup
    use BasicMM

    implicit none
! 
! This program tests all the matrices in Example 1 of our "Updating SVD Paper". 
! For random matrices with different dimensions, compute all the singular values 
! of the updated matrices by appending one row at the bottom. This routine uses 
! the usual HSS construction algorithm.
!
! ======================
!
! Written by Shengguo Li, on Aug. 25th, 2012
!
! Modified on Aug. 26th, 2012
! =====================
!
  INTEGER :: I, J, N, LWORK, ierr, info,N2, ni
  INTEGER :: DIM(5)
  DOUBLE PRECISION :: time, time1, temp, maxerr
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), D(:), W(:), V(:,:), V1(:,:), Work(:), X(:), &
                                   Z(:), Dt(:), Vt(:,:)
! A  : the inital N-by-N matrix;
! D  : the singular values of A;
! Dt : the singular vlaues of V1, used to check the orthogonality;
! W  : the singular values of updated matrix A1; 
! V  : the right singular vector matrix of A;
! V1 : the right singular vector matrix of A1; 
! X  : the appended random vector;
! Z  : the appended vector expressed as a combination of V;  Z=V^T*X
! Work: workspace for calling dgesvd.f;

  DOUBLE PRECISION :: ABS
  INTRINSIC        :: ABS
  
  data DIM /200, 600, 1000, 1400, 1800/
  
  DO I =2, 2
     N = DIM(I)
     print *, "The reuslts for Dim = ", N
     LWORK = (N/2)*(N+1)

     allocate(A(N,N), D(N), W(N), V(N,N), V1(N,N), Work(LWORK), STAT = ierr )
     IF(ierr /= 0) Then
        write(*,*) "Allocate Error1."
        stop
     END IF
     allocate(X(N), Z(N), Dt(N),Vt(N,N), STAT = ierr)
     IF(ierr /= 0) Then
        write(*,*) "Allocate Error2."
        stop
     END IF
     
     call random_number(A)
     call random_number(x)
     
     call dgesvd('N','A',N,N,A,N,D,A,N,V,N,work,lwork,info)  !svd of A. 
     IF(info /= 0) Then
        write(*,*) "dgesvd is not wrong. Info= ", info, "I= ", I
        stop
     END IF

! *****************************************************
!           Check the orthogonality of V              *
! *****************************************************
     call dlacpy('A',N,N,V,N,Vt,N)  ! Vt = V
     call dgesvd('N','N',N,N,Vt,N,Dt,A,N,Vt,N,work,lwork,info)  !svd of V1  V1 would be destroyed
     temp = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
     write(*,*) "Othogonality of exact: ", temp

     call dgemv('N',N,N,1d0,V,N,X,1,0d0,Z,1)  ! appended vector Z
     call dlasrt('I', N, D, info)      ! reverse D
     N2 = N/2
     Do J =1, N2                       ! reverse Z
        temp = Z(J)
        Z(J) = Z(N-J+1)
        Z(N-J+1) = temp
     END DO

     call cpu_time(time)
     call Cauchy_SVD1(N,W,D,A,Z,INFO)
     call cpu_time(time1)
     write(*,*) "Constructing Matrix A ", time1-time
     
! ************************************************************
!          Check the orthogonality of A                      *
! ************************************************************
     call dlacpy('A',N,N,A,N,Vt,N)  ! Vt = A
     call dgesvd('N','N',N,N,Vt,N,Dt,A,N,Vt,N,work,lwork,info)  !svd of Vt  Vt would be destroyed
     temp = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
     write(*,*) "Othogonality of Cauchy-like: ", temp

     call cpu_time(time)
     call dgemm('N','N',N,N,N,1d0,A,N,V,N,0d0,V1,N )   ! V1 is the CORRECT right singular vector matrix V^T
     call cpu_time(time1)
     time = time1 -time
     write(*,*) "mm time = ", time

! *********************************************************
!                   HSS method                            *
! *********************************************************
     ni = 50
     call cpu_time(time)
     call SvdHssup(A, N, V, 'R', ni )
     call cpu_time(time1)
     time = time1 - time
     write(*,*) "hss total time = ", time

! *********************************************************
!             Check the 1-norm error of HSS               *
! *********************************************************
     A= V1-V
     A = ABS(A)
     maxerr = maxval(A)
     print*, "The maximum error is ", maxerr
     A = ABS(A/V1)
     maxerr = maxval(A)
     print*, "The maximum relative error is ", maxerr
     
     ! Compute the orthogonality
     write(*,*) 'Norm of exact:'
     call testOrthAll('F',V1,N,N)
!     call testOrthAll('T',V1,N,N)
     write(*,*) 'Norm of HSS:'
     call testOrthAll('F',V,N,N)
!     call testOrthAll('T',V,N,N)

     Deallocate(A, D, W, V, V1, Work, X, Z, Dt, Vt) 

  END DO ! loop of matrices

  end program testAllHss

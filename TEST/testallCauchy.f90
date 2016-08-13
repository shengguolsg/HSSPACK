  program testAllCauchy
    use aux_hss
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
! Modified by S.G. Li, on Sep. 11th, 2013
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

  DOUBLE PRECISION ONE, ZERO, GG
  PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
  
  data DIM /299, 1000, 600, 1200, 1400, 1800, 3000, 5000/
  
  DO I =1, 1
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
     call random_number( A )
     call init_random_seed()
     call random_number( X )

     call dgesvd('N','A',N,N,A,N,F,A,N,V,N,work,lwork,info)  !svd of A. 
     IF(info /= 0) Then
        write(*,*) "dgesvd is not wrong. Info= ", info, "I= ", I
        stop
     END IF

! *****************************************************
!           Check the orthogonality of V              *
! *****************************************************
     call dgemv('N',N,N,1d0,V,N,X,1,0d0,Z,1)  ! appended vector Z
     ZZ = Z
     FF = F
     call dlasrt('I', N, FF, info)      ! reverse D
     N2 = N/2
     Do J =1, N2                       ! reverse Z
        temp = ZZ(J)
        ZZ(J) = ZZ(N-J+1)
        ZZ(N-J+1) = temp
     END DO

     time = comp_time()
     call Cauchy_SVD1(N,D,FF,A,ZZ,INFO)
     time1 = comp_time()
     write(*,*) "Constructing Matrix A ", time1-time

! ************************************************************
!          Matrix-Matrix Multiplication                      *
! ************************************************************
     write(*,*) "Othogonality of Cauchy-like: "
     call testorth(A,N,N)

     time = comp_time()   
     call dgemm('N','N',N,N,N,1d0,A,N,V,N,0d0,V1,N )   ! V1 is the CORRECT right singular vector matrix V^T
     time1 = comp_time()  
     time = time1 -time
     GG = 2*N**3 / time * 1e-9
     write(*,*) "mm time = ", time, GG, 'GB'

! *********************************************************
!                   HSS method                            *
! *********************************************************
     ni = 50
     tol = 1.0D-15
     time = comp_time()
     call SvdHssUpdat( D,FF,ZZ,V,N,TOL,ni )
     time1 = comp_time()
     time = time1 - time
     write(*,*) "hss total time = ", time 

! *********************************************************
!             Check the 1-norm error of HSS               *
! *********************************************************
     A= V1+V
     A = ABS(A)
     maxerr = maxval(A)
     print*, "The maximum error is ", maxerr
     A = ABS(A/V1)
     maxerr = maxval(A)
     print*, "The maximum relative error is ", maxerr

     ! Test orthogonality of V and V1
     write(*,*) 'Norm of exact:'
     call testOrthAll('F',V1,N,N)
!     call testOrthAll('T',V1,N,N)
     write(*,*) 'Norm of HSS:'
     call testOrthAll('F',V,N,N)
!     call testOrthAll('T',V,N,N)

     Deallocate( D,F,V,V1,Work,A,X,ZZ,Z,Dt,FF )

  END DO ! loop of matrices

end program testAllCauchy

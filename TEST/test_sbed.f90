  program Test_sbed

  implicit none
!
! This program tests correctness of mdsbedc which implements the 
! banded DC algorithm for symmetric banded matrices.

  integer(kind = 4) :: i,N,lwork,info,ierr,SB,LDAB,smlsiz,liwork
  integer           :: ii
  double precision  :: time,time1,one,zero,err
  parameter( one = 1.0E0, zero = 0.0E0 )
 
  integer, allocatable ::  IWORK(:)
  double precision, allocatable :: AB(:,:),Q(:,:),WORK(:),S(:),&
       D(:),Q1(:,:),AB0(:,:)
  
  integer,dimension(5) :: DIM= (/ 2000, 4000, 6000, 8000, 12000 /) 

 DOUBLE PRECISION   comp_Wtime
 EXTERNAL           comp_Wtime

  do ii = 3, 3 

  SB = 5
  N = DIM(II)
  LDAB = SB + 1
  smlsiz = 50

  write(*,*) 'The dimension of matrix A is', N  

  ! generate the diagonal and off-diagonal elements
  ALLOCATE( S(N),AB(LDAB,N),AB0(LDAB,N),Q1(N,N), stat=ierr )
  IF( ierr .ne. 0 ) THEN
     WRITE(*,*) 'Allocate failed in testdc.'
  END IF
  
  ! singualr vectors and workspaces
  lwork = 2*N*N + 7*N
  liwork = 6*N  ! 5*N+3
  ALLOCATE( Q(N,N),D(N),WORK(lwork),IWORK(liwork), stat=ierr )
  IF( ierr .ne. 0) THEN
     WRITE(*,*) 'Wrong while allocating Svecs.'
  END IF
  
  ! generate test matrix 
  call random_number( AB ) 
  AB0 = AB 
!  Call dlaset( 'L',SB,SB,ZERO,ZERO,AB(2,1),LDAB )

  ! Our algorithm 
   time = comp_Wtime()
!  call mdsbedc( 'V','U',N,SB,AB,LDAB,D,Q,N,WORK,LWORK,IWORK,LIWORK,INFO )
  call mdsbevd( 'I','U',N,SB,AB,LDAB,D,Q,N,WORK,LWORK,IWORK,LIWORK,INFO )
  time1 = comp_Wtime()
  time = time1-time
  write(*,*) 'Computing time of mdsbevd in HSSPACK is ', time

!  write(*,*) 'The eigenvalues computed by mdsbedc', D(1:N)
!  call TestOrth( Q,N,N,N )

!
  ! LAPACK, dsbevd
  lwork = 2*N*N + 7*N
  liwork = 5*N + 3
  time = comp_Wtime()
  call dsbevd( 'V','U',N,SB,AB0,LDAB,S,Q1,N,WORK,LWORK,IWORK,LIWORK,INFO )
  time1 = comp_Wtime()
  time = time1-time
  write(*,*) 'Computing time of dsbevd in MKL is ', time 

!  write(*,*) 'The eigenvalues computed by mdsbedc', S(1:N)
!
  DO i = 1, N
     D(i) = D(i)-S(i)
  END DO
  err = maxval( abs(D) )
  write(*,*) 'The max error is ', err
  err = maxval( abs(D)/abs(S) )
  write(*,*) 'The max relative error is', err
  write(*,*)
!     
  DEALLOCATE( AB,S,AB0,Q,Q1,D,WORK,IWORK )

 END DO
!
end program Test_sbed

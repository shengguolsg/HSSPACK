  program test_bdc1

  implicit none
!
! This program tests correctness of mdbdsvd which implements the 
! banded DC algorithm. 
!

  integer(kind = 4) :: i,j,k,N,lwork,info,ierr,SB,LDAB,smlsiz
  double precision  :: time,time1,one,zero,err
  parameter( one = 1.0E0, zero = 0.0E0 )
 
  integer, allocatable ::  IWORK(:)
  double precision, allocatable :: A(:),AB(:,:),U(:,:),VT(:,:),&
       WORK(:),S(:),D(:),C(:,:)

  DOUBLE PRECISION   DNRM2, comp_Wtime
  EXTERNAL           DNRM2, comp_Wtime


  SB = 5
  N = 5000
  LDAB = SB + 1
  smlsiz = 75
  
  write(*,*) 'DIM = ', N

  ! generate the diagonal and off-diagonal elements
  ALLOCATE( S(N),AB(LDAB,N),A(N*N),C(N,N), stat=ierr )
  IF( ierr .ne. 0 ) THEN
     WRITE(*,*) 'Allocate failed in testdc.'
  END IF
  
  ! singualr vectors and workspaces
  lwork = 5*N*N + 7*N
  ALLOCATE( U(N,N),D(N),VT(N,N),WORK(lwork),IWORK(8*N), stat=ierr )
  IF( ierr .ne. 0) THEN
     WRITE(*,*) 'Wrong while allocating Svecs.'
  END IF
  
  ! generate test matrix 
  call random_number( AB ) 
!  Call dlaset('L',SB,SB,ZERO,ZERO,AB(2,1),LDAB )

  ! Our algorithm 
  time  = comp_Wtime()
  call mdbdsvd( N,SB,AB,LDAB,A,S,U,N,VT,N,SMLSIZ,WORK,LWORK,IWORK,INFO )
  time1 = comp_Wtime()
  time  = time1-time
  write(*,*) 'Computing time of mdbdsvd in HSSPACK is ', time 

  ! LAPACK, dgesvd
  call dsb2flc( 'U',N,N,SB,AB,LDAB,C,N,INFO )
  write(*,*)   'Info is ', info, 'SB=', SB
  lwork = 4*N*N + 7*N
  time  = comp_Wtime()
  call dgesdd('A',N,N,C,N,D,U,N,VT,N,WORK,LWORK,IWORK,INFO)
  time1 = comp_Wtime()
  time  = time1-time
  write(*,*) 'Computing time of dgesdd in MKL is ', time 


  DO i    = 1, N
     D(i) = D(i)-S(i)
  END DO
  err = maxval( abs(D) )
  write(*,*) 'The max error is ', err
     
  DEALLOCATE( AB,S,A,U,VT,WORK,IWORK,C )

end program test_bdc1

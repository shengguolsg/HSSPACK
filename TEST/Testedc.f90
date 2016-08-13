program Testedc

  implicit none
!
! This program tests the correctness of mdlaed3 and CauchyQ. 

  integer(kind = 4) :: i, N, lwork, info, ierr, liwork, LDQ, wflag
  double precision  :: ONE, ZERO, err, time, time1, time2
  parameter( ONE = 1.0E0, ZERO=0.0E0 )
 
  double precision, allocatable :: E(:), D(:), A(:,:), WORK(:), &
                    Dt(:), Q(:,:), DD(:), EE(:)
  integer, allocatable :: IWORK(:)
!
  DOUBLE PRECISION   DNRM2, comp_Wtime
  EXTERNAL           DNRM2, comp_Wtime

  N = 5000
  LDQ = N
  ! generate the diagonal and off-diagonal elements
  lwork = N*(N+100)
  liwork = N + 5*N
  ALLOCATE( D(N),E(N),DD(N),EE(N),IWORK(LIWORK), stat=ierr )
  IF( ierr .ne. 0 ) THEN
     WRITE(*,*) 'Allocate failed in TestHss.'
  END IF

  ! singualr vectors and workspaces
  ALLOCATE( A(N,N),Q(N,N),WORK(lwork),Dt(N), stat=ierr )
  IF( ierr .ne. 0) THEN
     WRITE(*,*) 'Wrong while allocating.'
  END IF

  DO I =1, N
     D(I) = 2.0D0
  END DO
  E = ONE
  DD = D
  EE = E

  write(*,*) 'Dim = ', N

  time= comp_Wtime()
  call DSTEDC( 'I',N,DD,EE,Q,LDQ,WORK,LWORK,IWORK,LIWORK,INFO )
  time1 = comp_Wtime()
  time1 = time1 - time
  write(*,*) 'Lapack costs', time1

  time = comp_Wtime()
  call MDSTEDC( 'I',N,D,E,Q,LDQ,WORK,LWORK,IWORK,LIWORK,INFO )
  time2 = comp_Wtime()
  time2 = time2 - time
  write(*,*) 'Mine costs', time2, 'Speedup is', time1/time2

!!$  A = Q
!!$  Dt = 0.0D0
!!$  call dgesvd( 'N','N',N,N,A,N,Dt,A,N,A,N,work,lwork,info )
!!$  err = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
!!$  write(*,*) "Othogonality of Q: ", err

  D = D-DD
  err = maxval( abs(D) )
  write(*,*) 'Max error of computed eigenvalues:', err

  DEALLOCATE( E,D,A,Q,WORK,IWORK,Dt )

end program Testedc

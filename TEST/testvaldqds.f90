program TestValdqds

  implicit none
!
! This routine tests the correctness and speedups of our structured BiDC algorithm
! Its computed singular values are compared with those by LAPACK dbdsdc routine and 
! dqds routine. 
!
! ============
!
! Modified by Shengguo Li, on Sep. 12th, 2013
! ===========================================
!

  integer(kind = 4) :: i, j, k, N, lwork, info, ierr, DIM(6), itt
  double precision  :: Rerr, time, time1, one, Merr
  parameter( one = 1.0E0 )
 
  integer, allocatable ::  IWORK(:)
  double precision, allocatable :: D(:), E(:), U(:,:), VT(:,:), WORK(:), &
                    Dt(:), DD0(:), EE0(:),Dq(:),Eq(:),WORKQ(:)

  DOUBLE PRECISION   DNRM2, comp_Wtime
  EXTERNAL           DNRM2, comp_Wtime

  DATA DIM /10000, 4000, 5000, 8000, 9000, 10000/

  DO ITT = 2, 4
  ! Dim of matrix
     N = DIM(ITT)
     WRITE(*,*) 'DIM = ', DIM( ITT )
  
  ! generate the diagonal and off-diagonal elements
     ALLOCATE( D(N),E(N-1),DD0(N),EE0(N-1), stat=ierr )
     IF( ierr .ne. 0 ) THEN
        WRITE(*,*) 'Allocate failed in testdc.'
     END IF
  
  ! singualr vectors and workspaces
     lwork = 3*N*N + 4*N
     ALLOCATE( U(N,N),VT(N,N),WORK(lwork),IWORK(8*N), stat=ierr )
     IF( ierr .ne. 0) THEN
        WRITE(*,*) 'Wrong while allocating Svecs.'
     END IF

  ! workspace for testing 
     ALLOCATE( Dt(N),Dq(N),Eq(N),WORKQ(4*N), stat=ierr )
     IF( ierr .ne. 0) THEN
        WRITE(*,*) 'Wrong while allocating Svecs.'
     END IF

! Example 1
  call random_number(DD0)
  call random_number(EE0)

! Example 2
!     DD0 = 2.00
!     EE0 = 1.0

! Example 3
!!$     EE0 = 1.0
!!$     DO i = 1, N
!!$        DD0(i) = N+1-i
!!$     END DO

! Our algorithm 
     D = DD0
     E = EE0
     time = comp_Wtime()
     call mdbdsdc( 'U','I',N,D,E,U,N,VT,N,WORK,IWORK,WORK,IWORK,info )
     time1 = comp_Wtime()
     time = time1-time
     write(*,*) 'Computing time of HSS is ', time

! dqds algorithm
     Dq = DD0 
     Eq = EE0
     time = comp_Wtime()
     call dlasq1( N, Dq, Eq, WORKQ, info )
     time1 = comp_Wtime()
     time = time1-time
     write(*,*) 'Computing time of dqds is ', time

! LAPACK 
     time = comp_Wtime()
     call dbdsdc( 'U','I',N,DD0,EE0,U,N,VT,N,WORK,IWORK,WORK,IWORK,info )
     time1 = comp_Wtime()
     time = time1-time
     write(*,*) 'Computing time of LAPACK is ', time

! Compare the error
     Dt = D - DD0
     Merr = MAXVAL( ABS (Dt) )
     WRITE(*,*) 'Compare with dbdsdc: Max error ', Merr
     
     ! compare with dqds
     Dt = D - Dq
     Merr = MAXVAL( ABS (Dt) )
     WRITE(*,*) 'Ours compares with dqds: Max error ', Merr

     ! dbdsdc compares with dqds
!!$     Dt = DD0 - Dq
!!$     Merr = MAXVAL( ABS (Dt) )
!!$     Dt = ABS( Dt/Dq )
!!$     Rerr = MAXVAL( Dt )
!!$     WRITE(*,*) 'DBDSDC compares with dqds: Relative error ', Rerr, 'Max error ', Merr
  
  DEALLOCATE( D,E,DD0,EE0,U,VT,WORK,IWORK,Dt,Dq,Eq,WORKQ )

END DO
  
end program TestValdqds

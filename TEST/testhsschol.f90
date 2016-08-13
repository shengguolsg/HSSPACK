       program testhsschol
!         
         implicit none

       real(8), allocatable :: A0(:,:), A(:,:), R(:), x(:), b(:), b0(:)
       integer       :: prank,pp,lwork,i,ltr,n,ni,nit,lm, Maxit, info
       double precision   :: tol, s0, s1, time, tau

       double precision, external :: dnrm2

       !open(1, file = 'B.txt')
       open(1,file = '../DATA/A.dat')
       read(1,*) n
       allocate( A0(n,n), A(n,n) )
       do i = 1, n
          read(1,*) A0(i,1:n)
       enddo
       A(1:n,1:n) = A0(1:n,1:n) ! kept A0 for error checking

!!!!! Parameters
       prank  = 10
       tol   = 1d-8   ! tolerance
       tau   = 1d-13  ! accuracy parameter for pcg 
       lwork = n*n    ! work space size; for test only; could be improved
       ni    = 20     ! block row size
       Maxit = 30
       ltr   = 0
       pp = 5
       info = 0

!!!!! Solution test 
       allocate( x(n), b(n), b0(1:n) )
       do i = 1, n
          x(i) = 1d0
       enddo

       call DGEMV( 'N',n,n,1d0,A0,n,x,1,0d0,b,1 )
       b0( 1:n ) = b( 1:n ) 
       
       x( 1:n ) = b( 1:n )
       
       call HssCholSolver( 'R','N',A,n,x,ni,tol,prank,pp,tau,MaxIt,info )

       s0 = dnrm2(n,b,1)
       call DGEMV('N',n,n,-1d0,A0,n,x,1,1d0,b,1)
       s1 = dnrm2(n,b,1)
       print '(a38,e11.4)','relative residual: ',s1/s0

!!!!! PCG test
       call HssCholSolver( 'R','Y',A0,n,x,ni,tol,prank,pp,tau,MaxIt,info )
!
       deallocate( x,b,b0,A0,A )
!
     end program testhsschol

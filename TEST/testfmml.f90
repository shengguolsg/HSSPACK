   program testfmml
!
     implicit none

       real(8), allocatable :: A0(:,:), A(:,:), x1(:,:), y1(:,:), z1(:,:)
       real(8), allocatable :: x(:), y(:),z(:), v(:)
       real(8)       :: tol,s0,s1,time(2),error,time2,time1,sa,sb,gap
       integer       :: prank, pp, typ, i, j, tnk, lvl, pnh,ncol
       integer       :: n, ni, ltr, lm, info, n1, ldx, ntt

       double precision :: dnrm2
       external         :: dnrm2

! **************   initial the matrix ********************

       ! Example 1
       n = 2000
       ldx = n
       ncol = 100
       sa = 0.1
       sb = 9.0
       gap = (sb-sa) / (2*n)
       allocate (A0(n,n),A(n,n),x1(ldx,n),y1(ldx,n),z1(ldx,n) ) 
       allocate( x(n),y(n),z(n),v(n) )
       call random_number(z)
       call random_number(y1)
       call random_number(x1)
       call random_number(v)
       DO i = 1, n
          x(i) = sa + gap*(2*i-1)
          y(i) = x(i)+gap
       END DO

       call dlasrt('I', n, x, info)
       call dlasrt('I', n, y, info)
       call Cauchylike(A,x,y,z,v,n,n)
       A0(1:n,1:n) = A(1:n,1:n)   ! kept A0 for error checking
       z1(1:n,1:n) = x1(1:n,1:n)
 
! ********************   Construct HSS tree  ********************
       prank = 45   
       tol   = 1.0d-11   ! tolerance      
       ni    = 64        ! block row size
       call cpu_time(time1)
       call DHSSMML( A,N,N,N,x1,N,N,tol,prank,ni,ncol,time,info )
       call cpu_time(time2)
       time2 = time2 - time1
       print*, 'Total time is', time2 
       

!  ************************ Solution Test *******************************
       call cpu_time(time1)
       call dgemm('N','N',ldx,n,n,1d0,z1,n,A0,ldx,0d0,y1,ldx) 
       call cpu_time(time2)
       time2 = time2 - time1
       print*, 'LA time is', time2 

    
       z1(1:n,1:n) = y1(1:n,1:n)-x1(1:n,1:n)
       !     error = dnrm2(n,z1,1)
       error = maxval( abs(z1) )
       write(*,*) '2-norm absolute error is ', error, ldx

       deallocate(A0,A,x,y,z,x1,y1,z1)
            
     end program testfmml

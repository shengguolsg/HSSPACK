program testcauchyhss_serial

       implicit none
! ===========
! This procedure compare the accuracy and speedup of multiplying an Cauchy matrix with a general matrix.
! We assume that D and F are interlacing, whose entries are known exactly. 
! To show the efficiency of our algorithm, the size of matrix can be very large. 
! 
! ==============
! Written by Shengguo Li, on Nov. 25th, Changsha, 2012
! ==============
       integer dim(10)
       double precision  one, zero
       parameter( one=1.0D0, zero=0D0 )

       real(8), allocatable :: A(:,:), B(:,:),C(:,:)
       real(8), allocatable :: D(:),F(:),U(:),V(:)  
       real(8)       :: tol, s1, error, time1, time2, time3, gap, sa, bb
       real(8)       :: time(2)
       integer       :: i, j, lvl, pnh,tnk,ncol,ierr,nc,info,N0
       integer       :: n, ni, ltr, lm, lwork, it, nlast

       double precision :: dnrm2
       external         :: dnrm2 

! **************   initial the matrix ********************
     data dim /1000, 6000, 3000, 4000, 8000, 10000, 1800, 3000, 5000, 8000/

     do it = 1, 1
       N = dim(it)
       N0 = N
       lwork = N*30
       write(*,*) 'Dim = ', N

       allocate ( A(n,n),D(n),F(n),U(N),V(N),stat=ierr )
       if(ierr /= 0) Then
          write(*,*) "Allocate Error1."
          stop
       end if
       allocate( B(n,n),C(n,n), stat=ierr )  !,Tau(n),work(lwork)
       if(ierr /= 0) Then
          write(*,*) "Allocate Error2."
          stop
       end if
       
       sa = 0.1D0
       bb = 9.0D0
       gap = (bb-sa)/ (2*N)
       DO I = 1, N
          F(i) = sa + (2*I-1)*gap
          D(i) = F(i) + gap
       END DO
       call random_number(U)
       call random_number(V)
       call Cauchylike(A,D,F,U,V,N,N)

! ***********************************************
!             Choosing different matrix B       *
! ***********************************************
!!$       ! random matrix
       call random_number(B)

       ! diagonal matrix
!!$       B = zero      
!!$       do j = 1,n
!!$          B(j,j) = one
!!$       end do

!!$       call random_number(B)
!!$       call dgeqrf(N,N,B,N,Tau,Work,lwork,info)  
!!$       call dorgqr(N,N,N,B,N,Tau,Work,lwork,info)
       
! ***********************************************
!         Time of using dgemm                   *
! ***********************************************
       call cpu_time(time1)
       call dgemm('N','N',n,n,n,1d0,A,n,B,n,0d0,C,n) 
       call cpu_time(time2)
       time2 = time2-time1
       print *, 'MMM time: ',time2

!   ******************** HSS Matrix Approximation *********************
       tol  = 1.0d-17   ! tolerance      
       ni   = 120        ! block row size
       ncol = 128

       call cpu_time(time1)
       call cauchyhssmm(D,F,U,V,n,ni,tol,ncol,B,N,time )
       call cpu_time(time3)
       time3 = time3 - time1 
       print*, 'Const time', time(1),'multi time: ', time(2), 'total time ', time2
       print*, 'Speedup is', time2/time3

!  ************************ Solution Test *******************************
        A = C - B
        A = abs(A)
        error = maxval(A)
        write(*,*) 'Max absolute error is ', error

        deallocate(A,B,C,D,F,U,V ) !,Tau,work)
     end do

   end program testcauchyhss_serial

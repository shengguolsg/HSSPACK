program testcauchyhss
    use aux_hss
    use BasicMM
    use CauchyHss_VP2
    use ConstructHssd

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

!      HSS tree
       integer, allocatable :: tr(:), m(:)
       type(HTinfo)          :: TRE
       
       real(8), allocatable :: A(:,:), WW(:,:),SX(:),H(:), DD(:), B(:,:),C(:,:)
       real(8), allocatable :: D(:),F(:),U(:),V(:)  !,Tau(:),work(:)
       real(8)       :: tol, s1, time, error, time1, gap, sa, bb
       integer       :: i, j, lvl, pnh,tnk,ncol,ierr,nc,info,N0
       integer       :: n, ni, ltr, lm, n1,lwork, it
       type(hssmm)   :: PH

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

!!$       ! orthogonal matrix
!!$       call random_number(B)
!!$       call dgeqrf(N,N,B,N,Tau,Work,lwork,info)  
!!$       call dorgqr(N,N,N,B,N,Tau,Work,lwork,info)
       
! ***********************************************
!         Time of using dgemm                   *
! ***********************************************
       time = comp_time()
       call dgemm('N','N',n,n,n,1d0,A,n,B,n,0d0,C,n) 
       time1 = comp_time()
       time = time1-time
       print *, 'MMM time: ',time

! ***********************************************
!          Time of using HSS                    *
! ***********************************************

! ********************   Construct HSS tree  ********************
       tol   = 1.0d-17   ! tolerance      
       ni    = 50        ! block row size
       ltr   = 0
       allocate (m(n/ni+1), tr( 2*(n/ni+1)-1 ))
       m = 0
       call npart(n, ni, tr, ltr, m, lm)

!     Tree information
      lvl = ceiling( log2(ltr+1) ) -1
      n1 = 0.5* (ltr +1 )
      
      allocate(TRE%ttr(lvl,n1), TRE%lentr(lvl), TRE%lfx(ltr), TRE%ch(2,ltr) )
      call GetTreeInfo(TR, LTR, M, M, TRE, lvl)       

!   ******************** HSS Matrix Approximation *********************
!      Allocate space and initialize variables
       tnk = 0
       DO i = 1, lm
          tnk = tnk+ m(i)**2
       END DO

       allocate( PH%D(ltr), PH%U(2,ltr), PH%V(2,ltr), PH%B(2,ltr),PH%ph(4,ltr), PH%pd(2,ltr) )
       allocate( H(8*ni*n), DD(tnk) ) ! test only; can be optimized
       call hssexpmm0(PH, ltr)

!       HSS matrix
       pnh = 1
       time = comp_time()
!       call cauchy2hss_OMP( D,F,U,V,N,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh,TRE,info )
       call cauchy2hss( D,F,U,V,N,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh )
       time1 = comp_time()
       time = time1-time
       s1 = (1.0D0*pnh )/(ni*n)
       print*, 'construction time: ',time, s1
        
!   **************** HSS Matrix-Matrix Multiplication *******************
        ncol = 200
!        Allocate(WW(2*n,ncol), SX(2*lvl*n*ncol))

        nc = N0/ncol
        N = N0
        time = comp_time()
!$OMP PARALLEL DO
        Do j = 1, nc
           call fastHssmm_omp( DD,H,PH,m,B(1:N,(j-1)*ncol+1:j*ncol),tr,ltr,n,ncol,TRE,lvl )
        End Do
!$OMP END PARALLEL DO
        time1 = comp_time()
        time = time1 - time
        print*, 'multiplication time: ', time, 'ncol ', ncol

!  ************************ Solution Test *******************************
        A = C - B
        A = abs(A)
        error = maxval(A)
        write(*,*) 'Max absolute error is ', error

        deallocate(TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch, m, tr)
        deallocate(PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd)
        deallocate(A,B,C,D,F,U,V,H,DD ) ! ,WW,SX ) !,Tau,work)
     end do

   end program testcauchyhss

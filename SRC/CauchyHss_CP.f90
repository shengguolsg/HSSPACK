module CauchyHss_CP
  implicit none

contains

!!!!!!
    subroutine rrluCauchy(D,F,A,V,tol,U,W,PL,PU,M,N,Rk)
!
! .. Scalar parameter ..
    DOUBLE PRECISION  TOL
    INTEGER           M,N,Rk
! .. Array parameters ..
    DOUBLE PRECISION D(*),F(*),A(*),V(*),U(*),W(*)
    INTEGER          PL(*),PU(*)
!
! Purpose
! =========
! This routine computes a Rank Revealing Schur Complement, RRSC, for a Cauchy-like matrix
! via structured matrices. This matrix C has dimensions M-by-N, with generators D,F,A and V.
! C = (A(i)*V(j) / (D(i)-F(j)) )_{i,j}. 
!
! .. Parameters ..
! D  (INPUT/OUTPUT)  DOULBE PRECISION Array, DIMENSION(M). Row generators.
!    It will be modified and the first Rk elements are chosen ones. 
! 
! F (INPUT)  DOUBLE PRECISION Array, DIMENSION(N). Col generators.
!
! A (INPUT)  DOUBLE PRECISION Array, DIMENSION(M). Row generators
!   
! B (INPUT)  DOUBLE PRECISION Array, DIMENSION(M). Col generators
! 
! tol (INPUT) DOUBLE PRECISION Array
! 
! U (OUTPUT) DOUBLE PRECISION Array, DIMENSION(M). Row generators of Mat_N. 
!   Mat_N is a matrix in RRSC. The last M-Rk elements of U are useful.
!
! V (OUTPUT) DOUBLE PRECISION Array, DIMENSION( min(M,N) ). Row generators of Mat_N. 
!   Mat_N is a matrix in RRSC. The first Rk elements of V are useful.
!
! PL (OUTPUT) INTEGER Array, DIMENSION(M). Row permutation.
!
! PU (OUTPUT) INTEGER Array, DIMENSION(M). Col permutation.
!    Not used in our current problem.
! 
! M (INPUT) INTEGER, Row dimension of C
!
! N (INPUT) INTEGER, Col dimension of C
! 
! Rk (OUTPUT) INTEGER, rank of returned matrix.
!
! ============
! Written by Sheng-Guo Li, On Sept. 18th, 2012
! ============

! .. Local Scalars
    double precision rho, Amax, Ukk, zero, one, Nmax
    integer          nswap,j,k,mn,prd,flgl,ii,jj
    parameter        (zero = 0.0D0, one=1.0D0)

!  .. Intrinsic Functions ..
    intrinsic    max, abs, maxval,mod

    U(1:M)  = A(1:M)
    Rk = min(M,N)
    mn = Rk
    PL(1:M) = (/ (j, j=1,M) /)
    PU(1:N) = (/ (j, j=1,N) /)
    nswap = 0
    Amax = zero
    rho = 1.02D0
    prd = 10
    
    do k = 1, mn
       call CauchyPivtG_CP(D,F,U,V,PL,PU,k,A,M,N)

       Ukk = u(k)*V(k)/( D(k)-F(k))
       Amax = max(Amax,abs(Ukk) )
       if (abs(Ukk) .lt. tol*Amax ) then  ! converged
          Rk = k -1
!          call searchMax2(U(Rk+1),W,D(Rk+1),D,M-k,k,ii,jj,Nmax) 
!          if(Nmax .gt. one) then
!             write(*,*) 'Max N value is ', Nmax
!          end if

          exit
       end if
       
       U(k+1:M) = U(k+1:M)* (D(k+1:M)-D(k)) / (D(k+1:M)-F(k))
       V(k+1:N) = V(k+1:N)* (F(k+1:N)-F(k)) / (F(k+1:N)-D(k))
       do j = 1, k-1
          W(j) = W(j) * ( (F(k)-D(j))/ (D(k)-D(j)) )
       end do
       W(k) = (D(k)-F(k))/A(k)
       do j = 1, k-1
          W(k) = W(k) * (F(j)-D(k)) / ( D(j)-D(k) )
       end do

       ! swap
       flgl = mod(k,prd)
       do while(flgl .lt. 1)
          flgl = 1
          call searchMax2(U(k+1),W,D(k+1),D,M-k,k,ii,jj,Nmax)
          
          if(Nmax .gt. rho) then
             nswap = nswap + 1
             flgl = 0
             jj = jj + k 
             V(k+1:N)    = V(k+1:N) * ( (F(k+1:N)-D(ii)) / (F(k+1:N)-D(jj)) )
             U(k+1:jj-1) = U(k+1:jj-1) * ( (D(k+1:jj-1)-D(jj)) / (D(k+1:jj-1)-D(ii)) )
             U(jj+1:M)   = U(jj+1:M) * ( (D(jj+1:M)-D(jj)) / (D(jj+1:M)-D(ii)) )
             W(1:ii-1)   = W(1:ii-1) * ( (D(1:ii-1)-D(ii)) / (D(1:ii-1)-D(jj)) )
             W(ii+1:k)   = W(ii+1:k) * ( (D(ii+1:k)-D(ii)) / (D(ii+1:k)-D(jj)) )
             U(jj)       = A(ii) * ( (D(ii)-D(jj)) / (D(ii)-F(ii)) )  
             W(ii)       = (D(jj)-F(ii)) / A(jj)
             do j = 1, ii-1
                U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
             end do
             do j = ii+1, k
                U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
             end do
             call iswap(PL,ii,jj)
             call dswap(D,ii,jj)
             call dswap(A,ii,jj)
!             write(*,*) 'swap once '
          end if ! Nmax
          
       end do ! while

    end do ! main loop

  end subroutine rrluCauchy

!!!!!!
    subroutine searchMax2(U,W,D2,D1,LDU,LDW,ii,jj,Nmax)
! choose the maximum column by column

    double precision U(*),W(*),D2(*),D1(*)
    integer LDU,LDW,ii,jj

    double precision zero,junk,Nmax
    integer j,jjL
    parameter(zero = 0.0D0)
    
    call CauchyMax(D2,D1(1),u,LDU,junk,jjL)
    junk = junk*abs(w(1))       
    Nmax = junk
    ii = 1
    jj = jjL
    
    do j = 2, LDW   ! ii: col, jj: row
       call CauchyMax(D2,D1(j),u,LDU,junk,jjL)
       junk = junk*abs(w(j))       
       if(junk .gt. Nmax) then
          Nmax = junk
          ii = j
          jj = jjL
       end if
    end do

  end subroutine searchMax2

!!!!!!
  subroutine CauchyMax(D,F,u,N,junk,jj)
! returen the largest entry, junk, of a column of Cauchy matrix and its index, jj.

    double precision D(*),u(*),junk,F
    integer N,jj
    integer temp(1)

!  .. Intrinsic Functions ..
    intrinsic    Maxloc, ABS,MAXVAL

    double precision, allocatable :: LK(:) ! change it to LK(N)

    allocate( LK(N) )

    Lk = u(1:N)/ ( D(1:N)-F )
    Lk = abs(Lk)
    junk = maxval(Lk(1:N))
    temp = maxloc(Lk(1:N))
    jj = temp(1)
    
    deallocate(LK)
    
  end subroutine CauchyMax

!!!!!!
  subroutine CauchyPivtG(D,F,u,v,PL,PU,k,A,M,N)
!
    integer k,M,N
    double precision D(*),F(*),u(*),v(*),A(*)
    integer PL(*), PU(*)

!  .. Intrinsic Functions ..
    intrinsic    ABS

!   local parameters 
    integer jjL,jjU,flg,jj
    double precision junkL,junkU, Piv, pivot,zero,one    
    parameter    (zero = 0.0D0, one=1.0D0)

    call CauchyMax(D(k),F(k),u(k),M-k+1,junkL,jjL)
    call CauchyMax(F(k),D(k),v(k),N-k+1,junkU,jjU)
    junkL = junkL * abs(v(k))
    junkU = junkU * abs(u(k))
    Piv   = abs(u(k)*v(k)/ (D(k)-F(k)) )
    
    if (junkL .le. Piv .and. junkU .le. Piv ) then
       return
    end if
    
    pivot = zero
    flg = 0
    if(junkL > junkU) flg = 1
    
    do while (1 < 2)
       pivot = pivot +1
       if (flg == 1) then
          jj = jjL
          call dswap(D,k,jj+k-1)
          call dswap(u,k,jj+k-1)
          call dswap(A,k,jj+k-1)
          call iswap(PL,k,jj+k-1)
          call CauchyMax( F(k),D(k),v(k),N-k+1,junkU,jjU ) ! N-k+1
          if(jjU == 1) return
          
          flg = 0
          continue
       end if
       jj = jjU
       call dswap(F,k,jj+k-1)
       call dswap(v,k,jj+k-1)
       call iswap(PU,k,jj+k-1)
       call CauchyMax( D(k),F(k),u(k),M-k+1,junkL,jjL )
       if(jjL == 1) return 

       flg = 1
    end do

  end subroutine CauchyPivtG

!!!!!!!!
    subroutine CauchyPivtG_CP( D,F,u,v,PL,PU,k,A,M,N )
!
!   This routine chooses the largest entry for the k-th Schur Complement.
!   It uses complete pivoting. If this routine could be faster than Matrix-Matrix
!   multiplication, we can try to modify the codes for bidiagonal SVD too. 
! 
    integer k,M,N
    double precision D(*),F(*),u(*),v(*),A(*)
    integer PL(*), PU(*)

!  .. Intrinsic Functions ..
    intrinsic    ABS
! ..
! .. local parameters 
    integer ii,jj,j
    double precision junk,zero,one
    parameter    (zero = 0.0D0, one=1.0D0)
!   ii : the row index; jj : the column index;
! ..
! .. local array ..
    integer iiU(N-k+1)    ! the row index for each column
    double precision junkU(N-k+1)  ! the largest entry for each column
!
!$OMP PARALLEL DO    
    DO j =1, N-k+1
       call CauchyMax(D(k),F(k+j-1),u(k),M-k+1,junkU(j),iiU(j) )
       junkU(j) = abs( junkU(j)*V(k+j-1) ) 
    END DO
!$OMP END PARALLEL DO

    junk = junkU(1)
    ii = iiU(1)
    jj = 1
    DO j = 2, N-k+1
       IF(junkU(j) > junk ) THEN
          junk = junkU(j)
          ii = iiU(j)
          jj = j
       END IF
    END DO
        
    call dswap(D,k,ii+k-1)
    call dswap(u,k,ii+k-1)
    call dswap(A,k,ii+k-1)
    call iswap(PL,k,ii+k-1)
    
    call dswap(F,k,jj+k-1)
    call dswap(v,k,jj+k-1)
    call iswap(PU,k,jj+k-1)
    
  end subroutine CauchyPivtG_CP

!!!!!!
  subroutine dswap(D,mk,nk)
! double 1D array

    double precision D(*), tt
    integer mk,nk
    tt = D(mk)
    D(mk) = D(nk)
    D(nk) = tt
    
  end subroutine dswap

  subroutine iswap(D,mk,nk)
! integer 1D array

    integer D(*), tt, mk, nk
    tt = D(mk)
    D(mk) = D(nk)
    D(nk) = tt
    
  end subroutine iswap

!!!!!!
  subroutine cauchy2hss(D,F,U,V,LDA,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh)
       use aux_hss
!
!  .. Scalar Arguments ..
       INTEGER          :: LDA, LTR, lm, lvl
       DOUBLE PRECISION :: TOL
!
!  .. Array Arguments ..
       DOUBLE PRECISION :: D(*),F(*),U(*),V(*),H(*), DD(*)
       INTEGER          :: TR(*), M(*)
       TYPE(HSSMM)      :: PH
! 
!  Purpose
! =========
! Construct an HSS matrix for a Cauchy-like matrix A. The generators of A are D, F,U and V. 
! 
! Parameters
!========== 
! D,F,U,V (INPUT)  The generators of a Cauchy-like matrix, A_{ij}= U(i) V(j)/ (D(i)-F(j)).
!        They are all n-dimensional vectors. 
!
! LDA   (INPUT) INTEGER
! 
! TR    (INPUT) INTEGER array, DIMENSION(ltr), The HSS tree in post ordering.
!       
! LTR   (INPUT) INTEGER, length of HSS Tree
!
!  M    (INPUT) INTEGER array, DIMENSION(lm), the block partion of rows and columns
! 
!  LM   (INPUT) INTEGER, number of leaf nodes
!
!  PH   (OUTPUT) HSSMM TYPE, contains the information of generators in H and D
!         1) starting position of each generator, D, U, V, B
!         2) dimensions of each generator
!       Note that the generators are U, V^T, D and B. We get V^T not V. 
!
!  D    (OUTPUT) DOUBLE PRECISION array, DIMENSION(sum(m.^2)), the diagonal blocks
! 
!  H    (OUTPUT) DOUBLE PRECISION array, DIMENSION( LDA *K * alpha )
!         alpha is a constant; K is the HSS rank; 
! 
!  TOL  (INPUT) DOUBLE PRECISION, tolerance of low-rank approximation
!
!===================
!  Written by Shengguo Li, on Sept. 21st, 2012
!
!  Construct HSS matrix from Cauchy-like matrix by using Structured RRLU factorization
!  
!===================

!  .. Local Scalars ..
       INTEGER       :: ch1,ierr,info,nn,mi,ni,nt,Rk,it
       INTEGER       :: pnd,pnh,i,j,lt,ns
       DOUBLE PRECISION :: zero, one
       PARAMETER        (zero = 0.0D0, one=1.0D0)
!  .. Local Arrays ..
       integer, allocatable :: ch(:,:), l(:,:)
       integer, allocatable :: lsr(:), lsc(:),indx1(:),indx2(:)
       double precision, allocatable :: Di(:),Fi(:),Ui(:),Vi(:)
       ! lsr is the pointer of block row in stack A; lsc is the pointer of block column in stack A
!  ..
!  .. External Subroutines ..
      EXTERNAL  DGEMM, DLACPY
!  ..
!  .. External Functions .. 
      
!  .. Intrinsic Functions ..
      INTRINSIC    MIN, ABS,MAXVAL
              
       pnh = 1               ! 'pointer' in H
       pnd = 1               ! 'pointer' in D
       allocate( ch(2,ltr),l(ltr,2),lsc(lvl+2),lsr(lvl+2), stat=ierr )
       IF(ierr /= 0 ) Then
          write(*,*) "Allocate failed in cauchy2hss! "
          stop
       END IF

       allocate(indx1(lda),indx2(lda),Di(lda),Fi(lda),Ui(lda),Vi(lda))

       nn = ltr
       call child(tr, ch, ltr)

       l(1,1:2) = (/1,m(1)/)
       lt = 1
       it = 1
       do i = 1, ltr
          if (ch(1,i) == 0 .and. ch(2,i) == 0) then
             l(i,1:2) = (/lt,lt+m(it)-1/)
             lt = l(i,2)+1
             it = it+1
          else
             l(i,1:2) = (/l(ch(1,i),1), l(ch(2,i),2)/)
          endif
       enddo

       ns = 1                  ! (# of blocks in stack)+1
       lsr(ns) = 1             ! row start position at upper triangular
       lsc(ns) = 1             ! column start position at lower triangular
       
! *****************************************************
!                 MAIN LOOP                           *
!******************************************************
       do i = 1, nn
 
         ! leaf node
          if (ch(1,i) == 0 ) then 
             mi = l(i,2)-l(i,1)+1
             ni = lda-l(i,2)    !t^r, the right column block's length
             indx1(1:mi) = (/ (j, j=l(i,1),l(i,2) ) /)
             
             call Cauchylike(DD(pnd),D(l(i,1)),F(l(i,1)),U(l(i,1)),V(l(i,1)),mi,mi)
             PH%D(i) = mi
             PH%pd(1,i) = pnd
             PH%pd(2,i) = pnd + mi*mi -1
             pnd = pnd + mi*mi
             
             ! off-diag row compression
             if(ns .eq. 1 ) then
                indx2(1:ni) = (/ (j, j=l(i,2)+1,lda) /) ! col
                nt = ni
             else
                indx2(1:lsc(ns)-1) = (/ (j, j=1,lsc(ns)-1) /)
                indx2(lsc(ns):lsc(ns)+ni-1) = (/ (j, j=l(i,2)+1,lda) /)
                nt = ni + lsc(ns)-1
             end if
             
             Di(1:mi) = D(indx1(1:mi))
             Ui(1:mi) = U(indx1(1:mi))
             Fi(1:nt) = F(indx2(1:nt))
             Vi(1:nt) = V(indx2(1:nt))
             call comprcauchy('r',Di,Fi,Ui,Vi,tol,mi,nt,PH,pnh,H,i,Rk)

             lsr(ns+1) = lsr(ns)+Rk
             D(lsr(ns):lsr(ns+1)-1) = Di(1:Rk)
             U(lsr(ns):lsr(ns+1)-1) = Ui(1:Rk)

             ! off-diag col compression
             if(ns .eq. 1 ) then
                indx2(1:ni) = (/ (j, j=l(i,2)+1,lda) /) ! col
                nt = ni
             else
                indx2(1:lsr(ns)-1) = (/ (j, j=1,lsr(ns)-1) /)
                indx2(lsr(ns):lsr(ns)+ni-1) = (/ (j, j=l(i,2)+1,lda) /)
                nt = ni + lsr(ns)-1
             end if

             Di(1:nt) = D(indx2(1:nt))
             Ui(1:nt) = U(indx2(1:nt))
             Fi(1:mi) = F(indx1(1:mi))
             Vi(1:mi) = V(indx1(1:mi))
             call comprcauchy('c',Di,Fi,Ui,Vi,tol,mi,nt,PH,pnh,H,i,Rk)

             lsc(ns+1) = lsc(ns)+Rk
             F(lsc(ns):lsc(ns+1)-1) = Fi(1:Rk)
             V(lsc(ns):lsc(ns+1)-1) = Vi(1:Rk)
             ns = ns+1

          else ! parent nodes
             ch1 = ch(1,i)
             mi  = lsr(ns-1) - lsr(ns-2)
             ni  = lsc(ns)   - lsc(ns-1)
             call Cauchylike(H(pnh),D(lsr(ns-2)),F(lsc(ns-1)),U(lsr(ns-2)),V(lsc(ns-1)),mi,ni)
             call hssexpmm2('B',PH,pnh,mi,ni,ch1,info)  ! B{ch{tr{i}}(1)}

             ch1 = ch(2,i)
             mi  = lsr(ns) - lsr(ns-1)
             ni  = lsc(ns-1) - lsc(ns-2)
             call Cauchylike(H(pnh),D(lsr(ns-1)),F(lsc(ns-2)),U(lsr(ns-1)),V(lsc(ns-2)),mi,ni)
             call hssexpmm2('B',PH,pnh,mi,ni,ch1,info)  ! B{ch{tr{i}}(2)}

             if(i .eq. nn) exit
             
             ! off-diag row compression
             mi = lsr(ns) - lsr(ns-2)
             ni = lda-l(i,2)     !t^r
             indx1(1:mi) = (/ (j, j=lsr(ns-2),lsr(ns)-1) /)
             if(ns .eq. 3) then
                indx2(1:ni) = (/ (j, j=l(i,2)+1,lda) /) ! col
                nt = ni
             else
                indx2(1:lsc(ns-2)-1) = (/ (j, j=1,lsc(ns-2)-1) /)
                indx2(lsc(ns-2):lsc(ns-2)+ni-1) = (/ (j, j=l(i,2)+1,lda) /)
                nt = ni + lsc(ns-2)-1
             end if
             
             Di(1:mi) = D(indx1(1:mi))
             Ui(1:mi) = U(indx1(1:mi))
             Fi(1:nt) = F(indx2(1:nt))
             Vi(1:nt) = V(indx2(1:nt))
             call comprcauchy('r',Di,Fi,Ui,Vi,tol,mi,nt,PH,pnh,H,i,Rk)

             lsr(ns-1) = lsr(ns-2) + Rk
             D(lsr(ns-2):lsr(ns-1)-1) = Di(1:Rk)
             U(lsr(ns-2):lsr(ns-1)-1) = Ui(1:Rk)

             ! off-diag col compression
             mi = lsc(ns) - lsc(ns-2)
             ni = lda-l(i,2)    !t^r
             indx1(1:mi) = (/ (j, j=lsc(ns-2),lsc(ns)-1) /)
             if(ns .eq. 3 ) then
                indx2(1:ni) = (/ (j, j=l(i,2)+1,lda) /) ! col
                nt = ni
             else
                indx2(1:lsr(ns-2)-1) = (/ (j, j=1,lsr(ns-2)-1) /)
                indx2(lsr(ns-2):lsr(ns-2)+ni-1) = (/ (j, j=l(i,2)+1,lda) /)
                nt = ni + lsr(ns-2)-1
             end if

             Di(1:nt) = D(indx2(1:nt))
             Ui(1:nt) = U(indx2(1:nt))
             Fi(1:mi) = F(indx1(1:mi))
             Vi(1:mi) = V(indx1(1:mi))
             call comprcauchy('c',Di,Fi,Ui,Vi,tol,mi,nt,PH,pnh,H,i,Rk)

             lsc(ns-1) = lsc(ns-2)+Rk
             F(lsc(ns-2):lsc(ns-1)-1) = Fi(1:Rk)
             V(lsc(ns-2):lsc(ns-1)-1) = Vi(1:Rk)
             ns = ns-1

          end if
          
       end do ! main loop

       deallocate( ch,l,lsc,lsr,indx1,indx2,Di,Fi,Ui,Vi )

  end subroutine cauchy2hss

!!!!!!
  subroutine comprcauchy(rowcol,D,F,U,V,tol,M,N,PH,pnh,H,nodi,Rk)
    use aux_hss
!
! Scalar parameters
    integer M,N,pnh,nodi,Rk
    double precision tol
    character(len=1) rowcol
! Array parameters
    double precision D(*),F(*),U(*),V(*),H(*)
    TYPE(HSSMM)   :: PH
!
! Purpose
! ========
! This routine computes a low rank approximation of some off-diagonal blocks of a
! Cauchy-like matrix. 
! 
! D (input/output) double array, dimension(M), the first Rk entries is useful
!
! F (input)  double array, dimension(N)
! 
! U (input/output) double array, dimension(M), the first Rk entries is useful
! 
! V (input)  double array, dimension(N)
!
! tol (input) double precision
!
! M (input) integer, row dimension of A
!
! N (input) integer, col dimension of A
!
! PH (input/output) HSSMM type, stores the information of U,V and B
!
! pnh (input/output) the end of array H
!
! H (input/output) Double precision array, 1D, stores U,V and B.
!
! nodi (input) node i
!
! =========
! Written by S.G. Li, on Sept. 18th, 2012
! =========

! local scalars 
    integer mn,i,info
    logical  CR
    
! local arrays
    integer PL(M), PU(N)
    double precision zero, one
    parameter (zero = 0.0D0, one = 1.0D0)
    double precision, allocatable :: Q(:,:),Z(:),W(:)

!  .. Intrinsic Functions ..
    intrinsic    Min
    
!  .. External Functions
      LOGICAL     LSAME
      EXTERNAL    LSAME
      
      mn = min(M,N)
      allocate(Z(M),W(mn))
      Z = zero
      W = zero
      CR = lsame(rowcol,'r')
      if( CR )  then
         call rrluCauchy(D,F,U,V,tol,Z,W,PL,PU,M,N,Rk)
         call invp(PL,M)       ! PL --> InvPL
         allocate( Q(M,Rk) )
         if (Rk .lt. M) then            
            call Cauchylike2(Q,D(Rk+1),D,Z(Rk+1),W,M,Rk)    
            Q(1:M,1:Rk) = Q(PL,1:Rk)
            call dlacpy('A',M,Rk,Q,M,H(pnh),M)        ! copy U to H
            call hssexpmm2('U',PH,pnh,M,Rk,nodi,info) ! copy Q to generators            
         else
            ! copy identity matrix to generators
            Q(1:M,1:M) = zero
            do i = 1,M
               Q(i,i) = one
            end do
            Q(1:M,1:M) = Q(PL,1:M)
            call dlacpy('A',M,Rk,Q,M,H(pnh),M)        ! copy U to H
            call hssexpmm2('U',PH,pnh,M,Rk,nodi,info) ! copy Q to generators
         end if

      else ! block col
         
         call rrluCauchy(F,D,V,U,tol,Z,W,PL,PU,M,N,Rk)
         call invp(PL,M)  ! PL --> InvPL
         allocate( Q(Rk,M) )
         if (Rk .lt. M) then            
            call Cauchylike2(Q,F(Rk+1),F,Z(Rk+1),W,Rk,M)             
            Q(1:Rk,1:M) = Q(1:Rk,PL)
            call dlacpy('A',Rk,M,Q,Rk,H(pnh),Rk)      ! copy V to H
            call hssexpmm2('V',PH,pnh,Rk,M,nodi,info) ! copy Q to generators            
         else
            ! copy identity matrix to generators
            Q(1:M,1:M) = zero
            do i = 1,M
               Q(i,i) = one
            end do
            Q(1:M,1:M) = Q(1:M,PL)            
            call dlacpy('A',M,M,Q,M,H(pnh),M)         ! copy V to H
            call hssexpmm2('V',PH,pnh,Rk,M,nodi,info) ! copy Q to generators
         end if
         
      end if ! compr type

      deallocate(Z,W,Q)

    end subroutine comprcauchy

!!!!!!
    subroutine invp(P,N)
      integer N,P(*)
      integer i
      integer :: IP(N)

      IP = 0

      do i = 1, N
         IP( P(i) ) = i
      end do
      P(1:N) = IP(1:N)

    end subroutine invp

!!!!!!
  subroutine Cauchylike(A,D,F,u,v,M,N)
!   A(i,j) = u(i)v(j) / (D(i)-F(j))

    integer M,N,i,j
    double precision A(M,*), D(*),F(*),u(*),v(*)

    do j = 1,N
       do i = 1,M
          A(i,j) = u(i)*v(j) / (D(i)-F(j))
       end do
    end do
    
  end subroutine Cauchylike

!!!!!!
  subroutine Cauchylike2(A,D,F,u,v,M,N)
! A(i,j) = u(i)v(j) / (D(i)-F(j))

    integer M,N,i,j,Rk
    double precision A(M,*), D(*),F(*),u(*),v(*)
    double precision zero, one
    parameter(zero = 0.0D0, one = 1.0D0)
    
    Rk = M-N
    A(1:M,1:N) = zero
    
    if(Rk .eq. 0) then
       do j = 1, M
          A(j,j) = one
       end do
    else
       if(Rk .gt. 0) then
          do j = 1,N
             A(j,j) = one
             do i = N+1,M
                A(i,j) = u(i-N)*v(j) / (D(i-N)-F(j))
             end do
          end do
       else
          do j = 1, M
             A(j,j) = one
          end do
          do j = M+1, N
             do i = 1, M
                A(i,j) = v(i)*u(j-M) / (D(j-M)-F(i))
             end do
          end do
       end if
    end if
    
  end subroutine Cauchylike2


end module CauchyHss_CP

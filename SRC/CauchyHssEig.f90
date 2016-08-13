module CauchyHssEig
  implicit none

contains

!!!!!!!!
  SUBROUTINE rrluCauchyEigrow( F,D,A,V,DIFL,DIFR,tol,U,W,FF,PL,PU,M,N,Rk,PQ,nflops,nswap )
!
! ..Scalar parameter..
    INTEGER, INTENT(IN)          :: M,N
    INTEGER, INTENT(OUT)         :: Rk,nflops,nswap
    DOUBLE PRECISION, INTENT(IN) :: tol
! ..
! ..Array parameters..
    DOUBLE PRECISION, INTENT(IN)    :: DIFL(*),FF(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: A(*),V(*),U(*),W(*),D(*),F(*)
    INTEGER, INTENT(INOUT)          :: PL(*),PU(*),PQ(*)
!
! Purpose
! =========
! This routine computes a Rank Revealing Schur Complement, RRSC, for a Cauchy-like matrix
! via structured matrices. This matrix C has dimensions M-by-N, with generators D,F,A and V.
! C(i,j) = A(i)*V(j) / (F(i) -D(j) ) satisfies F*C - C*D = A*V.  
! Here F are the old eigenvalues and D are the updated ones. 
! F is the row generator, and D is the column generator. 
!
! It computes a low rank approximation of C, C(PQ,:) = [I; Mat_N]*C( PQ(1:Rk), :) where Z is also a 
! Cauchy-like matrix and satisfies F(Rk+1:M) * Mat_N - Mat_N * F(1:Rk) = U * W. 
! This routine is used to compress the HSS block row of eigenvector matrix.
!
! Parameters
! ==========
! F  (inout)  DOULBE PRECISION Array, DIMENSION(M). Row generators, the old eigenvalues.
!    It will be modified and the first Rk elements are the chosen ones. 
! 
! D  (inout)  DOUBLE PRECISION array, DIMENSION(N). Col generators, the updated eigenvalues.
!
! A (inout)  DOUBLE PRECISION array, DIMENSION(M). Row generators, updated z.
!   
! V (inout)  DOUBLE PRECISION array, DIMENSION(N). Col generators, alpha.
! 
! tol (in) DOUBLE PRECISION array
!
! DIFL (INPUT) DOUBLE PRECISION array, DIMENSION(LN) 
!      the distance between FF and D0, positive values, DIFL(i) = D0(i)-FF(i)
!
! DIFR (INPUT) DOUBLE PRECISION array, DIMENSION(LN) 
!      the distance between FF and D0, negative values, DIFL(i) = D0(i)-FF(i+1)
!
! U (OUTPUT) DOUBLE PRECISION array, DIMENSION(M). 
!   Row generators of Mat_N. Mat_N is a matrix in RRSC. 
!   The last M-Rk elements of U are useful. Initially, U = A. 
!
! W (OUTPUT) DOUBLE PRECISION array, DIMENSION( min(M,N) ). 
!   Row generators of Mat_N. Mat_N is a matrix in RRSC. 
!   The first Rk elements of W are useful.
!
! FF (input) double precision array, dimension( LN )
!    The orignal old eigenvalues.
!
! PL (OUTPUT) INTEGER array, DIMENSION(M). Row permutation. 
!    It relates with row generators, F and FF
!
! PU (OUTPUT) INTEGER array, DIMENSION(N). 
!    Col permutation. It relates with column generatros, D and D0
! 
! M (INPUT) INTEGER, Row dimension of C
!
! N (INPUT) INTEGER, Col dimension of C
! 
! Rk (OUTPUT) INTEGER, rank of returned matrix.
!
! PQ (INOUT)  INTEGER Array, DIMENSION( M )
!    Records the row permutations.
!
! ============
! Written by Sheng-Guo Li, On April 16th, 2013
! for the tridiagonal eigenvalue problem 
! ===========================================
!
! .. Local Parameters ..
    DOUBLE PRECISION  ONE, ZERO, NEGONE
    PARAMETER        (ZERO= 0.0E+0, ONE=1.0E+0, NEGONE= -1.0E+0 )
! ..
! .. Local Scalars ..
    INTEGER          j,k,mn,prd,flgl,ii,jj,lflops
    DOUBLE PRECISION rho,Amax,Ukk,Nmax,temp,temp2
! ..
! .. Intrinsic Functions ..
    INTRINSIC  MAX, MIN, ABS, MAXVAL, MOD

    U(1:M)  = A(1:M)
    PQ(1:M) = (/ (j, j=1,M) /)
    Rk = MIN(M,N)
    mn = Rk
    nswap = 0
    Amax = zero
    rho = 1.5E+0
    prd = 10
    nflops = 0
    
    DO k = 1, mn
       call CauchyPivtEig( F,D,U,V,PL,PU,k,A,M,N,PQ,lflops,FF,DIFL,DIFR )
       ! F is old; D is new; U is Z and V is alpha; PL is for F and PU is for D
       nflops = nflops + lflops

       temp = computfd( FF,k,k,PL,PU,DIFL,DIFR )
       Ukk  = U(k)*V(k) / temp
       Amax = MAX( Amax,ABS(Ukk) )
       IF( ABS(Ukk) .LT. tol*Amax ) THEN  ! first step converged
          ! Complete pivoting for Schur complement
          call CauchyPivtEig_CP( F,D,U,V,PL,PU,k,A,M,N,PQ,lflops,FF,DIFL,DIFR )
          temp = computfd( FF,k,k,PL,PU,DIFL,DIFR )
          Ukk  = U(k)*V(k) / temp
          Amax = MAX( Amax,ABS(Ukk) )
          IF( ABS(Ukk) .LT. tol*Amax ) THEN  ! final converged
             Rk = k-1
             EXIT
          END IF
       END IF
       
       DO j = k+1, M    !  U(k+1:M) = U(k+1:M)* (F(k+1:M)-F(k)) / (F(k+1:M)-D(k))
          temp = computfd( FF,j,k,PL,PU,DIFL,DIFR )  ! F(j) - D(k)
          Temp = (F(j)-F(k)) / Temp
          U(j) = U(j) * Temp
       END DO

       DO j = k+1, N   !  V(k+1:N) = V(k+1:N)* (D(k+1:N)-D(k)) / (D(k+1:N)-F(k))
          temp  = computdd( FF,j,k,PU,DIFL,DIFR )            ! D(j) - D(k)
          temp2 = NEGONE* computfd( FF,k,j,PL,PU,DIFL,DIFR ) ! D(j) - F(k)
          temp  = temp/ temp2
          V(j) = V(j) * Temp
       END DO
              
       DO j = 1, k-1    ! W(j) = W(j) * ( (D(k)-F(j))/ (F(k)-F(j)) )
          temp = NEGONE* computfd( FF,j,k,PL,PU,DIFL,DIFR )  ! D(k)-F(j)
          W(j) = W(j) *( Temp / (F(k)-F(j)) )
       END DO
       temp = computfd( FF,k,k,PL,PU,DIFL,DIFR )  ! F(k)-D(k)
       W(k) = temp / A(k)

       DO j = 1, k-1   ! W(k) = W(k) * (D(j)-F(k)) / ( F(j)-F(k) )
          temp = NEGONE* computfd( FF,k,j,PL,PU,DIFL,DIFR )   ! D(j)-F(k)
          W(k) = W(k) * ( Temp/ (F(j)-F(k) ) )
       END DO

       nflops = nflops + 5*(M-k)+7*(N-k)+2*6*(k-1)

       ! swap
       flgl = MOD(k,prd)
       DO WHILE(flgl .lt. 1)
          flgl = 1
          call searchMax2Eig(U(k+1),W,F(k+1),F,M-k,k,ii,jj,Nmax,lflops) 
          nflops = nflops + lflops

          IF(Nmax .gt. rho) THEN
!             write(*,*) 'Nmax=', Nmax, 'swap once row', nswap, 'k=', k
             nswap = nswap + 1
             flgl = 0
             jj = jj + k 

             DO j = k+1, N  ! V(k+1:N)    = V(k+1:N) * ( (D(k+1:N)-F(ii)) / (D(k+1:N)-F(jj)) )
                temp  = NEGONE* computfd( FF,ii,j,PL,PU,DIFL,DIFR )  ! D(j)-F(ii)
                temp2 = NEGONE* computfd( FF,jj,j,PL,PU,DIFL,DIFR )  ! D(j)-F(jj)
                temp  = temp / temp2
                V(j)  = V(j) * temp                
             END DO

             U(k+1:jj-1) = U(k+1:jj-1) * ( (F(k+1:jj-1)-F(jj)) / (F(k+1:jj-1)-F(ii)) )                            
             U(jj+1:M)   = U(jj+1:M) * ( (F(jj+1:M)-F(jj)) / (F(jj+1:M)-F(ii)) ) 
             W(1:ii-1)   = W(1:ii-1) * ( (F(1:ii-1)-F(ii)) / (F(1:ii-1)-F(jj)) )
             W(ii+1:k)   = W(ii+1:k) * ( (F(ii+1:k)-F(ii)) / (F(ii+1:k)-F(jj)) )

             ! U(jj)       = A(ii) * ( (F(ii)-F(jj)) / (F(ii)-D(ii)) )  
             temp  = computfd( FF,ii,ii,PL,PU,DIFL,DIFR )  ! D(ii)-F(ii)
             U(jj) = A(ii) * ( (F(ii)-F(jj)) / temp )

             ! W(ii)       = (F(jj)-D(ii)) / A(jj)
             temp  = computfd( FF,jj,ii,PL,PU,DIFL,DIFR )  ! F(jj)-D(ii)
             W(ii) = temp / A(jj)

             DO j = 1, ii-1
                ! U(jj) = U(jj) * ( (F(ii)-F(j)) / (F(ii)-D(j)) )
                temp  = computfd( FF,ii,j,PL,PU,DIFL,DIFR )  ! F(ii)-D(j)
                U(jj) = U(jj) * ( (F(ii)-F(j)) / temp )

                ! W(ii) = W(ii) * ( (D(j)-F(jj)) / (F(j)-F(jj)) )
                temp  = NEGONE* computfd( FF,jj,j,PL,PU,DIFL,DIFR )  ! D(j)-F(jj)
                W(ii) = W(ii) * ( temp / (F(j)-F(jj)) )
             END DO
             DO j = ii+1, k
                ! U(jj) = U(jj) * ( (F(ii)-F(j)) / (F(ii)-D(j)) )
                temp  = computfd( FF,ii,j,PL,PU,DIFL,DIFR )  ! F(ii)-D(j)
                U(jj) = U(jj) * ( (F(ii)-F(j)) / temp )

                ! W(ii) = W(ii) * ( (D(j)-F(jj)) / (F(j)-F(jj)) )
                temp  = NEGONE* computfd( FF,jj,j,PL,PU,DIFL,DIFR )  ! D(j)-F(jj)
                W(ii) = W(ii) * ( temp / (F(j)-F(jj)) )
             END DO
             call iswap(PL,ii,jj)
             call iswap(PQ,ii,jj)             
             call dswap(F,ii,jj)
             call dswap(A,ii,jj)

             nflops = nflops + 10*(N-k)+8*M+18*k
          END IF ! Nmax
          
       END DO ! while

    END DO ! main loop

  END SUBROUTINE rrluCauchyEigrow

!!!!!!!!!
    SUBROUTINE rrluCauchyEigcol( D,F,A,V,DIFL,DIFR,tol,U,W,FF,PL,PU,M,N,Rk,PQ,nflops,nswap )
!
! ..Scalar parameter..
    INTEGER, INTENT(IN)          :: M,N
    INTEGER, INTENT(OUT)         :: Rk,nflops,nswap
    DOUBLE PRECISION, INTENT(IN) :: tol
! ..
! ..Array parameters..
    DOUBLE PRECISION, INTENT(IN)    :: DIFL(*),FF(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: A(*),V(*),U(*),W(*),D(*),F(*)
    INTEGER, INTENT(INOUT)          :: PL(*),PU(*),PQ(*)
!
! Purpose
! =========
! This routine computes a Rank Revealing Schur Complement, RRSC, for a Cauchy-like matrix
! via structured matrix techniques. This matrix C has dimensions M-by-N, with generators D,F,A and V.
! C(i,j) = A(i)*V(j) / ( D(i) -F(j) ), i.e. D*C - C*F = A*V^T. 
!
! C is factorized as C(PQ,:) = [I; Z] * C( PQ(1:Rk), : ), and Z is a Cauchy matrix, defined as
! D_2 * Z - Z *D_1 = U(Rk+1:M)*W(1:Rk) with D_2 = D(Rk+1:M) and D_1 = D(1:Rk). 
! Z is an (M-Rk+1)-by-Rk matrix. Note that since D_1 and D_2 are the updated singular values,
! their differences can not be computed directly. D_1(i) may equal to D_2(j), but mathematically they
! are different. 
!
! ..Parameters..
! ===============
! D  (inout)  DOULBE PRECISION Array, DIMENSION( N )
!    Row generators, the updated eigenvalues and will be modified. F(i) < D(i) < F(i+1).
! 
! F  (inout)  DOUBLE PRECISION Array, DIMENSION( M ). 
!    Column generators, the old eigenvalues, which will be changed but are useless for 
!    current compression.
!
! A  (inout)  DOUBLE PRECISION Array, DIMENSION( M ). 
!    Row generators. It is the updated z, see our paper.
!   
! V (inout)  DOUBLE PRECISION Array, DIMENSION( N ). 
!    Col generators. It is the reciprocal of 2-norm of each column of the original Cauchy-like matrix.
! 
! tol (in) DOUBLE PRECISION, for low rank approximation
!
! DIFL (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the distance, DIFL(i) = D0(i)-FF(i), positive values. The entries of DIFL are 
!      referred by PL and PU. 
!
! DIFR (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the distance, DIFL(i) = D0(i)-FF(i+1), negative values. The entries of DIFR are 
!      referred by PL and PU. 
!
! U (out) DOUBLE PRECISION Array, DIMENSION(M) 
!    Row generators of Mat_N which is a matrix in RRSC. The last M-Rk elements of U are useful.
!
! W (out) DOUBLE PRECISION Array, DIMENSION( MIN( M,N ) )
!    Column generators of Mat_N and the first Rk elements of W are useful.
!
! FF (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!     The original old singular values F, and its entries are referred by PU.
!
! PL (inout) INTEGER Array, DIMENSION(M), Row permutation.
!    PL corresponds to the row generator, D or D0. 
!
! PU (inout) INTEGER Array, DIMENSION(N), Col permutation.
!    PU corresponds to the column generator FF
!    
! M (in) INTEGER, Row dimension of C
!
! N (in) INTEGER, Col dimension of C
! 
! Rk (out) INTEGER, rank of returned matrix.
!
! PQ (inout) INTEGER array, DIMENSION( M ) 
!    Records the local row permutations, which is used to compute the inverse 
!    permutations of rows.
! 
! nflops (out) INTEGER
!
! nswap (out) INTEGER, the total number of swap
! 
! ============
! Written by S.-G. Li, On April 17th, 2013
! for tridiagonal eigenvalue problems 
! ==================================================================
!
! .. Parameters ..
    DOUBLE PRECISION  ONE, ZERO, NEGONE
    PARAMETER         ( ZERO= 0.0E+0, ONE=1.0E+0, NEGONE= -1.0E+0 )
! ..
! .. Local Scalars
    INTEGER           j,k,mn,prd,flgl,ii,jj,lflops
    DOUBLE PRECISION  rho,Amax,Ukk,Nmax,temp,temp2
! ..
! .. Intrinsic Functions ..
    INTRINSIC         MAX, ABS, MAXVAL, MOD

    U(1:M)  = A(1:M)
    Rk = MIN( M,N )
    PQ(1:M) = (/ (j, j=1,M) /)
    mn = Rk
    nswap = 0
    Amax = ZERO
    rho = 1.5E+0
    prd = 10
    nflops = 0
    
    DO k = 1, mn
       call CauchyPivtEigm(D,F,U,V,PL,PU,k,A,M,N,PQ,lflops,FF,DIFL,DIFR )
       nflops = nflops + lflops

       temp = NEGONE* computfd( FF,k,k,PU,PL,DIFL,DIFR )  ! D(k)-F(k)
       Ukk = U(k)*V(k) / temp
       Amax = MAX( Amax,ABS(Ukk) )
       IF( ABS(Ukk) .LT. tol*Amax ) THEN  ! first step converged
          ! Complete pivoting for Schur complement
          call CauchyPivtEigm_CP( D,F,U,V,PL,PU,k,A,M,N,PQ,lflops,FF,DIFL,DIFR )
          temp = NEGONE* computfd( FF,k,k,PU,PL,DIFL,DIFR )  ! D(k)-F(k)
          Ukk = U(k)*V(k) / temp
          Amax = MAX( Amax,ABS(Ukk) )
          IF( ABS(Ukk) .LT. tol*Amax ) THEN  ! final converged
             Rk = k-1
!             write(*,*) 'Rk', rk, 'mn', mn
             EXIT
          END IF
       END IF
       
       DO j = k+1, M    !  U(k+1:M) = U(k+1:M)* (D(k+1:M)-D(k)) / (D(k+1:M)-F(k))

          temp  = computdd( FF,j,k,PL,DIFL,DIFR )            ! D(j)-D(k)
          temp2 = NEGONE* computfd( FF,k,j,PU,PL,DIFL,DIFR ) ! D(j)-F(k)
          temp =  temp/ temp2
          U(j) = U(j) * Temp
       END DO

       DO j = k+1, N   !  V(k+1:N) = V(k+1:N)* (F(k+1:N)-F(k)) / (F(k+1:N)-D(k))

          temp = computfd( FF,j,k,PU,PL,DIFL,DIFR ) ! F(j)-D(k)
          V(j) = V(j) *( (F(j)-F(k) ) /temp )
       END DO
              
       DO j = 1, k-1   ! W(j) = W(j) * ( (F(k)-D(j))/ (D(k)-D(j)) )

          temp  = computfd( FF,k,j,PU,PL,DIFL,DIFR )  ! F(k)-D(j)
          temp2 = computdd( FF,k,j,PL,DIFL,DIFR )     ! D(k)-D(j)
          temp  = temp / temp2
          W(j) = W(j) * Temp
       END DO

       temp = NEGONE* computfd( FF,k,k,PU,PL,DIFL,DIFR )  !  W(k) = (D(k)-F(k))/A(k)
       W(k) = temp / A(k)

       DO j = 1, k-1   ! W(k) = W(k) * (F(j)-D(k)) / ( D(j)-D(k) )
          temp  = computfd( FF,j,k,PU,PL,DIFL,DIFR )  !  F(j)-D(k)
          temp2 = computdd( FF,j,k,PL,DIFL,DIFR )     !  D(j)-D(k)
          temp  = temp / temp2
          W(k) = W(k) * Temp
       END DO

       nflops = nflops + 6*(M-k)+5*(N-k)+2*7*(k-1)

       ! swap
       flgl = mod(k,prd)
       DO while(flgl .lt. 1)
          flgl = 1
          call searchMax2Eigm( U(k+1),W,D,FF,PL,DIFL,DIFR,M-k,k,ii,jj,Nmax,lflops )
          nflops = nflops + lflops

          IF(Nmax .gt. rho) THEN
!             write(*,*) 'Nmax=', Nmax, 'swap once column', nswap, 'k=', k
             nswap = nswap + 1
             flgl = 0
             jj = jj + k 

             DO j = k+1, N   !V(k+1:N)    = V(k+1:N) * ( (F(k+1:N)-D(ii)) / (F(k+1:N)-D(jj)) )
                temp  = computfd( FF,j,ii,PU,PL,DIFL,DIFR )  ! F(j)-D(ii)
                temp2 = computfd( FF,j,jj,PU,PL,DIFL,DIFR )  ! F(j)-D(jj)
                temp  = temp / temp2
                V(j)  = V(j) * temp
             END DO

             DO j = k+1, jj-1  !U(k+1:jj-1) = U(k+1:jj-1) * ( (D(k+1:jj-1)-D(jj)) / (D(k+1:jj-1)-D(ii)) )
                temp  = computdd( FF,j,jj,PL,DIFL,DIFR )  !  D(j)-D(jj)
                temp2 = computdd( FF,j,ii,PL,DIFL,DIFR )  !  D(j)-D(ii)
                temp  = temp / temp2
                U(j)  = U(j) * temp
             END DO
             DO j = jj+1, M  !U(jj+1:M)   = U(jj+1:M) * ( (D(jj+1:M)-D(jj)) / (D(jj+1:M)-D(ii)) )
                temp  = computdd( FF,j,jj,PL,DIFL,DIFR )  !  D(j)-D(jj)
                temp2 = computdd( FF,j,ii,PL,DIFL,DIFR )  !  D(j)-D(ii)
                temp  = temp / temp2
                U(j)  = U(j) * temp
             END DO

             DO j = 1, ii-1  ! W(1:ii-1)   = W(1:ii-1) * ( (D(1:ii-1)-D(ii)) / (D(1:ii-1)-D(jj)) )
                temp  = computdd( FF,j,ii,PL,DIFL,DIFR )  !  D(j)-D(ii)
                temp2 = computdd( FF,j,jj,PL,DIFL,DIFR )  !  D(j)-D(jj)
                temp  = temp / temp2
                W(j)  = W(j) * temp
             END DO
             DO j = ii+1, k  ! W(ii+1:k)   = W(ii+1:k) * ( (D(ii+1:k)-D(ii)) / (D(ii+1:k)-D(jj)) )
                temp  = computdd( FF,j,ii,PL,DIFL,DIFR )  !  D(j)-D(ii)
                temp2 = computdd( FF,j,jj,PL,DIFL,DIFR )  !  D(j)-D(jj)
                temp  = temp / temp2
                W(j)  = W(j) * temp
             END DO

             ! U(jj)       = A(ii) * ( (D(ii)-D(jj)) / (D(ii)-F(ii)) )  
             temp  = computdd( FF,ii,jj,PL,DIFL,DIFR )    ! D(ii)-D(jj)
             temp2 = NEGONE* computfd( FF,ii,ii,PU,PL,DIFL,DIFR )  ! D(ii)-F(ii)
             temp = temp/ temp2
             U(jj) = A(ii) * temp

             ! W(ii)       = (D(jj)-F(ii)) / A(jj)
             temp  = NEGONE* computfd( FF,ii,jj,PU,PL,DIFL,DIFR )  ! D(jj)-F(ii)
             W(ii) =  temp / A(jj) 

             DO j = 1, ii-1
                ! U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                temp  = computdd( FF,ii,j,PL,DIFL,DIFR )    ! D(ii)-D(j)
                temp2 = NEGONE* computfd( FF,j,ii,PU,PL,DIFL,DIFR )  ! D(ii)-F(j)
                temp  = temp / temp2
                U(jj) = U(jj) *temp

                ! W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
                temp  = computfd( FF,j,jj,PU,PL,DIFL,DIFR )  ! F(j)-D(jj)
                temp2 = computdd( FF,j,jj,PL,DIFL,DIFR )     ! D(j)-D(jj)
                temp  = temp / temp2
                W(ii) = W(ii) * temp
             END DO
             DO j = ii+1, k
                ! U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                temp  = computdd( FF,ii,j,PL,DIFL,DIFR )    ! D(ii)-D(j)
                temp2 = NEGONE* computfd( FF,j,ii,PU,PL,DIFL,DIFR )  ! D(ii)-F(j)
                temp  = temp / temp2
                U(jj) = U(jj) * temp

                ! W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
                temp  = computfd( FF,j,jj,PU,PL,DIFL,DIFR )  ! F(j)-D(jj)
                temp2 = computdd( FF,j,jj,PL,DIFL,DIFR )     ! D(j)-D(jj)
                temp  = temp / temp2
                W(ii) = W(ii) * temp
             END DO
             call iswap(PL,ii,jj)
             call iswap(PQ,ii,jj)
             call dswap(D,ii,jj)
             call dswap(A,ii,jj)

             nflops = nflops + 10*(N-k)+12*M+(11*2)*k
          END IF ! Nmax
          
       END DO ! while

    END DO ! main loop

  END SUBROUTINE rrluCauchyEigcol

!!!!!!
  FUNCTION computdd( FF,j,k,PU,DIFL,DIFR )
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: j,k
! ..
! .. Array Arguments ..
    INTEGER, INTENT(IN) :: PU(*)
    DOUBLE PRECISION, INTENT(IN) :: FF(*), DIFL(*),DIFR(*)
! 
! Purpose
! =========
! Computes D(j)-D(k) and D is the updated eigenvalues. FF are the old eigenvalues,
! and D is the updated eigenvalues, D(j) = D0( PU(j) ). 
!
! ..Parameters..
! ==============
! FF (in) DOUBLE PRECISION array, DIMENSION( LN )
!     FF are the old eigenvalues. LN is the size of original problem and FF 
!     is referenced by calling PU. 
!
! j  (in) INTEGER
!
! k  (in) INTEGER
!
! PU (in) INTEGER array, DIMENSION(*)
!      Permutation, corresponding to D
!
! DIFL (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the distance, DIFL(i) = FF(i)-D0(i), negative values
!
! =============
! Written by S.-G. Li, on April 16th, 2013
! for tridiagonal eigenvalue problems
! ==========================================================
!
    INTEGER          :: PWJJ,PWKK
    DOUBLE PRECISION :: computdd

    PWJJ = PU(j)
    PWKK = PU(k)
    IF( PWJJ .lt. PWKK ) THEN
       computdd = (FF(PWJJ+1)-FF(PWKK) )+DIFR(PWJJ) - DIFL(PWKK)  ! j < k, d(j)-d(k)
    ELSE
       computdd = (FF(PWJJ)-FF(PWKK+1) )+DIFL(PWJJ)-DIFR(PWKK)  ! j > k
    END IF
    
  END FUNCTION computdd

!!!!!!
  FUNCTION computfd( FF,j,k,PL,PU,DIFL,DIFR ) 
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: j, k
! ..
! .. Array Arguments ..
    INTEGER, INTENT(IN) :: PL(*), PU(*)
    DOUBLE PRECISION, INTENT(IN) :: FF(*),DIFL(*),DIFR(*)
!
! Purpose
! ========
! This routine computes F(j)-D(k), where F contains the old eigenvalues and 
! D contains the new eigenvalues. This is a PURE function. 
!
! ..Parameters..
! ==============
! FF  (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!     FF are the old singular values. LN is the size of original problem and FF 
!     is referenced by calling PL and PU. 
!
! j   (in) INTEGER
!
! k   (in) INTEGER
!
! PL  (in) INTEGER array, DIMENSION(*)
!     Row permutation
! 
! PU  (in) INTEGER array, DIMENSION(*)
!     Permutation
!
! DIFL (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the distance, DIFL(i) = FF(i)-D0(i), negative values
!
! =============
! Written by S.-G. Li, on April 16th, 2013
! for tridiagonal eigenvalue problems
! =======================================================

    INTEGER PDJJ,PWKK
    DOUBLE PRECISION :: computfd

    PDJJ = PL(j)
    PWKK = PU(k)
    IF( PDJJ .LE. PWKK ) THEN
       computfd = ( FF(PDJJ)-FF(PWKK) ) - DIFL(PWKK)  ! j <= k, F(j) - D(k)
    ELSE
       computfd = ( FF(PDJJ)-FF(PWKK+1) )-DIFR(PWKK)  ! j > k
    END IF
    
  END FUNCTION computfd

!!!!!!
  SUBROUTINE CauchyPivtEig( F,D,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR )
!
! ..Scalar Arguments..
    INTEGER, INTENT(IN)  :: k, M, N
    INTEGER, INTENT(OUT) :: nflops
! ..
! ..Array Arguments..
    DOUBLE PRECISION, INTENT(INOUT) :: D(*),F(*),U(*),V(*),A(*)
    DOUBLE PRECISION, INTENT(IN)    :: FF(*),DIFL(*),DIFR(*)
    INTEGER, INTENT(INOUT) :: PL(*),PU(*),PQ(*)
!
! Purpose
! ========
! Choose the pivoting for the k-th Schur Complement factorization. This procedure
! implements a variant of complete pivoting of S by working on its generators, 
! F, D, U and V, where S is the k-th Schur complement and S(i,j)=U(i)*V(j) /(F(i)-D(j)).
! S is an M-by-N matrix. This routine is designed for rrluCauchyEigrow. 
! Here F are the old singular values and D are the updated singular values. 
! 
! Initially U == A, and all the parameters are only permuted. This pivoting strategy 
! is a variant of complete pivoting.
!   Search the largest entry of the first column, JunkL, and
!          the largest entry of the first row, JunkU.
!     IF both of them are less than the first entry,  STOP
!     IF JunkL > JunkU   FLAG = 1
!     WHILE (1)
!        IF FLAG == 1 THEN  
!           Permute the row of JunkL to the first, and then compute the largest entry of first row, JunkU
!           IF JunkU is the first entry, STOP
!        ELSE
!           Permute the column of JunkU to the first, and then compute the largest entry of first column, JunkL
!           IF JunkL is the first entry, STOP
!        END IF
!     END WHILE
!
!  Further details
!  This pivoting strategy is different from complete pivoting and it may not choose the largest entry.
! 
! ..Parameters..
! ==============
! F (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!   The row generator of S, which are the old eigenvalues, and can be used
!   directly;
! 
! D (inout) DOUBLE PRECISION array, DIMENSION( N ) 
!   The column generator of S, the updated eigenvalues;
! 
! U (inout) DOUBLE PRECISION, DIMENSION( M )
!   The row generators of S, and only the last M-k entries are used. 
!   S is an (M-k)-by-(N-k) matrix;
! 
! V (inout) DOUBLE PRECISION, DIMENSION( N )
!   The column generators of S, and only the last N-k entries are used.
!
! PL (inout) INTEGER, DIMENSION( M )
!    The row permutations corresponding the original generators, corresponding to F;
!
! PU (inout) INTEGER, DIMENSION( N )
!    The column permutations corresponding the original generators, corresponding to D;
!
! k  (in) INTEGER
!    The k-th Schur complement
! 
! A  (inout) DOUBLE PRECISION, DIMENSION( M )
!    The original generator of U.
! 
! M  (in) INTEGER, the row dimension of current HSS block;
! 
! N  (in) INTEGER, the column dimension of current HSS block;
!    
! PQ (inout) INTEGER array, DIMENSION( M )
!    Records the local row permutation, and PQ is used to computed the inverse of row permutations
!    for current HSS block;
! 
! nflops (out) INTEGER, the floating point operations during this pivoting;
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( * )
!    The original old singular values, and its entries are referenced by PL or PU;
!
! DIFL (in) DOUBLE PRECISION, DIMENSION( NN )
!    DIFL(i) = FF(i) - D0(i), negative
! 
! ==========
! Written by S.-G. Li, on April 16th, 2013, for tridiagonal
! =============================================================
!
! .. Local Parameters ..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! .. Local Scalars ..
    INTEGER           jjL,jjU,flg,jj,pivot
    DOUBLE PRECISION  junkL,junkU,Piv
! ..
! .. Intrinsic Functions ..
    INTRINSIC    ABS,MIN

    nflops = 0
    call CauchyMaxEigmc( FF,D(k),PL(k),PU(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! largest entry of the first column
    call CauchyMaxEigm(  FF,D(k),PL(k),PU(k),DIFL,DIFR,V( k ),N-k+1,junkU,jjU )  ! largest entry of the first row
    junkL = junkL * ABS( V(k) )
    junkU = junkU * ABS( U(k) )

    ! Compute the first entry of the K-th Schur Complement
    IF( PL(k) .LE. PU(k) ) THEN
       Piv = ABS( U(k)*V(k) / ( ( F(PL(k))-FF(PU(k)) )-DIFL(PU(k)) ) ) 
    ELSE
       Piv = ABS( U(k)*V(k) / ( ( F(PL(k))-FF(PU(k)+1) )-DIFR(PU(k)) ) ) 
    END IF
    nflops = nflops + 4* (M+N-2*k+3)
    
    IF (junkL .LE. Piv .AND. junkU .LE. Piv ) THEN
       RETURN
    END IF
    
    pivot = 0   ! do not return
    flg = 0
    IF(junkL > junkU) flg = 1
    
    DO WHILE( 1 < 2)
       pivot = pivot +1
       IF( flg == 1 ) THEN
          jj = jjL
          call dswap(F,k,jj+k-1)
          call dswap(U,k,jj+k-1)
          call dswap(A,k,jj+k-1)
          call iswap(PL,k,jj+k-1)
          call iswap(PQ,k,jj+k-1)
          call CauchyMaxEigm( FF, D(k),PL(k),PU(k),DIFL,DIFR,V( k ),N-k+1,junkU,jjU )  ! row
          nflops = nflops + 4* (N-k+1)
          IF(jjU == 1) RETURN
          
          flg = 0
          CONTINUE
       END IF
       jj = jjU
       call dswap(D,k,jj+k-1)
       call dswap(V,k,jj+k-1)
       call iswap(PU,k,jj+k-1)
       call CauchyMaxEigmc( FF,D(k),PL(k),PU(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! column
       nflops = nflops + 4* (M-k+1)
       IF(jjL == 1) RETURN

       flg = 1
    END DO

  END SUBROUTINE CauchyPivtEig

!!!!!!
  SUBROUTINE CauchyPivtEigm(D,F,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR )
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: k,M,N
    INTEGER, INTENT(OUT) :: nflops
! ..
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(INOUT) :: D(*),F(*),U(*),V(*),A(*)
    DOUBLE PRECISION, INTENT(IN)    :: FF(*),DIFL(*),DIFR(*)
    INTEGER, INTENT(INOUT) :: PL(*),PU(*),PQ(*)
!
! Purpose
! ========
! Choose the pivoting for the k-th Schur Complement factorization. This procedure
! implements a variant of complete pivoting of S by working on its generators, D, F, 
! U and V, where S is the k-th Schur complement and S(i,j)=U(i)*V(j) / (D(i)-F(j) ). 
! S is an M-by-N matrix. By row and col permutations, it makes S(1,1) the 'largest' one.
! 
! Initially U == A, and all the parameters are only permuted. This pivoting strategy is 
!   Search the largest entry of the first column, JunkL, and 
!          the largest entry of the first row, JunkU.
!     IF both of them are less than the first entry,    STOP
!     IF JunkL > JunkU   FLAG = 1
!     WHILE (1)
!        IF FLAG == 1 THEN  
!           Permute the row of JunkL to the first, and then compute the largest entry of first row, JunkU
!           IF JunkU is the first entry, STOP
!        ELSE
!           Permute the column of JunkU to the first, and then compute the largest entry of first column, JunkL
!           IF JunkL is the first entry, STOP
!        END IF
!     END WHILE
!
!  Further details
!  This pivoting strategy is different from complete pivoting and it may not choose the largest entry.
! 
! ..Parameters..
! =============
! D (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!   The row generator of S, and D contains the updated eigenvalues. S is an 
!   (M-k)-by-(N-k) matrix, and only D(k:M) are used.
! 
! F (inout) DOUBLE PRECISION array, DIMENSION( N ) 
!   The column generator of S, and F contains the old singular values. Only the 
!   last N-k elements of F are used. 
! 
! U (inout) DOUBLE PRECISION, DIMENSION( M )
!   The row generators of S, only the last M-k elements are used. 
! 
! V (inout) DOUBLE PRECISION, DIMENSION( N )
!   The column generators of S, only the last N-k elements are used. 
! 
! PL (inout) INTEGER, DIMENSION( M )
!    The row permutations corresponding the original generators, and 
!    D(1:M) = D0( PL(1:M) )
!
! PU (inout) INTEGER, DIMENSION( N )
!    The column permutations corresponding the original generators, and
!    F(1:N) = FF( PU(1:N) )
!
! k  (in) INTEGER
!    The k-th Schur complement
! 
! A  (inout) DOUBLE PRECISION, DIMENSION( M )
!    The original generator of U
! 
! M  (in) INTEGER, the row dimension of HSS block row, H
! 
! N  (in) INTEGER, the column dimension of HSS block row, H
!    
! PQ (inout) INTEGER array, DIMENSION( M )
!    Records the local row permutation, and PQ is used to computed the inverse of row permutations
!    for current HSS block
! 
! nflops (out) INTEGER, the floating point operations during this pivoting
! 
! FF (in) DOUBLE PRECISION, DIMENSION( LN )
!    The old singular values of original matrix, LN is the size of original problem.
! 
! DIFL (in) DOUBLE PRECISION, DIMENSION( LN )
!    DIFL(i) = D0(i) - FF(i), positive
! 
! DIFR (in) DOUBLE PRECISION, DIMENSION( LN )
!    DIFL(i) = D0(i) - FF(i+1), negative
! 
! ==========
! Written by S.-G. Li, on April 17th, 2013
! for tridiagonal eigenvalue problems
! ==========================================
!
! ..Local Parameters..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! ..Local Scalars..
    INTEGER           jjL,jjU,flg,jj,pivot
    DOUBLE PRECISION  junkL,junkU,Piv
! ..
! ..Intrinsic Functions..
    INTRINSIC    ABS,MIN

    nflops = 0
    call CauchyMaxEigm ( FF,D(k),PU(k),PL(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! largest entry of the first column
    call CauchyMaxEigmc( FF,D(k),PU(k),PL(k),DIFL,DIFR,V( k ),N-k+1,junkU,jjU )  ! largest entry of the first row
    junkL = junkL * ABS( V(k) )
    junkU = junkU * ABS( U(k) ) 

    ! compute the first entry of k-th Schur complement
    IF( PL(k) .LT. PU(k) ) THEN
       Piv = ABS( U(k)*V(k) / ( (FF(PL(k)+1)-FF(PU(k)) )+DIFR(PL(K)) ) )  ! D(k) - F(k)
    ELSE
       Piv = ABS( U(k)*V(k) / ( (FF(PL(k))-FF(PU(k)) ) + DIFL(PL(K)) ) )  ! D(k) - F(k)
    END IF
    nflops = nflops + 4* (M+N-2*k+3)
    
    IF (junkL .LE. Piv .AND. junkU .LE. Piv ) THEN
       RETURN
    END IF
    
    pivot = 0   ! do not return
    flg = 0
    IF(junkL > junkU) flg = 1
    
    DO WHILE (1 < 2)
       pivot = pivot +1
       IF (flg == 1) THEN
          jj = jjL
          call dswap(D,k,jj+k-1)
          call dswap(U,k,jj+k-1)
          call dswap(A,k,jj+k-1)
          call iswap(PL,k,jj+k-1)
          call iswap(PQ,k,jj+k-1)
          call CauchyMaxEigmc( FF,D(k),PU(k),PL(k),DIFL,DIFR,V( k ),N-k+1,junkU,jjU )  ! row
          nflops = nflops + 4* (N-k+1)
          IF(jjU == 1) RETURN
          
          flg = 0
          CONTINUE
       END IF
       jj = jjU
       call dswap(F,k,jj+k-1)
       call dswap(V,k,jj+k-1)
       call iswap(PU,k,jj+k-1)
       call CauchyMaxEigm( FF,D(k),PU(k),PL(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! column
       nflops = nflops + 4* (M-k+1)
       IF(jjL == 1) RETURN

       flg = 1
    END DO

  END SUBROUTINE CauchyPivtEigm

!!!!!!
  SUBROUTINE CauchyPivtEigm_CP( D,F,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR )
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: k,M,N
    INTEGER, INTENT(OUT) :: nflops
! ..
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(INOUT) :: D(*),F(*),U(*),V(*),A(*)
    DOUBLE PRECISION, INTENT(IN)    :: FF(*),DIFL(*),DIFR(*)
    INTEGER, INTENT(INOUT) :: PL(*),PU(*),PQ(*)
!
! Purpose
! ========
! Choose the pivoting for the k-th Schur Complement factorization. This procedure
! implements a complete pivoting of S by working on its generators, D, F, U and V,
! where S is the k-th Schur complement and S(i,j) = U(i)*V(j) / (D(i) - F(j) ). 
! S is an (M-k)-by-(N-k) matrix. 
! 
! Initially U == A, and all the parameters are only permuted. This pivoting strategy is 
! to choose the largest entry for each column and then choose the largest one. 
!   Search the largest entry of the first column, JunkL, with index jjL
!           Junk = JunkL and jj = jjL
!   Search the largest entry of j-th column, JunkR, and  index jjR
!     IF JunkR > JunkL  THEN
!          Junk = JunkR  and  jj = jjR
!   Permute the row and column generators
!
!  Further details
!  This pivoting strategy is complete pivoting different with Gu's variant strategy.
! 
! ..Parameters..
! ==============
! D (inout) DOUBLE PRECISION array, DIMENSION( N )
!   The row generator of S, the updated eigenvalues, and its relation with PL is
!   D(1:M) = D0( PL(1:M) );
! 
! F (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!   The column generator of S, which is useless and corresponding entries
!   are referenced by FF and PU, the old eigenvalues;
! 
! U (inout) DOUBLE PRECISION, DIMENSION( M )
!   The row generators of S, and only the last M-k entries are used;
! 
! V (inout) DOUBLE PRECISION, DIMENSION( N )
!   The column generators of S, and only the last N-k entries are used;
! 
! FF (in) DOUBLE PRECISION, DIMENSION( LN )
!    The old singular values of original matrix, and its entries are referenced by PL or PU.
! 
! DIFL (in) DOUBLE PRECISION, DIMENSION( LN )
!    DIFL(i) = D0(i)-FF(i), positive
! 
! DIFR (in) DOUBLE PRECISION, DIMENSION( LN )
!    DIFL(i) = D0(i)-FF(i+1), negative
! 
! PL (inout) INTEGER, DIMENSION( M )
!    The row permutations corresponding the original generators, corresponding to D;
!
! PU (inout) INTEGER, DIMENSION( N )
!    The column permutations corresponding the original generators, corresponding to F;
!
! k  (in) INTEGER
!    The k-th Schur complement
! 
! A  (inout) DOUBLE PRECISION, DIMENSION( M )
!    The original generator of U
! 
! M  (in) INTEGER, the row dimension of current HSS block row;
! 
! N  (in) INTEGER, the column dimension of current HSS block row;
!    
! PQ (inout) INTEGER array, DIMENSION( M )
!    Records the local row permutation, and PQ is used to computed the inverse of row permutations
!    for current HSS block row
! 
! nflops (out) INTEGER, the floating point operations during this pivoting
! 
! ==========
! Written by S.-G. Li, on April 17th, 2013
! for tridiagonal eigenvalue problems
! ===========================================
!
! ..Local Parameters..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! ..Local Scalars..
    INTEGER           ii,jj,j
    DOUBLE PRECISION  junk,junkT
! ..
! ..Local Arrays..
    INTEGER jjU(N-k+1)
    DOUBLE PRECISION junkU(N-k+1)
! ..
! ..Intrinsic Functions..
    INTRINSIC    ABS,MIN

    nflops = 0
    call CauchyMaxEigm( FF,D(k),PU(k),PL(k),DIFL,DIFR,U( k ),M-k+1,junkU(1),jjU(1) )  ! largest entry of the first column
    junk = ABS( junkU(1)*V(k) )
    ii = jjU(1)
    jj = 1

!$OMP PARALLEL DO PRIVATE(j)
    DO j = 2, N-k+1
       call CauchyMaxEigm( FF,D(k),PU(k+j-1),PL(k),DIFL,DIFR,U( k ),M-k+1,junkU(j),jjU(j) )  ! largest entry of the first column
    END DO
!$OMP END PARALLEL DO

    DO j = 2, N-k+1
       junkT = ABS( junkU(j)*V(k+j-1) )
       IF( junkT > junk ) THEN
          junk = junkT
          ii = jjU(j)  ! row
          jj = j       ! col
       END IF
    END DO

    call dswap(D,k,ii+k-1)
    call dswap(U,k,ii+k-1)
    call dswap(A,k,ii+k-1)
    call iswap(PL,k,ii+k-1)
    call iswap(PQ,k,ii+k-1)

    call dswap(F,k,jj+k-1)
    call dswap(V,k,jj+k-1)
    call iswap(PU,k,jj+k-1)
    
  END SUBROUTINE CauchyPivtEigm_CP

!!!!!!
  SUBROUTINE CauchyPivtEig_CP( F,D,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR )
!
! ..Scalar Arguments..
    INTEGER, INTENT(IN)  :: k, M, N
    INTEGER, INTENT(OUT) :: nflops
! ..
! ..Array Arguments..
    DOUBLE PRECISION, INTENT(INOUT) :: D(*),F(*),U(*),V(*),A(*)
    DOUBLE PRECISION, INTENT(IN)    :: FF(*),DIFL(*),DIFR(*)
    INTEGER, INTENT(INOUT) :: PL(*),PU(*),PQ(*)
!
! Purpose
! ========
! Choose the pivoting for the k-th Schur Complement factorization. This procedure
! implements a complete pivoting of S by working on its generators, D, F, U and V,
! where S is the k-th Schur complement and S(i,j) = U(i)*V(j) / (F(i) - D(j) ).
! S is an (M-k)-by-(N-k) matrix. This routine is designed for rrluCauchyEigrow.
! Here F are the old eigenvalues and D are the updated eigenvalues. F is the
! row generator and D is the column generator. 
! 
! Initially U == A, and all the parameters are only permuted. This pivoting strategy is 
! to choose the largest entry for each column and then choose the largest one. 
!   Search the largest entry of the first column, JunkL, with index jjL
!           Junk = JunkL and jj = jjL
!   Search the largest entry of j-th column, JunkR, and  index jjR
!     IF JunkR > JunkL  THEN
!          Junk = JunkR  and  jj = jjR
!   Permute the row and column generators
!
!  Further details
!  This pivoting strategy is complete pivoting instead of Gu's variant.
! 
! ..Parameters..
! ==============
! F (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!   The row generator of S, which are the old eigenvalues, and can be used
!   directly. 
! 
! D (inout) DOUBLE PRECISION array, DIMENSION( N ) 
!   The column generator of S, the updated eigenvalues.
! 
! U (inout) DOUBLE PRECISION, DIMENSION( M )
!   The row generators of S, only the last M-k entries are used, and 
!   S is an (M-k)-by-(N-k) matrix;
! 
! V (inout) DOUBLE PRECISION, DIMENSION( N )
!   The column generators of S, and only the last N-k entries are used;
!
! PL (inout) INTEGER, DIMENSION( M )
!    The row permutations corresponding the original generators, corresponding to F;
!
! PU (inout) INTEGER, DIMENSION( N )
!    The column permutations corresponding the original generators, corresponding to D;
!
! k  (in) INTEGER
!    The k-th Schur complement
! 
! A  (inout) DOUBLE PRECISION, DIMENSION( M )
!    The original generator of U
! 
! M  (in) INTEGER, the row dimension of current HSS block row H
! 
! N  (in) INTEGER, the column dimension of current HSS block row H
!    
! PQ (inout) INTEGER array, DIMENSION( M )
!    Records the local row permutation, and PQ is used to computed the inverse of row permutations
!    for current HSS blk
! 
! nflops (out) INTEGER, the floating point operations during this pivoting
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( * )
!    The original old singular values, and its entries are referenced by PL or PU;
!
! DIFL (in) DOUBLE PRECISION, DIMENSION( LN )
!    DIFL(i) = D0(i)-FF(i), positive
! 
! DIFR (in) DOUBLE PRECISION, DIMENSION( LN )
!    DIFL(i) = D0(i)-FF(i+1), negative
! 
! ==========
! Written by S.-G. Li, on April 16th, 2013
! for tridiagonal eigenvalue problems
! ===========================================================
!
! ..Local Parameters..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! ..Local Scalars..
    INTEGER           ii,jj,j
    DOUBLE PRECISION  junk, junkT
! ..
! ..Local Arrays..
    INTEGER jjU(N-k+1)
    DOUBLE PRECISION junkU(N-k+1)
! ..
! .. Intrinsic Functions ..
    INTRINSIC    ABS,MIN

    nflops = 0
    call CauchyMaxEigmc( FF,D(k),PL(k),PU(k),DIFL,DIFR,U( k ),M-k+1,junkU(1),jjU(1) )       ! largest entry of the first column
    junk = ABS( junkU(1)*V(k) )
    ii = jjU(1)
    jj = 1
    
!$OMP PARALLEL DO PRIVATE(j)
    DO j = 2, N-k+1
       call CauchyMaxEigmc( FF,D(k),PL(k),PU(k+j-1),DIFL,DIFR,U( k ),M-k+1,junkU(j),jjU(j) )  ! largest entry of the first column
    END DO
!$OMP END PARALLEL DO

    DO j = 2, N-k+1
       junkT = ABS( junkU(j)*V(k+j-1) )
       IF( junkT > junk ) THEN
          junk = junkT
          ii = jjU(j)   ! row
          jj = j        ! col
       END IF
    END DO

    call dswap(F,k,ii+k-1)
    call dswap(U,k,ii+k-1)
    call dswap(A,k,ii+k-1)
    call iswap(PL,k,ii+k-1)
    call iswap(PQ,k,ii+k-1)

    call dswap(D,k,jj+k-1)
    call dswap(V,k,jj+k-1)
    call iswap(PU,k,jj+k-1)

  END SUBROUTINE CauchyPivtEig_CP

!!!!!!
  SUBROUTINE searchMax2Eig( U,W,F2,F1,LDU,LDW,ii,jj,Nmax,nflops )
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDU, LDW
    INTEGER, INTENT(OUT) :: ii, jj, nflops
    DOUBLE PRECISION, INTENT(OUT) :: Nmax
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN) :: U(*),W(*),F2(*),F1(*)
!
! Purpose 
! ========
! Choose the maximum entry of Mat_N by searching column by column, and Mat_N is defined as
! F2* Mat_N -Mat_N*F1 = U*W, a Cauchy-like matrix. This routine is designed for rrluCauchyEigrow
! and here F1 and F2 are the old eigenvalues, which can be used directly. 
! F2 is a LDU-by-1 vector and F1 is a LDW-by-1 vector. 
! 
! ..Parameters..
! ================
! U  (in) DOUBLE PRECISION array, DIMENSION( LDU )
!    The row generators of Mat_N
! 
! W  (in) DOUBLE PRECISION array, DIMENSION( LDW )
!    The column generators of Mat_N
! 
! LDU (in) INTEGER, row dimension of Mat_N
! 
! LDW (in) INTEGER, column dimension of Mat_N
!
! ii (out) INTEGER, the column index of largest entry in Mat_N
!
! jj (out) INTEGER, the row index of largest entry in Mat_N
!
! Nmax (out) DOUBLE PRECISION, the largest entry in Mat_N
! 
! nflops (out) INTEGER 
! 
! ==========
! Written by S.-G. Li, on April 16th, 2013
! for tridiagonal eigenvalue problems
! ========================================
!
! .. Local Parameters ..
    DOUBLE PRECISION ZERO
    PARAMETER        ( ZERO = 0.0E+0 )
! .. Local Scalars ..
    INTEGER          :: j, jjL
    DOUBLE PRECISION :: junk

    nflops = 0    
    call CauchyMaxEig( F2,F1(1),U,LDU,junk,jjL )  ! first column
    junk = junk*ABS( W(1) )       
    Nmax = junk
    ii = 1
    jj = jjL
    nflops = nflops + 4* LDU*LDU
    
    DO j = 2, LDW   ! ii: col, jj: row
       call CauchyMaxEig( F2,F1(j),U,LDU,junk,jjL )  ! jth col
       junk = junk*ABS( W(j) )       
       IF(junk .GT. Nmax) THEN
          Nmax = junk
          ii = j
          jj = jjL
       END IF
    END DO

  END SUBROUTINE searchMax2Eig

!!!!!!
  SUBROUTINE searchMax2Eigm( U,W,D,FF,PL,DIFL,DIFR,LDU,LDW,ii,jj,Nmax,nflops )
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDU, LDW
    INTEGER, INTENT(OUT) :: ii, jj, nflops
    DOUBLE PRECISION, INTENT(OUT) :: Nmax
! ..
! .. Array Arguments ..
    INTEGER, INTENT(IN)  :: PL(*)
    DOUBLE PRECISION, INTENT(IN) :: U(*),W(*),D(*),FF(*),DIFL(*),DIFR(*)
!
! Purpose 
! ========
! Choose the maximum entry of Z by searching column by column, and Z is defined as
! D2* Z -Z*D1 = U*W, a Cauchy-like matrix. 
! D1 and D2 are updated eigenvalues, and therefore their difference should be
! computed via old eigenvalues and DIFL and PL.
! D2 is a LDU-by-1 vector and D1 is a LDW-by-1 vector. D2 and D1 is referred by DD and PL,
! D1(1)= DD( PL(1) ), D2(1) = DD( PL(k+1) ). 
!
! .. Parameters ..
! ================
! U  (in) DOUBLE PRECISION array, DIMENSION( LDU )
!    The row generators of Z
! 
! W  (in) DOUBLE PRECISION array, DIMENSION( LDW )
!    The column generators of Z
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( N )
!    The original old eigenvalues, and its entries are referred by PL.
!    It is used to compute  DD(i)-DD(j).
!
! PL (in) INTEGER array, DIMENSION( LDW+LDU )
!    The permutations of D0 for current node
!
! DIFL (in) DOUBLE PRECISION array, DIMENSION( N )
!      The original differences between updated singular values and old ones,
!      DIFL(i) = FF(i) - D0(i), negative
!
! LDU (in) INTEGER, row dimension of Z
! 
! LDW (in) INTEGER, column dimension of Z
!
! ii  (out) INTEGER 
!     The column index of chosen largest entry in Z
! 
! jj  (out) INTEGER 
!     The row index of chosen largest entry in Z
! 
! ==========
! Written by S.-G. Li, on April 17th, 2013
! for tridiagonal eigenvalue problems
! =========================================
!
! .. Local Parameters ..
    DOUBLE PRECISION ZERO
    PARAMETER        ( ZERO = 0.0E+0 )
! ..
! .. Local Scalars ..
    INTEGER          :: j, jjL
    DOUBLE PRECISION :: junk

    nflops = 0    
    call CauchyMaxEigmm(D,FF,PL,1,DIFL,DIFR,U,LDU,LDW,junk,jjL)  ! first col
    junk = junk*ABS( W(1) )       
    Nmax = junk
    ii = 1   ! col
    jj = jjL ! row
    nflops = nflops + 4* LDU*LDU
    
    DO j = 2, LDW   ! ii: col, jj: row
       call CauchyMaxEigmm(D,FF,PL,j,DIFL,DIFR,U,LDU,LDW,junk,jjL) ! jth col
       junk = junk*ABS( W(j) )       
       IF(junk .GT. Nmax) THEN
          Nmax = junk
          ii = j
          jj = jjL
       END IF
    END DO

  END SUBROUTINE searchMax2Eigm

!!!!!!
  SUBROUTINE CauchyMaxEig( D,F,U,N,junk,jj )
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(OUT) :: jj
    DOUBLE PRECISION, INTENT(IN)  :: F
    DOUBLE PRECISION, INTENT(OUT) :: junk
! 
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN) :: D(*)
    DOUBLE PRECISION, INTENT(IN) :: U(*)
! 
! Purpose
! =======
! This routine returns the largest entry in magnitude and its index in the first column of A,
! where A is a Cauchy-like matrix, defined as A(i,j) = U(i)*V(j) / (D(i) - F(j) ). 
! Here F = F(j), and D and F are only variables, both can be used directly.
! The leading dimension of A is N. 
!
! ===========
! Written by S.-G. Li, on April 16th, 2013, in Changsha China
! ==========================================================
!
! .. Local Scalars ..
    INTEGER temp(1)
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MAXLOC, ABS,MAXVAL
! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: LK(:)
    ALLOCATE( LK(N) )

    Lk = U(1:N)/ ( D(1:N)-F )
    Lk = ABS( Lk )
    junk = MAXVAL( Lk(1:N) )
    temp = MAXLOC( Lk(1:N) )
    jj = temp(1)
    
    DEALLOCATE( LK )
    
  END SUBROUTINE CauchyMaxEig

!!!!!!
  SUBROUTINE CauchyMaxEigmm( D,FF,PL,k,DIFL,DIFR,U,M,N,junk,jj )
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(OUT) :: jj
    INTEGER, INTENT(IN)  :: N, M, k
    DOUBLE PRECISION, INTENT(OUT) :: junk
! 
! .. Array Arguments ..
    INTEGER, INTENT(IN) :: PL(*)
    DOUBLE PRECISION, INTENT(IN) :: D(*),FF(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(IN) :: U(*)
! 
! Purpose
! =======
! This routine returns the largest entry in magnitude and its index of k-th column of matrix A,
! where A is a Cauchy-like matrix, defined as A(i,j) = U(i)*V(j) / (D2(i) - D1(j) ),
! an M-by-N matrix. D1(1:N) = D0(PL(1:N))=D(1:N), D2(1:M) = D0(PL(N+1:M+N))=D(N+1: M+N), 
! and their entries are updated eigenvalues, and FF are old eigenvalues which are 
! used to compute the difference between D1(i) and D2(j). 
!
! A(i,j) = U(i)*V(j) / (D(i+N)-D(j)), and D is not used. 
!
! ..Parameters..
! ==============
! D (in) DOUBLE PRECISION array, DIMENSION( M+N )
!    The updated singular values, and its entries are corresponding to PL; 
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( * )
!    The original old singular values, and its entries are referred by PL; 
! 
! PL (in) INTEGER, DIMENSION(N)
!    The permutation of updated singular values, the row generators of A
! 
! k  (in) INTEGER 
!    The index of computed column of A
! 
! DIFL (in) DOUBLE PRECISION array, DIMENSION( * )
!      The original differences between updated eigenvalues and old ones,
!      DIFL(i) = FF(i) - D0(i), negative
!
! U   (in) DOUBLE PRECISION array, DIMENSION( M ) 
!     The row generators of A
!
! M   (in) INTEGER
!     The row dimension of A
! 
! N   (in) INTEGER
!     The column dimension of A
! 
! Junk (out) DOUBLE
!      The absolute value of largest entry in the k-th column of A
! 
! jj   (out) INTEGER 
!      The row index of Junk
! ===========
! Written by S.-G. Li, on April 17th, 2013, in Changsha China
! for tridiagonal eigenvalue problems
! ==========================================================
!
! ..Local Scalars..
    INTEGER temp(1), I
    DOUBLE PRECISION TT
! ..
! ..Intrinsic Functions..
    INTRINSIC    MAXLOC, ABS,MAXVAL
! ..
! ..Local Arrays..
    DOUBLE PRECISION, ALLOCATABLE :: LK(:)
    ALLOCATE( LK(M) )

    DO i = 1, M
       TT = computdd(FF,i+N,k,PL,DIFL,DIFR)  ! D(i+N)-D(k)
       LK(i) = U(i) / TT
    END DO
    Lk = ABS( Lk )
    junk = MAXVAL( Lk(1:M) )
    temp = MAXLOC( Lk(1:M) )
    jj = temp(1)
    
    DEALLOCATE( LK )
    
  END SUBROUTINE CauchyMaxEigmm

!!!!!!
  SUBROUTINE CauchyMaxEigm( FF,D,Fk,PU,DIFL,DIFR,V,N,junk,jj )
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(IN)  :: Fk
    INTEGER, INTENT(OUT) :: jj
    DOUBLE PRECISION, INTENT(OUT) :: junk
! 
! .. Array Arguments ..
    INTEGER, INTENT(IN) :: PU(*)
    DOUBLE PRECISION, INTENT(IN) :: D(*),FF(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(IN) :: V(*)
! 
! Purpose
! =======
!  This routine returns the largest entry in magnitude and its index of one row of matrix A,
!  where A is a Cauchy-like matrix, defined as A(i,j) = U(i)*V(j) / (F(i) - D(j) ). 
!  Here Fk = F(i), and Fk is a scalar, one of the old singular values. D are 
!  the updated eigenvalues, and the relation between D and D0 is D(1:N) = D0( PU(1:N) ). 
!  The leading dimension of A is N. 
!
! ..Parameters..
! ==============
! FF (in) DOUBLE PRECISION array, DIMENSION( LN )
!    LN is the size of original problem. FF contains the original eigenvalues.
! 
! D  (in) DOUBLE PRECISION, DIMENSION( N )
!    The updated singular values, its entries have relationship with PU by 
!    D(1:N) = D0( PU(1:N) ). 
!
! PU (in) INTEGER, DIMENSION(N)
!    The permutation of updated eigenvalues, the column generators of A
! 
! Fk  (in) INTEGER 
!     The original position of F(j) in FF, i.e., FF(Fk) = F;
! 
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!     DIFL(i) = D0(i)-FF(i), the order of which will not change and its elements
!     are referenced by using PL or PU. Positive
!
! DIFR (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!     DIFR(i) = D0(i)-FF(i+1), the order of which will not change and its elements
!     are referenced by using PL or PU. Negative
!
! V  (in) DOUBLE PRECISION array, DIMENSION( N )
!     The row generator of A.
!      
! N  (in) INTEGER, column dimension of A.
! 
! junk (out) DOUBLE PRECISION, the computed largest entry of magnitude.
!
! jj (out) INTEGER, the column dimension of computed largest entry
! 
! ===========
!  Written by S.-G. Li, on April 16th, 2013, in Changsha China
!  for tridiagonal eigenvalue problem
! ==========================================================
!
! .. Local Scalars ..
    INTEGER temp(1), I, PUi
    DOUBLE PRECISION F
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MAXLOC, ABS,MAXVAL
! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: LK(:)
    ALLOCATE( LK(N) )

    F = FF(Fk)
    DO i = 1, N
       PUi = PU(i)
       IF( Fk .le. PUi ) THEN
          LK(i) = V(i) / ( (F -FF(PUi) )-DIFL(PUi) ) ! F(k) - D(i)
       ELSE
          LK(i) = V(i) / ( (F-FF(PUi+1) )-DIFR(PUi) ) ! F(k) - D(i)
       END IF
    END DO
    Lk = ABS( Lk )
    junk = MAXVAL( Lk(1:N) )
    temp = MAXLOC( Lk(1:N) )
    jj = temp(1)
    
    DEALLOCATE( LK )
    
  END SUBROUTINE CauchyMaxEigm

!!!!!!
  SUBROUTINE CauchyMaxEigmc( FF,D,PL,PU,DIFL,DIFR,U,N,junk,jj )
! 
! ..Scalar Arguments..
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(OUT) :: jj
    DOUBLE PRECISION, INTENT(IN)  :: D
    DOUBLE PRECISION, INTENT(OUT) :: junk
! 
! ..Array Arguments..
    INTEGER, INTENT(IN) :: PL(*),PU(*)
    DOUBLE PRECISION, INTENT(IN) :: FF(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(IN) :: U(*)
! 
! Purpose
! =======
! This routine returns the largest entry of magnitude and its index of one column of matrix A,
! where A is a Cauchy-like matrix, defined as A(i,j) = U(i)*V(j) / (F(i) - D(j) ). 
! Here D = D(j), F is referred by FF, F(1:N) = FF( PL(1:N) ). D is a scalar, one of the updated singular
! values, and FF are the old eigenvalues. 
!
! D is useless in this routine. 
!
! ..Parameters..
! ==============
! FF  (in) DOUBLE PRECISION, DIMENSION( N )
!    The original eigenvalues, its entries are referred by PL;
! 
! D (in) DOUBLE PRECISION scalar
!    One of the original updated singular values, the column generator of A,
!    and D = D0( PU ).
!
! PL (in) INTEGER array, DIMENSION( N )
!    The permutation of old eigenvalues, such that F = FF(PL);
! 
! PU (in) INTEGER 
!    The permutation of updated eigenvalues, the column generators of A;
! 
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!     DIFL(i) = FF(i)-D0(i), the order of which will not change and its elements
!     are referenced by using PL or PU, Positive.
!
! DIFR (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!     DIFL(i) = D0(i)-FF(i+1), the order of which will not change and its elements
!     are referenced by using PL or PU, Negative.
!
! U  (in) DOUBLE PRECISION array, DIMENSION( N )
!    The row generator of A;
! 
! junk (out) DOUBLE PRECISION, the computed largest entry of magnitude;
! 
! jj  (out) INTEGER, the column index of computed largest entry;
! 
! ===========
! Written by S.-G. Li, on April 16th, 2013, in Changsha China
! 
! ==========================================================
!
! ..Local Scalars..
    INTEGER temp(1), I, PLi, Dk
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MAXLOC, ABS,MAXVAL
! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: LK(:)
    ALLOCATE( LK(N) )

    Dk = PU(1)
    DO i = 1, N
       PLi = PL(i)
       IF( PLi .le. Dk ) THEN
          LK(i) = U(i) / ( (FF(PLi)-FF(Dk))-DIFL(Dk) )
       ELSE
          LK(i) = U(i) / ( (FF(PLi) - FF(Dk+1))-DIFR(Dk)  )
       END IF
    END DO
    Lk = ABS( Lk )
    junk = MAXVAL( Lk(1:N) )
    temp = MAXLOC( Lk(1:N) )
    jj = temp(1)
    
    DEALLOCATE( LK )
    
  END SUBROUTINE CauchyMaxEigmc

!!!!!!
  SUBROUTINE dswap(D,mk,nk)
! 
! .. Purpose 
! Switch the mk-th and nk-th entry of D, and D is a 1-D double precision array
! 
    INTEGER, INTENT(IN) :: mk,nk
    DOUBLE PRECISION, INTENT(INOUT) :: D(*)
!
    DOUBLE PRECISION tt

    tt = D(mk)
    D(mk) = D(nk)
    D(nk) = tt
    
  END SUBROUTINE dswap

  SUBROUTINE iswap(D,mk,nk)
! 
! .. Purpose 
! Switch the mk-th and nk-th entry of D, which is a 1-D integer array
! 
    INTEGER, INTENT(IN) :: mk,nk
    INTEGER, INTENT(INOUT) :: D(*)
!
    INTEGER tt

    tt = D(mk)
    D(mk) = D(nk)
    D(nk) = tt
    
  END SUBROUTINE iswap

!!!!!!!!
  SUBROUTINE showMaxMinVec( UVW, A, LENA, SK, EK )
! ..
! ..Scalar Arguments..
    CHARACTER(LEN=1) :: UVW
    INTEGER, INTENT(IN) :: LENA, SK, EK
! ..Array Arguments..
    DOUBLE PRECISION, INTENT(IN) :: A(*)
!
! Purpose
! =======
! Show the largest and smallest entry of A(SK: EK), and assume A at least has
! LENA entries, LENA=EK-SK+1. 
! 
! =========
! Written on April 16th, 2013
! ===========================
    DOUBLE PRECISION, PARAMETER :: STOL = 1.0D-12
    DOUBLE PRECISION, PARAMETER :: LTOL = 1.0D+10
!
    DOUBLE PRECISION junk1, junk2, W(LENA)
! ..
! .. INTRINSIC FUNCTIONS ..
    INTRINSIC  ABS, MAXVAL,MINVAL
!
    W(1:LENA) = ABS( A(SK : EK) )
    junk1 = maxval(W)
    junk2 = minval(W)

    IF( junk1 .ge. LTOL) THEN
       write(*,*) 'Max entry of ', UVW, ' is ', junk1
    END IF

    IF( junk2 .le. STOL) THEN
       write(*,*) 'Min entry of ', UVW, ' is ', junk2
    END IF
    
  END SUBROUTINE showMaxMinVec

!!!!!!
  SUBROUTINE invp(P,N)
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: N
! .. Array Arguments ..
    INTEGER, INTENT(INOUT) :: P(*)
! 
! Purpose
! =======
! This routine computes the inverse permutation of defined in P.
! 
! ======
! Modified on Dec. 24th, 2012
! ==========================
!
! .. Local Parameters ..
    INTEGER :: i
    INTEGER :: IP(N)
    
    IP = 0
    DO i = 1, N
       IP( P(i) ) = i
    END DO
    P(1:N) = IP(1:N)

  END SUBROUTINE invp

!!!!!!
  SUBROUTINE CauchylikeEig( A,F,D,U,V,DIFL,DIFR,PL,PU,M,N,nflops )
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  ::  M, N
    INTEGER, INTENT(OUT) ::  nflops
! .. Array Arguments ..
    INTEGER, INTENT(IN)  ::   PL(*),PU(*)
    DOUBLE PRECISION, INTENT(IN)  :: D(*),F(*),U(*),V(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(OUT) :: A(M,N)
!
! Purpose 
! =======
! It computes an Cauchy matrix A with dimension M-by-N, 
! A(i,j) = U(i)*V(j) / ( F(i) - D(j) ). The entries of F are the old eigenvalues,
! and D contains the updated eigenvalues. The computed eigenvector
! matrix satisfies F*A-A*D = U*V^T. 
!
! Parameters
! F   (input) DOUBLE PRECISION array, DIMENSION ( LN )
!     The original old eigenvalues which are used directly, the row generators.
!     Its entries are referenced by PL.
!
! D   (input)  DOUBLE PRECISION array, DIMENSION ( N )
!     D contains the updated eigenvalues for current HSS block,
!     D corresponds to PU, and their relation is D(i) = D0( PU(i) );
!
! U   (input) DOUBLE PRECISION array, DIMENSION ( M )
!     Row Generator of Cauchy matrix, which are referenced directly.
! 
! V   (input) DOUBLE PRECISON array, DIMENSION ( N )
!     Column Generator of Cauchy matrix, which are referenced directly.
!
! DIFL (input) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFL(i) = D0(i)-F(i), the order of which will not change and its elements
!      are also referenced by using PL or PU. Positive
!
! DIFR (input) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFR(i) = D0(i)-F(i+1), the order of which will not change and its elements
!      are also referenced by using PL or PU. Negative
!
! PL  (input) INTEGER array, DIMENSION ( N )
!     The row permutations of generator.  It relates with F.
!
! PU  (input) INTEGER array, DIMENSION ( N )
!     The column permutations of generator. It relates with D.
!
!  M  (input) INTEGER, the row dimension of A
!
!  N  (input) INTEGER, the column dimension of A
! 
! nflops (output) INTEGER 
!        Floating point operations of generating A.
!
! =============
!  Written by S.-G. Li, on April 17th, 2013
!  for tridiagonal eigenvalue problems
! =======================================
!
! .. Local Scalars ..
    INTEGER           i,j,PFII,PDJJ
    DOUBLE PRECISION  FJ,VJ

    nflops = 6*M*N
    DO j = 1,N
       PDJJ = PU(j)
       FJ = F(PDJJ)
       VJ = V(j)
       DO i = 1, M
          PFII = PL(i)
          IF( PFII .le. PDJJ ) THEN
             A(i,j) = U(i)*VJ / ( (F(PFII)-FJ )-DIFL(PDJJ) )
          ELSE
             A(i,j) = U(i)*VJ / ( ( F(PFII)-F(PDJJ+1) )-DIFR(PDJJ) )
          END IF
       END DO
    END DO
    
  END SUBROUTINE CauchylikeEig

!!!!!!
  SUBROUTINE CauchylikeEigQ2( A,D1,D2,W,Z,M,N,FF,DIFL,DIFR,PL,nflops )
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: M, N
    INTEGER, INTENT(OUT) :: nflops
! ..
! .. Array Arguments ..
    INTEGER, INTENT(IN)  :: PL(*)
    DOUBLE PRECISION, INTENT(IN)  :: D1(*),D2(*),W(*),Z(*),FF(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(OUT) :: A(M,N)
!
! Purpose
! ========
! It returns an interpolative matrix, A = [I MatN], where I is an M-by-M identity matrix and
! MatN is an M-by-(N-M) Cauchy-like matrix, satisfying D1*MatN - MatN*D2 = W*Z.
! MatN(i,j) = W(i)*Z(j) / (D1(i)-D2(j) ), D2 is the second part of updated svals, and 
! D1 is the first part of updated svals.  
! A would be a M-by-N fat matrix and the interpolative matrix for an HSS block column.
!
! More details: 
! 
! This routine is used to construct an interpolative matrix for an HSS block column, H, 
! such that F*H-H*D = U*V^T, and first H is first transposed, get B, 
! D*B-B*F = (-V)*U^T. B is compressed by rrluCauchyEigcol to
! B == PQ*[I; MatN^T] * B1, such that D1*B1 - B1*F= (-V1)*U^T, and
! D2 * MatN^T -MatN^T * D1 = (-Z2)*W, initially Z = U, and Z2 is the second part of Z,
! B1=B(PQ(1:Rk),:). And then, D1 * MatN - MatN * D2 = W*Z2. 
! The dummy variable Z is indeed Z2. 
! 
! Parameters
! ==========
! A  (out) DOUBLE PRECISION array, DIMENSION( M,N )
!    The interpolative matrix, A = [I MatN]. It is a fat matrix. 
!    N >= M.
!
! D1 (in) DOUBLE PRECISION array, DIMENSION( M )
!    The row generators of MatN, D_1
!
! D2 (in) DOUBLE PRECISION array, DIMENSION( N-M )
!    The column generators of MatN, D_2
! 
! W  (in) DOUBLE PRECISION array, DIMENSION( M )
!    The row generators of MatN
!
! Z  (in) DOUBLE PRECISION array, DIMENSION( N-M )
!    The column generators of MatN
! 
! M  (in) INTEGER, row dimension of A
!    
! N  (in) INTEGER, column dimension of A
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( LN )
!    The original F, which is referenced by calling PL.
! 
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN )
!      The distances between FF(i) and D0(i), negative
!
! PL (in) INTEGER array, DIMENSION( M )
!    It is used to compute the subtraction of D1(i)-D2(j) which must be computed by 
!    using FF and DIFL
!
! nflops (out) INTEGER
!
! =============
! Written by S.-G. Li, on April 17th, 2013
! ========================================
! 
! .. Local Parameters ..
    DOUBLE PRECISION  ZERO, ONE, NEGONE
    PARAMETER         ( ZERO = 0.0E+0, ONE=1.0E+0, NEGONE=-1.0E+0 ) 
! ..
! .. Local Scalars ..
    INTEGER i,j,PDJJ, PDII
    
    nflops = 0
    A(1:M,1:N) = ZERO
    
    IF(M .EQ. N) THEN
       DO j = 1, N
          A(j,j) = ONE
       END DO
    ELSE
       DO j = 1,M
          A(j,j) = ONE
       END DO
       DO j = M+1, N
          PDJJ = PL(j)
          DO i = 1, M
             PDII = PL(i)
             IF( PDII .LT. PDJJ ) THEN
                A(i,j) = NEGONE*W(i)*Z(j-M) / &
                     ( (FF(PDII+1)-FF(PDJJ) )+( DIFR(PDII)-DIFL(PDJJ) ) ) ! D1(i) - D2(j)
             ELSE
                A(i,j) = NEGONE*W(i)*Z(j-M) / &
                     ( (FF(PDII)-FF(PDJJ+1) )+( DIFL(PDII)-DIFR(PDJJ) ) ) ! D1(i) - D2(j)
             END IF
          END DO
       END DO
       nflops = 6*N*(M-N)           
    END IF
    
  END SUBROUTINE CauchylikeEigQ2

!!!!!!
  SUBROUTINE CauchylikeEigQ(A,D2,D1,U,V,M,N,nflops)
!
! ..Scalar Arguments..
    INTEGER, INTENT(IN)  :: M, N
    INTEGER, INTENT(OUT) :: nflops
! ..Array Arguments..
    DOUBLE PRECISION, INTENT(IN)  :: D2(*),D1(*),U(*),V(*)
    DOUBLE PRECISION, INTENT(OUT) :: A(M,N)
! 
! Purpose
! ========
! It returns an interpolative matrix, A=[I; Z] where I is an N-by-N identity matrix and Z
! is an (M-N)-by-N. Z is a Cauchy-like matrix defined as (D2) * Z - Z* (D1) = U*V.
! M >= N, and A would be a tall and skiny matrix. and Z is defined as 
! Z(i,j) = U(i)*V(j) / (D2(i) -D1(j) ).
!
! This routine is for the compression of HSS block row, and 
! both D2 and D1 are the original F which can be computed directly. 
! That is why we do not use DIFL in this routine. 
!
! Parameters
! ==========
! A  (out) DOUBLE PRECISION array, DIMENSION( M, N )
!    The interpolative matrix, a tall matrix, A = [I; Z ]
!
! D2  (in) DOUBLE PRECISION array, DIMENSION( M-N )
!    The row generators of Z, the second part of original F
! 
! D1  (in) DOUBLE PRECISION array, DIMENSION( N )
!    The column generators of Z, the fist part of original F
! 
! U  (in) DOUBLE PRECISION array, DIMENSION( M-N )
!    The row generators of Z
! 
! V  (in) DOUBLE PRECISION array, DIMENSION( N )
!    The column generators of Z
! 
! M  (in) INTEGER, row dimension of A
! 
! N  (in) INTEGER, column dimension of A, M >= N
! 
! nflops (in) INTEGER 
!
! ============
! Written by S.-G. Li, on April 17th, 2013
! =======================================
!
! .. Parameters ..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO = 0.0E+0, ONE=1.0E+0 ) 
! .. Local Variables    
    INTEGER  i,j
    
    nflops = 0
    A(1:M,1:N) = ZERO
    
    IF(M .EQ. N) THEN
       DO j = 1, N
          A(j,j) = ONE
       END DO
    ELSE
       DO j = 1, N
          A(j,j) = ONE
          DO i = N+1, M
             A(i,j) = U(i-N)*V(j) / (D2(i-N)-D1(j) )
          END DO
       END DO
       nflops = nflops + 5*M*(N-M)
    END IF
        
  END SUBROUTINE CauchylikeEigQ

!!!!!!!!!
  SUBROUTINE Cauchy2hssEig_omp(F,D0,U,V,DIFL,DIFR,LDU,TR,LTR,M,LM,PH,H,DD,TOL,TRE,lvl,pnh,nflops,nswap)
    USE aux_hss
!
!  .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDU, LTR, lm, lvl
    INTEGER, INTENT(OUT) :: pnh, nflops, nswap
    DOUBLE PRECISION, INTENT(IN)  :: TOL
!
!  .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN)    :: D0(*),F(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: U(*),V(*) 
    DOUBLE PRECISION, INTENT(INOUT) :: H(*), DD(*)
    INTEGER, INTENT(IN)             :: TR(*), M(*)
    TYPE(HSSMM), INTENT(INOUT)      :: PH
    type(HTinfo)                    :: TRE
!
! Purpose
! =========
! Construct an HSS matrix approximation to the eigenvectors of a diagonal matrix M with a rank-one modification.
! The diagonals of M are F, the old eigenvalues. The eigenvector matrix is  Cauchy-like matrix 
! A( i,j ) = U(i)*V(j) / (F(i) - D(j) ), i.e. F*A-A*D= u*v^T, 
! Note that F is the row generator and D is the column generator, the updated eigenvalues. 
! 
! This routine returns an HSS matrix approximation to A, stored in H and DD. 
! 
! Parameters
! ========== 
! F   (in) DOUBLE PRECISION array, DIMENSION( N )
!     Its entries are the old eigenvalues and the row generators of A. 
!
! D0   (in)  DOUBLE PRECISION array, DIMENSION( N ). 
!     Its entries are the newly computed eigenvalues and the column generators of A. 
!
! U   (inout) DOUBLE PRECISION array, DIMENSION( N )
!     The row generators of A, which are U(i) = z(i), see our paper. 
!
! V   (inout)  DOUBLE PRECISION array, DIMENSION( N )
!     The column generators of A, the reciprocal of the 2-norm of each column of A
!
! DIFL  (in) DOUBLE PRECISION array, DIMENSION( N )
!       DIFL(i) is the distance between F(i) and D(i), DIFL(i) = D(i)-F(i), a positive value.
!
! DIFR  (in) DOUBLE PRECISION array, DIMENSION( N )
!       DIFR(i) is the distance between F(i) and D(i), DIFR(i) = D(i)-F(i+1), a negative value.
!
! LDU   (in) INTEGER, Row dimension of A
!
! TR    (in) INTEGER array, DIMENSION( LTR ), The HSS tree in post ordering.
!       
! LTR   (in) INTEGER, length of HSS Tree
!
! M     (in) INTEGER array, DIMENSION( LM ), the block partion of rows and columns
! 
! LM    (in) INTEGER, number of leaf nodes
!
! PH    (inout) HSSMM TYPE, contains the information of generators in H and D
!         1) starting position of each generator, D, U, V, B
!         2) dimensions of each generator
!       Note that the generators are U, V^T, D and B. We get V^T not V. 
!
! DD    (out) DOUBLE PRECISION array, DIMENSION( sum(M.^2) ), the diagonal blocks
! 
! H     (out) DOUBLE PRECISION array, DIMENSION( LDA *K * alpha )
!        alpha is a constant; K is the HSS rank; 
! 
! TOL   (in) DOUBLE PRECISION, tolerance for low-rank approximation
!
! TRE   (in) TYPE HTinfo, stores the structure of HSS tree
!
! lvl   (in) INTEGER, the total level of HSS tree.
! 
! pnh   (out) INTEGER, the total storage of H which stores U, B and V^T.
!
! nflops (out) INTEGER, the total flops of constructing this HSS matrix 
!
! nswap  (out) INTEGER, the total swaps whiling constructing this HSS matrix
!
! ===================
!  Written by Shengguo Li, on April 17th, 2013
!  for tridiagonal eigenvalue problem
! =====================================================================================
!
! ..Parameters..
    DOUBLE PRECISION :: ZERO, ONE
    PARAMETER        ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! ..Local Scalars..
    INTEGER       :: ch1,ierr,info,nn,mi,ni,nt,Rk,it,n1,lca,lda,iinfo
    INTEGER       :: i,j,lt,ns,lflops,lswap,pdt,LL,nodi,LLt,nodj
    double precision :: temp
! ..
! ..Local Arrays..
    INTEGER, ALLOCATABLE :: ch(:,:),lr(:,:),lc(:,:),indx(:)
    INTEGER, ALLOCATABLE :: PL(:),PU(:),PLi(:),PUi(:),pnd(:)
    DOUBLE PRECISION, ALLOCATABLE :: Di(:),Fi(:),Ui(:),Vi(:),D(:)
! ..
! ..External Subroutines..
    EXTERNAL  DGEMM, DLACPY
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MIN,ABS,MAXVAL,SQRT,MINVAL
    
    lda = ldu
    lca = lda
    nflops = 0
    nswap = 0
    pnh = 1               ! 'pointer' in H
    ALLOCATE( D(lda),PL(lda),PU(lda),ch(2,ltr),lr(ltr,2),lc(ltr,2),pnd(ltr),stat=ierr )
    IF( ierr /= 0 ) THEN
       WRITE(*,*) "Allocate failed in cauchy2hss! "
       RETURN
    END IF

    ALLOCATE( indx(lda),Di(lda),Fi(lda),Ui(lda),Vi(lda),PLi(lda),PUi(lda),stat=ierr )
    IF( ierr /= 0 ) THEN
       WRITE(*,*) "Allocate failed in cauchy2hss! 2 "
       RETURN
    END IF

    PL(1:ldu) = (/ (j, j=1,ldu) /)
    PU(1:ldu) = (/ (j, j=1,ldu) /)
    D(1:lda) = D0(1:lda)

    call child(tr, ch, ltr)
    lr(1,1:2) = (/1,m(1)/)
    lt = 1
    it = 1
    pdt = 1
    do i = 1, ltr
       if( ch(1,i) == 0 .and. ch(2,i) == 0 ) then
          lr(i,1:2) = (/ lt,lt+m(it)-1 /)
          lt = lr(i,2)+1
          pnd(i) = pdt
          pdt = pdt + m(it)*m(it)
          it = it+1
       else
          lr(i,1:2) = (/ lr(ch(1,i),1), lr(ch(2,i),2) /)
       endif
    enddo
    lc(1:ltr, 1:2 ) = lr(1:ltr, 1:2)

! *****************************************************
!          Compute the main diagonals D_i             *
! *****************************************************
!$OMP PARALLEL DO PRIVATE( mi,i )
    do i = 1, ltr
       if ( ch(1,i) == 0 .and. ch(2,i) == 0 ) then
          mi = lr(i,2)-lr(i,1)+1
          call CauchylikeEig( DD(pnd(i)),F,D(lr(i,1)),U(lr(i,1)),V(lr(i,1)),&
               DIFL,DIFR,PL(lr(i,1)),PU(lr(i,1)),mi,mi,lflops )
          PH%D(i) = mi
          PH%pd(1,i) = pnd(i)
          PH%pd(2,i) = pnd(i) + mi*mi -1
       end if
    end do
!$OMP END PARALLEL DO

! *****************************************************
!                 Traverse  upwards                   *
!******************************************************
    DO LL = lvl, 1, -1     ! different level

! ******************* Row Compression ************************! 
       n1 = TRE%lentr(LL)
!$OMP PARALLEL PRIVATE(nodi,mi,ni,nt,Rk), FirstPrivate(indx,Fi,Vi,PUi,Di) 
!$OMP DO SCHEDULE(dynamic)
       DO i = 1, n1 
          nodi = TRE%ttr( LL, i )
          IF( ch( 1,nodi ) .ne. 0 ) THEN   ! parent node
             lr(nodi,1) = lr( ch(1,nodi),1 )
             lr(nodi,2) = lr( ch(2,nodi),2 )
             lc(nodi,1) = lc( ch(1,nodi),1 )
             lc(nodi,2) = lc( ch(2,nodi),2 )
          END IF

          mi = lr(nodi,2)-lr(nodi,1)+1
          ni = lc(nodi,2)-lc(nodi,1)+1
          nt = lca - ni
          if( lc(nodi,1) .eq. 1 ) then
             indx( 1:nt ) = (/ (j, j=lc(nodi,2)+1, lca ) /)
          else
             indx( 1:lc(nodi,1)-1 ) = (/ (j, j=1, lc(nodi,1)-1 ) /)
             indx( lc(nodi,1): nt ) = (/ (j, j=lc(nodi,2)+1, lca ) /)
          end if

          PUi(1:nt) = PU( indx(1:nt) )
          Vi(1:nt)  = V( indx(1:nt) )
          Di(1:nt)  = D( indx(1:nt) )
          Fi(1:mi)  = F( PL(lr(nodi,1):lr(nodi,2)) )
          call ComprcauchyEig( 'r',Fi,Di,U(lr(nodi,1)),Vi,DIFL,DIFR,F,PL(lr(nodi,1)), & 
               PUi,tol,mi,nt,PH,pnh,H,nodi,Rk,lflops,lswap )
       END DO ! (i)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! ****** reorder the remained row generators, U, PL ******
           lt = 1
           DO i = 1, n1
              nodi = TRE%ttr( LL, i )
              Rk = PH%U( 2, nodi )
              Ui( lt:lt+Rk-1 ) = U( lr(nodi,1): lr(nodi,1)+Rk-1 )
              PLi( lt:lt+Rk-1 ) = PL( lr(nodi,1): lr(nodi,1)+Rk-1 )
              lr( nodi,1 ) = lt
              lr( nodi,2 ) = lt + Rk -1
              lt = lt + Rk
           END DO

           DO LLt = LL-1, 1, -1       ! leaf nodes above the current level
              n1 = TRE%lentr( LLt )
              DO i = 1, n1
                 nodi = TRE%ttr( LLt,i )
                 IF( ch( 1,nodi ) == 0 .AND. ch( 2,nodi ) == 0 ) THEN
                    mi = lr(nodi,2)-lr(nodi,1)+1
                    Ui( lt:lt+mi-1 ) = U( lr(nodi,1):lr(nodi,2) )
                    PLi( lt:lt+mi-1 ) = PL( lr(nodi,1):lr(nodi,2) )
                    lr( nodi,1 ) = lt
                    lr( nodi,2 ) = lt + mi -1
                    lt = lt + mi
                 END IF
              END DO
           END DO
           U( 1:lt-1 ) = Ui( 1:lt-1 )
           PL( 1:lt-1 )= PLi( 1:lt-1 )
           LDA = lt-1

! ******************* column compression  ********************
           n1 = TRE%lentr( LL )
!$OMP PARALLEL PRIVATE(nodi,mi,ni,nt,Rk), FirstPrivate(indx,Ui,PLi,Fi)
!$OMP DO SCHEDULE(DYNAMIC)
           DO i = 1, n1
              nodi = TRE%ttr( LL,i )
              mi = lc(nodi,2)-lc(nodi,1)+1
              ni = lr(nodi,2)-lr(nodi,1)+1
              nt = lda - ni    
              if( lr(nodi,1) .eq. 1 ) then
                 indx( 1:nt ) = (/ (j, j=lr(nodi,2)+1, lda ) /)
              else
                 indx( 1:lr(nodi,1)-1 ) = (/ (j, j=1, lr(nodi,1)-1 ) /)
                 indx( lr(nodi,1): nt ) = (/ (j, j=lr(nodi,2)+1, lda ) /)
              end if

              Ui(1:nt)  = U( indx(1:nt) )
              PLi(1:nt) = PL( indx(1:nt) )
              Fi(1:nt)  = F( PL(indx(1:nt)) )
              call ComprcauchyEig('c',Fi,D(lc(nodi,1)),Ui,V(lc(nodi,1)),DIFL,&
                   DIFR,F,PLi,PU(lc(nodi,1)),tol,nt,mi,PH,pnh,H,nodi,Rk,lflops,lswap )
           END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

 ! ****** reorder the remained col generators, V, D, PUi ******
           lt = 1
           n1 = TRE%lentr( LL )
           DO i = 1, n1
              nodi = TRE%ttr( LL, i )
              Rk = PH%V(1, nodi) 
              Vi(lt : lt+Rk-1 ) = V( lc(nodi,1): lc(nodi,1)+Rk-1 )
              Di(lt : lt+Rk-1 ) = D( lc(nodi,1): lc(nodi,1)+Rk-1 )
              PUi(lt : lt+Rk-1 ) = PU( lc(nodi,1): lc(nodi,1)+Rk-1 )
              lc(nodi,1 ) = lt
              lc(nodi,2 ) = lt + Rk -1
              lt = lt + Rk
           END DO

           DO LLt = LL-1, 1, -1       ! leaf nodes above the bottom level
              n1 = TRE%lentr( LLt )
              DO i = 1, n1
                 nodi = TRE%ttr( LLt,i )
                 IF( ch( 1,nodi ) == 0 .AND. ch( 2,nodi ) == 0 ) THEN
                    mi = lc(nodi,2)-lc(nodi,1)+1
                    Vi( lt :lt+mi-1 ) = V( lc(nodi,1):lc(nodi,2) )
                    Di( lt :lt+mi-1 ) = D( lc(nodi,1):lc(nodi,2) )
                    PUi( lt :lt+mi-1 ) = PU( lc(nodi,1):lc(nodi,2) )
                    lc( nodi,1 ) = lt
                    lc( nodi,2 ) = lt + mi -1
                    lt = lt + mi
                 END IF
              END DO
           END DO
           V( 1:lt-1 ) = Vi( 1:lt-1 )
           D( 1:lt-1 ) = Di( 1:lt-1 )
           PU( 1:lt-1 ) = PUi( 1:lt-1 )
           LCA = lt-1

 ! ****** Store the generator B for nodes at the same level ******
           n1 = TRE%lentr( LL )
           DO i = 1, n1
              nodi = TRE%ttr( LL,i )
              mi  = lr(nodi,2)-lr(nodi,1) + 1

              if( nodi .eq. TRE%ch(1,tr(nodi) ) ) then 
                 nodj = TRE%ch(2, tr(nodi) )    ! adjacent of nodi
              else
                 nodj = TRE%ch(1, tr(nodi) )
              end if
              ni  = lc(nodj,2)-lc(nodj,1) + 1
              call CauchylikeEig(H(pnh),F,D(lc(nodj,1)),U(lr(nodi,1)),V(lc(nodj,1)),DIFL,&
                   DIFR,PL(lr(nodi,1)),PU(lc(nodj,1)),mi,ni,lflops )
              call hssexpmm2( 'B',PH,pnh,mi,ni,nodi,iinfo )  ! B{nodi}
           END DO

        END DO ! (LL)

    DEALLOCATE( ch,lr,lc,Di,Fi,Ui,Vi,PLi,PUi,D )

  END SUBROUTINE Cauchy2hssEig_omp

!!!!!!
  SUBROUTINE ComprcauchyEig(rowcol,F,D,U,V,DIFL,DIFR,FF,PL,PU,tol,&
       M,N,PH,pnh,H,nodi,Rk,nflops,nswap)
    USE aux_hss
!
! ..Scalar Arguments..
    INTEGER, INTENT(IN)    :: M,N,nodi
    INTEGER, INTENT(OUT)   :: Rk,nflops,nswap
    INTEGER, INTENT(INOUT) :: pnh
    DOUBLE PRECISION, INTENT(IN)  :: tol
    CHARACTER(LEN=1), INTENT(IN)  :: rowcol
! ..
! ..Array Arguments..
    DOUBLE PRECISION, INTENT(IN)    :: DIFL(*),FF(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: D(*),F(*),U(*),V(*),H(*)
    INTEGER, INTENT(INOUT)     :: PL(M),PU(N)
    TYPE(HSSMM), INTENT(INOUT) :: PH
!
! Purpose
! ========
! This routine is written for computing the eigenvector matrix of a diagonal matrix with rank-1 modification. 
! In this routine C is an M-by-N matrix, and defined as F*C-C*D= u v^T. 
! For row compression, C is factorized apprximately as 
! C = PQ*[I; Mat_N]* C(PQ(1:Rk),:), where PQ is a permutation and Mat_N is also a Cauchy-like matrix sastifying
! F(Rk+1:M) * Mat_N - Mat_N*F(1:Rk) = Z*W. It is computed via rrluCauchyEigrow
!
! For column compression, C is approximated as C = C(:,PQ(1:Rk) ) * [I Mat_N] * PQ, where
! Mat_N satisfies D2 *Mat_N- Mat_N * D1 = W(1:Rk)*Z(Rk+1:N) and D2=D(Rk+1:N) and D1=D(1:Rk). 
! 
! Parameters
! ==========
! rowcol (in) CHARACTER ( LEN = 1 )
!        = 'R', row compression
!        = 'C', column compresion
! 
! F  (in) DOUBLE PRECISION array, DIMENSION( M )
!    The entries of F are the old eigenvalues and row generators of A, and corresponds to PL.
!
! D  (in)  DOUBLE PRECISION array, DIMENSION( N )
!    The updated eigenvalues, column generators, and corresponding to PU. 
!
! U  (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!    U is the updated z in our old notation and it will be perumuted, only
!    the first Rk entries are useful for row compression. 
! 
! V  (inout)  DOUBLE PRECISION array, DIMENSION( N )
!    V is the column normalization scalars, Alpha. It will be permuted and
!    the first Rk entries are useful for column compression and V will not be used in 
!    row compression.
!
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!      DIFL(i) = D0(i)-FF(i), the order of which will not change and its elements
!      are also referenced by using PL or PU. Positive
! 
! DIFR (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!      DIFL(i) = D0(i)-FF(i+1), the order of which will not change and its elements
!      are also referenced by using PL or PU. Negative
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( LN ), original F
!    LN is the size of original problem. 
!
! PL (inout) INTEGER array, DIMENSION( M )
!    It stores the permutations for row generators of A. It corresponds to F for row compression and
!    it corresponds to F.
!
! PU (inout) INTEGER array, DIMENSION( N )
!    It stores the permutations for col generators of A. It corresponds to D for column compression and
!    it corresponds to D.
! 
! TOL (in) DOUBLE PRECISION, parameter for low rank approximation
!
! M  (in) INTEGER, row dimension of A
!     Now it only works when M == N.
!
! N (in) INTEGER, col dimension of A
!
! PH (inout) HSSMM type, stores the information of U,V and B
!
! pnh (inout) INTEGER, the end of array H
!
! H (inout) DOUBLE PRECISION array, DIMENSION( 7*LN*K )
!   It stores U,V and B.
!
! nodi (in) INTEGER, node i
!
! Rk   (out) INTEGER, the computed numerical rank
!
! nflops (out) INTEGER, the total flops of constructing this HSS matrix 
!
! nswap  (out) INTEGER, the total swaps whiling constructing this HSS matrix
!
! =========
! Written by S.-G. Li, on April 17th, 2013
! for tridiagonal eigenvalue problem
! ========================================
!
! .. Parameters ..
    DOUBLE PRECISION ZERO, ONE
    PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0)
! ..    
! .. Local scalars ..
    INTEGER          mn,i,info,lflops, ierr
    LOGICAL          CR
! ..    
! .. Local arrays ..
    INTEGER, ALLOCATABLE :: PQ(:)
    DOUBLE PRECISION, ALLOCATABLE :: Q(:,:),Z(:),W(:)
!   Q(:,:): the interpolative matrix, Z^{(k)}
!   Z(:)  : the row generators of Z^{(k)}
!   W(:)  : the column generators of Z^{(k)}
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MIN
! ..
! .. External Functions ..
    LOGICAL     LSAME
    EXTERNAL    LSAME
      
    MN = MIN( M,N )      
    nflops = 0
    CR = LSAME( rowcol,'r' )
    IF( CR )  THEN                                   ! block row
       ALLOCATE( Z(M),W(MN),PQ(M), stat=ierr )
       Z = ZERO
       W = ZERO
       call rrluCauchyEigrow( F,D,U,V,DIFL,DIFR,tol,Z,W,FF,PL,PU,M,N,Rk,PQ,lflops,nswap )
       nflops = nflops + lflops
       call invp(PQ,M)  ! PL --> InvPL
       IF( Rk .LT. MN ) THEN
          ALLOCATE( Q(M,Rk) )
          call CauchylikeEigQ( Q,F(Rk+1),F,Z(Rk+1),W,M,Rk,lflops )
          nflops = nflops + lflops
          Q(1:M,1:Rk) = Q(PQ,1:Rk)
!$OMP CRITICAL (Row2)
          call dlacpy('A',M,Rk,Q,M,H(pnh),M)      ! copy U to H
          call hssexpmm2('U',PH,pnh,M,Rk,nodi,info) ! copy Q to generators
!$OMP END CRITICAL (Row2)
       ELSE
          ! copy identity matrix to generators
          ALLOCATE( Q(M,M) )
          Rk = M
          Q(1:M,1:M) = ZERO
          DO i = 1,M
             Q(i,i) = ONE
          END DO
          Q(1:M,1:M) = Q(PQ,1:M)            
!$OMP CRITICAL (Row2)
          call dlacpy('A',M,M,Q,M,H(pnh),M)         ! copy U to H
          call hssexpmm2('U',PH,pnh,Rk,M,nodi,info) ! copy Q to generators
!$OMP END CRITICAL (Row2)
       END IF

    ELSE ! block column

       ALLOCATE( Z(N),W(MN),PQ(N), stat = ierr )    ! Check ??
       Z = ZERO
       W = ZERO
       call rrluCauchyEigcol( D,F,V,U,DIFL,DIFR,tol,Z,W,FF,PU,PL,N,M,Rk,PQ,lflops,nswap )
       nflops = nflops + lflops
       call invp(PQ,N)       ! PL --> InvPL
       IF( Rk .LT. MN ) THEN
          ALLOCATE( Q(Rk, N) )
          call CauchylikeEigQ2(Q,D,D(Rk+1),W,Z(Rk+1),Rk,N,FF,DIFL,DIFR,PU,lflops)
          nflops = nflops + lflops
          Q(1:Rk,1:N) = Q(1:Rk,PQ)            
!$OMP CRITICAL (Column2)
          call dlacpy('A',Rk,N,Q,Rk,H(pnh),Rk)       ! copy V to H
          call hssexpmm2('V',PH,pnh,Rk,N,nodi,info)  ! copy Q to generators
!$OMP END CRITICAL (Column2)
       ELSE
          ! copy identity matrix to generators
          allocate( Q(N,N) )
          Rk = N
          Q(1:N,1:N) = ZERO
          DO i = 1,N
             Q(i,i) = ONE
          END DO
          Q(1:N,1:N) = Q(1:N,PQ)
!$OMP CRITICAL (Column2)
          call dlacpy('A',N,Rk,Q,N,H(pnh),N)        ! copy V to H
          call hssexpmm2('V',PH,pnh,N,Rk,nodi,info) ! copy Q to generators
!$OMP END CRITICAL (Column2)
       END IF

    END IF ! compr type

    DEALLOCATE(Z,W,Q,PQ )

  END SUBROUTINE ComprcauchyEig


!!!!!!!!
  SUBROUTINE DHSSEVC( K,Ni,F,D,U,V,DIFL,DIFR,Q,LDQ,INFO )
     USE aux_hss
     USE ConstructHssd
     USE BasicMM
!
! ..Scalar Parameters..
    INTEGER, INTENT(IN) ::  K, LDQ, INFO, Ni
!
! ..Array Parameters..
    DOUBLE PRECISION, INTENT(IN) :: D(*),F(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: U(*),V(*),Q(LDQ,*)
!
! Purpose
! ========
! Firstly, this routine constructs an HSS approximation H to an orthogonal Cauchy-like 
! matrix, A(i,j) = U(i)*V(j) / (F(i)-D(j)), where F contains the old eigenvalues, and
! D contains the updated eigenvalues in ascending order. A is a K-by-K matrix. 
! Secondly, the first K rows of Q is multiplied by H from left via HSS matrix multiplication.
!
! There is no enough reason for choosing NCOL = 200
!
! Arguments  
! =========
! K    (in) INTEGER 
!      The size of secular equation, the dimension of that diagonal matrix
!
! Ni   (in) INTEGER
!      The block size of nodes at the bottom level
! 
! F    (in) DOUBLE PRECISION array, DIMENSION( K )
!      The old eigenvalues of previous tridiagonal matrices, row generators
!
! D    (in) DOUBLE PRECISION array, DIMENSION( K )
!      The updated eigenvalues of current level, column generators
!
! U    (inout) DOUBLE PRECISION array, DIMENSION( K )
!      The updated vector z, row generators
! 
! V    (inout) DOUBLE PRECISION array, DIMENSION( K )
!      The reciprocal of the 2-norm of column of A, column generators
!
! DIFL (in) DOUBLE PRECISION array, DIMENSION( K )
!      The difference between old eigenvalues and new eigenvalues
!      DIFL(J) = D(J) - F(J), positive value.
!
! DIFL (in) DOUBLE PRECISION array, DIMENSION( K )
!      The difference between old eigenvalues and new eigenvalues
!      DIFR(J) = D(J) - F(J+1), negative value.
!
! Q    (inout) DOUBLE PRECISION array, DIMENSION( K, LDQ )
!      The first K rows of right singular vectors that are needed
!      to be updated. 
! 
! LDQ (in) INTEGER
!      The row and column dimension of eigenvector matrix of this level
! 
! DD   (workspace) DOUBLE PRECISION array
!      Stores the diagonal blocks of HSS matrix, its size is about K*Ni
! 
! H    (workspace) DOUBLE PRECISION array
!      Stores the other generators of HSS matrix, U, V and B, its size is
!      about 7*K*Ni
! 
! WORK (workspace) DOUBLE PRECISION array, with size at least K*200
!      It is used in HSS matrix multiplication. 
! 
! INFO (output) INTEGER
!      = 0: successful exit.
!      otherwise: something is wrong.
! 
! Further Details
! ===============
! Written by S.-G. Li, on April 17th, 2013, in Changsha China.
! =============================================================
! 
! .. Parameters ..
    INTEGER           NCOL
    PARAMETER         ( NCOL = 200 )
    DOUBLE PRECISION  ONE, ZERO, TOL
    PARAMETER         ( ONE= 1.0D+0, ZERO= 0.0E+0, TOL= 1.0E-17 )
! ..
! .. Local Scalars ..
    INTEGER  i, nflops, nswap, nc, nlast,ierr,tnk, mr
    DOUBLE PRECISION  Maxerr
! ..
! ..Parameters for HSS tree..
    INTEGER              :: ltr,lm,n1,lvl,pnh
    type(HTinfo)         :: TRE
    type(hssmm)          :: PH
    INTEGER, ALLOCATABLE :: tr(:), M(:)
    DOUBLE PRECISION, ALLOCATABLE :: DD(:), HH(:)

! *********************************************************
!             Construct HSS tree                          *
! *********************************************************
    ltr   = 0
    mr = mod(K,Ni)
    IF( mr >= Ni/2 ) THEN
       lm = K/ni+1
    ELSE
       lm = K/Ni
    END IF
    ALLOCATE ( M(lm),tr( 2*lm-1 ), stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssevc'
    END IF

    call npart( K,Ni,tr,ltr,M,lm )

!   Tree information
    lvl = ceiling( log2(ltr+1) ) -1
    n1 = 0.5* (ltr +1 )

    ALLOCATE( TRE%ttr(lvl,n1),TRE%lentr(lvl),TRE%lfx(ltr),TRE%ch(2,ltr), stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssevc'
    END IF

    call GetTreeInfo( TR,LTR,M,M,TRE,lvl )

    ALLOCATE( PH%D(ltr),PH%U(2,ltr),PH%V(2,ltr),PH%B(2,ltr),PH%ph(4,ltr),PH%pd(2,ltr),stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssevc'
    END IF

    call hssexpmm0( PH, ltr )

    tnk = 0
    DO i = 1, lm
       tnk = tnk+ M(i)**2
    END DO
    ALLOCATE( HH( 10*Ni*K ),DD(tnk), stat=ierr ) 
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in DHSSEVC'
    END IF

    pnh = 1
    nswap = 0
    nflops = 0
    call cauchy2hssEig_omp( F,D,U,V,DIFL,DIFR,K,TR,LTR,M,LM,PH,HH,DD,TOL,TRE,lvl,pnh,nflops,nswap )
!
    NC = FLOOR( DBLE(LDQ) /NCOL )
    DO i = 1, NC
       call fastHssmmL_omp( DD,HH,PH,M,Q((i-1)*NCOL+1:i*NCOL,1:K),tr,ltr,NCOL,K,TRE,lvl )
    END DO
    nlast = LDQ - NC * NCOL
    IF( nlast .NE. 0 ) THEN
       call fastHssmmL_omp( DD,HH,PH,M,Q( NC*NCOL+1:LDQ,1:K ),tr,ltr,nlast,K,TRE,lvl )
    END IF

    DEALLOCATE( TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch )
    DEALLOCATE( PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd )
    DEALLOCATE( M, tr ) 

  END SUBROUTINE DHSSEVC


end module CauchyHssEig

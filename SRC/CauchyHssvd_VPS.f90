module CauchyHssvd_VPS
  use ConstructHssd
!
  implicit none

contains

!!!!!!!!!
    SUBROUTINE rrluCauchysvdrow( D,F,A,V,DIFL,DIFR,tol,U,W,FF,PL,PU,M,N,Rk,PQ,nflops,nswap )
!
! .. Scalar parameter ..
    INTEGER, INTENT(IN)          :: M,N
    INTEGER, INTENT(OUT)         :: Rk,nflops,nswap
    DOUBLE PRECISION, INTENT(IN) :: tol
! .. Array parameters ..
    DOUBLE PRECISION, INTENT(IN)    :: DIFL(*),DIFR(*),FF(*)
    DOUBLE PRECISION, INTENT(INOUT) :: A(*),V(*),U(*),W(*),D(*),F(*)
    INTEGER, INTENT(INOUT)          :: PL(*),PU(*),PQ(*)
!
! Purpose
! =========
! This routine computes a Rank Revealing Schur Complement, RRSC, for a Cauchy-like matrix
! via structured matrix techniques. This matrix C has dimensions M-by-N, with generators D,F,A and V.
! C(i,j) = A(i)*V(j) / ( D(i)**2 -F(j)**2 ). 
!
! C is factorized as C(PQ,:) = [I; Z] * C( PQ(1:Rk), : ), and Z is a Cauchy matrix, defined as
! D_2 * Z - Z *D_1 = U(Rk+1:M)*W(1:Rk) with D_2 = D(Rk+1:M)**2 and D_1 = D(1:Rk)**2. 
! Z is an (M-Rk+1)-by-Rk matrix. Note that since D_1 and D_2 are the updated singular values,
! their differences can not be computed directly. D_1(i) may equal to D_2(j), but mathematically they 
! are different. 
!
! ..Parameters..
! ===============
! D  (inout)  DOULBE PRECISION Array, DIMENSION( N ) 
!    Row generators, the updated singular values and will be modified. F(i) < D(i) < F(i+1).
! 
! F  (inout)  DOUBLE PRECISION Array, DIMENSION( M ). 
!    Column generators, the old singular values, which will be changed but are useless for 
!    row compression.
!
! A  (inout)  DOUBLE PRECISION Array, DIMENSION( M ). 
!    Row generators. It may the updated Z or the reciprocal of 2-norm of each row (col) of the original Cauchy-like matrix.
!   
! V (inout)  DOUBLE PRECISION Array, DIMENSION( N ). 
!    Col generators. It may the updated Z or the reciprocal of 2-norm of each row (col) of the original Cauchy-like matrix.
! 
! tol (in) DOUBLE PRECISION, for low rank approximation
!
! DIFL (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the left distance, DIFL(i) = D0(i)-F0(i), positive values. The entries of DIFL are 
!      referred by PL. 
!
! DIFR (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the right distance, DIFR(i) = D(i)-F(i+1), negative values. The entries of DIFR are 
!      referred by PU. 
! 
! U (out) DOUBLE PRECISION Array, DIMENSION(M) 
!    Row generators of Mat_N which is a matrix in RRSC. The last M-Rk elements of U are useful.
!
! W (out) DOUBLE PRECISION Array, DIMENSION( MIN( M,N ) )
!    Row generators of Mat_N and the first Rk elements of W are useful.
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
! Written by S.-G. Li, On Dec. 15th, 2012
! 
! Modified on April 9th, 2013
! 1) Remove DD, the original updated singular values, and use only D.
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
    rho = 1.1E+0
    prd = 10
    nflops = 0
    
    DO k = 1, mn
       call CauchyPivtGsvdm(D,F,U,V,PL,PU,k,A,M,N,PQ,lflops,FF,DIFL,DIFR)
       nflops = nflops + lflops

       temp = NEGONE* computdw( FF,k,k,PL,PU,DIFL,DIFR )  ! D(k)-F(k)
       Ukk = U(k)*V(k) / (D(k)+F(k)) / temp
       Amax = MAX( Amax,ABS(Ukk) )
       IF( ABS(Ukk) .LT. tol*Amax ) THEN  ! first step converged
          ! Complete pivoting for Schur complement
          call CauchyPivtGsvdm_CP(D,F,U,V,PL,PU,k,A,M,N,PQ,lflops,FF,DIFL,DIFR)
          temp = NEGONE* computdw( FF,k,k,PL,PU,DIFL,DIFR )  ! D(k)-F(k)
          Ukk = U(k)*V(k) / (D(k)+F(k)) / temp
          Amax = MAX( Amax,ABS(Ukk) )
          IF( ABS(Ukk) .LT. tol*Amax ) THEN  ! final converged
             Rk = k-1
             EXIT
          END IF
       END IF
       
       DO j = k+1, M    !  U(k+1:M) = U(k+1:M)* (D(k+1:M)-D(k)) / (D(k+1:M)-F(k))

          temp  = computww( FF,j,k,PL,DIFL,DIFR )            ! D(j)-D(k)
          temp2 = NEGONE* computdw( FF,k,j,PL,PU,DIFL,DIFR ) ! D(j)-F(k)
          temp =  temp/ temp2

          U(j) = U(j) * ( (D(j)+D(k)) / (D(j)+F(k)) ) * Temp
       END DO

       DO j = k+1, N   !  V(k+1:N) = V(k+1:N)* (F(k+1:N)-F(k)) / (F(k+1:N)-D(k))

          temp = computdw( FF,j,k,PL,PU,DIFL,DIFR ) ! F(j)-D(k)
          V(j) = V(j) *( (F(j)-F(k)) /temp )* ( ( F(j)+F(k) ) / (F(j)+D(k)) )
       END DO
              
       DO j = 1, k-1   ! W(j) = W(j) * ( (F(k)-D(j))/ (D(k)-D(j)) )

          temp  = computdw( FF,k,j,PL,PU,DIFL,DIFR )  ! F(k)-D(j)
          temp2 = computww( FF,k,j,PL,DIFL,DIFR )     ! D(k)-D(j)
          temp  = temp / temp2
          W(j) = W(j) *( (F(k)+D(j))/ (D(k)+D(j)) ) * Temp
       END DO

       temp = NEGONE* computdw( FF,k,k,PL,PU,DIFL,DIFR )  !  W(k) = (D(k)-F(k))/A(k)
       W(k) = (D(k)+F(k)) * temp / A(k)

       DO j = 1, k-1   ! W(k) = W(k) * (F(j)-D(k)) / ( D(j)-D(k) )
          temp  = computdw( FF,j,k,PL,PU,DIFL,DIFR )  !  F(j)-D(k)
          temp2 = computww( FF,j,k,PL,DIFL,DIFR )     !  D(j)-D(k)
          temp  = temp / temp2
          W(k) = W(k) * ( (F(j)+D(k))/ (D(j)+D(k)) ) * Temp
       END DO

       nflops = nflops + 6*(M-k)+5*(N-k)+2*7*(k-1)

       ! swap
       flgl = mod(k,prd)
       DO while(flgl .lt. 1)
          flgl = 1
          call searchMax2svdm(U(k+1),W,D,FF,PL,DIFL,DIFR,M-k,k,ii,jj,Nmax,lflops)
          nflops = nflops + lflops

          IF(Nmax .gt. rho) THEN
!             write(*,*) 'Nmax=', Nmax, 'swap once'
             nswap = nswap + 1
             flgl = 0
             jj = jj + k 

             DO j = k+1, N   !V(k+1:N)    = V(k+1:N) * ( (F(k+1:N)-D(ii)) / (F(k+1:N)-D(jj)) )
                temp  = computdw( FF,j,ii,PL,PU,DIFL,DIFR )  ! F(j)-D(ii)
                temp2 = computdw( FF,j,jj,PL,PU,DIFL,DIFR )  ! F(j)-D(jj)
                temp  = temp / temp2
                V(j)  = V(j)* ( (F(j)+D(ii) ) / (F(j)+D(jj)) ) * temp
             END DO

             DO j = k+1, jj-1  !U(k+1:jj-1) = U(k+1:jj-1) * ( (D(k+1:jj-1)-D(jj)) / (D(k+1:jj-1)-D(ii)) )
                temp  = computww( FF,j,jj,PL,DIFL,DIFR )  !  D(j)-D(jj)
                temp2 = computww( FF,j,ii,PL,DIFL,DIFR )  !  D(j)-D(ii)
                temp  = temp / temp2
                U(j)  = U(j)* ( (D(j)+D(jj))/ (D(j)+D(ii)) ) * temp
             END DO
             DO j = jj+1, M  !U(jj+1:M)   = U(jj+1:M) * ( (D(jj+1:M)-D(jj)) / (D(jj+1:M)-D(ii)) )
                temp  = computww( FF,j,jj,PL,DIFL,DIFR )  !  D(j)-D(jj)
                temp2 = computww( FF,j,ii,PL,DIFL,DIFR )  !  D(j)-D(ii)
                temp  = temp / temp2
                U(j)  = U(j)* ( (D(j)+D(jj))/ (D(j)+D(ii)) ) * temp
             END DO

             DO j = 1, ii-1  ! W(1:ii-1)   = W(1:ii-1) * ( (D(1:ii-1)-D(ii)) / (D(1:ii-1)-D(jj)) )
                temp  = computww( FF,j,ii,PL,DIFL,DIFR )  !  D(j)-D(ii)
                temp2 = computww( FF,j,jj,PL,DIFL,DIFR )  !  D(j)-D(jj)
                temp  = temp / temp2
                W(j)  = W(j)* ( (D(j)+D(ii))/ (D(j)+D(jj)) ) * temp
             END DO
             DO j = ii+1, k  ! W(ii+1:k)   = W(ii+1:k) * ( (D(ii+1:k)-D(ii)) / (D(ii+1:k)-D(jj)) )
                temp  = computww( FF,j,ii,PL,DIFL,DIFR )  !  D(j)-D(ii)
                temp2 = computww( FF,j,jj,PL,DIFL,DIFR )  !  D(j)-D(jj)
                temp  = temp / temp2
                W(j)  = W(j)* ( (D(j)+D(ii))/ (D(j)+D(jj)) ) * temp
             END DO

             ! U(jj)       = A(ii) * ( (D(ii)-D(jj)) / (D(ii)-F(ii)) )  
             temp  = computww( FF,ii,jj,PL,DIFL,DIFR )    ! D(ii)-D(jj)
             temp2 = NEGONE* computdw( FF,ii,ii,PL,PU,DIFL,DIFR )  ! D(ii)-F(ii)
             temp = temp/ temp2
             U(jj) = A(ii) * ( (D(ii)+D(jj))/ (D(ii)+F(ii)) ) * temp

             ! W(ii)       = (D(jj)-F(ii)) / A(jj)
             temp  = NEGONE* computdw( FF,ii,jj,PL,PU,DIFL,DIFR )  ! D(jj)-F(ii)
             W(ii) = ( (D(jj)+F(ii))/ A(jj) ) * temp 

             DO j = 1, ii-1
                ! U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                temp  = computww( FF,ii,j,PL,DIFL,DIFR )    ! D(ii)-D(j)
                temp2 = NEGONE* computdw( FF,j,ii,PL,PU,DIFL,DIFR )  ! D(ii)-F(j)
                temp  = temp / temp2
                U(jj) = U(jj) * ( (D(ii)+D(j)) / (D(ii)+F(j)) ) * temp

                ! W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
                temp  = computdw( FF,j,jj,PL,PU,DIFL,DIFR )  ! F(j)-D(jj)
                temp2 = computww( FF,j,jj,PL,DIFL,DIFR )     ! D(j)-D(jj)
                temp  = temp / temp2
                W(ii) = W(ii) * ( (F(j)+D(jj))/ (D(j)+D(jj)) ) * temp
             END DO
             DO j = ii+1, k
                ! U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                temp  = computww( FF,ii,j,PL,DIFL,DIFR )    ! D(ii)-D(j)
                temp2 = NEGONE* computdw( FF,j,ii,PL,PU,DIFL,DIFR )  ! D(ii)-F(j)
                temp  = temp / temp2
                U(jj) = U(jj) * ( (D(ii)+D(j)) / (D(ii)+F(j)) ) * temp

                ! W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
                temp  = computdw( FF,j,jj,PL,PU,DIFL,DIFR )  ! F(j)-D(jj)
                temp2 = computww( FF,j,jj,PL,DIFL,DIFR )     ! D(j)-D(jj)
                temp  = temp / temp2
                W(ii) = W(ii) * ( (F(j)+D(jj))/ (D(j)+D(jj)) ) * temp
             END DO
             call iswap(PL,ii,jj)
             call iswap(PQ,ii,jj)
             call dswap(D,ii,jj)
             call dswap(A,ii,jj)

             nflops = nflops + 10*(N-k)+12*M+(11*2)*k
          END IF ! Nmax
          
       END DO ! while

    END DO ! main loop

  END SUBROUTINE rrluCauchysvdrow

  SUBROUTINE rrluCauchysvdcol( F,D,A,V,DIFL,DIFR,tol,U,W,FF,PL,PU,M,N,Rk,PQ,nflops,nswap )
!
! ..Scalar parameter..
    INTEGER, INTENT(IN)          :: M,N
    INTEGER, INTENT(OUT)         :: Rk,nflops,nswap
    DOUBLE PRECISION, INTENT(IN) :: tol
! ..
! ..Array parameters..
    DOUBLE PRECISION, INTENT(IN)    :: DIFL(*),DIFR(*),FF(*)
    DOUBLE PRECISION, INTENT(INOUT) :: A(*),V(*),U(*),W(*),D(*),F(*)
    INTEGER, INTENT(INOUT)          :: PL(*),PU(*),PQ(*)
!
! Purpose
! =========
! This routine computes a Rank Revealing Schur Complement, RRSC, for a Cauchy-like matrix
! via structured matrices. This matrix C has dimensions M-by-N, with generators D,F,A and V.
! C(i,j) = A(i)*V(j) / (D(j)**2 -F(i)**2 ) satisfies F*C - C*D = -A*V.  
! In order to compress HSS block column, we transform it to HSS block row. 
! Here F are the old singular values and D are the updated ones. 
!
! It computes a low rank approximation of C, C(PQ,:) = [I; Z]*C( PQ(1:Rk), :) where Z is also a 
! Cauchy-like matrix and satisfies D(Rk+1:M) * Z - Z * D(1:Rk) = U * W. 
! This routine can be used to compress the off-diagonal HSS block row of left singular vector matrix 
! and the HSS block column of right singular vector matrix which is first transposed to HSS block row form. 
! Different from rrluCauchysvdrow, the row generator of C is F and their difference can be computed directly,
! and therefore this routine is a little simpler. 
!
! Now F is the row generator, and D is the column generator. We transform an HSS block column to an
! HSS block row, and consider it as an HSS block row.
!
! Parameters
! ==========
! F  (inout)  DOULBE PRECISION Array, DIMENSION(M). Row generators, the old singular values.
!    It will be modified and the first Rk elements are chosen ones. 
! 
! D  (inout)  DOUBLE PRECISION array, DIMENSION(N). Col generators, the updated singular values.
!
! A (inout)  DOUBLE PRECISION array, DIMENSION(M). Row generators, updated Z.
!   
! V (inout)  DOUBLE PRECISION array, DIMENSION(M). Col generators, alpha.
! 
! tol (in) DOUBLE PRECISION array
!
! DIFL (INPUT) DOUBLE PRECISION array, DIMENSION(N), the left distance, negative values
!
! DIFR (INPUT) DOUBLE PRECISION array, DIMENSION(N), the right distance, negative values
! 
! U (OUTPUT) DOUBLE PRECISION array, DIMENSION(M). Row generators of Mat_N. 
!   Mat_N is a matrix in RRSC. The last M-Rk elements of U are useful.
!
! W (OUTPUT) DOUBLE PRECISION array, DIMENSION( min(M,N) ). Row generators of Mat_N. 
!   Mat_N is a matrix in RRSC. The first Rk elements of V are useful.
!
! FF (input) double precision array, dimension( LN )
!    The orignal old singular values
!
! PL (OUTPUT) INTEGER array, DIMENSION(M). Row permutation. 
!    It relates with row generators, F and FF
!
! PU (OUTPUT) INTEGER array, DIMENSION(M). Col permutation.
!    It relates with column generatros, D and DD
! 
! M (INPUT) INTEGER, Row dimension of C
!
! N (INPUT) INTEGER, Col dimension of C
! 
! Rk (OUTPUT) INTEGER, rank of returned matrix.
!
! ============
! Written by Sheng-Guo Li, On Sept. 18th, 2012
! Modified on April 10th, 2013
! 1) Exchange D and F. Now D contains the updated singular values, and 
!    F contains old singular values.
! 2) Remove the use of DD, the original updated singular values, and
!    use D directly.
! ===========================================
!
! .. Local Parameters ..
    DOUBLE PRECISION ONE, ZERO, NEGONE
    PARAMETER        ( ZERO= 0.0E+0, ONE=1.0E+0, NEGONE= -1.0E+0 )
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
    rho = 1.1E+0
    prd = 10
    nflops = 0
    
    DO k = 1, mn
       call CauchyPivtGsvd(F,D,U,V,PL,PU,k,A,M,N,PQ,lflops,FF,DIFL,DIFR )
       ! D is new; F is old; U is Z and V is alpha; PL is for D and PU is for W
       nflops = nflops + lflops

       temp = computdw( FF,k,k,PU,PL,DIFL,DIFR )
       Ukk  = u(k)*V(k) / (F(k)+D(k)) / temp
       Amax = MAX( Amax,ABS(Ukk) )
       IF( ABS(Ukk) .LT. tol*Amax ) THEN  ! first step converged
          ! Complete pivoting for Schur complement
          call CauchyPivtGsvd_CP(F,D,U,V,PL,PU,k,A,M,N,PQ,lflops,FF,DIFL,DIFR)
          temp = computdw( FF,k,k,PU,PL,DIFL,DIFR )
          Ukk  = u(k)*V(k) / (D(k)+F(k)) / temp
          Amax = MAX( Amax,ABS(Ukk) )
          IF( ABS(Ukk) .LT. tol*Amax ) THEN  ! final converged
             Rk = k-1
             EXIT
          END IF
       END IF
       
       DO j = k+1, M    !  U(k+1:M) = U(k+1:M)* (D(k+1:M)-D(k)) / (D(k+1:M)-F(k))
          temp = computdw( FF,j,k,PU,PL,DIFL,DIFR )  ! D(j) - F(k)
          Temp = (F(j)-F(k)) / Temp
          U(j) = U(j) * ( (F(j)+F(k)) / (F(j)+D(k)) ) * Temp
       END DO

       DO j = k+1, N   !  V(k+1:N) = V(k+1:N)* (F(k+1:N)-F(k)) / (F(k+1:N)-D(k))
          temp  = computww( FF,j,k,PU,DIFL,DIFR )            ! F(j) - F(k)
          temp2 = NEGONE* computdw( FF,k,j,PU,PL,DIFL,DIFR ) ! F(j) - D(k)
          temp  = temp/ temp2
          V(j) = V(j) * ( (D(j)+D(k)) / (D(j)+F(k)) )* Temp
       END DO
              
       DO j = 1, k-1    ! W(j) = W(j) * ( (F(k)-D(j))/ (D(k)-D(j)) )
          temp = NEGONE* computdw( FF,j,k,PU,PL,DIFL,DIFR )  ! F(k)-D(j)
          W(j) = W(j) *( (D(k)+F(j))/ (F(k)+F(j)) ) * ( Temp / (F(k)-F(j)) )
       END DO
       temp = computdw( FF,k,k,PU,PL,DIFL,DIFR )
       W(k) = ( D(k)+F(k) ) * temp / A(k)

       DO j = 1, k-1   ! W(k) = W(k) * (F(j)-D(k)) / ( D(j)-D(k) )
          temp = NEGONE* computdw( FF,k,j,PU,PL,DIFL,DIFR )   ! F(j)-D(k)
          W(k) = W(k) * ( (D(j)+F(k))/ (F(j)+F(k)) ) * ( Temp/ (F(j)-F(k) ) )
       END DO

       nflops = nflops + 5*(M-k)+7*(N-k)+2*6*(k-1)

       ! swap
       flgl = MOD(k,prd)
       DO WHILE(flgl .lt. 1)
          flgl = 1
          call searchMax2svd(U(k+1),W,F(k+1),F,M-k,k,ii,jj,Nmax,lflops) 
          nflops = nflops + lflops

          IF(Nmax .gt. rho) THEN
!             write(*,*) 'Nmax=', Nmax, 'swap once'
             nswap = nswap + 1
             flgl = 0
             jj = jj + k 

             DO j = k+1, N  ! V(k+1:N)    = V(k+1:N) * ( (F(k+1:N)-D(ii)) / (F(k+1:N)-D(jj)) )
                temp  = NEGONE* computdw( FF,ii,j,PU,PL,DIFL,DIFR )  ! F(j)-D(ii)
                temp2 = NEGONE* computdw( FF,jj,j,PU,PL,DIFL,DIFR )  ! F(j)-D(jj)
                temp  = temp / temp2
                V(j)  = V(j)* ( (D(j)+F(ii) ) / (D(j)+F(jj)) ) * temp                
             END DO

             U(k+1:jj-1) = U(k+1:jj-1) * ( (F(k+1:jj-1)-F(jj)) / (F(k+1:jj-1)-F(ii)) )* &
                            ( (F(k+1:jj-1)+F(jj)) / (F(k+1:jj-1)+F(ii)) )
             U(jj+1:M)   = U(jj+1:M) * ( (F(jj+1:M)-F(jj)) / (F(jj+1:M)-F(ii)) ) * &
                            ( (F(jj+1:M)+F(jj)) / (F(jj+1:M)+F(ii)) )
             W(1:ii-1)   = W(1:ii-1) * ( (F(1:ii-1)-F(ii)) / (F(1:ii-1)-F(jj)) )* & 
                            ( (F(1:ii-1)+F(ii)) / (F(1:ii-1)+F(jj)) )
             W(ii+1:k)   = W(ii+1:k) * ( (F(ii+1:k)-F(ii)) / (F(ii+1:k)-F(jj)) )* & 
                            ( (F(ii+1:k)+F(ii)) / (F(ii+1:k)+F(jj)) )

             ! U(jj)       = A(ii) * ( (D(ii)-D(jj)) / (D(ii)-F(ii)) )  
             temp  = computdw( FF,ii,ii,PU,PL,DIFL,DIFR )  ! D(ii)-F(ii)
             U(jj) = A(ii) * ( (F(ii)+F(jj)) / (D(ii)+F(ii)) )* ( (F(ii)-F(jj)) / temp )

             ! W(ii)       = (D(jj)-F(ii)) / A(jj)
             temp  = computdw( FF,jj,ii,PU,PL,DIFL,DIFR )  ! D(jj)-F(ii)
             W(ii) = temp* (F(jj)+D(ii)) / A(jj)

             DO j = 1, ii-1
                ! U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                temp  = computdw( FF,ii,j,PU,PL,DIFL,DIFR )  ! D(ii)-F(j)
                U(jj) = U(jj) * ( (F(ii)+F(j)) / (F(ii)+D(j)) ) * ( (F(ii)-F(j)) / temp )

                ! W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
                temp  = NEGONE* computdw( FF,jj,j,PU,PL,DIFL,DIFR )  ! F(j)-D(jj)
                W(ii) = W(ii) * ( (D(j)+F(jj)) / (F(j)+F(jj)) ) * ( temp / (F(j)-F(jj)) )
             END DO
             DO j = ii+1, k
                ! U(jj) = U(jj) * ( (D(ii)-D(j)) / (D(ii)-F(j)) )
                temp  = computdw( FF,ii,j,PU,PL,DIFL,DIFR )  ! D(ii)-F(j)
                U(jj) = U(jj) * ( (F(ii)+F(j)) / (F(ii)+D(j)) ) * ( (F(ii)-F(j)) / temp )

                ! W(ii) = W(ii) * ( (F(j)-D(jj)) / (D(j)-D(jj)) )
                temp  = NEGONE* computdw( FF,jj,j,PU,PL,DIFL,DIFR )  ! F(j)-D(jj)
                W(ii) = W(ii) * ( (D(j)+F(jj)) / (F(j)+F(jj)) ) * ( temp / (F(j)-F(jj)) )
             END DO
             call iswap(PL,ii,jj)
             call iswap(PQ,ii,jj)             
             call dswap(F,ii,jj)
             call dswap(A,ii,jj)
!             write(*,*) 'Col swap once ', Nmax

             nflops = nflops + 10*(N-k)+8*M+18*k
          END IF ! Nmax
          
       END DO ! while

    END DO ! main loop

  END SUBROUTINE rrluCauchysvdcol

!!!!!!
  FUNCTION computww( FF,j,k,PL,DIFL,DIFR )
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: j,k
! ..
! .. Array Arguments ..
    INTEGER, INTENT(IN) :: PL(*)
    DOUBLE PRECISION, INTENT(IN) :: FF(*), DIFL(*), DIFR(*)
! 
! Purpose
! =========
! Computes D(j)-D(k) and D is the updated singular values. F ar the old singular 
! values, and D is the updated singular values, D(j) = D0( PL(j) ). 
!
! ..Parameters..
! ==============
! FF (in) DOUBLE PRECISION array, DIMENSION( LN )
!     FF are the old singular values. LN is the size of original problem and FF 
!     is referenced by calling PL. 
!
! j  (in) INTEGER
!
! k  (in) INTEGER
!
! PL (in) INTEGER array, DIMENSION(*)
!      Permutation, corresponding to D
!
! DIFL (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the left distance, DIFL(i) = D0(i)-FF(i), positive values
!
! DIFR (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the right distance, DIFR(i) = D0(i)-FF(i+1), negative values
!
! =============
! Modified by S.-G. Li, on Dec. 24th, 2012
! Modified on April 12th, 2013
! 1) Change DD to FF to be consistent with the other routines.
! ==========================================================
!
    INTEGER          :: PWJJ,PWKK
    DOUBLE PRECISION :: computww

    PWJJ = PL(j)
    PWKK = PL(k)
    IF( PWJJ .lt. PWKK ) THEN
       computww = (FF(PWJJ+1)-FF(PWKK) )-DIFL(PWKK) + DIFR(PWJJ)  ! j < k
    ELSE
       computww = (FF(PWJJ) -FF(PWKK+1) )-DIFR(PWKK)+DIFL(PWJJ)   ! j > k
    END IF

  END FUNCTION computww

!!!!!!
  FUNCTION computdw( FF,j,k,PL,PU,DIFL,DIFR) 
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: j, k
! .. Array Arguments ..
    INTEGER, INTENT(IN) :: PL(*), PU(*)
    DOUBLE PRECISION, INTENT(IN) :: FF(*),DIFL(*),DIFR(*)
!
! Purpose
! ========
! This routine computes F(j)-D(k), where F contains the old singular values and 
! D contains the new singular values. This is a PURE function. 
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
!     Permutation
! 
! PU  (in) INTEGER array, DIMENSION(*)
!     Permutation
!
! DIFL (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the left distance, DIFL(i) = D0(i)-FF(i), positive values
!
! DIFR (in) DOUBLE PRECISION Array, DIMENSION( LN )
!      the right distance, DIFR(i) = D0(i)-FF(i+1), negative values
!
! =============
! Modified by S.-G. Li, on Dec. 24th, 2012
! Modified on April 12th, 2013
! 1) Change DD to FF to be consistent with other routines.
! =======================================================

    INTEGER PDJJ,PWKK
    DOUBLE PRECISION :: computdw

    PDJJ = PU(j)
    PWKK = PL(k)
    IF( PDJJ .LE. PWKK ) THEN
       computdw = ( FF(PDJJ)-FF(PWKK) )- DIFL(PWKK)   ! j < k
    ELSE
       computdw = ( FF(PDJJ)-FF(PWKK+1) )-DIFR(PWKK)  ! j > k
    END IF

  END FUNCTION computdw

!!!!!!
  SUBROUTINE CauchyPivtGsvd( F,D,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR )
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
! where S is the k-th Schur complement and S(i,j) = U(i)*V(j) / (F(i)**2 - D(j)**2 ).
! S is an M-by-N matrix. This routine is designed for rrluCauchysvdcol. 
! Here F are the old singular values and D are the updated singular values. 
! 
! As in rrluCauchysvdcol, the row generators is F and the column generators is D.
! S is different with the HSS block in rrluCauchysvdcol with a minus sign.
! Since we are choosing the largest entry of magnitude, the sign dose not matter.
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
!   The row generator of S, which are the old singular values, and can be used
!   directly;
! 
! D (inout) DOUBLE PRECISION array, DIMENSION( N ) 
!   The column generator of S, the updated singular values;
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
!    The original generator of U
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
!    DIFL(i) = D0(i) - FF(i), positive
! 
! DIFR (in) DOUBLE PRECISION, DIMENSION( NN )
!    DIFR(i) = D0(i) - FF(i+1), negative
!
! ==========
! Modified by S.-G. Li, in Dec. 24th, 2012
! Modified on April 10th, 2013
! 1) Exchange D and F, and remove the use of DD, use D directly.
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
    call CauchyMaxsvdmc( D(k),FF,PU(k),PL(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! largest entry of the first column
    call CauchyMaxsvdm ( D(k),FF,PU(k),PL(k),DIFL,DIFR,V( k ),N-k+1,junkU,jjU )  ! largest entry of the first row
    junkL = junkL * ABS( V(k) )
    junkU = junkU * ABS( U(k) )

    ! Compute the first entry of the K-th Schur Complement
    ! Because of off-diagonal blocks, PL(k) can not equal to PU(k)
    IF( PL(k) .LT. PU(K) ) THEN
       Piv = ABS( U(k)*V(k) / ( ( F(k)-FF(PU(k)) )-DIFL(PU(k)) ) / ( D(k)+F(k) ) )
    ELSE
       Piv = ABS( U(k)*V(k) / ( ( FF(PU(k)+1)-F(k) )+DIFR(PU(k)) ) / ( D(k)+F(k) ) )
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
          call dswap(F,k,jj+k-1)
          call dswap(U,k,jj+k-1)
          call dswap(A,k,jj+k-1)
          call iswap(PL,k,jj+k-1)
          call iswap(PQ,k,jj+k-1)
!          call CauchyMaxsvd( F(k),D(k),V(k),N-k+1,junkU,jjU) ! row
          call CauchyMaxsvdm( D(k),FF,PU(k),PL(k),DIFL,DIFR,V( k ),N-k+1,junkU,jjU )  ! row
          nflops = nflops + 4* (N-k+1)
          IF(jjU == 1) RETURN
          
          flg = 0
          CONTINUE
       END IF
       jj = jjU
       call dswap(D,k,jj+k-1)
       call dswap(V,k,jj+k-1)
       call iswap(PU,k,jj+k-1)
!       call CauchyMaxsvd( D(k),F(k),U(k),M-k+1,junkL,jjL ) ! column
       call CauchyMaxsvdmc( D(k),FF,PU(k),PL(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! column
       nflops = nflops + 4* (M-k+1)
       IF(jjL == 1) RETURN

       flg = 1
    END DO

  END SUBROUTINE CauchyPivtGsvd

!!!!!!
  SUBROUTINE CauchyPivtGsvdm(D,F,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR)
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
! U and V, where S is the k-th Schur complement and S(i,j)=U(i)*V(j) / (D(i)**2-F(j)**2 ). 
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
!   The row generator of S, and D contains the updated singular values. S is an 
!   (M-k)-by-(N-k) matrix. 
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
!    DIFR(i) = D0(i) - FF(i+1), negative
! 
! ==========
! Modified by S.-G. Li, in Dec. 24th, 2012
! 
! Modified on April 9th, 2013
! 1) Remove the use of DD, and use D directly.
! ==========================================
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
    call CauchyMaxsvdm ( D(k),FF,PL(k),PU(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! largest entry of the first column
    call CauchyMaxsvdmc( D(k),FF,PL(k),PU(k),DIFL,DIFR,V( k ),N-k+1,junkU,jjU )  ! largest entry of the first row
    junkL = junkL * ABS( V(k) )
    junkU = junkU * ABS( U(k) ) 

    ! compute the first entry of k-th Schur complement
    ! Because of off-diagonal blocks, PL(k) can not equal to PU(k)
    IF( PL(k) .lt. PU(k) ) THEN
       Piv = ABS( U(k)*V(k) / ( (FF(PL(k)+1)-F(k))+DIFR(PL(K)) ) / ( D(k)+F(k) ) )
    ELSE
       Piv = ABS( U(k)*V(k) / ( (FF(PL(k))-F(k))+DIFL(PL(k)) ) / ( D(k)+F(k) ) ) 
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
          call CauchyMaxsvdmc( D(k),FF,PL(k),PU(k),DIFL,DIFR,V( k ),N-k+1,junkU,jjU )  ! row
!          call CauchyMaxsvdm( F(k),D(k),V(k),N-k+1,junkU,jjU) ! row
          nflops = nflops + 4* (N-k+1)
          IF(jjU == 1) RETURN
          
          flg = 0
          CONTINUE
       END IF
       jj = jjU
       call dswap(F,k,jj+k-1)
       call dswap(V,k,jj+k-1)
       call iswap(PU,k,jj+k-1)
       call CauchyMaxsvdm ( D(k),FF,PL(k),PU(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! column
!       call CauchyMaxsvdm( D(k),F(k),U(k),M-k+1,junkL,jjL ) ! column
       nflops = nflops + 4* (M-k+1)
       IF(jjL == 1) RETURN

       flg = 1
    END DO

  END SUBROUTINE CauchyPivtGsvdm

!!!!!!
  SUBROUTINE CauchyPivtGsvdm_VP( D,F,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR )
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
! where S is the k-th Schur complement and S(i,j) = U(i)*V(j) / (D(i)**2 - F(j)**2 ). 
! S is an M-by-N matrix. 
! 
! Initially U == A, and all the parameters are only permuted. This pivoting strategy is 
! a variant of complete pivoting. 
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
! ==============
! D (inout) DOUBLE PRECISION array, DIMENSION( N ) 
!   The row generator of S.
! 
! F (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!   The column generator of S.
! 
! U (inout) DOUBLE PRECISION, DIMENSION( M )
!   The row generators of S
! 
! V (inout) DOUBLE PRECISION, DIMENSION( N )
!   The column generators of S
! 
! DD (in) DOUBLE PRECISION, DIMENSION( NN )
!    The updated singular values of original matrix, the difference 
!    between DD(i) and DD(j) is computed via FF(k), DIFL and DIFR. 
!    DD(i) and DD(j) may be equal 'mistakely'. 
!
! FF (in) DOUBLE PRECISION, DIMENSION( NN )
!    The old singular values of original matrix
! 
! DIFL (in) DOUBLE PRECISION, DIMENSION( NN )
!    DIFL(i) = D(i) - F(i), positive
! 
! DIFR (in) DOUBLE PRECISION, DIMENSION( NN )
!    DIFR(i) = D(i) - F(i+1), negative
! 
! PL (inout) INTEGER, DIMENSION( M )
!    The row permutations corresponding the original generators
!
! PU (inout) INTEGER, DIMENSION( N )
!    The column permutations corresponding the original generators
!
! k  (in) INTEGER
!    The k-th Schur complement
! 
! A  (inout) DOUBLE PRECISION, DIMENSION( M )
!    The original generator of U
! 
! M  (in) INTEGER, the row dimension of S
! 
! N  (in) INTEGER, the column dimension of S
!    
! PQ (inout) INTEGER array, DIMENSION( M )
!    Records the local row permutation, and PQ is used to computed the inverse of row permutations
!    for current HSS blk
! 
! nflops (out) INTEGER, the floating point operations during this pivoting
! 
! ==========
! Modified by S.-G. Li, in Dec. 24th, 2012
! ========================================
!
! .. Local Parameters ..
    DOUBLE PRECISION  ZERO, ONE, RelCont
    PARAMETER         ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! .. Local Scalars ..
    INTEGER           j,jj,npivot,MN,ierr,ii,flg
    DOUBLE PRECISION  junk,Piv, Ratio
! ..
! .. Local Arrays ..
    integer :: iiU(10)
    double precision :: junkU(10)
! ..
! .. Intrinsic Functions ..
    INTRINSIC    ABS,MIN

    nflops = 0
    flg = 0
    npivot = 0
    MN = MIN(10, N-k+1) 
    RelCont = ONE - 1.0E-13
    
    DO WHILE( flg .ne. 1 )
       call CauchyPivtGsvdm(D,F,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR)

       ! maximun of first 10 columns
       DO j =1, MN
          call CauchyMaxsvdm( D(k),FF,PL(k),PU(k+j-1),DIFL,DIFR,U( k ),M-k+1,junkU(j),iiU(j) )  ! largest entry of the first column
          junkU(j) = abs( junkU(j)*V(k+j-1) ) 
       END DO

       junk = junkU(1)
       ii = iiU(1)
       jj = 1
       DO j = 2, MN
          IF(junkU(j) > junk ) THEN
             junk = junkU(j)
             ii = iiU(j)
             jj = j
          END IF
       END DO

       ! compute the first entry of k-th Schur complement
       ! Because of off-diagonal blocks, PL(k) can not equal to PU(k)
       IF( PL(k) .lt. PU(k) ) THEN
          Piv = ABS( U(k)*V(k) / ( (FF(PL(k)+1)-F(k))+DIFR(PL(K)) ) / ( D(k)+F(k) ) )
       ELSE
          Piv = ABS( U(k)*V(k) / ( (FF(PL(k))-F(k))+DIFL(PL(k)) ) / ( D(k)+F(k) ) ) 
       END IF
       nflops = nflops + 4* (M+N-2*k+3)

       Ratio = Piv/junk
       IF( Ratio .GE. RelCont ) THEN
          flg = 1
!          write(*,*) 'Pivot', npivot
          RETURN
       ELSE
          npivot = npivot +1
          call dswap(D,k,ii+k-1)
          call dswap(U,k,ii+k-1)
          call dswap(A,k,ii+k-1)
          call iswap(PL,k,ii+k-1)
          call iswap(PQ,k,ii+k-1)

          call dswap(F,k,jj+k-1)
          call dswap(V,k,jj+k-1)
          call iswap(PU,k,jj+k-1)
       END IF
    
    END DO

  END SUBROUTINE CauchyPivtGsvdm_VP

!!!!!!
  SUBROUTINE CauchyPivtGsvd_VP( F,D,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR )
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
! where S is the k-th Schur complement and S(i,j) = U(i)*V(j) / (F(i)**2 - D(j)**2 ).
! S is an M-by-N matrix. This routine is designed for rrluCauchysvdcol. 
! Here D contains the updated singular values, and F are the old singular values.
! As in rrluCauchysvdcol, F is the row generator of S, and D is the column generator
! of S. 
! 
! Initially U == A, and all the parameters are only permuted. This pivoting strategy 
! is a modification of Gu's strategy. The idea is
!   First use Gu's strategy to pivot a large entry, and then check whether it is the
!   largest of the first TEN columns; 
!   If it does, 
!       return; 
!   else  
!       permute the largest to the front and then use Gu's strategy again.
!   end 
!
!  Further details
!  This pivoting strategy is not good enough, and it is not quite useful. Junk. 
! 
! ..Parameters..
! ==============
! F (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!   The row generator of S, which are the old singular values, and can be used
!   directly. 
! 
! D (inout) DOUBLE PRECISION array, DIMENSION( N ) 
!   The column generator of S, the updated singular values.
! 
! U (inout) DOUBLE PRECISION, DIMENSION( M )
!   The row generators of S
! 
! V (inout) DOUBLE PRECISION, DIMENSION( N )
!   The column generators of S
!
! PL (inout) INTEGER, DIMENSION( M )
!    The row permutations corresponding the original generators
!
! PU (inout) INTEGER, DIMENSION( N )
!    The column permutations corresponding the original generators
!
! k  (in) INTEGER
!    The k-th Schur complement
! 
! A  (inout) DOUBLE PRECISION, DIMENSION( M )
!    The original generator of U
! 
! M  (in) INTEGER, the row dimension of S
! 
! N  (in) INTEGER, the column dimension of S
!    
! PQ (inout) INTEGER array, DIMENSION( M )
!    Records the local row permutation, and PQ is used to computed the inverse of row permutations
!    for current HSS blk
! 
! nflops (out) INTEGER, the floating point operations during this pivoting
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( * )
!    The original old singular values
!
! DD (in) DOUBLE PRECISION array, DIMENSION( * )
!    The original updated singular values. 
!
! DIFL (in) DOUBLE PRECISION, DIMENSION( NN )
!    DIFL(i) = FF(i) - DD(i), positive
! 
! DIFR (in) DOUBLE PRECISION, DIMENSION( NN )
!    DIFR(i) = FF(i) - DD(i+1), negative
!
! ==========
! Modified by S.-G. Li, in April 1st, 2013
! Modified on April 11th, 2013
! 1) Exchange D and F, remove the use of DD and use D directly.  
!    Now D are the updated svals
! =================================================
!
! .. Local Parameters ..
    DOUBLE PRECISION  ZERO, ONE, RelCont
    PARAMETER         ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! .. Local Scalars ..
    INTEGER           j,jj,npivot,MN,ii,flg
    DOUBLE PRECISION  junk,Piv, Ratio
! ..
! .. Local Arrays ..
    integer :: iiU(10)
    double precision  :: junkU(10)
! ..
! .. Intrinsic Functions ..
    INTRINSIC    ABS,MIN

    nflops = 0
    flg = 0
    npivot = 0
    MN = MIN(10, N-K+1)
    RelCont = ONE - 1.0E-13

    DO WHILE( flg .ne. 1 )
       call CauchyPivtGsvd(F,D,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR) 
       
       ! maximun of first 10 columns
       DO j =1, MN
       call CauchyMaxsvdmc( D(k),FF,PU(k+j-1),PL(k),DIFL,DIFR,U( k ),M-k+1,junkU(j),iiU(j) )   ! largest entry of the (k+j)-th column
          junkU(j) = abs( junkU(j)*V(k+j-1) ) 
       END DO

       junk = junkU(1)
       ii = iiU(1)
       jj = 1
       DO j = 2, MN
          IF(junkU(j) > junk ) THEN
             junk = junkU(j)
             ii = iiU(j)
             jj = j
          END IF
       END DO

       ! Compute the first entry of the K-th Schur Complement
       ! Because of off-diagonal blocks, PL(k) can not equal to PU(k)
       IF( PL(k) .LT. PU(K) ) THEN
          Piv = ABS( U(k)*V(k) / ( ( F(k)-FF(PU(k)) )-DIFL(PU(k)) ) / ( F(k)+D(k) ) )
       ELSE
          Piv = ABS( U(k)*V(k) / ( ( FF(PU(k)+1)-F(k) )+DIFR(PU(k)) ) / ( F(k)+D(k) ) )
       END IF
       nflops = nflops + 4* (M+N-2*k+3)

       Ratio = Piv/junk
       IF( Ratio .GE. RelCont  ) THEN
          flg = 1
!          write(*,*) 'Pivot', npivot
          RETURN
       ELSE
          flg = 0
          npivot = npivot + 1
          call dswap(F,k,ii+k-1)
          call dswap(U,k,ii+k-1)
          call dswap(A,k,ii+k-1)
          call iswap(PL,k,ii+k-1)
          call iswap(PQ,k,ii+k-1)

          call dswap(D,k,jj+k-1)
          call dswap(V,k,jj+k-1)
          call iswap(PU,k,jj+k-1)
       END IF

    END DO ! (WHILE)

  END SUBROUTINE CauchyPivtGsvd_VP

!!!!!!
  SUBROUTINE CauchyPivtGsvdm_CP(D,F,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR)
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
! where S is the k-th Schur complement and S(i,j) = U(i)*V(j) / (D(i)**2 - F(j)**2 ). 
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
!   The row generator of S, the updated singular values, and its relation with PL is
!   D(1:M) = D0( PL(1:M) );
! 
! F (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!   The column generator of S, which is useless, and the old singular values;
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
!    DIFL(i) = D0(i) - FF(i), positive
! 
! DIFR (in) DOUBLE PRECISION, DIMENSION( LN )
!    DIFR(i) = D0(i) - FF(i+1), negative
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
! Modified by S.-G. Li, in Dec. 24th, 2012
!
! Modified on April, 9th, 2013
! 1) Remove the use of DD, and use D directly.
! ===========================================
!
! .. Local Parameters ..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! .. Local Scalars ..
    INTEGER           jjL,jjU,ii,jj, j
    DOUBLE PRECISION  junkL,junkU,junk
! ..
! .. Intrinsic Functions ..
    INTRINSIC    ABS,MIN

    nflops = 0
    call CauchyMaxsvdm( D(k),FF,PL(k),PU(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! largest entry of the first column
    junkL = ABS( junkL*V(k) )
    junk = junkL
    ii = jjL
    jj = 1

    DO j = 1, N-k
       call CauchyMaxsvdm( D(k),FF,PL(k),PU(k+j),DIFL,DIFR,U( k ),M-k+1,junkU,jjU )  ! largest entry of the first column
       junkU = ABS( junkU*V(k+j) )
       IF( junkU > junk ) THEN
          junk = junkU
          ii = jjU
          jj = j+1
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
    
  END SUBROUTINE CauchyPivtGsvdm_CP

!!!!!!
  SUBROUTINE CauchyPivtGsvd_CP( F,D,U,V,PL,PU,k,A,M,N,PQ,nflops,FF,DIFL,DIFR )
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
! where S is the k-th Schur complement and S(i,j) = U(i)*V(j) / (F(i)**2 - D(j)**2 ).
! S is an M-by-N matrix. This routine is designed for rrluCauchysvdcol. 
! Here F are the old singular values and D are the updated singular values. F is the
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
!   The row generator of S, which are the old singular values, and can be used
!   directly. 
! 
! D (inout) DOUBLE PRECISION array, DIMENSION( N ) 
!   The column generator of S, the updated singular values.
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
! M  (in) INTEGER, the row dimension of S
! 
! N  (in) INTEGER, the column dimension of S
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
! DD (in) DOUBLE PRECISION array, DIMENSION( * )
!    The original updated singular values. 
!
! DIFL (in) DOUBLE PRECISION, DIMENSION( NN )
!    DIFL(i) = FF(i) - DD(i), positive
! 
! DIFR (in) DOUBLE PRECISION, DIMENSION( NN )
!    DIFR(i) = FF(i) - DD(i+1), negative
!
! ==========
! Modified by S.-G. Li, in Dec. 24th, 2012
! Modified on April 10th, 2013
! 1) Exchanged D and F, and remove the use DD, use D directly.
! ===========================================================
!
! .. Local Parameters ..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! .. Local Scalars ..
    INTEGER           jjL,jjU,ii,jj, j
    DOUBLE PRECISION  junkL,junkU,junk
! ..
! .. Intrinsic Functions ..
    INTRINSIC    ABS,MIN

    nflops = 0
    call CauchyMaxsvdmc( D(k),FF,PU(k),PL(k),DIFL,DIFR,U( k ),M-k+1,junkL,jjL )  ! largest entry of the first column
    junkL = ABS( junkL*V(k) )
    junk = junkL
    ii = jjL
    jj = 1
    
    DO j = 1, N-k
       call CauchyMaxsvdmc( D(k),FF,PU(k+j),PL(k),DIFL,DIFR,U( k ),M-k+1,junkU,jjU )  ! largest entry of the first column
       junkU = ABS( junkU*V(k+j) )
       IF( junkU > junk ) THEN
          junk = junkU
          ii = jjU
          jj = j+1
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

  END SUBROUTINE CauchyPivtGsvd_CP

!!!!!!
  SUBROUTINE searchMax2svd(U,W,D2,D1,LDU,LDW,ii,jj,Nmax,nflops)
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDU, LDW
    INTEGER, INTENT(OUT) :: ii, jj, nflops
    DOUBLE PRECISION, INTENT(OUT) :: Nmax
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN) :: U(*),W(*),D2(*),D1(*)
!
! Purpose 
! ========
! Choose the maximum entry of Z by searching column by column, and Z is defined as
! D2* Z -Z*D1 = U*W, a Cauchy-like matrix. This routine is designed for rrluCauchysvdcol,
! and here D1 and D2 are the old singular values, which can be used directly. 
! D2 is a LDU-by-1 vector and D1 is a LDW-by-1 vector. 
! 
! ..Parameters..
! ================
! U  (in) DOUBLE PRECISION array, DIMENSION( LDU )
!    The row generators of Z
! 
! W  (in) DOUBLE PRECISION array, DIMENSION( LDW )
!    The column generators of Z
! 
! LDU (in) INTEGER, row dimension of Z
! 
! LDW (in) INTEGER, column dimension of Z
!
! ==========
! Modified by S.-G. Li, on Dec. 24th, 2012
! ========================================
!
! .. Local Parameters ..
    DOUBLE PRECISION ZERO
    PARAMETER        ( ZERO = 0.0E+0 )
! .. Local Scalars ..
    INTEGER          :: j, jjL
    DOUBLE PRECISION :: junk

    nflops = 0    
    call CauchyMaxsvd( D2,D1(1),U,LDU,junk,jjL )
    junk = junk*ABS( W(1) )       
    Nmax = junk
    ii = 1
    jj = jjL
    nflops = nflops + 4* LDU*LDU
    
    DO j = 2, LDW   ! ii: col, jj: row
       call CauchyMaxsvd( D2,D1(j),U,LDU,junk,jjL )
       junk = junk*ABS( W(j) )       
       IF(junk .GT. Nmax) THEN
          Nmax = junk
          ii = j
          jj = jjL
       END IF
    END DO

  END SUBROUTINE searchMax2svd

!!!!!!
  SUBROUTINE searchMax2svdm(U,W,D,FF,PL,DIFL,DIFR,LDU,LDW,ii,jj,Nmax,nflops)
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
! D1 and D2 are updated singular values, and therefore their difference should be
! computed via old singular values and DIFL and DIFR. 
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
! DD (in) DOUBLE PRECISION array, DIMENSION( N )
!    The original updated singular values, and its entries are referred by PL.
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( N )
!    The original old singular values, and its entries are referred by PL.
!    It is used to compute  DD(i)-DD(j).
!
! PL (in) INTEGER array, DIMENSION( LDW+LDU )
!    The permutations of DD for current node
!
! DIFL (in) DOUBLE PRECISION array, DIMENSION( N )
!      The original differences between updated singular values and old ones,
!      DIFL(i) = DD(i) - FF(i), positive
!
! DIFR (in) DOUBLE PRECISION array, DIMENSION( N )
!      The original differences between updated singular values and old ones,
!      DIFR(i) = DD(i) - FF(i+1), negative
!
! LDU (in) INTEGER, row dimension of Z
! 
! LDW (in) INTEGER, column dimension of Z
!
! ii  (out) INTEGER 
!     The row index of chosen largest entry in Z
! 
! jj  (out) INTEGER 
!     The column index of chosen largest entry in Z
! 
! ==========
! Modified by S.-G. Li, on Mar. 13th, 2013
! Modified on April 12th, 2013
! 1) Remove the use of DD and change it to D
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
    call CauchyMaxsvdmm(D,FF,PL,1,DIFL,DIFR,U,LDU,LDW,junk,jjL)
    junk = junk*ABS( W(1) )       
    Nmax = junk
    ii = 1
    jj = jjL
    nflops = nflops + 4* LDU*LDU
    
    DO j = 2, LDW   ! ii: col, jj: row
!       call CauchyMaxsvd(D2,D1(j),U,LDU,junk,jjL)
       call CauchyMaxsvdmm(D,FF,PL,j,DIFL,DIFR,U,LDU,LDW,junk,jjL)
       junk = junk*ABS( W(j) )       
       IF(junk .GT. Nmax) THEN
          Nmax = junk
          ii = j
          jj = jjL
       END IF
    END DO

  END SUBROUTINE searchMax2svdm

!!!!!!
  SUBROUTINE CauchyMaxsvd( D,F,U,N,junk,jj )
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
! This routine returns the largest entry in magnitude and its index of one column of matrix A,
! where A is a Cauchy-like matrix, defined as A(i,j) = U(i)*V(j) / (D(i)**2 - F(j)**2 ). 
! Here F = F(j), and D and F are only variables, have no relation with new or old singular 
! values. 
!
! ===========
! Written by S.-G. Li, on Dec. 23th, 2012, in Changsha China
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

    Lk = U(1:N)/ ( (D(1:N)+F)*(D(1:N)-F) )
    Lk = ABS( Lk )
    junk = MAXVAL( Lk(1:N) )
    temp = MAXLOC( Lk(1:N) )
    jj = temp(1)
    
    DEALLOCATE( LK )
    
  END SUBROUTINE CauchyMaxsvd

!!!!!!
  SUBROUTINE CauchyMaxsvdmm( D,FF,PL,k,DIFL,DIFR,U,N,M,junk,jj )
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
! where A is a Cauchy-like matrix, defined as A(i,j) = U(i)*V(j) / (D2(i)**2 - D1(j)**2 ),
! an N-by-M matrix. D1(1:M) = D0(PL(1:M))=D(1:M), D2(1:N) = D0(PL(1+M:N+M))=D(1+M: M+N), 
! and their entries are updated singular values, and FF are old singular values which are 
! used to compute the difference between D1(i) and D2(j). 
!
! ..Parameters..
! ==============
! D (in) DOUBLE PRECISION array, DIMENSION( M+N )
!    The original updated singular values, and its entries are corresponding to PL; 
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
!      The original differences between updated singular values and old ones,
!      DIFL(i) = DD(i) - FF(i), positive
!
! DIFR (in) DOUBLE PRECISION array, DIMENSION( * )
!      The original differences between updated singular values and old ones,
!      DIFR(i) = DD(i) - FF(i+1), negative
!
! U   (in) DOUBLE PRECISION array, DIMENSION( N ) 
!     The row generators of A
!
! N   (in) INTEGER
!     The row dimension of A
! 
! M   (in) INTEGER
!     The column dimension of A
! 
! Junk (out) DOUBLE
!      The absolute value of largest entry in the k-th column of A
! 
! jj   (out) INTEGER 
!      The row index of Junk
! ===========
! Written by S.-G. Li, on Mar. 13th, 2013, in Changsha China
! Modified on April, 9th, 2013
! 1) Remove the use of DD, and use D directly.
! ==========================================================
!
! ..Local Scalars..
    INTEGER temp(1), I, K2, K1
! ..
! ..Intrinsic Functions..
    INTRINSIC    MAXLOC, ABS,MAXVAL
! ..
! ..Local Arrays..
    DOUBLE PRECISION, ALLOCATABLE :: LK(:)
    ALLOCATE( LK(N) )

    K1 = PL(k) 
    DO i = 1, N
       K2 = PL(i+M)
       IF( K2 .GT. K1 ) THEN
          IF( K2 .EQ. (K1+1) ) THEN
             LK(i) = U(i) / ( DIFL(K2)-DIFR(K1) ) / ( D(k)+D(i+M) )
          ELSE
             LK(i) = U(i) / ( (FF(K2)-FF(K1+1))+DIFL(K2)-DIFR(K1) ) / ( D(k)+D(i+M) )
          END IF
       ELSE
          IF( K1 .EQ. (K2+1) ) THEN
             LK(i) = U(i) / ( DIFL(K1)-DIFR(K2) ) / ( D(k)+D(i+M) )
          ELSE
             LK(i) = U(i) / ( (FF(K1)-FF(K2+1))+DIFL(K1)-DIFR(K2) ) / ( D(k)+D(i+M) )
          END IF
       END IF
    END DO
    Lk = ABS( Lk )
    junk = MAXVAL( Lk(1:N) )
    temp = MAXLOC( Lk(1:N) )
    jj = temp(1)
    
    DEALLOCATE( LK )
    
  END SUBROUTINE CauchyMaxsvdmm

!!!!!!
  SUBROUTINE CauchyMaxsvdm( D,FF,PL,Fk,DIFL,DIFR,U,N,junk,jj )
! 
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(IN)  :: Fk
    INTEGER, INTENT(OUT) :: jj
    DOUBLE PRECISION, INTENT(OUT) :: junk
! 
! .. Array Arguments ..
    INTEGER, INTENT(IN) :: PL(*)
    DOUBLE PRECISION, INTENT(IN) :: D(*),FF(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(IN) :: U(*)
! 
! Purpose
! =======
!  This routine returns the largest entry in magnitude and its index of one column of matrix A,
!  where A is a Cauchy-like matrix, defined as A(i,j) = U(i)*V(j) / (D(i)**2 - F(j)**2 ). 
!  Here F = F(j), and F is a scalar, one of the old singular values. D are 
!  the updated singular values, and the relation between D and D0 is D(1:N) = D0( PL(1:N) ). 
!  The leading dimension of A is N. 
!
! ..Parameters..
! ==============
! D  (in) DOUBLE PRECISION, DIMENSION( N )
!    The updated singular values, its entries have relationship with PL by 
!    D(1:N) = D0( PL(1:N) ). 
!
! FF (in) DOUBLE PRECISION array, DIMENSION( LN )
!    LN is the size of original problem. FF contains the original singular values.
! 
! PL (in) INTEGER, DIMENSION(N)
!    The permutation of updated singular values, the row generators of A
! 
! Fk  (in) INTEGER 
!     The original position of F(j) in FF, i.e., FF(Fk) = F;
! 
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!     DIFL(i) = D(i)-F(i), the order of which will not change and its elements
!     are referenced by using PL or PU. Positive
!
! DIFR (in) DOULBE PRECISION array, DIMENSION ( LN )
!     DIFR(i) = D(i)-F(i+1), the order of which will not change and its elements
!     are referenced by using PL or PU. Negative
!
! U  (in) DOUBLE PRECISION array, DIMENSION( N )
!     The row generator of A.
!      
! N  (in) INTEGER, leading dimension of A.
! 
! junk (out) DOUBLE PRECISION, the computed largest entry of magnitude.
!
! jj (out) INTEGER, the row dimension of computed largest entry
! 
! ===========
!  Written by S.-G. Li, on Mar. 13th, 2013, in Changsha China
!  Modified on April 9th, 2013
!  1) Remove the use of DD, and use D directly. 
! ==========================================================
!
! .. Local Scalars ..
    INTEGER temp(1), I, PLi
    DOUBLE PRECISION F
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MAXLOC, ABS,MAXVAL
! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: LK(:)
    ALLOCATE( LK(N) )

    F = FF(Fk)
    DO i = 1, N
       PLi = PL(i)
       IF( PLi .EQ. Fk ) THEN
          LK(i) = U(i) / DIFL( PLi ) / ( D(i)+F )
       ELSE
          IF( PLi .GT. Fk ) THEN
            LK(i) = U(i) / ( (FF(PLi)-F)+DIFL(PLi) ) / ( D(i)+F )
          ELSE
            LK(i) = U(i) / ( (FF(PLi+1)-F)+DIFR(PLi) ) / ( D(i)+F ) 
         END IF
      END IF
    END DO
    Lk = ABS( Lk )
    junk = MAXVAL( Lk(1:N) )
    temp = MAXLOC( Lk(1:N) )
    jj = temp(1)
    
    DEALLOCATE( LK )
    
  END SUBROUTINE CauchyMaxsvdm

!!!!!!
  SUBROUTINE CauchyMaxsvdmc( D,FF,PL,PU,DIFL,DIFR,V,N,junk,jj )
! 
! ..Scalar Arguments..
    INTEGER, INTENT(IN)  :: N, PL
    INTEGER, INTENT(OUT) :: jj
    DOUBLE PRECISION, INTENT(IN)  :: D
    DOUBLE PRECISION, INTENT(OUT) :: junk
! 
! ..Array Arguments..
    INTEGER, INTENT(IN) :: PU(*)
    DOUBLE PRECISION, INTENT(IN) :: FF(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(IN) :: V(*)
! 
! Purpose
! =======
! This routine returns the largest entry of magnitude and its index of one row of matrix A,
! where A is a Cauchy-like matrix, defined as A(i,j) = U(i)*V(j) / (D(i)**2 - F(j)**2 ). 
! Here D = D(i), F is referred by FF, F(1:N) = FF( PU(1:N) ). D is a scalar, one of the updated singular
! values, and FF are the updated singular values. 
!
! ..Parameters..
! ==============
! D (in) DOUBLE PRECISION scalar
!    One of the original updated singular values, the row of generator of one row of A,
!    and D = D0( PL ).
!
! FF  (in) DOUBLE PRECISION, DIMENSION( N )
!    The original singular values, its entries are referred by PU;
! 
! PL (in) INTEGER array
!    The permutation of updated singular values, such that D = D0(PL);
! 
! PU (in) INTEGER, DIMENSION( N )
!    The permutation of old singular values, the column generators of A;
! 
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!     DIFL(i) = D0(i)-FF(i), the order of which will not change and its elements
!     are referenced by using PL or PU. Positive
!
! DIFR (in) DOULBE PRECISION array, DIMENSION ( LN )
!     DIFR(i) = D0(i)-FF(i+1), the order of which will not change and its elements
!     are referenced by using PL or PU. Negative
! 
! V  (in) DOUBLE PRECISION array, DIMENSION( N )
!    The column generator of A;
! 
! junk (out) DOUBLE PRECISION, the computed largest entry of magnitude;
! 
! jj  (out) INTEGER, the column index of computed largest entry;
! 
! ===========
! Written by S.-G. Li, on Mar. 13th, 2013, in Changsha China
! 
! Modified on April 9th, 2013
! 1) Remove the use of DD, and use D directly. 
! ==========================================================
!
! ..Local Scalars..
    INTEGER temp(1), I, PUi, Dk
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MAXLOC, ABS,MAXVAL
! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: LK(:)
    ALLOCATE( LK(N) )

    Dk = PL
    DO i = 1, N
       PUi = PU(i)
       IF( PUi .EQ. Dk ) THEN
          LK(i) = V(i) / DIFL( PUi ) / ( D+FF(PUi) )
       ELSE
          IF( PUi .LT. Dk ) THEN
            LK(i) = V(i) / ( (FF(Dk)-FF(PUi))+DIFL(Dk) ) / ( D+FF(PUi) )
          ELSE
            LK(i) = V(i) / ( (FF(Dk+1)-FF(PUi))+DIFR(Dk) ) / ( D+FF(PUi) ) 
         END IF
      END IF
    END DO
    Lk = ABS( Lk )
    junk = MAXVAL( Lk(1:N) )
    temp = MAXLOC( Lk(1:N) )
    jj = temp(1)
    
    DEALLOCATE( LK )
    
  END SUBROUTINE CauchyMaxsvdmc

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

!!!!!!
  SUBROUTINE Comprcauchysvd(rowcol,D,F,U,V,DIFL,DIFR,FF,PL,PU,tol,M,N,PH,pnh,H,nodi,Rk,time,nflops,nswap)
    USE aux_hss
!
! Scalar parameters
    INTEGER, INTENT(IN)    :: M,N,nodi
    INTEGER, INTENT(OUT)   :: Rk,nflops,nswap
    INTEGER, INTENT(INOUT) :: pnh
    DOUBLE PRECISION, INTENT(IN)  :: tol
    DOUBLE PRECISION, INTENT(OUT) :: time
    CHARACTER(LEN=1), INTENT(IN)  :: rowcol
! Array parameters
    DOUBLE PRECISION, INTENT(IN)    :: DIFL(*),DIFR(*),FF(*)
    DOUBLE PRECISION, INTENT(INOUT) :: D(*),F(*),U(*),V(*),H(*)
    INTEGER, INTENT(INOUT)     :: PL(M),PU(N)
    TYPE(HSSMM), INTENT(INOUT) :: PH
!
! Purpose
! ========
! This routine computes a low rank approximation of some off-diagonal block of a
! Cauchy-like matrix C(i,j) = U(i)*V(j) / ( D(i)**2 - F(j)**2 ). 
! This one is written for computing the transpose of right singular vector matrix, VT. 
! The difference between this one and ComprcauchysvdU is that one is another's transpose. 
!
! Parameters
! ==========
! rowcol (in) CHARACTER ( LEN = 1)
!        = 'R', row compression
!        = 'C', column compresion
! 
! D  (in) DOUBLE PRECISION array, DIMENSION( M )
!    D are permuted in essential and only the first Rk entries are useful. 
!    But its permutation is recorded by PL and entries of D will not change. The entries of D 
!    are the new singular values and are the row generators of Cauchy-like matrix A. D and PL are in 
!    one group.
!
! F  (in)  DOUBLE PRECISION array, DIMENSION( N )
!    F will be permuted through PU in the process. Finally the permuted of F
!    will not be used since we are using interpolative decomposition and F corresponds
!    to the column. The entries of F are the old singular values. 
!
! U  (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!    U is the column normilzation scalars, ALPHA in our old notation and it will be perumuted, only 
!    the first Rk entries are useful for row compression. 
! 
! V  (inout)  DOUBLE PRECISION array, DIMENSION( N )
!    V is the Z in our old notation. It will be permuted and 
!    the first Rk entries are useful for column compression and V will not be used in 
!    row compression.
!
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!      DIFL(i) = D(i)-F(i), the order of which will not change and its elements
!      are also referenced by using PL or PU. Positive
! 
! DIFR (in) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFR(i) = D(i)-F(i+1), the order of which will not change and its elements
!      are also referenced by using PL or PU. Negative
!
! FF (in) DOUBLE PRECISION array, DIMENSION( LN ), original F
!    LN is the size of original problem. 
!
! PL (inout) INTEGER array, DIMENSION( M )
!    It stores the permutations for row generators of A. It corresponds to F for row compression and
!    it corresponds to D for column compression since 
!
! PU (inout) INTEGER array, DIMENSION( N )
!    It stores the permutations for col generators of A. It corresponds to D for row compression and
!    it corresponds to F for column compression since 
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
! time  (out) DOUBLE PRECISION, the accumulated time of calling rrluCauchy.
!
! nflops (out) INTEGER, the total flops of constructing this HSS matrix 
!
! nswap  (out) INTEGER, the total swaps whiling constructing this HSS matrix
!
! =========
! Written by S.G. Li, on Dec. 15th, 2012
! ======================================
!
! .. Parameters ..
    DOUBLE PRECISION ZERO, ONE
    PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0)
! ..    
! .. Local scalars ..
    INTEGER          mn,i,info,lflops, ierr
    DOUBLE PRECISION time1
    LOGICAL          CR
! ..    
! .. Local arrays ..
    INTEGER, ALLOCATABLE :: PLQ(:)
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
    IF( CR )  THEN
       ALLOCATE( Z(M),W(MN),PLQ(M), stat = ierr )    ! Check ??
       Z = ZERO
       W = ZERO
       call cpu_time(time)
       call rrluCauchysvdrow( D,F,U,V,DIFL,DIFR,tol,Z,W,FF,PL,PU,M,N,Rk,PLQ,lflops,nswap )
       nflops = nflops + lflops
       call cpu_time(time1)
       time = time1 - time
       call invp(PLQ,M)       ! PL --> InvPL
       IF( Rk .LT. MN ) THEN            
          ALLOCATE( Q(M,Rk) )            
          call Cauchylike2svd(Q,D(Rk+1),D,Z(Rk+1),W,M,Rk,FF,DIFL,DIFR,PL,lflops)
!!$            write(*,*) 'Row Rk is ', Rk, 'm is ', m, 'n = ',n, JunkN
!!$            call testaccury('r',D,F,U,V1,Q,M,N,Rk)
          nflops = nflops + lflops
          Q(1:M,1:Rk) = Q(PLQ,1:Rk)            
          call dlacpy('A',M,Rk,Q,M,H(pnh),M)        ! copy U to H
          call hssexpmm2('U',PH,pnh,M,Rk,nodi,info) ! copy Q to generators            
       ELSE
          ! copy identity matrix to generators
          allocate( Q(M,M) )
          Rk = M
          Q(1:M,1:M) = ZERO
          DO i = 1,M
             Q(i,i) = ONE
          END DO
          Q(1:M,1:M) = Q(PLQ,1:M)
          call dlacpy('A',M,Rk,Q,M,H(pnh),M)        ! copy U to H
          call hssexpmm2('U',PH,pnh,M,Rk,nodi,info) ! copy Q to generators
       END IF

    ELSE ! block col

       ALLOCATE( Z(N),W(MN),PLQ(N) )
       Z = ZERO
       W = ZERO
       call cpu_time(time)
       call rrluCauchysvdcol(F,D,V,U,DIFL,DIFR,tol,Z,W,FF,PU,PL,N,M,Rk,PLQ,lflops,nswap )
       nflops = nflops + lflops
       call cpu_time(time1)
       time = time1 - time
       call invp(PLQ,N)  ! PL --> InvPL
       IF( Rk .LT. MN ) THEN
          ALLOCATE( Q(Rk,N) )
          call Cauchylike2svdcol( Q,F(Rk+1),F,Z(Rk+1),W,Rk,N,lflops )
!!$            write(*,*) 'Col Rk is ', Rk, 'm is ', m, 'n = ',n, JunkN
!!$            call testaccury('c',D,F,V1,V,Q,Rk,M,N)            
          nflops = nflops + lflops
          Q(1:Rk,1:N) = Q(1:Rk,PLQ)
          call dlacpy('A',Rk,N,Q,Rk,H(pnh),Rk)      ! copy V to H
          call hssexpmm2('V',PH,pnh,Rk,N,nodi,info) ! copy Q to generators            
       ELSE
          ! copy identity matrix to generators
          ALLOCATE( Q(N,N) )
          Rk = N
          Q(1:N,1:N) = ZERO
          DO i = 1,N
             Q(i,i) = ONE
          END DO
          Q(1:N,1:N) = Q(1:N,PLQ)            
          call dlacpy('A',N,N,Q,N,H(pnh),N)         ! copy V to H
          call hssexpmm2('V',PH,pnh,Rk,N,nodi,info) ! copy Q to generators
       END IF

    END IF ! compr type

    DEALLOCATE(Z,W,Q,PLQ )

  END SUBROUTINE comprcauchysvd

!!!!!!
    SUBROUTINE testaccury(rowcol,D,F,U,V,Q,M,N,Rk)
      USE BasicMM

      integer M,N,Rk
      character(len=1) rowcol
      double precision D(*),F(*),U(*),V(*),Q(M,*)

      logical CR
      integer TM,TN,TRk
      double precision Relerr
      double precision, allocatable :: D1(:),F1(:),C(:,:),Err(:,:)

!  .. External Functions
      LOGICAL     LSAME
      EXTERNAL    LSAME

      CR = lsame(rowcol, 'r')
      IF( CR ) THEN
         allocate( D1(M),F1(N),C(M,N),Err(M,N) )
         D1 = D(1:M)*D(1:M)
         F1 = F(1:N)*F(1:N)
         call Cauchy2(C,D1,F1,U,V,M,N)
         Relerr = maxval( abs(Q(Rk+1:M,1:Rk)) )
         write(*,*) 'The max val of N is ', Relerr
         
         Err(1:M-Rk,1:N) = matmul(Q(Rk+1:M,1:Rk), C(1:RK,1:N))
         Err(1:M-Rk,1:N) = C(Rk+1:M,1:N) - Err(1:M-Rk,1:N)
         Err(1:M-Rk,1:N) = abs(Err(1:M-Rk,1:N)/C(Rk+1:M,1:N) )
         Relerr = maxval( Err(1:M-Rk,1:N) )
         write(*,*) 'The max ROW rel error is ', Relerr
         
      ELSE
         TM = N
         TN = Rk
         TRk = M 
         allocate( D1(TM),F1(TN),C(TM,TN),Err(TM,TN) )
         D1 = D(1:TM)*D(1:TM)
         F1 = F(1:TN)*F(1:TN)
         call Cauchy2(C,D1,F1,U,V,TM,TN)

         Err(1:TM,1:TN-TRk) = matmul(C(1:TM,1:TRk), Q(1:TRk,TRk+1:TN) )
         Err(1:TM,1:TN-TRk) = C(1:TM,TRk+1:TN) - Err(1:TM,1:TN-TRk)
         Err(1:TM,1:TN-TRk) = abs(Err(1:TM,1:TN-TRk)/C(1:TM,TRk+1:TN) )
         Relerr = maxval( Err(1:TM,1:TN-TRk) )
         write(*,*) 'The max COL rel error is ', Relerr         

      END IF
    
      deallocate( D1,F1,C,Err )
      
    END SUBROUTINE testaccury

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
  SUBROUTINE Cauchylikesvd(A,D,F,U,V,DIFL,DIFR,PL,PU,M,N,nflops)
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: M, N
    INTEGER, INTENT(OUT) :: nflops
! .. Array Arguments ..
    INTEGER, INTENT(IN) :: PL(*),PU(*)
    DOUBLE PRECISION, INTENT(IN)  :: D(*),F(*),U(*),V(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(OUT) :: A(M,N)
!
! Purpose 
! =======
! It computes an Cauchy matrix A with dimension M-by-N.
! A(i,j) = U(i)*V(j) / ( D(i)**2 - F(j)**2 ) 
!
! Parameters
! D   (input)  DOUBLE PRECISION array, DIMENSION ( LN ), LN is the size
!     of original problem. The reason is original arrays D and F are permuted and
!     they are referenced by using their index PL and PU. The new singular values 
!     which can not be used for subtraction. D(i)-D(j) is computed by using 
!     F(i), F(j), DIFL(*) and DIFR(*). 
!
! F   (input) DOUBLE PRECISION array, DIMENSION ( LN )
!     The old singular values which are used directly. 
!
! U   (input) DOUBLE PRECISION array, DIMENSION ( N )
!     Row Generator of Cauchy matrix, which are referenced directly.
! 
! V   (input) DOUBLE PRECISON array, DIMENSION ( N )
!     Column Generator of Cauchy matrix, which are referenced directly.
!
! DIFL (input) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFL(i) = D(i)-F(i), the order of which will not change and its elements
!      are also referenced by using PL or PU. 
!
! DIFR (input) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFR(i) = D(i)-F(i+1), the order of which will not change and its elements
!      are also referenced by using PL or PU. 
! 
! PL  (input) INTEGER array, DIMENSION ( N )
!     The row permutations of generator.
!
! PU  (input) INTEGER array, DIMENSION ( N )
!     The column permutations of generator.
!
! M,N (input) INTEGER, the dimension of A, M == N
! 
! nflops (output) INTEGER 
!        Floating point operations of generating A.
!
! =============
! Modified by S.-G. Li, on Dec. 13th, 2012
!  1) Add some comments
! =======================================
!
! .. Local Scalars ..
    INTEGER           i,j,PDII,PFJJ
    DOUBLE PRECISION  FJ

! .. External Functions ..
    DOUBLE PRECISION   DLAMC3
    EXTERNAL           DLAMC3
    
    nflops = 6*M*N
    DO j = 1,N
       PFJJ = PU(j)
       FJ = -F(PFJJ)
       DO i = 1, M
          PDII = PL(i)
          IF( PDII .lt. PFJJ ) THEN
             A(i,j) = U(i)*V(j) / ( DLAMC3(F(PDII+1),FJ)+DIFR(PDII) ) / ( D(PDII)+F(PFJJ) )
          ELSE
             A(i,j) = U(i)*V(j) / ( DLAMC3(F(PDII),FJ)+DIFL(PDII) ) / ( D(PDII)+F(PFJJ) )
          END IF
       END DO
    END DO
    
  END SUBROUTINE Cauchylikesvd

!!!!!!
  SUBROUTINE CauchylikesvdU( A,D,F,U,V,DIFL,DIFR,PL,PU,M,N,nflops )
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
! It computes an Cauchy matrix A with dimension M-by-N.
! A(i,j) = U(i)*V(j) / ( D(j)**2 - F(i)**2 ). This routine is different with 
! Cauchylikesvd(B, D, F, ...) and their relationship is B = A^T. 
!
! Parameters
! D   (input)  DOUBLE PRECISION array, DIMENSION ( LN ), LN is the size
!     of original problem. The reason is that original arrays D and F would be 
!     permuted and they are referenced by using their index PL and PU. 
!     The subtraction of new singular values can not be computed directly. 
!     D(i)-D(j) is computed by using F(i), F(j), DIFL(*) and DIFR(*). 
!
! F   (input) DOUBLE PRECISION array, DIMENSION ( LN )
!     The old singular values which are used directly. 
!
! U   (input) DOUBLE PRECISION array, DIMENSION ( N )
!     Row Generator of Cauchy matrix, which are referenced directly.
! 
! V   (input) DOUBLE PRECISON array, DIMENSION ( N )
!     Column Generator of Cauchy matrix, which are referenced directly.
!
! DIFL (input) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFL(i) = D(i)-F(i), the order of which will not change and its elements
!      are also referenced by using PL or PU. Positive
!
! DIFR (input) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFR(i) = D(i)-F(i+1), the order of which will not change and its elements
!      are also referenced by using PL or PU. Negative
! 
! PL  (input) INTEGER array, DIMENSION ( N )
!     The row permutations of generator.  It relates with F.
!
! PU  (input) INTEGER array, DIMENSION ( N )
!     The column permutations of generator. It relates with D.
!
! M,N (input) INTEGER, the dimension of A, M == N
! 
! nflops (output) INTEGER 
!        Floating point operations of generating A.
!
! =============
!  Written by S.-G. Li, on Dec. 23th, 2012
!  
!  Change assumed-shape dummmy array instead of assumed-size dummy array.
!
! =======================================
!
! .. Local Scalars ..
    INTEGER           i,j,PFII,PDJJ
    DOUBLE PRECISION  FJ

    nflops = 6*M*N
    DO j = 1,N
       PDJJ = PU(j)
       FJ = F(PDJJ)
       DO i = 1, M
          PFII = PL(i)
          IF( PFII .le. PDJJ ) THEN
             A(i,j) = U(i)*V(j) / ( (FJ-F(PFII) )+DIFL(PDJJ) ) / ( D(PDJJ)+F(PFII) )
          ELSE
             A(i,j) = U(i)*V(j) / ( ( F(PDJJ+1)-F(PFII) )+DIFR(PDJJ) ) / ( D(PDJJ)+F(PFII) )
          END IF
       END DO
    END DO
    
  END SUBROUTINE CauchylikesvdU

!!!!!!
  SUBROUTINE Cauchylike2svd(A,D,F,U,V,M,N,FF,DIFL,DIFR,PL,nflops)
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: M, N
    INTEGER, INTENT(OUT) :: nflops
! .. Array Arguments ..
    INTEGER, INTENT(IN)  :: PL(*)
    DOUBLE PRECISION, INTENT(IN)  :: D(*),F(*),U(*),V(*),FF(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(OUT) :: A(M,N)
!
! Purpose
! ========
! It returns an interpolative matrix, A = [I; Z], where I is an N-by-N identity matrix and
! Z is an (M-N)-by-N Cauchy-like matrix, satisfying (D**2)*Z-Z*(F**2) = U*V. 
! Z(i,j) = U(i)V(j) / (D(i)**2-F(j)**2), D**2 = D_2 and F**2 = D_1, parts of the updated svals. 
! 
! .. Parameters ..
! A  (out) DOUBLE PRECISION array, DIMENSION( M,N )
!    The interpolative matrix, A = [I; Z].
!
! D  (in) DOUBLE PRECISION array, DIMENSION( M-N )
!    The row generators of Z, D_2
!
! F  (in) DOUBLE PRECISION array, DIMENSION( N )
!    The column generators of Z, D_1
! 
! U  (in) DOUBLE PRECISION array, DIMENSION( M-N )
!    The row generators of Z
!
! V  (in) DOUBLE PRECISION array, DIMENSION( N )
!    The column generators of Z
! 
! M  (in) INTEGER, row dimension of A
!    
! N  (in) INTEGER, column dimension of A
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( LN )
!    The original F, which is referenced by calling PL.
! 
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN )
!      The distances between D(i) and F(i), positive numbers
!
! DIFR (in) DOUBLE PRECISION array, DIMENSION( LN )
!      The distances between D(i) and F(i+1), negative numbers
!
! nflops (out) INTEGER
!
! =============
! Modified by S.-G. Li, on Dec. 24th, 2012
! ========================================
! 
! .. Local Parameters ..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO = 0.0E+0, ONE=1.0E+0 ) 
! ..
! .. Local Scalars ..
    INTEGER i,j,Rk, PDJJ, PDII
    
    Rk = M-N
    nflops = 0
    A(1:M,1:N) = ZERO
    
    IF(Rk .EQ. 0) THEN
       DO j = 1, M
          A(j,j) = ONE
       END DO
    ELSE
       DO j = 1,N
          PDJJ = PL(j)
          A(j,j) = ONE
          DO i = N+1,M
             PDII = PL(i)
             IF(PDII .LT. PDJJ) THEN
                A(i,j) = U(i-N)*V(j) / (D(i-N)+F(j)) / ( (FF(PDII+1)-FF(PDJJ) ) + DIFR(PDII)-DIFL(PDJJ) )
             ELSE
                A(i,j) = U(i-N)*V(j) / (D(i-N)+F(j)) / ( (FF(PDII) -FF(PDJJ+1) ) - DIFR(PDJJ)+DIFL(PDII) )
             END IF
          END DO
       END DO
       nflops = 6*N*(M-N)           
    END IF
    
  END SUBROUTINE Cauchylike2svd

!!!!!!
  SUBROUTINE Cauchylike2svdcolU( A,D2,D1,U,V,M,N,FF,DIFL,DIFR,PL,nflops )
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: M, N
    INTEGER, INTENT(OUT) :: nflops
! .. Array Arguments ..
    INTEGER, INTENT(IN)  :: PL(*)
    DOUBLE PRECISION, INTENT(IN)  :: D2(*),D1(*),U(*),V(*),FF(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(OUT) :: A(M,N)
!
! Purpose
! ========
! It returns an interpolative matrix, A = [I Z], where I is an M-by-M identity matrix and
! Z is an M-by-(N-M) Cauchy-like matrix, satisfying Z* (D2**2) - (D1**2)*Z = V*U. 
! Z(i,j) = V(i)*U(j) / (D2(j)**2-D1(i)**2), D2 is the second part of updated svals, and 
! D1 is the first part of updated svals.  A would be a fat matrix and is the interpolative matrix
! for an HSS block column of the left singular vector matrix. 
!
! This routine is used to construct an interpolative matrix for an HSS block column, H, which is 
! first transposed to B. B is compressed by calling rrluCauchysvdrow 
! B = PQ*[I; Z^T]*B(PQ(1:Rk),:), where Z^T satisfies (D2**2)* Z^T - Z^T*(D1**2) = U*V. 
! 
! Parameters
! ==========
! A  (out) DOUBLE PRECISION array, DIMENSION( M,N )
!    The interpolative matrix, A = [I Z]. It is a fat matrix. 
!
! D2  (in) DOUBLE PRECISION array, DIMENSION( N-M )
!    The row generators of Z, D_2
!
! D1  (in) DOUBLE PRECISION array, DIMENSION( M )
!    The column generators of Z, D_1
! 
! U  (in) DOUBLE PRECISION array, DIMENSION( N-M )
!    The row generators of Z
!
! V  (in) DOUBLE PRECISION array, DIMENSION( M )
!    The column generators of Z
! 
! M  (in) INTEGER, row dimension of A
!    
! N  (in) INTEGER, column dimension of A
! 
! FF (in) DOUBLE PRECISION array, DIMENSION( LN )
!    The original F, which is referenced by calling PL.
! 
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN )
!      The distances between D(i) and F(i), positive numbers
!
! DIFR (in) DOUBLE PRECISION array, DIMENSION( LN )
!      The distances between D(i) and F(i+1), negative numbers
!
! PL (in) INTEGER array, DIMENSION( M )
!    It is used to compute the subtraction of D2(i)-D1(j) which must be computed by 
!    using FF, DIFL and DIFR. 
!
! nflops (out) INTEGER
!
! =============
! Modified by S.-G. Li, on Dec. 24th, 2012
! ========================================
! 
! .. Local Parameters ..
    DOUBLE PRECISION  ZERO, ONE, NEGONE
    PARAMETER         ( ZERO = 0.0E+0, ONE=1.0E+0, NEGONE=-1.0E+0 ) 
! ..
! .. Local Scalars ..
    INTEGER i,j,Rk, PDJJ, PDII
    
    Rk = N-M
    nflops = 0
    A(1:M,1:N) = ZERO
    
    IF(Rk .EQ. 0) THEN
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
                A(i,j) = NEGONE* V(i)*U(j-M) / ( D2(j-M)+D1(i) ) / ( (FF(PDII+1)-FF(PDJJ) ) + DIFR(PDII)-DIFL(PDJJ) )
             ELSE
                A(i,j) = NEGONE* V(i)*U(j-M) / ( D1(i)+D2(j-M) ) / ( (FF(PDII) -FF(PDJJ+1) ) - DIFR(PDJJ)+DIFL(PDII) )
             END IF
          END DO
       END DO
       nflops = 6*N*(M-N)           
    END IF
    
  END SUBROUTINE Cauchylike2svdcolU

!!!!!!
  SUBROUTINE Cauchylike2svdcol(A,D,F,U,V,M,N,nflops)
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: M, N
    INTEGER, INTENT(OUT) :: nflops
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN)  :: D(*),F(*),U(*),V(*)
    DOUBLE PRECISION, INTENT(OUT) :: A(M,N)
! 
! Purpose
! ========
! It returns an interpolative matrix, A=[I Z] where I is an N-by-N identity matrix and Z
! is an M-by-(N-M). Z is a Cauchy-like matrix defined as (F**2) * Z - Z* (D**2) = V*U.
! N >= M, A would be a fat matrix. 
! 
! The different with Cauchylike2svdrow is that both D and F are the original old singular values
! which can be computed directly. That is why we do not use DIFL and DIFR in this routine. 
! The column compression is derived by considering its transpose, (D**2) * Z^T - Z^T* (F**2) = U*V.
!
! Parameters
! ==========
! A  (out) DOUBLE PRECISION array, DIMENSION( M, N )
!    The interpolative matrix, a fat matrix, A = [I Z ]
!
! D  (in) DOUBLE PRECISION array, DIMENSION( N-M )
!    The column generators of Z, the second part of original F
! 
! F  (in) DOUBLE PRECISION array, DIMENSION( M )
!    The row generators of Z, the fist part of original F
! 
! U  (in) DOUBLE PRECISION array, DIMENSION( N-M )
!    The column generators of Z
! 
! V  (in) DOUBLE PRECISION array, DIMENSION( N-M )
!    The row generators of Z
! 
! M  (in) INTEGER, row dimension of A
! 
! N  (in) INTEGER, column dimension of A, N >= M
! 
! nflops (in) INTEGER 
!
! ============
! Modified by S.-G. Li, on Dec. 24th, 2012
! =======================================
!
! .. Parameters ..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO = 0.0E+0, ONE=1.0E+0 ) 
! .. Local Variables    
    INTEGER  i,j,Rk
    
    Rk = N - M
    nflops = 0
    A(1:M,1:N) = ZERO
    
    IF(Rk .EQ. 0) THEN
       DO j = 1, M
          A(j,j) = ONE
       END DO
    ELSE
       DO j = 1, M
          A(j,j) = ONE
       END DO
       DO j = M+1, N
          DO i = 1, M
             A(i,j) = V(i)*U(j-M) / (D(j-M)-F(i)) / ( D(j-M)+F(i) )
          END DO
       END DO

       nflops = nflops + 5*M*(N-M)
    END IF
        
  END SUBROUTINE Cauchylike2svdcol

!!!!!!
  SUBROUTINE Cauchylike2svdrowU(A,D2,D1,U,V,M,N,nflops)
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: M, N
    INTEGER, INTENT(OUT) :: nflops
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN)  :: D2(*),D1(*),U(*),V(*)
    DOUBLE PRECISION, INTENT(OUT) :: A(M,N)
! 
! Purpose
! ========
! It returns an interpolative matrix, A=[I; Z] where I is an N-by-N identity matrix and Z
! is an (M-N)-by-N. Z is a Cauchy-like matrix defined as (D2**2) * Z - Z* (D1**2) = U*V.
! M >= N, and A would be a tall and skiny matrix. 
!
! This routine is for compression of HSS block row of the left singular vector matrix. 
! The different with Cauchylike2svdrow is that both D2 and D1 are the original F 
! which can be computed directly. That is why we do not use DIFL and DIFR in this routine. 
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
! Modified by S.-G. Li, on Dec. 24th, 2012
! =======================================
!
! .. Parameters ..
    DOUBLE PRECISION  ZERO, ONE
    PARAMETER         ( ZERO = 0.0E+0, ONE=1.0E+0 ) 
! .. Local Variables    
    INTEGER  i,j,Rk
    
    Rk = N - M
    nflops = 0
    A(1:M,1:N) = ZERO
    
    IF(Rk .EQ. 0) THEN
       DO j = 1, N
          A(j,j) = ONE
       END DO
    ELSE
       DO j = 1, N
          A(j,j) = ONE
          DO i = N+1, M
             A(i,j) = U(i-N)*V(j) / (D2(i-N)-D1(j)) / ( D2(i-N)+D1(j) )
          END DO
       END DO
       nflops = nflops + 5*M*(N-M)
    END IF
        
  END SUBROUTINE Cauchylike2svdrowU

!!!!!!!!!
  SUBROUTINE Cauchy2hssvd( D,F,U,V,DIFL,DIFR,LDU,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh,time,nflops,nswap )
    USE aux_hss
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDU, LTR, lm, lvl
    INTEGER, INTENT(OUT) :: pnh, nflops
    INTEGER, INTENT(INOUT) :: nswap
    DOUBLE PRECISION, INTENT(IN)  :: TOL
    DOUBLE PRECISION, INTENT(OUT) :: time
!
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN)    :: D(*),F(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: V(*),U(*) 
    DOUBLE PRECISION, INTENT(INOUT) :: H(*), DD(*)
    INTEGER, INTENT(IN)             :: TR(*), M(*)
    TYPE(HSSMM), INTENT(INOUT)      :: PH
! 
! Purpose
! =========
! Construct an HSS matrix approximation to the right singular vectors of a broken-arrow matrix 
! M = [Z; DIAG(F)]. The diagonals of M are F, and the appended row is Z. The entries of D are 
! the singular values of M. The SVD of M is M = U * DIAG( D ) * VT. This routine returns an HSS 
! matrix approximation to VT. The matrix VT is defined as A(i,j) = U(i)*V(j)/(D(i)**2-F(j)**2), 
! where D is the updated singular values and F is the old singular values, D*A-A*F=u*v. 
!
! DIFL(i) is the distance between D(i) and F(i), DIFL(i) = D(i) - F(i), a positive value.
! DIFR(i) is the distance from D(i) to F(i+1), DIFR(i) = D(i)-F(i+1), a negative value. 
! This subroutine can be used for updating SVD problem.
! 
! Parameters
! ========== 
! D   (in)  DOUBLE PRECISION array, DIMENSION( N ). 
!     The newly computed singular values, and its entries can not be used directly.
!     The row generators of A
!
! F   (in) DOUBLE PRECISION array, DIMENSION( N )
!     The old singular values, and the column generator of A. 
!
! U   (inout)  DOUBLE PRECISION array, DIMENSION( N )
!     The row generators of A, the reciprocal of the 2-norm of each column of A
!
! V   (inout) DOUBLE PRECISION array, DIMENSION( N )
!     The column generators of A, which are V(i) = d(i)*z(i), see our paper. 
!
! LDU   (in) INTEGER, the leading dimension of A
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
! H     (out) DOUBLE PRECISION array, DIMENSION( LDA *Mi * alpha )
!        alpha is a constant(=7); Mi is size of bottom HSS block, in the order of 
!       the HSS rank; 
! 
! TOL   (in) DOUBLE PRECISION, tolerance for low-rank approximation
!
! lvl   (in) INTEGER, the total level of HSS tree.
! 
! pnh   (out) INTEGER, the total storage of H which stores U, B and V^T.
!
! time  (out) DOUBLE PRECISION, the accumulated time of calling rrluCauchy.
!
! nflops (out) INTEGER, the total flops of constructing this HSS matrix 
!
! nswap  (inout) INTEGER, the total swaps whiling constructing this HSS matrix
!        If the entry of nswap = -1, then this is used for updating SVD problem.
!
!===================
!  Written by Shengguo Li, on Dec. 15th, 2012
!  Construct HSS matrix from Cauchy-like matrix by using Structured RRLU factorization
!  for the right singular vector matrix. 
!  Modified on April 11th, 2013
!  1) Remove the use of D, and use Di directly.
!=====================================================================================
!
! .. Parameters ..
    DOUBLE PRECISION :: ZERO, ONE,time1
    PARAMETER        ( ZERO=0.0E+0, ONE=1.0E+0 )
! .. Local Scalars ..
    INTEGER       :: ch1,ierr,info,nn,mi,ni,nt,Rk,it
    INTEGER       :: pnd,i,j,lt,ns,lflops,lswap
! .. Local Arrays ..
    INTEGER, ALLOCATABLE :: ch(:,:),l(:,:),PLi(:),PUi(:)
    INTEGER, ALLOCATABLE :: lsr(:),lsc(:),PL(:),PU(:)
    DOUBLE PRECISION, ALLOCATABLE :: Di(:),Fi(:),Ui(:),Vi(:)
    ! lsr is the pointer of block row in stack A; lsc is the pointer of block column in stack A
! ..
! .. External Subroutines ..
    EXTERNAL  DGEMM, DLACPY
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MIN,ABS,MAXVAL,SQRT
              
    nflops = 0
    pnh = 1               ! 'pointer' in H
    pnd = 1               ! 'pointer' in D
!    ALLOCATE( PL(ldu),PU(ldu),ch(2,ltr),l(ltr,2),lsc(lvl+2),lsr(lvl+2),stat=ierr )
    ALLOCATE( PL(ldu),PU(ldu),ch(2,ltr),l(ltr,2),lsc(ltr),lsr(ltr),stat=ierr )
    IF(ierr /= 0 ) THEN
       WRITE(*,*) "Allocate failed in cauchy2hss! "
       RETURN
    END IF

    ALLOCATE( Di(ldu),Fi(ldu),Ui(ldu),Vi(ldu),PLi(ldu),PUi(ldu), stat=ierr )
    IF(ierr /= 0 ) THEN
       WRITE(*,*) "Allocate failed in cauchy2hss! 2 "
       RETURN
    END IF

    nn = ltr
    PL(1:ldu) = (/ (j, j=1,ldu) /)
    PU(1:ldu) = (/ (j, j=1,ldu) /)
    call child(tr, ch, ltr)

    l(1,1:2) = (/1,m(1)/)
    lt = 1
    it = 1
    DO i = 1, ltr
       IF (ch(1,i) == 0 .and. ch(2,i) == 0) THEN
          l(i,1:2) = (/lt,lt+m(it)-1/)
          lt = l(i,2)+1
          it = it+1
       ELSE
          l(i,1:2) = (/l(ch(1,i),1), l(ch(2,i),2)/)
       END IF
    END DO

!     U is alpha, D is the new svals, V is the updated Z, F is old svals
!     The following four lines are for updating SVD problems. Uncomment them
!     IF you are solving updating SVD problem.
    IF(nswap .eq. -1 ) THEN
       PL(1:ldu) = (/ (j, j=ldu,1,-1) /)
       PU(1:ldu) = (/ (j, j=ldu,1,-1) /)
       V(1:ldu) = V(PL)
       U(1:ldu) = U(PL)
    END IF

    nswap = 0
    pnh = 1               ! 'pointer' in H
    ns = 1                  ! (# of blocks in stack)+1
    lsr = 0
    lsc = 0
    lsr(ns) = 1             ! row start position at upper triangular
    lsc(ns) = 1             ! column start position at lower triangular
    time = zero

! *****************************************************
!                 MAIN LOOP                           *
!******************************************************
    DO i = 1, nn

       ! leaf node
       IF( ch(1,i) == 0 ) THEN
          mi = l(i,2)-l(i,1)+1
          ni = ldu-l(i,2)  !t^r, the right column block's length

          call Cauchylikesvd( DD(pnd),D,F,U(l(i,1)),V(l(i,1) ),DIFL,DIFR,PL(l(i,1)),PU(l(i,1)),mi,mi,lflops)
          PH%D(i) = mi
          PH%pd(1,i) = pnd
          PH%pd(2,i) = pnd + mi*mi -1
          pnd = pnd + mi*mi
          nflops = nflops + lflops

          ! off-diag row compression
          IF( ns .EQ. 1 ) THEN
             nt = ni
             Vi(1:nt ) = V(l(i,2)+1:ldu)
             PUi(1:nt) = PU(l(i,2)+1:ldu)
          ELSE
             nt = ni + lsc(ns)-1
             Vi(1:lsc(ns)-1) = V(1:lsc(ns)-1)
             Vi(lsc(ns):nt ) = V(l(i,2)+1:ldu) 
             PUi(1:lsc(ns)-1) = PU(1:lsc(ns)-1)
             PUi(lsc(ns):nt ) = PU(l(i,2)+1:ldu) 
          END IF

          Ui(1:mi) = U(l(i,1):l(i,2))
          PLi(1:mi)= PL(l(i,1):l(i,2))
          Di(1:mi) = D(PLi(1:mi))
          Fi(1:nt) = F(PUi(1:nt))
          call Comprcauchysvd('r',Di,Fi,Ui,Vi,DIFL,DIFR,F,PLi,PUi,tol,mi,nt,PH,pnh,H,i,Rk,time1,lflops,lswap)
          time = time1+time
          nflops = nflops + lflops
          nswap = nswap + lswap

          lsr(ns+1) = lsr(ns)+Rk
          U(lsr(ns):lsr(ns+1)-1)  = Ui(1:Rk)
          PL(lsr(ns):lsr(ns+1)-1) = PLi(1:Rk)

          ! off-diag col compression
          IF( ns .EQ. 1 ) THEN
             nt = ni
             Ui(1:nt ) = U(l(i,2)+1:ldu )
             PLi(1:nt) = PL(l(i,2)+1:ldu)
          ELSE
             nt = ni + lsr(ns)-1
             Ui(1:lsr(ns)-1) = U( 1:lsr(ns)-1 )
             Ui(lsr(ns):nt ) = U( l(i,2)+1:ldu)
             PLi(1:lsr(ns)-1) = PL( 1:lsr(ns)-1 )
             PLi(lsr(ns):nt ) = PL( l(i,2)+1:ldu)
          END IF

          Vi(1:mi) = V( l(i,1):l(i,2) )
          PUi(1:mi)= PU( l(i,1):l(i,2) )
          Di(1:nt) = D(PLi(1:nt))
          Fi(1:mi) = F(PUi(1:mi))
          call Comprcauchysvd('c',Di,Fi,Ui,Vi,DIFL,DIFR,F,PLi,PUi,tol,nt,mi,PH,pnh,H,i,Rk,time1,lflops,lswap)
          time = time1+time
          nflops = nflops + lflops
          nswap = nswap + lswap

          lsc(ns+1) = lsc(ns)+Rk
          V(lsc(ns):lsc(ns+1)-1)  = Vi(1:Rk)
          PU(lsc(ns):lsc(ns+1)-1) = PUi(1:Rk)
          ns = ns+1

       ELSE ! parent nodes
          ch1 = ch(1,i)
          mi  = lsr(ns-1) - lsr(ns-2)
          ni  = lsc(ns)   - lsc(ns-1)
          call Cauchylikesvd(H(pnh),D,F,U(lsr(ns-2)),V(lsc(ns-1)),DIFL,DIFR,PL(lsr(ns-2)),PU(lsc(ns-1)),mi,ni,lflops)
          call hssexpmm2('B',PH,pnh,mi,ni,ch1,info)  ! B{ch{tr{i}}(1)}
          nflops = nflops + lflops

          ch1 = ch(2,i)
          mi  = lsr(ns) - lsr(ns-1)
          ni  = lsc(ns-1) - lsc(ns-2)
          call Cauchylikesvd(H(pnh),D,F,U(lsr(ns-1)),V(lsc(ns-2)),DIFL,DIFR,PL(lsr(ns-1)),PU(lsc(ns-2)),mi,ni,lflops)
          call hssexpmm2('B',PH,pnh,mi,ni,ch1,info)  ! B{ch{tr{i}}(2)}
          nflops = nflops + lflops

          IF(i .EQ. nn) exit

          ! off-diag row compression
          mi = lsr(ns) - lsr(ns-2)
          ni = ldu-l(i,2)    !t^r
          IF(ns .EQ. 3) THEN
             nt = ni
             Vi( 1:nt ) = V( l(i,2)+1:ldu )  ! col
             PUi(1:nt ) = PU( l(i,2)+1:ldu )  
          ELSE
             nt = ni + lsc(ns-2)-1
             Vi(1:lsc(ns-2)-1 ) = V(1:lsc(ns-2)-1) 
             Vi( lsc(ns-2):nt ) = V( l(i,2)+1:ldu )
             PUi(1:lsc(ns-2)-1 ) = PU(1:lsc(ns-2)-1) 
             PUi( lsc(ns-2):nt ) = PU( l(i,2)+1:ldu )
          END IF

          Ui(1:mi) = U( lsr(ns-2):lsr(ns)-1 )
          PLi(1:mi)= PL( lsr(ns-2):lsr(ns)-1 )
          Di(1:mi) = D(PLi(1:mi))
          Fi(1:nt) = F(PUi(1:nt))
          call Comprcauchysvd('r',Di,Fi,Ui,Vi,DIFL,DIFR,F,PLi,PUi,tol,mi,nt,PH,pnh,H,i,Rk,time1,lflops,lswap) 
          time = time1+time
          nflops = nflops + lflops
          nswap = nswap + lswap

          lsr(ns-1) = lsr(ns-2) + Rk
          U(lsr(ns-2):lsr(ns-1)-1) = Ui(1:Rk)
          PL(lsr(ns-2):lsr(ns-1)-1) = PLi(1:Rk)

          ! off-diag col compression
          mi = lsc(ns) - lsc(ns-2)
          ni = ldu-l(i,2)     !t^r
          IF(ns .EQ. 3 ) THEN
             nt = ni
             Ui(1:nt ) = U( l(i,2)+1:ldu )
             PLi(1:nt) = PL(l(i,2)+1:ldu)
          ELSE
             nt = ni + lsr(ns-2)-1
             Ui(1:lsr(ns-2)-1 ) = U(1:lsr(ns-2)-1 )
             Ui(lsr(ns-2):nt ) = U(l(i,2)+1:ldu)
             PLi(1:lsr(ns-2)-1 ) = PL(1:lsr(ns-2)-1 )
             PLi(lsr(ns-2):nt ) = PL(l(i,2)+1:ldu)
          END IF

          Vi(1:mi) = V( lsc(ns-2):lsc(ns)-1 ) 
          PUi(1:mi)= PU( lsc(ns-2):lsc(ns)-1 ) 
          Di(1:nt) = D(PLi(1:nt))
          Fi(1:mi) = F(PUi(1:mi))
          call Comprcauchysvd('c',Di,Fi,Ui,Vi,DIFL,DIFR,F,PLi,PUi,tol,nt,mi,PH,pnh,H,i,Rk,time1,lflops,lswap)
          time = time1+time
          nflops = nflops + lflops
          nswap = nswap + lswap

          lsc(ns-1) = lsc(ns-2)+Rk
          V(lsc(ns-2):lsc(ns-1)-1) = Vi(1:Rk)
          PU(lsc(ns-2):lsc(ns-1)-1) = PUi(1:Rk)
          ns = ns-1
       END IF

    END DO ! main loop

    DEALLOCATE( ch,l,lsc,lsr,Di,Fi,Ui,Vi,PLi,PUi )

  END SUBROUTINE Cauchy2hssvd

!!!!!!!!!
  SUBROUTINE Cauchy2hssvdU(D,F,U,V,DIFL,DIFR,LDU,TR,LTR,M,LM,PH,H,DD,TOL,lvl,pnh,time,nflops,nswap)
    USE aux_hss
!
!  .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDU, LTR, lm, lvl
    INTEGER, INTENT(OUT) :: pnh, nflops
    INTEGER, INTENT(INOUT) :: nswap
    DOUBLE PRECISION, INTENT(IN)  :: TOL
    DOUBLE PRECISION, INTENT(OUT) :: time
!
!  .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN)    :: D(*),F(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: U(*),V(*) 
    DOUBLE PRECISION, INTENT(INOUT) :: H(*), DD(*)
    INTEGER, INTENT(IN)             :: TR(*), M(*)
    TYPE(HSSMM), INTENT(INOUT)      :: PH
!
! Purpose
! =========
! Construct an HSS matrix approximation to the left singular vectors of a broken-arrow matrix M = [Z; DIAG(F)].
! The diagonals of M are F, and the appended row is Z. The entries of D are the singular values of M.
! The SVD of M is M = U * DIAG( D ) * VT.
! 
! This routine returns an HSS approximation to Cauchy-like matrix A( i,j ) = U(i)*V(j) / (D(j)**2 - F(i)**2), 
! and is for structured divide and conquer algorithm. 
!
! This subroutine can be used for updating SVD problems, uncommment the four lines below. 
! 
! Parameters
! ========== 
! D   (in)  DOUBLE PRECISION array, DIMENSION( N ). 
!     Its entries are the newly computed singular values and the column generators of A. 
!
! F   (in) DOUBLE PRECISION array, DIMENSION( N )
!     Its entries are the old singular values and the row generators of A. 
!
! U   (inout) DOUBLE PRECISION array, DIMENSION( N )
!     The row generators of A, which are U(i) = d(i)*z(i), see our paper. 
!
! V   (inout)  DOUBLE PRECISION array, DIMENSION( N )
!     The column generators of A, the reciprocal of the 2-norm of each column of A
!
! DIFL  (in) DOUBLE PRECISION array, DIMENSION( N )
!       DIFL(i) is the distance between D(i) and F(i), DIFL(i) = D(i) - F(i), a positive value.
!
! DIFR  (in) DOUBLE PRECISION array, DIMENSION( N )
!       DIFR(i) is the distance from D(i) to F(i+1), DIFR(i) = D(i)-F(i+1), a negative value. 
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
! lvl   (in) INTEGER, the total level of HSS tree.
! 
! pnh   (out) INTEGER, the total storage of H which stores U, B and V^T.
!
! time  (out) DOUBLE PRECISION, the accumulated time of calling rrluCauchy.
!
! nflops (out) INTEGER, the total flops of constructing this HSS matrix 
!
! nswap  (inout) INTEGER, the total swaps whiling constructing this HSS matrix
!        If nswap = -1, it is for updating SVD problem
!
! ===================
!  Written by Shengguo Li, on Dec. 15th, 2012
!  Construct HSS matrix from Cauchy-like matrix by using Structured RRLU factorization
!  for the left singular vector matrix. 
!  
! =====================================================================================
!
! .. Parameters 
    DOUBLE PRECISION :: ZERO, ONE
    PARAMETER        ( ZERO=0.0E+0, ONE=1.0E+0 )
! ..
! .. Local Scalars ..
    INTEGER       :: ch1,ierr,info,nn,mi,ni,nt,Rk,it
    INTEGER       :: pnd,i,j,lt,ns,lflops,lswap
    DOUBLE PRECISION  time1
! ..
! .. Local Arrays ..
    INTEGER, ALLOCATABLE :: ch(:,:),l(:,:),PLi(:),PUi(:)
    INTEGER, ALLOCATABLE :: lsr(:),lsc(:),PL(:),PU(:)
    DOUBLE PRECISION, ALLOCATABLE :: Di(:),Fi(:),Ui(:),Vi(:)
    ! lsr is the pointer of block row in stack A; lsc is the pointer of block column in stack A
! ..
! .. External Subroutines ..
    EXTERNAL  DGEMM, DLACPY
! ..
! .. Intrinsic Functions ..
    INTRINSIC    MIN,ABS,MAXVAL,SQRT
    
    nflops = 0
    nswap = 0
    pnh = 1               ! 'pointer' in H
    pnd = 1               ! 'pointer' in D
!    ALLOCATE( PL(ldu),PU(ldu),ch(2,ltr),l(ltr,2),lsc(lvl+2),lsr(lvl+2),stat=ierr )
    ALLOCATE( PL(ldu),PU(ldu),ch(2,ltr),l(ltr,2),lsc(ltr),lsr(ltr),stat=ierr )
    IF(ierr /= 0 ) THEN
       WRITE(*,*) "Allocate failed in cauchy2hss! "
       RETURN
    END IF

    ALLOCATE( Di(ldu),Fi(ldu),Ui(ldu),Vi(ldu),PLi(ldu),PUi(ldu), stat=ierr )
    IF(ierr /= 0 ) THEN
       WRITE(*,*) "Allocate failed in cauchy2hss! 2 "
       RETURN
    END IF

    nn = ltr
    PL(1:ldu) = (/ (j, j=1,ldu) /)
    PU(1:ldu) = (/ (j, j=1,ldu) /)
    call child(tr, ch, ltr)

    l(1,1:2) = (/1,m(1)/)
    lt = 1
    it = 1
    DO i = 1, ltr
       IF (ch(1,i) == 0 .and. ch(2,i) == 0) THEN
          l(i,1:2) = (/lt,lt+m(it)-1/)
          lt = l(i,2)+1
          it = it+1
       ELSE
          l(i,1:2) = (/l(ch(1,i),1), l(ch(2,i),2)/)
       END IF
    END DO

!     U is alpha, D is the new svals, V is the updated Z, F is old svals
!     The following four lines are for updating SVD problems. Uncomment them
!     IF you are solving updating SVD problem.
    IF(nswap .eq. -1 ) THEN
       PL(1:ldu) = (/ (j, j=ldu,1,-1) /)
       PU(1:ldu) = (/ (j, j=ldu,1,-1) /)
       V(1:ldu) = V(PL)
       U(1:ldu) = U(PL)
    END IF

    ns = 1                  ! (# of blocks in stack)+1
    lsr = 0
    lsc = 0
    lsr(ns) = 1             ! row start position at upper triangular
    lsc(ns) = 1             ! column start position at lower triangular
    time = zero
       
! *****************************************************
!                 MAIN LOOP                           *
!******************************************************
    DO i = 1, nn

       ! leaf node
       IF( ch(1,i) == 0 ) THEN
          mi = l(i,2)-l(i,1)+1
          ni = ldu-l(i,2)  !t^r, the right column block's length

          call CauchylikesvdU( DD(pnd),D,F,U(l(i,1)),V(l(i,1)),DIFL,DIFR,PL(l(i,1)),PU(l(i,1)),mi,mi,lflops )
          PH%D(i) = mi
          PH%pd(1,i) = pnd
          PH%pd(2,i) = pnd + mi*mi -1
          pnd = pnd + mi*mi
          nflops = nflops + lflops

          ! off-diag row compression
          IF( ns .EQ. 1 ) THEN
             nt = ni
             Vi(1:nt ) = V(l(i,2)+1:ldu)
             PUi(1:nt) = PU(l(i,2)+1:ldu)
          ELSE
             nt = ni + lsc(ns)-1
             Vi(1:lsc(ns)-1) = V(1:lsc(ns)-1)
             Vi(lsc(ns):nt ) = V(l(i,2)+1:ldu)
             PUi(1:lsc(ns)-1) = PU(1:lsc(ns)-1)
             PUi(lsc(ns):nt ) = PU(l(i,2)+1:ldu)
          END IF

          Ui(1:mi) = U(l(i,1):l(i,2))
          PLi(1:mi)= PL(l(i,1):l(i,2))
          Fi(1:mi) = F(PLi(1:mi))
          Di(1:nt) = D(PUi(1:nt))
          call ComprcauchysvdU('r',Di,Fi,Ui,Vi,DIFL,DIFR,F,PLi,PUi,tol,mi,nt,PH,pnh,H,i,Rk,time1,lflops,lswap)
          time = time1+time
          nflops = nflops + lflops
          nswap = nswap + lswap

          lsr(ns+1) = lsr(ns)+Rk
          U(lsr(ns):lsr(ns+1)-1)  = Ui(1:Rk)
          PL(lsr(ns):lsr(ns+1)-1) = PLi(1:Rk)

          ! off-diag col compression
          IF( ns .EQ. 1 ) THEN
             nt = ni
             Ui(1:nt ) = U(l(i,2)+1:ldu )
             PLi(1:nt) = PL(l(i,2)+1:ldu)
          ELSE
             nt = ni + lsr(ns)-1
             Ui(1:lsr(ns)-1) = U( 1:lsr(ns)-1 )
             Ui(lsr(ns):nt ) = U( l(i,2)+1:ldu)
             PLi(1:lsr(ns)-1) = PL( 1:lsr(ns)-1 )
             PLi(lsr(ns):nt ) = PL( l(i,2)+1:ldu)
          END IF

          Vi(1:mi) = V( l(i,1):l(i,2) )
          PUi(1:mi)= PU( l(i,1):l(i,2) )
          Fi(1:nt) = F(PLi(1:nt))
          Di(1:mi) = D(PUi(1:mi))
          call ComprcauchysvdU('c',Di,Fi,Ui,Vi,DIFL,DIFR,F,PLi,PUi,tol,nt,mi,PH,pnh,H,i,Rk,time1,lflops,lswap)
          time = time1+time
          nflops = nflops + lflops
          nswap = nswap + lswap

          lsc(ns+1) = lsc(ns)+Rk
          V(lsc(ns):lsc(ns+1)-1)  = Vi(1:Rk)
          PU(lsc(ns):lsc(ns+1)-1) = PUi(1:Rk)
          ns = ns+1

       ELSE ! parent nodes
          ch1 = ch(1,i)
          mi  = lsr(ns-1) - lsr(ns-2)
          ni  = lsc(ns)   - lsc(ns-1)
          call CauchylikesvdU(H(pnh),D,F,U(lsr(ns-2)),V(lsc(ns-1)),DIFL,DIFR,PL(lsr(ns-2)),PU(lsc(ns-1)),mi,ni,lflops)
          call hssexpmm2('B',PH,pnh,mi,ni,ch1,info)  ! B{ch{tr{i}}(1)}
          nflops = nflops + lflops

          ch1 = ch(2,i)
          mi  = lsr(ns) - lsr(ns-1)
          ni  = lsc(ns-1) - lsc(ns-2)
          call CauchylikesvdU(H(pnh),D,F,U(lsr(ns-1)),V(lsc(ns-2)),DIFL,DIFR,PL(lsr(ns-1)),PU(lsc(ns-2)),mi,ni,lflops)
          call hssexpmm2('B',PH,pnh,mi,ni,ch1,info)  ! B{ch{tr{i}}(2)}
          nflops = nflops + lflops

          IF(i .EQ. nn) exit

          ! off-diag row compression
          mi = lsr(ns) - lsr(ns-2)
          ni = ldu-l(i,2)    !t^r
          IF(ns .EQ. 3) THEN
             nt = ni
             Vi( 1:nt ) = V( l(i,2)+1:ldu )  ! col
             PUi(1:nt ) = PU( l(i,2)+1:ldu )  
          ELSE
             nt = ni + lsc(ns-2)-1
             Vi(1:lsc(ns-2)-1 ) = V(1:lsc(ns-2)-1) 
             Vi( lsc(ns-2):nt ) = V( l(i,2)+1:ldu )
             PUi(1:lsc(ns-2)-1 ) = PU(1:lsc(ns-2)-1) 
             PUi( lsc(ns-2):nt ) = PU( l(i,2)+1:ldu )
          END IF

          Ui(1:mi) = U( lsr(ns-2):lsr(ns)-1 )
          PLi(1:mi)= PL( lsr(ns-2):lsr(ns)-1 )
          Fi(1:mi) = F(PLi(1:mi))
          Di(1:nt) = D(PUi(1:nt))
          call ComprcauchysvdU('r',Di,Fi,Ui,Vi,DIFL,DIFR,F,PLi,PUi,tol,mi,nt,PH,pnh,H,i,Rk,time1,lflops,lswap) 
          time = time1+time
          nflops = nflops + lflops
          nswap = nswap + lswap

          lsr(ns-1) = lsr(ns-2) + Rk
          U(lsr(ns-2):lsr(ns-1)-1) = Ui(1:Rk)
          PL(lsr(ns-2):lsr(ns-1)-1) = PLi(1:Rk)

          ! off-diag col compression
          mi = lsc(ns) - lsc(ns-2)
          ni = ldu-l(i,2)     !t^r
          IF(ns .EQ. 3 ) THEN
             nt = ni
             Ui(1:nt ) = U( l(i,2)+1:ldu )
             PLi(1:nt) = PL(l(i,2)+1:ldu)
          ELSE
             nt = ni + lsr(ns-2)-1
             Ui(1:lsr(ns-2)-1 ) = U(1:lsr(ns-2)-1 )
             Ui(lsr(ns-2):nt ) = U(l(i,2)+1:ldu)
             PLi(1:lsr(ns-2)-1 ) = PL(1:lsr(ns-2)-1 )
             PLi(lsr(ns-2):nt ) = PL(l(i,2)+1:ldu)
          END IF

          Vi(1:mi) = V( lsc(ns-2):lsc(ns)-1 ) 
          PUi(1:mi)= PU( lsc(ns-2):lsc(ns)-1 ) 
          Fi(1:nt) = F(PLi(1:nt))
          Di(1:mi) = D(PUi(1:mi))
          call ComprcauchysvdU('c',Di,Fi,Ui,Vi,DIFL,DIFR,F,PLi,PUi,tol,nt,mi,PH,pnh,H,i,Rk,time1,lflops,lswap)
          time = time1+time
          nflops = nflops + lflops
          nswap = nswap + lswap

          lsc(ns-1) = lsc(ns-2)+Rk
          V(lsc(ns-2):lsc(ns-1)-1) = Vi(1:Rk)
          PU(lsc(ns-2):lsc(ns-1)-1) = PUi(1:Rk)
          ns = ns-1
       END IF

    END DO ! main loop

    DEALLOCATE( ch,l,lsc,lsr,Di,Fi,Ui,Vi,PLi,PUi )

  END SUBROUTINE Cauchy2hssvdU

!!!!!!
  SUBROUTINE ComprcauchysvdU(rowcol,D,F,U,V,DIFL,DIFR,FF,PL,PU,tol,M,N,PH,pnh,H,nodi,Rk,time,nflops,nswap)
    USE aux_hss
!
! Scalar parameters
    INTEGER, INTENT(IN)    :: M,N,nodi
    INTEGER, INTENT(OUT)   :: Rk,nflops,nswap
    INTEGER, INTENT(INOUT) :: pnh
    DOUBLE PRECISION, INTENT(IN)  :: tol
    DOUBLE PRECISION, INTENT(OUT) :: time
    CHARACTER(LEN=1), INTENT(IN)  :: rowcol
! Array parameters
    DOUBLE PRECISION, INTENT(IN)    :: DIFL(*),DIFR(*),FF(*)
    DOUBLE PRECISION, INTENT(INOUT) :: D(*),F(*),U(*),V(*),H(*)
    INTEGER, INTENT(INOUT)     :: PL(M),PU(N)
    TYPE(HSSMM), INTENT(INOUT) :: PH
!
! Purpose
! ========
! This routine is written for computing the left singular vector matrix. Here we are assume that
! C is an M-by-N matrix. 
! For row compression, it computes a low rank approximation of some off-diagonal block of a
! Cauchy-like matrix C(i,j)= U(i)*V(j) / ( D(j)**2 - F(i)**2 ), where D are the updated singular values
! and F are the old ones. C is factorized as 
! C = PQ*[I; Mat_N]* C(PQ(1:Rk),:), where PQ is a permutation and Mat_N is also a Cauchy-like matrix sastifying
! F(Rk+1:M) * Mat_N - Mat_N*F(1:Rk) = Z*W. 
!
! For column compression, it computes a low-rank approximation of a Cauchy-like matrix 
! C(i,j) = U(i)*V(j) / (D(j)**2 - F(i)**2), where D are the updated singular values and F are the old ones. 
! We get a low-rank approximation of C as C = C(:,PQ(1:Rk) ) * [I Z] * PQ, where
! Z satisfies Z* (D2**2) - (D1**2) * Z = W(1:Rk)*Z(Rk+1:N) and D2=D(Rk+1:N) and D1=D(1:Rk). 
! 
!
! The difference between this one and Comprcauchysvd is that one is another's transpose. 
!
! Parameters
! ==========
! rowcol (in) CHARACTER ( LEN = 1 )
!        = 'R', row compression
!        = 'C', column compresion
! 
! F  (in) DOUBLE PRECISION array, DIMENSION( M )
!    F are permuted in essential and only the first Rk entries are useful. 
!    But its permutation is recorded by PL and entries of F will not change. The entries of F 
!    are the old singular values and are the row generators of Cauchy-like matrix A. F and PL are in 
!    one group.
!
! D  (in)  DOUBLE PRECISION array, DIMENSION( N )
!    D will be permuted through PU in the process. Finally the permuted of D
!    will not be used since we are using interpolative decomposition and D corresponds
!    to the column. The entries of D are the new singular values. 
!
! U  (inout) DOUBLE PRECISION array, DIMENSION( M ) 
!    U is the updated Z in our old notation and it will be perumuted, only 
!    the first Rk entries are useful for row compression. 
! 
! V  (inout)  DOUBLE PRECISION array, DIMENSION( N )
!    V is the column normalization scalars, Beta. It will be permuted and 
!    the first Rk entries are useful for column compression and V will not be used in 
!    row compression.
!
! DIFL (in) DOUBLE PRECISION array, DIMENSION( LN ) 
!      DIFL(i) = D(i)-F(i), the order of which will not change and its elements
!      are also referenced by using PL or PU. Positive
! 
! DIFR (in) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFR(i) = D(i)-F(i+1), the order of which will not change and its elements
!      are also referenced by using PL or PU. Negative
!
! FF (in) DOUBLE PRECISION array, DIMENSION( LN ), original F
!    LN is the size of original problem. 
!
! PL (inout) INTEGER array, DIMENSION( M )
!    It stores the permutations for row generators of A. It corresponds to F for row compression and
!    it corresponds to F for row compression.
!
! PU (inout) INTEGER array, DIMENSION( N )
!    It stores the permutations for col generators of A. It corresponds to D for column compression and
!    it corresponds to D for column compression.
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
! time  (out) DOUBLE PRECISION, the accumulated time of calling rrluCauchy.
!
! nflops (out) INTEGER, the total flops of constructing this HSS matrix 
!
! nswap  (out) INTEGER, the total swaps whiling constructing this HSS matrix
!
! =========
! Written by S.-G. Li, on Dec. 15th, 2012
! =========
!
! .. Parameters ..
    DOUBLE PRECISION ZERO, ONE
    PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0)
! ..    
! .. Local scalars ..
    INTEGER          mn,i,info,lflops, ierr
    DOUBLE PRECISION time1
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
       call cpu_time(time)
!       write(*,*) 'Row compression, Col M ', M, ' Row N ', N
       call rrluCauchysvdcol(F,D,U,V,DIFL,DIFR,tol,Z,W,FF,PL,PU,M,N,Rk,PQ,lflops,nswap )
       nflops = nflops + lflops
       call cpu_time(time1)
       time = time1 - time
       call invp(PQ,M)  ! PL --> InvPL
       IF( Rk .LT. MN ) THEN
          ALLOCATE( Q(M,Rk) )
          call Cauchylike2svdrowU( Q,F(Rk+1),F,Z(Rk+1),W,M,Rk,lflops )
          nflops = nflops + lflops
          Q(1:M,1:Rk) = Q(PQ,1:Rk)
          call dlacpy('A',M,Rk,Q,M,H(pnh),M)      ! copy U to H
          call hssexpmm2('U',PH,pnh,M,Rk,nodi,info) ! copy Q to generators
       ELSE
          ! copy identity matrix to generators
          ALLOCATE( Q(M,M) )
          Rk = M
          Q(1:M,1:M) = ZERO
          DO i = 1,M
             Q(i,i) = ONE
          END DO
          Q(1:M,1:M) = Q(PQ,1:M)            
          call dlacpy('A',M,M,Q,M,H(pnh),M)         ! copy U to H
          call hssexpmm2('U',PH,pnh,Rk,M,nodi,info) ! copy Q to generators
       END IF

    ELSE ! block column

       ALLOCATE( Z(N),W(MN),PQ(N), stat = ierr )    ! Check ??
       Z = ZERO
       W = ZERO
       call cpu_time(time)
!       write(*,*) 'Col compression, Col M ', M, ' Row N ', N
       call rrluCauchysvdrow( D,F,V,U,DIFL,DIFR,tol,Z,W,FF,PU,PL,N,M,Rk,PQ,lflops,nswap )
       nflops = nflops + lflops
       call cpu_time(time1)
       time = time1 - time
       call invp(PQ,N)       ! PL --> InvPL
       IF( Rk .LT. MN ) THEN
          ALLOCATE( Q(Rk, N) )
          call Cauchylike2svdcolU(Q,D(Rk+1),D,Z(Rk+1),W,Rk,N,FF,DIFL,DIFR,PU,lflops)
          nflops = nflops + lflops
          Q(1:Rk,1:N) = Q(1:Rk,PQ)            
          call dlacpy('A',Rk,N,Q,Rk,H(pnh),Rk)        ! copy V to H
          call hssexpmm2('V',PH,pnh,Rk,N,nodi,info) ! copy Q to generators            
       ELSE
          ! copy identity matrix to generators
          allocate( Q(N,N) )
          Rk = N
          Q(1:N,1:N) = ZERO
          DO i = 1,N
             Q(i,i) = ONE
          END DO
          Q(1:N,1:N) = Q(1:N,PQ)
          call dlacpy('A',N,Rk,Q,N,H(pnh),N)        ! copy V to H
          call hssexpmm2('V',PH,pnh,N,Rk,nodi,info) ! copy Q to generators
       END IF

    END IF ! compr type

    DEALLOCATE(Z,W,Q,PQ )

  END SUBROUTINE ComprcauchysvdU

!!!!!!
  SUBROUTINE DHSSVCS( K, Ni, D, Z, F, DIFL, DIFR, APHAU, APHAV, ZZ, &
                      UT,LCUT,VT,LCVT,WORK,INFO )
     USE aux_hss
     USE BasicMM
     USE ConstructHssd
!
! .. Scalar Parameters ..
    INTEGER, INTENT(IN) ::  K, LCUT, LCVT, INFO, Ni
!
! .. Array Parameters ..
    DOUBLE PRECISION, INTENT(IN) :: D(*),F(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: APHAU(*), APHAV(*),WORK(*),Z(*), & 
                                       UT(LCUT,K),VT(K,LCVT),ZZ(*)
!
! Purpose
! ========
! This routine computes the singular vectors of a broken arrow matrix, M =[Z; F] with
! Z is a row vector, and F is a diagonal matrix, its entries are old singular values in
! ascending order. 
! We use HSS matrices to approximate the left and right singular vectors without forming 
! them explicitly. The SVD of M is computed as M = UH * D * VHT, where UH and VHT are HSS
! matrices. Then, the left and right singular vectors of the upper level matrix are updated 
! by using matrix-matrix multiplication, U = U*UH and VT = VHT*VT, where U is an LDU-by-K
! matrix and VT is K-by-LDVT. 
!
! More details: 
! We first transpose U and then use fasthssmm to update the left singular vector matrix.
! Whether using fasthssmmL is faster is not known.  
!
! There is no enough reason for choosing NCOL = 200
!
! Arguments  
! =========
! 
! K    (in) INTEGER 
!      The size of secular equation, the dimension of that broken arrow matrix.
!
! Ni   (in) INTEGER
!      The block size of most nodes at the bottom level
! 
! D    (in) DOUBLE PRECISION array, DIMENSION( K )
!      The singular values of that broken arrow matrix are ordered increasingly. 
!
! Z    (in) DOUBLE PRECISION array, DIMENSION( K )
!      The appended row of that broken arrow matrix
!
! F    (in) DOUBLE PRECISION array, DIMENSION( K )
!      The diagonal elements of that broken arrow matrix, except 
!      F(1) = zero. F(2-K) are the old singular values of previous level.
!
! DIFL (in) DOUBLE PRECISION array, DIMENSION( K )
!      The difference between old singular values and new singular values
!      DIFL(J) = D(J) - F(J), positive value, the same as that in LAPACK. 
!
! DIFR (in) DOUBLE PRECISION array, DIMENSION( K )
!      The difference between old singular values and new singular values
!      DIFR(J) = D(J) - F(J+1), negative value, the same as LAPACK. DIFR(K)
!      will not be referenced. 
!
! APHAU (inout) DOUBLE PRECISION array, DIMENSION( K )
!      The normalization scalars of the right singular vectors of that 
!      broken arrow matrix, the reciprocal of the 2-norm of columns of UH. 
! 
! APHAV (inout) DOUBLE PRECISION array, DIMENSION( K )
!      The normalization scalars of the left singular vectors of that 
!      broken arrow matrix, the reciprocal of the 2-norm of rows of VHT. 
!
! ZZ   (inout) DOUBLE PRECISION array, DIMENSION( K )
!      One copy of Z
! 
! UT    (inout) DOUBLE PRECISION array, DIMENSION( K, LCUT )
!      The transpose of first K columns of left singular vectors that are needed
!      to be updated. 
!
! LCUT  (in) INTEGER
!      The row dimension of left singular vector of this level
! 
! VT   (inout) DOUBLE PRECISION array, DIMENSION( K, LCVT )
!      The first K rows of right singular vectors that are needed
!      to be updated. 
! 
! LCVT (in) INTEGER
!      The column (row) dimension of right singular vector of this level
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
! Written by S.-G. Li, on Dec. 7th, 2012, in Changsha China.
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
    DOUBLE PRECISION  timerrlu, Maxerr
! ..
! ..Parameters for HSS tree..
    INTEGER              :: ltr,lm,n1,lvl,pnh,IDD,IH,IWW
    type(HTinfo)         :: TRE
    type(hssmm)          :: PH
    INTEGER, ALLOCATABLE :: tr(:), M(:)
    DOUBLE PRECISION, ALLOCATABLE :: DD(:),HH(:)

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
       WRITE(*,*) 'Allocate failed in dhssvcs'
    END IF

    call npart( K,Ni,tr,ltr,M,lm )

!   Tree information
    lvl = ceiling( log2(ltr+1) ) -1
    n1 = 0.5* (ltr +1 )

    ALLOCATE( TRE%ttr(lvl,n1),TRE%lentr(lvl),TRE%lfx(ltr),TRE%ch(2,ltr), stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssvcs'
    END IF

    call GetTreeInfo( TR,LTR,M,M,TRE,lvl )

    ALLOCATE( PH%D(ltr),PH%U(2,ltr),PH%V(2,ltr),PH%B(2,ltr),PH%ph(4,ltr),PH%pd(2,ltr),stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssvcs'
    END IF

    call hssexpmm0( PH, ltr )

    tnk = 0
    DO i = 1, lm
       tnk = tnk+ M(i)**2
    END DO
    IDD = 1
    IH  = IDD + tnk
    IWW = IH + 10 * Ni *K
    ALLOCATE( HH( 10*Ni*K ), DD( tnk ), stat=ierr ) 
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in DHSSVCS'
    END IF

! ****************************************************************
!                  HSS Matrix Approximation for Right            *
! ****************************************************************
    pnh = 1
    nswap = 0
    nflops = 0
    ZZ(1:K) = Z(1:K)

!!$    call cauchy2hssvd( D,F,ZZ,APHAV,DIFL,DIFR,K,TR,LTR,M,LM,PH,WORK(IH),WORK(IDD),TOL,lvl,pnh,timerrlu,nflops,nswap )
    call cauchy2hssvd( D,F,APHAV,ZZ,DIFL,DIFR,K,TR,LTR,M,LM,PH,HH,DD,TOL,lvl,pnh,timerrlu,nflops,nswap ) 
    
    NC = FLOOR( DBLE(LCVT) /NCOL )
    DO i = 1, NC
       call fastHssmm_omp( DD,HH,PH,M,VT(1:K,(i-1)*NCOL+1:i*NCOL),tr,ltr,K,NCOL,TRE,lvl )
    END DO
    nlast = LCVT - NC * NCOL
    IF( nlast .NE. 0 ) THEN
       call fastHssmm_omp( DD,HH,PH,m,VT(1:K,NC*NCOL+1:LCVT ),tr,ltr,K,nlast,TRE,lvl )
    END IF

! ****************************************************************
!                  HSS Matrix Approximation for Left             *
! **************************************************************** 
!   We use HSS matrix to approximate hat{UT}, the first column of which is
!   zero and then which is a Cauchy matrix. 

!   Define new ZZ. ZZ(1) equal to zero. 
    ZZ(1:K) = F(1:K) * Z(1:K)
    Z(1:K)  = APHAU(1:K) ! To protect APHAU

!   HSS approximation to the left singular vector matrix 
    call hssexpmm0( PH, ltr )
    pnh = 1
    call cauchy2hssvdU( D,F,ZZ,Z,DIFL,DIFR,K,TR,LTR,M,LM,PH,WORK(IH),WORK(IDD),TOL,lvl,pnh,timerrlu,nflops,nswap )
    
!   Multiplication of hat(UT) and UT
    NC = FLOOR( DBLE(LCUT) /NCOL )
    Do i = 1, NC
       call fastHssmmL_omp( WORK(IDD),WORK(IH),PH,M,UT((i-1)*NCOL+1:i*NCOL,1:K),tr,ltr,NCOL,K,TRE,lvl )
    End Do ! ( NC )
    nlast = LCUT - NC * NCOL
    IF( nlast .NE. 0 ) THEN
       i = nc + 1
       call fastHssmmL_omp( WORK(IDD),WORK(IH),PH,M,UT((i-1)*NCOL+1:LCUT,1:K ),tr,ltr,nlast,K,TRE,lvl )
    END IF ! ( nlast )

    DEALLOCATE( TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch )
    DEALLOCATE( PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd )
    DEALLOCATE( M, tr ) 
  END SUBROUTINE DHSSVCS

!!!!!!
  SUBROUTINE DHSSVCS_BD( K, SB, Ni, D, Z, F, DIFL, DIFR, APHAU, APHAV, ZZ, &
                      UT,LCU,VT,LCVT,WORK,INFO )
     USE aux_hss
     USE BasicMM
!
! .. Scalar Parameters ..
    INTEGER, INTENT(IN) ::  K, LCU, LCVT, INFO, Ni, SB
!
! .. Array Parameters ..
    DOUBLE PRECISION, INTENT(IN) :: D(*),F(*),DIFL(*),DIFR(*)
    DOUBLE PRECISION, INTENT(INOUT) :: APHAU(*), APHAV(*),WORK(*),Z(LCVT,*), & 
                                       UT(LCU,K),VT(K,LCVT),ZZ(*)
!
! .. Parameters ..
    INTEGER           NCOL
    PARAMETER         ( NCOL = 200 )
    DOUBLE PRECISION  ONE, ZERO, TOL
    PARAMETER         ( ONE= 1.0D+0, ZERO= 0.0E+0, TOL= 1.0E-17 )
! ..
! .. Local Scalars ..
    INTEGER  i, nflops, nswap, nc, nlast,ierr,tnk, mr
    DOUBLE PRECISION  Maxerr, timerrlu
! ..
! ..Parameters for HSS tree..
    INTEGER              :: ltr,lm,n1,lvl,pnh,IDD,IH,IWW
    type(HTinfo)         :: TRE
    type(hssmm)          :: PH
    INTEGER, ALLOCATABLE :: tr(:), M(:)
    DOUBLE PRECISION, ALLOCATABLE :: TA(:,:),DD(:),HH(:),HA(:,:),ZT(:)
! SX(:,:), TA(:,:),DD(:),HH(:),HA(:,:)

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
    ALLOCATE ( M(lm),tr( 2*lm-1 ),ZT(K), stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssvcs'
    END IF

    call npart( K,Ni,tr,ltr,M,lm )

!   Tree information
    lvl = ceiling( log2(ltr+1) ) -1
    n1 = 0.5* (ltr +1 )

    ALLOCATE( TRE%ttr(lvl,n1),TRE%lentr(lvl),TRE%lfx(ltr),TRE%ch(2,ltr), stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssvcs'
    END IF

    call GetTreeInfo( TR,LTR,M,M,TRE,lvl )

    ALLOCATE( PH%D(ltr),PH%U(2,ltr),PH%V(2,ltr),PH%B(2,ltr),PH%ph(4,ltr),PH%pd(2,ltr),stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssvcs'
    END IF

    call hssexpmm0( PH, ltr )

    tnk = 0
    DO i = 1, lm
       tnk = tnk+ M(i)**2
    END DO
    IDD = 1
    IH  = IDD + tnk
    IWW = IH + 10 * Ni *K
    ALLOCATE( HH( 10*Ni*K ),DD(tnk), stat=ierr ) 
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in DHSSVCS'
    END IF

! ****************************************************************
!                  HSS Matrix Approximation for Right            *
! ****************************************************************
    pnh = 1
    nswap = 0
    nflops = 0
    ZZ(1:K) = Z(SB:SB+K-1,SB)

!    call cauchy2hssvd_omp( D,F,APHAV,ZZ,DIFL,DIFR,K,TR,LTR,M,LM,PH,HH,DD,TOL,TRE,lvl,pnh,nflops,nswap ) 
    call cauchy2hssvd( D,F,APHAV,ZZ,DIFL,DIFR,K,TR,LTR,M,LM,PH,HH,DD,TOL,lvl,pnh,timerrlu,nflops,nswap ) 
    
    NC = FLOOR( DBLE(LCVT) /NCOL )
    DO i = 1, NC
       call fastHssmm_omp( DD,HH,PH,M,VT(1:K,(i-1)*NCOL+1:i*NCOL),tr,ltr,K,NCOL,TRE,lvl )
    END DO
    nlast = LCVT - NC * NCOL
    IF( nlast .NE. 0 ) THEN
       call fastHssmm_omp( DD,HH,PH,m,VT(1:K,NC*NCOL+1:LCVT ),tr,ltr,K,nlast,TRE,lvl )
    END IF
!
!   Update the other columns by multiplying them by VT
    IF( SB .GT. 1 ) THEN
       call fastHssmm_omp( DD,HH,PH,M,Z(SB:SB+K-1,1:SB-1 ),tr,ltr,K,SB-1,TRE,lvl )
    END IF

! ****************************************************************
!                  HSS Matrix Approximation for Left             *
! **************************************************************** 
!   We use HSS matrix to approximate hat{UT}, the first column of which is
!   zero and then which is a Cauchy matrix. 

!   Define new ZZ. ZZ(1) equal to zero. 
    ZT(1:K) = Z( SB:SB+K-1,SB )
    ZZ(1:K) = F(1:K) * ZT(1:K)
    Z( SB:SB+K-1,SB )  = APHAU( 1:K ) ! To protect APHAU

!   HSS approximation to the left singular vector matrix 
    call hssexpmm0( PH, ltr )
    pnh = 1
!    call cauchy2hssvdU_omp( F,D,ZZ,ZT,DIFL,DIFR,K,TR,LTR,M,LM,PH,HH,DD,TOL,TRE,lvl,pnh,nflops,nswap )
    call cauchy2hssvdU( D,F,ZZ,ZT,DIFL,DIFR,K,TR,LTR,M,LM,PH,HH,DD,TOL,lvl,pnh,timerrlu,nflops,nswap )

!    write(*,*) 'Left construction is right'

    ALLOCATE( HA(K,K), stat=ierr )
    IF( ierr .ne. 0) THEN
       WRITE(*,*) 'Allocate failed in dhssvcs for fasthssmm'
    END IF

!   Multiplication of hat(UT) and UT
    NC = FLOOR( DBLE(LCU) /NCOL )
    Do i = 1, NC
       call fastHssmmL_omp( DD,HH,PH,M,UT((i-1)*NCOL+1:i*NCOL,1:K),tr,ltr,NCOL,K,TRE,lvl )
    End Do ! ( NC )
    nlast = LCU - NC * NCOL
    IF( nlast .NE. 0 ) THEN
       i = nc + 1
       call fastHssmmL_omp( DD, HH,PH,M,UT((i-1)*NCOL+1:LCU,1:K ),tr,ltr,nlast,K,TRE,lvl )
    END IF ! ( nlast )

    DEALLOCATE( TRE%ttr, TRE%lentr, TRE%lfx, TRE%ch )
    DEALLOCATE( PH%D, PH%U, PH%V, PH%B, PH%ph, PH%pd )
    DEALLOCATE( M, tr, ZT ) 

  END SUBROUTINE DHSSVCS_BD

end module CauchyHssvd_VPS

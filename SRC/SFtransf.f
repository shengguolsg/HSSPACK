      SUBROUTINE DSY2BUC( UPLO, N, B, AFULL, LDFULL, ABAND, LDBAND,
     $                    INFO )
*
      IMPLICIT NONE
* ..
* .. Scalars arguments..
        CHARACTER*1           UPLO
        INTEGER               N, B, LDFULL, LDBAND, INFO
* ..
* .. Array arguments ..
        DOUBLE PRECISION      AFULL( LDFULL,* ), ABAND( LDBAND,* )
*
* Purpose
* =======
*
*   The full routine contains sparse-full transformation routines.
*
*   dsy2buc copies a full symmetric banded A matrix from (upper or lower) (full) storage 
*   to upper banded storage.
*
*   Note that afull and aband must refer to different memory locations,
*   i.e., A may NOT be repacked within the same array.
*
* Arguments
* =========
*
*   UPLO    (in) character*1
*           Is the matrix A stored in the upper or lower triangle of the
*           array afull ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   N       (in) integer
*           The size of the matrix A.
*           n >= 0.
*
*   b       (in) integer
*           The (semi-)bandwidth of the matrix A.
*           0 <= b < n, if n >= 1.
*           b = 0     , if n = 0.
*
*   AFULL   (in) double precision array, dimension ( LDFULL,N )
*           The banded matrix A in upper or lower storage.
*
*   LDFULL  (in) integer
*           The leading dimension of the array afull.
*           LDFULL >= MAX( N, 1 ).
*
*   ABAND   (out) double precision array, dimension ( LDBAND, N )
*           The banded matrix A in upper banded storage (lower
*           banded storage is not supported).
*
*   LDBAND  (in)
*           The leading dimension of the array aband.
*           LDBAND >= b + 1.
*
*   INFO    (out) integer
*           On exit, info indicates the consistency of the arguments.
*           info = -1 : uplo was none of 'U', 'L' (upper/lower case).
*                = -2 : n was out of range (negative).
*                = -3 : b was out of range (negative or >= n ).
*                = -5 : ldfull was out of range ( < n or < 1 ).
*                = -7 : ldband was out of range ( < b + 1 ).
*           info =  1 : All arguments were ok.
*
* ===================
* Author: Shengguo Li, Nati. Univ. of Def. Tech., China
* email:  nudtlsg@gmail.com
* Date: May 20, 2013, Version: HSSPACK, V-1.0.0
* 
* Tested on May 20th, 2013
* ======================================================
*
* Local variables:
*
        LOGICAL               UPPER
        INTEGER               J, I
*
*   upper   is A given in upper symmetric storage ?
*
* Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*     --- check for errors in the parameters ---
*
        UPPER = LSAME( UPLO, 'U' )
*
        IF ( .NOT. ( UPPER .OR. ( LSAME( UPLO, 'L' ) ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( B .LT. 0 ) .OR.
     $           ( ( N .GT. 0 ) .AND. ( B .GE. N ) ) ) THEN
          INFO = - 3
        ELSEIF ( LDFULL .LT. MAX( N, 1 ) ) THEN
          INFO = - 5
        ELSEIF ( LDBAND .LT. ( B+1 ) ) THEN
          INFO = - 7
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 990
*
*       --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 990
*
*       --- non-trivial case ---
*
        IF ( UPPER ) THEN
*
*         --- "upper to upper" ---
          DO 110 J = 1, N
            DO 100 I = MAX( J-B, 1 ), J
              ABAND( B+1+I-J, J ) = AFULL( I, J )
  100       CONTINUE
  110     CONTINUE
        ELSE
*
*         --- "lower to upper" ---
          DO 210 J = 1, N
            DO 200 I = J, MIN( J+B, N )
              ABAND( B+1+J-I, I ) = AFULL( I, J )
  200       CONTINUE
  210     CONTINUE
        ENDIF
*
 990    RETURN
        END
*
      SUBROUTINE DSB2FLC( UPLO, N, M, B, ABAND, LDBAND, AFULL, LDFULL, 
     $                    INFO )
*
      IMPLICIT NONE
* ..
* .. Scalars arguments..
        CHARACTER*1           UPLO
        INTEGER               N, B, LDFULL, LDBAND, INFO, M
* ..
* .. Array arguments ..
        DOUBLE PRECISION      AFULL( LDFULL,* ), ABAND( LDBAND,* )
*
* Purpose
* =======
*
*   DSB2FLC copies an upper or lower N-by-M banded A matrix of compact 
*   storage form to its corresponding full storage form. 
*
*   Note that afull and aband must refer to different memory locations,
*   i.e., A may NOT be repacked within the same array.
*
* Arguments
* =========
*
*   UPLO    (in) character*1
*           Is the matrix A stored in the upper or lower triangle of the
*           array afull ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   N       (in) integer
*           The row size of the matrix A.
*           n >= 0.
*
*   M       (in) integer
*           The column size of the matrix A.
*           m > = 0, which is only used in the upper banded case. 
*
*   b       (in) integer
*           The (semi-)bandwidth of the matrix A.
*           0 <= b < n, if n >= 1.
*           b = 0     , if n = 0.
*
*   AFULL   (in) double precision array, dimension ( LDFULL,N )
*           The banded matrix A in upper or lower storage.
*
*   LDFULL  (in) integer
*           The leading dimension of the array afull.
*           LDFULL >= MAX( N, 1 ).
*
*   ABAND   (out) double precision array, dimension ( LDBAND, N )
*           The banded matrix A in upper banded storage (lower
*           banded storage is not supported).
*
*   LDBAND  (in)
*           The leading dimension of the array aband.
*           LDBAND >= b + 1.
*
*   INFO    (out) integer
*           On exit, info indicates the consistency of the arguments.
*           info = -1 : uplo was none of 'U', 'L' (upper/lower case).
*                = -2 : n was out of range (negative).
*                = -3 : b was out of range (negative or >= n ).
*                = -5 : ldfull was out of range ( < n or < 1 ).
*                = -7 : ldband was out of range ( < b + 1 ).
*           info =  1 : All arguments were ok.
*
* ===================
* Author: Shengguo Li, Nati. Univ. of Def. Tech., China
* email:  nudtlsg@gmail.com
* Date: May 20, 2013, Version: HSSPACK, V-1.0.0
* 
* Tested on May 20th, 2013
* ======================================================
*
* Local variables:
*
        LOGICAL               UPPER
        INTEGER               J, I
        DOUBLE PRECISION      ZERO
        PARAMETER             ( ZERO = 0.0D0 )
*
*   upper   is A given in upper symmetric storage ?
*
*   Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*     --- check for errors in the parameters ---
*
        UPPER = LSAME( UPLO, 'U' )
*
        IF ( .NOT. ( UPPER .OR. ( LSAME( UPLO, 'L' ) ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( B .LT. 0 ) .OR.
     $           ( ( N .GT. 0 ) .AND. ( B .GT. N ) ) ) THEN
          INFO = - 3
        ELSEIF ( LDFULL .LT. MAX( N,1 ) ) THEN
          INFO = - 5
        ELSEIF ( LDBAND .LT. ( B+1 ) ) THEN
          INFO = - 7
        ELSE
          INFO = 0
        ENDIF
*
        IF ( INFO .NE. 0 )     GOTO 991
*
*     --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 991
*
*     --- non-trivial case ---
*
        DO 300 I = 1, N
           DO 400 J = 1, M
              AFULL( I,J ) = ZERO
  400      CONTINUE
  300   CONTINUE
*
        IF( UPPER ) THEN
*
*         --- "upper to upper" ---
          DO 510 J = 1, N
            DO 500 I = MAX( J-B, 1 ), J
               AFULL( I,J ) = ABAND( B+1+I-J,J )
  500       CONTINUE
  510     CONTINUE
          IF( M .GT. N ) THEN
             DO 900 J = N+1, M
                DO 910 I = J-B, N
                   AFULL( I,J ) = ABAND( I-J+B+1,J )
  910           CONTINUE
  900        CONTINUE
          END IF
        ELSE
*
*         --- "lower to lower, this works for the cases of M <= N <= M+b " ---
          DO 610 J = 1, N
            DO 600 I = J, MIN( J+B, N )  ! ?? J+B 
               AFULL( I, J ) = ABAND( 1+I-J, J )
  600       CONTINUE

  610     CONTINUE
        ENDIF
*
 991    RETURN
        END
*
      SUBROUTINE DFL2SBC( UPLO, N, B, AFULL, LDFULL, ABAND, LDBAND,
     $                     INFO )
*
      IMPLICIT NONE
* ..
* .. Scalars arguments..
        CHARACTER*1           UPLO
        INTEGER               N, B, LDFULL, LDBAND, INFO
* ..
* .. Array arguments ..
        DOUBLE PRECISION      AFULL( LDFULL, * ), ABAND( LDBAND, * )
*
* Purpose
* =======
*
*   DFL2SBC copies a full banded A matrix from (upper or lower) (full) storage 
*   to its corresponding sparse storage form. 
*
*   Note that afull and aband must refer to different memory locations,
*   i.e., A may NOT be repacked within the same array.
*
* Arguments
* =========
*
*   UPLO    (in) character*1
*           Is the matrix A stored in the upper or lower triangle of the
*           array afull ?
*           uplo = 'U' : Upper triangle.
*                = 'L' : Lower triangle.
*
*   N       (in) integer
*           The size of the matrix A.
*           n >= 0.
*
*   b       (in) integer
*           The (semi-)bandwidth of the matrix A.
*           0 <= b < n, if n >= 1.
*           b = 0     , if n = 0.
*
*   AFULL   (in) double precision array, dimension ( LDFULL,N )
*           The banded matrix A in upper or lower storage.
*
*   LDFULL  (in) integer
*           The leading dimension of the array afull.
*           LDFULL >= MAX( N, 1 ).
*
*   ABAND   (out) double precision array, dimension ( LDBAND, N )
*           The banded matrix A in upper banded storage (lower
*           banded storage is not supported).
*
*   LDBAND  (in)
*           The leading dimension of the array aband.
*           LDBAND >= b + 1.
*
*   INFO    (out) integer
*           On exit, info indicates the consistency of the arguments.
*           info = -1 : uplo was none of 'U', 'L' (upper/lower case).
*                = -2 : n was out of range (negative).
*                = -3 : b was out of range (negative or >= n ).
*                = -5 : ldfull was out of range ( < n or < 1 ).
*                = -7 : ldband was out of range ( < b + 1 ).
*           info =  1 : All arguments were ok.
*
* ===================
* Author: Shengguo Li, Nati. Univ. of Def. Tech., China
* email:  nudtlsg@gmail.com
* Date: May 20, 2013, Version: HSSPACK, V-1.0.0
* 
* Tested on May 20th, 2013
* ======================================================
*
* Local variables:
*
        LOGICAL               UPPER
        INTEGER               J, I
*
*   upper   is A given in upper symmetric storage ?
*
*   Routines called:
*
        LOGICAL               LSAME
        EXTERNAL              LSAME
*
*   lsame   case-insensitive character matching (BLAS)
*
        INTRINSIC             MAX, MIN
*
* ----------------------------------------------------------------------
*
*       --- check for errors in the parameters ---
*
        UPPER = LSAME( UPLO, 'U' )
*
        IF ( .NOT. ( UPPER .OR. ( LSAME( UPLO, 'L' ) ) ) ) THEN
          INFO = - 1
        ELSEIF ( N .LT. 0 ) THEN
          INFO = - 2
        ELSEIF ( ( B .LT. 0 ) .OR.
     $           ( ( N .GT. 0 ) .AND. ( B .GE. N ) ) ) THEN
          INFO = - 3
        ELSEIF ( LDFULL .LT. MAX( N, 1 ) ) THEN
          INFO = - 5
        ELSEIF ( LDBAND .LT. ( B+1 ) ) THEN
          INFO = - 7
        ELSE
          INFO = 1
        ENDIF
*
        IF ( INFO .NE. 1 )     GOTO 992
*
*       --- check for quick return ---
*
        IF ( N .EQ. 0 )     GOTO 992
*
*          --- non-trivial case ---
*
        IF ( UPPER ) THEN
*
*         --- "upper to upper, this works for the cases of N <= M <= N+b " 
          DO 710 J = 1, N
            DO 700 I = MAX( J-B, 1 ), J
              ABAND( B+1+I-J, J ) = AFULL( I, J )
  700       CONTINUE
  710     CONTINUE
        ELSE
*
*         --- "lower to lower" ---
          DO 810 J = 1, N
            DO 800 I = J, MIN( J+B, N )
              ABAND( 1+I-J, J ) = AFULL( I, J )
  800       CONTINUE
  810     CONTINUE
        ENDIF
*
 992    RETURN
        END

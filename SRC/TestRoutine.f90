module TestRoutine

  implicit none

CONTAINS

!!!!!!!
  subroutine TestSVD( N,M,SB,AB,LDAB,A,U,LDU,S,VT,LDVT,WORK,info)
!
! .. Scalar parameter ..
    INTEGER, INTENT(IN)  :: M,N,LDAB,LDU,LDVT,info,SB
!
! .. Array parameters ..
    DOUBLE PRECISION, INTENT(IN) :: S(*),AB(LDAB,*),U(LDU,*),VT(LDVT,*)
    DOUBLE PRECISION, INTENT(INOUT) :: A(*), WORK(*)
!
! Purpose
! ==========
! This routine first copies a N-by-M submatrix of AB to A, and then 
! construct matrix C = U*S*VT. Finally, check whether C-A is small.
!
! ==========
! Written by Shengguo Li, on July 14, 2013
! ========================================
! 
    INTEGER IC, IE, IU, I, ITEMP, J
    DOUBLE PRECISION  ERR
    DOUBLE PRECISION, PARAMETER :: ONE = 1.0E0, ZERO = 0.0E0

    IU = 1
    IC = IU + N*N
    IE = IC + M*N
    CALL DSB2FLC( 'U',N,M,SB,AB,LDAB,A,M,INFO )
    
    CALL DLACPY( 'A', N,N,U,LDU,WORK(IU),N )
    DO I = 1, N
       ITEMP = (I-1)*N+1
       CALL DSCAL( N,S(I),WORK(ITEMP),1 )
    END DO
    CALL DGEMM('N','N',N,M,N,ONE,WORK(IU),N,VT,LDVT,ZERO,WORK(IC),N)
    
    DO I = 1, M
       DO J = 1, N
          WORK(IE+J-1+(I-1)*N) = WORK(IC+J-1+(I-1)*N) - A( J+(I-1)*M )
       END DO
    END DO
    
    ERR = maxval( abs( WORK(IE:IE+M*N-1) ) )
    write(*,*) 'Max error is ', Err

  end subroutine TestSVD

!!!!!!!
  subroutine TestSVD1( N,M,SB,AB,LDAB,A,U,LDU,S,VT,LDVT,ZZ,WORK,ORGNRM,info )
!
! .. Scalar parameter ..
    INTEGER, INTENT(IN)  :: M,N,LDAB,LDU,LDVT,info,SB
    DOUBLE PRECISION, INTENT(IN) :: ORGNRM
!
! .. Array parameters ..
    DOUBLE PRECISION, INTENT(IN) :: AB(LDAB,*),U(LDU,*),VT(LDVT,*),ZZ(M,*)
    DOUBLE PRECISION, INTENT(INOUT) :: A(*),WORK(*),S(*)
!
! Purpose
! ==========
! This routine first copies a N-by-M submatrix of AB to A, and then 
! construct matrix C = U*S*VT. Finally, check whether C-A is small.
!
! ==========
! Written by Shengguo Li, on July 14, 2013
! ========================================
! 
    INTEGER IVT, IE, IU, I, ITEMP, J, SBM1, IZ
    DOUBLE PRECISION  ERR
    DOUBLE PRECISION, PARAMETER :: ONE = 1.0E0, ZERO = 0.0E0

    IZ = 1 
    IVT = IZ + M*(SB-1)
    IU = IVT + M*N
    IE = IU + N*N
    SBM1 = SB - 1
    A(1:M*M) = ZERO
    CALL DSB2FLC( 'U',N,M,SB,AB,LDAB,A,M,INFO )

    call transp2d( ZZ,M,SB-1,WORK(IZ) )
    CALL DLACPY( 'A',N,M,VT,LDVT,WORK(IVT),N )

    CALL DLACPY( 'A',N,N,U,LDU,WORK(IU),N )
    DO I = SB, N
       ITEMP = (I-1)*N+IU
       CALL DSCAL( N,S(I),WORK(ITEMP),1 )
    END DO
    call DGEMM( 'N','N',SBM1,M,N,ONE,WORK(IZ),SBM1,VT,LDVT,ZERO,WORK(IVT),N )
    CALL DGEMM( 'N','N',N,M,N,ONE,WORK(IU),N,WORK(IVT),N,ZERO,WORK(IE),N )
    
    DO I = 1, M
       DO J = 1, N
          WORK(IVT+J-1+(I-1)*N) = WORK(IE+J-1+(I-1)*N)*ORGNRM - A( J+(I-1)*M )
       END DO
    END DO

!!$    open( unit=5, file='Work' )
!!$    write(5,*) WORK(IVT:IVT+M*N-1)
!!$    close(5)
    
    ERR = maxval( abs( WORK(IVT:IVT+M*N-1) ) )
    write(*,*) 'Max error of level 2 is ', Err

  end subroutine TestSVD1

!!!!!!!
  subroutine TestSVD3( N,M,SB,AB,LDAB,A,U,LDU,S,VT,LDVT,ZZ,WORK,ORGNRM,ZB,info )
!
! .. Scalar parameter ..
    INTEGER, INTENT(IN)  :: M,N,LDAB,LDU,LDVT,info,SB,ZB
    DOUBLE PRECISION, INTENT(IN) :: ORGNRM
!
! .. Array parameters ..
    DOUBLE PRECISION, INTENT(IN) :: AB(LDAB,*),U(LDU,*),VT(LDVT,*),ZZ(M,*)
    DOUBLE PRECISION, INTENT(INOUT) :: A(*),WORK(*),S(*)
!
! Purpose
! ==========
! This routine first copies a N-by-M submatrix of AB to A, and then 
! construct matrix C = U*S*VT. Finally, check whether C-A is small.
!
! ==========
! Written by Shengguo Li, on July 14, 2013
! ========================================
! 
    INTEGER IVT, IE, IU, I, ITEMP, J, SBM1, IZ
    DOUBLE PRECISION  ERR
    DOUBLE PRECISION, PARAMETER :: ONE = 1.0E0, ZERO = 0.0E0

    IZ = 1 
    IVT = IZ + M*(ZB-1)
    IU = IVT + M*N
    IE = IU + N*N
    SBM1 = ZB - 1
    A(1:M*M) = ZERO
    CALL DSB2FLC( 'U',N,M,SB,AB,LDAB,A,M,INFO )

    call transp2d( ZZ,M,ZB-1,WORK(IZ) )
    CALL DLACPY( 'A',N,M,VT,LDVT,WORK(IVT),N )

    CALL DLACPY( 'A',N,N,U,LDU,WORK(IU),N )
    DO I = ZB, N
       ITEMP = (I-1)*N+IU
       CALL DSCAL( N,S(I),WORK(ITEMP),1 )
    END DO
    IF(SBM1 .GE. 1 ) THEN
       call DGEMM( 'N','N',SBM1,M,N,ONE,WORK(IZ),SBM1,VT,LDVT,ZERO,WORK(IVT),N )
    END IF
    CALL DGEMM( 'N','N',N,M,N,ONE,WORK(IU),N,WORK(IVT),N,ZERO,WORK(IE),N )
    
    DO I = 1, M
       DO J = 1, N
          WORK(IVT+J-1+(I-1)*N) = WORK(IE+J-1+(I-1)*N)*ORGNRM - A( J+(I-1)*M )
       END DO
    END DO

    open( unit=5, file='Work' )
    write(5,*) WORK(IVT:IVT+M*N-1)
    close(5)

    ERR = maxval( abs( WORK(IVT:IVT+M*N-1) ) )
    write(*,*) 'Max error of level 2 is ', Err

  end subroutine TestSVD3

!!!!!!!
  subroutine TestSVD2( N,M,SB,AB,LDAB,A,U,LDU,S,VT,LDVT,ZZ,WORK,info)
!
! .. Scalar parameter ..
    INTEGER, INTENT(IN)  :: M,N,LDAB,LDU,LDVT,info,SB
!
! .. Array parameters ..
    DOUBLE PRECISION, INTENT(IN) :: AB(LDAB,*),U(LDU,*),VT(LDVT,*),ZZ(M,*)
    DOUBLE PRECISION, INTENT(INOUT) :: A(*), WORK(*),S(*)
!
! Purpose
! ==========
! This routine first copies a N-by-M submatrix of AB to A, and then 
! construct matrix C = U*S*VT. Finally, check whether C-A is small.
!
! ==========
! Written by Shengguo Li, on July 14, 2013
! ========================================
! 
    INTEGER IVT, IE, IU, I, ITEMP, J, SB1, IZ
    DOUBLE PRECISION  ERR
    DOUBLE PRECISION, PARAMETER :: ONE = 1.0E0, ZERO = 0.0E0

    IZ = 1 
    IVT = IZ + M*SB
    IU = IVT + M*N
    IE = IU + N*N
    SB1 = SB + 1
    A(1:M*M) = ZERO
    CALL DSB2FLC( 'U',N,M,SB,AB,LDAB,A,M,INFO )

    call transp2d( ZZ,M,SB,WORK(IZ) )
    CALL DLACPY( 'A',N,M,VT,LDVT,WORK(IVT),N )

    CALL DLACPY( 'A',N,N,U,LDU,WORK(IU),N )
    DO I = SB1, N
       ITEMP = (I-1)*N+IU
       CALL DSCAL( N,S(I),WORK(ITEMP),1 )
    END DO
    call DGEMM( 'N','N',SB,M,N,ONE,WORK(IZ),SB,VT,LDVT,ZERO,WORK(IVT),N )
    CALL DGEMM( 'N','N',N,M,N,ONE,WORK(IU),N,WORK(IVT),N,ZERO,WORK(IE),N )
    
    DO I = 1, M
       DO J = 1, N
          WORK(IVT+J-1+(I-1)*N) = WORK(IE+J-1+(I-1)*N) - A( J+(I-1)*M )
       END DO
    END DO
    
    ERR = maxval( abs( WORK(IVT:IVT+M*N-1) ) )
    write(*,*) 'Max error of level 2 is ', Err

  end subroutine TestSVD2

end module TestRoutine

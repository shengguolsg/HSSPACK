module ConstructHssd
! This procedure is random svd + column majority. 
  implicit none

contains

!!!!!!!!!
  subroutine mat2hssvdr(A, LDA, TR, LTR, M, LM, PH, H, DD, TOL, RK, T, PP, lvl,pnh)
       use aux_hss
!
!  .. Scalar Arguments ..
       INTEGER          :: LDA, PP, RK, LTR, lm, lvl
       DOUBLE PRECISION :: TOL
!
!  .. Array Arguments ..
       DOUBLE PRECISION :: A(LDA, *), H(*), DD(*), T(*)
       INTEGER          :: TR(*), M(*)
       TYPE(HSSMM)      :: PH
! 
!  Purpose
! =========
! Construct an HSS matrix for a rank-structured matrix A
! 
! Parameters
!========== 
!  A    (INPUT) DOUBLE PRECISION array, DIMENSION(lda, lda)
!       We assume A is given explicitly and a modified version without forming A will
!       be written later. A is used as partly workspace and its entries will be changed.
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
!  RK   (INPUT) INTEGER, presetted rank for fixed-rank approximation
! 
!  T    (INPUT) DOUBLE PRECISION, DIMENSION( log2(lm)*mi* LDA ),
!       mi is the average of diagonal block size; 
!  
!  PP  (INPUT) INTEGER, a parameter approximating HSS-rank + 10, used random algorithm
!
!===================
!  Written by Shengguo Li, on Aug. 8th, 2012
!
!  Finished on Aug. 14th, 2012
!  
!===================

!  .. Local Scalars ..
       INTEGER       :: CH1, IERR, INFO, LWORK, pmr,pmc, rkc, rkr,i,j
       INTEGER       :: lt,it,ns,pnd,ni,pnh,mm,ir1,ir2,nt1,nt2,tpn,mi,nn, &
                        ppa
!  .. Local Arrays ..
       INTEGER, ALLOCATABLE :: CH(:,:), L(:,:)
       INTEGER, ALLOCATABLE :: lsr(:), wsr(:), lsc(:), wsc(:)
       DOUBLE PRECISION, ALLOCATABLE :: WW(:), WORK(:), svs(:), Omega(:)
       ! WW is an extra workspace for random sampling. WORK is an extra workspace for DGESVD
!  ..
!  .. External Subroutines ..
      EXTERNAL  DGESVD, DGEMM, DLACPY
!  ..
!  .. External Functions .. 
      
!  .. Intrinsic Functions ..
      INTRINSIC    MIN, ABS,MAXVAL
              
       MM = MAXVAL( M(1:lm) )
       LWORK = LDA*MM*2 ! modify later
       pnh = 1              ! 'pointer' in H
       pnd = 1              ! 'pointer' in D
       ppa = 10             ! number of oversampling
       allocate( CH(2,LTR),L(LTR,2),lsr(LTR), wsr(LTR), WW(LDA*MM*lvl), WORK(LDA*MM*2), stat=ierr )
       IF(ierr /= 0 ) Then
          write(*,*) "Allocate failed in math2hssvdr! "
          stop
       END IF
       allocate( svs(LDA),lsc(ltr),wsc(ltr), Omega((lvl*MM+10)*LDA), STAT=iERR )  ! Omega may not have enough entries.
       IF(ierr /= 0 ) Then
          write(*,*) "Allocate failed in math2hssvdr! "
          stop
       END IF

       call random_number(Omega)

       nn = ltr - 1 
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
       wsr(ns) = l(1,2)+1      ! column width at upper triangular
       lsc(ns) = 1             ! column start position at lower triangular
       wsc(ns) = l(1,2)+1      ! row width at lower triangular
       
       do i = 1, nn
          ! factorize (leaf) diag & get off-diag block

          if (ch(1,i) == 0 .and. ch(2,i) == 0) then 
             mi = l(i,2)-l(i,1)+1
             ni = lda-l(i,2)      !t^r, the right column block's length

             call dlacpy('A', mi, mi, A(l(i,1),l(i,1)), lda, DD(pnd), mi)    ! copy D{i}
             PH%D(i) = mi
             PH%pd(1,i) = pnd
             PH%pd(2,i) = pnd + mi*mi -1
             pnd = pnd + mi*mi

             if (ns == 1) then
                nt1 = ni    ! lenght of whole HSS block column
                nt2 = ni    ! lenght of whole HSS block row
                call dlacpy('A', ni, mi, A(l(i,2)+1,l(i,1)), lda, T(1+mi*ni), ni) ! Construct T, the HSS block column
                ! Use Random Sampling when nt is big
                IF( (ni .gt. 50) .and. (ni .gt. mi+ppa) ) THEN
                   pp = mi + ppa
                   call dgemm('N','N',pp,mi,ni,1d0,Omega,pp,T(1+mi*ni),ni,0d0,ww,pp) !ww=Omega*T, starting from 1.
                ELSE
                   pp = ni
                   call dlacpy('A',pp,mi,T(1+mi*ni),pp,ww,pp)   ! copy T2 to ww
                END IF

                call dgesvd('N','A',pp,mi,ww,pp,svs,ww,pp,ww(1+mi*pp),mi,work,lwork,info)  !svd of ww. 
                rkc = count(svs(1:mi)>tol)
                ir2 = mi - rkc   ! the number of compressed columns
                call dgemm('N','T',ni,rkc,mi,1d0,T(1+mi*ni),ni,ww(1+mi*pp),mi,0d0,T(1),ni) ! T2 is stored in T(1)
                call hssexpmm('V',PH,pnh,H,rkc,mi,i,ww(1+mi*pp),mi,info) ! Vi^T is in ww(1+mi*pp), the first rkc rows
!
                call dlacpy('A', mi, ni, A(l(i,1),l(i,2)+1), lda, T(1+(rkc+mi)*ni), mi) ! T1, HSS block row
                ! Use Random Sampling when nt is big
                IF( (ni .gt. 50) .and. (ni .gt. mi+ppa) ) THEN
                   pp = mi + ppa
                   call dgemm('N','N',mi,pp,ni,1d0,T(1+(rkc+mi)*ni),mi,Omega,ni,0d0,ww,mi)  ! ww = T*omega
                ELSE
                   pp = ni
                   call dlacpy('A',mi,ni,T(1+(rkc+mi)*ni),mi,ww,mi)   ! copy T1 to ww
                END IF

                call dgesvd('A','N',mi,pp,ww,mi,svs,ww(1+mi*pp),mi,ww,mi,work,lwork,info) 
                rkr = count(svs(1:mi)>tol)
                ir1 = mi - rkr   ! the number of compressed rows
                call dgemm('T','N',rkr,ni,mi,1d0,ww(1+mi*pp),mi,T(1+(rkc+mi)*ni),mi,0d0,T(1+rkc*ni),rkr) 
                ! T1 is stored in T(1+rkc*nt)
                call hssexpmm('U',PH,pnh,H,mi,rkr,i,ww(1+mi*pp),mi,info) ! U is in ww(1+mi*pp), the first rkr columns
!
             else   ! leaf nodes except node 1
!
                nt2 = lsr(ns)-1+ni        ! lenght of whole HSS block column
                call dlacpy('A', lsr(ns)-1, mi, A(1,wsr(ns-1)), lda, T(1+mi*nt2), nt2)
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A',ni,mi,A(l(i,2)+1,l(i,1)),lda,T(1+mi*nt2+lsr(ns)-1),nt2)     ! T2
                END IF
                ! Use Random Sampling when nt is big
                IF( (nt2 .gt. 50) .and. (nt2 .gt. mi+ppa) ) THEN
                   pp = mi + ppa
                   call dgemm('N','N',pp,mi,nt2,1d0,Omega,pp,T(1+mi*nt2),nt2,0d0,ww,pp)     ! ww=Omega*T, starting from 1.
                 ELSE
                   pp = nt2
                   call dlacpy('A',pp,mi,T(1+mi*nt2),pp,ww,pp)   ! copy T2 to ww
                END IF

                call dgesvd('N','A',pp,mi,ww,pp,svs,ww,pp,ww(1+mi*pp),mi,work,lwork,info)   ! svd of ww. 
                rkc = count(svs(1:mi)>tol)
                ir2 = mi - rkc   ! the number of compressed columns
                call dgemm('N','T',nt2,rkc,mi,1d0,T(1+mi*nt2),nt2,ww(1+mi*pp),mi,0d0,T(1),nt2) ! T2 is stored in T(1)
                call hssexpmm('V',PH,pnh,H,rkc,mi,i,ww(1+mi*pp),mi,info) ! Vi^T is in ww(1+mi*pp), the first rkc rows
!
                nt1 = lsc(ns)-1+ni        ! length of whole HSS block row
                tpn = 1+rkc*nt2+mi*nt1
                call dlacpy('A', mi,lsc(ns)-1, A(wsc(ns-1),1),lda,T(tpn),mi)
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A', mi,ni,A(l(i,1),l(i,2)+1), lda,T( tpn+mi*(lsc(ns)-1) ), mi) ! T1
                END IF
                ! Use Random Sampling when nt is big
                IF( (nt1 .gt. 50) .and. (nt1 .gt. mi+ppa) ) THEN
                   pp = mi + ppa
                   call dgemm('N','N',mi,pp,nt1,1d0,T(tpn),mi,Omega,nt1,0.0d0,ww,mi)  ! ww = T*omega                   
                ELSE
                   pp = nt1
                   call dlacpy('A',mi,nt1,T(tpn),mi,ww,mi)   ! copy T1 to ww
                END IF
                
                call dgesvd('A','N',mi,pp,ww,mi,svs,ww(1+mi*pp),mi,ww,mi,work,lwork,info) 
                rkr = count(svs(1:mi)>tol)
                ir1 = mi - rkr   ! the number of compressed rows
                call dgemm('T','N',rkr,nt1,mi,1d0,ww(1+mi*pp),mi,T(tpn),mi,0d0,T(1+rkc*nt2),rkr) 
                ! T1 is stored in T(1+rkc*nt)
                call hssexpmm('U',PH,pnh,H,mi,rkr,i,ww(1+mi*pp),mi,info) ! U is in ww(1+mi*pp), the first rkr columns

             endif
             
               lsr(ns+1) = lsr(ns)+rkr
               wsr(ns)   = l(i,2)+1
               lsc(ns+1) = lsc(ns)+rkc
               wsc(ns)   = l(i,2)+1
               ns = ns +1
               IF( l(i,2)+1 .lt. lda ) THEN
                  call dlacpy('A', ni, rkc, T(lsr(ns-1)), nt2, A(l(i,2)+1,lsc(ns-1)), lda)    ! T2
                  call dlacpy('A', rkr, ni, T(1+nt2*rkc+rkr*(lsc(ns-1)-1) ), rkr,A(lsr(ns-1),l(i,2)+1),lda) ! T1 ??
               END IF
!               
          else  ! parent nodes
!
             pmr = lsr(ns)-lsr(ns-2)  ! row dim
             pmc = lsc(ns)-lsc(ns-2)  ! col dim
             ni = lda-l(i,2)

!            store B 
             ch1 = ch(1,i)
             lt  = lsr(ns-1)-lsr(ns-2)
             tpn = wsr(ns-1)-wsr(ns-2)
             call hssexpmm('B',PH,pnh,H,lt,tpn,ch1,A(lsr(ns-2),wsr(ns-2)), lda, info)  ! B{ch{tr{i}}(1)}
             ch1 = ch(2,i)
             lt  = wsc(ns-1)-wsc(ns-2)
             tpn = lsc(ns-1)-lsc(ns-2)
             call hssexpmm('B',PH,pnh,H,lt,tpn,ch1,A(wsc(ns-2),lsc(ns-2)), lda, info)  ! B{ch{tr{i}}(2)}

             if (ns == 3) then
                nt1 = ni    ! lenght of whole HSS block column
                nt2 = ni    ! lenght of whole HSS block row
                call dlacpy('A', ni, pmc, A(l(i,2)+1,lsc(ns-2)), lda, T(1+ni*pmc), ni) ! T2, blk col
                ! Use Random Sampling when nt is big
                IF( (nt2 .gt. 50) .and. (nt2 .gt. pmc+ppa) ) THEN
                   pp = pmc + ppa
                   call dgemm('N','N',pp,pmc,ni,1d0,Omega,pp,T(1+pmc*ni),ni,0d0,ww,pp) !ww=Omega*T, starting from 1.
                 ELSE
                   pp = nt2
                   call dlacpy('A',pp,pmc,T(1+pmc*nt2),pp,ww,pp)   ! copy T2 to ww
                END IF

                call dgesvd('N','A',pp,pmc,ww,pp,svs,ww,pp,ww(1+pmc*pp),pmc,work,lwork,info)  !svd of ww. 
                rkc = count(svs(1:pmc)>tol)
                ir2 = pmc - rkc   ! the number of compressed columns
                call dgemm('N','T',ni,rkc,pmc,1d0,T(1+pmc*ni),ni,ww(1+pmc*pp),pmc,0d0,T(1),ni) ! T2 is stored in T(1)
                call hssexpmm('V',PH,pnh,H,rkc,pmc,i,ww(1+pmc*pp),pmc,info) ! Vi^T is in ww(1+mi*pp), the first rkc rows
!
                tpn = 1+ni*pmc
                call dlacpy('A', pmr, ni, A(lsr(ns-2),l(i,2)+1), lda, T(tpn+ni*pmr), pmr) ! T1, blk row
                ! Use Random Sampling when nt is big
                IF( (ni .gt. 50) .and. (ni .gt. pmr+ppa) ) THEN
                   pp = pmr + ppa
                   call dgemm('N','N',pmr,pp,ni,1d0,T(tpn+ni*pmr),pmr,Omega,ni,0d0,ww,pmr)  ! ww = T*omega
                ELSE
                   pp = ni
                   call dlacpy('A',pmr,ni,T(tpn+ni*pmr),pmr,ww,pmr)   ! copy T1 to ww
                END IF

                call dgesvd('A','N',pmr,pp,ww,pmr,svs,ww(1+pmr*pp),pmr,ww,pp,work,lwork,info) 
                rkr = count(svs(1:pmr)>tol)
                ir1 = pmr - rkr   ! the number of compressed rows
                call dgemm('T','N',rkr,ni,pmr,1d0,ww(1+pmr*pp),pmr,T(tpn+ni*pmr),pmr,0d0,T(1+rkc*ni),rkr) 
                ! T1 is stored in T(1+rkc*nt)
                call hssexpmm('U',PH,pnh,H,pmr,rkr,i,ww(1+pmr*pp),pmr,info) ! U is in ww(1+mi*pp), the first rkr columns
!
             else  ! other parent nodes
!
                nt2 = lsr(ns-2)-1+lda-l(i,2)    ! HSS block column
                call dlacpy('A',lsr(ns-2)-1, pmc, A(1,wsr(ns-3)), lda, T(1+pmc*nt2),nt2)  ! T2
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A',ni,pmc, A(l(i,2)+1, lsc(ns-2)), lda, T(1+pmc*nt2+lsr(ns-2)-1), nt2)
                END IF
                ! Use Random Sampling when nt is big
                IF( (nt2 .gt. 50) .and. (nt2 .gt. pmc+ppa) ) THEN
                   pp = pmc + ppa
                   call dgemm('N','N',pp,pmc,nt2,1d0,Omega,pp,T(1+pmc*nt2),nt2,0d0,ww,pp)  ! ww=Omega*T, starting from 1.
                 ELSE
                   pp = nt2
                   call dlacpy('A',pp,pmc,T(1+pmc*nt2),pp,ww,pp)   ! copy T2 to ww
                END IF

                call dgesvd('N','A',pp,pmc,ww,pp,svs,ww,pp,ww(1+pmc*pp),pmc,work,lwork,info)   ! svd of ww. 
                rkc = count(svs(1:pmc)>tol)
                ir2 = pmc - rkc   ! the number of compressed columns
                call dgemm('N','T',nt2,rkc,pmc,1d0,T(1+pmc*nt2),nt2,ww(1+pmc*pp),pmc,0d0,T(1),nt2) ! T2 is stored in T(1)
                call hssexpmm('V',PH,pnh,H,rkc,pmc,i,ww(1+pmc*pp),pmc,info) ! Vi^T is in ww(1+mi*pp), the first rkc rows
                
!
                nt1 = lsc(ns-2)-1+lda-l(i,2)    ! HSS block column
                tpn = 1+rkc*nt2+pmr*nt1
                call dlacpy('A', pmr, lsc(ns-2)-1, A(wsc(ns-3), 1), lda, T(tpn), pmr) ! T1
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A', pmr, ni, A(lsr(ns-2),l(i,2)+1), lda, T( tpn+pmr*(lsc(ns-2)-1) ), pmr)
                END IF
                ! Use Random Sampling when nt is big
                IF( (nt1 .gt. 50) .and. (nt1 .gt. pmr+ppa) ) THEN
                   pp = pmr + ppa
                   call dgemm('N','N',pmr,pp,nt1,1d0,T(tpn),pmr,Omega,nt1,0d0,ww,pmr)  ! ww = T*omega
                ELSE
                   pp = nt1
                   call dlacpy('A',pmr,nt1,T(tpn),pmr,ww,pmr)   ! copy T1 to ww
                END IF

                call dgesvd('A','N',pmr,pp,ww,pmr,svs,ww(1+pmr*pp),pmr,ww,pmr,work,lwork,info) 
                rkr = count(svs(1:pmr)>tol)
                ir1 = pmr - rkr   ! the number of compressed rows
                call dgemm('T','N',rkr,nt1,pmr,1d0,ww(1+pmr*pp),pmr,T(tpn),pmr,0d0,T(1+rkc*nt2),rkr) 
                ! T1 is stored in T(1+rkc*nt)
                call hssexpmm('U',PH,pnh,H,pmr,rkr,i,ww(1+pmr*pp),pmr,info) ! U is in ww(1+pmr*pp), the first rkr columns

             endif  !(parent nodes)

                ns = ns-1  ! shrink ls
                lsr(ns) = lsr(ns-1)+rkr
                lsc(ns) = lsc(ns-1)+rkc
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A', ni, rkc, T(lsr(ns-1)), nt2, A(l(i,2)+1,lsc(ns-1)), lda)    ! T2 ??
                   call dlacpy('A', rkr, ni, T(1+nt2*rkc+rkr*(lsc(ns-1)-1) ), rkr,A(lsr(ns-1),l(i,2)+1),lda) ! T1 ??
                END IF
         endif
!
          do j = 1, ns-3
             A( lsr(j):lsr(j+1)-1,wsr(j)+ir2:wsr(ns-2)-1+ir2 ) = A( lsr(j):lsr(j+1)-1,wsr(j):wsr(ns-2)-1 )
             A( wsc(j)+ir1:wsc(ns-2)-1+ir1,lsc(j):lsc(j+1)-1 ) = A( wsc(j):wsc(ns-2)-1,lsc(j):lsc(j+1)-1 )
             wsr(j) = wsr(j) + ir2
             wsc(j) = wsc(j) + ir1
          enddo
          wsr(ns-1) = l(i,2)+1   ! Only for nonleaf nodes
          wsc(ns-1) = l(i,2)+1
          if (ns > 2 ) then
             call dlacpy('A',lsr(ns-1)-1,rkc,T(1),nt2,A(1,wsr(ns-2)+ir2), lda)          ! B{ch{tr(i)}(1)}
             call dlacpy('A',rkr,lsc(ns-1)-1,T(1+nt2*rkc),rkr, A(wsc(ns-2)+ir1,1), lda) ! B{ch{tr(i)}(2)}
             wsr(ns-2) = wsr(ns-2)+ir2
             wsc(ns-2) = wsc(ns-2)+ir1
          endif

       enddo

!   Root, Store B
       i=ltr
       ch1 = ch(1,i)
       lt  = lsr(ns-1)-lsr(ns-2)
       tpn = wsr(ns-1)-wsr(ns-2)
       call hssexpmm('B',PH,pnh,H,lt,tpn,ch1,A(lsr(ns-2),wsr(ns-2)), lda, info)  ! B{ch{tr{i}}(1)}
       ch1 = ch(2,i)
       lt  = wsc(ns-1)-wsc(ns-2)
       tpn = lsc(ns-1)-lsc(ns-2)
       call hssexpmm('B',PH,pnh,H,lt,tpn,ch1,A(wsc(ns-2),lsc(ns-2)), lda, info)  ! B{ch{tr{i}}(2)}

       deallocate (svs,wsr,wsc,lsc,ww,lsr,l,ch,work,omega)

     end subroutine mat2hssvdr

!!!!!!!
     subroutine mat2hssvd(A, LDA, TR, LTR, M, LM, PH, H, DD, TOL, RK, T, pnh )
       use aux_hss
!
!  .. Scalar Arguments ..
       INTEGER          :: LDA, RK, LTR, lm
       DOUBLE PRECISION :: TOL
!
!  .. Array Arguments ..
       DOUBLE PRECISION :: A(LDA, *), H(*), DD(*), T(*)
       INTEGER          :: TR(*), M(*)
       TYPE(HSSMM)      :: PH
! 
!  Purpose
! =========
! Construct an HSS matrix for a rank-structured matrix A
! 
! Parameters
!========== 
!  A    (INPUT) DOUBLE PRECISION array, DIMENSION(lda, lda)
!       We assume A is given explicitly and a modified version without forming A will
!       be written later. A is used as partly workspace and its entries will be changed.
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
!  RK   (INPUT) INTEGER, presetted rank for fixed-rank approximation
! 
!  T    (INPUT) DOUBLE PRECISION, DIMENSION( log2(lm)*mi* LDA ),
!       mi is the average of diagonal block size; 
!  
!===================
!  Written by Shengguo Li, on Aug. 8th, 2012
!
!  Modified on Aug. 14th, 2012
!  1) Use QR first when the matrix is tall and skiny
!===================

!  .. Local Scalars ..
       INTEGER       :: CH1, ERR, INFO, LWORK, pmr,pmc, rkc, rkr,i,j
       INTEGER       :: lt,it,ns,pnd,ni,pnh,mm,ir1,ir2,nt1,nt2,tpn,mi,nn,sz

!  .. Local Arrays ..
       INTEGER, ALLOCATABLE :: CH(:,:), L(:,:)
       INTEGER, ALLOCATABLE :: lsr(:), wsr(:), lsc(:), wsc(:)
       DOUBLE PRECISION, ALLOCATABLE :: WW(:), WORK(:), svs(:)
       ! WW is an extra workspace for random sampling. WORK is an extra workspace for DGESVD
!  ..
!  .. External Subroutines ..
      EXTERNAL  DGESVD, DGEMM, DLACPY
!  ..
!  .. External Functions .. 
      
!  .. Intrinsic Functions ..
      INTRINSIC    MIN, ABS,MAXVAL
              
       MM = MAXVAL( M(1:lm) )
       LWORK = LDA*MM*2 ! modify later
       pnh = 1              ! 'pointer' in H
       pnd = 1              ! 'pointer' in D
       allocate( CH(2,LTR),L(LTR,2),lsr(LTR), wsr(LTR), WW(LDA*LDA), WORK(LDA*MM*3) )
       allocate( svs(LDA),lsc(ltr),wsc(ltr), STAT=ERR ) 

       nn = ltr - 1 
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
       wsr(ns) = l(1,2)+1      ! column width at upper triangular
       lsc(ns) = 1             ! column start position at lower triangular
       wsc(ns) = l(1,2)+1      ! row width at lower triangular
       
       do i = 1, nn
          ! factorize (leaf) diag & get off-diag block

          if (ch(1,i) == 0 .and. ch(2,i) == 0) then 
             mi = l(i,2)-l(i,1)+1
             ni = lda-l(i,2)      !t^r, the right column block's length

             call dlacpy('A', mi, mi, A(l(i,1),l(i,1)), lda, DD(pnd), mi)    ! copy D{i}
             PH%D(i) = mi
             PH%pd(1,i) = pnd
             PH%pd(2,i) = pnd + mi*mi -1
             pnd = pnd + mi*mi

             if (ns == 1) then
                nt1 = ni    ! lenght of whole HSS block column
                nt2 = ni    ! lenght of whole HSS block row
                call dlacpy('A', ni, mi, A(l(i,2)+1,l(i,1)), lda, T(1+mi*ni), ni) ! Construct T, the HSS block column
                call dlacpy('A',ni,mi,T(1+mi*ni),ni,ww,ni)   ! copy T2 to ww
                
                IF( ni .gt. mi ) THEN
                   call dgeqrf(ni,mi,ww,ni,ww(1+mi*ni),work,lwork,info)
                   T(1:mi*mi) = 0
                   call dlacpy('U',mi,mi,ww,ni,T,mi)
                   call dgesvd('N','A',mi,mi,T,mi,svs,T,mi,ww(1+mi*ni),mi,work,lwork,info)
                ELSE
                   call dgesvd('N','A',ni,mi,ww,ni,svs,ww,ni,ww(1+mi*ni),mi,work,lwork,info)  !svd of ww. 
                END IF
                sz = min(mi, ni)  
                rkc = count(svs(1:mi)>tol)
                ir2 = mi - rkc   ! the number of compressed columns
                call dgemm('N','T',ni,rkc,mi,1d0,T(1+mi*ni),ni,ww(1+mi*ni),mi,0d0,T(1),ni) ! T2 is stored in T(1)
                call hssexpmm('V',PH,pnh,H,rkc,mi,i,ww(1+mi*ni),mi,info) ! Vi^T is in ww(1+mi*ni), the first rkc rows
!
                call dlacpy('A', mi, ni, A(l(i,1),l(i,2)+1), lda, T(1+(rkc+mi)*ni), mi) ! T1, HSS block row
                call dlacpy('A',mi,ni,T(1+(rkc+mi)*ni),mi,ww,mi)   ! copy T1 to ww
                
                call dgesvd('A','N',mi,ni,ww,mi,svs,ww(1+mi*ni),mi,ww,mi,work,lwork,info) 
                rkr = count(svs(1:sz)>tol)
                ir1 = mi - rkr   ! the number of compressed rows
                call dgemm('T','N',rkr,ni,mi,1d0,ww(1+mi*ni),mi,T(1+(rkc+mi)*ni),mi,0d0,T(1+rkc*ni),rkr) 
                ! T1 is stored in T(1+rkc*nt)
                call hssexpmm('U',PH,pnh,H,mi,rkr,i,ww(1+mi*ni),mi,info) ! U is in ww(1+mi*ni), the first rkr columns
!
             else   ! leaf nodes except node 1
!
                nt2 = lsr(ns)-1+ni        ! lenght of whole HSS block column
                call dlacpy('A', lsr(ns)-1, mi, A(1,wsr(ns-1)), lda, T(1+mi*nt2), nt2)
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A',ni,mi,A(l(i,2)+1,l(i,1)),lda,T(1+mi*nt2+lsr(ns)-1),nt2)     ! T2
                END IF
                call dlacpy('A',nt2,mi,T(1+mi*nt2),nt2,ww,nt2)   ! copy T2 to ww
                
                 IF( nt2 .gt. mi ) THEN
                   call dgeqrf(nt2,mi,ww,nt2,ww(1+mi*nt2),work,lwork,info)
                   T(1:mi*mi) = 0
                   call dlacpy('U',mi,mi,ww,nt2,T,mi)
                   call dgesvd('N','A',mi,mi,T,mi,svs,T,mi,ww(1+mi*nt2),mi,work,lwork,info)
                ELSE
                   call dgesvd('N','A',nt2,mi,ww,nt2,svs,ww,nt2,ww(1+mi*nt2),mi,work,lwork,info)   ! svd of ww. 
                END IF
                sz = min(mi, nt2)
                rkc = count(svs(1:sz)>tol)
                ir2 = mi - rkc   ! the number of compressed columns
                call dgemm('N','T',nt2,rkc,mi,1d0,T(1+mi*nt2),nt2,ww(1+mi*nt2),mi,0d0,T(1),nt2) ! T2 is stored in T(1)
                call hssexpmm('V',PH,pnh,H,rkc,mi,i,ww(1+mi*nt2),mi,info) ! Vi^T is in ww(1+mi*nt2), the first rkc rows
!
                nt1 = lsc(ns)-1+ni        ! length of whole HSS block row
                tpn = 1+rkc*nt2+mi*nt1
                call dlacpy('A', mi,lsc(ns)-1, A(wsc(ns-1),1),lda,T(tpn),mi)
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A', mi,ni,A(l(i,1),l(i,2)+1), lda,T( tpn+mi*(lsc(ns)-1) ), mi) ! T1
                END IF
                call dlacpy('A',mi,nt1,T(tpn),mi,ww,mi)   ! copy T1 to ww
                
                call dgesvd('A','N',mi,nt1,ww,mi,svs,ww(1+mi*nt1),mi,ww,mi,work,lwork,info) 
                sz = min(mi, nt1)
                rkr = count(svs(1:sz)>tol)
                ir1 = mi - rkr   ! the number of compressed rows
                call dgemm('T','N',rkr,nt1,mi,1d0,ww(1+mi*nt1),mi,T(tpn),mi,0d0,T(1+rkc*nt2),rkr) 
                ! T1 is stored in T(1+rkc*nt)
                call hssexpmm('U',PH,pnh,H,mi,rkr,i,ww(1+mi*nt1),mi,info) ! U is in ww(1+mi*nt1), the first rkr columns

             endif
             
               lsr(ns+1) = lsr(ns)+rkr
               wsr(ns)   = l(i,2)+1
               lsc(ns+1) = lsc(ns)+rkc
               wsc(ns)   = l(i,2)+1
               ns = ns +1
               IF( l(i,2)+1 .lt. lda ) THEN
                  call dlacpy('A', ni, rkc, T(lsr(ns-1)), nt2, A(l(i,2)+1,lsc(ns-1)), lda)    ! T2
                  call dlacpy('A', rkr, ni, T(1+nt2*rkc+rkr*(lsc(ns-1)-1) ), rkr,A(lsr(ns-1),l(i,2)+1),lda) ! T1 ??
               END IF
!               
          else  ! parent nodes
!
             ni = lda-l(i,2)
             pmr = lsr(ns)-lsr(ns-2)  ! row dim
             pmc = lsc(ns)-lsc(ns-2)  ! col dim

!            store B 
             ch1 = ch(1,i)
             lt  = lsr(ns-1)-lsr(ns-2)
             tpn = wsr(ns-1)-wsr(ns-2)
             call hssexpmm('B',PH,pnh,H,lt,tpn,ch1,A(lsr(ns-2),wsr(ns-2)), lda, info)  ! B{ch{tr{i}}(1)}
             ch1 = ch(2,i)
             lt  = wsc(ns-1)-wsc(ns-2)
             tpn = lsc(ns-1)-lsc(ns-2)
             call hssexpmm('B',PH,pnh,H,lt,tpn,ch1,A(wsc(ns-2),lsc(ns-2)), lda, info)  ! B{ch{tr{i}}(2)}

             if (ns == 3) then
                nt1 = ni    ! lenght of whole HSS block column
                nt2 = ni    ! lenght of whole HSS block row
                call dlacpy('A', ni, pmc, A(l(i,2)+1,lsc(ns-2)), lda, T(1+ni*pmc), ni) ! T2, blk col
                call dlacpy('A',ni,pmc,T(1+pmc*nt2),ni,ww,ni)   ! copy T2 to ww
                
                 IF( ni .gt. pmc ) THEN
                   call dgeqrf(ni,pmc,ww,ni,ww(1+pmc*ni),work,lwork,info)
                   T(1:pmc*pmc) = 0
                   call dlacpy('U',pmc,pmc,ww,ni,T,pmc)
                   call dgesvd('N','A',pmc,pmc,T,pmc,svs,T,pmc,ww(1+pmc*ni),pmc,work,lwork,info)
                ELSE
                   call dgesvd('N','A',ni,pmc,ww,ni,svs,ww,ni,ww(1+pmc*ni),pmc,work,lwork,info)  !svd of ww. 
                END IF
                sz = min(pmc,ni)
                rkc = count(svs(1:sz)>tol)
                ir2 = pmc - rkc   ! the number of compressed columns
                call dgemm('N','T',ni,rkc,pmc,1d0,T(1+pmc*ni),ni,ww(1+pmc*ni),pmc,0d0,T(1),ni) ! T2 is stored in T(1)
                call hssexpmm('V',PH,pnh,H,rkc,pmc,i,ww(1+pmc*ni),pmc,info) ! Vi^T is in ww(1+mi*ni), the first rkc rows
!
                tpn = 1+ni*pmc
                call dlacpy('A', pmr, ni, A(lsr(ns-2),l(i,2)+1), lda, T(tpn+ni*pmr), pmr) ! T1, blk row
                call dlacpy('A',pmr,ni,T(tpn+ni*pmr),pmr,ww,pmr)   ! copy T1 to ww

                call dgesvd('A','N',pmr,ni,ww,pmr,svs,ww(1+pmr*ni),pmr,ww,ni,work,lwork,info) 
                rkr = count(svs(1:sz)>tol)
                ir1 = pmr - rkr   ! the number of compressed rows
                call dgemm('T','N',rkr,ni,pmr,1d0,ww(1+pmr*ni),pmr,T(tpn+ni*pmr),pmr,0d0,T(1+rkc*ni),rkr) 
                ! T1 is stored in T(1+rkc*nt)
                call hssexpmm('U',PH,pnh,H,pmr,rkr,i,ww(1+pmr*ni),pmr,info) ! U is in ww(1+mi*ni), the first rkr columns
!
             else  ! other parent nodes
!
                nt2 = lsr(ns-2)-1+lda-l(i,2)    ! HSS block column
                call dlacpy('A',lsr(ns-2)-1, pmc, A(1,wsr(ns-3)), lda, T(1+pmc*nt2),nt2)  ! T2
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A',ni,pmc, A(l(i,2)+1, lsc(ns-2)), lda, T(1+pmc*nt2+lsr(ns-2)-1), nt2)
                END IF
                call dlacpy('A',nt2,pmc,T(1+pmc*nt2),nt2,ww,nt2)   ! copy T2 to ww
                
                IF( nt2 .gt. pmc ) THEN
                   call dgeqrf(nt2,pmc,ww,nt2,ww(1+pmc*nt2),work,lwork,info)
                   T(1:pmc*pmc) = 0
                   call dlacpy('U',pmc,pmc,ww,nt2,T,pmc)
                   call dgesvd('N','A',pmc,pmc,T,pmc,svs,T,pmc,ww(1+pmc*nt2),pmc,work,lwork,info)
                ELSE
                   call dgesvd('N','A',nt2,pmc,ww,nt2,svs,ww,nt2,ww(1+pmc*nt2),pmc,work,lwork,info)   ! svd of ww. 
                END IF
                sz = min(pmc,nt2)
                rkc = count(svs(1:sz)>tol)
                ir2 = pmc - rkc   ! the number of compressed columns
                call dgemm('N','T',nt2,rkc,pmc,1d0,T(1+pmc*nt2),nt2,ww(1+pmc*nt2),pmc,0d0,T(1),nt2) ! T2 is stored in T(1)
                call hssexpmm('V',PH,pnh,H,rkc,pmc,i,ww(1+pmc*nt2),pmc,info) ! Vi^T is in ww(1+mi*nt2), the first rkc rows
                
!
                nt1 = lsc(ns-2)-1+lda-l(i,2)    ! HSS block row
                tpn = 1+rkc*nt2+pmr*nt1
                call dlacpy('A', pmr, lsc(ns-2)-1, A(wsc(ns-3), 1), lda, T(tpn), pmr) ! T1
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A', pmr, ni, A(lsr(ns-2),l(i,2)+1), lda, T( tpn+pmr*(lsc(ns-2)-1) ), pmr)
                END IF
                call dlacpy('A',pmr,nt1,T(tpn),pmr,ww,pmr)   ! copy T1 to ww

                call dgesvd('A','N',pmr,nt1,ww,pmr,svs,ww(1+pmr*nt1),pmr,ww,pmr,work,lwork,info) 
                sz = min(pmr,nt1)
                rkr = count(svs(1:sz)>tol)
                ir1 = pmr - rkr   ! the number of compressed rows
                call dgemm('T','N',rkr,nt1,pmr,1d0,ww(1+pmr*nt1),pmr,T(tpn),pmr,0d0,T(1+rkc*nt2),rkr) 
                ! T1 is stored in T(1+rkc*nt)
                call hssexpmm('U',PH,pnh,H,pmr,rkr,i,ww(1+pmr*nt1),pmr,info) ! U is in ww(1+pmr*nt1), the first rkr columns

             endif  !(parent nodes)

                ns = ns-1  ! shrink ls
                lsr(ns) = lsr(ns-1)+rkr
                lsc(ns) = lsc(ns-1)+rkc
                IF( l(i,2)+1 .lt. lda ) THEN
                   call dlacpy('A', ni, rkc, T(lsr(ns-1)), nt2, A(l(i,2)+1,lsc(ns-1)), lda)    ! T2 ??
                   call dlacpy('A', rkr, ni, T(1+nt2*rkc+rkr*(lsc(ns-1)-1) ), rkr,A(lsr(ns-1),l(i,2)+1),lda) ! T1 ??
                END IF
         endif
!
          do j = 1, ns-3
             A( lsr(j):lsr(j+1)-1,wsr(j)+ir2:wsr(ns-2)-1+ir2 ) = A( lsr(j):lsr(j+1)-1,wsr(j):wsr(ns-2)-1 )
             A( wsc(j)+ir1:wsc(ns-2)-1+ir1,lsc(j):lsc(j+1)-1 ) = A( wsc(j):wsc(ns-2)-1,lsc(j):lsc(j+1)-1 )
             wsr(j) = wsr(j) + ir2
             wsc(j) = wsc(j) + ir1
          enddo
          wsr(ns-1) = l(i,2)+1   ! Only for nonleaf nodes
          wsc(ns-1) = l(i,2)+1
          if (ns > 2 ) then
             call dlacpy('A',lsr(ns-1)-1,rkc,T(1),nt2,A(1,wsr(ns-2)+ir2), lda)          ! B{ch{tr(i)}(1)}
             call dlacpy('A',rkr,lsc(ns-1)-1,T(1+nt2*rkc),rkr, A(wsc(ns-2)+ir1,1), lda) ! B{ch{tr(i)}(2)}
             wsr(ns-2) = wsr(ns-2)+ir2
             wsc(ns-2) = wsc(ns-2)+ir1
          endif

       enddo

!   Root, Store B
       i=ltr
       ch1 = ch(1,i)
       lt  = lsr(ns-1)-lsr(ns-2)
       tpn = wsr(ns-1)-wsr(ns-2)
       call hssexpmm('B',PH,pnh,H,lt,tpn,ch1,A(lsr(ns-2),wsr(ns-2)), lda, info)  ! B{ch{tr{i}}(1)}
       ch1 = ch(2,i)
       lt  = wsc(ns-1)-wsc(ns-2)
       tpn = lsc(ns-1)-lsc(ns-2)
       call hssexpmm('B',PH,pnh,H,lt,tpn,ch1,A(wsc(ns-2),lsc(ns-2)), lda, info)  ! B{ch{tr{i}}(2)}

       deallocate (svs,wsr,wsc,lsc,ww,lsr,l,ch, work)

     end subroutine mat2hssvd

!!!!!!!!
  subroutine fastHssmmL(DD,H,PH,M,X,TR,LTR,LDX,NCOL,TRE,LVL,WW,SX)
    use aux_hss
!
!  .. Scalar Arguments ..
       INTEGER           ::  LTR, LDX, NCOL, LVL

!  .. Array Arguments ..
       DOUBLE PRECISION  ::  DD(*), H(*), X(LDX,*), WW(LDX,*), SX(*)
       INTEGER           ::  M(*), TR(*)
       TYPE(HSSMM)       ::  PH  
       TYPE(HTinfo)      ::  TRE
!
!  Purpose    
!  ========
!   Compute the multiplication of a general matrix with an HSS matrix. This HSS matrix is in post ordering. 
!   This algorithm is a sequential code and based on the assumption that row and column partitions are same.
!   
! .. Argumrents ..
! DD    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators D's.
! 
! H    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators, U, V and B's. 
! 
! PH   (INPUT)  TYPE HSSMM, stores the information of all these generators.
! 
! M    (INPUT)  INTEGER array, DIMENSION(:), stores the sizes of diagonal blocks.
! 
! X    (INPUT/OUTPUT) DOUBLE PRECISION array, DIMENSION(:,:), the input matix.
!
! TR    (INPUT) INTEGER array, DIMENSION(ltr), The HSS tree in post ordering.
!
! LTR  (INPUT)  INTEGER, the length of HSS Tree, LTR.
! 
! LDX  (INPUT) INTEGER, the row dimension of X
!
! NCOL (INPUT) INTEGER, the column dimension of X
!
! TRE  (INPUT) TYPE HTinfo, contains some information about HSS tree
! 
! LVL  (INPUT) INTEGER, the total levels of HSS tree
! 
! ====================
! Written by Shengguo Li, on Nov 20th, 2012
!
! ====================

! .. Local Parameters ..
       INTEGER    :: n1, i, it, lt, nj, flag, LL, pww, ldw, &
                     pdx, strid_m, j, k, leafi, ni, pnx, nstart1, nstart2

! .. Array Parameters ..
       INTEGER    :: lx(ltr), l(ltr), px(ltr)

! .. Executable Statements ..
       n1 = 0.5*(ltr+1)
       ldw = ldx

!*********************************************************
!                      leaf nodes                        *
!*********************************************************
! lx(i): row partitions of x
       lt = 1
       it = 1
       lx = 0
       l  = 0
       n1 = TRE%lentr(lvl)
       do i = 1, n1
          lx(i) = lt
          lt = lt+m(i)

          ni = TRE%ttr(lvl,i)
          l(ni) = it
          it = it+PH%U(2,ni)
       end do

       pnx = 1
       do i = 1, n1
          ni = TRE%ttr(lvl,i)  ! can change it to column major
          call dgemm('N','N',ldx,PH%U(2,ni),m(i),1d0,X(1,lx(i)),ldx,H( PH%ph(1,ni) ),PH%U(1,ni), & 
               0d0,ww(1,pnx),ldw)          ! ww = Xi * U(i)
          pnx = pnx + PH%U(2,ni)
       end do

       pdx = 1
       do i  = 1, n1
          ni = TRE%ttr(lvl, i)
          
          if (ni .eq. TRE%ch(1,tr(ni) ) ) then 
             nj = TRE%ch(2, tr(ni) )    ! adjacent of ni
          else
             nj = TRE%ch(1, tr(ni) )
          end if
          
          ! xt{ni}: f_i in note
          call dgemm('N','N',ldw,PH%B(2,nj),PH%B(1,nj),1d0,ww(1,l(nj)),ldw,H( PH%ph(3,nj) ),&
               PH%B(1,nj),0d0, SX(pdx),ldx) 
          px(ni) = pdx
          pdx = pdx + PH%B(2,nj)*ldx
       end do

!***************************************************************
!                    Forward Traversal                         *
!***************************************************************
       flag = 0
       pww = 0
       DO LL = lvl-1, 1, -1
          n1 = TRE%lentr(LL)
          
          ! merge and determine the column dim
          lt = 1
          it = 1
          DO i =1, n1
             ni = TRE%ttr(ll, i)
             IF ( TRE%ch(1,ni) .ne. 0  ) THEN
                lx(i)  = lt
                lt = lt+PH%U( 2, TRE%ch(1,ni) )+PH%U( 2, TRE%ch(2,ni) )
             ELSE
                if( flag .eq. 0) then
                   flag  = 1
                   leafi = i   ! leafi is a leaf node, above is full binary tree
                end if
             END IF
             
             l(ni) = it    ! l(ni) records the starting column of ww(ni)
             it = it+PH%U(2,ni)
          END DO !(i)

          ! update ww i.e. g_i in our note
          Nstart1 = max( pww * ncol, 1)
          Nstart2 = max( (1-pww) * ncol, 1)
          IF ( flag .eq. 0 ) THEN          !there is no leaf node except bottom level
             DO i = 1, n1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',ldx,PH%U(2,ni),PH%U(1,ni),1d0,ww(1,Nstart1+lx(i)-1),ldw,H( PH%ph(1,ni) ), & 
                     PH%U(1,ni),0d0,ww(1,Nstart2+l(ni)-1),ldx)
             END DO
          ELSE
             DO i = 1, leafi-1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',ldx,PH%U(2,ni),PH%U(1,ni),1d0,ww(1,Nstart1+lx(i)-1),ldw,H( PH%ph(1,ni) ), & 
                     PH%U(1,ni),0d0,ww(1,Nstart2+l(ni)-1),ldx)
             END DO
             
             DO i = leafi, n1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',ldx,PH%U(2,ni),PH%U(1,ni),1d0,X(1,TRE%lfx(ni)),ldx,H( PH%ph(1,ni) ), & 
                     PH%U(1,ni), 0d0,ww(1,Nstart2+l(ni)-1), ldw )
             END DO
             flag = 0
          END IF !(flag)          

          do i = 1, n1
             ni = TRE%ttr(LL, i)
          
             if (ni .eq. TRE%ch(1,tr(ni) ) ) then 
                nj = TRE%ch(2, tr(ni) )    ! adjacent of ni
             else
                nj = TRE%ch(1, tr(ni) )
             end if
          
             ! xt{ni}
             call dgemm('N','N',ldx,PH%B(2,nj),PH%B(1,nj),1d0,ww(1,Nstart2+l(nj)-1),ldw,H( PH%ph(3,nj) ), &
                  PH%B(1,nj), 0d0,SX(pdx), ldx )      
             px(ni) = pdx
             pdx = pdx + PH%B(2,nj)*ldx
          end do
          pww = 1 - pww

       END DO ! (LL)

!*********************************************************
!                 Backward Traversal                     *
!*********************************************************
       ! First level
       LL = 1
       lt = 1
       DO i = 1, 2
          ni = TRE%ttr(LL, i)
          call dgemm('N','N',ldx,PH%V(2,ni),PH%V(1,ni),1d0,SX( px(ni) ),ldx, H( PH%ph(2,ni) ), & 
               PH%V(1,ni),0d0, ww(1,lt),ldw)
          
          IF ( TRE%ch(1,ni) .eq. 0 ) THEN   ! add it to X
             call dgemm('N','N',ldx,PH%D(ni),PH%D(ni),1d0,X(1,TRE%lfx(ni)),ldx,DD( PH%pd(1,ni) ),PH%D(ni), & 
                  1d0, ww(1,lt), ldw )
             call dlacpy('A', ldx,PH%D(ni),ww(1,lt),ldw,X(1,TRE%lfx(ni)), ldx)
          END IF
          lt = lt + PH%V(2,ni)
       END DO
       
       ! Traverse downward
       DO LL = 2, lvl
          n1 = TRE%lentr(LL)

          ! update the data of children
          lt =1
          DO i = 1, n1
             ni = TRE%ttr(LL,i)

             strid_m = PH%V(1,ni)
             it = px(ni)-1
             DO k = 1, strid_m
                DO j = 1, ldx         ! use 2D array there is a huge stride.
                   SX(it+j+(k-1)*ldx ) = SX(it+j+(k-1)*ldx ) + ww(j,lt+k-1)
                END DO
             END DO ! matrix-matrix addition
             lt = lt+strid_m
          END DO
!
          ! extend
          lt = 1
          DO i = 1, n1
             ni = TRE%ttr(LL, i)
             call dgemm('N','N',ldx,PH%V(2,ni),PH%V(1,ni),1d0, SX( px(ni) ), ldx, H( PH%ph(2,ni) ),& 
                  PH%V(1,ni),0d0, ww(1,lt),ldw )
             
             IF( TRE%ch(1,ni) .eq. 0 ) THEN   ! add it to X
                call dgemm('N','N',ldx,PH%D(ni),PH%D(ni),1d0,X(1,TRE%lfx(ni)),ldx,DD( PH%pd(1,ni) ),PH%D(ni), & 
                     1d0, ww(1,lt), ldw )
                call dlacpy('A', ldx,PH%D(ni),ww(1,lt),ldw,X(1,TRE%lfx(ni)), ldx)
             END IF             
             lt = lt + PH%V(2,ni)
          END DO

       END DO ! (LL)

     end subroutine fastHssmmL

!!!!!!!!
  subroutine fastHssmm(DD,H,PH,M,X,TR,LTR,LDX,NCOL,TRE,LVL,WW,SX)
    use aux_hss
!
!  .. Scalar Arguments ..
       INTEGER           ::  LTR, LDX, NCOL, LVL

!  .. Array Arguments ..
       DOUBLE PRECISION  ::  DD(*), H(*), X(LDX,*), WW(2*LDX,*), SX(*)
       INTEGER           ::  M(*), TR(*)
       TYPE(HSSMM)       ::  PH  
       TYPE(HTinfo)      ::  TRE
!
!  Purpose    
!  ========
!   Compute HSS matrix with a general matrix. This HSS matrix is in post ordering. 
!   
! .. Argumrents ..
! DD    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators D's.
! 
! H    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators, U, V and B's. 
! 
! PH   (INPUT)  TYPE HSSMM, stores the information of all these generators.
! 
! M    (INPUT)  INTEGER array, DIMENSION(:), stores the sizes of diagonal blocks.
! 
! X    (INPUT/OUTPUT) DOUBLE PRECISION array, DIMENSION(:,:), the input matix.
!
! TR    (INPUT) INTEGER array, DIMENSION(ltr), The HSS tree in post ordering.
!
! LTR  (INPUT)  INTEGER, the length of HSS Tree, LTR.
! 
! LDX  (INPUT) INTEGER, the row dimension of X
!
! NCOL (INPUT) INTEGER, the column dimension of X
!
! TRE  (INPUT) TYPE HTinfo, contains some information about HSS tree
! 
! LVL  (INPUT) INTEGER, the total levels of HSS tree
! 
! ====================
! Written by Shengguo Li, on Aug 9th, 2012
!
! Modified on Aug. 10th, 2012.
!   1) ww may be changed to 1D array, and the indexing would be much more complicated.
!   2) Some indexing l and lx may be done duing processing. Modify it later and make it work first.  
! ====================

! .. Local Parameters ..
       INTEGER    :: n1, i, it, lt, nj, flag, LL, pww, ldw, &
                     pdx, strid_m, j, k, leafi, ni, pnx, nstart1, nstart2

! .. Array Parameters ..
       INTEGER     :: lx(ltr), l(ltr), px(ltr)
!
! .. Executable Statements ..
       n1 = 0.5*(ltr+1)
       ldw = 2*ldx

!*********************************************************
!                      leaf nodes                        *
!*********************************************************
! lx(i): row partition of x
       lt = 1
       it = 1
       lx = 0
       l  = 0
       n1 = TRE%lentr(lvl)
       do i = 1, n1
          lx(i) = lt
          lt = lt+m(i)

          ni = TRE%ttr(lvl, i)
          l(ni) = it
          it = it+PH%V(1,ni)
       end do

       pnx = 1
       do i = 1, n1
          ni = TRE%ttr(lvl,i)  ! can change it in column major
          call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0,H( PH%ph(2,ni) ),PH%V(1,ni), & 
               X(lx(i),1), ldx, 0d0,ww(pnx,1), ldw)
          ! ww = V{i}'*Xi
          pnx = pnx + PH%V(1,ni)
       end do
!
       pdx = 1
       do i = 1, n1
          ni = TRE%ttr(lvl, i)
          
          if (ni .eq. TRE%ch(1,tr(ni) ) ) then 
             nj = TRE%ch(2, tr(ni) )    ! adjacent of ni
          else
             nj = TRE%ch(1, tr(ni) )
          end if
          
          ! xt{ni}
          call dgemm('N','N',PH%B(1,ni),ncol,PH%B(2,ni),1d0,H( PH%ph(3,ni) ),&
               PH%B(1,ni),ww(l(nj),1),ldw, 0d0, SX(pdx),PH%B(1,ni)) 
          px(ni) = pdx
          pdx = pdx + PH%B(1,ni)*ncol 
       end do

!***************************************************************
!                    Forward Traversal                         *
!***************************************************************
       flag = 0
       pww = 0
       DO LL = lvl-1, 1, -1
          n1 = TRE%lentr(LL)
          
          ! merge and determine the row dim
          lt = 1
          it = 1
          DO i =1, n1
             ni = TRE%ttr(ll, i)
             IF ( TRE%ch(1,ni) .ne. 0  ) THEN
                lx(i)  = lt
                lt = lt+PH%V( 1, TRE%ch(1,ni) )+PH%V( 1, TRE%ch(2,ni) )
             ELSE
                if( flag .eq. 0) then
                   flag  = 1
                   leafi = i   ! leafi is a leaf node, above is full binary tree
                end if
             END IF
             
             l(ni) = it    ! l(ni) records the starting row of ww(ni)
             it = it+PH%V(1,ni)
          END DO !(i)

          ! update ww
          Nstart1 = max( pww * ldx, 1)
          Nstart2 = max( (1-pww) * ldx, 1)
          IF ( flag .eq. 0 ) THEN
             DO i = 1, n1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0,H( PH%ph(2,ni) ), & 
                     PH%V(1,ni),ww(Nstart1+lx(i)-1,1),ldw,0d0,ww(Nstart2+l(ni)-1,1),ldw)
             END DO
          ELSE
             DO i = 1, leafi-1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0,H( PH%ph(2,ni) ), & 
                     PH%V(1,ni),ww(Nstart1+lx(i)-1,1),ldw,0d0,ww(Nstart2+l(ni)-1,1),ldw)
             END DO
             
             DO i = leafi, n1
                ni = TRE%ttr(LL, i)
                call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0, H( PH%ph(2,ni) ), & 
                     PH%V(1,ni), X(TRE%lfx(ni),1), ldx, 0d0,ww(Nstart2+l(ni)-1,1), ldw )
             END DO
             flag = 0
          END IF !(flag)          

          do i = 1, n1
             ni = TRE%ttr(LL, i)
          
             if (ni .eq. TRE%ch(1,tr(ni) ) ) then 
                nj = TRE%ch(2, tr(ni) )    ! adjacent of ni
             else
                nj = TRE%ch(1, tr(ni) )
             end if
          
             ! xt{ni}
             call dgemm('N','N',PH%B(1,ni),ncol,PH%B(2,ni),1d0,H( PH%ph(3,ni) ), &
                  PH%B(1,ni),ww(Nstart2+l(nj)-1,1), ldw, 0d0,SX(pdx), PH%B(1,ni) )      
             px(ni) = pdx
             pdx = pdx + PH%B(1,ni)*ncol 
          end do
          pww = 1 - pww

       END DO ! (LL)

!*********************************************************
!                 Backward Traversal                     *
!*********************************************************

       ! First level
       LL = 1
       lt = 1
       DO i = 1, 2
          ni = TRE%ttr(LL, i)
          call dgemm('N','N',PH%U(1,ni),ncol,PH%U(2,ni),1d0, H( PH%ph(1,ni) ), & 
               PH%U(1,ni),SX( px(ni) ), PH%U(2,ni),0d0, ww(lt,1),ldw)
          
          IF ( TRE%ch(1,ni) .eq. 0 ) THEN   ! add it to X
             call dgemm('N','N',PH%D(ni),ncol,PH%D(ni), 1d0, DD( PH%pd(1,ni) ), PH%D(ni), & 
                  X(TRE%lfx(ni), 1),ldx, 1d0, ww(lt,1), ldw )
             call dlacpy('A', PH%D(ni),ncol,ww(lt,1),ldw,X(TRE%lfx(ni),1), ldx)
          END IF
          lt = lt + PH%U(1,ni)

       END DO
       
       ! Traverse downward
       DO LL = 2, lvl
          n1 = TRE%lentr(LL)

          ! update the data of children
          lt =1
          DO i = 1, n1
             ni = TRE%ttr(LL,i)

             strid_m = PH%B(1,ni)
             it = px(ni)-1
             DO k = 1, ncol
                DO j = 1, strid_m         ! use 2D array there is a huge stride.
                   SX(it+j+(k-1)*strid_m ) = SX(it+j+(k-1)*strid_m ) + ww(lt+j-1,k)
                END DO
             END DO ! matrix-matrix addition
             lt = lt+strid_m
          END DO
!
          ! extend
          lt = 1
          DO i = 1, n1
             ni = TRE%ttr(LL, i)
             call dgemm('N','N',PH%U(1,ni),ncol,PH%U(2,ni),1d0, H( PH%ph(1,ni) ),& 
                  PH%U(1,ni),SX( px(ni) ), PH%U(2,ni),0d0, ww(lt,1),ldw )
             
             IF( TRE%ch(1,ni) .eq. 0 ) THEN   ! add it to X
                call dgemm('N','N',PH%D(ni),ncol,PH%D(ni), 1d0, DD( PH%pd(1,ni) ), PH%D(ni), & 
                     X(TRE%lfx(ni), 1),ldx, 1d0, ww(lt,1), ldw )
                call dlacpy('A', PH%D(ni),ncol,ww(lt,1),ldw,X(TRE%lfx(ni),1), ldx)
             END IF             
             lt = lt + PH%U(1,ni)
          END DO

       END DO ! (LL)

     end subroutine fastHssmm


!!!!!!!
     SUBROUTINE fastHssmm_omp( DD,H,PH,M,X,TR,LTR,LDX,NCOL,TRE,LVL )
     USE aux_hss
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: LTR, LDX, NCOL, LVL
!
! .. Array Arguments ..
    INTEGER, INTENT(IN) ::  M(*), TR(*)
    TYPE(HSSMM),  INTENT(IN) ::  PH  
    TYPE(HTinfo), INTENT(IN) ::  TRE
    DOUBLE PRECISION, INTENT(IN) :: DD(*), H(*)
    DOUBLE PRECISION, INTENT(INOUT) :: X(LDX,*)
!
!  Purpose    
!  ========
!  Compute HSS matrix with a general matrix. This HSS matrix is in post ordering. 
!  This subroutine is written in OpenMP. We use OpenMP to speed up its computation. 
!  
!
!*****************************************
! HSS Tree must have at least two levels *
!*****************************************
!
! ..Argumrents..
! ==============
! DD   (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators D's.
! 
! H    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators, U, V and B's. 
! 
! PH   (INPUT)  TYPE HSSMM, stores the information of all these generators.
! 
! M    (INPUT)  INTEGER array, DIMENSION(:), stores the sizes of diagonal blocks.
! 
! X    (INPUT/OUTPUT) DOUBLE PRECISION array, DIMENSION(:,:), the input matix.
!
! TR   (INPUT) INTEGER array, DIMENSION(ltr), The HSS tree in post ordering.
!
! LTR  (INPUT)  INTEGER, the length of HSS Tree, LTR.
! 
! LDX  (INPUT) INTEGER, the row dimension of X
!
! NCOL (INPUT) INTEGER, the column dimension of X
!
! TRE  (INPUT) TYPE HTinfo, contains some information about HSS tree
! 
! LVL  (INPUT) INTEGER, the total levels of HSS tree
! 
! ====================
! Written by Shengguo Li, on Mar. 18th, 2013
! ==========================================
!
! ..Local Parameters..
    INTEGER   n1, nc1, i, it, lt, nj, LL, ii, ierr, &
         pdx, j, k, ni, LDG, LDF,nxi
!
! ..Array Parameters..
    INTEGER    :: lx(ltr),lg(ltr),lf(ltr)
    DOUBLE PRECISION, PARAMETER   :: ZERO = 0.0D0
    DOUBLE PRECISION, allocatable ::  G(:,:), F(:,:)
! 
! Check input variables
    IF( lvl .lt. 3) Then
       write(*,*) "HSS tree must have at least 3 levels"
       return
    END IF

! .. Executable Statements ..
    n1 = 0.5*(ltr+1)

!*********************************************************
!                      partition x                       *
!*********************************************************
! lx(i): row partition of x, only for leaf nodes;
! lg(i): starting row position for each G_i;
! lf(i): starting row position for each F_i;
       lt = 1
       lx = 0
       i = 1
       n1 = TRE%lentr( lvl )
       DO ni = 1, ltr
          nc1 = TRE%ch( 1,ni )
          IF( nc1 .eq. 0 ) THEN
             lx(ni) = lt
             lt = lt + m(i)
             i = i + 1
          END IF
       END DO

!*********************************************************
!                   positions of G                       *
!*********************************************************
       it = 1
       lg = 0
       DO LL = lvl, 1, -1  ! each level
          n1 = TRE%lentr( LL )
          DO i = 1, n1
             ni = TRE%ttr( LL,i )
             lg( ni ) = it
             it = it+PH%V( 1,ni )
          END DO
       END DO
       LDG = it

!*********************************************************
!                   positions of F                       *
!*********************************************************
       it = 1
       lf = 0
       DO LL = 1, lvl  ! each level
          n1 = TRE%lentr( LL )
          DO i = 1, n1
             ni = TRE%ttr( LL,i )
             lf( ni ) = it
             it = it+PH%B( 1,ni )
          END DO
       END DO
       LDF = it

       ALLOCATE( G(LDG,NCOL), F(LDF,NCOL), stat=ierr )
       IF(ierr /= 0 ) THEN
          write(*,*) "Allocate failed in fasthssmm_omp "
       END IF
       G = zero
       F = zero

!*********************************************************
!                    Upsweep for G                       *
!*********************************************************
       ! bottom level
       n1 = TRE%lentr( lvl )

!$OMP PARALLEL PRIVATE(i,ni)
!$OMP DO SCHEDULE(dynamic)
       DO i = 1, n1
          ni = TRE%ttr( lvl,i )
          call dgemm('N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1.0d0,H( PH%ph(2,ni) ),PH%V(1,ni), & 
               X(lx(ni),1), LDX, 0.0d0,G(lg(ni),1), LDG )  ! G{i} = V{i}'*Xi
       END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

       DO LL = lvl-1, 1, -1
          n1 = TRE%lentr( LL )
!$OMP PARALLEL PRIVATE(i,nc1,ni)
!$OMP DO SCHEDULE(dynamic)
          DO i = 1, n1
             ni  = TRE%ttr( LL,i )
             nc1 = TRE%ch( 1,ni )
             IF( nc1 .ne. 0 ) THEN  ! parent node
                call dgemm( 'N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0,H( PH%ph(2,ni) ), & 
                    PH%V(1,ni),G(lg(nc1),1),LDG,0.0d0,G(lg(ni),1),LDG ) ! G{i}=V{i}'*[G{i1}; G{i2}]
             ELSE ! leaf node
                call dgemm( 'N','N',PH%V(1,ni),NCOL,PH%V(2,ni),1d0, H( PH%ph(2,ni) ), & 
                     PH%V(1,ni),X(lx(ni),1),LDX,0.0d0,G(lg(ni),1),LDG )   ! G{i} = V{i}'*Xi
             END IF
          END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       END DO

!*********************************************************
!                 Backward Traversal                     *
!*********************************************************
       ! the second top level
       LL = 1
       lt = 1
!!$OMP PARALLEL PRIVATE(i,ni,nj,nc1)
!!$OMP DO SCHEDULE(dynamic)
       DO i = 1, 2
          ni = TRE%ttr( LL,i )
          IF( ni .EQ. TRE%ch( 1,tr(ni) ) ) THEN
             nj = TRE%ch( 2,tr(ni) )    ! adjacent of ni
          ELSE
             nj = TRE%ch( 1,tr(ni) )
          END IF
          call dgemm('N','N',PH%B(1,ni),ncol,PH%B(2,ni),1.0d0,H( PH%ph(3,ni) ), &
               PH%B(1,ni),G(lg(nj),1), LDG, 0.0d0,F( lf(ni),1 ), LDF ) ! F{i}= B{i}*G{j}

          ! update F of its children
          nc1 = TRE%ch(1,ni)      ! left child
          IF( nc1 .ne. 0 ) THEN
             call dgemm( 'N','N',PH%U(1,ni),ncol,PH%U(2,ni),1d0,H( PH%ph(1,ni) ), &
                  PH%U(1,ni),F(lf(ni),1), LDF, 0.0d0,F(lf(nc1),1), LDF )
          END IF
       END DO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL
       
       ! Traverse downward
       DO LL = 2, lvl
          n1 = TRE%lentr( LL )

!$OMP PARALLEL PRIVATE(i,ni,nj,nc1)
!$OMP DO SCHEDULE(dynamic)
          DO i = 1, n1
             ni = TRE%ttr( LL,i )
             IF( ni .EQ. TRE%ch(1,tr(ni) ) ) THEN
                nj = TRE%ch( 2,tr(ni) )      ! adjacent of ni
             ELSE
                nj = TRE%ch( 1,tr(ni) )
             END IF
             call dgemm( 'N','N',PH%B(1,ni),ncol,PH%B(2,ni),1.0d0,H( PH%ph(3,ni) ), &
                  PH%B(1,ni),G(lg(nj),1), LDG, 1.0D0,F(lf(ni),1), LDF )
             
             nc1 = TRE%ch( 1,ni )    ! update F of its children
             IF( nc1 .ne. 0 ) THEN
                call dgemm( 'N','N',PH%U(1,ni),ncol,PH%U(2,ni),1d0,H( PH%ph(1,ni) ), &
                     PH%U(1,ni),F(lf(ni),1), LDF, 0.0d0,F(lf(nc1),1), LDF )
             END IF
          END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       END DO ! (LL)

!*********************************************************
!                   Compute X                            *
!*********************************************************
       DO ni = 1, ltr-1
          nc1 = TRE%ch( 1,ni )
          IF( nc1 .eq. 0 ) THEN
             call dgemm('N','N',PH%U(1,ni),ncol,PH%U(2,ni),1d0, H( PH%ph(1,ni) ),& 
                  PH%U(1,ni),F(lf(ni),1),LDF,0d0, G(1,1),LDG )
             
             call dgemm('N','N',PH%D(ni),ncol,PH%D(ni), 1d0, DD( PH%pd(1,ni) ), PH%D(ni), & 
                  X(lx(ni), 1),LDX, 1.0d0, G(1,1), LDG )
             call dlacpy( 'A', PH%D(ni),ncol,G(1,1),LDG,X(lx(ni),1), LDX )
          END IF
       END DO

       DEALLOCATE( G, F )

     END SUBROUTINE fastHssmm_omp


!!!!!!!!
     SUBROUTINE fastHssmmL_omp( DD,H,PH,M,X,TR,LTR,LDX,NCOL,TRE,LVL )
       USE aux_hss
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: LTR, LDX, NCOL, LVL
!
! .. Array Arguments ..
    INTEGER, INTENT(IN) ::  M(*), TR(*)
    TYPE(HSSMM),  INTENT(IN) ::  PH  
    TYPE(HTinfo), INTENT(IN) ::  TRE
    DOUBLE PRECISION, INTENT(IN) :: DD(*),H(*)
    DOUBLE PRECISION, INTENT(INOUT) :: X(LDX,NCOL)
!
!  Purpose    
!  ========
!  Multiplication a general matrix with an HSS matrix. This HSS matrix is in post ordering. 
!  This subroutine is written in OpenMP. We use OpenMP to speed up its computation. 
!
!*****************************************
! HSS Tree must have at least two levels *
!*****************************************
!
! ..Argumrents..
! ==============
! DD   (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators D's.
! 
! H    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the generators, U, V and B's. 
! 
! PH   (INPUT)  TYPE HSSMM, stores the information of all these generators.
! 
! M    (INPUT)  INTEGER array, DIMENSION(:), stores the sizes of diagonal blocks.
! 
! X    (INPUT/OUTPUT) DOUBLE PRECISION array, DIMENSION(:,:), the input matix.
!       X is usually a fat matrix.
!
! TR   (INPUT) INTEGER array, DIMENSION(ltr), The HSS tree in post ordering.
!
! LTR  (INPUT)  INTEGER, the length of HSS Tree, LTR.
! 
! LDX  (INPUT) INTEGER, the leading dimension of X
!
! NCOL (INPUT) INTEGER, the column dimension of X
!
! TRE  (INPUT) TYPE HTinfo, contains some information about HSS tree
! 
! LVL  (INPUT) INTEGER, the total levels of HSS tree
! 
! ====================
! Written by Shengguo Li, on Mar. 18th, 2013
! ==========================================
!
! ..Local Parameters..
    INTEGER   n1, nc1, i, it, lt, nj, LL, ii, ierr, &
         pdx, j, k, ni, LCG, LCF, nxi
!
! ..Array Parameters..
    INTEGER    :: lx(ltr),lg(ltr),lf(ltr)
    DOUBLE PRECISION, PARAMETER   :: ZERO = 0.0D0
    DOUBLE PRECISION, allocatable ::  G(:,:), F(:,:)
! 
! Check input variables
    IF( lvl .lt. 3) Then
       write(*,*) "HSS tree must have at least 3 levels"
       return
    END IF
!
! ..Executable Statements ..
    n1 = 0.5*(ltr+1)

!*********************************************************
!                      partition x                       *
!*********************************************************
! lx(i): column partition of x, only for leaf nodes;
! lg(i): starting column position for each G_i;
! lf(i): starting column position for each F_i;
       lt = 1
       lx = 0
       i = 1
       n1 = TRE%lentr( lvl )
       DO ni = 1, ltr
          nc1 = TRE%ch( 1,ni )
          IF( nc1 .eq. 0 ) THEN
             lx(ni) = lt
             lt = lt + m( i )
             i = i + 1
          END IF
       END DO

!*********************************************************
!                   positions of G                       *
!*********************************************************
       it = 1
       lg = 0
       DO LL = lvl, 1, -1  ! each level
          n1 = TRE%lentr( LL )
          DO i = 1, n1
             ni = TRE%ttr( LL,i )
             lg( ni ) = it
             it = it+PH%U( 2,ni )
          END DO
       END DO
       LCG = it

!*********************************************************
!                   positions of F                       *
!*********************************************************
       it = 1
       lf = 0
       DO LL = 1, lvl  ! each level
          n1 = TRE%lentr( LL )
          DO i = 1, n1
             ni = TRE%ttr( LL,i )

             IF( ni .EQ. TRE%ch( 1,tr(ni) ) ) THEN
                nj = TRE%ch( 2,tr(ni) )    ! adjacent of ni
             ELSE
                nj = TRE%ch( 1,tr(ni) )
             END IF
             lf( ni ) = it
             it = it+PH%B( 2,nj )  ! F{ni} = G{nj}*B{nj}
          END DO
       END DO
       LCF = it

       ALLOCATE( G(LDX,LCG), F(LDX,LCF), stat=ierr )
       IF(ierr /= 0 ) THEN
          write(*,*) "Allocate failed in fasthssmm_omp "
       END IF
       G = zero
       F = zero

!*********************************************************
!                    Upsweep for G                       *
!*********************************************************
       ! bottom level
       n1 = TRE%lentr( lvl )
!$OMP PARALLEL PRIVATE(i,ni)
!$OMP DO SCHEDULE(dynamic)
       DO i = 1, n1
          ni = TRE%ttr( lvl,i )
          call dgemm('N','N',LDX,PH%U(2,ni),PH%U(1,ni),1d0,X(1,lx(ni)),LDX,H( PH%ph(1,ni) ),PH%U(1,ni), & 
               0d0,G(1,lg(ni)),LDX )   ! G{i} = Xi * U(i)
       END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

       DO LL = lvl-1, 1, -1
          n1 = TRE%lentr( LL )
!$OMP PARALLEL PRIVATE(i,ni,nc1)
!$OMP DO SCHEDULE(dynamic)
          DO i = 1, n1
             ni  = TRE%ttr( LL,i )
             nc1 = TRE%ch( 1,ni )
             IF( nc1 .ne. 0 ) THEN  ! parent node
                call dgemm('N','N',LDX,PH%U(2,ni),PH%U(1,ni),1d0,G(1,lg(nc1)),LDX,H( PH%ph(1,ni) ), & 
                    PH%U(1,ni),0d0,G(1,lg(ni)),LDX )     ! G{i}=[G{i1} G{i2}]*U{i}
             ELSE ! leaf node
                call dgemm('N','N',LDX,PH%U(2,ni),PH%U(1,ni),1d0,X(1,lx(ni)),LDX,H( PH%ph(1,ni) ),& 
                     PH%U(1,ni), 0d0,G(1,lg(ni)),LDX )   ! G{i} = Xi * U(i)
             END IF
          END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       END DO

!*********************************************************
!                 Backward Traversal                     *
!*********************************************************
       ! the second top level
       LL = 1
       lt = 1
!$OMP PARALLEL PRIVATE(i,ni,nj,nc1)
!$OMP DO SCHEDULE(dynamic)
       DO i = 1, 2
          ni = TRE%ttr( LL,i )
          IF( ni .EQ. TRE%ch( 1,tr(ni) ) ) THEN
             nj = TRE%ch( 2,tr(ni) )    ! adjacent of ni
          ELSE
             nj = TRE%ch( 1,tr(ni) )
          END IF
          call dgemm('N','N',LDX,PH%B(2,nj),PH%B(1,nj),1d0,G(1,lg(nj)),LDX,H( PH%ph(3,nj) ), &
               PH%B(1,nj),0d0,F(1,lf(ni)),LDX )  ! F{i} = G{j}*B{j}

          ! update F of its children
          nc1 = TRE%ch(1,ni)      ! left child
          IF( nc1 .ne. 0 ) THEN
             call dgemm( 'N','N',LDX,PH%V(2,ni),PH%V(1,ni),1d0, F(1,lf(ni)), &
                  LDX, H( PH%ph(2,ni) ), PH%V(1,ni), 0.0d0,F(1,lf(nc1)), LDX )
          END IF
       END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

       ! Traverse downward
       DO LL = 2, lvl
          n1 = TRE%lentr( LL )
!$OMP PARALLEL PRIVATE(i,ni,nj,nc1)
!$OMP DO SCHEDULE(dynamic)
          DO i = 1, n1
             ni = TRE%ttr( LL,i )
             IF( ni .EQ. TRE%ch(1,tr(ni) ) ) THEN
                nj = TRE%ch( 2,tr(ni) )      ! adjacent of ni
             ELSE
                nj = TRE%ch( 1,tr(ni) )
             END IF
             call dgemm('N','N',LDX,PH%B(2,nj),PH%B(1,nj),1d0,G(1,lg(nj)),LDX,H( PH%ph(3,nj) ), &
                  PH%B(1,nj),1d0,F(1,lf(ni)),LDX )  ! F{i} = G{j}*B{j}

             ! update F of its children
             nc1 = TRE%ch(1,ni)      ! left child
             IF( nc1 .ne. 0 ) THEN
                call dgemm( 'N','N',LDX,PH%V(2,ni),PH%V(1,ni),1d0, F(1,lf(ni)), &
                     LDX, H( PH%ph(2,ni) ), PH%V(1,ni), 0.0d0,F(1,lf(nc1)), LDX )
             END IF
          END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       END DO ! (LL)

!*********************************************************
!                   Compute X                            *
!*********************************************************
       DO ni = 1, ltr-1
          nc1 = TRE%ch( 1,ni )
          IF( nc1 .eq. 0 ) THEN
             call dgemm('N','N',LDX,PH%V(2,ni),PH%V(1,ni),1d0,F(1,lf(ni)),LDX,H( PH%ph(2,ni) ),& 
                  PH%V(1,ni),0d0,G(1,1),LDX )
             
             call dgemm('N','N',LDX,PH%D(ni),PH%D(ni),1d0,X(1,lx(ni)),LDX,DD( PH%pd(1,ni) ), & 
                  PH%D(ni), 1d0, G(1,1), LDX )
             call dlacpy('A', LDX,PH%D(ni),G(1,1),LDX,X(1,lx(ni) ), LDX )
          END IF
       END DO

       DEALLOCATE( G, F )

     END SUBROUTINE fastHssmmL_omp

!!!!!!
  SUBROUTINE MDLASD32( K,D,Z,DSIGMA,DIFL,DIFR,ALPHA_L,ALPHA_R,INFO )
!
    implicit none

!     .. Scalar Arguments ..
      INTEGER            INFO, K
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D(*),Z(*),DSIGMA(*),DIFL(*),DIFR(*),ALPHA_L(*),ALPHA_R(*)
!     ..
!
!  Purpose
!  =======
!
!  MDLASD3 finds the square roots of the roots of the secular equation,
!  as defined by the values in DSIGMA and Z. It makes the appropriate
!  calls to DLASD4.
!
!  DLASD8 is called from DLASD6.
!
!  Arguments
!  =========
!
!  K       (input) INTEGER
!          The number of terms in the rational function to be solved
!          by DLASD4.  K >= 1.
!
!  D       (output) DOUBLE PRECISION array, dimension ( K )
!          On output, D contains the updated singular values.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension ( K )
!          On entry, the first K elements of this array contain the
!          components of the deflation-adjusted updating row vector.
!          On exit, Z is updated.
!
!  DSIGMA  (input/output) DOUBLE PRECISION array, dimension ( K )
!          On entry, the first K elements of this array contain the old
!          roots of the deflated updating problem.  These are the poles
!          of the secular equation.
!          On exit, the elements of DSIGMA may be very slightly altered
!          in value.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension at least 3 * K
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J
      DOUBLE PRECISION   RHO, TEMP,DIFLJ,DIFRJ,DJ,DSIGJ,DSIGJP
      DOUBLE PRECISION, allocatable :: WORK(:)
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASCL, DLASD4, DLASET, XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      Allocate( Work(3*K) )
!      
      IF( K.LT.1 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASD8', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
         DIFL( 1 ) = D( 1 )
         RETURN
      END IF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
      END DO
!
!     Book keeping.
!
      IWK1 = 1
      IWK2 = IWK1 + K
      IWK3 = IWK2 + K
      IWK2I = IWK2 - 1
      IWK3I = IWK3 - 1
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Initialize WORK(IWK3).
!
      CALL DLASET( 'A', K, 1, ONE, ONE, WORK( IWK3 ), K )
!
!     Compute the updated singular values, the arrays DIFL, DIFR,
!     and the updated Z.
!
      DO J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ), &
                     WORK( IWK2 ), INFO )
!
!        If the root finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'DLASD4', -INFO )
            RETURN
         END IF
         WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
         DIFL( J ) = -WORK( J )  
         DIFR( J ) = -WORK( J+1 )
         DO I = 1, J - 1
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )*      &
                             WORK( IWK2I+I ) / ( DSIGMA( I )- &
                             DSIGMA( J ) ) / ( DSIGMA( I )+   &
                             DSIGMA( J ) )
         END DO
         DO I = J + 1, K
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )*       &
                             WORK( IWK2I+I ) / ( DSIGMA( I )- &
                             DSIGMA( J ) ) / ( DSIGMA( I )+   &
                             DSIGMA( J ) )
         END DO
      END DO
!
!     Compute updated Z.
!
      DO I = 1, K
         Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
      END DO
!
!     Compute the left and right scalars
!
      DO J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         IF( J.LT.K ) THEN
            DIFRJ = -DIFR( J )
            DSIGJP = -DSIGMA( J+1 )
         END IF
         WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
         DO I = 1, J - 1
            WORK( I ) = Z( I ) / ( ( DSIGMA( I )+ DSIGJ )-DIFLJ ) / ( DSIGMA( I )+DJ )
         END DO
         DO I = J + 1, K
            WORK( I ) = Z( I ) / ( ( DSIGMA( I )+DSIGJP )+DIFRJ ) / ( DSIGMA( I )+DJ )
         END DO
         TEMP = DNRM2( K, WORK, 1 )
         ALPHA_R( J ) = one / TEMP
      END DO

      WORK(K+1)= -1.0D0
      DO J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         IF( J.LT.K ) THEN
            DIFRJ = -DIFR( J )
            DSIGJP = -DSIGMA( J+1 )
         END IF
         WORK( J ) = -DSIGMA(J)*Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
         DO I = 1, J - 1
            WORK( I ) = DSIGMA(I)*Z( I ) / ( ( DSIGMA( I )+ DSIGJ )-DIFLJ ) / ( DSIGMA( I )+DJ )
         END DO
         DO I = J + 1, K
            WORK( I ) = DSIGMA(I)*Z( I ) / ( ( DSIGMA( I )+DSIGJP )+DIFRJ ) / ( DSIGMA( I )+DJ )
         END DO
         TEMP = DNRM2( K+1, WORK, 1 )
         ALPHA_L( J ) = one / TEMP
      END DO

      deallocate( WORK )

!     End of MDLASD32
!
    END SUBROUTINE MDLASD32

!
! 
   end module ConstructHssd

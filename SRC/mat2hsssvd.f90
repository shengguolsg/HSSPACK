module mat2hsssvd
! This procedure is random svd + column majority. 

  implicit none

contains

!!!!!!!!!!! 
  subroutine mat2hsscholsvdr(A,lda,tr,ltr,m,lm,npv,RH,R,tol,ncol,T,RT,pp)
    use aux_hss
!
! ..
! .. Scalar parameter ..
    INTEGER, INTENT(IN)          :: lda,ltr,lm,npv,ncol,pp
    DOUBLE PRECISION, INTENT(IN) :: tol
! ..
! .. Arrary parameter ..
    INTEGER, INTENT(INOUT)  :: TR(*), M(*)
    DOUBLE PRECISION, INTENT(INOUT)  :: A(lda,*),R(*),T(*) ! T: workspace
    type(hss)     :: RH
    type(ddcell)  :: RT
!
! Purpose
! ==========
!    
! This routine constructs a Cholesky factorization for an HSS positive definite matrix.
! It uses randomized SVD to compute the low rank approximations of off-diagonal blocks. 
! This routine is for our paper:
!    "New efficient and robust HSS Cholesky factorization of Symmetric positive 
!     definite matrices, SIAM J. Matrix Anal. Appl., 33, 2012, pp 886--904."
! 
! This routine is written when I was not familiar with Fortran language. Therefore it
! may not be efficient enough, and readable. At least it works. We can make it better
! in later time. 
!
! .. Parameters .. 
! ===============
! 
! A   (inout)  DOUBLE PRECISION Array, DIMENSION(LDA, *)
!      The original matrix to be approximated. It is also used as a workspace. 
!      Its entries will be destroyed. 
! 
! LDA (in) INTEGER, the leading dimension of A
!
! TR  (inout) INTEGER Array, DIMENSION( ltr )
!     It stores the postordering HSS tree. 
! 
! LTR (in) INTEGER, the length of HSS tree
!
! M   (inout)  INTEGER Array, DIMENSION( LM )
!     The number of rows of each leaf nodes.
! 
! lm  (inout) INTEGER, number of lead nodes
!     Seems never used.
! 
! npv (in)  INTEGER, number of nodes, equals to LTR
!
! RH  (inout) DOUBLE PRECISION Array, DIMENSION( * )
!     hsschol type, stores all the generator of Cholesky factor
! 
! R  (inout) DOUBLE PRECISION Array, DIMENSION( * )
!    Used as workspace.
! 
! tol (in) DOUBLE PRECISION, tolerance for low-rank approximation
!     
! ncol (in) INTEGER, maximum rank allowed.
!
! T  (inout) DOUBLE PRECISION Array, DIMENSION( * )
!    Used as workspace.
! 
! RT (inout) DOUBLE PRECISION Array, DIMENSION( * )
!    ddcell type. 
! 
! pp (in) INTEGER, parameter for oversampling.
!     
! ========
! Written by Shengguo Li, Sep. 13th, 2013 
! =======================================
!
! .. local scalars ..
    integer  :: rk,ch1,pnum,err,pnr,lwork,pnd,ns,nt,i,info,ir,&
                it,j,lt,mi,ni,nn
! ..
! .. local arrays ..
    integer, allocatable, dimension(:,:) :: ch, l
    integer, allocatable, dimension(:) :: ls, ws
    double precision, allocatable, dimension(:) :: ww,DD,svs,work,Omega
       
       lwork=lda*lda
       pnr=1
       pnd=1
       allocate(ch(2,ltr),l(ltr,2),ls(ltr), ws(ltr),ww(lda*lda))
       allocate(DD(lda*lda), svs(lda),work(lda*lda), Omega(lda*pp), stat=err)

       nn = min(npv,ltr-1)
       call child(tr, ch, ltr)
       ! Omega should be a random matrix, pp*lda. 
       
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

       ns = 1 ! (# of blocks in stack)+1
       ls(ns) = 1
       ws(ns) = l(1,2)+1
       call random_number(Omega)
       
       do i = 1, nn
          ! factorize (leaf) diag & get off-diag block

          if (ch(1,i) == 0 .and. ch(2,i) == 0) then 
             mi = l(i,2)-l(i,1)+1
             ni = lda-l(i,2) 		!t^r, the right column block's length

             if (ns == 1) then
                call dlacpy('A', ni, mi, A(l(i,2)+1,l(i,1)), lda, R(pnr), ni) ! Construct T, the HSS column block
                nt = ni
             else
                nt = ls(ns)-1+ni        !the whole HSS block's length, which is store in R, starting from 'pnr'
                call transpose1d(mi, ls(ns)-1, A(ws(ns-1),1), lda, R(pnr), nt)
                if( l(i,2)+1 < lda ) then 
                   call dlacpy('A',ni,mi,A(l(i,2)+1,l(i,1)),lda,R(pnr+ls(ns)-1),nt) 
                end if
             endif
             
             if (mi >= 20 .AND. nt >= 20) then
		call dlacpy('A', mi, mi, A(l(i,1),l(i,1)), lda, T, mi)  !the front of T store A_{ii}
		call DPOTRF('L', mi, T, mi, info)        

                call dgemm('N','N',pp,mi,nt,1d0,Omega,pp,R(pnr),nt,0d0,ww,pp) !ww=Omega*T, starting from 1.
		call DTRSM('R', 'L', 'T', 'N', pp, mi, 1d0, T, mi, ww, pp) !ww (pp-by-mi) stores (Omega*T)*(L^{-T}) 
                call dgesvd('N','A',pp,mi,ww,pp,svs,ww,pp,ww(1+mi*pp),mi,work,lwork,info)  !svd of ww. 
                rk=min(ncol,count(svs(1:mi)>tol))

		ir =mi-rk
		RH%pp(i)=rk !remainig blocks
		if (ir > 0) then	   
                   call DTRSM('R','L','T','N',nt,mi,1d0,T,mi,R(pnr),nt) ! T=TL^{-T}, still stores in R(pnr)
                   call transpose1ds(ww(1+mi*pp),mi,mi)                 ! V**T --> V
                   call dgemm('N','N',nt,rk,mi,1d0,R(pnr),nt,ww(1+mi*pp),mi,0d0,T(1+mi*mi),nt) 

                   call dlacpy('A',mi,rk,ww(1+mi*pp),mi,ww(1+mi*pp+mi*mi),mi)
                   call dlacpy('A',mi,mi,ww(1+mi*pp+mi*rk),mi,R(pnr),mi) ! R(pnr) stores Q
		   call dtrmm('L','L','N','N',mi,mi,1d0,T,mi,R(pnr),mi)  !XQ=UL, R output
		   call dgeqlf( mi, mi, R(pnr), mi, svs, work, lwork, INFO)  !D and Q are togother, plus Tau=svs

		   !update T=L(r+1:end,r+1:end)*T		   
		   call dtrmm('R','L','T','N',nt,rk,1d0,R(pnr+ir*(mi+1)),mi,T(mi*mi+1),nt) 
		   call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .true.) 
                else
                   call dlacpy('A',mi,mi,A(l(i,1),l(i,1)),lda,DD(pnd),mi)  !DD{i} 
		   call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .false.)
		endif                
             else
                rk = mi 
		ir=0
                call dlacpy('A',mi,mi,A(l(i,1),l(i,1)),lda,DD(pnd),mi)  !DD{i} 
                call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .false.)
             endif
             
               ls(ns+1) = ls(ns)+rk
               ws(ns) = l(i,2)+1
               ns = ns+1
               if(l(i,2)+1 < lda ) then 
                  if (ir >0) then
                     call dlacpy('A', ni, rk, T(mi*mi+nt-ni+1), nt, A(l(i,2)+1,ls(ns-1)), lda)
                  else   
                     call dlacpy('A', ni, rk, A(l(i,2)+1,l(i,1)), lda, A(l(i,2)+1,ls(ns-1)), lda)		
                  endif
               end if

          else
             mi = ls(ns)-ls(ns-2)
             ni = lda-l(i,2)
	     ! Merge, use the front part of T to store Dt
	     ch1=ch(1,i)
             call dlacpy('A', RH%pp(ch1), RT%pmi(2,ch1), DD(RT%pmi(1,ch1)), RT%pmi(2,ch1),T,mi) !DD{1} 	     
	     pnum=1+RH%pp(ch1)*mi
             call dlacpy('A',RH%pp(ch(2,i)),RH%pp(ch1),A(ws(ns-2),ls(ns-2)),lda,T(1+RH%pp(ch1)),mi) ! B{i}
	     call transpose1d(RH%pp(ch(2,i)),RH%pp(ch1),A(ws(ns-2),ls(ns-2)), lda, T(pnum), mi) !B{i}'
             ch1=ch(2,i)
             pnum=RH%pp(ch1)
	     call dlacpy('A', pnum, pnum, DD(RT%pmi(1,ch1)), pnum,T((mi+1)*RH%pp(ch(1,i))+1),mi) !DD{2}
	     	     	     
             if (ns == 3) then
                call dlacpy('A', ni, mi, A(l(i,2)+1,ls(ns-2)), lda, R(pnr), ni)
                nt = ni
             else
		nt = ls(ns-2)-1+lda-l(i,2)
                call transpose1d(mi,ls(ns-2)-1, A(ws(ns-3),1), lda, R(pnr), nt)
                if(l(i,2)+1 < lda ) then 
                   call dlacpy('A', ni, mi, A(l(i,2)+1,ls(ns-2)), lda, R(pnr+ls(ns-2)-1), nt)
                end if
             endif
	
  	     if (mi >= 20 .AND. nt >= 20) then
		call DPOTRF('L', mi, T, mi, info) ! lower-triangular is chol
                ! pp = min(mi+10, nt)
                call dgemm('N','N',pp,mi,nt,1d0,Omega,pp,R(pnr),nt,0d0,ww,pp) !ww=Omega*T, starting from 1.
		call DTRSM('R', 'L', 'T', 'N', pp, mi, 1d0, T, mi, ww, pp) !ww (pp-by-mi) stores (Omega*T)*(L^{-T})
		call dgesvd('N','A',pp,mi,ww,pp,svs,ww,pp,ww(1+mi*pp),mi,work,lwork,info)  !svd of ww. 
        	rk=min(ncol,count(svs(1:mi)>tol))

		ir = mi-rk
		RH%pp(i)=rk
		if (ir > 0) then
                  call DTRSM('R','L','T','N',nt,mi,1d0,T,mi,R(pnr),nt) ! T=TL^{-T}, still stores in R(pnr)
                  call transpose1ds(ww(1+mi*pp),mi,mi)                 ! V**T --> V
                  call dgemm('N','N',nt,rk,mi,1d0,R(pnr),nt,ww(1+mi*pp),mi,0d0,T(1+mi*mi),nt) 

                  !call circshift2(ww(1+mi*pp),mi,mi,rk,2) ! Q=[V_2 V_1]
                  call dlacpy('A',mi,rk,ww(1+mi*pp),mi,ww(1+mi*pp+mi*mi),mi)
                  call dlacpy('A',mi,mi,ww(1+mi*pp+mi*rk),mi,R(pnr),mi) ! R(pnr) stores Q
		  call dtrmm('L','L','N','N',mi,mi,1d0,T,mi,R(pnr),mi)  !XQ=UL, R output
		  call dgeqlf( mi, mi, R(pnr), mi, svs, work, lwork, INFO)  

		  !update T=T*L(r+1:end,r+1:end)^T
		  call dtrmm('R','L','T','N',nt,rk,1d0,R(pnr+ir*(mi+1)),mi,T(mi*mi+1),nt) 
		  call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .true.) 
	 	else
		    call dtrmms('L','L','T','N',mi,1d0,T,mi) 
		    call dlacpy('A',mi,mi,T,mi,DD(pnd),mi) 
                    call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .false.)
               	endif                
             else
                call dlacpy('A',mi,mi,T,mi,DD(pnd),mi)  !copy Dt to DD 
                rk = mi
                ir=0
                call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .false.)                
             endif

        	ns = ns-1 ! shrink ls
                ls(ns) = ls(ns-1)+rk
                if (ir >0 .and. (l(i,2)+1 < lda) ) then
                   call dlacpy('A', ni, rk, T(mi*mi+nt-ni+1), nt, A(l(i,2)+1,ls(ns-1)), lda)
                endif
         endif

          do j = 1, ns-3
	     if (ir > 0) then
                A(ws(j)+ir:ws(ns-2)-1+ir,ls(j):ls(j+1)-1) = A(ws(j):ws(ns-2)-1,ls(j):ls(j+1)-1)
                ws(j) = ws(j)+ir
	     endif
          enddo
          ws(ns-1) = l(i,2)+1
          if (ns > 2 .and. ir >0) then
             call transpose1d(ls(ns-1)-1, rk, T(mi*mi+1), nt, A(ws(ns-2)+ir,1), lda)
             ws(ns-2) = ws(ns-2)+ir
          endif

       enddo
	    i=ltr
            mi = ls(ns)-ls(ns-2)         
	    ! Merge, use the front part of T to store Dt
      	    ch1=ch(1,i)
            call dlacpy('A', RT%pmi(2,ch1), RT%pmi(2,ch1), DD(RT%pmi(1,ch1)), RT%pmi(2,ch1),T,mi) !DD{1} 	     
            pnum=1+RH%pp(ch1)*mi
            call dlacpy('A',RH%pp(ch(2,i)),RH%pp(ch1), A(ws(ns-2),ls(ns-2)), lda,T(1+RH%pp(ch1)),mi) ! B{i}
            call transpose1d(RH%pp(ch(2,i)),RH%pp(ch1),A(ws(ns-2),ls(ns-2)), lda, T(pnum), mi) !B{i}'
	    ch1=ch(2,i)
            pnum=RH%pp(ch1)
	    call dlacpy('A', pnum, pnum, DD(RT%pmi(1,ch1)), pnum,T((mi+1)*RH%pp(ch(1,i))+1),mi) !DD{2}

            call DPOTRF('L',mi,T,mi,info)
            call dlacpy('A',mi,mi,T,mi,R(pnr),mi)
            RH%p(1,i)=pnr
            RH%p(2,i)=pnr+mi*mi
            RH%dmi(i)=mi

            deallocate( Omega,svs,DD,ww,ws,ls,l,ch )
!
       end subroutine mat2hsscholsvdr

!!!!!!!
    subroutine mat2hsscholsvd(A, lda, tr, ltr, m, lm, npv, RH, R, tol, ncol, T, RT)
       use aux_hss
! ..
! .. Scalar parameter ..
    INTEGER, INTENT(IN)          :: lda,ltr,lm,npv,ncol
    DOUBLE PRECISION, INTENT(IN) :: tol
! ..
! .. Arrary parameter ..
    INTEGER, INTENT(INOUT)  :: TR(*), M(*)
    DOUBLE PRECISION, INTENT(INOUT)  :: A(lda,*),R(*),T(*) ! T: workspace
    type(hss)     :: RH
    type(ddcell)  :: RT
!
! Purpose
! ==========
!    
! This routine constructs a Cholesky factorization for an HSS positive definite matrix.
! It uses truncated SVD to compute the low rank approximations of off-diagonal blocks. 
! This routine is for our paper:
!    "New efficient and robust HSS Cholesky factorization of Symmetric positive 
!     definite matrices, SIAM J. Matrix Anal. Appl., 33, 2012, pp 886--904."
! 
! This routine is written when I was not familiar with Fortran language. Therefore it
! may not be efficient enough, and readable. At least it works. We can make it better
! in later time. 
!
! .. Parameters .. 
! ===============
! 
! A   (inout)  DOUBLE PRECISION Array, DIMENSION(LDA, *)
!      The original matrix to be approximated. It is also used as a workspace. 
!      Its entries will be destroyed. 
! 
! LDA (in) INTEGER, the leading dimension of A
!
! TR  (inout) INTEGER Array, DIMENSION( ltr )
!     It stores the postordering HSS tree. 
! 
! LTR (in) INTEGER, the length of HSS tree
!
! M   (inout)  INTEGER Array, DIMENSION( LM )
!     The number of rows of each leaf nodes.
! 
! lm  (inout) INTEGER, number of lead nodes
!     Seems never used.
! 
! npv (in)  INTEGER, number of nodes, equals to LTR
!
! RH  (inout) DOUBLE PRECISION Array, DIMENSION( * )
!     hsschol type, stores all the generator of Cholesky factor
! 
! R  (inout) DOUBLE PRECISION Array, DIMENSION( * )
!    Used as workspace.
! 
! tol (in) DOUBLE PRECISION, tolerance for low-rank approximation
!     
! ncol (in) INTEGER, maximum rank allowed.
!
! T  (inout) DOUBLE PRECISION Array, DIMENSION( * )
!    Used as workspace.
! 
! RT (inout) DOUBLE PRECISION Array, DIMENSION( * )
!    ddcell type. 
! 
! ========
! Modified by Shengguo Li, Sep. 13th, 2013 
! =======================================
!
! .. local scalars ..
    integer       :: rk,ch1,pnum,pnr,lwork,pnd,ns,nt,i,info,ir,&
                     lt,mi,ni,nn,it,j
    double precision :: err
! ..
! .. local arrays ..
    integer, allocatable, dimension(:,:) :: ch, l
    integer, allocatable, dimension(:) :: ls, ws
    double precision, allocatable, dimension(:) :: ww,Tau,DD,svs,work
       
       lwork=lda*lda
       pnr=1
       pnd=1
       allocate (ch(2,ltr), l(ltr,2), ls(ltr), ws(ltr), ww(lda*lda), Tau(lda), DD(lda*lda), svs(3*lda),work(lda*lda))
       nn = min(npv,ltr-1)
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

       ns = 1 ! (# of blocks in stack)+1
       ls(ns) = 1
       ws(ns) = l(1,2)+1

       do i = 1, nn
          ! factorize (leaf) diag & get off-diag block
          if (ch(1,i) == 0 .and. ch(2,i) == 0) then 
             mi = l(i,2)-l(i,1)+1
             ni = lda-l(i,2) 		!remaining column block's length
              
             if (ns == 1) then
                call dlacpy('A', ni, mi, A(l(i,2)+1,l(i,1)), lda, R(pnr), ni) ! Construct T
                nt = ni
             else
                nt = ls(ns)-1+ni
                call transpose1d(mi, ls(ns)-1, A(ws(ns-1),1), lda, R(pnr), nt)
                if(l(i,2)+1<lda) then
                   call dlacpy('A',ni,mi,A(l(i,2)+1,l(i,1)),lda,R(pnr+ls(ns)-1),nt) 
                end if
             endif
             
             if (mi >= 20 .AND. nt >= 20) then
		call dlacpy('A', mi, mi, A(l(i,1),l(i,1)), lda, T, mi)
		call DPOTRF('L', mi, T, mi, info)         !lower-triangular is chol

                call dgeqrf(nt,mi,R(pnr),nt,Tau,work,lwork,info) !T=QR
                ww(1:mi*mi)=0
                call dlacpy('U',mi,mi,R(pnr),nt,ww,mi)
		call DTRSM('R', 'L', 'T', 'N', mi, mi, 1d0, T, mi, ww, mi) !ww stores R*(L^{-T}) ok
                call dgesvd('O','S',mi,mi,ww,mi,svs,ww,mi,ww(1+mi*mi),mi,work,lwork,info)  !T(1+mi*mi) ok. sign different
                rk=min(ncol,count(svs(1:mi)>tol))

		ir =mi-rk
		RH%pp(i)=rk
		if (ir > 0) then	   
                   T(1+mi*mi:mi*mi+nt*rk)=0
                   call Usigma(ww,mi,svs,rk,nt,T(1+mi*mi))
                   call dormqr('L','N',nt,rk,min(nt,mi),R(pnr),nt,Tau,T(mi*mi+1),nt,work,lwork,info) ! T is already T*Q_1^T

                   call circshift2(ww(1+mi*mi),mi,mi,rk,1) ! V**T
                   call transpose1d(mi,mi,ww(1+mi*mi),mi,R(pnr),mi) ! R(pnr) stores Q
		   call dtrmm('L','L','N','N',mi,mi,1d0,T,mi,R(pnr),mi)  !XQ=UL, R output
		   call dgeqlf( mi, mi, R(pnr), mi, TAU, work, lwork, INFO)  !D and Q are togother, plus Tau

		   !update T=L(r+1:end,r+1:end)*T		   
		   call dtrmm('R','L','T','N',nt,rk,1d0,R(pnr+ir*(mi+1)),mi,T(mi*mi+1),nt) 
		   call hssexpsvd(RH, R, pnr, i, Tau, ir, rk, mi, DD, RT, pnd, ww, .true.) 
                else
                   call dlacpy('A',mi,mi,A(l(i,1),l(i,1)),lda,DD(pnd),mi)  !DD{i} 
                   call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .false.)
		endif                
             else
                call dlacpy('A',mi,mi,A(l(i,1),l(i,1)),lda,DD(pnd),mi)  !DD{i} 
                rk = mi 
		ir=0
                call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .false.)
             endif
             
               ls(ns+1) = ls(ns)+rk
               ws(ns) = l(i,2)+1
               ns = ns+1
               if(l(i,2)+1 < lda ) then
                  if (ir >0) then
                     call dlacpy('A', ni, rk, T(mi*mi+nt-ni+1), nt, A(l(i,2)+1,ls(ns-1)), lda)
                  else   
                     call dlacpy('A', ni, rk, A(l(i,2)+1,l(i,1)), lda, A(l(i,2)+1,ls(ns-1)), lda)		
                  endif
               end if
          else
             mi = ls(ns)-ls(ns-2)
             ni = lda-l(i,2)
	     !Merge, use the front part of T to store Dt
	     ch1=ch(1,i)
             call dlacpy('A', RH%pp(ch1), RT%pmi(2,ch1), DD(RT%pmi(1,ch1)), RT%pmi(2,ch1),T,mi) !DD{1} 	     
	     pnum=1+RH%pp(ch1)*mi
             call dlacpy('A',RH%pp(ch(2,i)),RH%pp(ch1),A(ws(ns-2),ls(ns-2)),lda,T(1+RH%pp(ch1)),mi) ! B{i}
	     call transpose1d(RH%pp(ch(2,i)),RH%pp(ch1),A(ws(ns-2),ls(ns-2)), lda, T(pnum), mi) !B{i}'
             ch1=ch(2,i)
             pnum=RH%pp(ch1)
	     call dlacpy('A', pnum, pnum, DD(RT%pmi(1,ch1)), pnum,T((mi+1)*RH%pp(ch(1,i))+1),mi) !DD{2}
	     	     	     
             if (ns == 3) then
                call dlacpy('A', ni, mi, A(l(i,2)+1,ls(ns-2)), lda, R(pnr), ni)
                nt = ni
             else
		nt = ls(ns-2)-1+lda-l(i,2)
                call transpose1d(mi,ls(ns-2)-1, A(ws(ns-3),1), lda, R(pnr), nt)
                if( l(i,2)+1 < lda ) then
                   call dlacpy('A', ni, mi, A(l(i,2)+1,ls(ns-2)), lda, R(pnr+ls(ns-2)-1), nt)                
                end if
             endif
	
  	     if (mi >= 20 .AND. nt >= 20) then
		call DPOTRF('L', mi, T, mi, info) ! lower-triangular is chol
          	call dgeqrf(nt,mi,R(pnr),nt,Tau,work,lwork,info) !T=QR
        	ww(1:mi*mi)=0
        	call dlacpy('U',mi,mi,R(pnr),nt,ww,mi)
		call DTRSM('R', 'L', 'T', 'N', mi, mi, 1d0, T, mi, ww, mi) !ww stores R*(L^{-T})
        	call dgesvd('O','S',mi,mi,ww,mi,svs,ww,mi,ww(1+mi*mi),mi,work,lwork,info) 
        	rk=min(ncol,count(svs(1:mi)>tol))

		ir = mi-rk
		RH%pp(i)=rk
		if (ir > 0) then
                  T(1+mi*mi:mi*mi+nt*rk)=0
                  call Usigma(ww,mi,svs,rk,nt,T(1+mi*mi))
                  call dormqr('L','N',nt,rk,min(mi,nt),R(pnr),nt,Tau,T(mi*mi+1),nt,work,lwork,info) ! T is already T*Q_1^T

                  call circshift2(ww(1+mi*mi),mi,mi,rk,1) ! V**T
                  call transpose1d(mi,mi,ww(1+mi*mi),mi,R(pnr),mi)
		  call dtrmm('L','L','N','N',mi,mi,1d0,T,mi,R(pnr),mi)  !XQ=UL, R output
		  call dgeqlf( mi, mi, R(pnr), mi, TAU, work, lwork, INFO)  

		   !update T=T*L(r+1:end,r+1:end)^T
		   call dtrmm('R','L','T','N',nt,rk,1d0,R(pnr+ir*(mi+1)),mi,T(mi*mi+1),nt) 
		   call hssexpsvd(RH, R, pnr, i, Tau, ir, rk, mi, DD, RT, pnd, ww, .true.) 
	 	else
		    call dtrmms('L','L','T','N',mi,1d0,T,mi) 
		    call dlacpy('A',mi,mi,T,mi,DD(pnd),mi) !Dt to DD
                    call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .false.)
		endif                
             else
                call dlacpy('A',mi,mi,T,mi,DD(pnd),mi)  !copy Dt to DD 
                rk = mi
                ir=0		 		
                call hssexpsvd(RH, R, pnr, i, svs, ir, rk, mi, DD, RT, pnd, ww, .false.)
             endif
        	ns = ns-1 ! shrink ls
	       ls(ns) = ls(ns-1)+rk
             
             if (ir >0 .and. (l(i,2)+1 < lda ) ) then
                call dlacpy('A', ni, rk, T(mi*mi+nt-ni+1), nt, A(l(i,2)+1,ls(ns-1)), lda)
             endif

         endif

          do j = 1, ns-3
	     if (ir > 0) then
             	A(ws(j)+ir:ws(ns-2)-1+ir,ls(j):ls(j+1)-1) = A(ws(j):ws(ns-2)-1,ls(j):ls(j+1)-1)
             	ws(j) = ws(j)+ir
	     endif
          enddo
          ws(ns-1) = l(i,2)+1
          if (ns > 2 .and. ir >0) then
             call transpose1d(ls(ns-1)-1, rk, T(mi*mi+1), nt, A(ws(ns-2)+ir,1), lda)
             ws(ns-2) = ws(ns-2)+ir
          endif

       enddo
	    i=ltr
            mi = ls(ns)-ls(ns-2)         
	    !Merge, use the front part of T to store Dt
      	    ch1=ch(1,i)
            call dlacpy('A', RT%pmi(2,ch1), RT%pmi(2,ch1), DD(RT%pmi(1,ch1)), RT%pmi(2,ch1),T,mi) !DD{1} 	     
            pnum=1+RH%pp(ch1)*mi
            call dlacpy('A',RH%pp(ch(2,i)),RH%pp(ch1), A(ws(ns-2),ls(ns-2)), lda,T(1+RH%pp(ch1)),mi) ! B{i}
            call transpose1d(RH%pp(ch(2,i)),RH%pp(ch1),A(ws(ns-2),ls(ns-2)), lda, T(pnum), mi) !B{i}'
	    ch1=ch(2,i)
            pnum=RH%pp(ch1)
	     call dlacpy('A', pnum, pnum, DD(RT%pmi(1,ch1)), pnum,T((mi+1)*RH%pp(ch(1,i))+1),mi) !DD{2}

             call DPOTRF('L',mi,T,mi,info)
             call dlacpy('A',mi,mi,T,mi,R(pnr),mi)
	     RH%p(1,i)=pnr
             RH%p(2,i)=pnr+mi*mi
             RH%dmi(i)=mi

       deallocate (ch, l, ls, ws, ww, Tau, DD)
       end subroutine mat2hsscholsvd

!!!!!!
    subroutine solvecholsvd(RH, R, m, tr, ltr, b, lda, Work, lwork)
       use aux_hss
!
! .. Scalar parameters ..
       INTEGER, INTENT(IN)  :: ltr, lwork, lda
! ..
! .. Array parameters ..
       INTEGER, INTENT(IN) :: M(*), TR(*)
       DOUBLE PRECISION, INTENT(IN) :: R(*)
       DOUBLE PRECISION, INTENT(INOUT) :: WORK(*),b(*)
       type(hss)     :: RH
!
! ===========
! Modified by Shengguo, Li, on Sept, 13th, 2013
! =============================================
!
! .. local scalars ..
       integer    :: mi,rk,nn,i,lt,it,info,j
!
! .. local arrays ..
       integer, allocatable, dimension(:,:) :: ch,l
       double precision, allocatable, dimension(:) :: tb
       double precision, allocatable, dimension(:,:) :: Y
!
       nn=maxval(RH%dmi)
       allocate (ch(2,ltr), l(ltr,2), tb(lda), Y(nn,ltr)) 
       nn=ltr-1
       call child(tr, ch, ltr)

       l(1,1:2) = (/1,m(1)/)
       lt = 1
       it = 1
       do i = 1, ltr
          if (ch(1,i) == 0 .and. ch(2,i) == 0) then
             l(i,1:2) = (/lt,lt+m(it)-1/)
             Y(1:l(i,2)-l(i,1)+1,i)=b(l(i,1):l(i,2))
             lt = l(i,2)+1
             it = it+1
          else
             l(i,1:2) = (/l(ch(1,i),1), l(ch(2,i),2)/)
          endif
       enddo

	!Forward substitution
	nn=ltr-1
	do i=1, nn
	   if (ch(1,i) == 0 .and. ch(2,i) == 0) then
		if (RH%d(1,i) .NE. 0) then
		    mi=RH%dmi(i)
                    rk=mi-RH%pp(i)
                    call dormql('L','T',mi,1,mi,R(RH%p(1,i)),mi,R(RH%p(2,i)),Y(1:mi,i),mi,work,lwork,info)

                    if (RH%pp(i) /= mi) then
                       call dtrsv('L','N','N',rk,R(RH%p(1,i)),mi,Y(1,i),1)
                       call dgemv('N',RH%pp(i),rk,-1d0,R(RH%p(1,i)+rk),mi,Y(1:rk,i),1,1d0,Y(1+rk:mi,i),1)
                    endif
                 endif
	   else
		mi=RH%dmi(i)
                rk=mi-RH%pp(i)
                Y(1:mi,i)=[Y(RH%dmi(ch(1,i))-RH%pp(ch(1,i))+1:RH%dmi(ch(1,i)),ch(1,i)), & 
                    &Y(RH%dmi(ch(2,i))-RH%pp(ch(2,i))+1:RH%dmi(ch(2,i)),ch(2,i))]
               if (RH%d(1,i) .NE. 0) then
                    call dormql('L','T',mi,1,mi,R(RH%p(1,i)),mi,R(RH%p(2,i)),Y(1:mi,i),mi,work,lwork,info)
                    
                    call dtrsv('L','N','N',rk,R(RH%p(1,i)),mi,Y(1:mi,i),1)
                    call dgemv('N',RH%pp(i),rk,-1d0,R(RH%p(1,i)+rk),mi,Y(1:rk,i),1,1d0,Y(1+rk:RH%dmi(i),i),1)
               endif
	   endif
	enddo

 !the root
 i=ltr
 mi=RH%dmi(i)
 tb(1:mi)=[Y(RH%dmi(ch(1,i))-RH%pp(ch(1,i))+1:RH%dmi(ch(1,i)),ch(1,i)), &
                    & Y(RH%dmi(ch(2,i))-RH%pp(ch(2,i))+1:RH%dmi(ch(2,i)),ch(2,i))]
  call dtrsv('L','N','N',mi,R(RH%p(1,i)),mi,tb,1)
  call dtrsv('L','T','N',mi,R(RH%p(1,i)),mi,tb,1)

  Y(RH%dmi(ch(1,i))-RH%pp(ch(1,i))+1:RH%dmi(ch(1,i)),ch(1,i))=tb(1:RH%pp(ch(1,i)))
  Y(RH%dmi(ch(2,i))-RH%pp(ch(2,i))+1:RH%dmi(ch(2,i)),ch(2,i))=tb(RH%pp(ch(1,i))+1:mi)

  !backward substitution
  do i=nn,1,-1
     if (ch(1,i) /=0 .and. ch(2,i) /=0) then
        mi=RH%pp(ch(1,i))+RH%pp(ch(2,i))
        rk=mi-RH%pp(i)
        if (rk > 0) then
           call dgemv('T',RH%pp(i), rk, -1d0,R(RH%p(1,i)+rk), mi,Y(1+rk:mi,i),1,1d0,Y(1:rk,i),1)
           call dtrsv('L','T','N',rk,R(RH%p(1,i)),mi,Y(1:rk,i),1)
           tb(1:mi)=Y(1:mi,i)

           call dormql('L','N',mi,1,mi,R(RH%p(1,i)),mi,R(RH%p(2,i)),Y(1:mi,i),mi,work,lwork,info)
         endif
           Y(RH%dmi(ch(1,i))-RH%pp(ch(1,i))+1:RH%dmi(ch(1,i)),ch(1,i))=Y(1:RH%pp(ch(1,i)),i)
           Y(RH%dmi(ch(2,i))-RH%pp(ch(2,i))+1:RH%dmi(ch(2,i)),ch(2,i))=Y(RH%pp(ch(1,i))+1:mi,i)
     else
        if (RH%d(1,i) /= 0) then
           mi=RH%dmi(i)
           rk=mi-RH%pp(i)
           call dgemv('T',RH%pp(i), rk, -1d0,R(RH%p(1,i)+rk), mi,Y(1+rk:mi,i),1,1d0,Y(1:rk,i),1)
           call dtrsv('L','T','N',rk,R(RH%p(1,i)),mi,Y(1:rk,i),1)
           tb(1:mi)=Y(1:mi,i)

           call dormql('L','N',mi,1,mi,R(RH%p(1,i)),mi,R(RH%p(2,i)),Y(1:mi,i),mi,work,lwork,info)
        endif
     endif
  enddo

  do i=1,nn
     if (ch(1,i) == 0 .and. ch(2,i) == 0) then
        mi=RH%dmi(i)
        b(l(i,1):l(i,2))=Y(1:mi,i)
     endif
  enddo

 deallocate(ch,Y,l,tb)

  end subroutine solvecholsvd
     
!!!!!!!!
  subroutine cg(A, lda, n, b, tol, nit, RH, R, m, tr, ltr, Work, lwork, prec)

       use aux_hss
       integer :: lda, n, m(*), tr(*), ltr, nit, prec
       integer :: lwork, ni
       real(8) :: A(lda,*), b(*), R(*), Work(*)
       type(hss)     :: RH

       real(8) :: tol,dt,nu,mu, tb,res, time
       real(8), allocatable, dimension(:) :: rs,p,y,z,e

       allocate (rs(n), p(n), y(n), z(n), e(n))

       time = comp_time()

       ni = 0
       tb = dnorm2(n,b,1)
       p(1:n) = b(1:n)
       rs(1:n) = b(1:n)
       b(1:n) = 0
       if (prec == 1) call solvecholsvd(RH, R, m, tr, ltr, p, lda, Work, lwork)
       y(1:n) = p(1:n)

       ! e(1:n) = b(1:n)-x0(1:n)
       ! tt = dnorm2(n,x0,1)
       ! err = dnorm2(n,e,1)/tt
       res = dnorm2(n,rs,1)/tb
       
       do while (res >= tol .and. ni <= nit)
          ni = ni+1
          call DGEMV('N',n,n,1d0,A,n,p,1,0d0,z,1)
          dt = dot_product(y(1:n),rs(1:n))
          nu = dt/dot_product(p(1:n),z(1:n))
          b(1:n) = b(1:n)+nu*p(1:n)
          rs(1:n) = rs(1:n)-nu*z(1:n)
          y(1:n) = rs(1:n)
          if (prec == 1) call solvecholsvd(RH, R, m, tr, ltr, y, lda, Work, lwork)
          mu = dot_product(y(1:n),rs(1:n))/dt
          p(1:n) = y(1:n)+mu*p(1:n)

          !e(1:n) = b(1:n)-x0(1:n)
          !err = dnorm2(n,e,1)/tt
          res = dnorm2(n,rs,1)/tb
          !if (mod(ni,int(floor(nit/5.))) == 0) 
          print '(a10,i6,a22,e11.4)','CG step ', ni, ' relative residual: ',res
	enddo
       time = comp_time()-time

       print '(a10,i4,a30,e11.4)','>> After',ni,' CG steps, relative residual: ',res
       print*, 'PCG time:',time

       deallocate (rs, p, y, z, e)
     end subroutine cg
!
!
end module mat2hsssvd

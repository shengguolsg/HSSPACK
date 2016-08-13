module aux_hss
! 
! This module contains the basic structure of HSS matrices and construction procedure.
!  
  implicit none

! **************************************************
!        Types for HSS Cholesky Factorization      *
! **************************************************
       type hss
          integer,    allocatable, dimension(:) :: dmi !block size of each node
          integer,    allocatable, dimension(:) :: pp !remaining part, dimension: 1-by-nodes
          integer(8), allocatable, dimension(:,:) :: p
          ! p(1,:):   D, ..., p(2,:):   Q (Tau)
          integer,  allocatable, dimension(:,:) :: d ! 4xnn(i): ([rdim; cdim], i)
          ! d(1:2,:): D, ..., d(3:4,:): Q (Tau)
       end type

       type ddcell
          integer, allocatable, dimension(:,:) :: pmi
       end type ddcell

! **************************************************
!        Types for HSS Mat-Mat Multiplication      *
! **************************************************

       type HSSMM
!         for HSS matrix-matrix multiplication
          integer, allocatable, dimension(:)   :: D ! Diagonal block size of each node
          integer, allocatable, dimension(:,:) :: U ! (rdim, cdim) of U for each node, 2 x nn
          integer, allocatable, dimension(:,:) :: V ! (rdim, cdim) of V for each node, 2 x nn
          integer, allocatable, dimension(:,:) :: B ! (rdim, cdim) of B for each node, 2 x nn

          integer, allocatable, dimension(:,:) :: ph
          ! ph(1,:): U; ph(2,:):V; ph(3,:): B
          integer, allocatable, dimension(:,:) :: pd
          ! pd(:): D;
       end type HSSMM

       type HTinfo
          integer, allocatable :: ttr(:,:)  ! level tree
          integer, allocatable :: lentr(:)  ! length of each level
          integer, allocatable :: lfx(:)    ! row partition of X
          integer, allocatable :: ch(:,:)   ! children info
          integer  :: ln                    ! number of leaf nodes
          integer  :: pn                    ! number of parent nodes
          integer  :: lvl                   ! number of levels
       end type HTinfo

contains

!!!!!!!!!!! expand hss with zero pointers and sizes, a{i} has ni nodes
       subroutine hssexp0(a, b, ni)
       type(hss)    :: a
       type(ddcell) :: b
       integer      :: ni
!    
       a%dmi(1:ni) = 0
       a%pp(1:ni) = 0
       a%p(1:2,1:ni) = 0
       a%d(1:4,1:ni) = 0
       b%pmi(1:2,1:ni)=0

       end subroutine hssexp0


!!!!!!!!!!! for RRQR
subroutine hssexp(a, av, lav, j, b, ldb, mb, nb, copy)
!       use aux_basic
       type(hss) :: a
       real(8) :: av(*), b(*)
       integer :: lav,j,nb,mb,ldb,ip
       logical :: copy
       character :: durb

       ip = j 
       a%dmi(ip)=nb
       if (copy) then
	     a%p(1,ip) = lav !D(Q)'s starting point
	     lav = lav+nb*nb !#(Tau) = ld(av)
	     a%p(2,ip)=lav   !Tau's starting point
	     lav=lav+nb      !the length of Tau
	     call dlacpy('A', mb, nb, b, ldb, av(a%p(2,ip)), mb) !choose 'A' is better

	     a%d(1:2,ip) = (/nb, nb/) !sometime nb=0
	     a%d(3:4,ip) = (/ldb,nb/) !ldb=mb=1
	else     
	     a%p(1,ip)=lav
	     lav = lav+nb*nb !pnr increase, but no Tau
	     a%p(2,ip)=lav
	     
	endif       

       end subroutine hssexp

!!!!!!!!!!! 
       subroutine hssexpsvd(RH, R, pnr, i, Tau, pir, prk, mi, DD, RT, pnd, ww, copy)
       implicit none

       integer :: pir, prk, mi, info, err,i, ip, j,k
       type(hss) :: RH 
       type(ddcell)  :: RT
       real(8) :: R(*), Tau(*), DD(*), bb, ww(prk,*)  !R, Tau and DD
       !real(8), allocatable, dimension(:,:) :: Ttemp
       integer :: pnr,pnd,lwork
       logical :: copy

       !real(8) :: Ttemp(prk,prk)

       !allocate(Ttemp(prk,prk), stat=err)
       !print *, 'err=', err
       
       lwork=mi*mi
       ip = i 
       RH%dmi(ip)=mi
       if (copy) then
             RH%p(1,ip) = pnr  !D's starting point
             pnr = pnr+mi*mi  
             RH%p(2,ip)=pnr    !Tau's starting point
             pnr=pnr+mi
             call dlacpy('A',1,mi,Tau,1,R(RH%p(2,ip)),1)

             !compute DD{i}, ok
             ww(1:prk,1:prk)=0.0d+0
             call dlacpy('L',prk,prk,R(RH%p(1,ip)+(1+mi)*pir),mi,ww,prk) 
             !Ttemp=transpose(Ttemp)
             do k =1,prk
                do j=k,prk
                   bb=ww(j,k)
                   ww(j,k)=ww(k,j)
                   ww(k,j)=bb
                end do
            end do
             call dtrmm('L','L','N','N',prk,prk,1d0,R(RH%p(1,ip)+(1+mi)*pir),mi,ww,prk)
             call dlacpy('A',prk,prk,ww,prk,DD(pnd),prk)
             RT%pmi(1,ip)=pnd
             RT%pmi(2,ip)=prk
             pnd=pnd+prk*prk
             
             RH%d(1:2,ip) = (/mi, pir/) !sometime nb=0
             RH%d(3:4,ip) = (/1,mi/)    !useless
       else
            !not compressed
            RT%pmi(1,i)=pnd
            RT%pmi(2,i)=mi
            pnd=pnd+mi*mi
            RH%dmi(i)=mi
            RH%pp(i)=mi
       endif       
        
       !deallocate(Ttemp)
       end subroutine hssexpsvd

! **************************************************
!        Types for HSS Mat-Mat Multiplication      *
! **************************************************

!!!!! GetHssInfo()
  subroutine GetHssInfo(PH,M,ltr,TRE,ll, sh_rank ) 
!
! .. Scalar Arguments ..
    Character(Len=4) :: sh_rank
    INTEGER           :: LTR, ll

! .. Array Arguments ..
    TYPE(HSSMM)  :: PH
    TYPE(HTinfo) :: TRE
    INTEGER      :: M(*)
!
! Purpose 
! ========
! Print out the information of HSS matrix, the dim of generators, U, V and B.
! 
! .. Parameters ..
! PH  (INPUT)  TYPE(HSSMM)  contains all these generators
!
! M   (INPUT)  INTEGER Array, DIMENSION(LTR)
! 
! ltr (INPUT)  INTEGER, total number of nodes 
!
! TRE (INPUT)  TYPE(HTinfo), contains all the information of HSS tree
! 
! ll  (INPUT)  INTEGER, show all the information of the bottom ll levels
! 
! sh_rank (INPUT) CHARACTER STRING
! 
! ==============
! Written by Shengguo Li, on Aug. 15, 2012.
! ==============

! .. Scalar Parameters ..
      LOGICAL  ::  SHOW
      INTEGER  ::  STN, Ni, I, N, lvl, J
! .. Array Parameters ..
      INTEGER, ALLOCATABLE :: ViSet(:)

!  .. External Functions
      LOGICAL     LSAME
      EXTERNAL    LSAME
      
      allocate( ViSet(LTR) )
      SHOW = LSAME(sh_rank, 'SHOW')
      lvl = TRE%lvl

      J = 0
      DO I = 1, ltr
         IF( TRE%ch(1,I) .eq. 0) THEN
            ViSet(I) = J
            J = J+1
         ELSE
            J = J-1
            ViSet(I) = J
         END IF            
      END DO

      Ni = ltr -1
      STN = PH%ph(3,Ni) + PH%B(1,Ni)*PH%B(2,Ni)
      print *, 'The total storage of generators is ', STN
      print *, 'There are ', TRE%ln, 'leaf nodes', TRE%pn, 'parent nodes.'

      DO I = lvl, lvl-ll+1, -1
         N = TRE%lentr(I)
         print *, 'The row and column dim of U at level ', I
         DO J = 1, N
            print *, PH%U( 1:2, TRE%ttr(I,J) ), ViSet( TRE%ttr(I, J) ), TRE%ttr(I,J)
         END DO

!!$       print *, 'The row dims of VT at level ', I
!!$       print *, PH%V( 1,TRE%ttr(I,1:N) )
!!$
!!$       print *, 'The column dims of B at level ', I
!!$       print *, PH%B( 2,TRE%ttr(I,1:N) )
    END DO

  end subroutine GetHssInfo

!!!!! GetTreeInfo()
       subroutine GetTreeInfo(TR, LTR, RM, CM, TRE, LVL )
! 
! .. Scalar Parameters ..
     INTEGER  :: LTR, LVL

! .. Array Parameters ..
     INTEGER      :: TR(*), RM(*), CM(*)
     TYPE(HTinfo) :: TRE

! Purpose 
! ============
!    GetTreeInfo returns some information of HSS tree, TR. 
!    All these information are in Type HTinfo. 
!
! .. Parameters ..
! TR   (INPUT)  INTEGER Array, DIMENSION(LTR), post ordering HSS tree
!
! LTR  (INPUT)  INTEGER, the total number of nodes
!
! RM   (INPUT)  INTEGER Array, DIMENSION(N1), N1 is number of leaf nodes
!               row partition of matrix A
! CM   (INPUT)  INTEGER Array, DIMENSION(N1), N1 is number of leaf nodes
!               column partition of matrix A
! TRE  (OUTPUT) TYPE HTinfo, contains all the informaiton of TR.
! 
! ==================
! Written by Shengguo Li, on Aug. 14th, 2012
!===================
     
! .. Local Scalars ..
     INTEGER   :: I, J, It
     
! .. Local Arrays ..
     INTEGER   :: LV(ltr)

! .. Executive Statements ..
     
     call hsslevel(lv, tr, ltr)

     ! Tranform postordering to level-by-level
       TRE%ttr = 0
       do i =  lvl, 1, -1
          it = 1
          do j = 1, ltr-1
             if (lv(j) == i) then
                TRE%ttr(i, it ) = j
                it = it +1
             end if
          end do
          TRE%lentr(i) = it -1 
       end do

       ! Children
       call child(tr, TRE%ch, ltr)

       it = 1
       J = 1
       TRE%lfx = 0
       do i = 1, ltr
          if (TRE%ch(1,i) == 0 .and. TRE%ch(2,i) == 0 ) then
             TRE%lfx(i) = J
             J = CM(it) + J
             it = it +1
          end if
       end do
       TRE%ln = it -1       ! number of leaf nodes
       TRE%pn = ltr-(it-1)  ! number of parent nodes
       TRE%lvl = lvl        ! number of levels

       end subroutine GetTreeInfo

!!!!!!! hssblkrc()
       subroutine hssblkrc(D,W,Z,alpha,T,M,N,LDT)
!
! .. Scalar Arguments ..
         INTEGER    M,N,LDT
! .. Array Arguments ..
         DOUBLE PRECISION  D(*),W(*),Z(*),alpha(*),T(*)
!
! Purpose
! =========
! Construct the HSS block row or HSS block column
!
! .. Parameters ..
! D,W,Z,alpha  (INPUT)  DOUBLE PRECISION Array, DIMENSION 1D, the generators of Cauchy-like matrix
!              see mat2hssCauchy for details
! T   (INPUT/OUTPUT) DOUBLE PRECISION Array, DIMENSION 1D, contains the HSS block row or HSS block column
! 
! M   (INPUT)  INTEGER, the row dimension of off-diagonal block, part of T
! 
! N   (INPUT)  INTEGER, the col dimension of off-diagonal block, part of T
! 
! LDT (INPUT)  INTEGER, the leading Dim of T
!
! LK  (INPUT)  INTEGER, the starting position in T
! 
! =================
! Written by Shengguo Li, on Aug. 25th, 2012
! =================
         INTEGER   I,J
         DOUBLE PRECISION WW(M,N)
       
         DO I = 1, N
            DO J = 1, M
               WW(J,I) = alpha(J)* Z(I)/(( D(I)+W(j) )*( D(I)-W(j) ) )
            END DO
         END DO
         
         call dlacpy('A', M, N,WW,M, T, LDT)    

       end subroutine hssblkrc


!!!!!!! hssexpmm0, initial the HSSMM TYPE data
       subroutine hssexpmm0(a, ni)
       type(hssmm)    :: a
       integer      :: ni

       a%D(1:ni) = 0
       a%U(1:2,1:ni) = 0
       a%V(1:2,1:ni) = 0
       a%B(1:2,1:ni) = 0
       a%ph(1:3,1:ni)= 0
       a%pd(1:2,1:ni)= 0

     end subroutine hssexpmm0

!!!!!!! hssexpmmd()
     subroutine hssexpmmd(SD,DD,PH,pnd,D,W,Z,alpha,M,N,IP)
!
! .. Scalar Arguments ..
       CHARACTER      SD
       INTEGER        M, N, IP, pnd
! .. Array Arguments ..
       DOUBLE PRECISION  DD(*), D(*), W(*),Z(*),alpha(*)
       TYPE(HSSMM)       PH
! 
! Purpose 
! =========
! Store the information of main diagonal block generator, D. This subroutine is only for updating-SVD problem, 
! Cauchy-like matrices. 
! 
! .. Arguments ..
! SD    (INPUT)  CHARACTER, 'D' is the only choice
!
! PH    (INPUT/OUTPUT) HSSMM TYPE, stores all the information of generator D. It will be updated. 
!   
! PND   (INPUT/OUTPUT) INTEGER, points to the location of the end of array, H. 
!       Will be updated afterward.
!
! DD    (INPUT/OUTPUT) DOUBLE PRECISION array, DIMENSION(:), workspace which stores generator D
!
! M     (INPUT)  INTEGER, row dimension of this new generator
!
! N     (INPUT)  INTEGER, column dimension of this new generator
!
! IP    (INPUT)  INTEGER, the index of node 
!
! INFO  (INPUT/OUTPUT) INTEGER, check whether the inputs are correct
!
! ======================
! Written by Shengguo Li, on Aug. 25th, 2012
! ======================
       INTEGER  I, J
       DOUBLE PRECISION WW(M,N)
       
       DO I = 1, N
          DO J = 1, M
             WW(J,I) = alpha(J)* Z(I)/(( D(I)+W(j) )*( D(I)-W(j) ) )
          END DO
       END DO
       
       call dlacpy('A', M, N,WW,M, DD(pnd), M)    ! copy D{i}
       PH%D(IP) = M
       PH%pd(1,IP) = pnd
       PH%pd(2,IP) = pnd + M*N -1
       pnd = pnd + M*N

     end subroutine hssexpmmd

!!!!!!!  hssexp()
       subroutine hssexpmm2(UVB, PH, PNH, M, N, IP, INFO )
!
!  .. Scalar Arguments ..      
      CHARACTER         UVB
      INTEGER           PNH, M, N, IP, INFO
!  .. Array Arguments ..
      TYPE(HSSMM)       PH
!
!   Purpose
!   =======
!   Store the information of generators, U, V and B. D will be stored explicitly during the construction in mat2hss.
!   U, V and B are stored together and D is stored separated. This strategy is designed for
!   HSS matrix-matrix multiplication. 
!  
! .. Arguments ..
!   UVB  (INPUT)  CHARACTOR, 
!       'U','V','B': store the other generators, U, V and B.
! 
! PH    (INPUT/OUTPUT) HSSMM TYPE, stores all the information of generators. It will be updated. 
!   
! PNH   (INPUT/OUTPUT) INTEGER, points to the location of the end of array, H. 
!       Will be updated afterward.
!
! M     (INPUT)  INTEGER, row dimension of this new generator
!
! N     (INPUT)  INTEGER, column dimension of this new generator
!
! IP    (INPUT)  INTEGER, the index of node 
!
! INFO  (INPUT/OUTPUT) INTEGER, check whether the inputs are correct
!
! =============
! Written by Shengguo Li, on Aug 9th. 2012
! =============

!  .. Local Parameters
      LOGICAL     GU, GV, GB
!
!  .. External Functions
      LOGICAL     LSAME
      EXTERNAL    LSAME
!
! .. Executable Statements ..
!
      info = 0
      GU = LSAME( UVB, 'U')
      GV = LSAME( UVB, 'V')
      GB = LSAME( UVB, 'B')
      
      IF (.NOT. (GU .OR. GV .OR. GB ) ) THEN
         INFO = -1
      END IF
      
      IF ( INFO .EQ. 0) THEN
         IF ( GU ) THEN        ! store U
            PH%U(1,IP) = M      ! row dim
            PH%U(2,IP) = N      ! col dim
            PH%ph(1,IP) = pnh   ! location of U
            pnh = pnh + M*N    ! update the pointer 
         ELSE
            IF ( GV ) THEN
               PH%V(1,IP) = M      ! row dim
               PH%V(2,IP) = N      ! col dim
               PH%ph(2,IP) = pnh   ! location of V
               pnh = pnh + M*N    ! update the pointer 
            ELSE
               PH%B(1,IP) = M      ! row dim
               PH%B(2,IP) = N      ! col dim
               PH%ph(3,IP) = pnh   ! location of U
               pnh = pnh + M*N    ! update the pointer 
            END IF
         END IF
      END IF
      
       end subroutine hssexpmm2
     
!!!!!!!  hssexp()
       subroutine hssexpmm(UVB, PH, PNH, H, M, N, IP, WW, LDW, INFO)
!
!  .. Scalar Arguments ..      
      CHARACTER         UVB
      INTEGER           PNH, M, N, IP,INFO,LDW
!  .. Array Arguments ..
      DOUBLE PRECISION  H(*), WW(*)
      TYPE(HSSMM)       PH
!
!   Purpose
!   =======
!   Store the information of generators, U, V and B. D will be stored explicitly during the construction in mat2hss.
!   U, V and B are stored together and D is stored separated. This strategy is designed for
!   HSS matrix-matrix multiplication. 
!  
! .. Arguments ..
!   UVB  (INPUT)  CHARACTOR, 
!       'U','V','B': store the other generators, U, V and B.
! 
! PH    (INPUT/OUTPUT) HSSMM TYPE, stores all the information of generators. It will be updated. 
!   
! PNH   (INPUT/OUTPUT) INTEGER, points to the location of the end of array, H. 
!       Will be updated afterward.
!
! H     (INPUT/OUTPUT) DOUBLE PRECISION array, DIMENSION(:), workspace which stores all generators.
!
! M     (INPUT)  INTEGER, row dimension of this new generator
!
! N     (INPUT)  INTEGER, column dimension of this new generator
!
! IP    (INPUT)  INTEGER, the index of node 
!
! WW    (INPUT)  DOUBLE PRECISION array, DIMENSION(:), stores the new generator which will be copied to H. 
!
! LDW   (INPUT)  INTEGER, Leading dimension of WW
!
! INFO  (INPUT/OUTPUT) INTEGER, check whether the inputs are correct
!
! =============
! Written by Shengguo Li, on Aug 9th. 2012
! =============

!  .. Local Parameters
      LOGICAL     GU, GV, GB
!
!  .. External Functions
      LOGICAL     LSAME
      EXTERNAL    LSAME, DLACPY
!
! .. Executable Statements ..
!
      GU = LSAME( UVB, 'U')
      GV = LSAME( UVB, 'V')
      GB = LSAME( UVB, 'B')
      
      info = 0
      IF (.NOT. (GU .OR. GV .OR. GB ) ) THEN
         INFO = -1
      END IF
      
      IF ( INFO .EQ. 0) THEN
         IF ( GU ) THEN        ! store U
            PH%U(1,IP) = M      ! row dim
            PH%U(2,IP) = N      ! col dim
            PH%ph(1,IP) = pnh   ! location of U
            call dlacpy('A',M,N,WW,LDW,H(pnh),M)   ! copy U to H
            pnh = pnh + M*N    ! update the pointer 
         ELSE
            IF ( GV ) THEN
               PH%V(1,IP) = M      ! row dim
               PH%V(2,IP) = N      ! col dim
               PH%ph(2,IP) = pnh   ! location of V
               call dlacpy('A',M,N,WW,LDW,H(pnh),M)   ! copy V to H
               pnh = pnh + M*N    ! update the pointer 
            ELSE
               PH%B(1,IP) = M      ! row dim
               PH%B(2,IP) = N      ! col dim
               PH%ph(3,IP) = pnh   ! location of U
               call dlacpy('A',M,N,WW,LDW,H(pnh),M)   ! copy U to H
               pnh = pnh + M*N    ! update the pointer 
            END IF
         END IF
      END IF
      
       end subroutine hssexpmm

!!!!!!!   ConHssTree(x, y, n)  
       subroutine ConHssTree(x, y, ln, rm, tr, N, ltr)
!       Construct HSS tree for Cauchy matrix
!
!  .. Parameter
!   x  (input) DOUBLE PRECISION array, DIMENSION(N)
! 
!   y  (input) DOUBLE PRECISION array, DIMENSION(N)
! 
!  ln  (input) INTEGER, the number of leaf nodes
!      It may become smaller afterward. 
!
!  rm  (input/output)  INTEGER array, DIMENSION(ln)
!
!  tr  (input/output)  INTEGER array, DIMENSION(N)
!
!   N  (input)  INTEGER, the number of particles
!
!  ltr (input/output) INTEGER,  total number of nodes
!===============
!
!  Written by Shengguo Li, on Aug. 8th, 2012
!
!===============
         implicit none
!            
         double precision x(N), y(N)  ! (1/x_i - y_j)
         integer rm(ln), tr(*)
         integer  ln                ! number of leaf nodes
         integer     ltr               ! total number of nodes
!
         double precision a, b, gap, temp
         integer N                    ! length of x
         integer i, j1, msize, flag
         integer, allocatable :: cm(:)
         allocate( cm(ln) )

         a = minval( x(1:N) )
         temp = minval( y(1:N) )
         if ( temp < a ) a = temp
         b = maxval( x(1:N) )
         temp = maxval( y(1:N) )
         if (temp > b) b = temp
         gap = (b-a) / ln

         j1 = 1
         rm = 0
         cm = 0
         do i = 1, N

            j1 = ceiling( (x(i)-a)/gap )
            rm(j1) = rm(j1)+1

         end do !(particles i)
!
!  modify this tree if some rm(i) is too small, less than 20
         j1 = 1
         i = 1
         flag = 0
         msize = 20
         cm(1) = rm(1)
         do i = 1, ln-1
            if ( cm(j1) > msize ) then
               j1 = j1 + 1
               cm(j1) = rm(i+1)
            else
               cm(j1) = cm(j1)+rm(i+1)
               flag = 1
            end if
         end do
         
         if (flag == 1) then
            ln = j1
            rm(1:j1) = cm(1:j1)
         end if
         ltr = 2*ln -1

        call n2tree(2*ln-1, tr)

       end subroutine ConHssTree
         
!!!!!! hsslevel(lv, tr)
       recursive subroutine hsslevel(lv, tr, ncase)
!  Transform postordering HSS tree to level-by-level tree
!  
!  .. Parameter
!  lv   (input/output) INTEGER array, DIMENSION(nn)
!
!  tr   (input) INTEGER array, DIMENSION(nn)
!       nn is the number of nodes include the root node, nn = 2*ln-1 for binary tree.
!       ln   : number of leaf nodes
! 
! ch    (input) INTEGER array, DIMENSION(nn, 2)
!       the index of child of node i
!
! =====================
!
!  written by Shengguo Li, on Aug. 8th, 2012
!
!======================
         implicit none
!         
         integer tr(ncase), lv(ncase)
!
         integer ncase, lc, lt1, lt2
         integer, allocatable, dimension(:) :: lv1, lv2
         integer, allocatable, dimension(:,:) :: ch
         
         select case(ncase)
            case(0)  
               write(*,*) 'Error: Tree is NULL'
            case(1)
               lv(1) = 0
            case(2)
               lv(1:2)=(/1,0/)
            case(3)
               lv(1:3)=(/1,1,0/)
            case default
               allocate( ch(2, ncase) )
               call child(tr, ch, ncase)
               lc = ch(1, ncase)
               do while ( ch(1, lc) /= 0 .and. ch(2,lc) /= 0 )
                  lc = ch(1, lc)
               end do
!
               if ( ch(2, ncase) /= 0 ) then
                  lt1 = ch(1, ncase)-lc+1
                  lt2 = ch(2, ncase)-ch(1, ncase)
                  allocate( lv1(lt1),  lv2(lt2) )
                  call hsslevel(lv1, tr(lc:lc+lt1-1), lt1 )
                  call hsslevel(lv2, tr( ch(1, ncase)+1:ch(2,ncase) )-ch(1,ncase), lt2 )
                  lv( 1:lt1 ) = lv1( 1:lt1 ) +1 
                  lv( lt1+1:lt1+lt2 ) = lv2( 1:lt2 ) + 1
                  lv( ncase ) = 0

                  deallocate(ch, lv1, lv2)

               else
                  if (ch(1, ncase) /= 0 .and. ch(2, ncase) == 0 ) then  ! this must never happen
                     lt1 = ch(1, ncase)-lc+1
                     allocate( lv1(lt1) ) 
                     call hsslevel( lv1, tr( lc:ch(1, ncase) ), lt1 )
                     lv(1:lt1) = lv1(1:lt1) +1 
                     lv(ncase) = 0

                     deallocate(ch, lv1)

                  end if
               end if

            end select
!            
       end subroutine hsslevel
         
!!!!!!!!!!! n2tree(p, n)
       recursive subroutine n2tree(n, tr) ! not well balanced, need to improve
!   Construct post ordering binary tree
! 
       implicit none

       integer tr(*)
       integer n, n1

       if (mod(n, 2) == 0) stop "Input # of tree nodes must be odd"

       if (n == 1) then
          tr(1) = 0
          return
       else
          n1 = 2**floor(log(real(n))/log(2.))-1
          call btree(n1, tr)
          call n2tree(n-n1-1, tr(n1+1))
          tr(n1+1:n-1) = tr(n1+1:n-1)+n1
          tr(n1) = n
          tr(n-1:n) = (/n, 0/)
       endif

       end subroutine n2tree

!!!!!!!!!!! btree(p, n)
       recursive subroutine btree(n, tr)

       integer tr(*)
       integer m, n

       if (n.gt.3) then
          m = (n-1)/2
          call btree(m, tr)
          tr(m+1:m+m) = tr(1:m)+m
          tr(m) = n
          tr(m+m) = n
          tr(n) = 0
       else
          tr(1:3) = (/3,3,0/)
       end if

       end subroutine btree

!!!!!!!!!!! child(tr, ch, n)
       subroutine child(tr, ch, n)
         
       integer :: tr(*), ch(2,*), n
       integer i

       do i = 1,n
          ch(1,i) = 0
          ch(2,i) = 0
       end do

       !print '(20i5)',tr(1:n)
       do i = 1,n-1
          if (ch(1,tr(i)).eq.0) then
             ch(1,tr(i)) = i
          else
             ch(2,tr(i)) = i
          end if
       end do

       end subroutine child

!!!!!!!!!!!
       subroutine npart(n, ni, tr, ltr, m, lm)
       integer :: n, ni, tr(*), ltr, m(*), lm, mr
       

       if (ltr == 0) then
          lm = n/ni
          m(1:lm) = ni
          mr = mod(n,ni)
       
          if (mr >= ni/2) then
             lm = lm+1
             m(lm) = mr
          else
             m(lm) = m(lm)+mr
          endif
          ltr = 2*lm-1
       else
          lm = (ltr+1)/2
          ni = n/lm
          m(1:lm) = ni
          m(lm) = m(lm)+mod(n,ni);
       endif
       call n2tree(ltr, tr)

       end subroutine npart

!!!!!!!!!!!!
       subroutine compr(m,n,A,lda,R,ldr,rk,typ,tol,ncol,nflops)
       integer :: m, n, lda, ldr, rk, typ, ncol
       real(8) :: A(lda,*), R(ldr,*), tol
!       real(8),allocatable,dimension(:,:) :: a0
       integer :: nflops

       if (typ == 1) rk = ncol
       call mgsclpv(M,N,A,LDA,R,LDR,typ,rk,TOL,nflops)

       end subroutine compr

!!!!!!!!!!! 
       SUBROUTINE mgsclpv(M,N,A,LDA,R,LDR,typ,RANK,TOL,nflops)

       IMPLICIT NONE
!!!!
!!!!     .. Scalar Arguments ..
!!!!     ..
       INTEGER  :: LDA,LDR,M,N,iii
       integer  :: typ, RANK, nflops
       DOUBLE PRECISION  TOL
!!!!     ..
!!!!     .. Array Arguments ..
!!!!     ..
       INTEGER            PIV(N)
       DOUBLE PRECISION   A(LDA,*),R(LDR,*),VN(N)
!!!!     ..
!!!!     .. Local Scalars ..
!!!!     ..
       INTEGER            I,S,J,MN,IMAX
       DOUBLE PRECISION   VJ,VMAX, ZERO
       PARAMETER          ( ZERO = 0.0D+0)
!!!!     ..
!!!!     .. Intrinsic Functions ..
       INTRINSIC          SQRT
!!!!     ..
!!!!     .. Executable Statements ..
!!!!

       nflops = 0
       MN         = MIN(M,N)
       PIV(1:MN)  = 0
       R(1:MN,1:N)= ZERO
       DO I = 1,N
           VN(i) = dot_product(A(1:M,I),A(1:M,I))
           nflops = nflops+2*M
       END DO
!!!!
!!!!     Compute factorization.
!!!!
       DO I = 1,MN
           IMAX = I 
           VMAX = VN(I)
           DO J = I + 1, N
               IF (VN(J) > VMAX) THEN
                   VMAX = VN(J)
                   IMAX = J
               END IF
           END DO
           IF (IMAX > I) THEN
               DO S = 1, M
                   VJ        = A(S,IMAX)
                   A(S,IMAX) = A(S,I)
                   A(S,I)    = VJ
               END DO
               DO S = 1, I - 1
                   VJ        = R(S, I)
                   R(S,I)    = R(S, IMAX)
                   R(S,IMAX) = VJ
               END DO
               PIV(I)      = IMAX
               VN(IMAX)    = VN(I)
           END IF
           R(I,I) = dot_product(A(1:M,I),A(1:M,I))
           nflops = nflops+2*M

           R(I,I)     = SQRT(R(I,I))
           if (r(i,i) < 1d-15) then ! important, may need to improve the error
              rank = max(i-1,1)
              iii = i ! permutation already done for i, need to recover
              goto 111
           endif

           A(1:M,I) = A(1:M,I)/R(I,I)
           nflops = nflops+M
           VN(I) = R(I,I)
           DO J = I + 1, N
               R(I,J) = dot_product(A(1:M,I),A(1:M,J))

               A(1:M,J) = A(1:M,J)-R(I,J)*A(1:M,I)
               nflops = nflops+4*M+2
               VN(J)    = VN(J) - R(I,J) * R(I,J)
           END DO

           IF (typ == 0 .and. R(I,I) < (TOL * R(1,1)) &
            .or. typ == 1 .and. i == rank) THEN
               RANK = I
               iii = i
               GO TO 111
           END IF
       END DO

       RANK = MN
       iii = rank ! needed when rank=mn

 111   CONTINUE
       DO I = iii, 1, - 1
           IMAX      = PIV(I)
           DO S = 1, MIN(RANK,RANK * IMAX)
               VJ        = R(S,I)
               R(S,I)    = R(S,IMAX)
               R(S,IMAX) = VJ
           END DO
       END DO

       END SUBROUTINE mgsclpv
         
!!!!!!!!!!
       SUBROUTINE DTRMMS(SIDE,UPLO,TRANSA,DIAG,M,ALPHA,A,LDA)
!      .. Scalar Arguments ..
         DOUBLE PRECISION ALPHA
         INTEGER LDA,M,N
         CHARACTER DIAG,SIDE,TRANSA,UPLO
!      .. Array Arguments ..
         DOUBLE PRECISION A(LDA,*)
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: B

!  Purpose
!  =======
!  DTRMMS performs one of the matrix-matrix operations
!
!  A:=alpha*op(A)*A or A:=alpha*A*op(A)
!
!  where alpha is a scalar, A is a unit, or non-unit, upper or lower triangular, 
!  m by m matrix and op(A) is one of  
!   
!    op(A) = A or op(A) = A'
!  =======================

!   .. External Functions ..
    LOGICAL LSAME
    EXTERNAL LSAME
!   ..
!   .. External Subroutines ..
    EXTERNAL XERBLA
!   ..
!   .. Intrinsic Functions ..
    INTRINSIC MAX
!   ..
!   .. Local Scalars ..
    LOGICAL LSIDE
!   ..
!   .. Parameters ..
    DOUBLE PRECISION ONE, ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!   ..
!   
!   Test the input parameters.
!

!    compute routine
     ALLOCATE(B(M,M))

     call DLACPY('L',M,M,A,LDA,B,M)
     call dtrmm('R','L','T','N',M,M,1d0,B,M,A,M)

     DEALLOCATE(B)

    END SUBROUTINE DTRMMS

subroutine transpose1d(m, n, A, lda, B, ldb)
       integer m, n, i,j, lda ,ldb
       real(8) :: A(lda,*), B(ldb,*)

       do i = 1, m
          do j = 1, n
             B(j,i) = A(i,j)
          end do
       end do
       end subroutine transpose1d

!!!!!!!!       
       subroutine transpose1ds(A,m,n)
         integer ::  m,n, i,j
         real(8) :: A(m,m), b

         do i =1,m
            do j=i,n
               b=A(j,i)
               A(j,i)=A(i,j)
               A(i,j)=b
             end do
          end do

       end subroutine transpose1ds

!!!!!!!!!!
       subroutine circshift2(A,lda,ldb,m,dn)

         integer :: lda,m,dn, ldb
         real(8) :: A(lda,*)

         A(1:lda,1:ldb)=cshift(A(1:lda,1:ldb),m,dn)
       end subroutine circshift2

!!!!!!!!!!
       subroutine Usigma(A,lda,S,k,nt,B)
         integer :: k,i, lda, nt
         real(8) :: A(lda,*), B(nt,k), S(*)

         do i=1,k
            A(1:lda,i)=A(1:lda,i)*S(i)
         end do
         
         call dlacpy('A',lda,k,A,lda,B,nt)
       end subroutine Usigma

!!!!!!
       subroutine Usigma2(A,lda,S,k)
         integer :: k,i, lda
         real(8) :: A(lda,*), S(*)

         do i=1,k
            A(1:lda,i)=A(1:lda,i)*S(i)
         end do
       end subroutine Usigma2

!!!!!!!!!!!
       function flops(typ, transA, mA, nA, transB, mB, nB)
       character :: typ, transA, transB
       integer(8) :: flops, t
       integer mA, nA, mB, nB

       select case (typ)
          case ('s') ! sum
             flops = mA*nA
          case ('c') ! chol
             t = ma
             flops = t**3/3-t*t+t*2/3
          case ('p') ! prod
             t = mA*nA
             if (transB == 'n') then
                flops = t*nB*2
             else
                flops = t*mB*2
             endif
          case ('t') ! tri sol
             t = mB*mB
             if (transA == 'n') then
                flops = t*mA
             else
                flops = t*nA
             endif
       end select

       end function flops

!!!!!!!!!!!
       function comp_time()
	
         real(8) :: comp_time
         integer :: time0(8)
         call date_and_time(values=time0)
         comp_time = time0(5) * 3600 + time0(6)*60 + time0(7) +0.001 * time0(8)
	
       end function comp_time

!!!!!!!!!!!
       function dnorm2(n,x,incx) result(d)
         integer :: n,incx, i
         real(8) :: x(*),d

         d = 0.0D+0
         do i = 1,n,incx
            d = d+x(i)*x(i)
         enddo
         d = sqrt(d)

       end function dnorm2

!!!!!!!!
!   Calculate log2 (n)
       function log2 (n) result (resultValue)
         implicit none
         integer, intent (in) :: n
         DOUBLE PRECISION :: resultValue
    
         resultValue = real (LOG (n+0.0) / LOG (2.0))
  
       end function log2

     end module aux_hss

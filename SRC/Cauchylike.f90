  subroutine Cauchylike(A,D,F,u,v,M,N)
!   A(i,j) = u(i)v(j) / (D(i)-F(j))

    integer M,N,i,j
    double precision A(M,*), D(*),F(*),u(*),v(*)

    do j = 1,N
       do i = 1,M
          A(i,j) = u(i)*v(j)/( D(i)-F(j) )
       end do
    end do
    
  end subroutine Cauchylike


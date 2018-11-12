program seidel
  implicit none
  integer,parameter:: n=2
  real*8 a(n,n),b(n),l(n,n),u(n,n),l_inv(n,n)
  real*8 t(n,n),c(n),x(n)
  integer*8 i,j
 
  ! a matrix
  write(6,*) "insert values of a matrix"
  do i=1,n
     do j=1,n
        read(5,*) a(i,j)
     end do
  end do
 
  ! b matrix
  write(6,*) "insert values of b matrix"
  do i=1,n
     read(5,*) b(i)
  end do
 
  call lower_component(a,n,l)
  call strict_upper(a,l,n,u)
  call inverse(l,l_inv,n)
  call tc_seidel(t,l_inv,u,n,c,b)
  x(1)=1.1
  x(2)=2.3
  do i=1,100
     x = matmul(t,x)+c
     print*, i,x
  end do
  
end program seidel
!**************************************************************************
subroutine lower_component(a,n,l)
  implicit none
  integer*4 n
  integer*8 i,j
  real*8 a(n,n),l(n,n)

  do i=1,n
     do j=1,n
        if (i >= j) then
           l(i,j) = a(i,j)
        else
           l(i,j) = 0
        end if
     end do
  end do
  !do i=1,n
   !  do j=1,n
    !    write(6,*) l(i,j)
    ! end do
 ! end do
end subroutine lower_component
!***************************************************************************
subroutine strict_upper(a,l,n,u)
  implicit none
  integer*4 n
  integer*8 i,j
  real*8 a(n,n),l(n,n),u(n,n)
  u = a-l
  !do i=1,n
   !  do j=1,n
    !    write(6,*) u(i,j)
    ! end do
  !end do
end subroutine strict_upper
!**************************************************************************
subroutine inverse(a,c,n)
implicit none 
integer n
real*8 a(n,n), c(n,n)
real*8 L(n,n), U(n,n), b(n), d(n), x(n)
real*8 coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do

 !do i=1,n
    ! do j=1,n
    !   write(6,*) c(i,j)
    ! end do
 ! end do
  
end subroutine inverse
!*******************************************************************************
subroutine tc_seidel(t,l_inv,u,n,c,b)
  implicit none
  integer*4 n
  integer*8 i,j
  real*8 t(n,n),l_inv(n,n),u(n,n),c(n),b(n)
  t = - matmul(l_inv,u)
  c = matmul(l_inv,b)

   do i=1,n
     do j=1,n
       write(6,*) t(i,j)
     end do
  end do
  
  do i=1,n
       write(6,*) c(i) 
  end do
end subroutine tc_seidel
!************************************************************************************




  
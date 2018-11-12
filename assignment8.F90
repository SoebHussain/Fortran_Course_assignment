  program asd
  implicit none
  real hp,Ta,h,m,u(0:100),l(0:100),c(0:100),T(0:100),d(0:100)
  integer n,i,j
  hp = 0.01
  Ta = 20
  h=1
  n=9
  T(0) = 40
  T(10) = 200

  do i=1,n
     d(i) = hp*h*h + 2
     u(i) = -1
     l(i) = -1
     c(i) = hp*h*h*Ta
  end do
  c(1)  = hp * h*h*Ta + T(0)
  c(n) = hp*h*h*Ta + T(10);

  do i=1,n
     m=-l(i+1)/d(i)
     d(i+1)=d(i+1)+m*u(i)
        c(i+1)=c(i+1)+m*c(i)
        l(i+1)=0
     end do
        
     T(n)=c(n)/d(n)
     do i=1,n-1
        j=n-i
        T(j)=(c(j)-u(j)*T(j+1))/d(j)
        
     end do

     do i=0,n+1
        write(6,*) "T[",i,"]",T(i)
     end do


   end program
program summation
  implicit none
  integer,parameter:: n= 200
  real*8 temp,x,sum,sum1
  integer i
  temp = 1
  sum=0
  print*, "write the value of x"
  read(5,*) x
  do i=1,n
     temp =(temp*x)/i
     sum1 = sum
     sum = sum + temp*(i+2)
     if(abs(sum-sum1).le.10e-12) then
        write(6,*) sum,i
     exit
     end if
  end do
end program summation


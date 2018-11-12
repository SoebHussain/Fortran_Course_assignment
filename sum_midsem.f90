program suma
  implicit none
  integer i,j
  integer , parameter:: n=100
  real*8 temp,sum,x,sum1
  sum=1
  temp=1
  sum1=0
  x=0.5
  j=0
  do i=1,n,2
     temp = (temp*i*(-1)*x)/(i+1)
     j=j+1
     sum = sum + temp
     print*,sum,sum1
     
     
   if (abs(sum-sum1).lt.10e-7) then
        write(6,*) j,sum,sum-sum1
     exit
    end if
      sum1 = sum
  end do

  write(6,*) sum

  
     
end program suma


program ass10
implicit none

integer n,i,seed,nrand
real h,a,b,f,ans1,f0,fn,sum,ans2,x,ans3,aa,c,m,ans4,rem

a=0
b=0.8
n=2
h=(b-a)/n

! single trapezoidal
call fun(a,f0)
call fun(b,fn)
ans1=(b-a)*(f0+fn)/2
write(6,*)ans1

! multiple trapezoidal
sum=f0
do i=1,n-1
   call fun(a+i*h,f)
   sum=sum+2*f
end do
sum=sum+fn
ans2=h*sum/2
write(6,*)ans2


aa=0.1
c=1
m=b-a
nrand=100

!sample mean method using modulo generator
sum=0
x=a
do i=1,nrand
   x=mod(aa*x+c,m)
   call fun(x,f)
   sum=sum+f
end do
ans3=(b-a)*sum/nrand
write(6,*)ans3


! using rand function
call srand(seed)
sum=0
do i=1,nrand
   x=a+(b-a)*rand()
   call fun(x,f)
   sum=sum+f
end do
ans4=(b-a)*sum/nrand
write(6,*)ans4

end program ass10


subroutine fun(x,f)
implicit none
real x,f
f=0.2+25*x-200*x**2+675*x**3-900*x**4+400*x**5
end subroutine fun
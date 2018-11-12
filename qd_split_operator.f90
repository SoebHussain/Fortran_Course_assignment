program split
  implicit none
  integer,parameter:: n=400 
  integer i,t,j,p
  real*8 xst,xend,dx,x(n),v(n),dt,temp
  complex*16 z,expov(n),a,psi(n),g(n),b,gd(n),pxp,php,pxhp,pvp
  dt=0.1
  z=(0,1)
  xst=-10
  xend=30
  a= -z*dt
  b= -a/2

  
  !****grid prepare
  call grid(xst,xend,dx,x,n)
  !print*, x

  !*****calculation of psi
  do i=1,n
     psi(i) = exp((-x(i)**2)/2 + z*x(i))
     !write(6,*) psi(i)
  end do

  !*****normalization of psi
  call normalise(psi,psi,dx,n)
  do i=1,n
     !write(6,*) psi(i)
  end do

  !****loop of time
  t=20
 p=20/dt
 do j=1,p
 
  !****calculation of v
  do i=1,n
     v(i)= 0!-(x(i)**2)/2
     !write(6,*) x(i),v(i)
  end do

  !****exp(av/2)
  do i=1,n
     expov(i) = 1 ! exp((a*v(i))/2)
     !write(6,*) expov(i)
  end do
  
  !***multiplication of expov and psi
  do i=1,n
     g(i)= expov(i) * psi(i)
    ! write(6,*) expov(i),psi(i),g(i)
  end do
  
  call differentiation_fourier(g,n,gd,x,b)
  !print*, gd
  
  !****final value of psi
  do i=1,n
     psi(i) = expov(i)*gd(i)
      !write(6,*) psi(i)
  end do
  temp = j*dt
  call integration(psi,x,dx,pxp,n)  !!calculation of <psi|x|psi>
  !print*, dx
  !write(27,*) temp,abs(pxp)
  call psihpsi(psi,n,x,dx,php)   !calculation of <psi|h|psi>
  !write(6,*) temp,php

  call psixhpsi(psi,n,x,dx,pxhp) ! calculation of <psi|xh|psi>
  !write(6,*) temp,abs(pxhp)

  call  psivpsi(psi,v,x,n,dx,pvp) ! calculation of <psi|v|psi>
  !write(6,*)  temp, pvp
  write(6,*) temp,abs(pvp+pxp) ! calculation of total energy
  
  
 end do

 

end program split
!********************subroutines**********************************************8

subroutine grid(qst,qend,dq,q,n)
  implicit none
  integer n,i
  real*8  qst,qend,q(n),dq
  dq=(qend-qst)/n
  do i=1,n
     q(i)=qst+(i-1)*dq
  enddo
end subroutine grid

!*******************

subroutine integral_c(integrand,n,dx,out1)
  implicit none
  integer n,i
  real*8   dx
  complex*16 integrand(n),out1

  out1=0
  do i=1,n
     out1=out1+integrand(i)*dx
  enddo
  
end subroutine integral_c

!*******************

subroutine normalise(input,output,dx,n)
implicit none
integer n,i
real*8 dx,normal
complex*16 integrand(n),input(n),output(n),out1
do i=1,n
integrand(i) = input(i)*conjg(input(i))
end do
call integral_c(integrand,n,dx,out1)
normal = 1/sqrt(out1)
do i=1,n
   output(i) = normal*input(i)
end do
end subroutine normalise 
!*******************
subroutine fourier(q1,q2,q1g,q2g,n,sign)
  implicit none
  integer n,i1,i2,sign
  complex*16 q1(n),q2(n),integrand(n),zi,out
  real*8 q1g(n),q2g(n),pi
  pi=dacos(-1.d0)
  zi=(0,1)
  do i2=1,n

     integrand=q1(:)*exp(sign*zi*q1g(:)*q2g(i2))
      call integral_c(integrand,n,q1g(2)-q1g(1),out)
      q2(i2)=out/sqrt(2*pi)
   enddo
   end subroutine fourier

!*******************
subroutine differentiation_fourier(f,n,fd,x,b)
  implicit none
  integer n,i
  complex*16 f(n),fd(n),g(n),z,b
  real*8 x(n),k(n),pi,dx,dk
  z=(0,1)
  dx=x(2)-x(1)
  pi=dacos(-1.d0)
   call grid(-pi/dx,pi/dx,dk,k,n)
  
   call fourier(f,g,x,k,n,-1)
   g=g*exp(b*((z*k)**2))
  call fourier(g,fd,k,x,n,+1)
!!do i=1,n ; write(8,*)x(i),real(f(i)) ; enddo
!! do i=1,n ; write(9,*)k(i),real(g(i)) ; enddo   
end subroutine differentiation_fourier
!*******************
subroutine integration(psi,x,dx,energy,n)
  implicit none
  integer i,n
  real*8 x(n),dx
  complex*16 integrand,psi(n),energy
  integrand =0
  energy= integrand
  !print*, dx
   do i=1,n
     integrand =integrand  + psi(i)*x(i)*conjg(psi(i))*dx
  end do
  energy= integrand
!  call integral_c(integrand,n,dx,energy)
end subroutine integration
!********************

subroutine differ_fourier(f,n,fd,x,order)!** with order
  implicit none
  integer n,order,i
  complex*16 f(n),fd(n),g(n),zi
  real*8 x(n),k(n),pi,dx,dk
  zi=(0,1)
  dx=x(2)-x(1)
  pi=dacos(-1.d0)
   call grid(-pi/dx,pi/dx,dk,k,n)
  
   call fourier(f,g,x,k,n,-1)
   g=g*(zi*k)**order
  call fourier(g,fd,k,x,n,+1)
do i=1,n ; write(8,*)x(i),real(f(i)) ; enddo
   do i=1,n ; write(9,*)k(i),real(g(i)) ; enddo   
end subroutine differ_fourier


!********************

subroutine psihpsi(psi,n,x,dx,out2)
  implicit none
  integer i,n
  real*8 x(n),dx
  complex*16 psi(n),psid(n),out2,out1,integrand
  
  call differ_fourier(psi,n,psid,x,+2)
  integrand = 0
  do i=1,n
     integrand = integrand + conjg(psi(i))*psid(i)*dx
  end do
  out1 = integrand

  out2= out1/2
end subroutine psihpsi

!******************
subroutine psixhpsi(psi,n,x,dx,out2)
  implicit none
  integer i,n
  real*8 x(n),dx
  complex*16 psi(n),psid(n),out2,out1,integrand
  
  call differ_fourier(psi,n,psid,x,+2)
  integrand = 0
  do i=1,n
     integrand = integrand + conjg(psi(i))*x(i)*psid(i)*dx
  end do
  out1 = integrand

  out2= out1/2
end subroutine psixhpsi
!************************
subroutine psivpsi(psi,v,x,n,dx,out2)
  implicit none
  integer i,n
  real*8 x(n),dx,v(n)
  complex*16 integrand,psi(n),out2
  integrand=0
  do i=1,n
     integrand = integrand + psi(i)*v(i)*conjg(psi(i))*dx
  end do

  out2 =integrand
end subroutine psivpsi

!**********************
  

     

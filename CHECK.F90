program main
  implicit none
  integer,parameter:: nx=500,np=25
  real*8  xst,xend,x(nx),dx,pi,hp(nx,0:np),integrand(nx),out,factorial,a,fr(nx),frd(nx)
  integer i,j
  complex*16 f(nx),f2d(nx),integrandc(nx),outc

  pi=dacos(-1.d0)
  xst=-25; xend=25;  call grid_prepare(xst,xend,x,nx)!;  write(6,*)x
  call hermite(hp,nx,np,x)
  do i=0,np
     do j=0,np
        integrand(:)=hp(:,i)*hp(:,j)*exp(-x(:)**2)
        call  integral_r(integrand,nx,x(2)-x(1),out)
        call check_hermite_ortho(out,i,j)
     enddo
  enddo

   fr(:)=hp(:,0)*exp(-x(:)**2/2)
   call diff(fr,nx,frd,x(2)-x(1))
   do i=1,nx
      write(10,*)x(i),frd(i)
write(11,*)x(i),(x(i)**2-1)*exp(-x(i)**2/2)
   enddo
   stop
  
  do i=0,np
     
     f=hp(:,i)*exp(-x(:)**2/2)
     call differentiation_fourier(f,nx,f2d,x,2)
     do j=1,nx
        write(7,*)x(j),real(f(j))*(x(j)**2-1),real(f2d(j))
     enddo
     integrandc=f*f2d
     call integral_c(integrandc,nx,x(2)-x(1),outc)
     integrand(:)=real(f**2)
     call  integral_r(integrand,nx,x(2)-x(1),out)
     write(6,*)i,outc/out
  enddo
end program main

subroutine check_hermite_ortho(out,i,j)
  implicit none
  integer i,j
  real*8  out,a,factorial,pi
  pi=dacos(-1.d0)
  if(i.ne.j)then
     if(abs(out).gt.1e-8)then
        write(6,*)i,j,'error'
     endif
  else
     a=sqrt(pi)*2**i*factorial(i)
     if(abs(out-a).gt.1e-8)then
        write(6,*)i,j,'error'
     endif
  endif
end subroutine check_hermite_ortho

function factorial(n)
  implicit none
  integer n,i
  real*8 factorial
  factorial=1
  do i=1,n
     factorial=factorial*i
  enddo
end function factorial
  
subroutine grid_prepare(qst,qend,q,n)
  implicit none
  integer n,i
  real*8  qst,qend,q(n),dq
  dq=(qend-qst)/n
  do i=1,n
     q(i)=qst+(i-1)*dq
  enddo
end subroutine grid_prepare


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

subroutine hermite(hp,nx,np,xg)
  implicit none
  integer nx,i,j,np
  real*8 hp(nx,0:np),xg(nx)
  do i=1,nx
     hp(i,0)=1
     hp(i,1)=2.d0*xg(i)
     do j=2,np
        hp(i,j)=2*xg(i)*hp(i,j-1)-2*(j-1)*hp(i,j-2)
     enddo
  enddo
end subroutine hermite

subroutine integral_r(integrand,n,dx,out)
  implicit none
  integer n,i
  real*8 integrand(n),dx,out
  out=0
  do i=1,n
     out=out+integrand(i)*dx
  enddo
  
end subroutine integral_r

subroutine integral_c(integrand,n,dx,out)
  implicit none
  integer n,i
  real*8   dx
  complex*16 integrand(n),out

  out=0
  do i=1,n
     out=out+integrand(i)*dx
  enddo
  
end subroutine integral_c



subroutine differentiation_fourier(f,n,fd,x,order)
  implicit none
  integer n,order,i
  complex*16 f(n),fd(n),g(n),zi
  real*8 x(n),k(n),pi,dx
  zi=(0,1)
  dx=x(2)-x(1)
  pi=dacos(-1.d0)
   call grid_prepare(-pi/dx,pi/dx,k,n)
  
   call fourier(f,g,x,k,n,-1)
   g=g*(zi*k)**order
  call fourier(g,fd,k,x,n,+1)
do i=1,n ; write(8,*)x(i),real(f(i)) ; enddo
   do i=1,n ; write(9,*)k(i),real(g(i)) ; enddo   
end subroutine differentiation_fourier





subroutine ea(a,u,n)
  implicit none
  integer n,i
  real*8 a(n,n),u(n,n),temp(n,n)
  u=0
  do i=1,n
     u(i,i)=1
  enddo
  temp=u
  do i=1,30
     temp=matmul(temp,a)/i
     u=u+temp
  enddo
end subroutine ea



subroutine diff(f,n,df,dx)
  implicit none
  integer n,i
  real*8 F(n),DF(n)
  real*8 dx
  df=0
  do i=2,n-1
     df(i)=(f(i+1)+f(i-1)-2*f(i))/dx**2
  enddo
end subroutine diff
     

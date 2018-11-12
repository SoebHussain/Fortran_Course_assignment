program main2
  call q3
end program main2



subroutine q2
  implicit none
  integer,parameter:: nx=500
  real*8  xst,xend,x(nx),pi
  complex*16 f(nx),f2d(nx),integrandc(nx),outc,zi,f1d(nx),outc1,outc2,outc3
  zi=(0,1)
  pi=dacos(-1.d0)
  xst=-25; xend=25;  call grid_prepare(xst,xend,x,nx)
  f=exp(-x**2/2+2*zi*x)
  
  integrandc=conjg(f)*x*f
  call integral_c(integrandc,nx,x(2)-x(1),outc)
  write(6,*)outc
  
  integrandc=conjg(f)*x*x*f
  call integral_c(integrandc,nx,x(2)-x(1),outc1)
  write(6,*)outc1
  
  call d1dx_4pt(f,f1d,nx,x(2)-x(1))
  integrandc=-zi*conjg(f)*f1d
  call integral_c(integrandc,nx,x(2)-x(1),outc2)
  write(6,*)outc2

  call d1dx_4pt(f1d,f2d,nx,x(2)-x(1))
  integrandc=-1*conjg(f)*f2d
  call integral_c(integrandc,nx,x(2)-x(1),outc3)
  write(6,*)outc3
  write(6,*)sqrt((outc**2-outc1)*(outc2**2-outc3))
  
  
  
  
end subroutine q2

subroutine d1dx_4pt(f,f1d,nx,dx)
  implicit none
  integer nx,i
  complex*16 f(nx),f1d(nx)
  real*8 dx
  f1d=0
  do i=3,nx-2
     f1d(i)=(-f(i+2)+8*f(i+1)  +f(i-2)-8*f(i-1))/(12*dx)
  enddo
end subroutine d1dx_4pt




subroutine Q1_1
  
  call Q1(.1d0,4)
  call Q1(.1d0,10)
  call Q1(.1d0,40000)
  
  call Q1(.5d0,4)
  call Q1(.5d0,10)
  call Q1(.5d0,40000)
  
  call Q1(.98d0,4)
  call Q1(.98d0,10)
  call Q1(.98d0,40000)
end subroutine Q1_1
  
Subroutine Q1(a,n)
  implicit none
  integer i,n
  real*8 a,u,temp
  u=1
  
  temp=u
  do i=1,n
     temp=-temp*a*(2*i-1)/2/i
     u=u+temp
     if (abs(temp) .lt. 1e-6)        goto 10
  enddo
  10 write(6,*)a,i,u
end Subroutine Q1




subroutine  q3
  implicit none
  integer,parameter:: nx=200,np=20
  real*8  xst,xend,x(nx),dx,pi,hp(nx,0:np),integrand(nx),out&
       ,factorial,a,fr(nx),frd(nx),h0mat(np+1,np+1),&
       xmat(np+1,np+1),f1(nx),t,psi(nx),dt
  integer i,j
  complex*16 f(nx),f2d(nx),integrandc(nx),outc,c(np+1),zi,hmat(np+1,np+1),u(np+1,np+1)
  zi=(0,1)
  dt=0.1d0
  pi=dacos(-1.d0)
  xst=-10; xend=10;  call grid_prepare(xst,xend,x,nx)!;  write(6,*)x
  psi=exp(-x**2/2)
  integrand=psi**2
  call  integral_r(integrand,nx,x(2)-x(1),out)
  psi=psi/sqrt(out)
  
  call hermite(hp,nx,np,x)
  
  do j=0,np
     i=j
     integrand(:)=hp(:,i)*hp(:,j)*exp(-x(:)**2)
     call  integral_r(integrand,nx,x(2)-x(1),out)
     !call check_hermite_ortho(out,i,j)
     hp(:,i)=hp(:,i)/sqrt(out)
     
     integrand(:)=psi*hp(:,j)*exp(-x(:)**2/2)
     call  integral_r(integrand,nx,x(2)-x(1),out)
     c(i+1)=out
  enddo
  !write(6,*)c; stop
  
  do i=0,np
     f=hp(:,i)*exp(-x(:)**2/2)
     call differentiation_fourier(f,nx,f2d,x,2)
     do j=0,np
        f1=hp(:,j)*exp(-x(:)**2/2)
        integrandc=-0.5d0*f1*f2d +0.5d0*f*f1*x**2
        call integral_c(integrandc,nx,x(2)-x(1),outc)
        h0mat(j+1,i+1)=outc
        !if(abs(outc).gt.1e-8)write(6,*)i,j,outc,'check'
        
        integrand(:)=hp(:,i)*hp(:,j)*exp(-x(:)**2)*x
        call  integral_r(integrand,nx,x(2)-x(1),out)
        xmat(j+1,i+1)=out
     enddo
  enddo
  do i=1,200
     t=i*dt
     hmat=-zi*dt*(h0mat+0.01d0*xmat*cos(t))
     call ea_complex(hmat,u,np+1)
     c=matmul(u,c)
     outc=sum(conjg(c)*matmul(h0mat,c))
     write(16,*)t,real(outc),aimag(outc)
     
  enddo


  
end subroutine q3

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

subroutine ea_complex(a,u,n)
  implicit none
  integer n,i
  complex*16 a(n,n),u(n,n),temp(n,n)
  u=0
  do i=1,n
     u(i,i)=1
  enddo
  temp=u
  do i=1,30
     temp=matmul(temp,a)/i
     u=u+temp
  enddo
end subroutine ea_complex

subroutine ex(a,u,n)
  implicit none
  integer i,n
  real*8 a,u,temp
  u=1

  temp=u
  do i=1,30
     temp=temp*a/i
     u=u+temp
  enddo
end subroutine ex

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
     
subroutine q5
  integer, parameter:: n=5
  integer i,j
  complex*16 a(n,n),ainv(n,n),tempa(n,n),zi

  zi=(0.d0,1.d0)
  
  
  do i=1,n
     do j=1,n
        if(i.eq.j)then
           a(i,j)=i+j
        elseif(i.gt.j)then
           a(i,j)=(i+j)*zi
        elseif(i.lt.j)then
           a(i,j)=(-zi*(i+j))
        endif
        
     enddo
  enddo


  tempa=a
call mat_inverse(a,ainv,n)

do i=1,n
   do j=1,n

write(6,*)i,j,a(i,j),ainv(i,j)
 enddo
enddo

tempa=matmul(tempa,ainv)
do i=1,n
   do j=1,n

if (abs( tempa(i,j)) .gt. 1e-10)         write(6,*)i,j,tempa(i,j)
 enddo
enddo

  
end subroutine q5

subroutine mat_inverse(a,ainv,n)
integer  n
integer i,j
complex*16 a(n,n),ainv(n,n),tempa(n,n),int,zi
zi=(0.d0,1.d0)


do i=1,n
   do j=1,n
      ainv(i,j)=0.d0
   enddo
   ainv(i,i)=1.d0
enddo



do i=1,n
   ainv(i,:)= ainv(i,:)/a(i,i)
   a(i,:)= a(i,:)/a(i,i)
   !write(*,*)a,ainv
   ! read(*,*)
   
   do j=1,n
      if( i.ne.j)then
         ainv(j,:)=ainv(j,:)-a(j,i)*ainv(i,:)
         a(j,:)=a(j,:)-a(j,i)*a(i,:)
      endif
      !write(*,*)a,ainv
      !read(*,*)
      
   enddo
enddo
!stop



end subroutine mat_inverse

subroutine q4
  implicit none
  integer,parameter:: n=2
  real*8 y(n),dt,t,y_old(n),t0,pi
  integer i,iy,iy2,it0
  
  dt=0.1d0
  do it0=1,50
     t0=20*dacos(-1.d0)/50*it0
     
     y=0.d0
     do i=1,1000
        t=i*dt
        y_old=y
        call rk4_f90(y,t,dt,n,t0)
        if (t.gt. .5d0 .and. y (1)*y_old(1) .lt. 0.d0) then
           write(6,*)t0,0.5d0*y(2)**2
           goto 10
        endif
        
     enddo
10   continue   
  enddo
end subroutine q4


subroutine ff(y,t,f,n,t0)
  implicit none
  integer n
  real*8 y(n),f(n),t,k,t0
!  k=1 ;f(1)=y(2); f(2)=-k*y(1) !single harmonic oscillator
!f(1)=y(3);f(2)=y(4);f(3)=-y(1);f(4)=-y(2) ! problem 1
!f(1)=y(3);f(2)=y(4);f(3)=-y(1)/sqrt(y(1)**2+y(2)**2);f(4)=-y(2)/sqrt(y(1)**2+y(2)**2) ! problem 2
!f(1)=y(3);f(2)=y(4);f(3)=-(y(1)+2*y(1)*y(2));f(4)=-(y(2)+y(1)**2-y(2)**2) ! problem 3

  f(1)=y(2); f(2)=cos(0.1*(t-t0))
end subroutine ff


subroutine rk4_f90(y,t,dt,n,t0)
  implicit none
  integer n,i
  real*8 y(n),f(n),t,k,dt,temp(n),k1(n),k2(n),k3(n),k4(n),t0
  call  ff(y,t,f,n,t0)
  k1=dt*f
  temp=y+0.5d0*k1
  call  ff(temp,t+0.5d0*dt,f,n,t0)
  
  k2=dt*f
  temp=y+0.5d0*k2
  call  ff(temp,t+0.5d0*dt,f,n,t0)
  
  k3=dt*f
  temp=y+k3
  call  ff(temp,t+dt,f,n,t0)
  
  k4=dt*f
  y(1:n)=y(1:n)+(k1(1:n)+2.d0*k2(1:n)+2.d0*k3(1:n)+k4(1:n))/6.d0
  
end subroutine rk4_f90

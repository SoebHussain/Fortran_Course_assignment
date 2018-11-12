program QDbySplit

  implicit none
  integer,parameter:: nx=400
  integer i,j
  real*8 dx,xst,xend,x(nx),dt,energy
  complex*8 psi(nx),npsi(nx),ev(nx),v(nx),g(nx),gd(nx)
  
  
  dt =0.1
  xst=-10
  xend=10
  do i=1,100
  call grid(xst,xend,nx,dx,x)
  call calpsi(psi,x,nx)
  call normalize(psi,nx,dx,npsi)
  call assignV(x,v,nx)
  call expov(x,v,nx,ev,dt)
  call gx(ev,npsi,g,nx)

  call differentiation_fourier(g,nx,gd,x)
  call gx(gd,ev,g,nx)
  energy=0
 
  write(7,*) x(i) , abs(psi(i)) 
end do

  
end program QDbySplit
!**********************DIFFerenciation by fourier********************
  subroutine grid_preparef(qst,qend,q,n)
  implicit none
  integer n,i
  real*8  qst,qend,q(n),dq
  dq=(qend-qst)/n
  do i=1,n
     q(i)=qst+(i-1)*dq
  enddo
end subroutine grid_preparef


subroutine fourierf(q1,q2,q1g,q2g,n,sign)
  implicit none
  integer n,i1,i2,sign
  complex*8 q1(n),q2(n),integrand(n),zi,out
  real*8 q1g(n),q2g(n),pi
  pi=dacos(-1.d0)
  zi=(0,1)
  do i2=1,n

     integrand=q1(:)*exp(sign*zi*q1g(:)*q2g(i2))
      call integral_cf(integrand,n,q1g(2)-q1g(1),out)
      q2(i2)=out/sqrt(2*pi)
   enddo
   end subroutine fourierf

subroutine integral_rf(integrand,n,dx,out)
  implicit none
  integer n,i
  real*8 integrand(n),dx,out
  out=0
  do i=1,n
     out=out+integrand(i)*dx
  enddo
  
end subroutine integral_rf

subroutine integral_cf(integrand,n,dx,out)
  implicit none
  integer n,i
  real*8   dx
  complex*8 integrand(n),out

  out=0
  do i=1,n
     out=out+integrand(i)*dx
  enddo
  
end subroutine integral_cf



subroutine differentiation_fourier(f,n,fd,x)
  implicit none
  integer n,order,i
  complex*8 f(n),fd(n),g(n),zi,a
  real*8 x(n),k(n),pi,dx,dt
  dt=0.1
  zi=(0,1)
    a= zi * dt
  dx=x(2)-x(1)
  pi=dacos(-1.d0)
   call grid_preparef(-pi/dx,pi/dx,k,n)
  
   call fourierf(f,g,x,k,n,-1)
   g=g*exp(a*((k*zi)**2))
  call fourierf(g,fd,k,x,n,+1)
!do i=1,n ; write(8,*)x(i),real(f(i)) ; enddo
!   do i=1,n ; write(9,*)k(i),real(g(i)) ; enddo   
end subroutine differentiation_fourier
!**********************gx=ev*psi*************************************
subroutine gx(ev,psi,g,n)
  implicit none
  integer n,i
  complex*8 ev(n),g(n),psi(n)
  do i=1,n
     g(i)=psi(i)*ev(i)
   !  print*,ev(i),psi(i),g(i)
  end do
  
end subroutine gx

!**********************expo(V)***************************************
subroutine expov(x,v,nx,ev,dt)
  implicit none
  integer nx,i
  complex*8 v(nx),ev(nx),a,zi
  real*8 x(nx),dt
  zi = (0,-1)
  a= zi * dt
  do i=1,nx
     ev(i) = 1! exp(a*v(i)/2)
   !  print*,ev(i)
  end do
end subroutine expov
!*********************ASSIGNING V*************************************
subroutine assignV(x,v,n)
  implicit none
  integer i,n
  complex*8 v(n)
  real*8 x(n)

  do i=1,n
     v(i) = -(x(i)**2)/2
    ! print*,v(i)
  end do
end subroutine assignV
!*********************NORMALISED PSI *********************************
subroutine normalize(psi,n1,dx,npsi)
  implicit none
  integer i,j,n,n1
  complex*8 psi(n1),npsi(n1)
   real*8 ni,sum1,dx
   do i=1,n1
        sum1=sum1+(psi(i)*conjg(psi(i))*dx)
     end do
     ni=sqrt(1/sum1)
  
  do i=1,n1
     npsi(i)=ni*psi(i)
    ! print*,npsi(i)
  end do
end subroutine normalize
!**********************CALCULATION OF PSI******************************
subroutine calPSI(psi,x,n)
  implicit none
  integer i,n
  real*8 x(n)
  complex*8 psi(n),zi
  zi=(0,1)
  do i=1,n
     psi(i) = exp(((x(i)**2)/(-2)) + (zi*x(i))  )
    ! print*,psi(i)
  end do
end subroutine calPSI
!**************************GRID PREPARE*****************************
subroutine grid(xi,xf,n,dx,x)
  implicit none
  integer::i,n
  real*8::xi,xf,dx,x(n)
     dx=(xf-xi)/n
  
 do i=1,n
    x(i)=xi+(i-1)*dx
 end do
end subroutine grid
  

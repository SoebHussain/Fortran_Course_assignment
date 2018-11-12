program main
  implicit none
  integer,parameter::nx=400,nt=200,np=21
  integer i,j,k
  real*8 x1,x2,dx,t1,t2,dt,x(nx),t(nt),phi(np,nx),V(nx),fx(nx),ft(nt),her(np,nx)
  
  complex*16 zi, expH(nt),expX(nt),psi(nt,nx),g(nx),fn(nx),psi0(nx),expH2(nt),out
  zi=(0.d0,1.d0)
  x1=-10.0; x2=30.0; t1=0.0; t2=20.0
  dx=(x2-x1)/nx; dt=(t2-t1)/nt
  
  call creategrid(x1,x2,nx,x)
  call creategrid(t1,t2,nt,t)
  call hermite(x1,x2,nx,np-1,her)

  do i=1,np
     phi(i,:)=her(i,:)*exp(-0.5*x**2)
  end do

  psi0=exp(-0.5*(x-2)**2)
  V=x
  fx=0.0
  ft=0.0
  call QdpBasis(x1,x2,t2,nx,nt,np,phi,psi0,V,fx,ft,psi,expH)

  
  do i=1,nt
     g=psi(i,:)
     fn=conjg(g)*x*g
     expX(i)=sum(fn)*dx
  end do

   do i=1,nt           ! calculation of exp. of H
     fn=psi(i,:)
     call compsciexp(x1,x2,nx,conjg(fn),fn,V,out)
     expH2(i)=out
  end do
  
  do i=1,nt
     write(6,*)t(i),expH(i),expH2(i)
     write(9,*)t(i),real(expX(i))
 end do
  
end program main





subroutine QdpBasis(x1,x2,t2,nx,nt,np,phi,psi0,V,fx,ft,psi,expH0)   ! phi is real
  implicit none
   integer i,j,k,nx,nt,np
   real*8 x1,x2,x(nx),dx,pi,H0mat(np,np),Hfmat(np,np),integrand(nx),integral
   
   real*8 t1,t2,t(nt),dt,V(nx),fx(nx),ft(nt),phi(np,nx),dphi(nx)
  
   complex*16 U(np,np),zi,C(np),Hmat(np,np),expH0(nt),psi(nt,nx),psi0Integrand(nx),psi0Int,psi0(nx)
  
   zi=(0.d0,1.d0)
   dx=(x2-x1)/nx 
   pi=dacos(-1.d0)
     
   t1=0.0; dt=(t2-t1)/nt  
   call creategrid(t1,t2,nt,t)
   call creategrid(x1,x2,nx,x)

   do i=1,np                  ! normalisation of phi
      integrand=abs(phi(i,:)**2)
      integral=sum(integrand)*dx
      phi(i,:)=phi(i,:)/sqrt(integral)
   end do
   
   integrand=abs(psi0**2)   ! normalisation of psi0
   integral=sum(integrand)*dx
   psi0=psi0/sqrt(integral)
   
   do j=1,np                  ! calculation of H mat
      call fourdiff(x1,x2,nx,phi(j,:),dphi)
      do i=1,np
         integrand=-0.5*phi(i,:)*dphi+phi(i,:)*V*phi(j,:)
         integral=sum(integrand)*dx
         H0mat(i,j)=integral

         integrand=phi(i,:)*phi(j,:)*fx
         integral=sum(integrand)*dx
         Hfmat(i,j)=integral       
      end do
   end do

   do i=1,np            ! calculation of C at t=0
      psi0Integrand=phi(i,:)*psi0
      psi0Int=sum(psi0Integrand)*dx
      C(i)=psi0Int
   end do

   do i=1,nt
      Hmat=-zi*dt*(H0mat+ft(i)*Hfmat)
      call expmat(Hmat,np,U)      ! calculation of U
      do j=1,nx
         psi(i,j)=sum(c*phi(:,j))     ! calculation of psi
      end do     
      expH0(i)=sum(conjg(c)*matmul(H0mat,c))          ! calculation of exp. of H=T+V
      C=matmul(U,C)
   end do
   
 end subroutine QdpBasis
 

!****************************************************************

! subroutine=1: take a real matrix from user
subroutine takemat(n1,n2,x)
  implicit none
  integer i,j,n1,n2
  real*8 x(n1,n2)
  
  print*,"enter the element of a matrix ",n1,"X",n2
  do i=1,n1
     do j=1,n2
        read(5,*)x(i,j)
     end do
     print*,""
  end do
  print*,"entered matrix hase been stored " 
end subroutine takemat


!subroutiine=2: exp(A) A is any complex square matrix
subroutine expmat(A,n,out)
  implicit none 
  integer i,j,n
  complex*16  A(n,n),out(n,n),temp(n,n)
  temp=0
  do i=1,n
     temp(i,i)=1
  end do
  out=temp
  
  do i=1,500
     temp= matmul(temp,A)/i
     out=out+temp
  end do
  
end subroutine expmat


!subroutine=3-a: fourier differentiation of real function   
subroutine fourdiff(x1,x2,n,f,ff)
  implicit none
  integer i,n
  real*8 x(n),x1,x2,dx,k(n),k1,k2,dk,pi,f(n),ff(n)
  complex*16 g(n),zi,integral,fn(n)
 
  zi = (0.d0,1.d0)
  pi = dacos(-1.d0)
  dx=(x2-x1)/n
  k1=(-pi)/dx; k2=pi/dx; dk=(k2-k1)/n
  call creategrid(x1,x2,n,x)
  call creategrid(k1,k2,n,k)
  
  do i=1,n
     fn=f*exp(-zi*k(i)*x)
     integral=sum(fn*dx)
     g(i)=integral/sqrt(2*pi)
  end do
  
  g=g*(zi*k)**2  ! change it only for higher differentiation
  
 do i=1,n
    fn=g*exp(zi*k*x(i))
    integral=sum(fn*dk)
    ff(i)=real(integral/sqrt(2*pi))
 end do
 
end subroutine fourdiff

!subroutine=3-b: to calculate double differentiation of any coplex number
subroutine compfourdiff(x1,x2,n,f,ff)
  implicit none
  integer i,n
  real*8 x(n),x1,x2,dx,k(n),k1,k2,dk,pi
  complex*16 g(n),zi,integral,fn(n),f(n),ff(n)
 
  zi = (0.d0,1.d0)
  pi = dacos(-1.d0)
  dx=(x2-x1)/n
  k1=(-pi)/dx; k2=pi/dx; dk=(k2-k1)/n
  call creategrid(x1,x2,n,x)
  call creategrid(k1,k2,n,k)
  
  do i=1,n
     fn=f*exp(-zi*k(i)*x)
     integral=sum(fn*dx)
     g(i)=integral/sqrt(2*pi)
  end do
  
  g=g*(zi*k)**2  ! change it only for higher differentiation
  
 do i=1,n
    fn=g*exp(zi*k*x(i))
    integral=sum(fn*dk)
    ff(i)=integral/sqrt(2*pi)
 end do
 
end subroutine compfourdiff

!subroutine=4: to calculate hermite polynomials
subroutine hermite(x1,x2,nx,np,h) 
   implicit none
   integer nx,np,i,j
   real*8  x1,x2,dx,x(nx),h(0:np,nx)
   call creategrid(x1,x2,nx,x)
   h(0,:)=1
   h(1,:)=2*x

  do j=1,np-1
     h(j+1,:)=(2*x*h(j,:))-2*j*h(j-1,:)
  end do
  
end subroutine hermite


! subroutine=5-a: to calculate expectation value of T+V(x) for any real sci
subroutine sciexp(x1,x2,n,sci1,sci2,V,out)
  implicit none
  integer n,i
  real*8 x1,x2,sci1(n),sci2(n),dsci2(n),V(n),h(n),out,dx
  dx=(x2-x1)/n
  call fourdiff(x1,x2,n,sci2,h)
  dsci2=-0.5*h
  h=sci1*dsci2+sci1*V*sci2
  out=sum(h*dx)
  
end subroutine sciexp


!subroutine=5-b:to calculate the expectation value T+V(x) for any complex sci
subroutine compsciexp(x1,x2,n,sci1,sci2,V,out)
  implicit none
  integer n,i
  real*8 x1,x2,V(n),dx
  complex*16 sci1(n),sci2(n),dsci2(n),h(n),out
  
  dx=(x2-x1)/n
  call compfourdiff(x1,x2,n,sci2,h)
  dsci2=-0.5*h
  h=sci1*dsci2+sci1*V*sci2
  out=sum(h*dx)
  
end subroutine compsciexp


!subroutine=6: create grid
subroutine creategrid(q1,q2,n,q)
  implicit none
  integer i,n
  real*8 q1,q2,q(n),dq
  dq=(q2-q1)/n
  do i=1,n
     q(i)=q1+(i-1)*dq
  end do
end subroutine creategrid

!subroutine=7: to calculate the grid number for any x
subroutine gridnum(x1,x2,n,gn)
  implicit none
  integer n,gn
  real*8 x1,x2,x
!  print*,"enter any x b/w ",x1,"and",x2,"at which you want to calculate the value"
  read(5,*)x
  gn=(x-x1)*n/(x2-x1)+1
end subroutine gridnum

  


!############   LESS USE #########################


!subroutine=2: real matrix multiplication
subroutine matmult(a,b,c,n1,n2,n3)
  implicit none
  integer i,j,k, n1,n2,n3
  real*8 a(n1,n2),b(n2,n3),c(n1,n3)
  do i=1,n1
     do j=1,n3
        c(i,j)=0.0
        do k=1,n2
           c(i,j)=c(i,j)+a(i,k)*b(k,j)
        end do
     end do
  end do
end subroutine matmult


!subroutine=7-a: integration of any real valued function
subroutine integration(l1,l2,n,integrand,out)
  implicit none
  integer j,n
  real*8 l1,l2,dl,out,integrand(n)
  dl=(l2-l1)/n
  out=0.0
  do j=1,n
     out=out+integrand(j)*dl
  end do    
end subroutine integration

!subroutine=7-b: integration of any complex  function  
subroutine compintegration(l1,l2,n,integrand,out) 
  implicit none
  integer j,n
  real*8 l1,l2,dl
  complex*16 out,integrand(n)
  dl=(l2-l1)/n
  out=0.0
  do j=1,n
     out=out+integrand(j)*dl
  end do    
end subroutine compintegration

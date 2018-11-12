program pr18
 implicit none
 integer*8,parameter::n=400,n1=25,lwork=10*n1,n2=1
 real*8,parameter::pi=dacos(-1.d0)
 integer*8::i,j,info,k
 real*8::phi(n1,n),h(n1,n),xi,xf,dx,x(n),rwork(lwork),w(n1),sum1,dt,ti,tf,t(n),l,ev(n1,n1),work(lwork)
 complex*16::zi,a(n1,n1),a1(n),u(n1,n1),c(n1,n),xt(n),diag(n1,n1),psi(n),sum,a2(n)


 zi=(0,1)
 xi=-10.0;xf=10.0;l=20.0
 call grid(xi,xf,n,dx,x)
 !call hermite(xi,xf,n1,n,h)
 do i=1,n1
    do j=1,n
       phi(i,j)=sqrt(2/l)*Sin(i*pi*(x(j)-xi)/l)
    end do
 end do
 call normalize(phi,n1,n,dx,phi)

 call kengfun(xi,xf,n1,n,phi,a) !creating h(i,j) matrix in a



 do i=1,n1
     do j=1,n1
     if(abs(a(i,j)-a(j,i)).gt.1e-10) write(6,*)a(i,j),a(j,i)
        if(i==j) then
       ev(i,j)=real(a(i,j))+(i*i*pi*pi)/(2*l*l)
     else
       ev(i,j)=real(a(i,j))
     end if
     end do
  end do





call dsyev('V','L',n1,ev,n1,w,work,lwork,info)
do i=1,n1
   write(6,*) w(i)
end do



           do j=1,n
              sum1=0
                do i=1,n1
                   sum1=sum1+ev(i,1)*phi(i,j)
                end do
                psi(j)=sum1
           end do
           do j=1,n
            write(7,*) x(j),real(psi(j))
          end do   
         
end program pr18


!************************************subrotuine

!subroutine for hermite
subroutine hermite(xi,xf,n1,n,h)
  implicit none
  integer*8::i,j,n1,n
  real*8::xi,xf,dx,x(n),h(n1,n)
  call grid(xi,xf,n,dx,x)
  do j=1,n
     h(1,j)=1
     h(2,j)=2*x(j)
  end do
  do i=3,n1
     do j=1,n
        h(i,j)=2*x(j)*h(i-1,j)-2*(i-2)*h(i-2,j) !recursion relation
     end do
  end do
end subroutine hermite

!subroutine for normalize
subroutine normalize(psi,n1,n,dx,npsi)
  implicit none
   integer*8::i,j,n,n1
   real*8::psi(n1,n),ni(n1),sum1,dx,npsi(n1,n)
   do i=1,n1
     sum1=0
     do j=1,n
        sum1=sum1+(psi(i,j))**2*dx
     end do
     ni(i)=sqrt(1/sum1)
  end do
 
   do i=1,n1
     do j=1,n
        npsi(i,j)=ni(i)*psi(i,j)
     end do
  end do
end subroutine normalize

!subroutine for grid
 subroutine grid(xi,xf,n,dx,x)
   implicit none
    integer*8::i,n
     real*8::xi,xf,dx,x(n)
      dx=(xf-xi)/n
  do i=1,n
     x(i)=xi+(i-1)*dx
  end do
end subroutine grid

!subroutine for psi_h_psi
subroutine kengfun(xi,xf,n1,n,psi,a)
 implicit none
  integer*8::i,m,j,k,n,n1
  real*8::x(n),xi,xf,dx,psi(n1,n),al
  complex*16::dfm(n1,n),psi1(n1,n),sum,g(n),a(n1,n1)
  m=2
  call grid(xi,xf,n,dx,x)
  call diffour(xi,xf,n1,n,m,psi,dfm)
  do i=1,n1
     do j=1,n
        psi1(i,j)=-(dfm(i,j)/2)+(x(j)**2*psi(i,j)/2)
     end do
      do k=1,n1
          do j=1,n
             g(j)=psi(k,j)*psi1(i,j)
          end do
         call integratec(xi,xf,n,g,sum)
         a(i,k)=sum
      end do
   end do
 end subroutine kengfun

 subroutine diffour(xi,xf,n1,n,m,f,dfm)
  implicit none
  real*8,parameter::pi=dacos(-1.d0)
  integer*8::i,m,j,n,n1,l
  real*8::x(n),xi,xf,dx,k(n),dk,ki,kf,h,f(n1,n)
  complex*16::zi,g(n1,n),sum1(n),dfm(n1,n),f1(n1,n)

  zi=(0,1)
  !grid point
  call grid(xi,xf,n,dx,x)
  ki=-pi/dx ;kf=pi/dx
  call grid(ki,kf,n,dk,k)
  !integration of gk
  do l=1,n1
     do i=1,n
         sum1(i)=0
       do j=1,n
          sum1(i)=sum1(i)+f(l,j)*exp(-zi*x(j)*k(i))*dx
       end do
        g(l,i)=sum1(i)/sqrt(2*pi) 
    end do
 end do
 !inverse integration
 do l=1,n1
  do i=1,n
     sum1(i)=0
     do j=1,n
         sum1(i)=sum1(i)+(-zi*k(j))**m*g(l,j)*exp(zi*x(i)*k(j))*dk
         f1(l,i)=sum1(i)/sqrt(2*pi) 
     end do
      dfm(l,i)=f1(l,i)
  end do
 end do
  return  
end subroutine diffour

!integration of complex function
subroutine integratec(xi,xf,n,f,sum)
  implicit none
  integer*8::i,n
  real*8::xi,xf,dx,x(n)
  complex*16::f(n),sum
  call grid(xi,xf,n,dx,x)
  sum=0
  do i=1,n
     sum=sum+f(i)*dx
  end do
end subroutine integratec

subroutine kengfunc(xi,xf,n1,n,psi,a)
 implicit none
  integer*8::i,m,j,k,n,n1
  real*8::x(n),xi,xf,dx,al
  complex*16::dfm(n1,n),psi1(n1,n),sum,g(n),a(n1),psi(n1,n)
  m=2
  call grid(xi,xf,n,dx,x)

  call diffourc(xi,xf,n1,n,m,psi,dfm)
  do i=1,n1
     sum=0
     do j=1,n
        psi1(i,j)=-(dfm(i,j)/2)+(x(j)**2*psi(i,j)/2)
     end do
    
          do j=1,n
             g(j)=conjg(psi(i,j))*psi1(i,j)
          end do
          call integratec(xi,xf,n,g,sum)
          a(i)=sum
      end do
 end subroutine kengfunc

 subroutine diffourc(xi,xf,n1,n,m,f,dfm)
  implicit none
  real*8,parameter::pi=dacos(-1.d0)
  integer*8::i,m,j,n,n1,l
  real*8::x(n),xi,xf,dx,k(n),dk,ki,kf,h
  complex*16::zi,g(n1,n),sum1(n),dfm(n1,n),f1(n1,n),f(n1,n)

  zi=(0,1)
  !grid point
  call grid(xi,xf,n,dx,x)
  ki=-pi/dx ;kf=pi/dx
  call grid(ki,kf,n,dk,k)
  !integration of gk
  do l=1,n1
     do i=1,n
         sum1(i)=0
       do j=1,n
          sum1(i)=sum1(i)+f(l,j)*exp(-zi*x(j)*k(i))*dx
       end do
        g(l,i)=sum1(i)/sqrt(2*pi) 
    end do
 end do
 !inverse integration
 do l=1,n1
  do i=1,n
     sum1(i)=0
     do j=1,n
         sum1(i)=sum1(i)+(-zi*k(j))**m*g(l,j)*exp(zi*x(i)*k(j))*dk
         f1(l,i)=sum1(i)/sqrt(2*pi) 
     end do
      dfm(l,i)=f1(l,i)
  end do
 end do
  return  
end subroutine diffourc

subroutine kengfunc1(xi,xf,n1,n,psi,a)
 implicit none
  integer*8::i,m,j,k,n,n1
  real*8::x(n),xi,xf,dx,al
  complex*16::dfm(n1,n),psi1(n1,n),sum,g(n),a(n1),psi(n1,n),psi2(n1,n)
  m=2
  call grid(xi,xf,n,dx,x)

  call diffourc(xi,xf,n1,n,m,psi,dfm)
  do i=1,n1
    
   do j=1,n
        psi1(i,j)=-(dfm(i,j)/2)+(x(j)**2*psi(i,j)/2)
     end do
  end do
  
  call diffourc(xi,xf,n1,n,m,psi1,dfm)
    do i=1,n1
      sum=0
     do j=1,n
        psi2(i,j)=-(dfm(i,j)/2)+(x(j)**2*psi1(i,j)/2)
     end do
    
          do j=1,n
             g(j)=conjg(psi1(i,j))*psi2(i,j)
          end do
          call integratec(xi,xf,n,g,sum)
          a(i)=sum
      end do
 end subroutine kengfunc1
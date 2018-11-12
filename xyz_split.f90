
program asd
  implicit none 
  integer ,parameter :: n=400
  real x(n),dx,xi,xf,dt,pi,time,temp
  complex ev(n),zi,a,fev(n),ff(n),fevg(n),f(n),V(n),answer
  integer i,p,j
  pi=dacos(-1.d0)
  dt=0.1
  f=0
  zi=(0,1)
  a=-zi*dt
 ! print*,zi
  
  xf=30
  xi=-10
  dx= (xf-xi)/n

  do i=1,n
  x(i) = xi + (i-1)*dx
  f(i) =exp((( (x(i)**2)/(-2))) + zi*x(i) )           !!! defined initial psi
end do

call noram(dx,f,n)                                     ! normalised psi
!CALL ier(f,n,v,x,2)  ! 


time = 20           !print*,f(1:10)
p = time/dt
do j=1,p

                  !   print*,f(1:10)
do i=1,n
   ev(i) =1 !exp( zi * dt * (x(i) **2) /(- 4))
end do
                          ! print*,ev(1:10)

call fevR(ev,f,fev,n)
                        !!print*,fev(1:10) 


call  differentiation_fourier(fev,n,ff,x,1)
!print*, ff(1:10)

call fevC(ev,ff,fevg,n)


f = fevg

call interg(f,answer,dx,v,n)
temp = j*dt

write(17,*),temp,answer


end do








end program asd


 
subroutine noram(dx,f,n)
  implicit none
  real dx
  complex f(n),g(n),sas
  integer i,n
  sas=0
  
   do i=1,n
      sas = sas + (f(i)*conjg(f(i))*dx)
   end do

   sas = abs(sas)
   sas = 1/sqrt(sas)
     ! print*,sas

   f = f*sas
  
 end subroutine noram

 subroutine fevR(ev,f,fev,n)
   complex ev(n),fev(n),f(n)
   integer n

   do i=1,n
      fev(i) = f(i) * ev(i)
   end do
 end subroutine fevR
 subroutine fevC(ev,f,fev,n)
   complex ev(n),fev(n)
   complex f(n)
   integer n

   do i=1,n
      fev(i) = f(i) * ev(i)
   end do
 end subroutine fevC

 
 subroutine difn_fourier(f,n,fd,x,order)
  implicit none
  integer n,order,i
  complex f(n),fd(n),g(n),zi
  real x(n),k(n),pi,dx
  zi=(0,1)
  dx=x(2)-x(1)
  pi=dacos(-1.d0)
  call grid_prepare(-pi/dx,pi/dx,k,n)
  
   call fourier(f,g,x,k,n,1)
  
    g=g*(zi*k)**order
    call fourier(g,fd,k,x,n,-1)
    
    
   
  end subroutine difn_fourier
  
  subroutine ier(f,n,fd,x,order)
  implicit none
  integer n,order,i
  complex f(n),fd(n),g(n),zi
  real x(n),k(n),pi,dx
  zi=(0,1)
  dx=x(2)-x(1)
  pi=dacos(-1.d0)
   call grid_prepare(-pi/dx,pi/dx,k,n)
  
   call fourier(f,g,x,k,n,-1)
   g=g*(zi*k)**order
  call fourier(g,fd,k,x,n,+1)
!do i=1,n ; write(8,*)x(i),real(f(i)) ; enddo
  ! do i=1,n ; write(9,*)k(i),real(g(i)) ; enddo   
end subroutine ier


   


subroutine differentiation_fourier(f,n,fd,x,order)
  implicit none
  integer n,order,i
  complex f(n),fd(n),g(n),zi
  real x(n),k(n),pi,dx
  zi=(0,1)
  dx=x(2)-x(1)
  pi=dacos(-1.d0)
  call grid_prepare(-pi/dx,pi/dx,k,n)
  
   call fourier(f,g,x,k,n,-1)
   g=g*exp(((zi*k)**2)*zi*0.05)
  call fourier(g,fd,k,x,n,+1)
   
    end subroutine differentiation_fourier

subroutine grid_prepare(qst,qend,q,n)
  implicit none
  integer n,i
  real qst,qend,q(n),dq
  dq=(qend-qst)/n
  do i=1,n
     q(i)=qst+(i-1)*dq
  enddo
end subroutine grid_prepare

subroutine interg(fevg,answer,dx,v,n)
  implicit none
  complex fevg(n),answer,v(n)
  
  real dx!,v(n)
  integer n,i
   answer =0
  do i=1,n
     answer = answer + v(i)*fevg(i)*dx*conjg(fevg(i))
  end do

  !!write(6,*)
end subroutine interg



   subroutine integral_c(integrand,n,dx,out)
  implicit none
  integer n,i
  real  dx
  complex integrand(n),out

  out=0
  do i=1,n
     out=out+integrand(i)*dx
    
  enddo
  
end subroutine integral_c

  subroutine fourier(q1,q2,q1g,q2g,n,sign)
  implicit none
  integer n,i1,i2,sign
  complex q1(n),q2(n),integrand(n),zi,out
  real q1g(n),q2g(n),pi
  pi=dacos(-1.d0)
  zi=(0,1)
  do i2=1,n

     integrand(:)=q1(:)*exp(sign*zi*q1g(:)*q2g(i2))
   !! print*,integrand
      call integral_c(integrand,n,q1g(2)-q1g(1),out)
      q2(i2)=out/sqrt(2*pi)
   enddo
   end subroutine fourier
  

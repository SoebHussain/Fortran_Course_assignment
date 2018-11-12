
Program main 
  implicit none
  integer,parameter:: n=5,nw=n
  real*8 y(n),dt,t,tend,tol ,  w(nw,9), c(50+n)
  integer  ierr,ind
  integer i
  external ff
  ind=1
  y=0.d0
  
  y(1)=0!2.0d8  
  y(2)=1.0d16
  y(3)=0!7.0d11
  y(4)=0
  y(5)=0
  dt=0.1d0
  tol=1d-8
t=0.d0
tend=0.d0


  do i=0,10**8*2
!dt=dt+0.00000001
     t=tend
     tend=t+dt
!     t=(i-1)*dt
     !call rk1(y,t,dt,n)
     ! call rk2(y,t,dt,n)
     call rk4_f90(y,t,dt,n)
   !  call dverk_old_modified (n,ff,t,y,tend,tol,ind,c,nw,w)
if(mod(i,1000).eq.0)     write(8,*)t,y(1),y(2),y(3)

enddo
tend=0
y(4)=1.0d8
y(5)=5.0d4
  do i=0,10**8
!dt=dt+0.00000001
     t=tend
     tend=t+dt
!     t=(i-1)*dt
     !call rk1(y,t,dt,n)
     ! call rk2(y,t,dt,n)
     call rk4_f90(y,t,dt,n)
   !  call dverk_old_modified (n,ff,t,y,tend,tol,ind,c,nw,w)
if(mod(i,1000).eq.0)     write(9,*)t,y(1),y(2),y(3)

  enddo
end Program main




!subroutine ff(f,x,t,nvar)
subroutine ff(x,t,f,nvar) !with rk4
!subroutine ff(nvar,t,x,f) !with dverk

  implicit none
  integer nvar
  double precision f(nvar),x(nvar),t,k,m
  double precision ak1,ak2,ak3,ak4,c1,c2,c3,c4,k5,k6
  !c4=9.0d17
  ak1=3.0d-12	
  ak2=1.22d-33
  ak3=5.5d-4
  ak4=6.9d-16
  k5=4.82d-11
  k6=10.39d-12
  c1=x(1)
  c2=x(2)
  c3=x(3)
  c4=c2/.22d0
  

  f(1)=-(ak2)*(c4)*(c1)*(c2)+(ak3)*(c3)-(ak4)*(c1)*(c3)+2*ak1*(c2)      &
	 -k5*x(1)*x(4)
  f(2)=-ak1*(c2)-ak2*(c4)*(c1)*c2+ak3*(c3)+2*(ak4)*c1*c3    		&
	+k5*x(1)*x(4)+k6*x(5)*x(3)
  !f(2)=0.d0
  f(3)=ak2*(c4)*(c1)*(c2)-ak3*c3-ak4*c1*c3     			&
	-k6*x(5)*x(3)
  f(4)=-k5*x(1)*x(4)+k6*x(5)*x(3)
  f(5)=+k5*x(1)*x(4)-k6*x(5)*x(3)
  !      write(6,*)c1,c2,c3,c4,ak1,ak2,ak3,ak4,f(1),f(2),f(3)
  !stop
end subroutine ff


subroutine rk4_f90(y,t,dt,n)
  implicit none
  integer n,i
  real*8 y(n),f(n),t,k,dt,temp(n),k1(n),k2(n),k3(n),k4(n)
  call  ff(y,t,f,n)
  k1=dt*f
  temp=y+0.5d0*k1
  call  ff(temp,t+0.5d0*dt,f,n)
  
  k2=dt*f
  temp=y+0.5d0*k2
  call  ff(temp,t+0.5d0*dt,f,n)
  
  k3=dt*f
  temp=y+k3
  call  ff(temp,t+dt,f,n)
  
  k4=dt*f
  y(1:n)=y(1:n)+(k1(1:n)+2.d0*k2(1:n)+2.d0*k3(1:n)+k4(1:n))/6.d0
  
end subroutine rk4_f90


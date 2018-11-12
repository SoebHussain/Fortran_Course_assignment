program split
  implicit none
  integer, parameter:: nx= 400
  
  integer i
  real*8 xst,xend,x(nx),dx,normal,psi(nx),integrand(nx),out,v(nx),dt
  complex*16 a,z,expo_v(nx)
  xst=-10
  xend=10
  dt =0.1
  z= (0,1)
  a = -z*dt
  call grid_prepare(xst,xend,x,nx,dx)
! write(6,*) dx
  do i=1,nx
     psi(i)= (1+x(i))*exp((x(i)**2)/(-2))
    ! write(6,*) psi(i)
  end do

  do i=1,nx
     integrand(i) = psi(i)*psi(i)
  end do

  !  write(6,*) integrand (250)

   call integral_r(integrand,nx,dx,out)
   normal = 1/sqrt(out)
   ! write(6,*) normal
   do i= 1,nx
      psi(i)= psi(i)/normal
     ! write(6,*) psi(i)
   end do

   do i=1,nx
      v(i) = x(i)**2/2
      !write(6,*) v(i)
   end do

   do i=1,10
      expo_v(i) = exp((a*v(i))/2)
      write(6,*) expo_v(i)
   end do

   
   

  
 
 
end program split

!*********************************

subroutine grid_prepare(qst,qend,q,n,dq)
  implicit none
  integer n,i
  real*8  qst,qend,q(n),dq
  dq=(qend-qst)/n
  do i=1,n
     q(i)=qst+(i-1)*dq
  enddo
end subroutine grid_prepare
!********************************

subroutine integral_r(integrand,n,dx,out)
  implicit none
  integer n,i
  real*8 integrand(n),dx,out
  ! print*,dx
  out=0
  do i=1,n
     out=out+integrand(i)*dx
    
  end do

 ! write(6,* ) out
  
end subroutine integral_r
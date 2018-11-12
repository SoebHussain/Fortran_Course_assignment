program rk
  implicit none

  integer,parameter :: n=5
  real*8 y(n),t,h,yf(n)
  integer i

  h=0.1
  y(1)=207683979.20418355  !!o
  y(2)=9999685902506188.0   !! o2
  y(3)=209332675446.10599      !!ozone
  y(4)=1e8  !! clo
  y(5)=5*1e4 !!cl 
  
  do i=1,1000*200000
     t=i*h
     call rk4(t,y,h,n,yf)
     y(1:n)=yf(1:n)  
    if(mod(i,1000).eq.0)  write(11,*)t,yf(1),yf(2),yf(3),yf(4),yf(5)
  end do

end program rk

!******************************************************************************

  subroutine rk4(t,y,h,n,yf)
    implicit none

    real*8 h,t,y(n),yf(n),yy(n)
    real*8 k1(n),k2(n),k3(n),k4(n),temp1(n),temp2(n),temp3(n)
    integer i,n

    do i=1,n
       call fnctn(t,y,yy)
       k1(i)=h*yy(i)
       temp1(i)=y(i)+0.5*k1(i)
    end do

    do i=1,n
       call fnctn(t+0.5*h,temp1,yy)
       k2(i)=h*yy(i)
       temp2(i)= y(i)+0.5*k2(i)
    end do

    do i=1,n
       call fnctn(t+0.5*h,temp2,yy)
        k3(i)=h*yy(i)
       temp3(i)= y(i)+k3(i)
    end do

    do i=1,n
       call  fnctn(t+h,temp3,yy)
       k4(i)=h*yy(i)
       yf(i)=y(i)+(k1(i)+2*k2(i)+2*k3(i)+k4(i))/6.d0
    end do

  end subroutine rk4

  
!********************************************************************
  subroutine fnctn(t,temp,ans)
    implicit none
    real*8 t,temp(5),ans(5),k1,k2,k3,k4,k5,k6,M

   
    k1=3.0*1e-12
    k2=1.2*1e-33
    k3=5.5*1e-4
    k4=6.9*1e-16
    k5=4.82*1e-11
    k6=10.39*1e-12
    M=temp(2)/.22d0

       ans(1)=2*k1*temp(2)-k2*M*temp(1)*temp(2)+k3*temp(3)-k4*temp(1)*temp(3)-k5*temp(1)*temp(4)!!o
       ans(2)=-temp(2)*k1-k2*M*temp(1)*temp(2)+k3*temp(3)+2*k4*temp(1)*temp(3)+k5*temp(1)*temp(4)+k6*temp(3)*temp(5)!!02
       ans(3)=k2*M*temp(1)*temp(2)-k3*temp(3)-k4*temp(1)*temp(3)-k6*temp(5)*temp(3)!!o3
       ans(4)=-k5*temp(1)*temp(4)+k6*temp(5)*temp(3)!! clo
       ans(5)=k5*temp(1)*temp(4)-k6*temp(5)*temp(3)!! cl
       
     end subroutine fnctn


     

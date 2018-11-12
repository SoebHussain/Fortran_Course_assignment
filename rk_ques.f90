program ranga
  implicit none
  integer i
  integer, parameter :: n=4
  real*8 t,y(n),yf(n),h
  h = 0.1
  y(1)=2
  y(2)=2
  y(3)=2
  y(4)=2
 ! call rk4(t,y,h,n,yf)
  do i=1,10
     t=i*h
     call rk4(t,y,h,n,yf)
     y(1:n)=yf(1:n)  
     write(6,*)t,yf(1),yf(2),yf(3),yf(4)
  end do

  
end program ranga
!**************************************************************************************
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
  !**********************************************************************************
  
  subroutine fnctn(t,temp,ans)
    implicit none
    real*8 t,temp(4),ans(4)
    ans(1)= temp(3)
    ans(2)= temp(4)
    ans(3)=-temp(1)/sqrt((temp(1)**2 + temp(2)**2))
    ans(4)=-temp(2)/sqrt((temp(1)**2 + temp(2)**2))
     end subroutine fnctn
     
!***************************************************************
  
 

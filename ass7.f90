program ass7
  implicit none

  integer,parameter :: n=2
  real*8 y(n),x,h,yf(n),xi,xf,ye(n),yeuler(n),yh1(n),yh2(n),yheun(n),yavg(n),yfake(n)
  integer i

  xi=0
  xf=2
  h=(xf-xi)/25
  y(1)=0
  y(2)=0

  write(6,*)'RK4'
  do i=1,25
     x=i*h
     call rk4(x,y,h,n,yf)
     y(1:n)=yf(1:n)
     write(*,*)x,yf
  end do

  yeuler(1)=0
  yeuler(2)=0
  x=0

  write(6,*)'euler'
  do i=1,25
     call fnctn(x,yeuler,ye)
     yeuler(1:n)=yeuler(1:n)+ye(1:n)*h
     x=i*h
     write(6,*)x,yeuler
  end do

  yheun(1)=0
  yheun(2)=0
  x=0

  write(6,*)'heun'
  do i=1,25
     call fnctn(x,yheun,yh1)
     yfake(1:n)=yheun(1:n)+yh1(1:n)*h
     x=i*h
     call fnctn(x,yfake,yh2)
     yavg(1:n)=(yh1(1:n)+yh2(1:n))/2
     yheun(1:n)= yheun(1:n)+yavg(1:n)*h
     write(6,*)x,yheun
  end do
     
end program ass7



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

  subroutine fnctn(t,temp,ans)
    implicit none
    real*8 t,temp(2),ans(2)
       ans(1)=temp(1)*temp(2)+cos(t)-0.5*sin(2*t)
       ans(2)=temp(1)*temp(1)+ temp(2)*temp(2)-(1+sin(t))
     end subroutine fnctn
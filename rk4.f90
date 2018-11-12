!********************************
!name: rupal saxena
!rollno. 150122031
!***********************************
Program Ass7

  integer, parameter :: n = 10000

  real*8 X(n),Y(n),Z(n)
  real*8 xl,x0,h
  integer k,kl,l

  x0=0.d0       
  xl=2.d0       
  kl=100         
  h=(xl-x0)/kl  
  X(0)=x0
  Y(0)=0.d0
  Z(0)=0.d0
  
  do k=0, kl-1
    call RK4(X(k),Y(k),Z(k),h,X(k+1),Y(k+1),Z(k+1))
  end do
  
    do k=0, kl
    write(6,*) X(k), Y(k), Z(k)
  end do
 
END Program Ass7
!***********************************************************************************
real*8 Function F(x,y,z)
        real*8 x,y,z
        F=y*z+cos(x)-0.5d0*sin(2.d0*x)
        
      End Function F
      
real*8 Function G(x,y,z)
        real*8 x,y,z
        G=y*y+z*z-(1.d0+sin(x))
               
              End Function G
!************************************************************* subroutine for rk4            
Subroutine RK4(x,y,z,h,x1,y1,z1)
      real*8 x,y,z,h,x1,y1,z1
      real*8 c1,c2,c3,c4,d1,d2,d3,d4,h2,F,G
        c1=F(x,y,z)
        d1=G(x,y,z)
        h2=h/2.d0
        c2=F(x+h2,y+h2*c1,z+h2*d1)
        d2=G(x+h2,y+h2*c1,z+h2*d1)
        c3=F(x+h2,y+h2*c2,z+h2*d2)
        d3=G(x+h2,y+h2*c2,z+h2*d2)
        c4=F(x+h,y+h*c3,z+h*d3)
        d4=G(x+h,y+h*c3,z+h*d3)
        x1=x+h
        y1=y+h*(c1+2.d0*c2+2.d0*c3+c4)/6.d0
        z1=z+h*(d1+2.d0*d2+2.d0*d3+d4)/6.d0
              
              End Subroutine RK4
!****************************************************************************
              

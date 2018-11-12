program ass8
implicit none
real x,T,z,f1,f2,H0,Ta,x0,T0,xn,Tn,z0,g1,g2,g3,b,b1,b2,b3,e,h,shoot
integer p
p=0
H0=0.01
Ta =20
x0=0.0
T0=40.0
xn=10.0
Tn=200.0
h=1


write(6,*) "write your first guess"
read(*,*) g1
b=Tn
z0=g1
b1=shoot(x0,T0,z0,xn,h,1,h0,Ta)
write(6,*) "g1 is", b1
if(abs(b1-b)<0.0000001) then
write(6,*) "value of x and respective T is"
e=shoot(x0,T0,z0,xn,h,1,h0,Ta)
else
write(6,*) "Enter the next guess"
read(*,*) g2
z0=g2
b2= shoot(x0,T0,z0,xn,h,1,h0,Ta)
write(6,*) "g2 is", b2
end if

if(abs(b2-b)<0.0000001) then
write(6,*) "value of x and respective T is"
e=shoot(x0,T0,z0,xn,h,1,h0,Ta)
else 
write(6,*) "g1", g1,"and g2", g2
g3=g2+(((g2-g1)*(b-b2))/(1.0*(b2-b1)))
if(b1-b2==0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!
end program ass8

real function f1(x,T,z,Ta,h0)
real x,T,z
f1=z
end function

real function f2(x,T,z,Ta,H0)
real x,T,z,H0,Ta
f2=H0*(T-Ta)
end function

real function shoot(x0,T0,z0,xn,h,p,h0,Ta)
real x0,T0,z0,xn,h,x,T,z,kt1,kt2,kt3,kt4,kz1,kz2,kz3,kz4,kt,kz,x1,T1,z1,Ta,h0
integer p
x=x0
z=z0
T=T0
do while (x<xn)
	kt1=f1(x,T,z,Ta,H0)
        kz1=f2(x,T,z,Ta,H0)
        
        kt2=f1(x+h/2.0,T+kt1/2.0,z+kz1/2.0,Ta,H0)
        kz2=f2(x+h/2.0,T+kt1/2.0,z+kz1/2.0,Ta,H0)
        
        kt3=f1(x+h/2.0,T+kt2/2.0,z+kz2/2.0,Ta,H0)
        kz3=f2(x+h/2.0,T+kt2/2.0,z+kz2/2.0,Ta,H0)
        
        kt4=f1(x+h,T+kt3,z+kz3,Ta,H0)
        kz4=f2(x+h,T+kt3,z+kz3,Ta,H0)

        kz=1/6.0*(kz1+2*kz2+2*kz3+kz4)*h
        kt=1/6.0*(kt1+2*kt2+2*kt3+kt4)*h

        T1=T+kt
        x1=x+h
        z1=z+kz

	x=x1;
        T=T1;
        z=z1;
	if (p==1) then
	write(6,*) x,T
	end if
end do
shoot = T
end function


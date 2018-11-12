program hartree
  implicit none
  integer ,parameter :: nb=2
  integer i,j,k,l
  real*8 t(2,2),v1(2,2),v2(2,2),hcore(2,2),s(2,2),zita,a(22,22),c(2,2)
  real*8 temp(2,2),f(2,2),sum1(2,2),sout(2,2),p
  real*8 fp(2,2),cp(2,2)
  complex*16 z
  
  zita = 1.24 !slater function

  !define kinetic energy
  t(1,1) = 0.7600
  t(1,2) = 0.2365
  t(2,1) = 0.2365
  t(2,2) = 0.7600
  
  !potencial energy of first electron
  v1(1,1)=-1.2266
  v1(1,2)=-0.5974
  v1(2,1)=-0.5974
  v1(2,2)=-0.6538
  
  !potencial energy of second electron
  v2(1,1)=-0.6538
  v2(1,2)=-0.5974
  v2(2,1)=-0.5974
  v2(2,2)=-1.2266

  !core-harmitian 
  hcore(:,:) = t(:,:) + v1(:,:) + v2(:,:)
  !write(6,*) hcore

  
  s(1,1)=1.0
  s(1,2)=0.6593
  s(2,1)=0.6593
  s(2,2)=1.0
  p=-0.5
  call s_inv(s,p,sout) !calculation of s**(-1/2)
 
 
a=0    !(psi1psi2|....) 

a(11,11)=0.7746
a(22,22)=0.7746
a(11,22)=0.5697
a(22,11)=0.5697
a(21,11)=0.4441
a(12,11)=0.4441
a(11,21)=0.4441
a(11,12)=0.4441
a(22,21)=0.4441
a(22,12)=0.4441
a(21,22)=0.4441
a(12,22)=0.4441
a(21,21)=0.2970
a(12,12)=0.2970
a(12,21)=0.2970
a(21,12)=0.2970

c=0
do i=1,20
call funct(a,s,hcore,c,f) ! calculation of f for the first time
call fprime(sout,f,fp)   !calculation of f prime
call cprime(fp,cp)       !calculation of c prime
call citerate(cp,sout,c) ! calculation of c from c prime
end do

end program hartree


!********************************************************************************************************************************

subroutine funct(a,s,hcore,c,f)
  integer i,j,k,l
  real*8 sum1(2,2),temp(2,2),c(2,2),a(22,22),hcore(2,2),f(2,2),s(2,2)
  
sum1 =0
temp=0


do i=1,2
do j=1,2
do k=1,2
do l=1,2
   temp(i,j) = c(k,1)*c(l,1)*(2*a((10*i+j),(10*l+k))-a((10*i+k),(10*l+j)))
end do
end do
   sum1(i,j) = sum1(i,j) + temp(i,j)
end do
end do

!print*, sum1

f(:,:) = hcore(:,:) + sum1(:,:)
!print*, f
end subroutine funct

!**********************************************************************************************************************

subroutine s_inv(s,p,sout)
  implicit none
  integer,parameter :: nl =2, lwork = 10**nl
  integer info,i
  complex*16 work(lwork),z
  real*8 s(2,2), sout(2,2),p,w(nl),ww(nl,nl)
  z = (0,1)
  call dsyev('V','U',nl,s,nl,w,work,lwork,info)
  ww=0
  do i=1,nl
     ww(i,i) = w(i)**(p)
  end do

  
  sout=matmul(s,matmul(ww,transpose((s))))
 ! write(6,*) sout
end subroutine s_inv
!****************************************************************************************************************

subroutine fprime(sout,f,fp)
  implicit none
  integer i
  real*8 fp(2,2),f(2,2),sout(2,2)

  fp =matmul(sout, matmul(f,sout))

! write(6,*) fp
end subroutine fprime
!*****************************************************************************************************************

subroutine cprime(fp,cp)
  implicit none
  integer , parameter:: nl = 2, lwork= 10**nl
  integer info , i
  real*8 fp(2,2), cp(2,2),w(nl),ww(nl,nl)
  complex*16 work(lwork),z
  z = (0,1)
  call dsyev('V','U',nl,fp,nl,w,work,lwork,info)
  cp = fp
 ! write(6,*) cp
end subroutine cprime

!***************************************************************************************************************

subroutine citerate(cp,sout,c)
  implicit none
  real*8 cp(2,2),sout(2,2), c(2,2)
  c = matmul (sout,cp)
  write(6,*) c
end subroutine citerate

!**************************************************************************************************************

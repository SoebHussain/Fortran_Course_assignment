program lapack
  implicit none
  integer,parameter :: nl = 2, lwork= 10**nl
  integer info,i
  real*8 w(nl),ww(nl,nl),rwork(lwork)
  complex*16 h(2,2),z,work(lwork),f(nl,nl)
  z = (0,1)

  h(1,1)=0
  h(1,2)=z
  h(2,1)=-z
  h(2,2)=2
  !write(6,*)h
  call zheev ('V','U',nl,h,nl,w,work,lwork,rwork,info)

  ww = 0
  do i=1,nl
     ww(i,i) = sin(w(i))
  end do

  f= matmul(h,matmul(ww,transpose(conjg(h))))
  write(6,*) f
end program lapack


  

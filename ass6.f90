  !*******************
  !Name: Rupal Saxena
  !Roll number: 150122031
  !*******************
program assignment6
  implicit none
 integer i
  real f,y_eu(1000),y_heu(1000),h,x(1000),y_ex(1000),error_heu(1000),error_eu(1000)
 h= 0.5
 x(1) = 0
 y_ex(1) = 1
 y_eu(1)= 1
 y_heu(1)=1
 do i=1,8
    x(i+1)= x(i)+h
    y_ex(i+1)= -0.5*(x(i+1)**4)+4*(x(i+1)**3)-10*(x(i+1)**2)+8.5*x(i+1)+1
    f=-2*(x(i)**3)+12*(x(i)**2)-20*x(i)+8.5 !! dy/dx
    y_eu(i+1) = y_eu(i)+h*f
    f= (f+(-2*(x(i+1)**3)+12*(x(i+1)**2)-20*x(i+1)+8.5))/2
    y_heu(i+1)=y_heu(i)+h*f
    error_eu(i+1)=abs((y_ex(i+1)-y_eu(i+1))/y_ex(i+1))*100
    error_heu(i+1)=abs((y_ex(i+1)-y_heu(i+1))/y_heu(i+1))*100
 end do
 do i=1,9
    write(7,*) x(i),y_ex(i),y_heu(i),error_heu(i),y_eu(i),error_eu(i)
 end do 
end program assignment6







program sum_quiz
  implicit none
  integer i,x,p
  real*8  temp
  x=11
  
  p= mod(x,2)
  if (p.eq.0) then
     temp=1
     do i=0,x-2,2
        temp = temp*(i+2)
        
     end do
     write(6,*) temp

  else
    
     temp = 1
     do i=1,x-2,2
        temp = temp*(i+2)
         write(6,*) temp
     end do
    
  end if
   
end program sum_quiz




  
     
